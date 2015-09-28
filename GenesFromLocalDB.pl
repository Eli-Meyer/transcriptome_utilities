#!/usr/bin/perl
# written by E Meyer, eli.meyer@science.oregonstate.edu
# distributed without any guarantees or restrictions

$scriptname=$0; $scriptname =~ s/.+\///g;

# -- program description and required arguments
unless ($#ARGV == 6)
        {print "\nAssigns gene names to a set of DNA sequences based on sequence similarity\n";
	print "with other genes of known function.  This script relies on both a local\n";
	print "sequence database (formatted for blast) (UniProt is recommended) and a local\n";
	print "definitions file in which each sequence ID is associated with a gene name.\n";
	print "e.g. \"Accession1\tGene name or description\", produced using GetGeneNames.sh\n";
        print "Output:\t a table of best matches and a fasta file of annotated sequences.\n";
        print "Usage:\t $scriptname -i=seqs -b=exclude -n/p=db -a=defs\n";
        print "Arguments:\n";
        print "\t -i=seqs\t fasta file of sequences to be annotated\n";
        print "\t -t=threads\t number of threads to use in blast search\n";
        print "\t -b=exclude\t a text file of \'bad words\' to be excluded (e.g. \'uncharacterized\')\n";
        print "\t -p/n=db1\t The sequence database to search. p=protein, n=nucleotide\n";
        print "\t -a=defs\t Tab delimited file of gene names associated with that db (ID, Name).\n";
        print "\t -e=evalue\t critical e-value for blast search.\n";
	print "\t -o=outopt:\t\"yes\" = initially write blast report to /data/custom_directory, then copy to\n";
	print "\t\t\t working directory and remove /data/custom_directory. (choose this at CGRB)\n";
	print "\t\t\t\"no\"= write to working directory (choose this option on most systems)\n\n";
        print "\n"; exit;
        }

# -- use statements and check for dependencies
$mod1="File::Which";
unless(eval("require $mod1")) {print "$mod1 not found. Exiting\n"; exit;}
use File::Which;
$mod2="Bio::SeqIO";
unless(eval("require $mod2")) {print "$mod2 not found. Exiting\n"; exit;}
use Bio::SeqIO;
$mod3="Bio::SearchIO";
unless(eval("require $mod3")) {print "$mod3 not found. Exiting\n"; exit;}
use Bio::SearchIO;
$dep1 = "blastx";
unless (defined(which($dep1))) {print $dep1, " not found. Exiting.\n"; exit;}

# -- hard coded default settings
my $nprog = "tblastx";
my $pprog = "blastx";
my $hno = 20;			# number of hits to consider
my $crit = 0.0001;		# maximum e-value to consider a true match

# -- data input
my @unk; my @dbs; my @ann; my %out; my %uhs;
my @plist; my $exopt; my $run; my %dhs; my @tlist; my $outdir;
system("date");
foreach $argi (0 .. $#ARGV)
	{$name = $ARGV[$argi]; chomp ($name);
	@flag = split(/=/, $name);
	if ($flag[0] eq "-i") {$qfil = $flag[1];}
	if ($flag[0] eq "-b") {$bfil = $flag[1];}
	if ($flag[0] eq "-p") {$prog = $pprog; $db = $flag[1]}
	if ($flag[0] eq "-n") {$prog = $nprog; $db = $flag[1]}
	if ($flag[0] eq "-a") {$def = $flag[1]}
	if ($flag[0] eq "-t") {$tno = $flag[1]}
	if ($flag[0] eq "-e") {$crit = $flag[1]}
	if ($flag[0] eq "-o") {$outopt = $flag[1]}
	}

my $qseq = new Bio::SeqIO (-file=>$qfil, -format=>'fasta');
system ("cp $qfil tq.fasta");
open (BIN, $bfil); @bword = <BIN>; close(BIN);
print "Avoided terms: "; for (@bword) {chomp ($_); print $_, ", ";}
print "\n\n";

print "Loading sequences...\n";
while (my $seq = $qseq->next_seq) 
	{$nom = $seq->display_id; 
	$dfl = $seq->description;
	$ss = $seq->seq;
	push @unk, $nom;
	$uhs{$nom} = $ss;
	$dhs{$nom} = $dfl;
	}
my %shs = %uhs; 
my @orig = keys(%shs);
print "Done.\n";

# -- load gene name database
print "Loading definitions file...\n";
open (TAB, $def);
my %tabh;
while(<TAB>)
	{
	chomp;
	@cols = split("\t", $_);
	$tabh{$cols[0]} = $cols[1];
	}
system("date");
print "Finished loading definitions file.\n\n";

# -- run BLAST search
print $db." database looks good.  Blasting...\n";
$eset = $crit*10;
if ($outopt =~ /y/)
	{
	@chars = ("A".."Z", "a".."z");
	$string .= $chars[rand @chars] for 1..8;
	print "Writing blast report to temporary directory /data/$string\n";
	system("mkdir /data/$string");
	system("$prog -db $db -query $qfil -evalue $eset -num_descriptions $hno -num_alignments $hno -out /data/$string/out.br -num_threads $tno");
	system("cp /data/$string/out.br .");
	system("rm -rf /data/$string");
	print "Moving blast report to ./out.br\n";
	}
elsif ($outopt =~ /n/)
	{
	print "Writing blast report to ./out.br\n";
	system("$prog -db $db -query $qfil -evalue $eset -num_descriptions $hno -num_alignments $hno -out out.br -num_threads $tno");
	}
system("date");
print "Finished blasting ".$db.".\n\n";

# -- parse out top hits for each query and identify the gene name for that hit
print "Parsing blast report from ".$db."...\n";
my $br = new Bio::SearchIO (-file=>"out.br", format=>'blast');
RESULTS: while (my $result = $br->next_result)
	{$qid = $result->query_accession;
	$hcount = 0; 
	HITS: while (my $hit = $result->next_hit)
		{
		$check = 0;
		if ($hcount>0) {next RESULTS;}
		my $hobs = $hit->significance;
		if ($hobs > $crit) {next RESULTS;}
#		print $qid, "\t", $hit->accession, "\n";
		HSPS: while (my $hsp = $hit->next_hsp)
			{$eobs = $hsp->evalue;
			if ($hcount>0) {next RESULTS;}
			if ($eobs > $crit) {next HITS;}
			if ($eobs <= $crit)
				{$hid = $hit->accession;
				if ($hid =~ /\|.*\|/)
					{
					$hid =~ s/^\w+\|//;
					$hid =~ s/\|.+//;
					}
				if (!exists($tabh{$hid})) {next HITS;}
				$thisname = $tabh{$hid};
				foreach $bw (@bword)
					{
					chomp($bw);
					if ($thisname =~ /$bw/i) {$check++; next HITS;}
					}
				if ($check == 0)
					{
					$thisname =~ s/ /_/g;				
					$out{$qid}{"hit"} = $hid;
					$out{$qid}{"name"} = $thisname;
					$out{$qid}{"e"} = $eobs;
					delete($uhs{$qid});
					$hcount++; next RESULTS;
					}
				}
			}
		}
	}
system("date");
print "Finished parsing report from ".$db.".\n\n";

# -- output both annotated and non-annotated sequences with whatever information
# -- is available for the sequence.  original descriptions are retained.
my $outseqs = new Bio::SeqIO(-file=>">gene_annotated.fasta", -format=>'fasta');
foreach $q (@orig)
	{if (defined($out{$q})) 
		{
		$so = new Bio::Seq(-display_id=>$q,
			-seq=>$shs{$q});
		$od = $dhs{$q};
		$nd = $od." Match_Acc=".$out{$q}{"hit"}." Gene=".$out{$q}{"name"};
		$so->description($nd);
		$outseqs->write_seq($so);
		}
	if (!defined($out{$q})) 
		{
		$so = new Bio::Seq(-display_id=>$q,
			-seq=>$shs{$q},
			-description=>$dhs{$q});
		$outseqs->write_seq($so);
		}
	}

# -- output summary information
my @un = keys (%uhs); my $uno = @un;
print $uno." sequences remained un-annotated.\n";
my @an = keys (%out); my $ano = @an;
print $ano." sequences successfully annotated.\n\n";
foreach $q (keys(%out))
	{
	print $q, "\t";
	print $out{$q}{"hit"}, "\t";
	print $out{$q}{"e"}, "\t";
	print $out{$q}{"name"}, "\t";
	print "\n";}
print "\n";
#system("rm error*.log");
system("rm tq.fasta");
#system("rm tqi.fasta");

system("cp $outdir/out.br .");
system("date");
print "\n";

