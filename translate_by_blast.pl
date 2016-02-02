#!/usr/bin/perl
# written by E Meyer, eli.meyer@science.oregonstate.edu
# distributed without any guarantees or restrictions

# -- check arguments and print usage statement
$scriptname=$0; $scriptname =~ s/.+\///g;
$usage = <<USAGE;
Translates a series of transcript sequences (DNA) into proteins
in the frame determined from blast matches in a user-specified database.
Output: a fasta file containing the longest ORF for each transcript.
Usage: $scriptname database queries threshold
Where:
        database:       protein sequences, whether already formatted for BLAST or not.
        queries:        nucleotide sequences of transcripts to be translated.
        threshold:      minimum bit-score required to count a hit as significant.
USAGE
if ($#ARGV != 2 || $ARGV[0] eq "-h") {print "\n", "-"x60, "\n", $scriptname, "\n", $usage, "-"x60, "\n\n"; exit;}

# -- module and executable dependencies
# use this block if checking for executable dependencies
# copy the block and edit to check for additional Perl modules required by the script
$mod1="File::Which";
unless(eval("require $mod1")) {print "$mod1 not found. Exiting\n"; exit;}
use File::Which;
$mod2="Bio::SeqIO";
unless(eval("require $mod2")) {print "$mod2 not found. Exiting\n"; exit;}
use Bio::SeqIO;
$mod3="Bio::SearchIO";
unless(eval("require $mod3")) {print "$mod3 not found. Exiting\n"; exit;}
use Bio::SearchIO;

# use this block and edit to check for executables required by the script
$dep1 = "blastall";
unless (defined(which($dep1))) {print $dep1, " not found. Exiting.\n"; exit;}

use Bio::SeqIO;
use Bio::SearchIO;

my $dbfile = $ARGV[0];
my $qfile = $ARGV[1];
my $critscore = $ARGV[2];
my $crite = 1;			# max e value to report
my $cpuno = 1;			# number of processors for blast search
my $hno = 3;			# max number of hits to consider

print STDERR "\nBeginning comparison of translated proteins with BLAST matches.\n";
$datei = `date`;
print STDERR $datei, "\n";

# -- format database if needed
if (glob("$dbfile*.pin"))
	{
	print STDERR "$dbfile is already formatted.\n";
	}
else	{
	print STDERR "Formatting database $dbfile ...\n";
	system("formatdb -i $dbfile -p T");
	}

# -- begin BLAST section
# -- check if BLAST has been completed already
if (-e "out.br") 
	{
	print STDERR "\nBLAST search already completed. (out.br exists)\n";
	print STDERR "If this is not correct, delete out.br and re-run.\n";
	print STDERR "Continuing using the existing out.br ...\n\n";
	}
else
	{
	print STDERR "\n";
	print STDERR "Blasting $qfile against $dbfile.\n";
	print STDERR `date`;

# -- run BLAST search if needed

# -- note that this section complies with the CGRB-specific need to write BLAST output to /data directory
# -- this may or may not be possible or useful on other clusters. If it doesnt work, delete this section 
# -- except the indicated lines and instead define $outfile as "./out.br"
my @chars = ("A".."Z", "a".."z");
$string .= $chars[rand @chars] for 1..8;
$outdir = "/data/$string";
system("mkdir $outdir");
$outfile = "$outdir/out.br";											# keep this
system("blastall -d $dbfile -i $qfile -a $cpuno -e $crite -v $hno -b $hno -p blastx -o /data/$string/out.br");	# keep this too
system("cp $outfile .");
system("rm -rf $outdir");
print STDERR "Done.\n";
print STDERR `date`;
	}

# -- find best frame for each query in BLAST report
my $report = new Bio::SearchIO(-file=>"out.br", -format=>"blast");
while($result = $report->next_result)
	{
	$qid = $result->query_accession;
	while ($hit = $result->next_hit)
		{
		$hid = $hit->accession;
		while ($hsp = $hit->next_hsp)
			{
			if(exists($rh{$qid}{"blasttrans"})) {next;}
			if ($hsp->bits>=$critscore)
				{
				$rh{$qid}{"blasttrans"} = $hsp->query_string;
				$rh{$qid}{"frame"} = $hsp->query->frame;
				$rh{$qid}{"strand"} = $hsp->query->strand;
				$rh{$qid}{"hit"} = $hid;
#				print $qid, "\t", $rh{$qid}{"blasttrans"}, "\t";
#				print $rh{$qid}{"frame"}, "\t";
#				print $rh{$qid}{"strand"}, "\n"; 
				}	
			}
		}
	}

# -- translate each query with a BLAST match into a protein, in the frame determined from BLAST
open(TO1, ">all_trans_alignments.aln");
open(TO2, ">all_hsp_alignments.aln");

my $dnaobj = new Bio::SeqIO(-file=>$qfile, -format=>"fasta");
my $proobj = new Bio::SeqIO(-file=>">$qfile.proteins.fasta", -format=>"fasta");
while($seq = $dnaobj->next_seq)
	{
	$sid = $seq->display_id;
	$slh{$sid}++;
	if(!exists($rh{$sid}{"strand"})) {next;}
	if ($rh{$sid}{"strand"} eq -1) {$tseq = $seq->revcom;}
	else	{$tseq = $seq;}
	$pobj = $tseq->translate(-frame => $rh{$sid}{"frame"});
        $pstr = $pobj->seq;
	@poa = split(/\*/, $pstr);
	%orfh = ();
	foreach $p (@poa) {$orfh{$p} = length($p);}
	@spoa = sort{$orfh{$b}<=>$orfh{$a}}(keys(%orfh));
	$pobj->seq($spoa[0]);
	$proobj->write_seq($pobj);
	$th{$sid} = $pobj->seq;
#	print $sid, "\n", $pobj->seq, "\n";
	}

# --clean up after yourself
print STDERR "\nCompleted.\n";
$datei = `date`;
print STDERR $datei, "\n\n";

