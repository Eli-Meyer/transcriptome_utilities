#!/usr/bin/perl
# written by E Meyer, eli.meyer@science.oregonstate.edu
# distributed without any guarantees or restrictions

$scriptname=$0; $scriptname =~ s/.+\///g;
if ($#ARGV != 3 || $ARGV[0] eq "-h") 
	{
	print "\nFinds candidate matches for a gene of interest in an annotated transcriptome assembly.\n";
	print "This script assumes the transcriptome was annotated using GenesFromLocalDB.pl.\n";
	print "Usage: $scriptname list db transcriptome label\n"; 
	print "Where:\tlist:\t\ta list of one or more accession numbers for your gene of interest.\n";
	print "\t\t\t(Note that these must match sequences in your BLAST db)\n";
	print "\tdb:\t\ta protein FASTA file, previously used to annotate the transcriptome (e.g. uniprot.fasta)\n";
	print "\ttranscriptome:\tthe FASTA formatted, annotated transcriptome assembly.\n";
	print "\tlabel:\t\ta short code describing the transcriptome or gene (useful when concatenating \n";
	print "\t\t\toutputs from multiple searches)\n\n";
	exit;
	}
#!/usr/bin/perl
# -- check for dependencies

$mod1="File::Which";
unless(eval("require $mod1")) {print "$mod1 not found. Exiting\n"; exit;} 
use File::Which;
$mod2="Bio::SeqIO";
unless(eval("require $mod2")) {print "$mod2 not found. Exiting\n"; exit;} 
use Bio::SeqIO;
$mod3="Bio::SearchIO";
unless(eval("require $mod3")) {print "$mod3 not found. Exiting\n"; exit;} 
use Bio::SearchIO;
$dep1 = "formatdb";
unless (defined(which($dep1))) {print $dep1, " not found. Exiting.\n"; exit;}
$dep2 = "blastall";
unless (defined(which($dep2))) {print $dep2, " not found. Exiting.\n"; exit;}


unless ($#ARGV==3) {print "usage: script list db ESTs prefix\n"; exit;}

my $lfile = $ARGV[0];		# list of protein accession numbers
my $dbfile = $ARGV[1];		# database containing those proteins
my $efile = $ARGV[2];		# annotated est file to search for matches
my $prefix = $ARGV[3];		# prefix for the output table (name of est source)
my $tnthd = 45;			# bitscore threshold for tblastn

print STDERR "Reading list of model proteins ...\n";
open(LF, $lfile);
while(<LF>)
	{
	chomp;
	$lh{$_}++;
	}
close(LF);
print STDERR "Done.\n";

print STDERR "Extracting model proteins from $dbfile ...\n";
if (-e "tmp_seqs.fasta") {system("rm tmp_seqs.fasta");}
foreach $l (sort(keys(%lh)))
	{
	system("grep $l $dbfile -m 1 -A 50 >>tmp_seqs.fasta");
	}

my $tseqs = new Bio::SeqIO(-file=>"tmp_seqs.fasta", -format=>"fasta");
my $tq = new Bio::SeqIO(-file=>">tq.fasta", -format=>"fasta");
while($seq = $tseqs->next_seq)
	{
	$si = $seq->display_id;
	$si =~ s/^\w+\|//;
	$si =~ s/\|.+//;
	if(!exists($lh{$si})) {next;}
#	print STDERR $si, "\n";
	$tq->write_seq($seq);
	}
print STDERR "Done.\n";

unless (-e "$efile.nhr") 
	{
	print STDERR "Formatting $efile for tblastn ...\n";
	system("formatdb -i $efile -p F");
	print STDERR "Done.\n";
	}

print STDERR "Comparing protein(s) against annotated ESTs...\n";
system("blastall -d $efile -p tblastn -a 1 -v 100 -b 100 -e 1 -i tq.fasta -o qve.br -F F");
my $report = new Bio::SearchIO(-file=>"qve.br", -format=>"blast");
while($result = $report->next_result)
	{
	while ($hit = $result->next_hit)
		{
		while ($hsp = $hit->next_hsp)
			{
			if ($hsp->bits < $tnthd) {next;}
			$hname = $hit->accession;
			$hdesc = $hit->description;
			$qname = $result->query_accession;
#			print $hname, "\t", $hdesc, "\t", $qname, "\n";
#			if ($hdesc !~ /Match_Acc/ || $hdesc !~ /Gene=/) {next;}
			@ha = split(" ", $hdesc);
			foreach $h (@ha)
				{
				if ($h =~ /Match_Acc=/)
					{
					$h =~ s/Match_Acc=//;
					$mai = $h;
					}
				if ($h =~ /Gene=/)
					{
					$h =~ s/Gene=//;
					$gni = $h;
					}
				}
#			$rmh{$mai} = $gni;
#			$hnh{$mai} = $hname;
#			$qnh{$mai} = $qname;
			$rmh{$hname} = $gni;
			$hnh{$hname} = $mai;
			$qnh{$hname} = $qname;
#			print $mai, "\t", $gni, "\t", $hname, "\t", qname, "\n";
			}
		}
	}
print STDERR "Done.\n";

print STDERR "Extracting gene annotation for each reciprocal match...\n";
foreach $ma (sort(keys(%qnh)))
	{
	system("grep \"$ma\" $efile >tmp.txt");
	system("perl -pi -e \"s/- /-/g\" tmp.txt");
	open(IN, "tmp.txt");
	while (<IN>) {chomp; $spatt = $_;} 
	close(IN);
	@bits = split(" ", $spatt);
	$gni = "";
	foreach $b (@bits) 
		{
		if ($b =~ /Gene=/)
			{
			$b =~ s/Gene=//;
			$gni = $b;
			}
		}
#	print $gni, "\n";
	$gah{$ma} = $gni;
	}
print STDERR "Done.\n";

print "Protein\tEST_name\tMatch\tGeneAnnotation\n";
foreach $ma (sort(keys(%rmh)))
	{
	print $prefix, "\t", $qnh{$ma}, "\t", $ma, "\t", $hnh{$ma}, "\t", $gah{$ma}, "\n";;
	}
print STDERR "Done.\n";
system("rm tmp_seqs.fasta; rm tq.fasta; rm qve.br; rm tmp.txt");
