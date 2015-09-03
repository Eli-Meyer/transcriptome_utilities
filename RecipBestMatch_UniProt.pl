#!/usr/bin/perl
# written by E Meyer, eli.meyer@science.oregonstate.edu
# distributed without any guarantees or restrictions

$scriptname=$0; $scriptname =~ s/.+\///g;
if ($#ARGV != 5 || $ARGV[0] eq "-h") 
	{
	print "\nIdentifies candidates for reciprocal best matches to a list of genes from UniProt\n";
	print "in a set of transcript sequences. Please note -- it is up to the user to examine these\n";
	print "candidate matches and decide whether they really count as the same gene or not.\n";
	print "Usage: $scriptname list db gene_names badwords transcripts bitscore\n"; 
	print "Where:\tlist:\t\ta list of UniProt accession numbers for the genes of interest.\n";
	print "\tdb:\t\tcomplete path to the local BLAST-formatted UniProt database.\n";
	print "\tgene_names:\ttab delimited file associating each accession in the DB with a gene name\n";
	print "\t\t\t(the output from GetGeneNames.sh)\n";
	print "\tbadwords:\ta list of words to exclude; annotations matching these strings considered uninformative and skipped.\n";
	print "\tbitscore:\tcritical (minimum) bit-score required to count as a valid match.\n\n";
	exit;
	}

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
$dep1 = "blastall";
unless (defined(which($dep1))) {print $dep1, " not found. Exiting.\n"; exit;}
$dep2 = "formatdb";
unless (defined(which($dep2))) {print $dep2, " not found. Exiting.\n"; exit;}

my $lfile = $ARGV[0];		# list of protein accession numbers from UniProt	
my $dbfile = $ARGV[1];		# uniprot db, complete path
my $gnfile = $ARGV[2];		# tab delimited file of gene annotations
my $badwords = $ARGV[3];	# file of badwords to exclude
my $efile = $ARGV[4];		# file of cDNA sequences to be searched
my $tnthd = $ARGV[5];		# bitscore threshold for blast
my $cores = 8;			# hard coded number of cores to use for BLAST

# read in list of badwords
open(BWD, $badwords); 
@bwords = ();
while (<BWD>) {chomp; push @bwords, $_;}
close(BWD);
foreach $b (@bwords) {print STDERR $b, " ";} print "\n";

# -- read in list of model protein accessions
print STDERR "Reading list of model proteins ...\n";
open(LF, $lfile);
while(<LF>)
	{
	chomp;
	$lh{$_}++;
	}
close(LF);
print STDERR "Done.\n";

# -- extract those protein sequences from UniProt
print STDERR "Extracting model proteins from $dbfile ...\n";
if (-e "tmp_seqs.fasta") {system("rm tmp_seqs.fasta");}
foreach $l (sort(keys(%lh)))
	{
	print STDERR $l, "...\n";
	system("grep $l $dbfile -m 1 -A 50 >>tmp_seqs.fasta");
	}

my $tseqs = new Bio::SeqIO(-file=>"tmp_seqs.fasta", -format=>"fasta");
my $tq = new Bio::SeqIO(-file=>">tq.fasta", -format=>"fasta");
while($seq = $tseqs->next_seq)
	{
	$si = $seq->display_id;
	$si =~ s/^[a-z]+\|//ig;	
	$si =~ s/\|.+//;
#	print STDERR $si, "\n";
	if(!exists($lh{$si})) {next;}
	if(exists($dh{$si})) {next;}
	$tq->write_seq($seq);
	$dh{$si}++;
	}
print STDERR "Done.\n";

# -- format for blast comparisons
unless (-e "$efile.nhr") 
	{
	print STDERR "Formatting $efile for tblastn ...\n";
	system("formatdb -i $efile -p F");
	print STDERR "Done.\n";
	}

# -- find all significant matches for each protein in the EST dataset using BLASTX
print STDERR "Comparing protein(s) against ESTs...\n";
system("blastall -d $efile -p tblastn -a $cores -v 100 -b 100 -e 1 -F F -i tq.fasta -o qve.br");
my $report = new Bio::SearchIO(-file=>"qve.br", -format=>"blast");
while($result = $report->next_result)
	{
	while ($hit = $result->next_hit)
		{
		$qname = $result->query_accession;
		$qname =~ s/^\w+\|//;
		$qname =~ s/\|.*//;
		while ($hsp = $hit->next_hsp)
			{
			$hname = $hit->accession;
			if ($hsp->bits < $tnthd) {next;}
			unless (exists($hh{$hname})) {$hh{$hname}=$qname;}
			if (exists($qveh{$qname}{$hname})) {next;}
			$qveh{$qname}{$hname}++;
			print STDERR $qname, "\t", $hname, "\n";
			}
		}
	}
print STDERR "Done.\n";

# -- extracting ESTs matched in that search for reciprocal blast
print STDERR "Extracting ESTs for reciprocal BLAST...\n";
my $eseqs = new Bio::SeqIO(-file=>$efile, -format=>"fasta");
my $trseqs = new Bio::SeqIO(-file=>">tmp_rq.fasta", -format=>"fasta");
while($seq = $eseqs->next_seq)
	{
	if (exists($hh{$seq->display_id}))
		{
		$trseqs->write_seq($seq);
		}
	}
print STDERR "Done.\n";

# -- find best annotated match for each in uniprot db using BLASTX
print STDERR "Comparing ESTs against Uniprot...\n";
system("blastall -d $dbfile -p blastx -a $cores -v 100 -b 100 -e 1 -F F -i tmp_rq.fasta -o evq.br");
my $report = new Bio::SearchIO(-file=>"evq.br", -format=>"blast");
while($result = $report->next_result)
	{
	$hitcount = 0;
#	$prevhit = 0;
	$qname = $result->query_accession;
	if (exists($evqh{$qname})) {next;}
	while ($hit = $result->next_hit)
		{
		if ($hitcount>0 && exists($evqh{$qname})) {last;}
		while ($hsp = $hit->next_hsp)
			{
			$thishit = $hsp->bits;
			if ($thishit < $tnthd) {last;}
#			if ($thishit < $prevhit) {last;}
			if ($thishit >= $tnthd) {$hitcount++;}
#			$prevhit = $thishit;
			$hname = $hit->accession;
			$hname =~ s/^\w+\|//;
			$hname =~ s/\|.*//;
			$annoti = `grep $hname $gnfile`;
			chomp($annoti);
			@anna = split("\t", $annoti);
			$badcount = 0;
			foreach $b (@bwords) {if ($anna[1] =~ /$b/i) {$badcount++;}}
			if ($badcount==0)
				{
				$evqh{$qname} = $hname;
				$evqa{$qname} = $anna[1];
				$eh{$hname}++;
				print STDERR $qname, "\t", $hname, "\t", $anna[1], "\t", $badcount, "\n";
				last;
				}
			}
		}
	if ($hitcount>0 && !exists($evqh{$qname}))
		{
		$evqh{$qname} = $hname;
		$evqa{$qname} = $anna[1];
		$eh{$hname}++;
		print STDERR $qname, "\t", $hname, "\t", $anna[1], "\t", $badcount, "\n";
		}
	}
print STDERR "Done.\n";

# -- write output
print STDERR "Writing output...\n";
print "cDNA\tFoundBy\tRecipMatch\tGeneAnnotation\n";
foreach $em (sort(keys(%evqh)))
	{
	print $em, "\t", $hh{$em}, "\t", $evqh{$em}, "\t", $evqa{$em}, "\n";
	}
print STDERR "Done.\n";
