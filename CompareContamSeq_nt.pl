#!/usr/bin/perl
# written by E Meyer, eli.meyer@science.oregonstate.edu
# distributed without any guarantees or restrictions

# -- check arguments and print usage statement
$scriptname=$0; $scriptname =~ s/.+\///g;
$usage = <<USAGE;
Identifies most likely origin of each sequence by comparing against a
nucleotide sequence database from a single close relative, and a database
of likely comtaminants. e.g. for corals, A. digitifera would be a
good choice for a close relative, and Symbiodinium would be a likely contaminant.

Please note that since this is a nucleotide comparison, this approach is only reasonable
if you have sequence information for extremely close relatives (preferably the same species)
for both the target and contaminant. 

Each sequence is assigned to the source which it matches best, or to neither
if it matches neither database.

Usage: $scriptname -q queries -s score -t target_db -c contam_db
Required arguments:
        queries:	FASTA-formatted file of DNA sequences (e.g. transcripts)
        score:		bit-score threshold (matches this high or higher count as a match)
        target_db:	BLAST-formatted nucleotide database of a closely related species
        contam_db:	BLAST-formatted nucleotide database of an expected contaminant
USAGE

# -- module and executable dependencies
# use this block if checking for executable dependencies
# copy the block and edit to check for additional Perl modules required by the script
$mod1="File::Which";
unless(eval("require $mod1")) {print "$mod1 not found. Exiting\n"; exit;}
use File::Which;
$mod1="Getopt::Std";
unless(eval("require $mod1")) {print "$mod1 not found. Exiting\n"; exit;}
use Getopt::Std;
$mod1="Bio::SearchIO";
unless(eval("require $mod1")) {print "$mod1 not found. Exiting\n"; exit;}
use Bio::SearchIO;;
$mod1="Bio::SeqIO";
unless(eval("require $mod1")) {print "$mod1 not found. Exiting\n"; exit;}
use Bio::SeqIO;;

# use this block and edit to check for executables required by the script
$dep1 = "blastn";
unless (defined(which($dep1))) {print $dep1, " not found. Exiting.\n"; exit;}

# get variables from input
getopts('q:s:t:c:h');	# in this example a is required, b is optional, h is help
if (!$opt_q ||!$opt_s ||!$opt_t ||!$opt_c || $opt_h) {print "\n", "-"x60, "\n", $scriptname, "\n", $usage, "-"x60, "\n\n"; exit;}
$qfile = $opt_q;
$score = $opt_s;
$tdb = $opt_t;
$tfname = $tdb; 
$tfname =~ s/.+\///g;
#print $tfname, "\n\n";
$hno = 3;
$cpu = 1;
$cdb = $opt_c;

@dba = ();
push @dba, $tdb;
push @dba, $cdb;

# BLAST search section
# for each database
foreach $d (@dba)
	{
# run blast search
	print "Comparing $qfile against $d...\n";
	$ofn = $d;
	$ofn =~ s/.+\///g;
	$ofni = $ofn;
	$ofn =~ s/\.[a-z]+$//g;
	$ofn = "$ofn.br";
	system("blastn -db $d -query $qfile -out $ofn -num_descriptions $hno -num_alignments $hno -num_threads $cpu");
#	system("mv $d.br .");
	print "Done.\n";

# extract and store information about best matches
#	$d =~ s/.+\///g;
#	$fni = $d.".br"; 
	$report = new Bio::SearchIO(-file=>$ofn, -format=>"blast");
	while ($result = $report->next_result)
		{
		$hitno = 0;
		while ($hit = $result->next_hit)
			{
			if ($hit->bits < $score) {next;}
			$hitno++;
			if ($hitno > 1) {next;}
#			print $result->query_accession, "\t";
#			print $hit->accession, "\t";
#			print $hit->bits, "\n";
			$hh{$result->query_accession}{$d}{"hit"} = $hit->accession;
			$hh{$result->query_accession}{$d}{"score"} = $hit->bits;
			}
		if ($hitno eq 0) 
			{
			
			}
		}
	}

# decisions and output section
# for each query sequence
open(TAB, ">origin_summary.tab");
foreach $qs (sort(keys(%hh)))
	{
	$constat = 0; $conid = "";
# compare best matches across databases
	%qh = %{$hh{$qs}};
	@sda = sort{$qh{$b}{"score"}<=>$qh{$a}{"score"}}(keys(%qh));
	$nsda = @sda;

# decide whether each sequence is a contaminant
	if ($nsda<1) {$constat = "NA";}
	elsif ($sda[0] eq $tdb) {$constat = 0;}
	else
		{
		$constat = 1;
		$conid = $sda[0];
		}
	if ($constat == 0) {$idh{$qs} = "target";}
	else {$idh{$qs} = $conid;}
	print $qs, "\t", "@sda", "\t", $constat, "\t", $idh{$qs}, "\n";

	}

# write out sequences and summary output to the appropriate files
$inseq = new Bio::SeqIO(-file=>$qfile, -format=>"fasta");
while ($seq = $inseq->next_seq)
	{
	$seqcount++;
	$qs = $seq->display_id;
	if(exists($idh{$qs}))
		{
		$qi = $idh{$qs};
		$qi =~ s/.+\///g;
		$qi =~ s/\.\w+$//;
		$destfile = $qi.".screened.fasta";
		print TAB $qs, "\t", $qi, "\n";		
		if ($idh{$qs} eq "target") {$tarcount++;}
		else {$concount++;}
		}
	else	{
		$destfile = "nomatch.screened.fasta";
		print TAB $qs, "\t", "no match", "\n";		
		$unkcount++;
		}
	$outseq = new Bio::SeqIO(-file=>">>$destfile", -format=>"fasta");
	$outseq->write_seq($seq);
	}
print "\n", $seqcount, " sequences input.\n";
print $tarcount, " of these matched ", $tfname, " more closely than any contaminants.\n";
print $concount, " matched contaminants more closely than ", $tfname, ".\n";
print $unkcount, " matched none of the supplied DB (nomatch.screened.fasta).\n\n";
