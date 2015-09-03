#!/usr/bin/perl
# written by E Meyer, eli.meyer@science.oregonstate.edu
# distributed without any guarantees or restrictions

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
$dep1 = "blastx";
unless (defined(which($dep1))) {print $dep1, " not found. Exiting.\n"; exit;}

# take user input
$scriptname=$0; $scriptname =~ s/.+\///g;
unless ($#ARGV>1)
	{
	print "\nIdentifies the most likely origin of each sequence in a transcriptome assembly\n";
	print "by comparison with a protein DB from a single close relative, and one or more \n";
	print "databases of likely contaminants. e.g. for corals, A. digitifera would be a\n";
	print "good target and Symbiodinium would be a likely contaminant.\n";
	print "Each sequence is assigned to the source which it matches best.\n";
	print "\nUsage: $scriptname queries score target_db contam_db1 ... contam_dbN\n";
	print "\tqueries:\tFASTA-formatted file of DNA sequences (e.g. transcripts)\n";
	print "\tscore:\t\tbit-score threshold (matches this high or higher count as a match)\n";
	print "\ttarget_db:\tBLAST-formatted protein database of a closely related species\n";
	print "\tcontam_db1:\tfirst BLAST-formatted protein database of an expected contaminant\n";
	print "\tcontam_dbN:\tlast BLAST-formatted protein database of an expected contaminant\n";
	print "\n";
	exit;
	}
$qfile = $ARGV[0];
$score = $ARGV[1];
$tdb = $ARGV[2];
$tfname = $tdb; 
$tfname =~ s/.+\///g;
#print $tfname, "\n\n";
$hno = 3;
$cpu = 1;

@dba = ();
for ($a=2;$a<=$#ARGV;$a++)
	{
	push @dba, $ARGV[$a];
	}

# BLAST search section
# for each database
foreach $d (@dba)
	{
# run blast search
	print "Comparing $qfile against $d...\n";
	$ofn = $d;
	$ofn =~ s/.+\///g;
	$ofni = $ofn;
	$ofn =~ s/\..+//g;
	$ofn = "$ofn.br";
	system("blastx -db $d -query $qfile -out $ofn -num_descriptions $hno -num_alignments $hno -num_threads $cpu");
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
	if ($nsda<1) {$constat = 0;}
	elsif ($sda[0] eq $tfname) {$constat = 0;}
	else
		{
		$constat = 1;
		$conid = $sda[0];
		}
	if ($constat == 0) {$idh{$qs} = "target";}
	else {$idh{$qs} = $conid;}
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
