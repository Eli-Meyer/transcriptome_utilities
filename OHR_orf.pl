#!/usr/bin/perl
# written by E Meyer, eli.meyer@science.oregonstate.edu
# distributed without any guarantees or restrictions

# -- check arguments and print usage statement
$scriptname=$0; $scriptname =~ s/.+\///g;
$usage = <<USAGE;
Calculates ortholog hit ratios based on longest ORF in BLAST frame
Output:  tab delimited text: sequence, OHR.
Usage: $scriptname blast_report threshold sequences
Arguments:
         blast_report    blast report from BLASTX against gene models
         threshold       minimum bit score
         sequences       fasta file of nucleotide sequence (BLAST queries)
USAGE
if ($#ARGV != 2 || $ARGV[0] eq "-h") {print "\n", "-"x60, "\n", $scriptname, "\n", $usage, "-"x60, "\n\n"; exit;}

# -- module and executable dependencies
# use this block if checking for executable dependencies
# copy the block and edit to check for additional Perl modules required by the script
$mod1="Bio::SearchIO";
unless(eval("require $mod1")) {print "$mod1 not found. Exiting\n"; exit;}
use Bio::SearchIO;
$mod2="Bio::SeqIO";
unless(eval("require $mod2")) {print "$mod2 not found. Exiting\n"; exit;}
use Bio::SeqIO;


# define variables from user input
my $br = $ARGV[0];	# blast report from blasting transcripts against 
			# gene models from a sequenced relative
my $thd = $ARGV[1];	# bit-score threshold
my $seqfile = $ARGV[2];	# fasta file of queries

# read in BLAST report and build a hash of query, hsps, ranges
# and a hash of query lengths
my $report = new Bio::SearchIO(-file=>$br, -format=>"blast");
while ($result = $report->next_result)
	{
	$qid = $result->query_accession;
	$hcount = 0;
	$qlh{$qid} = $result->query_length;
	while ($hit = $result->next_hit)
		{
		$scount = 0; $initstd = 0;
		if ($hit->bits < $thd) {next;}
		$hcount++; if ($hcount>1) {next;}
		$hid = $hit->accession;
		while ($hsp = $hit->next_hsp)
			{
			$scount++;
			if ($scount==1)
				{
				$initstd = $hsp->strand('query');
				}
			else
				{
				if ($hsp->strand('query') != $initstd) {next;}
				}
			$qsh{$qid} = $hid;
# for each position in the hsp
			for ($a=$hsp->start('hit'); $a<=$hsp->end('hit'); $a++)
				{
# record the position in hit coordinates (not query)
				$hsph{$qid}{$a} = $scount;
				}
			$hlh{$hid} = $hit->length;
			$frame = $hsp->query->frame;
			$strand = $hsp->query->strand;
			$string = $hsp->query_string;
			if(exists($tsh{$qid}{$hid})) {next;}
			$tsh{$qid}{"hit"} = $hid;
			$tsh{$qid}{"strand"} = $strand;
			$tsh{$qid}{"frame"} = $frame;
			$tsh{$qid}{"string"} = $string;
			}
		}
	}
# translate into blast-defined frame, get longest ORF, and print out OHR
$seqs = new Bio::SeqIO(-file=>$seqfile, -format=>"fasta");
while ($seq = $seqs->next_seq)
	{

	$sid = $seq->display_id;
	if(!exists($qsh{$sid})) {next;}

        if ($tsh{$sid}{"strand"} eq -1) {$tseq = $seq->revcom();}
        else    {$tseq = $seq;}
        $pobj = $tseq->translate(-frame => $tsh{$sid}{"frame"});
        $pstr = $pobj->seq;
	@poa = split(/\*/, $pstr);
	%orfh = ();
	foreach $p (@poa) {$orfh{$p} = length($p);}
	@spoa = sort{$orfh{$b}<=>$orfh{$a}}(keys(%orfh));
	$maxleni = $orfh{$spoa[0]};
	$hleni = $hlh{$qsh{$sid}};
	$ohri = $maxleni/$hleni;
	print $sid, "\t", $qsh{$sid}, "\t", $maxleni, "\t", $hleni, "\t", $ohri, "\n";
	}

