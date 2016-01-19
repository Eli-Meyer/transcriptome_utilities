#!/usr/bin/perl
# written by E Meyer, eli.meyer@science.oregonstate.edu
# distributed without any guarantees or restrictions

# -- check arguments and print usage statement
$scriptname=$0; $scriptname =~ s/.+\///g;
$usage = <<USAGE;
Calculates ortholog hit ratios based on the sum of Hsp lengths
Output:  tab delimited text: sequence, OHR.
Usage:    $scriptname blast_report threshold
Arguments:
         blast_report    blast report from BLASTX against gene models
         threshold       minimum bit score
USAGE
if ($#ARGV != 1 || $ARGV[0] eq "-h") {print "\n", "-"x60, "\n", $scriptname, "\n", $usage, "-"x60, "\n\n"; exit;}

# -- module and executable dependencies
$mod1="Bio::SearchIO";
unless(eval("require $mod1")) {print "$mod1 not found. Exiting\n"; exit;}
use Bio::SearchIO;

# -- build variables from user input
my $br = $ARGV[0];	# blast report from blasting transcripts against 
			# gene models from a sequenced relative
my $thd = $ARGV[1];	# bit-score threshold

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
			$hsl{$qid}{$scount}{"qi"} = $hsp->start('query');
			$hsl{$qid}{$scount}{"qf"} = $hsp->end('query');
			$hsl{$qid}{$scount}{"hi"} = $hsp->start('hit');
			$hsl{$qid}{$scount}{"hf"} = $hsp->end('hit');
			}
		}
	}

# loop through all queries
@qa = sort(keys(%qlh));
foreach $q (@qa)
	{
	if(!exists($qsh{$q})) {next;}
# for each query, sort best-hit hsps by position in the hit
	%subh = %{$hsl{$q}};
	@hspa = sort{$subh{$a}{"hi"} <=> $subh{$b}{"hi"}}(keys(%subh));
	$qleni = $qlh{$q};
	$hleni = $hlh{$qsh{$q}};
# for each hsp, for each position in the reference gene, ask: 
# is this position in an HSP? 
	%stath = ();
	foreach $h (@hspa)
		{
		$ii = $subh{$h}{"hi"};
		$fi = $subh{$h}{"hf"};
		for ($a=$ii; $a<=$fi; $a++)
			{
			$stath{1}{$a} = 1;
			$hith{$qsh{$q}}{$a} = $hleni;
			}
		}
	%oneh = %{$stath{1}};
	@ka = keys(%oneh); $nka = @ka;
	$ohr = $nka / $hleni;
	}

# output each hit and its ohr
foreach $h (sort(keys(%hith)))
	{
	%hhh = %{$hith{$h}};
	@ka = keys(%hhh); $nka = @ka;
	$ohr = $nka / $hhh{$ka[0]};
	print $h, "\t", $ohr, "\n";
	}

