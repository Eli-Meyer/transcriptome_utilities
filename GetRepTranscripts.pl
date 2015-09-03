#!/usr/bin/perl
# written by E Meyer, eli.meyer@science.oregonstate.edu
# distributed without any guarantees or restrictions

# -- check for dependencies and arguments
$mod1="Bio::SeqIO";
unless(eval("require $mod1")) {print "$mod1 not found. Exiting\n"; exit;} 
use Bio::SeqIO;

$scriptname=$0; $scriptname =~ s/.+\///g;
unless ($#ARGV==3)
        {
        print "\nSelects the longest representative for each component or\n";
	print "subcomponent in a Trinity transcriptome assembly.\n";
        print "Usage: $scriptname assembly.fasta option output.tab output.fasta\n";
	print "Where:\n";
	print "\tassembly.fasta:\tthe input file, assembled by Trinity\n";
	print "\toption:\t\tc for component, or s for subcomponent\n";
	print "\toutput.tab:\ta name for the output summary file\n";
	print "\toutput.fasta:\ta name for the output file of representative transcripts\n\n";
        exit;
        }

# define variables
$inseq = $ARGV[0];
$opt = $ARGV[1];
$outtab = $ARGV[2];
$outseq = $ARGV[3];
$iseqs = new Bio::SeqIO(-file=>$inseq, -format=>"fasta");
$oseqs = new Bio::SeqIO(-file=>">$outseq", -format=>"fasta");
open(OUT, ">$outtab");

# record length, component, and subcomponent for each transcript
while ($seq = $iseqs->next_seq)
	{
	$noseq++;
	$rsid = $seq->display_id;
	@rsa = split("_", $rsid);
	$nh{$rsid}{"c"} = $rsa[0];
	$nh{$rsid}{"s"} = $rsa[0]."_".$rsa[1];
	$clh{$rsa[0]}{$rsid} = $seq->length;
	$slh{$rsa[0]."_".$rsa[1]}{$rsid} = $seq->length;
	$nrch{$rsa[0]}++;
	$nrsh{$rsa[0]."_".$rsa[1]}++;
	}

# build a hash of selected transcripts
if ($opt eq "c")
	{
	print OUT "comp\tmembers\tlongest\n";
	foreach $c (sort(keys(%clh)))
		{
		%subch = %{$clh{$c}};
		@ta = sort{$subch{$b}<=>$subch{$a}}(keys(%subch));
		$nmc = @ta;
		$ti = $ta[0];
		$gh{$ti}++;
		print OUT $c, "\t", $nmc, "\t", $ti, "\n";
		}
	}
elsif ($opt eq "s")
	{
	print OUT "subcomp\tmembers\tlongest\n";
	foreach $s (sort(keys(%slh)))
		{
		%subsh = %{$slh{$s}};
		@ta = sort{$subsh{$b}<=>$subsh{$a}}(keys(%subsh));
		$nms = @ta;
		$ti = $ta[0];
		$gh{$ti}++;
		print OUT $s, "\t", $nms, "\t", $ti, "\n";
		}
	}

# write out selected transcripts to output
$iseqs = new Bio::SeqIO(-file=>$inseq, -format=>"fasta");
while ($seq = $iseqs->next_seq)
	{
	$rsid = $seq->display_id;
	if(exists($gh{$rsid}))
		{
		$oseqs->write_seq($seq);
		}
	}

# print summary output
@nrca = keys(%nrch);
@nrsa = keys(%nrsh);
$nnrc = @nrca;
$nnrs = @nrsa;
print STDERR $noseq, " sequences in input file.\n";
print STDERR $nnrc, " components.\n";
print STDERR $nnrs, " subcomponents.\n";

