#!/usr/bin/perl
# written by E Meyer, eli.meyer@science.oregonstate.edu
# distributed without any guarantees or restrictions

# -- program description and required arguments
$scriptname=$0; $scriptname =~ s/.+\///g;
if ($#ARGV < 4 || $ARGV[0] eq "-h") 
	{
	print "\nSearches DNA sequences against a series of contaminant sets (rRNA, Mt, etc), and\n";
	print "removes sequences matching one or more contaminants as they are identified.\n";
	print "Please note that these searches are conducted in series. If a sequence is eliminated\n";
	print "because it matched a record in the first contaminant file, that sequence is not searched\n";
	print "against subsequent contaminant files.\n";
        print "Output:\t table showing the reason for removing discarded reads, and fasta file of reads that passed\n";
        print "Usage:\t script type= score= reads= contam= table= passed=\n"; 
        print "Arguments:\n";
        print "\t type= \t\t type of blast (type=tblastx or blastn) \n";
        print "\t score= \t critcal bit-score threshold \n";
        print "\t reads= \t reads to search for contaminants, fasta file \n";
        print "\t contam= \t contaminant sequences, entered as contam=label,file \n";
        print "\t\t\t\t where label is a name for that set of sequences (e.g., rRNA)\n";
        print "\t\t\t\t and file is the path to those sequences (e.g., /path/rrna.fasta)\n";
        print "\t\t\t you must enter contam=label,file once for each set of contaminants to search \n";
        print "\t table= \t name for the summary output table \n";
        print "\t passed= \t name for the passed sequences fasta output \n";
        print "\n"; exit;
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
$dep2 = "ExcludeFasta.pl";
unless (defined(which($dep2))) {print $dep2, " not found. Exiting.\n"; exit;}

# -- data input
my @cons = ();
foreach $argi (0 .. $#ARGV)
        {$name = $ARGV[$argi]; chomp ($name);
        @flag = split(/=/, $name);
        if ($flag[0] eq "score") {$critscore = $flag[1];}
        if ($flag[0] eq "type") {$myprog = $flag[1];}
        if ($flag[0] eq "reads") {$rfil = $flag[1];}
        if ($flag[0] eq "contam") {push @cons, $flag[1];}
        if ($flag[0] eq "table") {$tfil = $flag[1];}
        if ($flag[0] eq "passed") {$ofil = $flag[1];}
        }

# -- count the sequences in initial sequence file
my $iseqs = new Bio::SeqIO(-file=>$rfil, -format=>'fasta');
my $iseqn = 0; while ($iseqs->next_seq)	{$iseqn++;}
my $maxhits = $iseqn*3;

# -- blast each set of contaminant sequences against the reads in series
# -- eliminating from the database at each round the sequences found
# -- in the previous round as contaminants 
system("date");
system("cp $rfil tq.fasta");
my %discard; my %cth;
foreach $c (@cons)
	{@tk = keys(%discard); $ntk = @tk;
	if ($ntk>0) 
		{open(TL, ">discard.list");
		for (@tk) {print TL $_, "\n";}
		close(TL);
		system("ExcludeFasta.pl discard.list $rfil >tq.fasta");
		}
	($conlabel, $condb) = split(",", $c);
	$cth{$conlabel} = 0;
	print "Blasting ", $conlabel, " against sequences...\n";
	system("formatdb -i tq.fasta -p F");
	system ("blastall -i $condb -d tq.fasta -p $myprog -a 16 -e 1 -F F -v $maxhits -b $maxhits -o tq.br");
	print "Finished blasting ", $conlabel, " against sequences.\n";
	system("date");
	my $report = new Bio::SearchIO(-file=>"tq.br", -format=>'blast');
	while ($result = $report->next_result)
		{
		while ($hit = $result->next_hit)
			{
			while ($hsp = $hit->next_hsp)
				{chomp($hsp);
				if ($hit->bits < $critscore) {next;}
				$discard{$hit->accession} = $conlabel;
				}
			}
		}
	print "Finished parsing blast report.\n";
	system("date");
	print "\n";
	}

# -- write out discard with contaminant types to a table
open(TOUT, ">$tfil");
@dkl = sort(keys(%discard));
$nallcons = 0;
for (@dkl) 
	{print TOUT $_, "\t", $discard{$_}, "\n";
	$cth{$discard{$_}}++;	
	$nallcons++;
	}

# -- write out sequences that passed to a fasta file
@tk = keys(%discard); $ntk = @tk;
if ($ntk>0) 
	{open(TL, ">discard.list");
	for (@tk) {print TL $_, "\n";}
	close(TL);
	system("invert_seq_ext.pl discard.list $rfil >tq.fasta");
	}
system("cp tq.fasta $ofil");

# -- print out summary of the process
print "\n";
print $iseqn, " sequences in input file\n";
print $nallcons, " sequences look like contaminants\n";
for (keys(%cth))
	{
	print "\t", $_, "\t", $cth{$_}, "\n";
	}
my $oseqs = new Bio::SeqIO(-file=>$ofil, -format=>'fasta');
my $oseqn = 0; while ($oseqs->next_seq)	{$oseqn++;}
print $oseqn, " sequences passed all tests\n";
print "\n";

# -- clean up after yourself!
#my $fdbname = $rfil.".x";
#system("rm $fdbname\*");
system("rm tq.f*");
system("rm error.log") if (-e "error.log");
system("rm formatdb.log") if (-e "formatdb.log");

