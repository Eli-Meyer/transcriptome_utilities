#!/usr/bin/perl
# written by E Meyer, eli.meyer@science.oregonstate.edu
# distributed without any guarantees or restrictions

$scriptname=$0; $scriptname =~ s/.+\///g;

# -- program description and required arguments
unless ($#ARGV == 2)
        {print "\nAssigns GO terms to a set of sequences already annotated with gene names based on UniProt.\n";
        print "Output:\t a fasta file of annotated sequences.\n";
        print "Usage:\t $scriptname input annotations output\n";
        print "Arguments:\n";
        print "\t input\t\tfasta file of sequences to be annotated\n";
        print "\t annotations\tfile of Uniprot GO associations, produced using GOAnnotTable.pl\n";
        print "\t output\t\ta name for the output file\n";
        print "\n"; exit;
        }

my $seqfile = $ARGV[0];
my $dbfile = $ARGV[1];
my $outfile = $ARGV[2];

open(DB, $dbfile);
while(<DB>)
	{
	chomp;
	@cols = split("\t", $_);
	$cols[1] =~ s/ /\,/g;
	$dbh{$cols[0]} = $cols[1];
	}
close(DB);

open(IN, $seqfile);
open(OUT, ">$outfile");
while(<IN>)
	{
	chomp;
	$goi = $mi = $igi = "";
	unless ($_ =~ />/) {print OUT $_, "\n"; next;}
	$_ =~ s/>//;
	@bits = split(" ", $_);
	$iti = $bits[0];
	$acci = "";
	foreach $b (@bits)
		{
		$bcount++;
		if ($b =~ /Match_Acc=/) {$acci = $b; $acci =~ s/Match_Acc=//;}
		}
	if (exists($dbh{$acci})) {$goi = $dbh{$acci};}
	else {$goi = "";}
	print OUT ">";
	foreach $b (@bits)
		{
		if ($b =~ /GOMatch=/) {next;}
		if ($b =~ /GOTerms=/) {next;}
		print OUT $b, " ";
		}	
	if ($goi ne "") {print OUT "GO=", $goi, " ";}
	print OUT "\n";
	}
