#!/usr/bin/perl
# written by E Meyer, eli.meyer@science.oregonstate.edu
# distributed without any guarantees or restrictions

$scriptname=$0; $scriptname =~ s/.+\///g;
if ($#ARGV != 0 || $ARGV[0] eq "-h") 
	{
	print "\nConverts the ontology file from Gene Ontology website (http://www.geneontology.org/ontology/gene_ontology.obo)\n";
	print "into a tab-delimited annotations file suitable for use with GOFromGeneAnnotation.pl\n"; 
	print "Usage: $scriptname input > output\n\n"; 
	exit;
	}

open (IN, $ARGV[0]);

while (<IN>)
	{
	if ($_ =~ /^\!/) {next;}
	chomp;
	@cols = split("\t", $_);
	$gh{$cols[1]}{$cols[4]}++;
	}
close(IN);

foreach $a (sort(keys(%gh)))
	{
	%ah = %{$gh{$a}};
	print $a, "\t";
	foreach $g (sort(keys(%ah)))
		{
		print $g, " ";
		}
	print "\n";
	}
