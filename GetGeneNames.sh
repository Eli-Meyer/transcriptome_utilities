#!/bin/bash

# define variables from input
INPUT=$1

# check input and print usage
if [[ -z $INPUT ]]
	then
	echo ""
	echo "Prepares a tab delimited text file of gene names from the UniProt fasta file";
	echo "Usage: script_name input"; 
	echo "	input:	FASTA file containing all UniProt records (Swiss-Prot + TrEMBL)";
	echo ""
	exit
fi

grep ">" $INPUT > working.tab
perl -pi -e "s/^>[a-z]+\|//" working.tab
perl -pi -e "s/\|\w+ / /" working.tab
perl -pi -e "s/OS=.*//" working.tab
grep -v "ncharacterized protein" working.tab >out1.tab
mv out1.tab gene_names.tab
rm working.tab
perl -pi -e "s/ /\t/" gene_names.tab
