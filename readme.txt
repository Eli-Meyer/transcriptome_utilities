-------------------------
CompareContamSeq_nt.pl
-------------------------

------------------------------------------------------------
CompareContamSeq_nt.pl
Identifies most likely origin of each sequence by comparing against a
nucleotide sequence database from a single close relative, and a database
of likely comtaminants. e.g. for corals, A. digitifera would be a
good choice for a close relative, and Symbiodinium would be a likely contaminant.

Please note that since this is a nucleotide comparison, this approach is only reasonable
if you have sequence information for extremely close relatives (preferably the same species)
for both the target and contaminant. 

Each sequence is assigned to the source which it matches best, or to neither
if it matches neither database.

Usage: CompareContamSeq_nt.pl -q queries -s score -t target_db -c contam_db
Required arguments:
        queries:	FASTA-formatted file of DNA sequences (e.g. transcripts)
        score:		bit-score threshold (matches this high or higher count as a match)
        target_db:	BLAST-formatted nucleotide database of a closely related species
        contam_db:	BLAST-formatted nucleotide database of an expected contaminant
------------------------------------------------------------

-------------------------
CompareContamSeq.pl
-------------------------

------------------------------------------------------------
CompareContamSeq.pl
Identifies most likely origin of each sequence by comparing against a
protein sequence database from a single close relative, and a protein database
from a likely comtaminants. e.g. for corals, A. digitifera would be a
good choice for a close relative, and Symbiodinium would be a likely contaminant.

Each sequence is assigned to the source which it matches best, or to neither
if it matches neither database.

Usage: CompareContamSeq.pl -q queries -s score -t target_db -c contam_db
Required arguments:
        queries:	FASTA-formatted file of DNA sequences (e.g. transcripts)
        score:		bit-score threshold (matches this high or higher count as a match)
        target_db:	BLAST-formatted protein database of a closely related species
        contam_db:	BLAST-formatted protein database of an expected contaminant
------------------------------------------------------------

-------------------------
GeneMiner.pl
-------------------------

Finds candidate matches for a gene of interest in an annotated transcriptome assembly.
This script assumes the transcriptome was annotated using GenesFromLocalDB.pl.
Usage: GeneMiner.pl list db transcriptome label
Where:	list:		a list of one or more accession numbers for your gene of interest.
			(Note that these must match sequences in your BLAST db)
	db:		a protein FASTA file, previously used to annotate the transcriptome (e.g. uniprot.fasta)
	transcriptome:	the FASTA formatted, annotated transcriptome assembly.
	label:		a short code describing the transcriptome or gene (useful when concatenating 
			outputs from multiple searches)

-------------------------
GenesFromLocalDB.pl
-------------------------

Assigns gene names to a set of DNA sequences based on sequence similarity
with other genes of known function.  This script relies on both a local
sequence database (formatted for blast) (UniProt is recommended) and a local
definitions file in which each sequence ID is associated with a gene name.
e.g. "Accession1	Gene name or description", produced using GetGeneNames.sh
Output:	 a table of best matches and a fasta file of annotated sequences.
Usage:	 GenesFromLocalDB.pl -i=seqs -b=exclude -n/p=db -a=defs
Arguments:
	 -i=seqs	 fasta file of sequences to be annotated
	 -t=threads	 number of threads to use in blast search
	 -b=exclude	 a text file of 'bad words' to be excluded (e.g. 'uncharacterized')
	 -p/n=db1	 The sequence database to search. p=protein, n=nucleotide
	 -a=defs	 Tab delimited file of gene names associated with that db (ID, Name).
	 -e=evalue	 critical e-value for blast search.
	 -o=outopt:	"yes" = initially write blast report to /data/custom_directory, then copy to
			 working directory and remove /data/custom_directory. (choose this at CGRB)
			"no"= write to working directory (choose this option on most systems)


-------------------------
GetRepTranscripts.pl
-------------------------

Selects the longest representative for each component or
subcomponent in a Trinity transcriptome assembly.
Usage: GetRepTranscripts.pl assembly.fasta option output.tab output.fasta
Where:
	assembly.fasta:	the input file, assembled by Trinity
	option:		c for component, or s for subcomponent
	output.tab:	a name for the output summary file
	output.fasta:	a name for the output file of representative transcripts

-------------------------
GOAnnotTable.pl
-------------------------

Converts the ontology file from Gene Ontology website (http://www.geneontology.org/ontology/gene_ontology.obo)
into a tab-delimited annotations file suitable for use with GOFromGeneAnnotation.pl
Usage: GOAnnotTable.pl input > output

-------------------------
GOFromGeneAnnotation.pl
-------------------------

Assigns GO terms to a set of sequences already annotated with gene names based on UniProt.
Output:	 a fasta file of annotated sequences.
Usage:	 GOFromGeneAnnotation.pl input annotations output
Arguments:
	 input		fasta file of sequences to be annotated
	 annotations	file of Uniprot GO associations, produced using GOAnnotTable.pl
	 output		a name for the output file

-------------------------
OHR_hits.pl
-------------------------

------------------------------------------------------------
OHR_hits.pl
Calculates ortholog hit ratios based on the sum of Hsp lengths
Output:  tab delimited text: sequence, OHR.
Usage:    OHR_hits.pl blast_report threshold
Arguments:
         blast_report    blast report from BLASTX against gene models
         threshold       minimum bit score
------------------------------------------------------------

-------------------------
OHR_orf.pl
-------------------------

------------------------------------------------------------
OHR_orf.pl
Calculates ortholog hit ratios based on longest ORF in BLAST frame
Output:  tab delimited text: sequence, OHR.
Usage: OHR_orf.pl blast_report threshold sequences
Arguments:
         blast_report    blast report from BLASTX against gene models
         threshold       minimum bit score
         sequences       fasta file of nucleotide sequence (BLAST queries)
------------------------------------------------------------

-------------------------
RecipBestMatch_UniProt.pl
-------------------------

Identifies candidates for reciprocal best matches to a list of genes from UniProt
in a set of transcript sequences. Please note -- it is up to the user to examine these
candidate matches and decide whether they really count as the same gene or not.
Usage: RecipBestMatch_UniProt.pl list db gene_names badwords transcripts bitscore
Where:	list:		a list of UniProt accession numbers for the genes of interest.
	db:		complete path to the local BLAST-formatted UniProt database.
	gene_names:	tab delimited file associating each accession in the DB with a gene name
			(the output from GetGeneNames.sh)
	badwords:	a list of words to exclude; annotations matching these strings considered uninformative and skipped.
	bitscore:	critical (minimum) bit-score required to count as a valid match.

-------------------------
RemoveContamSeq.pl
-------------------------

Searches DNA sequences against a series of contaminant sets (rRNA, Mt, etc), and
removes sequences matching one or more contaminants as they are identified.
Please note that these searches are conducted in series. If a sequence is eliminated
because it matched a record in the first contaminant file, that sequence is not searched
against subsequent contaminant files.
Output:	 table showing the reason for removing discarded reads, and fasta file of reads that passed
Usage:	 script type= score= reads= contam= table= passed=
Arguments:
	 type= 		 type of blast (type=tblastx or blastn) 
	 score= 	 critcal bit-score threshold 
	 reads= 	 reads to search for contaminants, fasta file 
	 contam= 	 contaminant sequences, entered as contam=label,file 
				 where label is a name for that set of sequences (e.g., rRNA)
				 and file is the path to those sequences (e.g., /path/rrna.fasta)
			 you must enter contam=label,file once for each set of contaminants to search 
	 table= 	 name for the summary output table 
	 passed= 	 name for the passed sequences fasta output 

-------------------------
TaxaOriginByBlast.pl
-------------------------

Determines the most likely source for each sequence in a transcriptome assembly
based on the taxonomic ID of each sequence's best match in a DB containing
sequences from diverse taxa (NCBI's nr)
Usage: TaxaOriginByBlast.pl query db db_dir threshold level blastopt
	query:		fasta file of cDNA sequences
	db:		path to local nr database, preformatted with taxonomic information
	db_dir:		path to the directory containing db and taxonomy files from NCBI:
			(nodes.dmp, names.dmp, taxadb.bti, taxadb.btd in db_dir,
			id2names, names2id, nodes, parents in db_dir/bin)
	threshold:	bit-score threshold for blast matches
	level:		taxonomic level (kingdom, phylum, class, etc.), as defined in NCBI taxonomy
	blastopt:	Run BLAST? first time, enter TRUE. If re-running, enter FALSE
	outopt:		"yes" = initially write blast report to /data/custom_directory, then copy to
			 working directory and remove /data/custom_directory. (choose this at CGRB)
			"no"= write to working directory (choose this option on most systems)

-------------------------
TranslateByBlast.pl
-------------------------

------------------------------------------------------------
TranslateByBlast.pl
Translates a series of transcript sequences (DNA) into proteins
in the frame determined from blast matches in a user-specified database.
Output: a fasta file containing the longest ORF for each transcript.
Usage: TranslateByBlast.pl database queries threshold
Where:
        database:       protein sequences, whether already formatted for BLAST or not.
        queries:        nucleotide sequences of transcripts to be translated.
        threshold:      minimum bit-score required to count a hit as significant.
------------------------------------------------------------

-------------------------
trinity_reps.pl
-------------------------

Selects the longest representative for each component or
subcomponent in a Trinity transcriptome assembly.
Usage: trinity_reps.pl assembly.fasta option output.tab output.fasta
Where:
	assembly.fasta:	the input file, assembled by Trinity
	option:		c for component, or s for subcomponent
	output.tab:	a name for the output summary file
	output.fasta:	a name for the output file of representative transcripts

-------------------------
GetGeneNames.sh
-------------------------

Prepares a tab delimited text file of gene names from the UniProt fasta file
Usage: script_name input
	input:	FASTA file containing all UniProt records (Swiss-Prot + TrEMBL)

