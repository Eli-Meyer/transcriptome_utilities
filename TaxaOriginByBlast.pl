#!/usr/bin/perl
# written by E Meyer, eli.meyer@science.oregonstate.edu
# distributed without guarantees or restrictions

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
$mod4="Bio::DB::Taxonomy";
unless(eval("require $mod4")) {print "$mod4 not found. Exiting\n"; exit;} 
use Bio::DB::Taxonomy;
$dep1 = "blastall";
unless (defined(which($dep1))) {print $dep1, " not found. Exiting.\n"; exit;}
$dep2 = "blastdbcmd";
unless (defined(which($dep2))) {print $dep2, " not found. Exiting.\n"; exit;}

# build variables from user input
my $qfile = $ARGV[0]; 
my $db = $ARGV[1]; 
my $db_dir = $ARGV[2];
my $crit = $ARGV[3]; 
my $ranklev = $ARGV[4]; 
my $blastopt = $ARGV[5];
my $outopt = $ARGV[6];

$scriptname=$0; $scriptname =~ s/.+\///g;
if ($#ARGV != 6 || $ARGV[0] eq "-h")
	{
	print "\nDetermines the most likely source for each sequence in a transcriptome assembly\n";
	print "based on the taxonomic ID of each sequence's best match in a DB containing\n";
	print "sequences from diverse taxa (NCBI's nr)\n";
	print "Usage: $scriptname query db db_dir threshold level blastopt\n";
	print "\tquery:\t\tfasta file of cDNA sequences\n";
	print "\tdb:\t\tpath to local nr database, preformatted with taxonomic information\n";
	print "\tdb_dir:\t\tpath to the directory containing db and taxonomy files from NCBI:\n";
	print "\t\t\t(nodes.dmp, names.dmp, taxadb.bti, taxadb.btd in db_dir,\n";
	print "\t\t\tid2names, names2id, nodes, parents in db_dir/bin)\n";
	print "\tthreshold:\tbit-score threshold for blast matches\n";
	print "\tlevel:\t\ttaxonomic level (kingdom, phylum, class, etc.), as defined in NCBI taxonomy\n";
	print "\tblastopt:\tRun BLAST? first time, enter TRUE. If re-running, enter FALSE\n";
	print "\toutopt:\t\t\"yes\" = initially write blast report to /data/custom_directory, then copy to\n";
	print "\t\t\t working directory and remove /data/custom_directory. (choose this at CGRB)\n";
	print "\t\t\t\"no\"= write to working directory (choose this option on most systems)\n\n";
	exit;
	}

# build hard coded variables
my @chars = ("A".."Z", "a".."z");

$scriptname=$0; $scriptname =~ s/.+\///g;
my $hno = 10;
my $pathi = $db_dir;
#my $pathi = "/nfs/ZOOLOGY/Meyer_Lab/DB/nr";	# path to database directory
# this must contain up to date versions of:
# database files preformatted with taxonomic information from NCBI (e.g. nr, nt)
# taxadb files from NCBI (taxadb.bti, taxadb.btd)
# nodes.dmp and names.dmp from NCBI
# executables from NCBI taxonomy required by bioperl 
# 	(id2names, names2id, nodes, parents) in $pathi/bin

# BLAST search section
if ($blastopt =~ /t/i || $blastopt !~ /f/i) 
	{
# run BLAST search
	print STDERR "\n";
	print STDERR "Blasting $qfile against $db.\n";
	print STDERR `system("date")`;
	if ($outopt =~ /y/)
		{
		$string .= $chars[rand @chars] for 1..8;
		system("mkdir /data/$string");
		system("blastall -d $db -p blastx -a 4 -e $crit -v $hno -b $hno -i $qfile -o /data/$string/blast.br");
		system("cp /data/$string/blast.br .");
		system("rm -rf /data/$string");
		}
	elsif ($outopt =~ /n/)
		{
		system("blastall -d $db -p blastx -a 4 -e $crit -v $hno -b $hno -i $qfile -o blast.br");
		}
	print STDERR "Blasting completed.\n";
	print STDERR `system("date")`;
	}
else	{
	print STDERR "BLAST search already completed. Using existing BLAST reports for next steps\n";
	print STDERR `system("date")`;
	}

# Parsing taxonomic info from best match
print STDERR "Parsing blast report.\n";
print STDERR "Extracting taxonomic identities from top hit records.\n";
print STDERR `system("date")`;
print "query\tmatch\tspecies\ttaxonomic_level_$taxlev\n";

my %th;
@allta = qw{kingdom phylum class family order genus species};
foreach $t (@allta) {$jcount++; $ath{$t}=$jcount; $oth{$jcount}=$t;}

my $report = new Bio::SearchIO (-file=>"blast.br", -format=>'blast');
while ($result = $report->next_result)
	{$status = 0;
	while ($hit = $result->next_hit)
		{
		if ($status > 0) {next;}
		$eobs = $hit->bits;
		$qid = $result->query_accession;
		$sid = $hit->accession;
		if ($eobs >= $crit)
			{

# user input - path to taxonomy files, accession number, level
		$acc = $sid;
#		$level = $taxlev;

# identify species from accession
		system("blastdbcmd -db nr -entry $acc -outfmt %T > tmp.log");	
		$taxonid = `head -n 1 tmp.log`; chomp($taxonid);
		system("rm tmp.log");
#		print $taxonid, "\n";

# extract taxonomy for that species
		my $idx_dir = "$pathi/bin";
		my $new_dir = "./tmpdir";
		unless (-e "./tmpdir")
			{
			system("mkdir ./tmpdir");
			system("cp $idx_dir/* $new_dir/");
			}
		my ($nodefile,$namesfile) = ("$pathi/nodes.dmp","$pathi/names.dmp");
		my $db = new Bio::DB::Taxonomy(-source => 'flatfile',
			-nodesfile => $nodefile, -namesfile => $namesfile,
			-directory => $new_dir);
#		print $idx_dir, "\t", $nodefile, "\t", $namesfile, "\n";
		@ta = ();
		if (defined $db->get_Taxonomy_Node(-taxonid => $taxonid)) # is your species in the database?
			{
		 	my $node = $db->get_Taxonomy_Node(-taxonid => $taxonid);
			my $kingdom = $node;	
			push @ta, $kingdom->scientific_name;
			for (1..25) 
				{
				if(!defined($db->get_Taxonomy_Node(-taxonid => $kingdom->parent_id))) {last;}
				$kingdom = $db->get_Taxonomy_Node(-taxonid => $kingdom->parent_id);
				push @ta, $kingdom->scientific_name;
		    		}
			}

# identify classification at the specified level
		$ranki=$ranklev;
		%trh = (); 
		$nexlev = $oth{$ath{$ranki}+1};
		foreach $t (@ta)
		        {
			if (!defined($db->get_taxon(-name=>$t))) {next;}
		        $taxon = $db->get_taxon(-name=>$t);
		        $tid = $t;
		        $tname = $taxon->scientific_name;
		        $trank = $taxon->rank;
		        $trh{$trank} = $tname;
		        }
		$nata = @allta;
		$mylab = "";
		for ($a=$ath{$ranki};$a<=$nata;$a++)
			{
			if (!exists($trh{$oth{$a}})) {next;}
			if ($a eq $ath{$ranki}) 
				{
				$mylab=$trh{$oth{$a}}; 
#				print $mylab, "\n"; 
				last;
				}
			else
				{
				$mylab="Uknown $ranki, $oth{$a}: $trh{$oth{$a}}";
#				print $mylab, "\n";
				last;
				}
			}

#$ss = "@ta[0]";
#if ($ss eq $exclude) {next;}			
#@rta = reverse(@ta);
#print "@rta\n";				### TESTING
#print $rta[$level-1], "\n";
#print "@ta[0]", "\n";


			print $qid, "\t", $sid, "\t";
#			print $ss, "\t";
#			$seltax = $rta[$level-1];
			$seltax=$mylab;
			print $seltax;
			print "\n";
			if ($th{$seltax}) {$th{$seltax}++;}
			if (!$th{$seltax}) {$th{$seltax} = 1;}
			$status++;
			}
		}
	}
print STDERR "Parsing completed.\n";
print STDERR `system("date")`;

print "\n";
print "Taxon\tNumber\n";
@taxlist = keys(%th); @taxlist = sort(@taxlist);
for $t (@taxlist) {print $t, "\t", $th{$t}, "\n";}
print "\n";

