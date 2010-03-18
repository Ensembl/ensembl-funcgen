#!/usr/local/bin/perl

=head1 NAME

peaks_report.pl -- generate R plots for generic properties of feature sets

=head1 SYNOPSIS

peaks_report.pl [-dbhost ($DB_HOST) -dbuser ($DB_USER) -dbport ($DB_PORT) -dbpass ("") -dbname ($DB_NAME) -species -R -nodump -feature_sets fsetA,fsetB,... ]

=head1 DESCRIPTION

take a set of hits in bed files and generate R plots of the lengths

=head1 OPTIONS

=over

=item B<help>

Give short help

=item B<species>

Mandatory: Species of the database (e.g. homo_sapiens)

=item B<dbhost>

Host where the database is (defaults to $DB_HOST)

=item B<dbuser>

User of the database (defaults to $DB_READ_USER)

=item B<dbpass>

Password for the database user (defaults to "")

=item B<dbport>

Port of the host where the database is (defaults to $DB_PORT)

=item B<dbname>

Name of the database (defaults to $DB_NAME)

=item B<dbhost>

Host of the specific core database to use (defaults to $DNADB_HOST)

=item B<dbuser>

User of the specific core database (defaults to $DNADB_USER)

=item B<dbpass>

Password for the specific core database user (defaults to "")

=item B<dbport>

Port of the host where the specific core database to use is (defaults to $DNADB_PORT)

=item B<dbname>

Name of the specific core database to use

=item B<R>

When specified runs the R code to generate plots of the data.

=item B<nodump>

When specified it will skip dumping the data (e.g. for reusing old dumps)

=item B<compare>

Will generate graphs comparing all datasets (if not specified, only set-specific graphs will be generated)

=item B<regstats>

When specified it generate and plot statistics regarding the regulatory build present in the database
Requires the presence of a table rf_stats in the database

=item B<feature_sets>

List of space separated FeatureSet names to be considered. Names with spaces must be quoted.
When not specified all features are considered

=item B<-all_seq_regions>

Processes all available seq_regions, by default only uses chromosomes.

=back


=head1 SEE ALSO

ensembl-functgenomics/scripts/environments/peaks.env PeaksReport function


=head1 LICENSE

  Copyright (c) 1999-2010 The European Bioinformatics Institute and
  Genome Research Limited.  All rights reserved.

  This software is distributed under a modified Apache license.
  For license details, please see

    http://www.ensembl.org/info/about/code_licence.html


=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <ensembl-dev@ebi.ac.uk>.

  Questions may also be sent to the Ensembl help desk at
  <helpdesk@ensembl.org>.



=cut


# To do
# 1 Move reg stats code in here and allow reg_stats mode only
# 2 Tighten up parameter checking
# 3 Change -feature_sets all sets behaviour?
# 4 Add outdir
# 5 Catch absent files (sugest removal of -no_dump?)
# 6 Handle output better lsf/R out? Use optional run name for file naming?
# 7 Use RMySQL to pull the data directly into R rather than dumping
# 8 Move to sub dir
# 9 DONE Restrict plots to main chromosomes (do NT contigs in a seaprate plot)
# 10 Show all chr names in plot axis
# 11 DONE Correct sql to prevent pruoduct with nr seq_region entries

use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor;
use Data::Dumper;

my ($species, $dnadbhost, $dnadbuser, $dnadbport, $dnadbpass, $dnadbname, $host, $user, $port, $pass, $dbname, $help, $R, $nodump, $compare, $regstats, $all_seq_regions);

$host = $ENV{DB_HOST};
$port = $ENV{DB_PORT};
$user = $ENV{DB_READ_USER};
$pass='';
$dbname = $ENV{DB_NAME};
$dnadbhost = $ENV{DNADB_HOST};
$dnadbport = $ENV{DNADB_PORT};
$dnadbuser = $ENV{DNADB_USER};

#get command line options

my (@fset_names);

print "peaks_report.pl @ARGV\n";

GetOptions (
			'species=s'          => \$species,
			'dnadbhost=s'        => \$dnadbhost,
			'dnadbuser=s'        => \$dnadbuser,
			'dnadbport=i'        => \$dnadbport,
			'dnadbpass=s'        => \$dnadbpass,
			'dnadbname=s'        => \$dnadbname,
			'dbhost=s'           => \$host,
			'dbuser=s'           => \$user,
			'dbport=i'           => \$port,
			'dbpass=s'           => \$pass,
			'dbname=s'           => \$dbname,
			"help|h"             => \$help,
			"R"                  => \$R,
			"nodump"             => \$nodump,
			"compare"            => \$compare,
			"regstats"           => \$regstats,
			"all_seq_regions"    => \$all_seq_regions,
			"feature_sets=s{,}"  => \@fset_names,
		   )  or pod2usage( -exitval => 1 ); #Catch unknown opts

pod2usage(1) if ($help);


# TO DO Need to validate vars here
# and fail nicely



#Check database connections
my $coredba;
if($dnadbname){
  my $coredba = Bio::EnsEMBL::DBSQL::DBAdaptor->new
    (
     -host => $dnadbhost,
     -port => $dnadbport,
     -user => $dnadbuser,
     -pass => $dnadbpass,
     -dbname => $dnadbname,
     -species => $species,
     -group   => 'core',
    );
}

my $efgdba = Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor->new
  (
   -host    => $host,
   -port    => $port,
   -user    => $user,
   -pass    => $pass,
   -dbname  => $dbname,
   -species => $species,
   -group   => 'funcgen',
   -dnadb   => $coredba, #Assumes that new will accept undef as parameter for this...
  );

my $fsa = $efgdba->get_FeatureSetAdaptor();

if(scalar(@fset_names)==0){
   # @fset_names = "all feature sets"
   #this returns the result of a generic query, which is a list inside a list... so the actual list is in [0]
   #TODO -> CHANGE THE API code?
   #TODO pass type as parameter??
   my @fsets = $fsa->fetch_all_by_type('annotated');
   foreach my $fset (@{$fsets[0]}){ push(@fset_names, $fset->name);  } 
}

#print the data of each set to individual files (maybe put all in one file??)
#give other saving options (e.g. clean-up, backup?)


my @sr_types = ('chromosome');
push @sr_types, 'non_chromosome' if $all_seq_regions;


my %sr_type_clauses = (
					'chromosome' => , "='chromosome'",
				   'non_chromosome' => , "!='chromosome'",
				   );

if(!$nodump){
  print "::Dumping Datasets\n";
  foreach my $fset (@fset_names){

	#This was not accounting for nr sr_ids
	foreach my $sr_type(@sr_types){

	  my $query ="SELECT s.name as 'region', (f.seq_region_end - f.seq_region_start) as 'length' FROM annotated_feature f, (select distinct(seq_region_id), sr.name from seq_region sr, coord_system cs where sr.coord_system_id=cs.coord_system_id and cs.name".$sr_type_clauses{$sr_type}.") s, feature_set fs WHERE f.feature_set_id=fs.feature_set_id AND f.seq_region_id=s.seq_region_id AND fs.name='$fset';";

	  #warn $query;

	  my $cmd = "mysql -e \"".$query."\" -quick -h$host -P$port -u$user ".(($pass)? "-p$pass" : "")." $dbname > ${fset}_data.${sr_type}.txt";
	  #print $cmd."\n";
	  system($cmd);
	}
  }

  if($regstats){
    my $query = "SELECT * from rf_stats";
    system("mysql -e \"".$query."\" -quick -h$host -P$port -u$user ".(($pass)? "-p$pass" : "")." $dbname > regstats.txt");
  }


}

#Todo add names without spaces and strange characters...
#print R file with analysis
#Check which parameters one might want...
if (defined $R) {
  print "::Generating the R plots\n";
  open(FO,">peaks_report.R");
  print FO "require(graphics)\n";
  print FO "pdf(file=\"peaks_report.pdf\")\n";
  print FO "par(cex=0.7,cex.main=1.5)\n";

  if($regstats){ 
     print FO "require(gplots);\n";
     print FO "regstats <- read.table(\"regstats.txt\",row.names=1,header=TRUE);\n"; 
     print FO "textplot(regstats,halign=\"left\")\n";
  }

  my %safenames;
  foreach my $fset (@fset_names){
    my $safename = $fset;
    $safename =~ s/\W/_/g;
    $safenames{$fset}=$safename;
  }
  
  foreach my $fset (@fset_names){
    #This will not work in R if the name has characters like ':'
    my $safename = $safenames{$fset};

	foreach my $sr_type(@sr_types){
	  my $data_name = "${safename}_${sr_type}";
		
	  print FO "$data_name <- read.table(\"${fset}_data..${sr_type}.txt\",header=TRUE)\n";
	  #the if is to avoid errors printing empty sets...
	  print FO "if (length(".$data_name."\$length) > 0) barplot(unlist(lapply(split(".$data_name.",".$data_name."\$region),function(x) length(x\$length))),main=\"".$fset."\",xlab=\"region\",ylab=\"number of peaks\")\n";
	  print FO "if (length(".$data_name."\$length) > 0) boxplot(lapply(split(".$data_name.",".$data_name."\$region),function(x) x\$length),main=\"".$fset."\",xlab=\"region\",ylab=\"Peak length\")\n";
	}
  }

  if($compare){ 
	#Maybe do some more analysis here?
	my @safe; 
	
	#For brevity only compare true chromosomes
	map { push @safe, $safenames{$_}.'_chromosome' } @fset_names;
	
	print FO "barplot(c(length(".join("\$length),length(",@safe)."\$length)),main=\"Comparison of Annotated Features\",xlab=\"Feature set\",ylab=\"number of peaks\",col=rainbow(".scalar(@safe)."))\n";
	print FO "legend(\"topright\",legend=c(\"".join("\",\"",@fset_names)."\"),fill=rainbow(".scalar(@fset_names)."))\n";
	print FO "boxplot(".join("\$length,",@safe)."\$length,main=\"Comparison of Annotated Features\",xlab=\"Feature set\",ylab=\"Peak length\",col=rainbow(".scalar(@safe)."))\n";
	print FO "legend(\"topright\",legend=c(\"".join("\",\"",@fset_names)."\"),fill=rainbow(".scalar(@fset_names)."))\n";
	
  }


  close FO;
  #This submits to yesterday by default
  system "R CMD BATCH --slave peaks_report.R";
}

__END__




