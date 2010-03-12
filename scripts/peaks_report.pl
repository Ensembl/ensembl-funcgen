#!/usr/local/bin/perl

=head1 NAME

peaks_report.pl -- generate R plots for generic properties of feature sets

=head1 SYNOPSIS

peaks_report.pl -species [-dbhost ($DB_HOST) -dbuser ($DB_USER) -dbport ($DB_PORT) -dbpass ("") -dbname ($DB_NAME) -R -nodump -feature_sets fsetA,fsetB,... ]

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

Set of feature names to be considered  (separated by ',')
When not specified all features are considered

=back

=head1 LICENCE

This code is distributed under an Apache style licence. Please see
http://www.ensembl.org/info/about/code_licence.html for details.

=head1 AUTHOR

Ensembl Functional Genomics

=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <ensembl-dev@ebi.ac.uk>.

  Questions may also be sent to the Ensembl help desk at
  <helpdesk@ensembl.org>.

=cut

use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor;
use Data::Dumper;

my ($species, $dnadbhost, $dnadbuser, $dnadbport, $dnadbpass, $dnadbname, $host, $user, $port, $pass, $dbname, $help, $R, $nodump, $compare, $regstats);

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
	    "feature_sets=s{,}"  => \@fset_names,
             );


pod2usage(1) if ($help);

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
if(!$nodump){
  print "::Dumping Datasets\n";
  foreach my $fset (@fset_names){
    my $query ="SELECT s.name as 'region', (f.seq_region_end - f.seq_region_start) as 'length' FROM annotated_feature f, seq_region s, feature_set fs WHERE f.feature_set_id=fs.feature_set_id AND f.seq_region_id=s.seq_region_id AND fs.name='$fset';";
    my $cmd = "mysql -e \"".$query."\" -quick -h$host -P$port -u$user ".(($pass)? "-p$pass" : "")." $dbname > ${fset}_data.txt";
    #print $cmd."\n";
    system($cmd);
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
    print FO "$safename <- read.table(\"${fset}_data.txt\",header=TRUE)\n";
    #the if is to avoid errors printing empty sets...
    print FO "if (length(".$safename."\$length) > 0) barplot(unlist(lapply(split(".$safename.",".$safename."\$region),function(x) length(x\$length))),main=\"".$fset."\",xlab=\"region\",ylab=\"number of peaks\")\n";
    print FO "if (length(".$safename."\$length) > 0) boxplot(lapply(split(".$safename.",".$safename."\$region),function(x) x\$length),main=\"".$fset."\",xlab=\"region\",ylab=\"Peak length\")\n";
  }

  if($compare){
    #Maybe do some more analysis here?
    my @safe; foreach my $fset (@fset_names){ push(@safe, $safenames{$fset}); }
    print FO "barplot(c(length(".join("\$length),length(",@safe)."\$length)),main=\"Comparison of Annotated Features\",xlab=\"Feature set\",ylab=\"number of peaks\",col=rainbow(".scalar(@safe)."))\n";
    print FO "legend(\"topright\",legend=c(\"".join("\",\"",@fset_names)."\"),fill=rainbow(".scalar(@fset_names)."))\n";
    print FO "boxplot(".join("\$length,",@safe)."\$length,main=\"Comparison of Annotated Features\",xlab=\"Feature set\",ylab=\"Peak length\",col=rainbow(".scalar(@safe)."))\n";
    print FO "legend(\"topright\",legend=c(\"".join("\",\"",@fset_names)."\"),fill=rainbow(".scalar(@fset_names)."))\n";
  }


  close FO;
  system "R CMD BATCH --slave peaks_report.R";
}

__END__




