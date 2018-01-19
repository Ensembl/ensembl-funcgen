#!/usr/bin/env perl

=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2018] EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

  Questions may also be sent to the Ensembl help desk at
  <http://www.ensembl.org/Help/Contact>.

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

=item B<-species>

Species of the database (e.g. homo_sapiens)

=item B<-dbhost>

Host where the database is

=item B<-dbuser>

User of the database 

=item B<-dbpass>

Password for the database user 

=item B<-dbport>

Port of the host where the database is 

=item B<-dbname>

Name of the database 

=item B<-dnadb_host>

Host of the specific core database to use 

=item B<-dnadb_user>

User of the specific core database 

=item B<-dnadbdb_pass>

Password for the specific core database user 

=item B<-dnadb_port>

Port of the host where the specific core database to use is

=item B<-dnadb_name>

Name of the specific core database to use

=item B<-R>

When specified runs the R code to generate plots of the data.

=item B<-nodump>

When specified it will skip dumping the data 

=item B<-compare>

Will generate graphs comparing all datasets (if not specified, only set-specific graphs will be generated)

=item B<-regstats>

When specified it generate and plot statistics regarding the regulatory build present in the database
Requires the presence of a table rf_stats in the database

=item B<-feature_sets>

List of space separated FeatureSet names to be considered. Names with spaces must be quoted.
When not specified all features are considered

=item B<-feature_table>

Type of the feature: annotated or regulatory (default is annotated)

=item B<-no_outliers>

When specified, the outliers are not drawn in the plots (by default, outliers are drawn)

=item B<-all_seq_regions>

Processes all available seq_regions, by default only uses chromosomes.

=item B<-name>

Name for this report, default is 'peaks_report'

=item B<-inset_main>

Space given for the legend for chromosomes/regions. Default is 0.1 (works well for ~20 chromosomes)
Depends on number of regions.

=item B<-inset_compare>

Space given for the legend for comparing cell types. Default is 0.5 (works for ~10 cell-types)
Depends on number of cell types.

=item B<-outdir>

Base folder for output. Defaults to $WORK_DIR/output/$DB_NAME

=back


=head1 SEE ALSO

ensembl-funcgen/scripts/environments/peaks.env PeaksReport function

=cut


# To do
# 1 Move reg stats code in here and allow reg_stats mode only
# 2 Tighten up parameter checking
# 3 Change -feature_sets all sets behaviour?
# 4 Add outdir
# 5 Catch absent files (sugest removal of -no_dump?)
# 6 Handle output better lsf/R out? Use optional run name for file naming?
# 7 Use RMySQL to pull the data directly into R rather than dumping(use optional r workspace save to save data within R?)
# 8 Move to sub dir
# 9 DONE Restrict plots to main chromosomes (do NT contigs in a seaprate plot)..much quicker to view now
# 10 Show all chr names in plot axis
# 11 DONE Correct sql to prevent pruoduct with nr seq_region entries
# 12 Add -farm option or just do this in the run script as is only needed when dumping?
# 13 Order chromosome numerically rather than lexically
# 14 CHange axis labels to real numbers
# 15 Move comparison legend outside of plot as it obscures view of data
# 16 DONE Implement report name
# 17 DONE These are rather large, can we compress the data some how by using different plots styles? -no_outliers
# 18 FeatureSet names/colour do not appear on last graph axis/plot
# 19 Add support for comparison of sets between DBs. Dump from both DBs | cat | sort. Append dbname to fset names in plots?

use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor;
use Cwd;


my %feature_tables = (
					  annotated  => 1,#values?
					  regulatory => 1,
					  #external   => 1,
					 );


#Best not use these as default as they are easily forgotten
#my $host = $ENV{DB_HOST};
#my $port = $ENV{DB_PORT};
#my $user = $ENV{DB_READ_USER};
#my $pass = $ENV{DB_PASS};
#my $dbname = $ENV{DB_NAME};
#my $dnadbhost = $ENV{DNADB_HOST};
#my $dnadbport = $ENV{DNADB_PORT};
#my $dnadbuser = $ENV{DNADB_USER};
#my $dnadbname =  $ENV{DNADB_NAME};
#my $dnadbpass =  $ENV{DNADB_PASS};

my ($species, $help, $R, $nodump, $compare, $regstats, $all_seq_regions, $no_outliers, $name, $outdir);
my ($feature_table, $host, $port, $user, $pass, $dbname, $dnadbhost, $dnadbport, $dnadbuser, $dnadbname, $dnadbpass);
my ($inset_main, $inset_compare);

#Default values
$inset_main=0.1;
$inset_compare=0.5;
$user = 'ensro';
$name = 'peaks_report_'.$$;#Add PID to avoid overwriting previous reports
$outdir = $ENV{'EFG_DATA'}."/output/".$ENV{DB_NAME};

#get command line options

my (@fset_names_tmp);

print "peaks_report.pl @ARGV\n";

GetOptions (
			'species=s'          => \$species,
			'dnadb_host=s'       => \$dnadbhost,
			'dnadb_user=s'       => \$dnadbuser,
			'dnadb_port=i'       => \$dnadbport,
			'dnadb_pass=s'       => \$dnadbpass,
			'dnadb_name=s'       => \$dnadbname,
			'dbhost=s'           => \$host,
			'dbuser=s'           => \$user,
			'dbport=i'           => \$port,
			'dbpass=s'           => \$pass,
			'dbname=s'           => \$dbname,
			'outdir=s'           => \$outdir,
			"help|h"             => \$help,
			"R"                  => \$R,
			"nodump"             => \$nodump,
			"no_outliers"        => \$no_outliers,
			"compare"            => \$compare,
			"regstats"           => \$regstats,
			"all_seq_regions"    => \$all_seq_regions,
			"feature_sets=s{,}"  => \@fset_names_tmp,
	                "feature_table=s",   => \$feature_table,
	                "inset_main=s",      => \$inset_main,
	                "inset_compare=s",   => \$inset_compare,
			"name=s"             => \$name,
		   )  or pod2usage( -exitval => 1 ); #Catch unknown opts

pod2usage(1) if ($help);

#Reset to undef so we don't try with empty string
#$pass      ||= undef;
#$dnadbpass ||= undef;

# Sould be failing a little nicer now... 
if(!$feature_table) { print "Missing Type of Feature: annotated or regulatory (use -h for help)\n"; exit 0; }
if(!$host || !$port || !$user || !$dbname )  {  print "Missing connection parameters (use -h for help)\n"; exit 0; }
if(!$outdir )  {  print "\$EFG_DATA not defined and -outdir not specified\n"; exit 0; }

if(! $feature_tables{$feature_table}){
  die("You have specified an invalid -feature_table. Must be one of:\t".join("\t", (keys %feature_tables)));
}

if(! -d $outdir){
  die("Error: $outdir is not a valid output folder");
} else{
  
  print "Setting default output directory to:\t".$outdir;
  
  if(! -d $outdir){
    system("mkdir -p $outdir") == 0 or 
      die("Could not create output directory:\t".$outdir);
  }
}

#warn "dbpass is x${pass}x";
#warn "dnadbpass is x${dnadbpass}";

#Check database connections
my ($coredba, $efgdba);
if($dnadbname){
  
  my $coredba = Bio::EnsEMBL::DBSQL::DBAdaptor->new
    (
     -host => $dnadbhost,
     -port => $dnadbport,
     -user => $dnadbuser,
     -dbname => $dnadbname,
     -species => $species,
     -group   => 'core',
     -pass    => $dnadbpass,
    );
}

 
$efgdba = Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor->new
  (
   -host    => $host,
   -port    => $port,
   -user    => $user,
   -dbname  => $dbname,
   -species => $species, #Not strictly necessary
   -dnadb   => $coredba, 
   -pass    => $pass,
  );


#Test connections
$efgdba->dbc->db_handle;
$efgdba->dnadb->dbc->db_handle;

my $fsa = $efgdba->get_FeatureSetAdaptor();

my @fset_names;

if(scalar(@fset_names_tmp)==0){
  map { push @fset_names, $_->name } @{$fsa->fetch_all_by_type($feature_table)};
}
else{#Validate fset names
  
  foreach my $fsname(@fset_names_tmp){
	my $fset = $fsa->fetch_by_name($fsname);
	
	if($fset){
	  push @fset_names, $fset->name;
	} else {
	  warn("Could not fetch FeatureSet:\t$fsname ... Skipping");
	}
  }
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
  print "\n\n:: Dumping Datasets\n";
  
  #This was not accounting for nr sr_ids
  foreach my $sr_type(@sr_types){
    
    #Save to only one file... though it may be big...
    my $query ="SELECT fs.name as 'name', s.name as 'region', (f.seq_region_end - f.seq_region_start) as 'length' FROM ${feature_table}_feature f, (select distinct(seq_region_id), sr.name from seq_region sr, coord_system cs where sr.coord_system_id=cs.coord_system_id and cs.name".$sr_type_clauses{$sr_type}." and cs.is_current is TRUE) s, feature_set fs WHERE f.feature_set_id=fs.feature_set_id AND f.seq_region_id=s.seq_region_id AND fs.name IN ('".join("','",@fset_names)."');";
    
    my $cmd = "mysql -e \"".$query."\" -quick -h$host -P$port -u$user ".(($pass)? "-p$pass" : "")." $dbname >${outdir}/${name}.data.${sr_type}.txt";
    print $cmd."\n";
    system($cmd) == 0 || die('Failed to dump data');
  }
  
  if($regstats){
    my $query = "SELECT * from rf_stats";
    system("mysql -e \"".$query."\" -quick -h$host -P$port -u$user ".(($pass)? "-p$pass" : "")." $dbname > ${outdir}/regstats.txt");
  }

}

my $dir = getcwd;
#print R file with analysis
#Check which parameters one might want...
if (defined $R) {
  
  print "::Generating the R plots\n";
  open(FO,">${outdir}/${name}.R");
  print FO "require(graphics)\n";
  print FO "pdf(file=\"${outdir}/${name}.pdf\")\n";
  print FO "par(cex=0.7,cex.main=1)\n";
  
  if($regstats){ 
     print FO "require(gplots);\n";
     print FO "regstats <- read.table(\"${outdir}/regstats.txt\",row.names=1,header=TRUE);\n"; 
     print FO "textplot(regstats,halign=\"left\")\n";
  }
  
  foreach my $sr_type(@sr_types){
    
    #Load the data
    print FO "data_${sr_type} <- read.table(\"${outdir}/${name}.data.${sr_type}.txt\",header=TRUE,sep=\"\\t\")\n";
    
    #Give a little space for the legend... outside the graph... (test how much space... and the size of text in lengend)
    print FO "par(xpd=T, mar=par()\$mar+c(0,0,0,5))\n";
      
    #Print individual graphs
    print FO "for (subset in split(data_${sr_type},data_${sr_type}\$name)){\n";
    print FO "    subdata <- lapply(split(subset, subset\$region),function(x) x\$length)\n";
    print FO "    barplot(sapply(subdata, function(x) length(x)), main=subset\$name[1], xlab=\"region\",ylab=\"Number of Peaks\", col=rainbow(length(subdata)))\n";
    print FO "    legend('topright', inset=c(-".$inset_main.",0), legend=levels(subset\$region),fill=rainbow(length(subdata)), cex=0.8)\n";
    print FO "    boxplot(subdata,main=subset\$name[1],xlab='region',ylab='Peak length', col=rainbow(length(subdata))";
    if($no_outliers){ print FO ",outline=FALSE"; }
    print FO ")\n"; 
    print FO "    legend('topright',inset=c(-".$inset_main.",0),legend=levels(subset\$region), fill=rainbow(length(subdata)), cex=0.8)\n";    
    print FO "}\n";

  }

  if($compare){ 
    
    #For brevity only compare true chromosomes
    print FO "par(xpd=T, mar=par()\$mar+c(0,0,0,10))\n";

    #Global overview comparison
    print FO "data_region <-lapply(split(data_chromosome, data_chromosome\$name), function(x) x\$length)\n";
    print FO "barplot(sapply(data_region, function(x) length(x)),main='Number of Peaks per Set', xlab='Set',ylab='Number of Peaks', col=rainbow(length(data_region)), xaxt='n')\n";
    print FO "legend('topright', inset=c(-".$inset_compare.",0),legend=levels(data_chromosome\$name),fill=rainbow(length(data_region)), cex=0.8)\n";  
    #print FO "axis(1, labels=FALSE, at=1:length(data_region), tick=TRUE)\n";
    print FO "boxplot(data_region,main='Peak Length per Dataset',xlab=\"Set\",ylab=\"Peaks Length\", col=rainbow(length(data_region)), xaxt='n'";
    if($no_outliers){ print FO ",outline=FALSE"; }
    print FO ")\n";
    print FO "legend('topright',inset=c(-".$inset_compare.",0),legend=levels(data_chromosome\$name),fill=rainbow(length(data_region)), cex=0.8)\n";    
    print FO "axis(1, labels=FALSE, at=1:length(data_region), tick=TRUE)\n";

    #Print Comparative graphs by Region
    print FO "for (subset in split(data_chromosome,data_chromosome\$region)){\n";
    print FO "    subdata <- lapply(split(subset, subset\$name),function(x) x\$length)\n";
    print FO "    barplot(unlist(lapply(subdata, function(x) length(x))), main=subset\$region[1], xlab=\"region\",ylab=\"Number of Peaks\", col=rainbow(length(subdata)), xaxt='n')\n";
    print FO "    legend('topright',inset=c(-".$inset_compare.",0),legend=levels(subset\$name),fill=rainbow(length(subdata)), cex=0.8)\n";
    #print FO "    axis(1, labels=FALSE, at=1:length(subdata), tick=TRUE)\n";

    print FO "    boxplot(subdata,main=subset\$region[1],xlab='Set',ylab='Peak length', col=rainbow(length(subdata)), xaxt='n'";
    if($no_outliers){ print FO ",outline=FALSE"; }

    print FO ")\n"; 
    print FO "    legend('topright',inset=c(-".$inset_compare.",0),legend=levels(subset\$name),fill=rainbow(length(subdata)), cex=0.8)\n";    

    print FO "    axis(1, labels=FALSE, at=1:length(subdata), tick=TRUE)\n";
    print FO "}\n";
    
  }
    
   
  close FO;
  chdir($outdir);
  #This submits to the yesterday queue by default
  system "R CMD BATCH --slave ${name}.R";
}

chdir($dir);

__END__




