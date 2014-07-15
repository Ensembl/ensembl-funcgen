#!/usr/bin/env perl

=head1 LICENSE

Copyright [1999-2014] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute

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

run_pwm_mapping_pipeline.pl -- gets all relevant data and runs the pwm mapping

=head1 DESCRIPTION

fetches all TFs for which there is data in the given funcgen database and runs the mapping pipeline with the matrices corresponding to those TFs

Requires 
/usr/local/ensembl/bin to be in $PATH, as well as 
$EFG_SRC/scripts and 
ensembl-personal/dkeefe/perl/scripts
TODO: Need to eliminate most of these dependencies...

=head1 OPTIONS

=over

=item B<help>

Shows this

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

=item B<-dnadb_pass>

Password for the specific core database user 

=item B<-dnadb_port>

Port of the host where the specific core database to use is

=item B<-dnadb_name>

Name of the specific core database to use

=item B<-workdir>

Folder where the data is found

=item B<-outputdir>

folder where the output of the pipeline will go

=item B<-species>

Species e.g. mus_musculus

=item B<-assembly>

Assembly e.g. GRCh37_58_37c

=item B<-schema>

Schema Build e.g. 61_37n

=item B<-feature_type_list>

Space separated list of feature types e.g. Max cMyb .
If not specified all FeatureTypes with class 'Transcription Factor' will be used

=back

=cut

use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Funcgen::Utils::EFGUtils qw( open_file );
use Bio::EnsEMBL::Utils::Exception qw(throw warning stack_trace_dump);

my ($host, $port, $user, $pass, $dbname);
my ($dnadb_host, $dnadb_port, $dnadb_user, $dnadb_pass, $dnadb_name);
my ($help, $workdir, $outputdir, $species, $assembly, $schema);


#get command line options

my (@fts, @fset_names);

print "run_pwm_mapping_pipeline.pl @ARGV\n";

#Maybe pass these to our bin folder?
# /usr/local/ensembl/bin/fastaexplode & fastaclean
GetOptions (
	    'dnadb_host=s'        => \$dnadb_host,
	    'dnadb_user=s'        => \$dnadb_user,
	    'dnadb_port=i'        => \$dnadb_port,
	    'dnadb_pass=s'        => \$dnadb_pass,
	    'dnadb_name=s'        => \$dnadb_name,
	    'host=s'           => \$host,
	    'user=s'           => \$user,
	    'port=i'           => \$port,
	    'pass=s'           => \$pass,
	    'dbname=s'           => \$dbname,
	    'workdir=s'          => \$workdir,
	    'outputdir=s'        => \$outputdir,
	    'species=s'          => \$species,
	    'assembly=s'         => \$assembly,
	    'schema=s'           => \$schema,
	    "feature_type_list=s{,}"  => \@fts,
	    "help|h"             => \$help,
	   )  or pod2usage( -exitval => 1 ); #Catch unknown opts

pod2usage(1) if ($help);

# Should be failing a little nicer now... 
if(!$host || !$port || !$user || !$dbname ) {  print "Missing connection parameters for funcgen db\n"; pod2usage(1); }
if(!$workdir || !$outputdir) {  print "Missing working folder(s)\n"; pod2usage(0); }
if(!$species || !$assembly || !$schema) {  print "Need species, assembly and schema information\n"; pod2usage(1); }


#Check database connections
my $coredba;
if($dnadb_name){
  
  my $dnadb_pass = {-pass => $dnadb_pass};
	  
  my $coredba = Bio::EnsEMBL::DBSQL::DBAdaptor->new
    (
     -host => $dnadb_host,
     -port => $dnadb_port,
     -user => $dnadb_user,
     -dbname => $dnadb_name,
     -species => $species,
     -group   => 'core',
     %$dnadb_pass
    );
}

my ($efgdba, $apass);

if($dbname){

  if ($pass){
    $apass = {-pass => $pass};
    $pass = "-p $pass";
  }


  $efgdba = Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor->new
    (
     -host    => $host,
     -port    => $port,
     -user    => $user,
     -dbname  => $dbname,
     -species => $species,
     -dnadb   => $coredba, #Assumes that new will accept undef as parameter for this...
     %$apass
    );
}

#Test connections
$efgdba->dbc->db_handle;
$efgdba->dnadb->dbc->db_handle;

my $fta = $efgdba->get_FeatureTypeAdaptor();
my $fsa = $efgdba->get_FeatureSetAdaptor();
my $bma = $efgdba->get_BindingMatrixAdaptor();
my $dbea = $efgdba->get_DBEntryAdaptor();


#Separate this into a separate step/script

=pod

if(! -d "${outputdir}/matrices"){
  system("mkdir ${outputdir}/matrices") && die "Error creating matrices folder";
}
open(FO,">".$outputdir."/matrix_list");

my @tfs;
if(scalar(@fts) > 0){
  
  foreach my $ft (@fts){
    my $tf = $fta->fetch_by_name($ft); 
  
    if(!$tf){
      throw "Could not find Transcription Factor $ft";
    } else {
      push(@tfs, $tf);
    } 
  }
}
else {
  #If none is specified, get all...
  #warn "Using all transcription factors";
  @tfs = @{$fta->fetch_all_by_class('Transcription Factor')};
  
  #Change this to use th bm ftypes?
  #but it appears that the matrix list uses the fset ftypes
  
}

# Just dump out the individual pfm files here from the new merged file
# as thats what the pipeline expects


my $pfms_file = "${workdir}/binding_matrices/Jaspar_5.0/JASPAR_CORE/pfm/nonredundant/pfm_all.txt";

my $pfms_fh = open_file($pfms_file);
my ($line, $pfm_fh, $record, $mname);

while(($line = $pfms_fh->getline) && defined $line){
    
  if($line =~ /^>([^ ]+)/){
    
    
    if(defined $pfm_fh){
      print $pfm_fh $record;
      $pfm_fh->close;
    } 
    
    $mname  = $1;
    chomp $mname;
    $pfm_fh = open_file("${outputdir}/matrices/${mname}.pfm", '>');
    $record = '';#$line;
  }
  else{
    $record .= $line;    
  }
}

print $pfm_fh $record;
$pfm_fh->close;


foreach my $tf (sort { $a->name cmp $b->name} @tfs){

  # if there is displayable data for it, get its matrix(ces) and print them.
  #my @fs = @{$fsa->fetch_all_by_FeatureType($tf,'displayable')};
  my @fs = @{$fsa->fetch_all_by_FeatureType($tf)};

  if(scalar(@fs)>0){
    #if tf is NOT a complex ft, but s involved in complexes, get their matrices
    # to see if the tf is a complex, just try getting an associated gene.
    my @matrices = @{$bma->fetch_all_by_FeatureType($tf)};
    
    
    #? There may be many other DBEtries for a FeatureType
    #Just get the associated ftypes anyway
    
    #In fact why are we even bothering with the TF
    #why not just access the bms directly?
    
    #keep this as is for now to get it running
   
    
    #if(scalar(@{$dbea->fetch_all_by_FeatureType($tf)})>0){
      
      
    foreach my $atf (@{$fta->fetch_all_by_association($tf)}){
    	push @matrices, @{$bma->fetch_all_by_FeatureType($atf)};
    }
    #}
    
    my @names;
    
    foreach my $matrix (@matrices){
      push @names, $matrix->name;
      #system("cp ${workdir}/binding_matrices/Jaspar_5.0/all_data/FlatFileDir/".$matrix->name.".pfm ${outputdir}/matrices") && die "could not copy matrix ".$matrix->name."";
      
    }
    
    if(scalar(@names)>0){
      my $namestr=join(";",@names);
      print FO $tf->name."\t".$namestr."\n";
    
    }
  }
}
close FO;

=cut

#matrix_list no longer exists as we use the DB
#where is this being used??!!
#system("cp ${workdir}/binding_matrices/Jaspar/matrix_list.txt ${outputdir}/matrices") && die "could not find fasta file";

#print "processing fasta\n";
#WHAT? Why are we gunzipping and copying here, why is this not used in situ?
#system("gunzip -dc ${workdir}/fasta/${species}/${species}_male_${assembly}_unmasked.fasta.gz > ${outputdir}/fasta.fas") && die "could not unzip fasta file";

my $fasta_file = "${workdir}/fasta/${species}/${species}_male_${assembly}_unmasked.fasta";
$assembly =~ s/_.*$//;
print $assembly."\n";

my $map_file = "${outputdir}/all_mappings.tab";
#This is currently hardcoded in pwm_filter_map!?

#Could do with a recover or force mode here

if(! -f $map_file){


  my $cmd = "pwm_genome_map.pl -g ${outputdir}/fasta.fas -a ${assembly} -o $map_file ".
    "-p 0.001 -w ${outputdir}/tmp_results/ ${outputdir}/matrices/*.pfm";
  print $cmd."\n";
  system($cmd) && die "error running pwm_genome_map"; 
}

#stashed
#  #There is a danger that the workdir will be removed by this script
#  #so append tmp_results in the script not here, and make mandatory
#
#  my $cmd = "pwm_genome_map.pl -g $fasta_file -a ${assembly} -o $map_file ".
#    "-p 0.001 -w $outputdir ${outputdir}/matrices/*.pfm";
#  print $cmd."\n";
#  #warn "skipping pwm_genome_map\n";
#  system($cmd) && die "error running pwm_genome_map"; 
#}

#print "Creating fasta header file...\n";
#system("grep '>' $fasta_file > ${outputdir}/fasta.id_lines") && die "error processing fasta file";




#The output thresholds need to be stored in the DB for future reference, and used as a factor in defining the permissive threshold above
system("pwm_filter_mappings.pl -t $map_file -i ${outputdir}/matrix_list -e ".$dbname." -H ".$host." -u ".$user." -P ".$port." -o ${outputdir}/thresholds -g ${outputdir}/fasta.id_lines -s ${schema}") && die "could not run pwm_filter_mappings";
#stashed
#system("pwm_filter_mappings.pl -t $map_file -i ${outputdir}/matrix_list -e ".$dbname." -H ".$host." -u ".$user." -P ".$port." $pass -o ${outputdir}/thresholds -g ${outputdir}/fasta.id_lines -s ${schema}") && die "could not run pwm_filter_mappings";


system("mkdir -p ${outputdir}/filtered") && die "could not create filtered folder";

system("pwm_thresh_filter.pl -i ${outputdir}/thresholds -o ${outputdir}/filtered") && die "could not run pwm_thresh_filter";

