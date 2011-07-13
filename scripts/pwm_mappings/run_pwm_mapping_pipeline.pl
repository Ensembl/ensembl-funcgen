#!/usr/bin/env perl

=head1 LICENSE

  Copyright (c) 1999-2011 The European Bioinformatics Institute and
  Genome Research Limited.  All rights reserved.

  This software is distributed under a modified Apache license.
  For license details, please see

    http://www.ensembl.org/info/about/code_licence.html

=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <ensembl-dev@ebi.ac.uk>.

  Questions may also be sent to the Ensembl help desk at
  <helpdesk@ensembl.org>.

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

=item B<-dnadbhost>

Host of the specific core database to use 

=item B<-dnadbuser>

User of the specific core database 

=item B<-dnadbpass>

Password for the specific core database user 

=item B<-dnadbport>

Port of the host where the specific core database to use is

=item B<-dnadbname>

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
use Bio::EnsEMBL::Utils::Exception qw(throw warning stack_trace_dump);

my ($host, $port, $user, $pass, $dbname);
my ($dnadbhost, $dnadbport, $dnadbuser, $dnadbpass, $dnadbname);
my ($help, $workdir, $outputdir, $species, $assembly, $schema);


#get command line options

my (@fts, @fset_names);

print "run_pwm_mapping_pipeline.pl @ARGV\n";

#Maybe pass these to our bin folder?
# /usr/local/ensembl/bin/fastaexplode & fastaclean
GetOptions (
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
if(!$host || !$port || !$user || !$dbname ) {  print "Missing connection parameters for funcgen db\n"; pod2usage(0); }
if(!$workdir || !$outputdir) {  print "Missing working folder(s)\n"; pod2usage(0); }
if(!$species || !$assembly || !$schema) {  print "Need species and assembly information\n"; pod2usage(0); }


#Check database connections
my $coredba;
if($dnadbname){
  
  my $dnadbpass = {-pass => $dnadbpass};
	  
  my $coredba = Bio::EnsEMBL::DBSQL::DBAdaptor->new
    (
     -host => $dnadbhost,
     -port => $dnadbport,
     -user => $dnadbuser,
     -dbname => $dnadbname,
     -species => $species,
     -group   => 'core',
     %$dnadbpass
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

system("mkdir ${outputdir}/matrices") && die "Error creating matrices folder";
open(FO,">".$outputdir."/matrix_list");

my @tfs;
if(scalar(@fts)>0){
  foreach my $ft (@fts){
    my $tf = $fta->fetch_by_name($ft); 
    if(!$tf){
      throw "Could not find Transcription Factor $ft";
    } else {
      push(@tfs, $tf);
    } 
  } 
} else {
  #If none is specified, get all...
  warn "Using all transcription factors";
  @tfs = @{$fta->fetch_all_by_class('Transcription Factor')};
}

foreach my $tf (sort { $a->name cmp $b->name} @tfs){

  # if there is displayable data for it, get its matrix(ces) and print them.
  #my @fs = @{$fsa->fetch_all_by_FeatureType($tf,'displayable')};
  my @fs = @{$fsa->fetch_all_by_FeatureType($tf)};

  if(scalar(@fs)>0){
    #if tf is NOT a complex ft, but s involved in complexes, get their matrices
    # to see if the tf is a complex, just try getting an associated gene.
    my @matrices = @{$bma->fetch_all_by_FeatureType($tf)};
    if(scalar(@{$dbea->fetch_all_by_FeatureType($tf)})>0){
      foreach my $atf (@{$fta->fetch_all_by_association($tf)}){
    	push @matrices, @{$bma->fetch_all_by_FeatureType($atf)};
      }
    }
    my @names;
    foreach my $matrix (@matrices){
      push @names, $matrix->name;
      system("cp ${workdir}/binding_matrices/Jaspar/all_data/FlatFileDir/".$matrix->name.".pfm ${outputdir}/matrices") && die "could not copy matrix ".$matrix->name."";
      
    }
    if(scalar(@names)>0){
      my $namestr=join(";",@names);
      print FO $tf->name."\t".$namestr."\n";
    
    }
  }
}
close FO;

print "processing fasta\n";
system("cp ${workdir}/binding_matrices/Jaspar/matrix_list.txt ${outputdir}/matrices") && die "could not find fasta file";

system("gunzip -dc ${workdir}/fasta/${species}/${species}_male_${assembly}_unmasked.fasta.gz > ${outputdir}/fasta.fas") && die "could not unzip fasta file";

$assembly =~ s/_.*$//;
print $assembly."\n";
system("pwm_genome_map.pl -g ${outputdir}/fasta.fas -a ${assembly} -o ${outputdir}/all_mappings.tab -p 0.001 -w ${outputdir}/tmp_results/ ${outputdir}/matrices/*.pfm") && die "error running pwm_genome_map"; 

system("grep '>' ${outputdir}/fasta.fas > ${outputdir}/fasta.id_lines") && die "error processing fasta file";

system("pwm_filter_mappings.pl -i ${outputdir}/matrix_list -e ".$dbname." -H ".$host." -u ".$user." -P ".$port." -o ${outputdir}/thresholds -g ${outputdir}/fasta.id_lines -s ${schema}") && die "could not run pwm_filter_mappings";

system("mkdir -p ${outputdir}/filtered") && die "could not create filtered folder";

system("pwm_thresh_filter.pl -i ${outputdir}/thresholds -o ${outputdir}/filtered") && die "could not run pwm_thresh_filter";

