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

Assembly e.g. GRCh38

=item B<-feature_type_list>

Space separated list of feature types e.g. Max cMyb .
If not specified all FeatureTypes with class 'Transcription Factor' will be used

=back

=cut



#TODO 
# 1 Pull this apart and put it into hive, some of the following may disappear when  this is done
# 2 Allow skip matrix extract in place of a matrix dir.
# 3 Make assembly optional, default to dnadb default


use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Funcgen::Sequencing::MotifTools qw( read_matrix_file
                                                      reverse_complement_matrix 
                                                      write_matrix_file
                                                      get_revcomp_file_path          );
use Bio::EnsEMBL::Funcgen::Utils::EFGUtils        qw( open_file run_system_cmd run_backtick_cmd );
use Bio::EnsEMBL::Utils::Exception                qw( throw warning stack_trace_dump );

my ($host, $port, $user, $pass, $dbname);
my ($dnadb_host, $dnadb_port, $dnadb_user, $dnadb_pass, $dnadb_name);
my ($help, $workdir, $outputdir, $species, $assembly);
my (@fts, @fset_names);

print "run_pwm_mapping_pipeline.pl @ARGV\n";

GetOptions (
 'dnadb_host=s' => \$dnadb_host,
 'dnadb_user=s' => \$dnadb_user,
 'dnadb_port=i' => \$dnadb_port,
 'dnadb_pass=s' => \$dnadb_pass,
 'dnadb_name=s' => \$dnadb_name,
 'host=s'     => \$host,
 'user=s'     => \$user,
 'port=i'     => \$port,
 'pass=s'     => \$pass,
 'dbname=s'     => \$dbname,
 'workdir=s'    => \$workdir,
 'outputdir=s'  => \$outputdir,
 'species=s'    => \$species,
 'assembly=s'   => \$assembly,
 "feature_type_list=s{,}"  => \@fts,
 "help|h"             => \$help,
) or pod2usage( -exitval => 1 ); #Catch unknown opts
  
pod2usage(1) if ($help);

# Should be failing a little nicer now... 
if(!$host || !$port || !$user || !$dbname ) {  print "Missing connection parameters for funcgen db\n"; pod2usage(1); }
if(!$workdir || !$outputdir) {  print "Missing working folder(s)\n"; pod2usage(0); }
if(!$species || !$assembly ) {  print "Need species and assembly information\n"; pod2usage(1); }


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


### Identify FeatureSets with associated BindingMatrices ###
#Separate this into a separate step/script
my $fta  = $efgdba->get_FeatureTypeAdaptor();
my $fsa  = $efgdba->get_FeatureSetAdaptor();
my $bma  = $efgdba->get_BindingMatrixAdaptor();
#my $dbea = $efgdba->get_DBEntryAdaptor();
my @tfs;

#Could remove this bit is we implement FeatureSetAdaptor::_constrain_feature_type_class
#Although this woudl not allow validation os specified FeatureTypes

#@fts = ('RXRA');


if(scalar(@fts) > 0){
  
  foreach my $ft (@fts){
    my $tf = $fta->fetch_by_name($ft); 
  
    if(! $tf){
      throw "Could not find Transcription Factor $ft";
    } 
    else {
      push @tfs, $tf;
    } 
  }
}
else {  #If none is specified, get all...
  @tfs = @{$fta->fetch_all_by_class('Transcription Factor')};
  #Change this to use th bm ftypes?
  #but it appears that the matrix list uses the fset ftypes
}


my %matrix_tf;
#If we ever support more than one source (i.e. something other than Jaspar)
#Then we need to restrict this list based on the source DB


foreach my $tf (sort { $a->name cmp $b->name} @tfs){
  # if there is displayable data for it, get matrices and cache print them.
  #my @fs = @{$fsa->fetch_all_by_FeatureType($tf,'displayable')};


  if(@{$fsa->fetch_all_by_FeatureType($tf)}){
    #Get direct TF binding_matrix associations     

    map { $matrix_tf{$_->name} = $tf->name } @{$bma->fetch_all_by_FeatureType($tf)};

    #And indirect associations through TF complexes
    foreach my $atf (@{$fta->fetch_all_by_association($tf)}){
      map { $matrix_tf{$_->name} = $tf->name } @{$bma->fetch_all_by_FeatureType($atf)};
    }
  }
}

if(! -d "${outputdir}/matrices"){
  system("mkdir ${outputdir}/matrices") && die "Error creating matrices folder";
}



### Create Individual pfm files ###

# Moods expects individual pfm files, but the new Jaspar release only 
# has a single mutli-record file.

# Currently also need collections which have ID prefix PB (PBM) PL (PBM_LHL) and PH (PBM_HOMEO)
# These matrices need to be grabbed form the previous release (JASPAR2010),a nd some even from the
# release before that!!
# They are not present in the current release as they have not been updated in Jaspar 5(???), even though 
# they are they are present in DB, the ARCHIVE sql directories and on the website(????!)
# This needs to be done manually before running this script
# TODO Also need to update the association script to use these collections.

# Have to use all these files, as there are omissions from each, even though they are in the 
# current live web site.
# These are in order of preference

my @source_pfm_files = (
  "${workdir}/binding_matrices/Jaspar_5.0/JASPAR_CORE/pfm/nonredundant/pfm_all.txt",
  "${workdir}/binding_matrices/Jaspar_5.0/JASPAR_CORE/pfm/redundant/pfm_all.txt",
  "${workdir}/binding_matrices/JASPAR2010/all_data/matrix_only/matrix.txt",
  "${workdir}/binding_matrices/Jaspar_pre_2010/all_data/matrix_only/matrix_only.txt",
  );

#Why are we now not getting these?
#PL0007.1(Max) is associated to a FeatureSet but is not present in the source pfm file
#PH0144.1(POU2F2) is associated to a FeatureSet but is not present in the source pfm file

#PH0144.1 is in 2010 and pre 2010, why are we not parsing it?


#Could even fork prior to the matrix file read, and submit batch jobs
# based on the matrix/tf name pair. This would however mean multiple parallel reads
# from the same matrix file.
# Also need to consider the ability to run a pfm without the dependancy on a Feature set.
# Basically we need to keep the Featture set identification, the pwm mapping and the theshold 
#generation/filtering steps separate.

my ($pfm_file, %pfm_files);
my @matrix_ids    = keys(%matrix_tf);
my $num_matrices = scalar(@matrix_ids);
#We could check for the pfm and recvomp pfm file in the loop above, to prevent
#reading and recreateing the pfm files everytime.

print "Extracting matrix files...\n";

PFM_FILE: foreach my $pfms_file(@source_pfm_files){
  print "Reading matrices from:\t$pfms_file\n";
  #restrict read matrices here, not below
  my $matrix_info = read_matrix_file($pfms_file, \@matrix_ids);

  foreach my $matrix_id(keys %$matrix_info){
    #Ensure we haven't see this matrix previously i.e. in the newer file

    #if(! exists $matrix_tf{$matrix_id}){
    #  die("$matrix_id is not in original filter list");
    #}

    if(! exists $pfm_files{$matrix_id}){
      #Original matrix
      #warn "Processing $matrix_id";
      $pfm_file = "${outputdir}/matrices/${matrix_id}.pfm";
      write_matrix_file([$matrix_info->{$matrix_id}], 
                        $pfm_file);
      $pfm_files{$matrix_id} = $pfm_file;

      #Revcomp matrix!
      $matrix_info->{$matrix_id}{matrix} = 
       reverse_complement_matrix($matrix_info->{$matrix_id}{matrix});
      write_matrix_file([$matrix_info->{$matrix_id}], 
                        get_revcomp_file_path($pfm_file));


      $num_matrices--;
      last PFM_FILE if ! $num_matrices;   
    }
  }
}


#This needs to be passed as an arg!
#Keep all the input file path architecture in the caller

my $fasta_file = "${workdir}/fasta/${species}/${species}_male_${assembly}_unmasked.fasta";


$assembly =~ s/_.*$//;
print "Using assembly:\t$assembly\n";


#Could do with a recover or force mode here


#TO DO
# Move all the functionality from these scripts to MotifTools
# So we never have to system call a separate script as this
# loses error output (at least without IPC::Cmd/Run which isn't in 5.10)


### Validate we have found all matrices registered in the DB
my $absent = 0;

foreach my $mname(keys(%matrix_tf)){

  if(! exists $pfm_files{$mname}){
    warn $mname.'('.$matrix_tf{$mname}.") is associated to a FeatureSet but is not present in the source pfm file\n";
    $absent = 1;
  }
}

if($absent){
  die("Identified some absent matrices in the input pfm files:\n\t".
    join("\n\t", @source_pfm_files));
}


my $cmd = "pwm_genome_map.pl -g $fasta_file -a ${assembly} ".
 "-p 0.001 -o ${outputdir} ".join(' ', values(%pfm_files));
#This bsubs a job array of run_moods.pl jobs for each matrix (inc rev comp) vs chr fasta files
#Then waits for them to finish before merging
run_system_cmd($cmd); 

my $fasta_headers_file = "${outputdir}/fasta.headers.txt";

#This is currently broken as the Y chr header was wrong

if(! -f $fasta_headers_file){
  print "Creating fasta header file...\n";
  run_system_cmd("grep '>' $fasta_file > ${fasta_headers_file}.tmp");
  run_system_cmd("mv ${fasta_headers_file}.tmp ${fasta_headers_file}");
}

run_system_cmd("mkdir -p ${outputdir}/filtered");


#Let's just pass one map file and the TF from the matrix_list
#Map tab files are generated
#Or just add everything to pwm_genome_map?
#By moving code to MotifTools?

#The output thresholds need to be stored in the DB for future reference, and used as a factor in defining the permissive threshold above

my $thresholds_file = "${outputdir}/thresholds.txt";
my $lsf_out_dir     = "${outputdir}/filtered/lsf_out/";
mkdir($lsf_out_dir); #This likely already exists from pwm_genome_map.pl Actually, this may be removed if all is successful.


#Touch filter file here, so we don't get and failure when the pwm_file_mappings.pl script tries to append
system("touch $thresholds_file");

MATRIX: foreach my $matrix (keys %matrix_tf){
  #Defined/written in pwm_genome_map.pl process_and_align_matrix
  my $map_file        = "${outputdir}/${matrix}.unfiltered.bed";
  print "Preparing to post-process:\t$map_file\n";
  
  #if(-e $thresholds_file){
    my $grep_cmd    = "grep -E '^$matrix\[\[:space:\]\]' $thresholds_file";
    my $thresh_line = run_backtick_cmd($grep_cmd, 1);#grep returns exit code 1 if not present
  
    if($thresh_line){
      warn "Skipping $matrix as is already present in output:\n$thresh_line";
      next MATRIX;  
    }
  #}

  my $job_cmd = "pwm_filter_mappings.pl -m $matrix -t $map_file -T ".$matrix_tf{$matrix}." -e ".$dbname.
   " -H ".$host." -u ".$user." -P ".$port.' '.$pass.
   " -o ${outputdir} -g $fasta_headers_file -a $assembly";

  my $bjob_name = "FILTER_${matrix}";
  my $bsub_params = '-M 2000 -R"select[mem>2000] rusage[mem=2000]" '.
   "-o $lsf_out_dir/$bjob_name.\%J.out -e $lsf_out_dir/$bjob_name.\%J.err";
  my $cmd         = "bsub -J $bjob_name $bsub_params '$job_cmd'";
  print "$cmd\n";
  my $bsub_output = run_backtick_cmd($cmd);

  #print $bsub_output."\n";

  if($bsub_output =~ /Job \<([0-9]+)\> is submitted to /){
    $bsub_output = $1;
  }
  else{  
    #This is slightly brittle, we could do some retry here?
    throw("Failed to submit job:\n$cmd\n$bsub_output");
  }


#run_system_cmd("pwm_filter_mappings.pl -t $map_file -i ${outputdir}/matrix_list -e ".$dbname.
#  " -H ".$host." -u ".$user." -P ".$port.
#  " -o ${outputdir}/thresholds -g ${outputdir}/fasta.id_lines -s ${schema}");

#run_system_cmd("pwm_thresh_filter.pl -i ${outputdir}/thresholds -o ${outputdir}/filtered");
  
}