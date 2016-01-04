#!/usr/bin/env perl

=head1 LICENSE

Copyright [1999-2016] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute

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

seed_hive.pl

=head1 SYNOPSIS

  seed_hive.pl -url <String> -hive_script_dir <String> -input_id <String> \
   -logic_name <IdentifyAlignInputSubsets | IdentifyReplicateResultSets | IdentifyMergedResultSets> \
   [ -no_write -help -man ]

=head1 PARAMETERS

  Mandatory:
    -logic_name       <String>      IdentifySetInputs logic name e.g. IdentifyAlignInputSubsets, 
                                    IdentifyReplicateResultSets or IdentifyMergedResultSets
    -input_id         <String>      Input id for the given logic name. This is a hash containing Set 
                                    filter key value pairs and other options i.e.
                                    '{#define lists as comma separated string
                                      feature_types       => "ftype_name1,ftype_name2",   
                                      #or as listref           
                                      cell_types          => ["celltype_name1", "celltype_name2"], 
                                      experimental_groups => "egroup_1",
                                      analyses            => ["set_analysis_1", "set_analysis_2"],
                                      #Review these states, uses OR logic here as opposed to AND for above
                                      states              => ["RERUN", "TO_RUN"],            

                                      #Or direct specificy sets names or dbIDs directly (cannot use comma separated string) e.g.
                                      set_names           => ["set_name_1", "set_name_2"],
                                      #OR
                                      set_ids             => [1,2],
                                      #OR 
                                      experiments         => ["CTYPE_FTYPE_EXPGROUP(_ARCHIVEID)", "CTYPE_%_EXPGROUP(_ARCHIVEID)"],

                                      #OPTIONAL DATAFLOW PARAMS TO THE NEXT LOGICAL ANALYSIS 
                                      #e.g. DefineOutputSet
                                      dataflow_params     => {recover=>1} }'
    WARNING: ONLY USE QUOTING STYLE AS ABOVE                                    
    WARNING: If calling from a bash script then omit all spaces to avoid word splitting                                  
    -url              <String>      Hive DB URL i.e. mysql://DB_USER:DB_PASS@DB_HOST:DB_PORT/DB_NAME
                                    (Normally defined by environment)
    -hive_script_dir  <String>      Hive scripts directory (Normally defined by environment)
                                      
  Optional:
    -no_write         Seeds pipeline and runs Identify analysis in -no_write mode to
                      list Sets which would be identified by given input_id
    -debug    [1-3]   Sets the debug output level for -no_write mode
    -help             Prints a helpful message
    -man              Prints the man page

=head1 DESCRIPTION

This is a wrapper script for the main hive seed_pipeline.pl script. It validates the inputs required
by the Ensembl regulation sequencing pipeline for the various IdentifySetInputs analyses, as well as 
handling a -no_write mode which can also run the given analysis to show which sets would be identified
from the given input_id. 

=cut

use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;

use Bio::EnsEMBL::Hive::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Funcgen::Utils::EFGUtils qw( run_system_cmd run_backtick_cmd );

my @valid_lnames = qw(IdentifyAlignInputSubsets IdentifyReplicateResultSets IdentifyMergedResultSets);

&main();

#todo
#1 Allow to run with no arguments, which will identify the seedable analyses???
#  Arguably done need this wrapper for that use, and this wrapper also lists the analyses
#  in the pod

sub main{
  my @tmp_args = @ARGV;
  my ($hive_url, $hive_script_dir, $lname, $no_write, $input_id);               
  my $debug = '';
  
  GetOptions 
   (#Mandatory
    'logic_name=s'      => \$lname,
    'url=s'             => \$hive_url,
    'hive_script_dir=s' => \$hive_script_dir,
    'input_id=s'        => \$input_id,
  
    #Optional
    'no_write'         => \$no_write,
    'debug=i'          => \$debug,
    'help'             => sub { pod2usage(-exitval => 0); }, 
    #removed ~? frm here as we don't want to exit with 0 for ?
    'man|m'            => sub { pod2usage(-exitval => 0, -verbose => 2); },
   ) or pod2usage(-exitval => 1, -message => "Specified parameters are:\t@tmp_args"); 
             
 
  
  ### VALIDATE PARAMETERS ###
  
  if(! ((defined $lname) && (defined $hive_url) && 
        (defined $hive_script_dir) && (defined $input_id)) ){
    pod2usage(-exitval => 1, 
              -message => "Mandatory parameters not met:\n\t-logic_name\t\t".($lname || ''). 
                            "\n\t-input_id\t\t".($input_id || ''). 
                            "\n\t-url\t\t\t".($hive_url || ''). 
                            "\n\t-hive_script_dir\t".($hive_script_dir|| '')."\n");      
  }
  
  if( ! grep(/^${lname}$/, @valid_lnames) ){
    pod2usage(-exitval => 1, 
              -message => $lname.' is not a valid logic_name');
  }
  
  
  
  if($no_write){
    $input_id =~ s/[,]*}$/,no_write=>1}/;
  }
   
  $debug = '-debug '.$debug if $debug; 
   
  warn "input_id is $input_id" if $debug; 
   
  ### SEED THE PIPELINE ### 
  
  #Need to preprocess $input_id to escape single quotes
  
  my $cmd = "perl $hive_script_dir/seed_pipeline.pl -url $hive_url ".
              "-logic_name $lname -input_id '$input_id' 2>&1";
  #We don't get an exit status here if the job_id already exists even though it is not seeded. 
  #Catch STDERR here has the transformed input_id is printed to STDERR
              
  warn "Running $cmd" if $debug;
  my $job_info = run_backtick_cmd($cmd);  # This get's polluted with STDERR sometimes
  # TODO Use IPC run to get just STDOUT.
  my $exit_status = $?;    
  warn "exit status is $?\njob_info is $job_info\n" if $debug;
  
  
  ### GET THE JOB ID ###                                    
  my $job_id;
  my $id_stored = 0;
  
  if( $job_info =~ /Could not create job/){
    # This will have been reformated as it appears in the DB
    # Can we replace the code which actually reformat the input_id?
    (my $stored_input_id = $job_info) =~ s/^.*Could not create job '//s;
    $stored_input_id =~ s/' \(it may have.*(?:\n.*)*$//;
    #optional multi line subsitution as it may have been suffixed with another 
 
    #Now try and fetch job_id from the job or analysis_data table;  
    my $hive_db = Bio::EnsEMBL::Hive::DBSQL::DBAdaptor->new(-url => $hive_url);
    
    my $sql_helper = $hive_db->dbc->sql_helper;
    my $sql = "select job_id from job where input_id='$stored_input_id'";
    #warn "\n\nSQL $sql";
    $job_id = $sql_helper->execute_single_result($sql, undef, undef, undef, 1);#no throw flag
    
    if(! defined $job_id){
      #todo, this need to join to the analysis table also!
      
      
      $sql = 'select job_id from job j, analysis_data ad where ad.data='.
        "'$stored_input_id' and j.input_id =concat('_extended_data_id ', ad.analysis_data_id)"; 
      #warn "\n\nSQL $sql";  
      $job_id = $sql_helper->execute_single_result($sql, undef, undef, undef, 1);#no throw flag
    }
    
    if(! defined $job_id){
      die("Failed to identify existing job_id for input_id:\n$input_id\nUsing SQL:\t$sql");  
    } 
    
    $id_stored = 1;    
  }
  elsif($exit_status){
    die("Failed to execute:\n$cmd\n$job_info");
  }
  else{
    #Handles multi-line output
    ($job_id = $job_info) =~ s/(?:.*\n)*.*Job ([0-9]+) \[ .*/$1/;
    
    if($job_id !~ /^[0-9]+$/){
      my $err = "Seeded job, but failed to parse job id from script output:\n$job_info";
      $err .= "\nSkipping -no_write run for $lname" if $no_write;
      die($err);
    }      
  }
  
  ### RUN IdentifySetInputs -no_write 
  
  if( $no_write ){
    $cmd = "perl $hive_script_dir/runWorker.pl -url $hive_url -no_write -job_id $job_id -force 1 $debug";
    #warn "Running $cmd";
    run_system_cmd($cmd);
  }
  else{
    #would be nice to get some output here, by actually running IdetifyInputSubsets
    
    if($id_stored){
      print "Job ($job_id) already previously seeded for input_id:\n$input_id\n";
    }
    else{
      print "Job ($job_id) successfully seeded:\n$input_id\n";  
    }
  }
  
  return;
}# end of main


### TA DAA! ###

1;
