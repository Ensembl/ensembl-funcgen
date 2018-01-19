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

run_moods.pl -- Runs the Moods find_pssm_dna pwm mapping algorithm 

=head1 DESCRIPTION




=head1 OPTIONS

=over

Mandatory:

=item B<-mapper>

Mapper name or full path

=item B<-target>

Target fasta file

=item B<-queries>

List of query files

=item B<-output>

Output file name

Optional:

=item B<-format>

Output format required. This is normally provided via a post processing method in MotifTools.

=item B<-help|man>

Shows this message

=back

=cut

use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use Bio::EnsEMBL::Funcgen::Utils::EFGUtils        qw( open_file );
use Bio::EnsEMBL::Funcgen::Sequencing::MotifTools qw( run_moods );
#use Bio::EnsEMBL::Funcgen::Utils::DBAdaptorHelper qw( get_DB_options_config
#                                                      create_DBConnection_from_options );

#In case we need to pass mapper specific params
#Getopt::Long::Configure("pass_through"); #Allows unknown options to pass on to @ARGV


my ($mapper, $job_array, @target_fastas, @query_files, 
    $p_thresh, $out_prefix, $format, $debug_level);
#my $db_opts = get_DB_options_config(['db']); #default db option names

#Keep this simple and moods specific for now
#rather than following the generic implementation of SeqTool::run_aligner


sub main(){
#  my @tmp_args = @ARGV;
  print "run_moods.pl @ARGV\n";

  GetOptions (
    #Mandatory
    #'mapper=s'     => \$mapper,
    'targets=s{,}' => \@target_fastas,
    'queries=s{,}' => \@query_files,
    'out_prefix=s' => \$out_prefix,
    #Optional
    'job_array'   => \$job_array,
    'format=s'    => \$format,
    'p_thresh=s'  => \$p_thresh,
    #%$db_opts,#dbname, host, user, pass, port
    'debug=i'   => \$debug_level,
    'help'      => sub { pod2usage(-exitval => 0); }, 
    'man|m'     => sub { pod2usage(-exitval => 0, -verbose => 2); },
    )  or pod2usage(-exitval => 1); 
  #can't do this if we enable pass_through 
  #as we expect the unexpected! 
  #  my @extra_params = @ARGV;

  
  #Handle job array

  my $target_fasta;

  if((scalar(@target_fastas) > 1) && ! $job_array){
    die('');
  }
  elsif($job_array){

    if(! defined $ENV{LSB_JOBINDEX}){
      die('-job_array has been specified but $ENV{LSB_JOBINDEX}'.
        ' is not available. Was this bsub\'d as an array?');
    }

    $target_fasta = $target_fastas[$ENV{LSB_JOBINDEX} -1];
  }
  else{ #scalar(@target_fastas) < 1
    $target_fasta = $target_fastas[0];
  }

  #Define output

  #We may have passed many queries
  #so we do need to pass an out_prefix, which we will 
  #append with the target base name sufficed with out or tab
  #will have to recreate this logic in caller anyway
  

  #All mandatory param check done in run_mapper
  #This somewhat hides the cmdline interface of the mapper
  #but allows for calling of different mappers with generic params
  #and optional pass_through params

  #Not setting { 'RaiseError' => 1 } here

  #my $jdbc = (defined $format) ?  
  # create_DBConnection_from_options($db_opts, 'db') : undef;

  #How are we going to handle output prefix here?
  #With another array?
  #No these are really intermediate files, so the name is not that important
  #Do it based on the output 


  #Handle out file recovery here on in run_moods?

  run_moods(-target     => $target_fasta,
            -queries    => \@query_files,
            -out_prefix => $out_prefix,
            -format     => $format,
            -p_thresh   => $p_thresh, #Default in run_moods is 0.001
            -debug      => $debug_level,
            #-jaspar_dbc => $jdbc,
            #-mapper_params => @extra_params
           );
}

#This avoids executing run_aligner when including using require 
#i.e. as a library for general usage or when unit testing
#Slightly redundant as there is nothing much to test in here really?
main() unless(caller);
1;