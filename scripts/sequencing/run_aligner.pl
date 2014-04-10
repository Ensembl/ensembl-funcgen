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

run_aligner.pl

=head1 SYNOPSIS

  run_aligner.pl  -aligner <String> -program_file <String> -query_file <String> \
    -target_file <String> [ -parameters '<String>' -output_dir <String> -debug <Int> -help -man ]
    ( + Aligner specific parameters ) 

=head1 PARAMETERS

  Mandatory:
    -aligner      <String>  Namespace of Bio::EnsEMBL::Funcgen::Hive::Aligner subclass e.g BWA
                            If the full namespace is not specified it will default to:
                              Bio::EnsEMBL::Funcgen::Hive::$aligner
    -program_file <String>  Path or name of program file e.g. bwa 
    -query_file   <String>  Query fastq file  e.g. my_reads.fastq 
    -target_file  <String>  Target fasta file e.g. homo_sapiens_female_GRCh37_unmasked.fasta
  
  Additional parameters may be required given the aligner specified e.g for BWA:
    -sam_ref_fai  <String>  Path to sam reference (target_file) fai index file.  
      
  Optional:
    -parameters <String>  A quoted string of command line parameters to pass to the aligner program
    -output_dir <String>  Output directory path
    -debug      <Int>     Debug level to pass on to Aligner.pm
    -help                 Prints the POD SYNOPSIS 
    -man                  Prints this full man page
 
=head1 DESCRIPTION

This is a simple wrapper script for the Bio::EnsEMBL::Funcgen::Hive::Aligner module. It will run
the aligner module specified. It has been written such that this script can be incoporated 
into an existing module/script using require e.g.

  require "$ENV{EFG_SRC}/scripts/sequencing/run_aligner.pl";
  my $output_file = run_aligner(@aligner_params);

This is most useful when wanting to run a quick analysis and the calling cantext does not have 
access to an Aligner Analysis object from a Funcgen DB. 

Alternatively, the ReadAlignment pipeline is a much more powerful way to run alignments and is 
also optimised to chunk the inputs to run parallel jobs.

=cut


use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;

use Bio::EnsEMBL::Funcgen::Utils::EFGUtils qw( run_system_cmd );
Getopt::Long::Configure("pass_through"); #Allows unknown options to pass on to @ARGV

&run_aligner();
#run_aligner rather than main for clarity when required into a module

sub run_aligner{
  #Assigning to @ARGV allows this script to be required and run from a module 
  #without having to fork via and external system call
  @ARGV = @_ if (! @ARGV) && @_;
  my $aln_pkg;

  GetOptions 
   (#Mandatory
    'aligner=s' => \$aln_pkg,  
    #Optional  
    'help'             => sub { pod2usage(-exitval => 0); }, 
    'man|m'            => sub { pod2usage(-exitval => 0, -verbose => 2); },
   );

  my @aligner_params = @ARGV;

#can't do this as we expect the unexpected! <Cue silouette of crazy naked dancing lady> <- 10 points for that reference
# or pod2usage(-exitval => 1, -message => "Specified parameters are:\t@tmp_args"); 
             
  if(! defined $aln_pkg){
    pod2usage(-exitval => 1, -message => '-aligner parameter must be defined');
  }
  
  if($aln_pkg !~ /::/){
    warn "Full $aln_pkg namespace was not specified, defaulting to:\tBio::EnsEMBL::Funcgen::Hive::$aln_pkg\n";
    $aln_pkg = "Bio::EnsEMBL::Funcgen::Hive::$aln_pkg";
  }
  
  #Quote so eval treats $aln_pkg as BAREWORD and converts :: to /
  if( ! eval "{require $aln_pkg; 1}"){
    die("Failed to require $aln_pkg\n$@");
  }

  if(! @aligner_params){
    pod2usage(-exitval => 1, 
              -message => "There are no aligner parameters specified\n".
              "For more information please run:\tperldoc $aln_pkg");
  }
 
  my $aligner = $aln_pkg->new(@aligner_params);

  if(! defined $aligner){
    die("Failed to create $aln_pkg from params:\n@aligner_params");
  }

  #This should return the output file
  return $aligner->run;
}

1;


__END__

