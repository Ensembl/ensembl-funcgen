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
    OR 1 or both of these
    -split                  
    -merge                          
                              
                              
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

# To do
# Move more functionality in here? Split and merge? Pass the helper for debug?
# No not here, but called from here, functionality should go in a Sequencing::Utils module
# Move All Aligners and PeakCallers to Sequencing
# Have code return info necessary for hive data flow
# Test run as library mode? Will @ARGV be inherited from a calling script?
# Leave dataflow and bsubing fan jobs to the caller
# Add main back to handle run modes and param parsing from cmdline
# Then other subs do not have to handle GetOptions at all.
# Hmm, is this a good idea to have these method requireable?
# These will use the global namespace i.e. $main::run_aligner
# This will redefine any other global of the same name
# Really, it is much safer to put as much of this as possible in a module
# Can we put all of the logic below in Aligner? 
# Or RunAligner, this will then delegate to the actual aligner as required.
# This RunAligner could then be used in the Hive::RunAligner
# This script could then either use RunAligner or SeqUtils as required.  
# The aim is to remove the need to require this script at all.
# Maybe this control class should be AlignTools?
# Which can be used as an object or a package

use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use Bio::EnsEMBL::Funcgen::Sequencing::SeqTools;# qw( run_aligner );
use Bio::EnsEMBL::Funcgen::Sequencing::SeqQC; #in case we just want to run QC

Getopt::Long::Configure("pass_through"); #Allows unknown options to pass on to @ARGV

#As pod2usage depends on trying to open $0 to access the pod and $0 is set to 
#the parent path, then it will fail when included from: 
#   perl -e 'require "run_aligner.pl"' 
#with:
#   Can't open -e: No such file or directory
#Will likely list the POD of the caller script when included into a defined script/module.
#Hence resetting $0 = __FILE__ which is definitely the real path of run_aligner.pl
#This is done prior to each pod2usage call to reduce the risk of
#adverse effects, but may still cause problems if this is run under eval
#or if there or DESTROY methods in calling objects which access $0.
#Will this affect stack traces???
#Other option here is to test caller and use throw or pod2usage as required.
#($0, the Cwd module and the FindBin module) fail under mod_perl

#my $real_path = __FILE__;

main() unless(caller);
#This avoids executing run_aligner when including as a libaray via require or do
#for general usage or when unit testing
#Although this is now slightly redundant now we have moved things to AlignTools


#process_sam_bam filters unmapped, sorts and rmdups
#merge only rmdups

#todo use helper here, then pass that through for debugging?
#this will prevent having to handle debug in SeqTools, and then in Aligner?


sub main{
  my ($aligner, $split, $merge, $outfile, $skip_qc, $no_rmdups, $qc_type, $out_dir, $batch_job,
      $sam_ref_fai, @files, $debug_level, $write_chk, $out_prefix, @qc_methods);

  GetOptions 
   (#Run modes
    'split'       => \$split,
    'aligner=s'   => \$aligner,  
    #'filter'      => \$filter, #process_sam_bam
    'merge'       => \$merge,
    
    #'peak_caller=s'
    
    
    'out_dir=s' => \$out_dir,

    #Split/Merge params
    'files=s{,}'   => \@files,
    'out_prefix=s' => \$out_prefix, 
    
    #Merge params
    'outfile=s'       => \$outfile,
    'sam_ref_fai=s'   => \$sam_ref_fai, #also required for BWA aligner
    'write_checksum'  => \$write_chk,
    'no_rmdups'       => \$no_rmdups,
    'qc_methods=s{,}' => \@qc_methods,
    'qc_type=s'       => \$qc_type,
    
    #these should return a hash ref of qc method keys
    #and qc results values
    #These will be used to push into the tracking DB
    #just store the pass/fail state and point to the output file for now?
    
    
     
    #Optional
    'skip_inline_qc' => \$skip_qc,
    'debug=i'   => \$debug_level,
    'help'      => sub { pod2usage(-exitval => 0); }, 
    'man|m'     => sub { pod2usage(-exitval => 0, -verbose => 2); },
   );
  # or pod2usage(-exitval => 1, -message => "Specified parameters are:\t@tmp_args"); 
  #can't do this as we expect the unexpected! <Cue silouette of crazy naked dancing lady> <- 10 points for that reference
  my @extra_params = @ARGV;
  
  #Passing ARGV like this will not handle flag, as they will have no value
  #so all unsupported flags must have boolean values
  #supported flags tend to be run modes
  #so we have no need for the constructor at all
  #if we are handling run modes here
  #This is the same for arrays. These will have to be passed as quoted string
  #and handled in AlignTools i.e. if they contain spaces, then they will be split where appropriate
 
 
  if(! ($aligner || $merge || $split || $qc_type || @qc_methods )){
    pod2usage(-exitval => 1, 
              -message => "Run mode has not been specified.\n".
                "Please specify one of:\t".join(' ', (qw(-aligner -split -merge -qc_type -qc_methods))).
                "\nOr -split and -merge together");  
  }
  
  
  #Need to handle multiple modes here
  
  if($aligner &&
     ($merge && $split)){
    pod2usage(-exitval => 1, 
              -message => "You have defined mutually exclusive run modes.\n".
                "Please specify -aligner, -split or -merge individually or -split and -merge together");          
  }
  
  
  #Here we have a problem, that some of these method may take hashrefs of params
  #(or other objects) which cannot be passed on the command line
  #hence will need to add option support for those params
 
  #Change these to named subs or despatch table?
 
  
  if($aligner){
    #we may have clash between -files and -query_file
    
    if(@files){
       pod2usage(-exitval => 1, 
                 -message => "The -files parameter is not valid for the -aligner mode.\n".
                  "Please specify -query_file");
    }
    
    #sam_ref_fai is BWA specific here
    #need to make this generic somehow

   warn "aligner params are @extra_params";


    run_aligner(-aligner        => $aligner, 
                -skip_qc         => $skip_qc,
                -aligner_params => [@extra_params,
                                    -sam_ref_fai, $sam_ref_fai]);
  }
  elsif($split){
    my ($fastq_chunks, $qc_results) = 
      split_fastqs(-files      => \@files, 
                   -out_prefix => $out_prefix,
                   -out_dir    => $out_dir,
                   -merge      => $merge, 
                   -skip_qc    => $skip_qc);
                   
                   
    #Do something for QC here?
                 
    #output cmdline for run_aligner and merge
    #enable batch job run?
    #This would be unsafe to automatically use LSB_INDEX to generate suffix, as this
    #may be part of a batch, unless done explcitly with a -batch_job flag in this script
    #which then would change the query file accordingly
    #This would be tricky, as we don't capture -query_file in here
    #Just pass -batch_job to run_aligner
    
    my $batch_size = scalar(@$fastq_chunks);
    (my $prefix    = $fastq_chunks->[0]) =~ s/_[0-9]{4}$//;
    (my $job_name  = $prefix)            =~ s/.*\///;
    
    print "To run aligner jobs use the following cmdline:\n\t".
      "bsub -o ${job_name}_alignments.\%J.\%I.out -e ${job_name}_alignments.\%J.\%I.err -J ${job_name}_alignments[1-$batch_size] -q normal -R ".
      '$EFG_SRC/scripts/sequencing/run_seq_tools.pl -aligner BWA -program_file bwa -batch_job 1 -sam_ref_fai $SAM_REF_FAI -target_file $TARGET_FILE '.
      " -query_file $prefix\n";
     
    #-w is lsf job dependancy on align batch
    #need to change the file to bam files
        
    print "To run post alignment merge job, use the following cmdline:\n\t".
      "bsub -w \"'done(${job_name}_alignments)'\" \$EFG_SRC/scripts/sequencing/run_seq_tools.pl -merge -sam_ref_fai \$SAM_REF_FAI ".
      "-outfile ${out_prefix}.bam -files ".join(' ', @$fastq_chunks);   
  }
  else{#must be $merge
    #
    #param hash debug, rm_dups, write_checksum
    
    #This is now directly calling merge_bams!
    #Maybe this should move to AlignTools
    #The idea is to move all functionality to AlignTools
    #so this is the one stop shop for Seq/Align functions
    #SeqUtils?
    
    #This script will be the cmdline interface to all those 
    #high level utils
    
    #SeqUtils need to include EFGUtils
    #The Aligners and PeakCaller cannot include SeqUtils
    #other wise there will be a circular dependancy
    #
    
    
    merge_bams($outfile, $sam_ref_fai, \@files, 
               {debug          => $debug_level,
                no_rmdups      => $no_rmdups,
                write_checksum => $write_chk});
  }
   
}




1;


__END__

