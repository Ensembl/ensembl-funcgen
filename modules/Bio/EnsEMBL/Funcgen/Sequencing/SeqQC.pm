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

=cut


=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

  Questions may also be sent to the Ensembl help desk at
  <http://www.ensembl.org/Help/Contact>.

=head1 NAME

Bio::EnsEMBL::Funcgen::Sequencing::SeqQC

=head1 DESCRIPTION

This module collates a variety of quality control utilities for sequence alignments
and peaks calling


=head1 SYNOPSIS

use Bio::EnsEMBL::Funcgen::Utils::Sequencing qw(list your required methods here);



=cut

###############################################################################

package Bio::EnsEMBL::Funcgen::Sequencing::SeqQC;

use warnings;
use strict;

use File::Find; 
use Bio::EnsEMBL::Utils::Argument  qw( rearrange );
use Bio::EnsEMBL::Utils::Exception qw( throw      );
use Bio::EnsEMBL::Utils::Scalar    qw( assert_ref );
use Bio::EnsEMBL::Funcgen::Utils::EFGUtils qw( run_system_cmd open_file );
#use File::Path                     qw( make_path  );
#use File::Basename                 qw( dirname fileparse );
#use File::Spec;
#use FileHandle;

use base qw( Exporter );
use vars qw( @EXPORT_OK );

@EXPORT_OK = qw( 
  run_QC
  run_spp_phantom
  );



#todo
#implement so infrastructure to for running all default QC for given stage
#will likely need a hash with some label keys:
# fastq, raw_bam, filtered_bam, peaks etc.
#Write wrapper methods for these if required, such that we don't need to handle option
#setting

#Return these methods as a list, such that they can either be run in series or submitted to the farm to run in parallel

#Input params should be standardised across methods with same qc type/label


sub run_QC{
    
  
  
}







=head2 run_spp_phantom

  Arg [1]     : String  - input.bam file path
  Arg [2]     : Hashref - Params: 
                  {-out_prefix => <String>,   #Output file prefix. Default = input.bam
                   -quality    => <Int>,      #Acceptable quality threshold [-2..2]
                   -out_dir    => <String>,   #Output directory(-odir), default is 
                                              directory of input bam file. 
                   -no_dups    => <Boolean>,  #For bam files where duplicates have been removed.
                                              Runs run_spp_nodups.R instead of run_spp.R.
                   -force      => <Boolean>}  #Force flag. Over-writes previous output
                                              by spcifying -rf
  Arg [3]     : Hashref - Raw SPP params. These can be any valid SPP command line
                options, with optional values. These will over-ride any of the above 
                params. The '=' does not need to be specified for key value pairs 
                e.g.
                  -odir=/your/out/put/directory
                  Should be specified as follows:
                  {'-odir' => '/your/out/put/directory'}             
  Description : A wrapper method to run the phantompeakqualtools of the SPP package.
                Defaults outputs will be:
                  $output_dir/$output_prefix.spp.txt
                  $output_dir/$output_prefix.spp.pdf
                
                Note: SPP phantom uses quite a lot of tmp space, roughly the same 
                amount of the input bam file. Hence it may be necessary to specify this 
                in any resource requirements. 
                  e.g. For LSF: bsub 
                  
                For more details about command line options and requirements see:
                  https://code.google.com/p/phantompeakqualtools/
  Returntype  : Int - SPP Quality tag [-2..2]
  Example     : my $qc_passed_results = run_spp_phantom($bam_file, {quality => 1}); 
  Exceptions  : Throws if input.bam not valid.
                Throws if acceptable quality threshold is defined and not met.
  Status      : At risk

=cut


#unix,bash,R-2.10 and above,awk,samtools,boost C++ libraries
#R packages: SPP, caTools, snow
#NOTE: The current package does not run on a MAC or WINDOWS.

#run post filtering with a threshold of 1 or higher
#if we get failures then we can always go back and look at the unfitlered bams if need be
#to look at pure library quality.


#8GB? Probably a lot less than that.

#Max 
#2MB for 459M

#We may want to define an RScript path to override the default
#This can be done outside of perl by appropriately order the $PATH
#But this may not always be possible

sub run_spp_phantom{
  my $in_bam     = shift;
  my $params     = shift || {};
  my $raw_params = shift || {}; #overrides $params
  my $result = 1;
  
  #VALIDATE PARAMS
  if(! defined $in_bam ){
    throw('Bam file argument is not defined');      
  }
  elsif(! -f $in_bam){
    throw("Bam file argument is not a file:\t".$in_bam);
  }
  
  assert_ref($params, 'HASH');
  assert_ref($raw_params, 'HASH');
  
  #PROCESS PARAMS
  #validate all param key start with -, 
  #otherwise rearrange will simply return them (hangover from CGI origins)
  #reported to core team
  
  if(my @non_hyphenated = grep {/^[^-]/} keys %$params){
    throw("Params must be prefixed with a hyphen:\t@non_hyphenated");  
  }
  
  my ($out_pfx, $req_qual, $out_dir, $no_dups, $force) = 
    rearrange([qw(out_prefix quality out_dir no_dups force )], %$params);
  
  #Generate defaults, spp would use these anyway, 
  #but we need them in the error output for clarity
  ($out_pfx = $in_bam) =~ s/.*\///o     if ! defined $out_pfx;
  ($out_dir = $in_bam) =~ s/$out_pfx//  if ! defined $out_dir;
  
  #SET DEFAULT OPTIONS
  if(! keys(%$raw_params)){ #Set default save option
    $raw_params->{'-savp'} = "${out_dir}/${out_pfx}.spp.pdf"; 
  }
  
  if(! exists $raw_params->{'-odir'}){
    #This maybe an empty string if we have passed a bam file in the cwd
    $raw_params->{'-odir'} = $out_dir if $out_dir;   
  }

  if( ! exists $raw_params->{'-out'}){
    $raw_params->{'-out'} = "${out_dir}/${out_pfx}.spp.txt";      
  }
  
  $raw_params->{'-rf'} = undef if $force;
    
  #BUILD THE CMD 
  my $cmd = $no_dups ? 'run_spp_nodups.R' : 'run_spp.R';
  $cmd    = _which_path($cmd);
  $cmd    = "Rscript $cmd -c=${in_bam}";
  
  foreach my $opt(keys(%$raw_params)){
    my $param = ' '.$opt;
    $param   .= '='.$raw_params->{$opt} if defined $raw_params->{$opt};
    $cmd     .= $param;
  }
  
  
  #It's a good idea to test the outputs will be writeable before wasting time
  # 
  
  
  #RUN THE CMD RETURN RESULTS
  #warn "cmd is $cmd";
  run_system_cmd($cmd);
  #could really do with passing error string back if this fails
  
  
  my $fh = open_file($raw_params->{'-out'});
  my $output = <$fh>;
  $fh->close;   
  chomp $output;
  my @output = split("\t", $output);
  #$output[0] is filename  
  my %results = (numReads         => $output[1],
                 estFragLen       => $output[2],
                 corr_estFragLen  => $output[3],
                 phantomPeak      => $output[4],
                 corr_phantomPeak => $output[5],
                 argmin_corr      => $output[6],
                 min_corr         => $output[7],
                 NSC              => $output[8],
                 RSC              => $output[9],
                 QualityTag       => $output[10]);
  
  $out_dir = $raw_params->{'-out'}.'/' || '';
  
  #todo change this so we don't throw, but we return a fail status
  
           
  if(defined $req_qual && 
     ($results{QualityTag} < $req_qual)){
    $result = 0;   
    warn("SPP:PhantomPeakQualTools failed on QualityTag:\tRequired $req_qual\tObserved ".
      $results{QualityTag}."\nInput file:\t$in_bam\nResults:\t${out_dir}${out_pfx}.spp.txt\n".
      "Plot:\t${out_dir}/${out_pfx}.spp.pdf\n");     
  }
  
  #return like as in the scalar context the last value of the array is returned
  #meaning we can do something like this in the caller
  #if(! run_spp_phantom()){ #do something }
  
  return \%results, $result;
}


#-f seems to be optional, and only acts to over-ride auto-detection

#fastqc online docs seem focused on gui usage only, 
#cmdline options are not available on line only via man page?
#

sub run_fastqc{
  my $in_files = shift;
  my $params   = shift || {};
  
  if(! defined $in_files){
    throw('Input file(s) argument not defined, mus be a single file path or an array ref of file paths');  
  }
  elsif(ref($in_files)){
    assert_ref($in_files, 'ARRAY', 'fastqc input files');  
  }
  else{
    $in_files = [$in_files];  
  }
  
  assert_ref($params, 'HASH');

  if(my @non_hyphenated = grep {/^[^-]/} keys %$params){
    throw("Params must be prefixed with a hyphen:\t@non_hyphenated");  
  }
  
  
  my ($out_dir) = 
    rearrange([qw(out_dir format)], %$params);
 
  #fastqc isn't normally installed into the top of the path 
  #a symlink is required, so let's do some whiching to see whether we have it in the $PATH
  #else call _which_path?
    
  #my $cmd = ;
  #warn "DEACTIVATED FASTQC FOR NOW:\nfastqc -f fastq -o ".$self->output_dir." @fastqs";
  #run_system_cmd('fastqc -o '.$self->output_dir." @fastqs"); 
  
  
  return;
}

#This is to get around the problem that some scripts are not available in the top level of the $PATH dirs
#and creating a toplevel link is not appropriate as we need there full path.
#move this to EFGUtils?
sub _which_path{
  my $filename  = shift;
  my @env_paths = split(/:/, $ENV{PATH});
  my $path;
   
  find(sub{ if ($_ eq $filename){$path=$File::Find::name; return; }}, @env_paths);
  
  if(! defined $path){
    throw("Unable to find file in \$PATH:\t$filename\n\$PATH contains:\t".join("\t", @env_paths));
  }
  
  return $path;
}

1;