=pod

=head1 NAME

Bio::EnsEMBL::Funcgen::Hive::BWA

=head1 DESCRIPTION


=cut

package Bio::EnsEMBL::Funcgen::Hive::BWA;

use warnings;
use strict;

#qw the methods even if they are EXPORTED, so we know where they come from
use Bio::EnsEMBL::Utils::Argument          qw( rearrange );
use Bio::EnsEMBL::Utils::Exception         qw( throw );
use Bio::EnsEMBL::Funcgen::Utils::EFGUtils qw( run_system_cmd );

use base qw( Bio::EnsEMBL::Funcgen::Hive::Aligner );#Does not import


#todo use Rob Davies hardened bwa https://github.com/daviesrob/bwa
#The following branches are available:
#0.5.10-mt_fixes    0.5.10 with multithreaded samse/sampe
#0.6.2_fixes        0.6.2 released version
#0.6.2-mt           0.6.2 with multithreaded samse/sampe
#master_fixes       The current bwa master version
#or see if Heng has incoprorated these into the main bwa github rep?

#todo look into multi-threading to speed thing up -t #nodes
#will this just reduce our submission rate as we will require an entire node
#to be free before we get a slot
#maybe try -t 2?

#sub new {
#  my $caller = shift;
#  my $class = ref($caller) || $caller;
#  my $self = $class->SUPER::new(@_);  
#  return;
#}



#todo handle bwa dir ans samtools dir better
#we could take a bin_dir and a program name/file
#todo handle samse/sampe? Default to samse
#this isn't really a parameter, and there are several commands to provide params for


#although output format is sam
#we always need to sort and convert to bam for merge

sub run {
  my $self = shift;

  my $query_file  = $self->query_file;
  my $input_dir   = $self->input_dir;
  my $output_dir  = $self->output_dir;
  my $ref_file    = $self->reference_file;
  my $bwa_bin     = $self->program_file;
  (my $bin_dir = $bwa_bin) =~ s/(.*\/)[^\/].*/$1/go;
  (my $outfile_prefix = $query_file) =~ s/\.fastq$//;
  
  #assume samtools is in the same dir as bwa
  
 

  #TODO Pass the location of the binary to be sure we'll be running the right version?
#  my $bwa_cmd = "$bwa_bin aln $bwa_index $input_file";
  #Allow this to work with paired reads?? Maybe not for the moment...
  #in that case pass bwa algorithm as parameter...
  #If using -q make sure we've got the correct Sanger quality scores...
#  $bwa_cmd .= " | $bwa_bin samse $bwa_index - $input_file";
#  $bwa_cmd .= " | samtools view -uS - ";
#  $bwa_cmd .= " | samtools sort - ${input_file}.sorted";
#  if(system($bwa_cmd) != 0){ throw "Problems running $bwa_cmd";  }


  #Piping all these together was causing uncaught errors

  ### FIND SUFFIX ARRAY COORDS OF SEQS
  my $bwa_cmd = "$bwa_bin aln $ref_file ${input_dir}/${query_file} > ".
    "${output_dir}/${outfile_prefix}.sai";
  run_system_cmd($bwa_cmd);
  

  ### GENERATE SINGLE END READ ALIGNMETNS
  #todo add -P option here to load index into memory which should speed execution
  #but will require ~5GB of memory
  #THIS ONLY WORKS FOR SAMPE?
  $bwa_cmd = "$bwa_bin samse $ref_file ${output_dir}/${outfile_prefix}.sai ${input_dir}/${query_file}".
    " > ${output_dir}/${outfile_prefix}.samse.sam";
  run_system_cmd($bwa_cmd);
  run_system_cmd("rm -f ${output_dir}/${outfile_prefix}.sai"); #No fail flag here?
  run_system_cmd("rm -f ${input_dir}/${query_file}"); #No fail flag here?
  
  
  #Now we need to change this to unfiltered? But do this in the merge
  #also, need utilise EFGUtils filtering and conversion?
  #filtering should be done here, such that we don't get parellel processes trying to do this for controls
  
  ### CONVERT TO BAM
  #-S input is sam, with header (else -t is required)
  #-u output is uncompressed bam, prefered for piping to other samtools cmds, although this has been shown to
  #be fragile, hence we have split the cmds up here
  $bwa_cmd = "${bin_dir}/samtools view -uS ${output_dir}/${outfile_prefix}.samse.sam > ${output_dir}/${outfile_prefix}.samse.bam.unsorted"; 
  run_system_cmd($bwa_cmd);
  run_system_cmd("rm -f ${output_dir}/${outfile_prefix}.samse.sam");

  ### SORT BAM
  $bwa_cmd = "${bin_dir}/samtools sort  ${output_dir}/${outfile_prefix}.samse.bam.unsorted ${output_dir}/${outfile_prefix}.bam";
  run_system_cmd($bwa_cmd);
  run_system_cmd("rm -f ${output_dir}/${outfile_prefix}.samse.bam.unsorted");
  
  #how are we going to handle samse here?
  #is this generic to all read aligners? nomenclature certainly isn't
  #merge job will need to know what the output format is.
  #so we need to make that available
  #although we are doing this here, as it is a pre-req for the merge step
  #just drop samse from file name for now, as we only do single end reads at present
  
  return;
}




1;
