package Bio::EnsEMBL::Funcgen::Hive::RefBuildFileLocator;

use warnings;
use strict;
use Bio::EnsEMBL::Utils::Argument  qw( rearrange );

sub new {
  my $caller = shift;
  my $class  = ref($caller) || $caller;
  my $self   = bless {}, $class;
  return $self;
}

sub locate {
  my $self  = shift;
  my $param = shift;
  
  my $species          = $param->{species};
  my $assembly         = $param->{assembly};
  my $epigenome_gender = $param->{epigenome_gender};
  my $file_type        = $param->{file_type};
  
  my $directory_species_gender;
  
  if ($species eq 'homo_sapiens') {
    $directory_species_gender = $species . '_' . $self->epigenome_gender_to_directory_species_gender($epigenome_gender) . '/' . $species;
  } else {
    $directory_species_gender = $species;
  }
  
  if ($file_type eq 'bwa_index') {
    return $self->bwa_index_by_species_assembly($directory_species_gender, $assembly);
  }
  if ($file_type eq 'samtools_fasta_index') {
    return $self->samtools_fasta_index_by_species_assembly($directory_species_gender, $assembly);
  }
  if ($file_type eq 'chromosome_lengths_by_species_assembly') {
    return $self->chromosome_lengths_by_species_assembly($directory_species_gender, $assembly);
  }
  die;
}

sub chromosome_lengths_by_species_assembly {

  my $self     = shift;
  my $species  = shift;
  my $assembly = shift;

  return $species . '/' . $assembly . '/genome_fasta/' . $assembly . '.sizes';
}

sub bwa_index_by_species_assembly {

  my $self     = shift;
  my $species  = shift;
  my $assembly = shift;

  return $species . '/' . $assembly . '/genome_index/bwa/' . $assembly;
}

sub samtools_fasta_index_by_species_assembly {

  my $self     = shift;
  my $species  = shift;
  my $assembly = shift;

  return $species . '/' . $assembly . '/genome_fasta/' . $assembly . '.fa.fai';
}

=head2 epigenome_gender_to_directory_species_gender

   Currently only used for human.
   
=cut
sub epigenome_gender_to_directory_species_gender {

  my $self             = shift;
  my $epigenome_gender = shift;
  
  if ($epigenome_gender eq 'male') {
    return 'male';
  }
  return 'female';
}

1;
