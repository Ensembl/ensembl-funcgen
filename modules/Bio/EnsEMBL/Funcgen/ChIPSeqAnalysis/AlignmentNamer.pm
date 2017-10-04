package Bio::EnsEMBL::Funcgen::ChIPSeqAnalysis::AlignmentNamer;

use strict;
use File::Spec;

use Role::Tiny::With;
with 'Bio::EnsEMBL::Funcgen::GenericConstructor';

sub _constructor_parameters {
  return {
    directory_name_builder      => 'directory_name_builder',
    experiment                  => 'experiment',
    biological_replicate_number => 'biological_replicate_number',
    technical_replicate_number  => 'technical_replicate_number',
  };
}

use Bio::EnsEMBL::Funcgen::GenericGetSetFunctionality qw(
  _generic_get_or_set
);

sub directory_name_builder      { return shift->_generic_get_or_set('directory_name_builder',      @_); }
sub biological_replicate_number { return shift->_generic_get_or_set('biological_replicate_number', @_); }
sub technical_replicate_number  { return shift->_generic_get_or_set('technical_replicate_number',  @_); }
sub experiment                  { return shift->_generic_get_or_set('experiment',                  @_); }
#sub name                        { return shift->_generic_get_or_set('name',                        @_); }

sub unset_biological_replicate_number {
  shift->{biological_replicate_number} = undef;
}

sub unset_technical_replicate_number {
  shift->{technical_replicate_number} = undef;
}

sub name {

  my $self       = shift;
  
  my $experiment                  = $self->experiment;
  my $biological_replicate_number = $self->biological_replicate_number;
  my $technical_replicate_number  = $self->technical_replicate_number;
  
  my $epigenome          = $experiment->epigenome;
  my $feature_type       = $experiment->feature_type;
  my $experimental_group = $experiment->experimental_group;

  my @name_components = (
    $epigenome->production_name,
    $feature_type->name,
    $experimental_group->name,
  );
  
  if ($biological_replicate_number) {
    push @name_components, 'BR' . $biological_replicate_number;
  }
  if ($technical_replicate_number) {
    push @name_components, 'TR' . $technical_replicate_number;
  }
  
  my $alignment_name = join '_', @name_components;
  
  use Bio::EnsEMBL::Funcgen::Utils::GoodUtils qw( create_production_name );
  return create_production_name($alignment_name);
}

sub base_name_with_duplicates {
  my $self = shift;
  return $self->name . '_with_duplicates';
}

sub base_name_no_duplicates {
  my $self = shift;
  return $self->name . '_no_duplicates';
}

sub bam_file_name_with_duplicates {
  my $self = shift;
  return $self->base_name_with_duplicates . '.bam';
}

sub bam_file_name_no_duplicates {
  my $self = shift;
  return $self->base_name_no_duplicates . '.bam';
}

sub bigwig_file_name {
  my $self = shift;
  return $self->base_name_no_duplicates . '.bw';
}

sub bam_file_with_duplicates {
  my $self = shift;
  return File::Spec->catfile(
    $self->directory_name_builder->bam_output_dir,
    $self->bam_file_name_with_duplicates
   );
}
sub bam_file_with_duplicates_stored {
  my $self = shift;
  return File::Spec->catfile(
    $self->directory_name_builder->bam_output_dir_stored,
    $self->bam_file_name_with_duplicates
   );
}
sub bam_file_no_duplicates {
  my $self = shift;
  return File::Spec->catfile(
    $self->directory_name_builder->bam_output_dir,
    $self->bam_file_name_no_duplicates
   );
}
sub bam_file_no_duplicates_stored {
  my $self = shift;
  return File::Spec->catfile(
    $self->directory_name_builder->bam_output_dir_stored,
    $self->bam_file_name_no_duplicates
   );
}
sub bigwig_file_no_duplicates {
  my $self = shift;
  return File::Spec->catfile(
    $self->directory_name_builder->bigwig_output_dir,
    $self->bigwig_file_name
   );
}
sub bigwig_file_no_duplicates_stored {
  my $self = shift;
  return File::Spec->catfile(
    $self->directory_name_builder->bigwig_output_dir_stored,
    $self->bigwig_file_name
   );
}

1;

