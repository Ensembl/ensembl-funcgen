package Bio::EnsEMBL::Funcgen::PeakCallingPlan::AlignmentNamer;

use strict;
use File::Spec;
use Bio::EnsEMBL::Utils::Exception qw( throw );

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

sub init {
  
  my $self = shift;
  
  use Carp;
  if (! defined $self->experiment) {
    confess("Experiment is a mandatory parameter!");
  }
  return;
}

use Bio::EnsEMBL::Funcgen::GenericGetSetFunctionality qw(
  _generic_get_or_set
);

sub directory_name_builder      { return shift->_generic_get_or_set('directory_name_builder',      @_); }
sub biological_replicate_number { return shift->_generic_get_or_set('biological_replicate_number', @_); }
sub technical_replicate_number  { return shift->_generic_get_or_set('technical_replicate_number',  @_); }
sub experiment                  { return shift->_generic_get_or_set('experiment',                  @_); }

sub unset_biological_replicate_number {
  shift->{biological_replicate_number} = undef;
}

sub unset_technical_replicate_number {
  shift->{technical_replicate_number} = undef;
}

sub name {

  my $self       = shift;
  my $experiment = $self->experiment;
  
  if (! defined $experiment) {
    throw("No experiment passed!");
  }
  
  my $biological_replicate_number = $self->biological_replicate_number;
  my $technical_replicate_number  = $self->technical_replicate_number;
  
  my $epigenome          = $experiment->epigenome;
  my $feature_type       = $experiment->feature_type;
  my $experimental_group = $experiment->experimental_group;

  my @name_components = (
    $experiment->dbID,
    $epigenome->production_name,
    $feature_type->name,
    $experimental_group->name,
  );
  
  my $alignment_name = join '_', @name_components;
  
  use Bio::EnsEMBL::Funcgen::Utils::GoodUtils qw( create_production_name );
  my $production_name = create_production_name($alignment_name);

  if ($biological_replicate_number) {
    $production_name .= '_BR' . $biological_replicate_number
  }
  if ($technical_replicate_number) {
    $production_name .= '_TR' . $technical_replicate_number
  }
  return $production_name;
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

sub bed_file_name {
  my $self = shift;
  return $self->base_name_no_duplicates . '.bed';
}

sub bam_file_with_duplicates {
  my $self = shift;
  return File::Spec->catfile(
    $self
      ->directory_name_builder
      ->bam_output_dir(
        $self->experiment
      ),
    $self->bam_file_name_with_duplicates
   );
}
sub bam_file_with_duplicates_stored {
  my $self = shift;
  return File::Spec->catfile(
    $self
      ->directory_name_builder
      ->bam_output_dir_stored(
        $self->experiment
      ),
      $self->bam_file_name_with_duplicates
   );
}
sub bam_file_no_duplicates {
  my $self = shift;
  return File::Spec->catfile(
    $self
      ->directory_name_builder
      ->bam_output_dir(
        $self->experiment
      ),
    $self->bam_file_name_no_duplicates
   );
}
sub bam_file_no_duplicates_stored {
  my $self = shift;
  return File::Spec->catfile(
    $self
      ->directory_name_builder
      ->bam_output_dir_stored(
        $self->experiment
      ),
    $self->bam_file_name_no_duplicates
   );
}
sub bigwig_file_no_duplicates {
  my $self = shift;
  return File::Spec->catfile(
    $self
      ->directory_name_builder
      ->bigwig_output_dir(
        $self->experiment
      ),
    $self->bigwig_file_name
   );
}
sub bed_file_no_duplicates {
  my $self = shift;
  return File::Spec->catfile(
    $self
      ->directory_name_builder
      ->bed_output_dir(
        $self->experiment
      ),
    $self->bed_file_name
   );
}
sub bigwig_file_no_duplicates_stored {
  my $self = shift;
  return File::Spec->catfile(
    $self
      ->directory_name_builder
      ->bigwig_output_dir_stored(
        $self->experiment
      ),
    $self->bigwig_file_name
   );
}

1;

