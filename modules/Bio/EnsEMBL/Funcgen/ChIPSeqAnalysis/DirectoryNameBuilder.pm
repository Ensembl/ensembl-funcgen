package Bio::EnsEMBL::Funcgen::ChIPSeqAnalysis::DirectoryNameBuilder;

use strict;
use Data::Dumper;

use Role::Tiny::With;
with 'Bio::EnsEMBL::Funcgen::GenericConstructor';

sub _constructor_parameters {
  return {
    root_dir                => 'root_dir',
    species                 => 'species',
    assembly                => 'assembly',
    ensembl_release_version => 'ensembl_release_version',
  };
}

use Bio::EnsEMBL::Funcgen::GenericGetSetFunctionality qw(
  _generic_get_or_set
);

sub root_dir                { return shift->_generic_get_or_set('root_dir',                @_); }
sub species                 { return shift->_generic_get_or_set('species',                 @_); }
sub assembly                { return shift->_generic_get_or_set('assembly',                @_); }
sub ensembl_release_version { return shift->_generic_get_or_set('ensembl_release_version', @_); }

=head1 regulation_directory

  The directory in which the regulation files go. This is shared with the 
  other teams, so best to put it in its own directory.
  
  Deactivated for now, use own directory after the development break when
  we move from result_set to alignment.

=cut
sub regulation_directory {
  #return undef;
  return 'funcgen';
}

=head1 version_directory

  The Ensembl release version when the files were generated. Padded with a 
  zero so it looks nice when we hit release 100.

=cut
sub version_directory {
  my $self = shift;
  return $self->ensembl_release_version;
}

sub db_output_dir {
  my $self = shift;
  return File::Spec->catfile(
    $self->root_dir,
    $self->species,
    $self->assembly,
  );
}

=head1 default_directory_by_table_and_file_type

  Returns the names of directories in which the files should be stored by 
  default. This is:
  
  <root directory> / <the directory forregulation specific files> / <the name of the table that represents it in the database> / < the release version > / ersa_signal / <the file type>

=cut
sub default_directory_by_table_and_file_type {
  my $self   = shift;
  my @params = @_;
  
  return File::Spec->catfile(
    $self->db_output_dir,
    $self->default_directory_by_table_and_file_type_stored(@params),
  );
}

sub default_directory_by_table_and_file_type_stored {
  my $self      = shift;
  my $table     = shift;
  my $file_type = shift;
  
  return '/' . File::Spec->catfile(
    $self->regulation_directory,
    $table, 
    $self->version_directory, 
    'ersa_signal', 
    $file_type
  );
}

sub peaks_output_dir         { return shift->default_directory_by_table_and_file_type        ('peak',      'peaks');  }
sub bigwig_output_dir        { return shift->default_directory_by_table_and_file_type        ('alignment', 'bigwig'); }
sub bigwig_output_dir_stored { return shift->default_directory_by_table_and_file_type_stored ('alignment', 'bigwig'); }
sub bam_output_dir           { return shift->default_directory_by_table_and_file_type        ('alignment', 'bam');    }
sub bam_output_dir_stored    { return shift->default_directory_by_table_and_file_type_stored ('alignment', 'bam');    }

=head1 quality_check_output_dir

  The root directory for the output of all quality checks.

=cut
sub quality_check_output_dir {
  my $self = shift;
  return File::Spec->catfile(
    $self->db_output_dir, 
    'quality_checks', 
    $self->version_directory
  )
}

sub fastqc_output_dir                       { File::Spec->catfile( shift->quality_check_output_dir, 'fastqc')                       }
sub flagstats_output_dir                    { File::Spec->catfile( shift->quality_check_output_dir, 'flagstats')                    }
sub phantom_peaks_output_dir                { File::Spec->catfile( shift->quality_check_output_dir, 'phantom_peaks')                }
sub proportion_of_reads_in_peaks_output_dir { File::Spec->catfile( shift->quality_check_output_dir, 'proportion_of_reads_in_peaks') }
sub chance_output_dir                       { File::Spec->catfile( shift->quality_check_output_dir, 'chance')                       }

1;
