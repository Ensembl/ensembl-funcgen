package Bio::EnsEMBL::Funcgen::HiveConfig::Ftp::ExportDataFilesBase;

use strict;
use warnings;
use base 'Bio::EnsEMBL::Funcgen::HiveConfig::Ftp::Base';
use Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf;

sub default_options {
    my $self = shift;
    
    return {
        %{$self->SUPER::default_options},
        'dbfile_registry_path' => [],
    }
}

sub create_dbfile_registry_path_parameter {
  my $self = shift;
    my $dbfile_registry_path_parameter = 'Not set';
    
    my $dbfile_registry_path;
    if ($self->o('dbfile_registry_path') !~ /^#/) {
      if (ref $self->o('dbfile_registry_path') eq 'ARRAY') {
        $dbfile_registry_path = $self->o('dbfile_registry_path');
      }
    }
    
    if ($dbfile_registry_path) {
      my $dbfile_registry_path_list = $dbfile_registry_path;
      my @dbfile_registry_path_parameters;
      foreach my $current_dbfile_registry_path (@$dbfile_registry_path_list) {
        push @dbfile_registry_path_parameters, "  --dbfile_registry_path $current_dbfile_registry_path/#species#/#assembly#";
      }
      $dbfile_registry_path_parameter = join "\\\n", @dbfile_registry_path_parameters;
    }
    return $dbfile_registry_path_parameter;
}

1;
