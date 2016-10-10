package Bio::EnsEMBL::Funcgen::HiveConfig::Ftp::ExportDataFiles;

use strict;
use warnings;
use base 'Bio::EnsEMBL::Funcgen::HiveConfig::Ftp::Base';
use Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf;

sub default_options {
    my ($self) = @_;
    
    return {
        %{$self->SUPER::default_options},
        'species_assembly_data_file_base_path' => [],
    }
}

sub pipeline_analyses {
    my $self = shift;

    my $data_freeze_date = $self->o('data_freeze_date');
    
    my $species_assembly_data_file_base_path_parameter = 'Not set';
    
    if (ref $self->o('species_assembly_data_file_base_path') eq 'ARRAY') {
      my $species_assembly_data_file_base_path_list = $self->o('species_assembly_data_file_base_path');
      my @species_assembly_data_file_base_path_parameters;
      foreach my $current_species_assembly_data_file_base_path (@$species_assembly_data_file_base_path_list) {
        push @species_assembly_data_file_base_path_parameters, "  --species_assembly_data_file_base_path $current_species_assembly_data_file_base_path";
      }
      $species_assembly_data_file_base_path_parameter = join "\\\n", @species_assembly_data_file_base_path_parameters;
    }
    
    return [
        {   -logic_name  => 'start_export',
            -module      => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
            -flow_into   => {
               MAIN => 'fetch_parameters_from_db'
            },
        },
        {   -logic_name  => 'fetch_parameters_from_db',
            -module      => 'Bio::EnsEMBL::Funcgen::Hive::Ftp::FetchParametersFromDb',
            -flow_into   => {
               MAIN => 'export_data_files_bam'
            },
        },
        {   -logic_name  => 'export_data_files_bam',
            -module      => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -parameters  => {
              cmd => '
                export_data_files.pl  \
                  --destination_root_path #ftp_base_dir#/#species#/Alignments  \
                  --file_type bam  \
                  --assembly #assembly# \
                ' . $species_assembly_data_file_base_path_parameter . ' \
                  --registry #reg_conf#  \
                  --die_if_source_files_missing 0 \
                  --data_freeze_date #data_freeze_date# \
                  --species #species#
                ',
            },
            -flow_into   => {
               MAIN => 'export_data_files_bigwig'
            },
        },
        {   -logic_name  => 'export_data_files_bigwig',
            -module      => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -parameters  => {
              cmd => '
                export_data_files.pl  \
                  --destination_root_path #ftp_base_dir#/#species#/Alignments  \
                  --file_type bigwig  \
                  --assembly #assembly# \
                ' . $species_assembly_data_file_base_path_parameter . ' \
                  --registry #reg_conf#  \
                  --die_if_source_files_missing 0 \
                  --data_freeze_date #data_freeze_date# \
                  --species #species#
                ',
            },
            -flow_into   => {
               MAIN => 'export_data_files_segmentation'
            },
        },
        {   -logic_name  => 'export_data_files_segmentation',
            -module      => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -parameters  => {
              cmd => '
                export_segmentation_files.pl            \
                  --registry #reg_conf#                 \
                  --species #species#                   \
                  --assembly #assembly#                 \
                  --data_freeze_date #data_freeze_date# \
                ' . $species_assembly_data_file_base_path_parameter . ' \
                  --destination_root_path #ftp_base_dir#/#species#/Segmentation  \
                ',
            },
        },
    ]
}

1;
