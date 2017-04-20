=pod 

=head1 NAME

Bio::EnsEMBL::Hive::RunnableDB::Funcgen::CollectionWriter

=head1 DESCRIPTION

=cut

package Bio::EnsEMBL::Funcgen::Hive::CollectionWriter;

use base ('Bio::EnsEMBL::Funcgen::Hive::BaseDB');

use strict;
use Bio::EnsEMBL::Funcgen::Utils::EFGUtils qw( run_system_cmd );

sub run {
  my $self = shift; 

  my $db = $self->param_required('out_db');
  my $result_set_adaptor = $db->get_ResultSetAdaptor;
  $result_set_adaptor->{file_type} = 'BAM';

  my $result_set = $self->fetch_Set_input('ResultSet');
  my $bam_file = $self->db_output_dir . '/' . $result_set->dbfile_path;
  
  my @file_to_delete;

  if(! -e $bam_file) {
    use Carp;
    confess("Can't find bam file $bam_file!");
  }

  my $expected_bam_index_file = $bam_file . '.bai';
  
  if(! -e $expected_bam_index_file) {
    my $cmd = 'samtools index ' . $bam_file;
    run_system_cmd($cmd);
  }
  
  push @file_to_delete, $expected_bam_index_file;
  
  my $feature_set = $self->FeatureSet;
  
  foreach my $current_file (@file_to_delete) {
    $self->dataflow_output_id( {
        file_to_delete => $current_file,
    }, 7);
  }
  
  my $feature_set_analysis_logic_name = $feature_set->analysis->logic_name;

  my $output_id = {
    %{$self->batch_params}, 
    # These are already param_required by fetch_Set_input
    dbID         => $self->param('dbID'),
    set_name     => $self->param('set_name'),  # mainly for readability
    set_type     => $self->param('set_type'),
    filter_from_format => undef,
    logic_name => $feature_set_analysis_logic_name,
  };
  $self->dataflow_output_id( $output_id, 2 );
  return;
}

1;
