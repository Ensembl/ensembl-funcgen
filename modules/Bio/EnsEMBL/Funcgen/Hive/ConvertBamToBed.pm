=pod 

=head1 NAME

Bio::EnsEMBL::Hive::RunnableDB::Funcgen::ConvertBamToBed

=head1 DESCRIPTION

=cut

package Bio::EnsEMBL::Funcgen::Hive::ConvertBamToBed;

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

  my $bed_file = $bam_file;
  $bed_file =~ s/.bam$/.bed/;
  
  if (! -e $bed_file) {
    my $cmd = qq(bamToBed -i $bam_file > ${bed_file}.part);
    $self->hive_run_system_cmd($cmd);
    $cmd = qq(mv ${bed_file}.part $bed_file);
    $self->hive_run_system_cmd($cmd);
    push @file_to_delete, $bed_file;
  }
  
  push @file_to_delete, $bed_file;
  
  my $control_result_set = $result_set->get_ControlResultSet;
  my $control_bam_file = $self->db_output_dir . '/' . $control_result_set->dbfile_path;
  
  if (! -e $control_bam_file) {
    die("$control_bam_file doesn't exist!");
  }
  my $control_bed_file = $control_bam_file;
  $control_bed_file =~ s/.bam$/.bed/;

  if (! -e $control_bed_file) {
    my $cmd = qq(bamToBed -i $control_bam_file > ${control_bed_file}.part);
    $self->hive_run_system_cmd($cmd);
    $cmd = qq(mv ${control_bed_file}.part $control_bed_file);
    $self->hive_run_system_cmd($cmd);
    push @file_to_delete, $control_bed_file;
  }
  
  foreach my $current_file (@file_to_delete) {
    $self->dataflow_output_id( {
        file_to_delete => $current_file,
    }, 7);
  }
}

1;
