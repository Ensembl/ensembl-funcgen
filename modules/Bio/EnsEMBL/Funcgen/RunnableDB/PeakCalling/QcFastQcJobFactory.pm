=pod
=head1 NAME

Bio::EnsEMBL::Funcgen::RunnableDB::PeakCalling::QcFastQcJobFactory

=head1 DESCRIPTION
=cut

package Bio::EnsEMBL::Funcgen::RunnableDB::PeakCalling::QcFastQcJobFactory;

use warnings;
use strict;
use Data::Dumper;

use base 'Bio::EnsEMBL::Hive::Process';

sub run {
  my $self = shift;
  my $species         = $self->param('species');
  my $tempdir         = $self->param('tempdir');
  my $experiment_name = $self->param('experiment');
  
  my $read_file_adaptor = Bio::EnsEMBL::Registry
    ->get_adaptor(
        $species, 
        'funcgen', 
        'ReadFile'
    );
  my $experiment_adaptor = Bio::EnsEMBL::Registry
    ->get_adaptor(
        $species, 
        'funcgen', 
        'Experiment'
    );
  my $experiment = $experiment_adaptor->fetch_by_name($experiment_name);
  if (!defined $experiment) {
    $self->throw("Couldn't find experiment $experiment_name!");
  }
  
  my $all_read_files = $read_file_adaptor->fetch_all_by_Experiment($experiment);

  my $funcgen_dbc = $read_file_adaptor->db->dbc;
  
  foreach my $current_read_file (@$all_read_files) {
  
    my $read_file_name = $current_read_file->name;
  
    my $this_runs_temp_dir = "$tempdir/$read_file_name";
  
    my $input_id = {
        fastqc_tempdir => $this_runs_temp_dir,
        read_file_id   => $current_read_file->dbID,
        read_file_name => $read_file_name,
        read_file      => $current_read_file->file,
        
        species        => $species,
        
        # Connection details for the db to which the results will be written
        tracking_db_user   => $funcgen_dbc->user,
        tracking_db_pass   => $funcgen_dbc->password,
        tracking_db_host   => $funcgen_dbc->host,
        tracking_db_name   => $funcgen_dbc->dbname,
        tracking_db_port   => $funcgen_dbc->port,
    };
    $self->dataflow_output_id($input_id, 2);
  }
  return;
}

1;
