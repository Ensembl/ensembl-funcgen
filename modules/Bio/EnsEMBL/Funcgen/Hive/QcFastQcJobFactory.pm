=pod
=head1 NAME

Bio::EnsEMBL::Funcgen::Hive::QcFastQcJobFactory

=head1 DESCRIPTION
=cut

package Bio::EnsEMBL::Funcgen::Hive::QcFastQcJobFactory;

use warnings;
use strict;

use base qw( Bio::EnsEMBL::Funcgen::Hive::BaseDB );

sub run {
  my $self = shift;
  my $input_subset_id = $self->param('input_subset_id');
  my $input_id = $self->create_input_id($input_subset_id);
  
  $self->dataflow_output_id($input_id, 2);
  return;
}

sub create_input_id {

  my $self = shift;
  my $input_subset_id = shift;
  my $work_dir = $self->fastqc_output_dir;
  
  
  my $out_db = $self->param('out_db');
  my $input_subset = $out_db->get_InputSubsetAdaptor->fetch_by_dbID($input_subset_id);
  my $epigenome_production_name = $input_subset->epigenome->production_name;

  my $temp_dir = "$work_dir/$epigenome_production_name/$input_subset_id";

  my $input_id = {
      tempdir               => $temp_dir,
      input_subset_id       => $input_subset_id,
      
      # Connection details for the db to which the results will be written
      tracking_db_user   => $out_db->dbc->user,
      tracking_db_pass   => $out_db->dbc->password,
      tracking_db_host   => $out_db->dbc->host,
      tracking_db_name   => $out_db->dbc->dbname,
      tracking_db_port   => $out_db->dbc->port,
  };
  return $input_id;
}

1;
