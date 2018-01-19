
=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2018] EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

  Questions may also be sent to the Ensembl help desk at
  <http://www.ensembl.org/Help/Contact>.

=head1 NAME

Bio::EnsEMBL::Funcgen::Hive::RunIDR

=head1 DESCRIPTION

=cut

package Bio::EnsEMBL::Funcgen::Hive::RunIDR;

use warnings;
use strict;

use Bio::EnsEMBL::Utils::Scalar                 qw( assert_ref ); 
use Bio::EnsEMBL::Utils::Exception              qw( throw );
use Bio::EnsEMBL::Funcgen::Sequencing::SeqTools qw( run_IDR );

use base qw( Bio::EnsEMBL::Funcgen::Hive::BaseDB );

sub fetch_input {
  my $self = shift;
  $self->SUPER::fetch_input;
 
  my $result_set_ids = $self->get_param_method('bed_files',  'required');
  assert_ref($result_set_ids, 'ARRAY', 'ResultSet dbIDs');

  $self->get_param_method('idr_threshold', 'required');
  $self->get_param_method('output_prefix', 'required');
  $self->get_param_method('batch_name',    'required');
  $self->get_output_work_dir_methods;
  return;
}

sub run {
  my $self = shift;

  my $num_peaks;
  my $png_file;
  
  eval {
    ($num_peaks, $png_file) = run_IDR(
      -out_dir       => $self->output_dir, 
      -output_prefix => $self->output_prefix, 
      -threshold     => $self->idr_threshold, 
      -bed_files     => $self->bed_files,
      -batch_name    => $self->batch_name, 
    );
  };

  if($@) {
    $self->throw_no_retry("Failed to run_IDR ".$self->output_prefix."\n$@");
  }
  
  warn("Num peaks is $num_peaks png is $png_file");
  
  $self->set_param_method('num_peaks', $num_peaks);
  return;
}

sub write_output {
  my $self = shift;
  $self->dataflow_output_id( {'idr_peak_counts'   => $self->num_peaks}, 2);
}

1;