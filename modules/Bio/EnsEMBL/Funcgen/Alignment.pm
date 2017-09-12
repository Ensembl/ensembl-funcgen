=head1 LICENSE

Copyright [1999-2016] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute

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

=head1 SYNOPSIS

=head1 DESCRIPTION

=cut

package Bio::EnsEMBL::Funcgen::Alignment;

use strict;
use warnings;

use base 'Bio::EnsEMBL::Funcgen::GenericGetSetFunctionality';

sub _constructor_parameters {
  return {
    name           => 'name',
    analysis_id    => 'analysis_id',
    bam_file_id    => 'bam_file_id',
    bigwig_file_id => 'bigwig_file_id',
    read_file_ids  => 'read_file_ids',
  };
}

sub _simple_accessors {
  return [
    { method_name => 'name',              hash_key    => '_name',           },
    { method_name => 'analysis_id',       hash_key    => '_analysis_id',    },
    { method_name => 'bam_file_id',       hash_key    => '_bam_file_id',    },
    { method_name => 'bigwig_file_id',    hash_key    => '_bigwig_file_id', },
    { method_name => 'read_file_ids',     hash_key    => '_read_file_ids',  },
  ]
}

sub _fetch_methods {
  return [
    {
      method_name             => 'fetch_Analysis',
      hash_key                => '_analysis',
      get_adaptor_method_name => 'get_AnalysisAdaptor',
      dbID_method             => 'analysis_id',
    },
  ]
}

sub _set_methods {
  return [
    {
      method_name   => 'set_Analysis',
      expected_type => 'Bio::EnsEMBL::Analysis',
      hash_key      => '_analysis',
    },
  ]
}

sub _fetch_DataFile {

  my $self         = shift;
  my $data_file_id = shift;
  
  my $data_file_adaptor = $self->db->db->get_DataFileAdaptor;
  if (! defined $data_file_adaptor) {
    throw("Couldn't get a DataFileAdaptor!");
  }
  my $data_file = $data_file_adaptor->fetch_by_dbID($data_file_id);
  return $data_file;
}

sub has_bam_DataFile {

  my $self = shift;
  return defined $self->bam_file_id;
}

sub has_bigwig_DataFile {

  my $self = shift;
  return defined $self->bigwig_file_id;
}

sub fetch_bam_DataFile {

  my $self        = shift;
  my $bam_file_id = $self->bam_file_id;
  
  return $self->_fetch_DataFile($bam_file_id);
}

sub fetch_bigwig_DataFile {

  my $self           = shift;
  my $bigwig_file_id = $self->bigwig_file_id;
  
  return $self->_fetch_DataFile($bigwig_file_id);
}

sub fetch_BamFile {
  my $self = shift;
  return $self->fetch_bam_DataFile;
}

sub fetch_BigWigFile {
  my $self           = shift;
  return $self->fetch_bigwig_DataFile;
}

sub fetch_all_ReadFiles {

  my $self = shift;
  
  my $read_file_id = $self->read_file_ids;
  my $read_file_adaptor = $self->db->db->get_ReadFileAdaptor;
  
  if (! defined $read_file_adaptor) {
    throw("Couldn't get a ReadFileAdaptor!");
  }
  
  my @all_read_files;
  
  foreach my $current_read_file_id (@$read_file_id) {
    my $current_read_file = $read_file_adaptor->fetch_by_dbID($current_read_file_id);
    push @all_read_files, $current_read_file;
  }
  return \@all_read_files;
}

1;

