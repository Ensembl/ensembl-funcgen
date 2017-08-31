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

package Bio::EnsEMBL::Funcgen::PeakCalling;

use strict;
use warnings;

use base 'Bio::EnsEMBL::Funcgen::GenericGetSetFunctionality';

sub _constructor_parameters {
  return {
    feature_type_id => 'feature_type_id',
    analysis_id     => 'analysis_id',
    alignment_id    => 'alignment_id',
  };
}

sub _simple_accessors {
  return [
    { method_name => 'feature_type_id', hash_key    => '_feature_type_id', },
    { method_name => 'analysis_id',     hash_key    => '_analysis_id',     },
    { method_name => 'alignment_id',    hash_key    => '_alignment_id',    },
  ]
}

sub _fetch_methods {
  return [
    {
      method_name             => 'fetch_FeatureType',
      hash_key                => '_feature_type',
      get_adaptor_method_name => 'get_FeatureTypeAdaptor',
      dbID_method             => 'feature_type_id',
    },
    {
      method_name             => 'fetch_Analysis',
      hash_key                => '_analysis',
      get_adaptor_method_name => 'get_AnalysisAdaptor',
      dbID_method             => 'analysis_id',
    },
    {
      method_name             => 'fetch_Alignment',
      hash_key                => '_alignment',
      get_adaptor_method_name => 'get_AlignmentAdaptor',
      dbID_method             => 'alignment_id',
    },
  ]
}

1;
