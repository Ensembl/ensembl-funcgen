#
# Ensembl module for Bio::EnsEMBL::DBSQL::Funcgen::AnnotatedFeatureAdaptor
#

=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2025] EMBL-European Bioinformatics Institute

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

package Bio::EnsEMBL::Funcgen::DBSQL::EpigenomeTrackAdaptor;

use strict;
use base 'Bio::EnsEMBL::Funcgen::DBSQL::GenericAdaptor';
use Bio::EnsEMBL::Utils::Exception qw( throw );

sub object_class {
    return 'Bio::EnsEMBL::Funcgen::EpigenomeTrack';
}

sub _tables {
  return ['epigenome_track', 'et']
}

sub fetch_all_by_Epigenome {
    my $self      = shift;
    my $epigenome = shift;
    
    if (! defined $epigenome) {
      throw("Epigenome was undefined");
    }
    my $constraint = 'epigenome_id  = ' . $epigenome->dbID;
    return $self->fetch_all($constraint);
}

sub fetch_all_by_FeatureType {
    my $self      = shift;
    my $feature_type = shift;

    if (! defined $feature_type) {
      throw("FeatureType was undefined");
    }
    my $constraint = 'feature_type_id  = ' . $feature_type->dbID;
    return $self->fetch_all($constraint);
}

1;
