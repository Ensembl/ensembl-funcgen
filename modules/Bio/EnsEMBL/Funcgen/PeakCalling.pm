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

use Bio::EnsEMBL::Funcgen::GenericGetSetFunctionality qw(
  _generic_get_or_set
  _generic_set
  _generic_get
  _generic_fetch
);

use Role::Tiny::With;
with 'Bio::EnsEMBL::Funcgen::GenericConstructor';

sub _constructor_parameters {
  return {
    dbID            => 'dbID',
    db              => 'db',
    feature_type_id => 'feature_type_id',
    analysis_id     => 'analysis_id',
    alignment_id    => 'alignment_id',
    name            => 'name',
    display_label   => 'display_label',
  };
}

sub dbID            { return shift->_generic_get_or_set('dbID',            @_); }
sub db              { return shift->_generic_get_or_set('db',              @_); }
sub name            { return shift->_generic_get_or_set('name',            @_); }
sub feature_type_id { return shift->_generic_get_or_set('feature_type_id', @_); }
sub analysis_id     { return shift->_generic_get_or_set('analysis_id',     @_); }
sub alignment_id    { return shift->_generic_get_or_set('alignment_id',    @_); }
sub display_label   { return shift->_generic_get_or_set('display_label',   @_); }

sub fetch_FeatureType { 
  return shift->_generic_fetch('feature_type', 'get_FeatureTypeAdaptor', 'feature_type_id');
}
sub fetch_Analysis {
  return shift->_generic_fetch('analysis', 'get_AnalysisAdaptor', 'analysis_id');
}
sub fetch_Alignment {
  return shift->_generic_fetch('alignment', 'get_AlignmentAdaptor', 'alignment_id');
}

sub _fetch_experiment_id {
    my $self = shift;
    my $db = $self->db;
    my $experiment_id = $db->sql_helper->execute_single_result(
      -SQL    => '
        select 
            distinct read_file_experimental_configuration.experiment_id
        from 
            peak_calling
            join alignment_read_file using (alignment_id)
            join read_file_experimental_configuration using (read_file_id)
            join experiment using (experiment_id)
        where
            experiment.is_control = False
            and peak_calling.peak_calling_id = ?
      ',
      -PARAMS => [ $self->dbID ],
    );
    return $experiment_id;
}

sub fetch_Epigenome {
    my $self = shift;
    my $db = $self->db;
    my $experiment_id = $self->_fetch_experiment_id;
    my $epigenome_id  = $db->sql_helper->execute_single_result(
      -SQL    => '
        select 
            distinct experiment.epigenome_id
        from 
            experiment
        where
            experiment_id = ?
      ',
      -PARAMS => [ $experiment_id ],
    );
    my $epigenome = $db->db->get_EpigenomeAdaptor->fetch_by_dbID($epigenome_id);
    return $epigenome;
}

sub fetch_source_label {
    my $self = shift;
    my $db = $self->db;
    my $experiment_id = $self->_fetch_experiment_id;
    my $source_label  = $db->sql_helper->execute_single_result(
      -SQL    => '
        select 
          experimental_group.name 
        from 
          experiment
          left join experimental_group using (experimental_group_id) 
        where
            experiment_id = ?
      ',
      -PARAMS => [ $experiment_id ],
    );
    return $source_label;
}

1;
