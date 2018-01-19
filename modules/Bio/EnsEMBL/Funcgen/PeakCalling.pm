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

Bio::EnsEMBL::Funcgen::PeakCalling - Object that represents a peak calling

=head1 SYNOPSIS

  use strict;
  use warnings;
  use Bio::EnsEMBL::Registry;
  use List::Util qw( min );
  use Data::Dumper;

  my $registry = 'Bio::EnsEMBL::Registry';

  $registry->load_registry_from_db(
      -host => 'ensembldb.ensembl.org',
      -user => 'anonymous'
  );

  my $peak_calling_adaptor = Bio::EnsEMBL::Registry->get_adaptor('homo_sapiens', 'funcgen', 'PeakCalling');

  my $all_peak_callings = $peak_calling_adaptor->fetch_all;

  my $number_of_peak_callings_available = @$all_peak_callings;

  print "There are $number_of_peak_callings_available peak callings available for querying:\n";

  # Print the first ten
  my $max_features_to_print = 10;

  for my $i ( 1.. min($max_features_to_print, $number_of_peak_callings_available) ) {

    my $current_peak_calling = $all_peak_callings->[$i];
    print "  - " . $current_peak_calling->display_label . "\n";
    
  }

=head1 DESCRIPTION

This object represents a peak calling from a ChIP-seq or other high-throughput
assay. It links to 

  - the set of Peaks that were generated and 
  - the alignment that the peak calling was done on and
  - the peak caller that was used via the analysis.

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
    experiment_id   => '_experiment_id',
  };
}

sub dbID            { return shift->_generic_get_or_set('dbID',            @_); }
sub db              { return shift->_generic_get_or_set('db',              @_); }

=head2 name

  Example    : my $name = $peak_calling->name;
  Description: Accessor for the name of the peak calling.
  Returntype : String
  Exceptions : None
  Caller     : general
  Status     : Stable

=cut
sub name            { return shift->_generic_get_or_set('name',            @_); }
sub feature_type_id { return shift->_generic_get_or_set('feature_type_id', @_); }
sub analysis_id     { return shift->_generic_get_or_set('analysis_id',     @_); }
sub alignment_id    { return shift->_generic_get_or_set('alignment_id',    @_); }
sub _experiment_id  { return shift->_generic_get_or_set('_experiment_id',  @_); }

=head2 display_label

  Example    : my $display_label = $peak_calling->display_label;
  Description: Accessor for the display_label of the peak calling. This is 
               used as the name displayed in the web browser.
  Returntype : String
  Exceptions : None
  Caller     : general
  Status     : Stable

=cut
sub display_label   { return shift->_generic_get_or_set('display_label',   @_); }

=head2 fetch_FeatureType

  Example    : my $feature_type = $peak_calling->fetch_FeatureType;
  Description: Fetches the feature type of the peak calling. This is the 
               type of feature the experiment was assaying for.
  Returntype : Bio::EnsEMBL::Funcgen::FeatureType
  Exceptions : None
  Caller     : general
  Status     : Stable

=cut
sub fetch_FeatureType { 
  return shift->_generic_fetch('feature_type', 'get_FeatureTypeAdaptor', 'feature_type_id');
}

=head2 fetch_Analysis

  Example    : my $analysis = $peak_calling->fetch_Analysis;
  Description: Fetches the analysis of the peak calling. This is the analysis
               representing the peak caller that was used to analyse the 
               alignment.
  Returntype : Bio::EnsEMBL::Analysis
  Exceptions : None
  Caller     : general
  Status     : Stable

=cut
sub fetch_Analysis {
  return shift->_generic_fetch('analysis', 'get_AnalysisAdaptor', 'analysis_id');
}

=head2 fetch_Alignment

  Example    : my $alignment = $peak_calling->fetch_Alignment;
  Description: Fetches the alignment on which the peak calling was done.
               alignment.
  Returntype : Bio::EnsEMBL::Funcgen::Alignment
  Exceptions : None
  Caller     : general
  Status     : Stable

=cut
sub fetch_Alignment {
  return shift->_generic_fetch('alignment', 'get_AlignmentAdaptor', 'alignment_id');
}

=head2 fetch_Epigenome

  Example    : my $epigenome = $peak_calling->fetch_Epigenome;
  Description: Fetches the epigenome that was used in the assay.
  Returntype : Bio::EnsEMBL::Funcgen::Epigenome
  Exceptions : None
  Caller     : general
  Status     : Stable

=cut
sub fetch_Epigenome {
    my $self = shift;
    my $db = $self->db;
    my $experiment_id = $self->_experiment_id;
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

sub fetch_Experiment {
  return shift->_generic_fetch('experiment', 'get_ExperimentAdaptor', '_experiment_id');
}

=head2 fetch_source_label

  Example    : my $source_label = $peak_calling->fetch_source_label;
  Description: Fetches the name of the experimental group of the experiment 
               that led to this peak calling.
  Returntype : String
  Exceptions : None
  Caller     : general
  Status     : Stable

=cut
sub fetch_source_label {
    my $self = shift;
    my $db = $self->db;
    my $experiment_id = $self->_experiment_id;
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
