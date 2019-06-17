=head1 LICENSE

    Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
    Copyright [2016-2019] EMBL-European Bioinformatics Institute

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

=cut

package MirnaTargetFeatureTest;

use strict;
use warnings;

use Test::More;
use Test::Exception;

use parent qw(Bio::EnsEMBL::Funcgen::Test);

sub parameters :Test(setup) {
    my $self = shift;

    # my $gene_adaptor = $self->{core_db}->get_adaptor('Gene');
    my $feature_type_adaptor = $self->{funcgen_db}->get_adaptor('FeatureType');
    my $analysis_adaptor     = $self->{funcgen_db}->get_adaptor('Analysis');

    my $feature_type =
        $feature_type_adaptor->fetch_by_name('TarBase miRNA target');
    my $analysis = $analysis_adaptor->fetch_by_logic_name('TarBase');

    my %mandatory_constructor_parameters = (
        '-feature_type'           => $feature_type,
        '-analysis'               => $analysis,
        '-gene_stable_id'         => 'ENSG00000000000',
        '-accession'              => 'MIMAT0000',
        '-evidence'               => 'dummy_evidence',
        '-method'                 => 'dummy_method',
        '-supporting_information' => 'dummy_information',
        '-display_label'          => 'dummy_label'
    );

    $self->{mandatory_constructor_parameters} =
        \%mandatory_constructor_parameters;

    my %optional_constructor_parameters = ();

    my %constructor_parameters = (%mandatory_constructor_parameters,
                                  %optional_constructor_parameters);

    $self->{constructor_parameters} = \%constructor_parameters;
}

sub define_expected :Test(setup) {
    my $self = shift;

    my $feature_type_adaptor = $self->{funcgen_db}->get_adaptor('FeatureType');
    my $feature_type = $feature_type_adaptor->fetch_by_dbID(179077);

    $self->{expected} = {
        'display_label'          => 'hsa-let-7a-5p',
        'accession'              => 'MIMAT0000062',
        'evidence'               => 'Experimental',
        'method'                 => 'HITS-CLIP',
        'gene_stable_id'         => 'ENSG00000002834',
        'supporting_information' => 'BS1; spliced_no',
        'feature_type'           => $feature_type
    };
}

sub dbIDs_to_fetch {return [1];}

sub getters {
    return [ 'display_label', 'accession', 'evidence', 'method',
             'gene_stable_id', 'supporting_information' ];
}

sub get_FeatureType :Test(1) {
    my $self = shift;

    is_deeply($self->{fetched}->[0]->get_FeatureType,
              $self->{expected}->{feature_type},
              'get_FeatureType works'
    );
}


1;