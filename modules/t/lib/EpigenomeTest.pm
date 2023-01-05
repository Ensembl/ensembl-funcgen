=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2023] EMBL-European Bioinformatics Institute

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

package EpigenomeTest;

use strict;
use warnings;

use Test::More;
use Test::Exception;

use parent qw(Bio::EnsEMBL::Funcgen::Test);

sub parameters :Test(setup) {
    my $self = shift;

    my %mandatory_constructor_parameters = (
        -name => 'dummy_epigenome_name'
    );

    $self->{mandatory_constructor_parameters} =
        \%mandatory_constructor_parameters;

    my %optional_constructor_parameters = (
        -short_name      => 'dummy_short_name',
        -description     => 'dummy_description',
        -production_name => 'dummy_production_name',
        -search_terms    => 'dummy_terms',
        -full_name       => 'dummy_full_name'
    );

    my %constructor_parameters = (%mandatory_constructor_parameters,
                                  %optional_constructor_parameters);

    $self->{constructor_parameters} = \%constructor_parameters;

    my %parameters = (
    );

    $self->{parameters} = \%parameters;
}

sub define_expected :Test(setup) {
    my $self = shift;

    my $description = 'Roadmap Epigenomics Mapping Consortium Epigenome'
        . ' (Class 2) for Fetal Intestine Large using'
        . ' donors/samples:H-24595; H-23769';
    my $search_terms = [ "IHECRE00000987", "large intestine",
                         "large intestine male embryo (composite)",
                         "Roadmap Epigenomics Mapping Consortium Epigenome (Class 2) for Fetal Intestine Large using donors/samples:H-24595; H-23769",
                         "ENCSR020SIU", "Roadmap", "intestine",
                         "large intestine", "endoderm",
                         "Entire large intestine", "intestinum crassum",
                         "large intestine", "BiosampleType", "Item",
                         "large intestine",
                         "A subdivision of the digestive tract that connects the small intestine to the cloaca or anus. Lacks or has few villi[Kardong].",
                         "Entire large intestine", "intestinum crassum" ];

    $self->{expected} = {
        'name'             => 'IHECRE00000987',
        'production_name'  => 'large_intestine',
        'gender'           => 'male',
        'description'      => $description,
        'short_name'       => 'large intestine',
        'search_terms'     => $search_terms,
        'full_name'        => 'large intestine male embryo (composite)',
        'efo_accession'    => 'EFO:0000000',
        'encode_accession' => 'ENCSR020SIU',
        'epirr_accession'  => 'IHECRE00000987.1'
    };

    my %summary = (
        'name'             => $self->{expected}->{name},
        'gender'           => $self->{expected}->{gender},
        'description'      => $self->{expected}->{description},
        'short_name'       => $self->{expected}->{short_name},
        'search_terms'     => $search_terms,
        'efo_accession'    => $self->{expected}->{efo_accession},
        'epirr_accession'  => $self->{expected}->{epirr_accession},
        'encode_accession' => $self->{expected}->{encode_accession},
        'full_name'        => $self->{expected}->{full_name}
    );

    $self->{expected}->{summary} = \%summary;
}

sub dbIDs_to_fetch {return [ 151 ];}

sub getters {
    return [ 'name', 'production_name', 'gender', 'description',
             'short_name', 'search_terms', 'full_name',
             'efo_accession', 'encode_accession', 'epirr_accession'];
}

sub reset_relational_attributes :Test(2) {
    my $self = shift;

    my $epigenome = $self->{fetched}->[0];
    $epigenome->reset_relational_attributes;

    is($epigenome->dbID, undef,
       'reset_relational_attributes() resets the dbID');
    is($epigenome->adaptor, undef, '... and also resets the adaptor');
}

1;
