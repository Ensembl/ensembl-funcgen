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
=cut

package MotifFeatureTest;

use strict;
use warnings;

use Test::More;
use Test::Exception;

use parent qw(Bio::EnsEMBL::Funcgen::Test);

sub new {
    my $class = shift;
    my $self = $class->SUPER::new(@_);

    return $self;
}

sub parameters :Test(setup) {
    my $self = shift;

    my $binding_matrix =
        $self->{funcgen_db}->get_adaptor('BindingMatrix')->fetch_by_dbID(403);

    my %mandatory_constructor_parameters = (
        -score          => 100,
        -binding_matrix => $binding_matrix
    );

    $self->{mandatory_constructor_parameters} =
        \%mandatory_constructor_parameters;

    my $slice =
        $self->{core_db}->get_adaptor('Slice')->fetch_by_seq_region_id(131560);

    my %optional_constructor_parameters = (
        -slice     => $slice,
        -stable_id => 'ENSM00000000000',
        -start     => 1,
        -end       => 15,
        -strand    => 1
    );

    my %constructor_parameters = (%mandatory_constructor_parameters,
                                  %optional_constructor_parameters);

    $self->{constructor_parameters} = \%constructor_parameters;

    my $epigenome =
        $self->{funcgen_db}->get_adaptor('Epigenome')->fetch_by_dbID(126);
  
    my $regulatory_feature =
        $self->{funcgen_db}->get_adaptor('RegulatoryFeature')->fetch_by_dbID(15179);  
  
    my %parameters = (
        'epigenome' => $epigenome,
        'regulatory_feature' => $regulatory_feature,
        'position'  => 4
    );

    $self->{parameters} = \%parameters;
}

sub define_expected :Test(setup) {
    my $self = shift;

    my $binding_matrix =
        $self->{funcgen_db}->get_adaptor('BindingMatrix')->fetch_by_dbID(403);

    my $peak_calling =
        $self->{funcgen_db}->get_adaptor('PeakCalling')->fetch_by_dbID(965);
    
    my $epigenomes =
        $self->{funcgen_db}->get_adaptor('Epigenome')
             ->fetch_all_by_dbID_list([ 126, 127, 131, 132, 153 ]);

    $self->{expected} = {
        'binding_matrix'          => $binding_matrix,
        'feature_so_acc'          => 'SO:0000235',
        'feature_so_term'         => 'TF_binding_site',
        'score'                   => -3.16613049708,
        'stable_id'               => 'ENSM00522226970',
        'number_of_peak_callings' => 1,
        'peak_callings'           => [ $peak_calling ],
        'epigenomes'              => $epigenomes,
        'is_position_informative' => 0,
        'is_verified'             => 1,
    };

    my %summary = (
        'binding_matrix_stable_id'              => $binding_matrix->stable_id(),
        'start'                                 => 26327498,
        'end'                                   => 26327520,
        'strand'                                => -1,
        'seq_region_name'                       => 8,
        'stable_id'                             => $self->{expected}->{stable_id},
        'score'                                 => $self->{expected}->{score},
        'transcription_factor_complex'          => '',
        'epigenomes_with_experimental_evidence' => 'HeLa-S3,A549,GM12878,K562,HepG2'
    );

    $self->{expected}->{summary} = \%summary;
}

sub dbIDs_to_fetch {return [5];}

sub getters {
    return [ 'score', 'stable_id' ];
}

sub get_BindingMatrix :Test(1) {
    my $self = shift;

    is_deeply($self->{fetched}->[0]->get_BindingMatrix,
              $self->{expected}->{binding_matrix},
              'get_BindingMatrix works'
    );
}

sub feature_so_acc :Test(1) {
    my $self = shift;

    is($self->{fetched}->[0]->feature_so_acc,
       $self->{expected}->{feature_so_acc},
       'feature_so_acc works');
}

sub feature_so_term :Test(1) {
    my $self = shift;

    is($self->{fetched}->[0] ->feature_so_term,
       $self->{expected}->{feature_so_term},
       'feature_so_term works');
}

sub get_all_overlapping_Peak_Callings_by_Epigenome_and_Regulatory_Feature :Test(2) {
    my $self = shift;

    my $epigenome = $self->{parameters}->{epigenome};
    my $regulatory_feature = $self->{parameters}->{regulatory_feature};

    my @peak_callings = @{$self->{fetched}->[0]->get_overlapping_Peak_Callings_by_Epigenome_and_Regulatory_Feature($epigenome, $regulatory_feature)};

    for my $peak_calling (@peak_callings) {
        isa_ok($peak_calling, 'Bio::EnsEMBL::Funcgen::PeakCalling');
    }

    is(scalar @peak_callings,
       $self->{expected}->{number_of_peak_callings},
       'get_all_overlapping_Peak_Callings_by_Epigenome_and_Regulatory_Feature works');
}

sub get_all_overlapping_Peak_Callings_by_Epigenome :Test(1) {
    my $self        = shift;

    my $epigenome = $self->{parameters}->{epigenome};
    my $peak_calling = $self->{fetched}->[0]->
        get_overlapping_Peak_Callings_by_Epigenome($epigenome);

    is_deeply($peak_calling,
              $self->{expected}->{peak_callings},
              'get_all_overlapping_Peak_Callings_by_Epigenome() works'
    );
}

sub is_experimentally_verified_in_Epigenome :Test(1) {
    my $self = shift;

    my $epigenome = $self->{parameters}->{epigenome};

    my $is_verified = $self->{fetched}->[0]->
        is_experimentally_verified_in_Epigenome($epigenome);

    is ($is_verified,
        $self->{expected}->{is_verified},
        'is_experimentally_verified_in_Epigenome() works'
    );

    # TODO test false
}
  
sub get_all_Epigenomes_with_experimental_evidence :Test(1) {
    my $self        = shift;

    my $epigenomes  = $self->{fetched}->[0]
                           ->get_all_Epigenomes_with_experimental_evidence; 
    is_deeply($epigenomes,
              $self->{expected}->{epigenomes},
              'get_all_Epigenomes_with_experimental_evidence() works');
}
    
sub is_position_informative :Test(1){
    my $self = shift;
    
    my $pos         = $self->{parameters}->{position};
    is($self->{fetched}->[0]->is_position_informative($pos),
       $self->{expected}->{is_position_informative},
       'is_position_informative() works'
    );
}
   
sub infer_variation_consequence :Test(0){
    #TODO write tests
}

1;
