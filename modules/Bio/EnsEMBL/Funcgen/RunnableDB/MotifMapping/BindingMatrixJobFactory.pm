=pod

=head1 NAME

    Bio::EnsEMBL::Funcgen::RunnableDB::MotifMapping::BindingMatrixJobFactory

=head1 DESCRIPTION



=head1 LICENSE

    Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
    Copyright [2016-2023] EMBL-European Bioinformatics Institute

    Licensed under the Apache License, Version 2.0 (the "License"); you may not use this file except in compliance with the License.
    You may obtain a copy of the License at

         http://www.apache.org/licenses/LICENSE-2.0

    Unless required by applicable law or agreed to in writing, software distributed under the License
    is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
    See the License for the specific language governing permissions and limitations under the License.

=head1 CONTACT

    ensembl-funcgen@ebi.ac.uk

=cut

package Bio::EnsEMBL::Funcgen::RunnableDB::MotifMapping::BindingMatrixJobFactory;

use base 'Bio::EnsEMBL::Hive::RunnableDB::JobFactory';

# sub param_defaults {
#     my $self = shift @_;
#
#     # $self->fetch_input();
#     # my $species = $self->param_required('species');
#     # my $species = $self->o('species');
#     # use Data::Printer; p $self->{'_param_hash'};
#
#     # my @binding_matrix_stable_ids = keys %{$self->param_required('bm_2_tf_ep')};
#     my $param_defaults = $self->SUPER::param_defaults; use Data::Printer;
#     p $param_defaults;
#     $param_defaults->{inputlist} = [1,2,3];
#     $param_defaults->{column_names} = ['matrix'];
#
#     return $param_defaults;
# }

sub run {
    my $self = shift;
    # use Data::Printer;
    # p $self;
    $self->param('column_names', ['matrix']);
    my $species = $self->param_required('species');
    my $param = $self->_get_mapping($species);
    use Data::Printer; p $param;
    $self->db->hive_pipeline->add_new_or_update('PipelineWideParameters',
                                                'param_name'  => 'bm_2_tf_ep',
                                                'param_value' => $param);
    p $self;
    p $self->db;
    p $self->db->hive_pipeline;

    my @binding_matrix_stable_ids = keys %{$self->param_required('bm_2_tf_ep')};
    if($binding_matrix_stable_ids == 0){die;}
    $self->param('inputlist', \@binding_matrix_stable_ids);
    p @binding_matrix_stable_ids;
    # use Data::Printer; p $self->{_input_job};

    $self->SUPER::run();

}

sub _get_mapping {
    my $species = shift;

    my $bm_2_tf_ep = {};

    my $funcgen_adaptor
                               = Bio::EnsEMBL::Registry->get_DBAdaptor($species, 'funcgen');
    my $binding_matrix_adaptor = $funcgen_adaptor->get_adaptor('BindingMatrix');
    my $all_binding_matrices = $binding_matrix_adaptor->fetch_all();

    for my $binding_matrix (@{$all_binding_matrices}) {
        my $bm_stable_id = $binding_matrix->stable_id();
        my $peak_callings = $binding_matrix->get_all_PeakCallings();
        my @pairs = ();

        for my $peak_calling (@{$peak_callings}){
            my $tf_name = uc($peak_calling->get_FeatureType->name());
            my $epigenome_name =
                $peak_calling->get_Epigenome()->production_name();
            push @pairs, [$tf_name, $epigenome_name];
        }

        $bm_2_tf_ep->{$bm_stable_id} = \@pairs;
    }

    return $bm_2_tf_ep;
}



1;