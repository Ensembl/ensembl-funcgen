=pod

=head1 NAME

    Bio::EnsEMBL::Funcgen::RunnableDB::MotifMapping::DumpBindingMatrixFrequenciesToDisk

=head1 DESCRIPTION

    

=head1 LICENSE

    Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
    Copyright [2016-2021] EMBL-European Bioinformatics Institute

    Licensed under the Apache License, Version 2.0 (the "License"); you may not use this file except in compliance with the License.
    You may obtain a copy of the License at

         http://www.apache.org/licenses/LICENSE-2.0

    Unless required by applicable law or agreed to in writing, software distributed under the License
    is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
    See the License for the specific language governing permissions and limitations under the License.

=head1 CONTACT

    ensembl-funcgen@ebi.ac.uk

=cut

package Bio::EnsEMBL::Funcgen::RunnableDB::MotifMapping::DumpBindingMatrixFrequenciesToDisk;

use strict;
use base ('Bio::EnsEMBL::Hive::Process');

sub run {
    my $self = shift;

    my $species            = $self->param('species');
    my $binding_matrix_dir = $self->param_required('matrices_dir');

    my $funcgen_adaptor =
        Bio::EnsEMBL::Registry->get_DBAdaptor($species, 'funcgen');
    my $binding_matrix_adaptor = $funcgen_adaptor->get_adaptor('BindingMatrix');

    my $binding_matrices = $binding_matrix_adaptor->fetch_all();

    for my $binding_matrix (@{$binding_matrices}) {
        open my $bm_file, '>',
             $binding_matrix_dir . '/' . $binding_matrix->stable_id() . '.pfm';
        print $bm_file $binding_matrix->get_elements_as_string;
        close $bm_file;
    }

}


1;
