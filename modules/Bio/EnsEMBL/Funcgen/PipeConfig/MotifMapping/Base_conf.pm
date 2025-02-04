=pod

=head1 NAME

    Bio::EnsEMBL::Funcgen::PipeConfig::MotifMapping::Base_conf

=head1 DESCRIPTION

    

=head1 LICENSE

    Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
    Copyright [2016-2025] EMBL-European Bioinformatics Institute

    Licensed under the Apache License, Version 2.0 (the "License"); you may not use this file except in compliance with the License.
    You may obtain a copy of the License at

         http://www.apache.org/licenses/LICENSE-2.0

    Unless required by applicable law or agreed to in writing, software distributed under the License
    is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
    See the License for the specific language governing permissions and limitations under the License.

=head1 CONTACT

    ensembl-funcgen@ebi.ac.uk

=cut

package Bio::EnsEMBL::Funcgen::PipeConfig::MotifMapping::Base_conf;

use strict;
use warnings;

use Bio::EnsEMBL::Hive::PipeConfig::EnsemblGeneric_conf;
use base ('Bio::EnsEMBL::Hive::PipeConfig::EnsemblGeneric_conf');

sub pipeline_wide_parameters {
    my $self = shift;

    my $species = $self->o('species');
    my $assembly = $self->o('assembly');

    $assembly = $self->check_assembly($species, $assembly);

    return {
        %{$self->SUPER::pipeline_wide_parameters},
        species                => $species,
        workdir                => $self->o('workdir'),
        registry               => $self->o('registry'),
        selex_dir              => $self->o('selex_dir'),
        genome                 => $self->o('genome'),
        chr_sizes              => $self->o('chr_sizes'),
        auto_sql               => $self->o('auto_sql'),
        release                => $self->o('release'),
        assembly               => $assembly,
        # bm_2_tf_ep             => {},
        MOODS_dir              => '#workdir#/#species#/MOODS/',
        matrices_dir           => '#MOODS_dir#/binding_matrices/',
        MOODS_raw_dir          => '#MOODS_dir#/raw/',
        MOODS_bed_dir          => '#MOODS_dir#/bed/',
        Peaks_dir              => '#workdir#/#species#/Peaks/',
        RF_dir                 => '#workdir#/#species#/RegulatoryFeatures/',
        Peaks_intersection_dir => '#workdir#/#species#/Peaks_intersection/',
        RF_intersection_dir    => '#workdir#/#species#/RF_intersection/',
        stable_id_mapping_dir  => '#workdir#/#species#/stable_id_mapping',
        temp_dir               => '#workdir#/#species#/temp',
        web_bed_dir            => '#workdir#/#species#/web_bed',
        ftp_dumps_dir          => '#workdir#/#species#/ftp_dumps',
        output_dir             => '#workdir#/#species#/output'
    };
}

sub check_assembly {
    my ($self, $species, $assembly) = @_;

    if ($species eq 'homo_sapiens') {
        if ($assembly =~ /h38/i){
            $assembly = 'GRCh38';
        }
        elsif ($assembly =~ /h37/i){
            $assembly = 'GRCh37';
        }
    }
    elsif ($species eq 'mus_musculus') {
        $assembly = 'GRCm38';
    }

    return $assembly;
}

sub default_options {
    my ($self) = @_;
    return {
        %{$self->SUPER::default_options()},

        workdir       => './',
        pipeline_name => 'MotifMapping',
    };
}

sub resource_classes {
    my ($self) = @_;
    return {
        %{$self->SUPER::resource_classes},

        'default'   => { 'LSF' => '-q production' },
        '250Mb_job' => {
            'LSF' =>
            '-q production -M250   -R"select[mem>250]   rusage[mem=250]"'
        },
        '500Mb_job' => {
            'LSF' =>
            '-q production -M500   -R"select[mem>500]   rusage[mem=500]"'
        },
        '1Gb_job'   => {
            'LSF' =>
            '-q production -M1000  -R"select[mem>1000]  rusage[mem=1000]"'
        },
        '2Gb_job'   => {
            'LSF' =>
            '-q production -M2000  -R"select[mem>2000]  rusage[mem=2000]"'
        },
        '4Gb_job'   => {
            'LSF' =>
            '-q production -M4000  -R"select[mem>4000]  rusage[mem=4000]"'
        },
        '6Gb_job'   => {
            'LSF' =>
            '-q production -M6000  -R"select[mem>6000]  rusage[mem=6000]"'
        },
        '8Gb_job'   => {
            'LSF' =>
            '-q production -M8000  -R"select[mem>8000]  rusage[mem=8000]"'
        },
        '16Gb_job'  => {
            'LSF' =>
            '-q production -M16000 -R"select[mem>16000] rusage[mem=16000]"'
        },
        '24Gb_job'  => {
            'LSF' =>
            '-q production -M24000 -R"select[mem>24000] rusage[mem=24000]"'
        },
        '32Gb_job'  => {
            'LSF' =>
            '-q production -M32000 -R"select[mem>32000] rusage[mem=32000]"'
        },
        '16c32Gb_job'  => {
            'LSF' =>
            '-q production -n 16 -M32000 -R"select[mem>32000] rusage[mem=32000]"'
        },
        '48Gb_job'  => {
            'LSF' =>
            '-q production -M48000 -R"select[mem>48000] rusage[mem=48000]"'
        },
        '64Gb_job'  => {
            'LSF' =>
            '-q production -M64000 -R"select[mem>64000] rusage[mem=64000]"'
        },
    };
}

1;
