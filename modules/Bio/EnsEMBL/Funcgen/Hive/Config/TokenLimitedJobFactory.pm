=pod 

=head1 NAME

    Bio::EnsEMBL::Funcgen::Hive::Config::RollingJobFactoryDemo2

=head1 SYNOPSIS

=head1 DESCRIPTION

=head1 LICENSE

    Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
    Copyright [2016-2018] EMBL-European Bioinformatics Institute

    Licensed under the Apache License, Version 2.0 (the "License"); you may not use this file except in compliance with the License.
    You may obtain a copy of the License at

         http://www.apache.org/licenses/LICENSE-2.0

    Unless required by applicable law or agreed to in writing, software distributed under the License
    is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
    See the License for the specific language governing permissions and limitations under the License.

=head1 CONTACT

=cut
package Bio::EnsEMBL::Funcgen::Hive::Config::TokenLimitedJobFactory;

use strict;
use warnings;

use base ('Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf');

sub pipeline_analyses {
    my ($self) = @_;
    return [
        {   -logic_name => 'JobPool',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
            -input_ids => [ 
	      { id => 1, },
	      { id => 2, },
	      { id => 3, },
	      { id => 4, },
	      { id => 5, },
            ],
            -wait_for    => 'JobPool',
            -meadow_type => 'LOCAL',
            -flow_into => {
	      '1' => [ 'TokenLimitedJobFactory' ],
            },
        },
        {   -logic_name => 'TokenLimitedJobFactory',
            -module     => 'Bio::EnsEMBL::Funcgen::Hive::TokenLimitedJobFactory',
            -meadow_type=> 'LOCAL',
            -input_ids => [ 
	      { x => 1 },
	      { x => 2 },
            ],
            -flow_into => {
	      '2->A' => [ 'DoPipelineWork' ],
	      'A->1' => [ 'TokenLimitedJobFactory' ],
            },
        },
        {   -logic_name => 'DoPipelineWork',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
            -meadow_type=> 'LOCAL',
            -flow_into => {
	      '1' => [ 'EndOfPipelineWork' ],
            },
        },
        {   -logic_name => 'EndOfPipelineWork',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
            -meadow_type=> 'LOCAL',
        },
    ];
}


1;




