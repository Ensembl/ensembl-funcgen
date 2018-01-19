=pod 

=head1 NAME

    Bio::EnsEMBL::Funcgen::Hive::Config::QC_Fastqc

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

package Bio::EnsEMBL::Funcgen::Hive::Config::QC_Fastqc;

use strict;
use warnings;

use Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf;
use base ('Bio::EnsEMBL::Funcgen::Hive::Config::Base');

sub pipeline_analyses {
    my ($self) = @_;
    return [
	{
	  -logic_name => 'create_job_batches',
	  -flow_into => {
	    2 => 'fastqc_start',
	  },
	},
        {   -logic_name => 'fastqc_start',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
            -flow_into => { 
              'MAIN->A' => 'QcFastQcInputIdsFromInputSet',
              'A->MAIN' => 'fastqc_done',
            },
        },
        {   -logic_name => 'QcFastQcInputIdsFromInputSet',
            -module     => 'Bio::EnsEMBL::Funcgen::Hive::QcFastQcInputIdsFromInputSet',
            -flow_into => { 
	      2 => 'QcFastQcJobFactory',
            },
        },
        {   -logic_name => 'QcFastQcJobFactory',
            -module     => 'Bio::EnsEMBL::Funcgen::Hive::QcFastQcJobFactory',
            -flow_into => { 
	      2 => 'MkFastQcTempDir',
            },
        },
        {   -logic_name => 'MkFastQcTempDir',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -parameters => { 
		  cmd => qq!mkdir -p #tempdir#!,
            },
            -flow_into => { 
	      MAIN => 'JobFactoryFastQC',
            },
        },
        {   -logic_name => 'JobFactoryFastQC',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::JobFactory',
            -parameters => { 
		  inputquery => qq(select local_url, "#tempdir#" as tempdir from input_subset_tracking where input_subset_id = #input_subset_id#),
		  db_conn    => 'mysql://#tracking_db_user#:#tracking_db_pass#@#tracking_db_host#:#tracking_db_port#/#tracking_db_name#'
# 		  db_conn    => 'funcgen:#species#'
            },
            -flow_into => {
                '2->A' => 'RunFastQC',
                'A->1' => 'QcFastQcLoaderJobFactory',
            },
        },
        {   -logic_name => 'RunFastQC',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -parameters => { 
		  cmd => qq(fastqc --extract -o #tempdir# #local_url#),
            },
            -rc_name => 'normal_2GB',
            # For some reason this analysis is very attractive to workers. 
            # So to prevent them from focusing on this at the beginning of the
            # pipeline, I'm limiting this.
            #
            -analysis_capacity => 5,
        },
        {   -logic_name => 'QcFastQcLoaderJobFactory',
            -module     => 'Bio::EnsEMBL::Funcgen::Hive::QcFastQcLoaderJobFactory',
            -flow_into => {
                2 => 'QcLoadFastQcResults',
            },
        },
        {   -logic_name        => 'QcLoadFastQcResults',
            -module            => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -parameters => {
                  use_bash_pipefail => 1,
		  cmd => qq(load_fastqc_summary_file.pl        )
		    . qq( --input_subset_id #input_subset_id#  )
		    . qq( --summary_file #fastqc_summary_file# )
		    . qq( --work_dir #tempdir#                 )
		    . qq( | mysql )
		    . qq( --host #tracking_db_host#  )
		    . qq( --port #tracking_db_port#  )
		    . qq( --user #tracking_db_user#  )
		    . qq( -p#tracking_db_pass#       )
		    . qq( #tracking_db_name#         ),

            },
        },
        {   -logic_name => 'fastqc_done',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
        },
    ];
}

1;
