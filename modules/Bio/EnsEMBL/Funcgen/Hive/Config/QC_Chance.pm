=pod 

=head1 NAME

    Bio::EnsEMBL::Funcgen::Hive::Config::QC_Chance

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

package Bio::EnsEMBL::Funcgen::Hive::Config::QC_Chance;

use strict;
use warnings;

use Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf;
use base ('Bio::EnsEMBL::Funcgen::Hive::Config::Base');

sub pipeline_analyses {
    my ($self) = @_;
    return [
	{
	    -logic_name => 'index_bam_files',
	    -module     => 'Bio::EnsEMBL::Funcgen::Hive::CollectionWriter',
	    -flow_into => {
		2  => 'qc_chance_start',
	    },
	},
        {   -logic_name => 'qc_chance_start',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
            -flow_into => { 
              'MAIN->A' => 'QcChanceJobFactory',
              'A->MAIN' => 'qc_chance_done',
            },
        },
	{
	  -logic_name    => 'QcChanceJobFactory',
	  -module        => 'Bio::EnsEMBL::Funcgen::Hive::QcChanceJobFactory',
# 	  -parameters => {
# 		chromosome_file => $self->o('chromosome_file')
# 	  },
	  -flow_into => { 2 => 'MkTempDir', },
	},
        {   -logic_name => 'MkTempDir',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -parameters => {
		  cmd => qq!mkdir -p #tempdir#!,
            },
            -flow_into => { 
	      '1->A' => [ 'JobFactoryArgenrich', 'CreateChanceBins'],
	      'A->1' => 'RunArgenrich',
            },
        },
	{   -logic_name => 'CreateChanceBins',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
	    -parameters => {
                cmd => qq(create_chance_bins.pl )
		  . qq(  --species                 #species#                     )
		  . qq(  --epigenome_gender        #epigenome_gender#            )
		  . qq(  --assembly                #assembly#                    )
		  . qq(  --outputfile              #tempdir#/#chance_bin_file#   )
  		  . qq(  --reference_data_root_dir #reference_data_root_dir#     )
            },
        },
        {   -logic_name => 'JobFactoryArgenrich',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::JobFactory',
            -flow_into => {
                2 => 'CpToTemp',
            },
        },
        {   -logic_name => 'CpToTemp',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -parameters => { 
                  cmd => qq!rm -f #tempdir#/#file# ; ln -s #sourcedir#/#file# #tempdir#!,
            },
            -flow_into => { MAIN => 'IndexBam' },
        },
        {   -logic_name => 'IndexBam',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -parameters => { 
		  cmd => qq!samtools index #tempdir#/#file#!,
            },
            -flow_into => { MAIN => 'CountReads' },
        },
        {   -logic_name => 'CountReads',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::JobFactory',
            -parameters => { 
		  inputcmd        => "samtools view -c #tempdir#/#file#",
		  column_names    => [ 'read_count' ],
            },
            -flow_into => {
                2 => [
		  ':////accu?read_count={kind}',
		  ':////accu?file={kind}',
                ],
             },
        },
        {   -logic_name => 'RunArgenrich',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
	    -parameters => {
                cmd => qq(argenrich_with_labels_and_rerunnable.R --args plot=TRUE outdir=#tempdir# )
		  . qq(    ipsz=#expr( #read_count#->{"signal"}            )expr# )
		  . qq( inputsz=#expr( #read_count#->{"control"}           )expr# )
		  . qq(      ip=#tempdir#/#expr( #file#->{"signal"}        )expr# )
		  . qq(   input=#tempdir#/#expr( #file#->{"control"}       )expr# )
		  . qq(    bins=#tempdir#/#chance_bin_file#                       )
		  
		  # This ends up in #tempdir#, because of the parameter 
		  # "outdir=#tempdir#" set further above.
		  #
		  . qq( outfile=#argenrich_outfile#)
            },
            -flow_into => { MAIN => 'LoadChanceToDB', },
            -rc_name   => 'normal_monitored_16GB',
        },
        {   -logic_name => 'LoadChanceToDB',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
	    -parameters => {
                cmd => qq(load_argenrich_qc_file.pl   )
		  . qq( --argenrich_file        #tempdir#/#argenrich_outfile#     )
		  . qq( --signal_result_set_id  #signal_result_set_id#            )
		  . qq( --user   #tracking_db_user#   )
                  . qq( --pass   #tracking_db_pass#   )
                  . qq( --port   #tracking_db_port#   )
		  . qq( --host   #tracking_db_host#   )
		  . qq( --dbname #tracking_db_name#   )
		  . qq( --work_dir #tempdir#  )
            },
        },
        {   -logic_name => 'qc_chance_done',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
        },
    ];
}

1;
