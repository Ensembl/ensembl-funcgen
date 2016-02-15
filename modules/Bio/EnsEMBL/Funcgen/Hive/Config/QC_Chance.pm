=pod 

=head1 NAME

    Bio::EnsEMBL::Funcgen::Hive::Config::QC_Chance

=head1 SYNOPSIS

=head1 DESCRIPTION

=head1 LICENSE

    Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute

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

use base ('Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf');  # All Hive databases configuration files should inherit from HiveGeneric, directly or indirectly


sub pipeline_wide_parameters {
    my ($self) = @_;
    return {
        %{$self->SUPER::pipeline_wide_parameters},          # here we inherit anything from the base class
    };
}

sub pipeline_analyses {
    my ($self) = @_;
    return [
#         {   -logic_name => 'ArgenrichJobDefinition',
#             -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
#             -meadow_type=> 'LOCAL',
#             -input_ids => [
#             {
# 		
#                 column_names => [ 'kind', 'file', 'sourcedir', 'tempdir' ],
# 		inputlist    => [
# 		  [ 'signal',  'F36P:hist:BR2_H3K27me3_3526_bwa_samse_1_2_3.bam', "#sourcedir#", "#tempdir#" ],
# 		  [ 'control', 'F36P:hist:BR2_WCE_3526_bwa_samse_1.bam',          "#sourcedir#", "#tempdir#" ],
# 		],
# 		# Directory in which the bam files are
# 		sourcedir             => '/warehouse/ensembl10/funcgen/alignments/homo_sapiens/GRCh38/3526',
# 		
# 		# Directory into which the bam files will be copied
# 		tempdir               => '/lustre/scratch109/ensembl/funcgen/mn1/ersa/debug/F36P:hist:BR2_H3K27me3_3526',
# 		
# 		# Name of the output file that argenrich creates. This will be in #tempdir#
# 		argenrich_outfile     => 'argenrich_outfile.txt',
# 		
# 		# result_set_id of the control to which this will be linked
# 		control_result_set_id => 1,
# 		
# 		# result_set_id of the signal to which this will be linked
# 		signal_result_set_id  => 2,
# 		
# 		# The file with chromosome lengths.
# 		chrlenfile            => '/lustre/scratch109/ensembl/funcgen/mn1/ersa/faang/reference_files/CCAT/homo_sapiens_.CCAT_chr_lengths.txt',
# 		
# 		# #chrlenfile# needs to be sorted first.  This is done in the 
# 		# SortChrLenFile analysis. #chrlenfilesorted# is the name of 
# 		# the file to which the sorted file is written.
# 		#
# 		chrlenfilesorted      => '/lustre/scratch109/ensembl/funcgen/mn1/ersa/faang/reference_files/CCAT/homo_sapiens_.CCAT_chr_lengths.chrlenfilesorted.txt',
# 		
# 		# Connection details for the db to which the results will be written
# 		tracking_db_user   => 'ensadmin',
# 		tracking_db_pass   => 'xxx',
# 		tracking_db_host   => 'ens-genomics2',
# 		tracking_db_name   => 'mn1_faang_tracking_homo_sapiens_funcgen_81_38',
# 		
#             }
#             ],
#             -flow_into => { 1 => 'MkTempDir', },
#         },
	{
	    -logic_name => 'PreprocessAlignments',
	    -module     => 'Bio::EnsEMBL::Funcgen::Hive::CollectionWriter',
	    -flow_into => {
		'3'   => [  'QcChanceJobFactory' ],
		'4'   => [  'QcChanceJobFactory' ],
		'5'   => [  'QcChanceJobFactory' ],
		'6'   => [  'QcChanceJobFactory' ],
		'100' => [  'QcChanceJobFactory' ],
	    },
	},
	{
	  -logic_name    => 'QcChanceJobFactory',
	  -module        => 'Bio::EnsEMBL::Funcgen::Hive::QcChanceJobFactory',
	  -meadow        => 'LOCAL',
	  -parameters => {
		#chromosome_file => $self->o('data_root_dir'). '/reference_files/CCAT/'.$self->o('species').'_'.$self->o('assembly').'.CCAT_chr_lengths.txt'
		chromosome_file => $self->o('chromosome_file')
	  },
	  -flow_into => { 2 => 'MkTempDir', },
	},
        {   -logic_name => 'MkTempDir',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -meadow_type=> 'LOCAL',
            -parameters => {
		  cmd => qq!mkdir -p #tempdir#!,
            },
            -flow_into => { 
	      '1->A' => [ 'JobFactoryArgenrich', 'SortChrLenFile' ],
	      'A->1' => [ 'RunArgenrich' ],
            },
        },
        {   -logic_name => 'SortChrLenFile',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -meadow_type=> 'LOCAL',
            -parameters => { 
		  cmd => qq!sort -k1,1 #chrlenfile# > #chrlenfilesorted#!,
            },
            -flow_into => { 1 => 'argenrichformregions', },
        },
        {   -logic_name => 'argenrichformregions',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -meadow_type=> 'LOCAL',
            -parameters => { 
		  cmd => qq!/software/ensembl/funcgen/argenrichformregions.pl #chrlenfilesorted#!,
            },
        },
        {   -logic_name => 'JobFactoryArgenrich',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::JobFactory',
            -meadow_type=> 'LOCAL',
            -flow_into => {
		# Skipping copy, we can work on the files directly.
                2 => [ 'CpToTemp' ],
                #2 => [ 'IndexBam' ],
            },
        },
        {   -logic_name => 'CpToTemp',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -meadow_type=> 'LOCAL',
            -parameters => { 
		  #cmd => qq!cp #sourcedir#/#file# #tempdir#!,
		  cmd => qq!ln -s #sourcedir#/#file# #tempdir#!,
            },
            -flow_into => { 1 => 'IndexBam' },
        },
        {   -logic_name => 'IndexBam',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -meadow_type=> 'LSF',
            -parameters => { 
		  cmd => qq!samtools index #tempdir#/#file#!,
            },
            -flow_into => { 1 => 'CountReads' },
        },
        {   -logic_name => 'CountReads',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::JobFactory',
            -meadow_type=> 'LSF',
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
            -meadow_type=> 'LSF',
	    -parameters => {
                'cmd' => qq(/nfs/users/nfs_m/mn1/work_dir_faang/argenrich.R --args plot=TRUE outdir=#tempdir# )
		  . qq(    ipsz=#expr( #read_count#->{"signal"}            )expr# )
		  . qq( inputsz=#expr( #read_count#->{"control"}           )expr# )
		  . qq(      ip=#tempdir#/#expr( #file#->{"signal"}        )expr# )
		  . qq(   input=#tempdir#/#expr( #file#->{"control"}       )expr# )
		  . qq( outfile=#argenrich_outfile#)
            },
            -flow_into => { 1 => 'LoadChanceToDB', },
            -rc_name   => 'normal_monitored_16GB',
        },
        {   -logic_name => 'LoadChanceToDB',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -meadow_type=> 'LOCAL',
	    -parameters => {
                'cmd' => qq(load_argenrich_qc_file.pl        )
		  . qq( --argenrich_file        #tempdir#/#argenrich_outfile#     )
		  #. qq( --control_result_set_id #control_result_set_id#           )
		  . qq( --signal_result_set_id  #signal_result_set_id#            )
		  . qq( --user #tracking_db_user# --pass #tracking_db_pass# --host #tracking_db_host# --dbname #tracking_db_name# )
            },
        },
    ];
}

1;



