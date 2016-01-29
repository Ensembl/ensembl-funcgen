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
        {   -logic_name => 'ArgenrichJobDefinition',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
            -meadow_type=> 'LOCAL',
            -input_ids => [
            {
                column_names => [ 'kind', 'file', 'sourcedir', 'tempdir' ],
		inputlist    => [
		  [ 'signal',  'F36P:hist:BR2_H3K27me3_3526_bwa_samse_1_2_3.bam', "#sourcedir#", "#tempdir#" ],
		  [ 'control', 'F36P:hist:BR2_WCE_3526_bwa_samse_1.bam',          "#sourcedir#", "#tempdir#" ],
		],
		sourcedir             => '/warehouse/ensembl10/funcgen/alignments/homo_sapiens/GRCh38/3526',
		tempdir               => '/lustre/scratch109/ensembl/funcgen/mn1/ersa/debug/F36P:hist:BR2_H3K27me3_3526',
		argenrich_outfile     => 'argenrich_outfile.txt',
		control_result_set_id => 1,
		signal_result_set_id  => 2,
		chrlenfile            => '/lustre/scratch109/ensembl/funcgen/mn1/ersa/faang/reference_files/CCAT/homo_sapiens_.CCAT_chr_lengths.txt',
		chrlenfilesorted      => '/lustre/scratch109/ensembl/funcgen/mn1/ersa/faang/reference_files/CCAT/homo_sapiens_.CCAT_chr_lengths.chrlenfilesorted.txt',
		
		tracking_db_user   => 'ensadmin',
		tracking_db_pass   => 'xxx',
		tracking_db_host   => 'ens-genomics1',
		tracking_db_name   => 'mn1_faang_tracking_homo_sapiens_funcgen_81_38',
		
            }
            ],
            -flow_into => { 1 => 'MkTempDir', },
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
                2 => [ 'CpToTemp' ],
            },
        },
        {   -logic_name => 'CpToTemp',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -meadow_type=> 'LOCAL',
            -parameters => { 
		  cmd => qq!cp #sourcedir#/#file# #tempdir#!,
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
            -flow_into => { 1 => 'LoadToDB', },
            -rc_name   => 'normal_monitored_16GB',
        },
        {   -logic_name => 'LoadToDB',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -meadow_type=> 'LOCAL',
	    -parameters => {
                'cmd' => qq(./scripts/sequencing/load_argenrich_qc_file.pl        )
		  . qq( --argenrich_file        #tempdir#/#argenrich_outfile#     )
		  . qq( --control_result_set_id #control_result_set_id#           )
		  . qq( --signal_result_set_id  #signal_result_set_id#            )
		  . qq( --user #tracking_db_user# --pass #tracking_db_pass# --host #tracking_db_host# --dbname #tracking_db_name# )
            },
        },
    ];
}

sub resource_classes {
  my $self = shift;
  return 
    {#todo add in lsf group spec here to top LSF warning output
     #todo pass DB_HOST_LSFNAME as param
     
     default                 => { 'LSF' => '' },    
     #urgent                  => { 'LSF' => '-q yesterday' },
     #Should never use this in the pipline, best to bswitch after submission if required
     normal_2GB              => { 'LSF' => ' -M2000 -R"select[mem>2000] rusage[mem=2000]"' },
     normal_monitored        => { 'LSF' => "" },
     normal_high_mem         => { 'LSF' => ' -M5000 -R"select[mem>5000] rusage[mem=5000]"' },
     normal_high_mem_2cpu    => { 'LSF' => ' -n2 -M5000 -R"select[mem>5000] rusage[mem=5000] span[hosts=1]"' },
     normal_monitored_2GB    => {'LSF' => " -M2000 -R\"select[mem>2000]".
                                                " rusage[mem=2000]\"" },
     normal_monitored_4GB    => {'LSF' => " -M4000 -R\"select[mem>4000] rusage[mem=4000]\"" },  
     normal_monitored_8GB    => {'LSF' => " -M8000 -R\"select[mem>8000] rusage[mem=8000]\"" },   
     normal_monitored_16GB   => {'LSF' => " -M16000 -R\"select[mem>16000] rusage[mem=16000]\"" }, 
     normal_16GB_2cpu        => {'LSF' => ' -n2 -M16000 -R"select[mem>16000] rusage[mem=16000] span[hosts=1]"' },
     normal_20GB_2cpu        => {'LSF' => ' -n2 -M20000 -R"select[mem>20000] rusage[mem=20000] span[hosts=1]"' }, 
     normal_25GB_2cpu        => {'LSF' => ' -n2 -M25000 -R"select[mem>25000] rusage[mem=25000] span[hosts=1]"' }, 
     normal_30GB_2cpu        => {'LSF' => ' -n2 -M30000 -R"select[mem>30000] rusage[mem=30000] span[hosts=1]"' },      
     normal_10gb_monitored   => {'LSF' => " -M10000 -R\"select[mem>10000] rusage[mem=10000]\"" },
     normal_5GB_2cpu_monitored => {'LSF' => " -n2 -M5000 -R\"select[mem>5000] rusage[mem=5000] span[hosts=1]\"" },
     normal_10gb             => { 'LSF' => ' -M10000 -R"select[mem>10000] rusage[mem=10000]"' },
     long_monitored          => { 'LSF' => "-q long " },
     long_high_mem           => { 'LSF' => '-q long -M4000 -R"select[mem>4000] rusage[mem=4000]"' },
     long_monitored_high_mem => { 'LSF' => "-q long -M4000 -R\"select[mem>4000] rusage[mem=4000]\"" },
    };
}

1;



