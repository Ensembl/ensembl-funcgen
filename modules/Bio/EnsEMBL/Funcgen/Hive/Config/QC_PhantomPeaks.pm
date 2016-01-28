=pod 

=head1 NAME

    Bio::EnsEMBL::Funcgen::Hive::Config::QC_PhantomPeaks

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
package Bio::EnsEMBL::Funcgen::Hive::Config::QC_PhantomPeaks;

use strict;
use warnings;

use base ('Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf');  # All Hive databases configuration files should inherit from HiveGeneric, directly or indirectly

sub pipeline_analyses {
    my ($self) = @_;
    return [
        {   -logic_name => 'CreateFakeInputId',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
            -meadow_type=> 'LOCAL',
            -input_ids => [
	      {
		"dbID" => 2291,
		"fastq_files" => [],
		"output_dir" => "/lustre/scratch109/ensembl/funcgen/mn1/ersa/faang/alignments/homo_sapiens/GRCh38/3526",
		"result_set_groups" => {
		  "MEG01:hist:BR1_H3K27ac_3526_bwa_samse" => {
		    "dbIDs" => [2291,2292,2293,2294,2295],"set_names" => ["MEG01:hist:BR1_H3K27ac_3526_bwa_samse_TR3","MEG01:hist:BR1_H3K27ac_3526_bwa_samse_TR1","MEG01:hist:BR1_H3K27ac_3526_bwa_samse_TR2","MEG01:hist:BR1_H3K27ac_3526_bwa_samse_TR5","MEG01:hist:BR1_H3K27ac_3526_bwa_samse_TR4"]
		  },
		  "MEG01:hist:BR1_H3K4me3_3526_bwa_samse" => {
		    "dbIDs" => [2286,2287,2288,2289,2290],
		    "set_names" => [ 
		      "MEG01:hist:BR1_H3K4me3_3526_bwa_samse_TR3","MEG01:hist:BR1_H3K4me3_3526_bwa_samse_TR1","MEG01:hist:BR1_H3K4me3_3526_bwa_samse_TR4","MEG01:hist:BR1_H3K4me3_3526_bwa_samse_TR2","MEG01:hist:BR1_H3K4me3_3526_bwa_samse_TR5"
		    ]
		  },
		  "merged" => {
		    "dbIDs" => [2296],"set_names" => ["MEG01:hist:BR1_H3K27me3_3526_bwa_samse"]
		    }
		},
		"set_name" => "MEG01:hist:BR1_H3K27ac_3526_bwa_samse_TR3",
		"set_prefix" => "MEG01:hist:BR1_WCE_3526_bwa_samse",
		"set_type" => "ResultSet",

		chromosome_file => "/lustre/scratch109/ensembl/funcgen/mn1/ersa/faang/reference_files/CCAT/homo_sapiens_.CCAT_chr_lengths.txt", 
		batch_param_names => [		    
		  "no_write","rollback","full_delete","slices","skip_slices","result_set_only","result_set_mode",
		  "recover","alignment_analysis","peak_analysis","permissive_peaks","control_feature","no_idr",
		  "indexed_ref_fasta","idr_analysis","max_peaks","checksum_optional"
		],
		use_tracking_db => 1, 
		dnadb => {"-dnadb_host" => "ens-livemirror","-dnadb_name" => "homo_sapiens_core_82_38","-dnadb_pass" => "","-dnadb_port" => 3306,"-dnadb_user" => "ensro"}, 
		out_db => {"-dbname" => "mn1_faang2_tracking_homo_sapiens_funcgen_81_38","-host" => "ens-genomics1","-pass" => "ensembl","-port" => 3306,"-user" => "ensadmin"}, 
		pipeline_name => "blah", 
		data_root_dir => "/lustre/scratch109/ensembl/funcgen/mn1/ersa/faang/", 
		"alignment_analysis" => "bwa_samse",
		"checksum_optional" => 0,
		"dbID" => "2270",
		"set_name" => "F36P:hist:BR1_H3K4me3_3526_bwa_samse",
		
	      },
            ],
            -flow_into => { 
	      '1' => [ 'QcPhantomPeaksJobFactory' ],
            },
        },
# 	{
# 	    -logic_name => 'PreprocessAlignments',
# 	    -module     => 'Bio::EnsEMBL::Funcgen::Hive::CollectionWriter',
# 	    -flow_into => {
# 		'3'   => [  'QcPhantomPeaksJobFactory' ],
# 		'4'   => [  'QcPhantomPeaksJobFactory' ],
# 		'5'   => [  'QcPhantomPeaksJobFactory' ],
# 		'6'   => [  'QcPhantomPeaksJobFactory' ],
# 		'100' => [  'QcPhantomPeaksJobFactory' ],
# 	    },
# 	},
        {   -logic_name => 'QcPhantomPeaksJobFactory',
            -module     => 'Bio::EnsEMBL::Funcgen::Hive::QcPhantomPeaksJobFactory',
            -meadow_type=> 'LOCAL',
            -flow_into => { 
	      '2' => [ 'QcRunPhantomPeaks' ],
            },
        },
        {   -logic_name => 'QcRunPhantomPeaks',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -meadow_type=> 'LSF',
            -parameters => { 
		  cmd => 
		    qq( Rscript /software/ensembl/funcgen/spp_package/run_spp.R )
		  . qq(    -c=#bam_file# )
		  . qq(    -savp -out=#phantom_peak_out_file# )
            },
            -flow_into => { 
	      '1' => [ 'QCLoadPhantomPeaksToDB' ],
            },
        },
        {   -logic_name => 'QCLoadPhantomPeaksToDB',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -meadow_type=> 'LOCAL',
	    -parameters => {
                'cmd' =>
		    qq( load_phantom_peak_file.pl )
		  . qq(    --result_set_id #result_set_id# )
		  . qq(    --result_file #phantom_peak_out_file# )
		  . qq(    --user #tracking_db_user# --pass #tracking_db_pass# --host #tracking_db_host# --dbname #tracking_db_name# )
            },
        },
    ];
}

sub resource_classes {
  my $self = shift;
  return 
    {     
     default                 => { 'LSF' => '' },    
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



