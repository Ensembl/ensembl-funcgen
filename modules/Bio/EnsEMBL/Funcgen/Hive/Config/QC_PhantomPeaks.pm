=pod 

=head1 NAME

    Bio::EnsEMBL::Funcgen::Hive::Config::QC_PhantomPeaks

=head1 SYNOPSIS

=head1 DESCRIPTION

=head1 LICENSE

    Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
    Copyright [2016-2017] EMBL-European Bioinformatics Institute

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

use Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf;
use base ('Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf');

sub pipeline_analyses {
    my ($self) = @_;
    return [
	{
	  -logic_name => 'BamFileQc',
	  -flow_into => {
	    MAIN => WHEN(
	      '#has_duplicates# eq "no"'  => 'QcPhantomPeaksJobFactory',
	      '!defined #has_duplicates#' => 'QcPhantomPeaksJobFactory',
	    ),
	  },
	},
        {   -logic_name => 'QcPhantomPeaksJobFactory',
            -module     => 'Bio::EnsEMBL::Funcgen::Hive::QcPhantomPeaksJobFactory',
#             -meadow_type=> 'LOCAL',
            -flow_into => { 
	      2 => 'QcRunPhantomPeaks4GB',
            },
        },
        {   -logic_name  => 'QcRunPhantomPeaks4GB',
            -module      => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -parameters  => { 
		  cmd => 
		    # Rscript does not search the path, so we use "which" to 
		    # do that. Also using single quotes to avoid interpolation
		    # of the dollar sign.
		    #
		    q( Rscript /software/ensembl/funcgen/spp_package/run_spp_nodups.R          )
		    # Overwrite plotfile, if one already exists
		  . qq(    -rf                             )
		  . qq(    -c=#bam_file#                   )
		  . qq(    -savp                           )
		  . qq(    -out=#phantom_peak_out_file#    )
		  . qq(    -odir=#tempdir#                 )
		  . qq(    -tmpdir=#tempdir#               )
		  # In case the job gets terminated for memlimit, this 
		  # ensures that the worker also dies. (or so we hope)
		  . qq(    ; sleep 30 )
            },
	    -rc_name    => 'normal_4GB_2cpu',
            -flow_into  => { 
	      MAIN     => 'QCLoadPhantomPeaksToDB',
	      MEMLIMIT => 'QcRunPhantomPeaks30GB',
            },
        },
        {   -logic_name  => 'QcRunPhantomPeaks30GB',
            -module      => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -parameters  => { 
		  cmd => 
		    # Rscript does not search the path, so we use "which" to 
		    # do that. Also using single quotes to avoid interpolation
		    # of the dollar sign.
		    #
		    q( Rscript /software/ensembl/funcgen/spp_package/run_spp_nodups.R )
		    # Overwrite plotfile, if one already exists
		  . qq(    -rf                             )
		  . qq(    -c=#bam_file#                   )
		  . qq(    -savp                           )
                  . qq(    -out=#phantom_peak_out_file#    )
                  . qq(    -odir=#tempdir#                 )
                  . qq(    -tmpdir=#tempdir#               )
		  # In case the job gets terminated for memlimit, this 
		  # ensures that the worker also dies. (or so we hope)
		  . qq(    ; sleep 30 )
            },
	    -rc_name    => 'normal_30GB_2cpu',
            -flow_into  => { 
	      MAIN => 'QCLoadPhantomPeaksToDB',
            },
        },
        {   -logic_name => 'QCLoadPhantomPeaksToDB',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
#             -meadow_type=> 'LOCAL',
	    -parameters => {
                cmd =>
		    qq( load_phantom_peak_file.pl                )
		  . qq(    --result_set_id #result_set_id#       )
		  . qq(    --result_file #phantom_peak_out_file# )
		  . qq(    --user   #tracking_db_user#   )
		  . qq(    --pass   #tracking_db_pass#   )
		  . qq(    --host   #tracking_db_host#   )
		  . qq(    --dbname #tracking_db_name#   )
		  . qq(    --work_dir #tempdir#          )
		  . qq(    --bam_file #bam_file#         )

            },
        },
    ];
}

1;
