=pod 

=head1 NAME

    Bio::EnsEMBL::Funcgen::Hive::Config::QC_PhantomPeaks

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
package Bio::EnsEMBL::Funcgen::Hive::Config::QC_PhantomPeaks;

use strict;
use warnings;

use Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf;
use base ('Bio::EnsEMBL::Funcgen::Hive::Config::Base');

sub pipeline_analyses {
    my ($self) = @_;
    return [
	{
	  -logic_name => 'BamFileQc',
	  -flow_into => {
	    MAIN => WHEN(
	      '#has_duplicates# eq "no"'  => 'qc_phantom_peaks_start',
	      '!defined #has_duplicates#' => 'qc_phantom_peaks_start',
	    ),
	  },
	},
        {   -logic_name => 'qc_phantom_peaks_start',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
            -flow_into => { 
              'MAIN->A' => 'QcPhantomPeaksJobFactory',
              'A->MAIN' => 'qc_phantom_peaks_done',
            },
        },
        {   -logic_name => 'QcPhantomPeaksJobFactory',
            -module     => 'Bio::EnsEMBL::Funcgen::Hive::QcPhantomPeaksJobFactory',
            -flow_into => { 
	      2 => 'QcRunPhantomPeaks4GB',
            },
        },
        {   -logic_name  => 'QcRunPhantomPeaks4GB',
            -module      => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -parameters  => {
                  use_bash_pipefail => 1,
		  cmd => 
		    # Rscript does not search the path, so we use "which" to 
		    # do that. Also using single quotes to avoid interpolation
		    # of the dollar sign.
		    #
		    q( Rscript $(which run_spp_nodups.R)   )
		    # Overwrite plotfile, if one already exists
		  . qq(    -rf                             )
		  . qq(    -c=#bam_file#                   )
		  . qq(    -savp                           )
		  . qq(    -out=#phantom_peak_out_file#    )
		  . qq(    -odir=#tempdir#                 )
		  . qq(    -tmpdir=#tempdir#               )

                  # Avoid the "Error: ignoring SIGPIPE signal error"
                  # from the EBI cluster
                  . qq(    &&                              )
		  # In case the job gets terminated for memlimit, this 
		  # ensures that the worker also dies. (or so we hope)
                  . qq(    sleep 30                        )
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
                  use_bash_pipefail => 1,
		  cmd => 
		    # Rscript does not search the path, so we use "which" to 
		    # do that. Also using single quotes to avoid interpolation
		    # of the dollar sign.
		    #
		    q( Rscript $(which run_spp_nodups.R)   )
		    # Overwrite plotfile, if one already exists
		  . qq(    -rf                             )
		  . qq(    -c=#bam_file#                   )
		  . qq(    -savp                           )
                  . qq(    -out=#phantom_peak_out_file#    )
                  . qq(    -odir=#tempdir#                 )
                  . qq(    -tmpdir=#tempdir#               )

                  # Avoid the "Error: ignoring SIGPIPE signal error"
                  # from the EBI cluster
                  . qq(    &&                              )
		  # In case the job gets terminated for memlimit, this 
		  # ensures that the worker also dies. (or so we hope)
		  . qq(    sleep 30                        )
            },
	    -rc_name    => 'normal_30GB_2cpu',
            -flow_into  => { 
	      MAIN => 'QCLoadPhantomPeaksToDB',
            },
        },
        {   -logic_name => 'QCLoadPhantomPeaksToDB',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
	    -parameters => {
                cmd =>
		    qq( load_phantom_peak_file.pl                )
		  . qq(    --result_set_id #result_set_id#       )
		  . qq(    --result_file #phantom_peak_out_file# )
		  . qq(    --user   #tracking_db_user#   )
                  . qq(    --pass   #tracking_db_pass#   )
                  . qq(    --port   #tracking_db_port#   )
		  . qq(    --host   #tracking_db_host#   )
		  . qq(    --dbname #tracking_db_name#   )
		  . qq(    --work_dir #tempdir#          )
		  . qq(    --bam_file #bam_file#         )

            },
        },
        {   -logic_name => 'qc_phantom_peaks_done',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
        },
    ];
}

1;
