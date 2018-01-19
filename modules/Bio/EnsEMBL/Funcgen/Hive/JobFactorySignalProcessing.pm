=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2018] EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

  Questions may also be sent to the Ensembl help desk at
  <http://www.ensembl.org/Help/Contact>.

=cut

package Bio::EnsEMBL::Funcgen::Hive::JobFactorySignalProcessing;

use strict;

use base qw( Bio::EnsEMBL::Funcgen::Hive::BaseDB );

sub fetch_input {
  my $self = shift;
  
  # This sets out_db which is needed to get a ResultSetAdapter in the next 
  # command.
  #
  $self->SUPER::fetch_input();
  
  $self->fetch_Set_input('ResultSet');
  $self->get_param_method('result_set_groups', 'required');
  
  return;
}

sub run {
  my $self = shift;

  my $result_set   = $self->ResultSet;
  my %batch_params = %{$self->batch_params};

  # Run signal fastqs
  my $result_set_groups               = $self->result_set_groups;
  my $result_set_analysis_logic_name  = $result_set->analysis->logic_name;
  
  if (exists $result_set_groups->{'merged'}) {

    my @merged_jobs;
    my $current_result_set_group = $result_set_groups->{'merged'};

    for my $i(0...$#{$current_result_set_group->{dbIDs}}) {
      push @merged_jobs, {
	%batch_params,
	set_type    => 'ResultSet',
	set_name    => $current_result_set_group->{set_names}->[$i],
	dbID        => $current_result_set_group->{dbIDs}->[$i]
      };
    }
    # This flows to Preprocess_bwa_samse_merged.
    #
    $self->branch_job_group(
      #'Preprocess_'.$result_set_analysis_logic_name.'_merged', 
      10,
      \@merged_jobs
    );
=head2

This is what is being seeded in the pipeline:

Branch 10 is Preprocess_bwa_samse_merged 

[
  10,
  [
    {
      'dbID' => 2323,
      'set_type' => 'ResultSet',
      'checksum_optional' => 0,
      'set_name' => 'F36P:hist:BR1_H3K27me3_3526_bwa_samse',
      'garbage' => [
		      '/lustre/scratch109/ensembl/funcgen/mn1/ersa/mn1_refactorMerge2_tracking_homo_sapiens_funcgen_81_38/output/ersa/alignments/homo_sapiens/GRCh38/3526/F36P:hist:BR1_WCE_3526_bwa_samse.0000.bam',
		      '/lustre/scratch109/ensembl/funcgen/mn1/ersa/mn1_refactorMerge2_tracking_homo_sapiens_funcgen_81_38/output/ersa/alignments/homo_sapiens/GRCh38/3526/F36P:hist:BR1_WCE_3526_bwa_samse.0001.bam',
		      '/lustre/scratch109/ensembl/funcgen/mn1/ersa/mn1_refactorMerge2_tracking_homo_sapiens_funcgen_81_38/output/ersa/alignments/homo_sapiens/GRCh38/3526/F36P:hist:BR1_WCE_3526_bwa_samse.0002.bam',
		      '/lustre/scratch109/ensembl/funcgen/mn1/ersa/mn1_refactorMerge2_tracking_homo_sapiens_funcgen_81_38/output/ersa/alignments/homo_sapiens/GRCh38/3526/F36P:hist:BR1_WCE_3526_bwa_samse.0003.bam',
		      '/lustre/scratch109/ensembl/funcgen/mn1/ersa/mn1_refactorMerge2_tracking_homo_sapiens_funcgen_81_38/output/ersa/alignments/homo_sapiens/GRCh38/3526/F36P:hist:BR1_WCE_3526_bwa_samse.0004.bam',
		      '/lustre/scratch109/ensembl/funcgen/mn1/ersa/mn1_refactorMerge2_tracking_homo_sapiens_funcgen_81_38/output/ersa/alignments/homo_sapiens/GRCh38/3526/F36P:hist:BR1_WCE_3526_bwa_samse.0005.bam',
		      '/lustre/scratch109/ensembl/funcgen/mn1/ersa/mn1_refactorMerge2_tracking_homo_sapiens_funcgen_81_38/output/ersa/alignments/homo_sapiens/GRCh38/3526/F36P:hist:BR1_WCE_3526_bwa_samse.0006.bam',
		      '/lustre/scratch109/ensembl/funcgen/mn1/ersa/mn1_refactorMerge2_tracking_homo_sapiens_funcgen_81_38/output/ersa/alignments/homo_sapiens/GRCh38/3526/F36P:hist:BR1_WCE_3526_bwa_samse.0007.bam',
		      '/lustre/scratch109/ensembl/funcgen/mn1/ersa/mn1_refactorMerge2_tracking_homo_sapiens_funcgen_81_38/output/ersa/alignments/homo_sapiens/GRCh38/3526/F36P:hist:BR1_WCE_3526_bwa_samse.0008.bam',
		      '/lustre/scratch109/ensembl/funcgen/mn1/ersa/mn1_refactorMerge2_tracking_homo_sapiens_funcgen_81_38/output/ersa/alignments/homo_sapiens/GRCh38/3526/F36P:hist:BR1_WCE_3526_bwa_samse.0009.bam',
		      '/lustre/scratch109/ensembl/funcgen/mn1/ersa/mn1_refactorMerge2_tracking_homo_sapiens_funcgen_81_38/output/ersa/alignments/homo_sapiens/GRCh38/3526/F36P:hist:BR1_WCE_3526_bwa_samse.0010.bam',
		      '/lustre/scratch109/ensembl/funcgen/mn1/ersa/mn1_refactorMerge2_tracking_homo_sapiens_funcgen_81_38/output/ersa/alignments/homo_sapiens/GRCh38/3526/F36P:hist:BR1_WCE_3526_bwa_samse.0011.bam'
		    ],
      'alignment_analysis' => 'bwa_samse'
    }
  ]
],

=cut
    delete $result_set_groups->{'merged'};
  }
    
  foreach my $current_result_set_group_name (keys %{$result_set_groups}) {
  
    # This:
    #
    # print "current_result_set_group_name = $current_result_set_group_name\n";
    #
    # would print this:
    # 
    # current_result_set_group_name = F36P:hist:BR1_H3K4me3_3526_bwa_samse
    # current_result_set_group_name = F36P:hist:BR1_H3K27ac_3526_bwa_samse
    # 
    
    my @replicate_jobs;
    my $current_result_set_group = $result_set_groups->{$current_result_set_group_name};

    # This builds the jobs for the signal fastqs. The job descriptions are
    # taken from the "result_set_groups". This is a hash that holds the
    # information to build the next sets of jobs.
    #
    for my $i(0...$#{$current_result_set_group->{dbIDs}}) {
      push @replicate_jobs, {
	%batch_params,
	set_type    => 'ResultSet',
	set_name    => $current_result_set_group->{set_names}->[$i],
	dbID        => $current_result_set_group->{dbIDs}->[$i]
      };
    }

    # This is run when flowing from the merging of controls. The fan
    # job goes to Preprocess_bwa_samse_replicate, the funnel job to
    # PreprocessIDR.
    #
    $self->branch_job_group(
    
      # Fan job
      #'Preprocess_'.$result_set_analysis_logic_name.'_replicate',
      
      11,
      \@replicate_jobs,
      
      # Funnel job
      #'PreprocessIDR', 
      3,
      [
	{
	  # These are result set ids
	  #
	  dbIDs     => $current_result_set_group->{dbIDs},
	  
	  # These are the names of the results sets. This is a redundant 
	  # piece of information.
	  #
	  set_names => $current_result_set_group->{set_names},
	  set_type  => 'ResultSet',
	  %batch_params,
	}
      ]
    );

=head2

This is what is being seeded in the pipeline:

Branch 11 is Preprocess_bwa_samse_replicate 
Branch  3 is PreprocessIDR

[
            11,
            [
              {
                'dbID' => 2303,
                'set_type' => 'ResultSet',
                'checksum_optional' => 0,
                'set_name' => 'F36P:hist:BR1_H3K4me3_3526_bwa_samse_TR5',
                'garbage' => $VAR1->[0][1][0]{'garbage'},
                'alignment_analysis' => 'bwa_samse'
              },
              {
                'dbID' => 2304,
                'set_type' => 'ResultSet',
                'checksum_optional' => 0,
                'set_name' => 'F36P:hist:BR1_H3K4me3_3526_bwa_samse_TR6',
                'garbage' => $VAR1->[0][1][0]{'garbage'},
                'alignment_analysis' => 'bwa_samse'
              },
              {
                'dbID' => 2305,
                'set_type' => 'ResultSet',
                'checksum_optional' => 0,
                'set_name' => 'F36P:hist:BR1_H3K4me3_3526_bwa_samse_TR1',
                'garbage' => $VAR1->[0][1][0]{'garbage'},
                'alignment_analysis' => 'bwa_samse'
              },
              {
                'dbID' => 2306,
                'set_type' => 'ResultSet',
                'checksum_optional' => 0,
                'set_name' => 'F36P:hist:BR1_H3K4me3_3526_bwa_samse_TR3',
                'garbage' => $VAR1->[0][1][0]{'garbage'},
                'alignment_analysis' => 'bwa_samse'
              },
              {
                'dbID' => 2307,
                'set_type' => 'ResultSet',
                'checksum_optional' => 0,
                'set_name' => 'F36P:hist:BR1_H3K4me3_3526_bwa_samse_TR9',
                'garbage' => $VAR1->[0][1][0]{'garbage'},
                'alignment_analysis' => 'bwa_samse'
              },
              {
                'dbID' => 2308,
                'set_type' => 'ResultSet',
                'checksum_optional' => 0,
                'set_name' => 'F36P:hist:BR1_H3K4me3_3526_bwa_samse_TR2',
                'garbage' => $VAR1->[0][1][0]{'garbage'},
                'alignment_analysis' => 'bwa_samse'
              },
              {
                'dbID' => 2309,
                'set_type' => 'ResultSet',
                'checksum_optional' => 0,
                'set_name' => 'F36P:hist:BR1_H3K4me3_3526_bwa_samse_TR8',
                'garbage' => $VAR1->[0][1][0]{'garbage'},
                'alignment_analysis' => 'bwa_samse'
              },
              {
                'dbID' => 2310,
                'set_type' => 'ResultSet',
                'checksum_optional' => 0,
                'set_name' => 'F36P:hist:BR1_H3K4me3_3526_bwa_samse_TR4',
                'garbage' => $VAR1->[0][1][0]{'garbage'},
                'alignment_analysis' => 'bwa_samse'
              },
              {
                'dbID' => 2311,
                'set_type' => 'ResultSet',
                'checksum_optional' => 0,
                'set_name' => 'F36P:hist:BR1_H3K4me3_3526_bwa_samse_TR10',
                'garbage' => $VAR1->[0][1][0]{'garbage'},
                'alignment_analysis' => 'bwa_samse'
              },
              {
                'dbID' => 2312,
                'set_type' => 'ResultSet',
                'checksum_optional' => 0,
                'set_name' => 'F36P:hist:BR1_H3K4me3_3526_bwa_samse_TR7',
                'garbage' => $VAR1->[0][1][0]{'garbage'},
                'alignment_analysis' => 'bwa_samse'
              }
            ],
            3,
            [
              {
                'dbIDs' => [
                             2303,
                             2304,
                             2305,
                             2306,
                             2307,
                             2308,
                             2309,
                             2310,
                             2311,
                             2312
                           ],
                'set_names' => [
                                 'F36P:hist:BR1_H3K4me3_3526_bwa_samse_TR5',
                                 'F36P:hist:BR1_H3K4me3_3526_bwa_samse_TR6',
                                 'F36P:hist:BR1_H3K4me3_3526_bwa_samse_TR1',
                                 'F36P:hist:BR1_H3K4me3_3526_bwa_samse_TR3',
                                 'F36P:hist:BR1_H3K4me3_3526_bwa_samse_TR9',
                                 'F36P:hist:BR1_H3K4me3_3526_bwa_samse_TR2',
                                 'F36P:hist:BR1_H3K4me3_3526_bwa_samse_TR8',
                                 'F36P:hist:BR1_H3K4me3_3526_bwa_samse_TR4',
                                 'F36P:hist:BR1_H3K4me3_3526_bwa_samse_TR10',
                                 'F36P:hist:BR1_H3K4me3_3526_bwa_samse_TR7'
                               ],
                'checksum_optional' => 0,
                'set_type' => 'ResultSet',
                'alignment_analysis' => 'bwa_samse'
              }
            ]
          ]

=cut

  }
}

sub write_output {
  my $self = shift;
  $self->dataflow_job_groups;
  return;
}

1;
