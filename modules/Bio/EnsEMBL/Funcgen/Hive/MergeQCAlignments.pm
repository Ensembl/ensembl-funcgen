=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute

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

=head1 NAME

Bio::EnsEMBL::Hive::Funcgen::MergeQCAlignments

=head1 DESCRIPTION

Merges bam alignments from split (replicate or merged) fastq files.

=cut

package Bio::EnsEMBL::Funcgen::Hive::MergeQCAlignments;

use warnings;
use strict;

use Bio::EnsEMBL::Utils::Exception         qw( throw );
use Bio::EnsEMBL::Funcgen::Utils::EFGUtils qw( run_system_cmd );
use Bio::EnsEMBL::Funcgen::Sequencing::SeqTools;# merge_bams 

use base qw( Bio::EnsEMBL::Funcgen::Hive::BaseDB );

my %valid_flow_modes = (
  replicate => undef,
  merged    => undef,
  signal    => undef
); 

sub fetch_input {  
  my $self = shift;
  
  $self->SUPER::fetch_input();
  my $result_set = $self->fetch_Set_input('ResultSet');
  
  $self->get_param_method('bam_files',  'silent');
  $self->get_param_method('fastq_files',  'silent');
    
  if((! $self->bam_files ) && $self->fastq_files) {
    my $bam_file;
    $self->bam_files([ map {($bam_file = $_) =~ s/\.fastq_([0-9]+)$/.$1.bam/o; $bam_file} @{$self->fastq_files} ]);
  } elsif(! $self->bam_files) {
    $self->throw_no_retry('No bam_files or fastq_files have been defined');    
  }
  
  $self->get_param_method('set_prefix', 'required');  #This is control specific
  my $flow_mode = $self->get_param_method('flow_mode',  'required');
  $self->set_param_method('run_controls', 0); 
  
  if(! exists $valid_flow_modes{$flow_mode}) {
    throw("$flow_mode is now a valid flow_mode parameter, please specify one of:\t".
      join(' ', keys %valid_flow_modes));  
  }
  
  if(($flow_mode ne 'signal') && $self->get_param_method('result_set_groups', 'silent')) {
    throw("The $flow_mode flow mode is not valid for use with result_set_groups"); 
  }

  # flow mode is set to signal as a pipeline wide parameter in the merge 
  # analysis that processes the controls.
  #
  if($flow_mode eq 'signal') {
  
    # result_set_groups are present in the jobs of the control processing
    # analyses 
    #
    # - Preprocess_bwa_samse_control and
    # - MergeControlAlignments_and_QC 
    #
    # It looks like this:
    # 
    # "result_set_groups" => {
    #   "KU812:hist:BR2_H3K27ac_3526_bwa_samse" => {"dbIDs" => [2283,2284,2285,2286,2287],"set_names" => ["KU812:hist:BR2_H3K27ac_3526_bwa_samse_TR2","KU812:hist:BR2_H3K27ac_3526_bwa_samse_TR4","KU812:hist:BR2_H3K27ac_3526_bwa_samse_TR5","KU812:hist:BR2_H3K27ac_3526_bwa_samse_TR1","KU812:hist:BR2_H3K27ac_3526_bwa_samse_TR3"]},
    #   "KU812:hist:BR2_H3K4me3_3526_bwa_samse" => {"dbIDs" => [2288,2289,2290,2291,2292],"set_names" => ["KU812:hist:BR2_H3K4me3_3526_bwa_samse_TR5","KU812:hist:BR2_H3K4me3_3526_bwa_samse_TR4","KU812:hist:BR2_H3K4me3_3526_bwa_samse_TR1","KU812:hist:BR2_H3K4me3_3526_bwa_samse_TR2","KU812:hist:BR2_H3K4me3_3526_bwa_samse_TR3"]},
    #   "merged" => {"dbIDs" => [2293],"set_names" => ["KU812:hist:BR2_H3K27me3_3526_bwa_samse"]}
    # },
    # 
    $self->get_param_method('result_set_groups', 'required');
    $self->run_controls(1); 
  }

  if($flow_mode eq 'replicate'){
    $self->get_param_method('permissive_peaks', 'required');
  }
 
  $self->init_branching_by_analysis;
  return;
}


sub run {
  my $self       = shift;
  my $result_set = $self->ResultSet;
  my $cmd;

  return;
  
  ### CLEAN FASTQS ###
  if($self->fastq_files){
    #Run with no exit flag so we don't fail on retry
    $cmd = 'rm -f '.join(' ', @{$self->fastq_files});
    $self->helper->debug(3, "Removing fastq chunks:\n$cmd");
    run_system_cmd($cmd, 1);
  }
  
  ### MERGE BAMS ###
  my $file_prefix  = $self->get_alignment_path_prefix_by_ResultSet($result_set, $self->run_controls); 
  my $unfiltered_bam     = $file_prefix.'.bam';
  
  $self->helper->debug(1, "Merging bams to:\t".$unfiltered_bam); 
  
  merge_bams({
    input_bams => $self->bam_files, 
    output_bam => $unfiltered_bam, 
    debug => $self->debug,
  });
  
  my $tmp_bam = "${unfiltered_bam}.nodups.bam";
  
  remove_duplicates_from_bam({
    input_bam  => $unfiltered_bam,
    output_bam => $tmp_bam, 
    debug      => $self->debug,
  });
  
  unlink($unfiltered_bam);
  run_system_cmd("mv $tmp_bam $unfiltered_bam");

  return;
}

sub write_output {

  my $self = shift;

  my $result_set   = $self->ResultSet;
  my %batch_params = %{$self->batch_params};
  my $flow_mode    = $self->flow_mode;
  
  if ($flow_mode eq 'merged') {
  
    my %output_id = (
      set_type    => 'ResultSet',
      set_name    => $self->ResultSet->name,
      dbID        => $self->ResultSet->dbID,
      garbage     => $self->bam_files,
    );

    $self->branch_job_group('DefineMergedDataSet', [{%batch_params, %output_id}]);
  }
  
  if ($flow_mode eq 'replicate') {
  
    my %output_id = (
      set_type      => 'ResultSet',
      set_name      => $self->ResultSet->name,
      dbID          => $self->ResultSet->dbID,
      garbage       => $self->bam_files,
      
      # If flowing to replicate processing, set the permissive_peaks parameter.
      # permissive_peaks is always set to "SWEmbl_R0005" in our runs.
      #
      peak_analysis => $self->permissive_peaks,
    );

    $self->branch_job_group('run_'.$self->permissive_peaks.'_replicate', [{%batch_params, %output_id}]);
  }
  
  if ($flow_mode eq 'signal') {
    # Run signal fastqs
    my $result_set_groups               = $self->result_set_groups;
    my $result_set_analysis_logic_name  = $result_set->analysis->logic_name;
    
    if (exists $result_set_groups->{'merged'}) {

      my @merged_jobs;
      my $current_result_set_group = $result_set_groups->{'merged'};

      for my $i(0...$#{$current_result_set_group->{dbIDs}}) {
        push @merged_jobs, {
	  %batch_params,
	  garbage     => $self->bam_files, 
	  set_type    => 'ResultSet',
	  set_name    => $current_result_set_group->{set_names}->[$i],
	  dbID        => $current_result_set_group->{dbIDs}->[$i]
	};
      }
      # This flows to Preprocess_bwa_samse_merged.
      #
      $self->branch_job_group(
	'Preprocess_'.$result_set_analysis_logic_name.'_merged', 
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
	  garbage     => $self->bam_files, 
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
	'Preprocess_'.$result_set_analysis_logic_name.'_replicate',
	\@replicate_jobs,
	
	# Funnel job
	'PreprocessIDR', [
	  {
	    dbIDs     => $current_result_set_group->{dbIDs},
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
  
  $self->dataflow_job_groups;
  return;
}

1;
