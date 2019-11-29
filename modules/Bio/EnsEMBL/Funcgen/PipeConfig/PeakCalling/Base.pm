package Bio::EnsEMBL::Funcgen::PipeConfig::PeakCalling::Base;

use strict;
use warnings;
use base 'Bio::EnsEMBL::Funcgen::PipeConfig::PeakCalling::ResourceClasses';
use Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf;


sub beekeeper_extra_cmdline_options {
  my $self = shift;
  return '-reg_conf ' . $self->o('reg_conf') . ' -keep_alive -can_respecialize 1 -sleep 0.5';
}

sub default_options {
  my $self = shift;
  
  return {
      %{$self->SUPER::default_options},
      pipeline_name => 'chip_seq_analysis',
   };
}

sub pipeline_wide_parameters {
  my $self = shift;
  
  return {
    %{$self->SUPER::pipeline_wide_parameters},
    pipeline_name           => $self->o('pipeline_name'),
    tempdir                 => $self->o('tempdir'),
    tempdir_peak_calling    => $self->o('tempdir') . '/peak_calling',
    data_root_dir           => $self->o('data_root_dir'),
    reports_dir             => $self->o('reports_dir'),
    ensembl_release_version => $self->o('ensembl_release_version'),
    reference_data_root_dir => $self->o('reference_data_root_dir'),
    reg_conf                => $self->o('reg_conf'),
  };
}

sub generate_parallel_alignment_analyses {
    my $self  = shift;
    my $param = shift;
    
    my $start     = $param->{start};
    my $prefix    = $param->{prefix};
    my $suffix    = $param->{suffix};
    my $flow_into = $param->{flow_into};
    my $after     = $param->{after};

    my $surround = sub {
        my $name = shift;
        return $prefix . $name . $suffix
    };
    
    if (! defined $start) {
        $start = $surround->('start_align');
    }
    if (! defined $flow_into) {
        $flow_into = {};
    }

    return [
        {   -logic_name  => $start,
            -module      => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
            -flow_into   => {
               'MAIN->A' => $surround->('split_experiment'),
               'A->MAIN' => $surround->('done_align'),
            },
        },
        {   -logic_name  => $surround->('split_experiment'),
            -module     => 'Bio::EnsEMBL::Funcgen::RunnableDB::PeakCalling::SplitExperimentReads',
            -parameters => {
              tempdir => '#tempdir_peak_calling#/#species#/alignments'
            },
            -flow_into   => {
               '3->A' => $surround->('split_fastq_files'),
               'A->2' => $surround->('merge'),
            },
        },
        {   -logic_name  => $surround->('split_fastq_files'),
            -module     => 'Bio::EnsEMBL::Funcgen::RunnableDB::PeakCalling::SplitFastq',
            -rc_name    => '4Gb_job',
            -batch_size => 3,
            -parameters => {
              tempdir => '#tempdir_peak_calling#/#species#/alignments'
            },
            -flow_into   => {
               '2->A' => $surround->('align_quick'),
               'A->3' => $surround->('merge_chunks'),
            },
        },
        {   -logic_name  => $surround->('align_quick'),
            -module     => 'Bio::EnsEMBL::Funcgen::RunnableDB::PeakCalling::AlignFastqFile',
            -priority   => 10,
            -batch_size => 10,
            -rc_name    => '8Gb_job_2h',
            -flow_into  => {
              RUNLIMIT => $surround->('align_slow'),
            }
        },
        {   -logic_name  => $surround->('align_slow'),
            -module     => 'Bio::EnsEMBL::Funcgen::RunnableDB::PeakCalling::AlignFastqFile',
            -priority   => 10,
            -rc_name    => '8Gb_job_8h',
        },
        {   -logic_name  => $surround->('merge_chunks'),
            -module      => 'Bio::EnsEMBL::Funcgen::RunnableDB::PeakCalling::MergeBamFiles',
            -priority   => 20,
        },
        {   -logic_name  => $surround->('merge'),
            -module     => 'Bio::EnsEMBL::Funcgen::RunnableDB::PeakCalling::MergeBamFiles',
            -priority   => 30,
            -flow_into   => {
               MAIN     => $surround->('remove_duplicates'),
               MEMLIMIT => $surround->('merge_himem'),
            },
        },
        {   -logic_name  => $surround->('merge_himem'),
            -module     => 'Bio::EnsEMBL::Funcgen::RunnableDB::PeakCalling::MergeBamFiles',
            -priority   => 30,
            -rc_name    => '8Gb_job',
            -flow_into   => {
               MAIN => $surround->('remove_duplicates'),
            },
        },
        {   -logic_name  => $surround->('remove_duplicates'),
            -module     => 'Bio::EnsEMBL::Funcgen::RunnableDB::PeakCalling::RemoveDuplicates',
            -priority   => 40,
            -flow_into   => {
               MAIN     => $surround->('check_alignments_to_y_chromosome'),
               MEMLIMIT => $surround->('remove_duplicates_himem'),
            },
        },
        {   -logic_name  => $surround->('remove_duplicates_himem'),
            -module     => 'Bio::EnsEMBL::Funcgen::RunnableDB::PeakCalling::RemoveDuplicates',
            -priority   => 40,
            -rc_name    => '8Gb_job',
            -flow_into   => {
               MAIN => $surround->('check_alignments_to_y_chromosome'),
            },
        },
        {   -logic_name  => $surround->('check_alignments_to_y_chromosome'),
            -module     => 'Bio::EnsEMBL::Funcgen::RunnableDB::PeakCalling::CheckAlignmentsToYChromosome',
            -flow_into   => {
               MAIN => $surround->('register_alignment'),
            },
        },
        {   -logic_name  => $surround->('register_alignment'),
            -module     => 'Bio::EnsEMBL::Funcgen::RunnableDB::PeakCalling::RegisterAlignment',
            -flow_into   => {
               MAIN => $surround->('check_alignments'),
            },
        },
        {   -logic_name  => $surround->('check_alignments'),
            -module           => 'Bio::EnsEMBL::DataCheck::Pipeline::RunDataChecks',
            -max_retry_count  => 3,
            -parameters => {
              registry_file    => '#reg_conf#',
              species          => '#species#',
              group            => 'funcgen',
              datacheck_names  => [ 'ControlAlignmentNamingConvention' ],
            },
        },
        {   -logic_name  => $surround->('done_align'),
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
            -priority   => 50,
            -flow_into  => $flow_into,
        },
    ]
}

1;

