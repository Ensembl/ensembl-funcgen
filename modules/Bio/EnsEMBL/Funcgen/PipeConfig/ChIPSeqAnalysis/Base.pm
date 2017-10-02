package Bio::EnsEMBL::Funcgen::PipeConfig::ChIPSeqAnalysis::Base;

use strict;
use warnings;
use base 'Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf';
use Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf;

sub beekeeper_extra_cmdline_options {
  my $self = shift;
  return '-reg_conf ' . $self->o('reg_conf') . ' -keep_alive -can_respecialize 1 -sleep 0.1';
}

sub default_options {
  my $self = shift;
  
  return {
      %{$self->SUPER::default_options},
      pipeline_name => 'chip_seq_analysis',
      in_test_mode  => 0,
   };
}

sub pipeline_wide_parameters {
  my $self = shift;
  
  return {
    %{$self->SUPER::pipeline_wide_parameters},
    pipeline_name           => $self->o('pipeline_name'),
    tempdir                 => $self->o('tempdir'),
    data_root_dir           => $self->o('data_root_dir'),
    ensembl_release_version => $self->o('ensembl_release_version'),
    reference_data_root_dir => $self->o('reference_data_root_dir'),
    in_test_mode            => $self->o('in_test_mode'),
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
               'MAIN->A' => $surround->('split'),
               'A->MAIN' => $surround->('done_align'),
            },
        },
        {   -logic_name  => $surround->('split'),
            -module     => 'Bio::EnsEMBL::Funcgen::RunnableDB::ChIPSeq::SplitFastq',
            -flow_into   => {
               '2->A' => $surround->('align'),
               'A->3' => $surround->('merge'),
            },
        },
        {   -logic_name  => $surround->('align'),
            -module     => 'Bio::EnsEMBL::Funcgen::RunnableDB::ChIPSeq::AlignFastqFile',
            -rc_name    => '8Gb_job',
        },
        {   -logic_name  => $surround->('merge'),
            -module     => 'Bio::EnsEMBL::Funcgen::RunnableDB::ChIPSeq::MergeBamFiles',
            -flow_into   => {
               MAIN => $surround->('remove_duplicates'),
            },
        },
        {   -logic_name  => $surround->('remove_duplicates'),
            -module     => 'Bio::EnsEMBL::Funcgen::RunnableDB::ChIPSeq::RemoveDuplicates',
            -flow_into   => {
               MAIN => $surround->('register_alignment'),
            },
        },
        {   -logic_name  => $surround->('register_alignment'),
            -module     => 'Bio::EnsEMBL::Funcgen::RunnableDB::ChIPSeq::RegisterAlignment',
        },
        {   -logic_name  => $surround->('done_align'),
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
            -flow_into  => $flow_into,
        },
    ]
}

sub resource_classes {
    my ($self) = @_;
    return {
        %{$self->SUPER::resource_classes},  # inherit 'default' from the parent class

         '250Mb_job'    => {'LSF' => '-M250   -R"select[mem>250]   rusage[mem=250]"'   },
         '500Mb_job'    => {'LSF' => '-M500   -R"select[mem>500]   rusage[mem=500]"'   },
         '1Gb_job'      => {'LSF' => '-M1000  -R"select[mem>1000]  rusage[mem=1000]"'  },
         '2Gb_job'      => {'LSF' => '-M2000  -R"select[mem>2000]  rusage[mem=2000]"'  },
         '4Gb_job'      => {'LSF' => '-M4000  -R"select[mem>4000]  rusage[mem=4000]"'  },
         '8Gb_job'      => {'LSF' => '-M8000  -R"select[mem>8000]  rusage[mem=8000]"'  },
         '16Gb_job'     => {'LSF' => '-M16000 -R"select[mem>16000] rusage[mem=16000]"' },
         '24Gb_job'     => {'LSF' => '-M24000 -R"select[mem>24000] rusage[mem=24000]"' },
         '32Gb_job'     => {'LSF' => '-M32000 -R"select[mem>32000] rusage[mem=32000]"' },
         '48Gb_job'     => {'LSF' => '-M48000 -R"select[mem>48000] rusage[mem=48000]"' },
         '64Gb_job'     => {'LSF' => '-M64000 -R"select[mem>64000] rusage[mem=64000]"' },
    };
}

1;

