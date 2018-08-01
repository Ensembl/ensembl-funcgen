package Bio::EnsEMBL::Funcgen::PipeConfig::PeakCalling::ResourceClasses;

use strict;
use warnings;
use base 'Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf';
use Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf;


sub resource_classes {
    my ($self) = @_;
    
    my $default_resource_classes = {
        %{$self->SUPER::resource_classes},  # inherit 'default' from the parent class

         '250Mb_job'      => {'LSF' => '-M250   -R"select[mem>250]   rusage[mem=250]"'   }     ,
         '500Mb_job'      => {'LSF' => '-M500   -R"select[mem>500]   rusage[mem=500]"'        },
         '1Gb_job'        => {'LSF' => '-M1000  -R"select[mem>1000]  rusage[mem=1000]"'       },
         '2Gb_job'        => {'LSF' => '-M2000  -R"select[mem>2000]  rusage[mem=2000]"'       },
         '4Gb_job'        => {'LSF' => '-M4000  -R"select[mem>4000]  rusage[mem=4000]"'       },
         '4Gb_job_2cpus'  => {'LSF' => '-M4000  -R"select[mem>4000]  rusage[mem=4000]" -n 2'  },
         '8Gb_job'        => {'LSF' => '-M8000  -R"select[mem>8000]  rusage[mem=8000]"'       },
         '16Gb_job'       => {'LSF' => '-M16000 -R"select[mem>16000] rusage[mem=16000]"'      },
         '24Gb_job'       => {'LSF' => '-M24000 -R"select[mem>24000] rusage[mem=24000]"'      },
         '32Gb_job'       => {'LSF' => '-M32000 -R"select[mem>32000] rusage[mem=32000]"'      },
         '32Gb_job_2cpus' => {'LSF' => '-M32000 -R"select[mem>32000] rusage[mem=32000]" -n 2' },
         '32Gb_job_3cpus' => {'LSF' => '-M32000 -R"select[mem>32000] rusage[mem=32000]" -n 3' },
         '48Gb_job'       => {'LSF' => '-M48000 -R"select[mem>48000] rusage[mem=48000]"'      },
         '64Gb_job'       => {'LSF' => '-M64000 -R"select[mem>64000] rusage[mem=64000]"'      },
         '64Gb_job_3cpus' => {'LSF' => '-M64000 -R"select[mem>64000] rusage[mem=64000]" -n 3' },

         'binarization' => {'LSF' => '-M16000 -R"select[mem>16000] rusage[mem=16000]" -n 5'                },
         'learn_model'  => {'LSF' => '-M31000 -R"span[hosts=1] select[mem>31000] rusage[mem=31000]" -n 12' },
    };
    
    my $production_resource_classes = {
        %{$self->SUPER::resource_classes},  # inherit 'default' from the parent class

         'default'        => {'LSF' => '-q production-rh7'   }     ,
         '250Mb_job'      => {'LSF' => '-q production-rh7 -M250   -R"select[mem>250]   rusage[mem=250]"'   }     ,
         '500Mb_job'      => {'LSF' => '-q production-rh7 -M500   -R"select[mem>500]   rusage[mem=500]"'        },
         '1Gb_job'        => {'LSF' => '-q production-rh7 -M1000  -R"select[mem>1000]  rusage[mem=1000]"'       },
         '2Gb_job'        => {'LSF' => '-q production-rh7 -M2000  -R"select[mem>2000]  rusage[mem=2000]"'       },
         '4Gb_job'        => {'LSF' => '-q production-rh7 -M4000  -R"select[mem>4000]  rusage[mem=4000]"'       },
         '4Gb_job_2cpus'  => {'LSF' => '-q production-rh7 -M4000  -R"select[mem>4000]  rusage[mem=4000]" -n 2'  },
         '8Gb_job'        => {'LSF' => '-q production-rh7 -M8000  -R"select[mem>8000]  rusage[mem=8000]"'       },
         '16Gb_job'       => {'LSF' => '-q production-rh7 -M16000 -R"select[mem>16000] rusage[mem=16000]"'      },
         '24Gb_job'       => {'LSF' => '-q production-rh7 -M24000 -R"select[mem>24000] rusage[mem=24000]"'      },
         '32Gb_job'       => {'LSF' => '-q production-rh7 -M32000 -R"select[mem>32000] rusage[mem=32000]"'      },
         '32Gb_job_2cpus' => {'LSF' => '-q production-rh7 -M32000 -R"select[mem>32000] rusage[mem=32000]" -n 2' },
         '32Gb_job_3cpus' => {'LSF' => '-q production-rh7 -M32000 -R"select[mem>32000] rusage[mem=32000]" -n 3' },
         '48Gb_job'       => {'LSF' => '-q production-rh7 -M48000 -R"select[mem>48000] rusage[mem=48000]"'      },
         '64Gb_job'       => {'LSF' => '-q production-rh7 -M64000 -R"select[mem>64000] rusage[mem=64000]"'      },
         '64Gb_job_3cpus' => {'LSF' => '-q production-rh7 -M64000 -R"select[mem>64000] rusage[mem=64000]" -n 3' },

         'binarization' => {'LSF' => '-q production-rh7 -M16000 -R"select[mem>16000] rusage[mem=16000]" -n 5'                },
         'learn_model'  => {'LSF' => '-q production-rh7 -M31000 -R"span[hosts=1] select[mem>31000] rusage[mem=31000]" -n 12' },

    };
    
    return $production_resource_classes;
}

1;

