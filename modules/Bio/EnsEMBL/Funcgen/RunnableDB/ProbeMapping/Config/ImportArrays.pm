# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2020] EMBL-European Bioinformatics Institute
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::Analysis::Config::ImportArrays

=head1 SYNOPSIS

    use Bio::EnsEMBL::Analysis:Config::ImportAarrys;

=head1 DESCRIPTION

  This contains the configuration for importing arrays from flat files.

=head1 CONTACT

=cut

package Bio::EnsEMBL::Funcgen::RunnableDB::ProbeMapping::Config::ImportArrays;

use warnings ;
use strict;
use vars qw( %Config );

# Hash containing config info
# -- one hashnode per logic name, with a 'DEFAULT' logic name provided
#

%Config =
  (
   ARRAY_CONFIG => {
    DEFAULT => {

     # Regular expression for parsing file headers
     #
     IIDREGEXP =>  '^>probe:(\S+):(\S+):(\S+:\S+;).*$',

     IFIELDORDER  => {},
     ARRAY_PARAMS => {},
     ARRAYS_WITH_DEFAULT_PARAMS => [],
    },



    # AFFY UTR

    IMPORT_AFFY_UTR_ARRAYS => {

       IIDREGEXP => '^>probe:(\S+?):(\S+?):(\S+;).*$',
       #IIDREGEXP => '^>probe:(\S+?):(\S+?):(\S+:\S+;).*$',

      IFIELDORDER => {
        -name       => 2,
        -array_chip => 0,
        -array      => 0,
        -probe_set  => 1
      },

      ARRAY_PARAMS => {

      Default =>  {
        -vendor  => 'AFFY',
        -format  => 'EXPRESSION',
        -type    => 'OLIGO',
        -class   => 'AFFY_UTR',
        
        -is_probeset_array       => 1,
        -is_linked_array         => 0,
        -has_sense_interrogation => 0,

      },

      platypus_exon => {
        -name    => 'platypus_exon',
        -vendor  => 'AFFY',
        -format  => 'EXPRESSION',
        -type    => 'OLIGO',
        -class   => 'AFFY_UTR',
        -is_probeset_array       => 1,
        -is_linked_array         => 0,
        -has_sense_interrogation => 0,
      },
     },

     INPUT_FORMAT => 'FASTA',

     ARRAYS_WITH_DEFAULT_PARAMS => [
        'Porcine',
        'X_tropicalis',
        'Canine_2',
        'Rhesus',
        'CINT06a520380F',
        'Bovine',
        'Chicken',
        'C_elegans',
        'Zebrafish',
        'RAE230A',
        'RAE230B',
        'Rat230_2',
        'RG-U34A',
        'RG-U34B',
        'RG-U34C',
        'RN-U34',
        'RT-U34',
        'RG-U34B',
        'PrimeView',
        'HC-G110',
        'U133_X3P',
        'HuGeneFL',
        'HG_U95A',
        'HG-U95E',
        'HG-U95D',
        'HG-U95C',
        'HG-U95B',
        'HG_U95Av2',
        'HG-U133_Plus_2',
        'HG-U133B',
        'HG-U133A',
        'HG-U133A_2',
        'HG-Focus',
        'MG-U74Cv2',
        'MG-U74A',
        'MG-U74Av2',
        'MG-U74B',
        'MG-U74Bv2',
        'MG-U74C',
        'MOE430A',
        'MOE430B',
        'Mouse430A_2',
        'Mouse430_2',
        'Mu11KsubA',
        'Mu11KsubB',
        'DrosGenome1',
        'Drosophila_2',
        'Yeast_2',
        'YG-S98',
        'E_coli_2',
        'E_coli_Antisense',
        'S_aureus',

        ## Plant arrays
        'ATH1-121501',
        'Barley1',
        'Cotton',		## I think GPL8672
        'Medicago',		## I think GPL4652
        'Poplar',		## I think GPL4359
        'Rice',			## I think GPL2025
        'Sorghum',
        'Soybean',
        'Tomato',		## 
        'Vitis_Vinifera',	## I think GPL1320
        'maize',		## I think GPL4032
        'GPL16720',             ## custom maize array
        'wheat',		## I think GPL3802

      ],

    },



    IMPORT_AFFY_ST_ARRAYS => {

     IIDREGEXP => '^>probe:(\S+?):(\S+?);.*[TranscriptCluster|ProbeSet]ID=(\S+);',
     

     #Can't use ProbeID=([0-9]+) as control probes only have there ProbeID in the concat'd full name string
     #Hence the match will fail.

     IFIELDORDER => {
        -name       => 1,
        -array_chip => 0,
        -array      => 0,
        -probe_set   => 2,
      },

     ARRAY_PARAMS => {
        'Default' => {
          -vendor => 'AFFY',
          -format  => 'EXPRESSION',
          -type    => 'OLIGO',
          -class   => 'AFFY_ST',
          
          -is_probeset_array       => 1,
          -is_linked_array         => 1,
          -has_sense_interrogation => 1,

        },

      'HTA-2_0' => {
        -name => 'HTA-2_0',
        -vendor => 'AFFY',
        -format  => 'EXPRESSION',
        -type    => 'OLIGO',
        -description => 'Human Transcriptome Array 2.0',
        -class   => 'AFFY_ST',

          -is_probeset_array       => 1,
          -is_linked_array         => 1,
          -has_sense_interrogation => 1,
      },
     },

     INPUT_FORMAT => 'FASTA',

     ARRAYS_WITH_DEFAULT_PARAMS => [
      'RaEx-1_0-st-v1',
      'RaGene-1_0-st-v1',
      'RaGene-2_1-st-v1',
      'HuGene-1_0-st-v1',
      'HuGene-2_0-st-v1',
      'HuEx-1_0-st-v2',
      'MoGene-1_0-st-v1',
      'MoEx-1_0-st-v1',
      'MoGene-2_1-st-v1',
      'CynGene-1_0-st-v1',
      'CyRGene-1_0-st-v1',
      'RheGene-1_0-st-v1',
      'RheGene-1_1-st-v1',
      'GPL19230',
      'ZebGene-1_1-st-v1',
      'ZebGene-1_0-st-v1',
      'EquGene-1_0-st-v1',
      'EquGene-1_1-st-v1',

      ## Plant arrays
      'AraGene-1_1-st-v1',
      'MedGene-1_0-st-v1',	## I think GPL18240
      'RCnGene-1_1-st-v1',
      'RJpGene-1_1-st-v1',
      'SoyGene-1_1-st-v1',
      'TomGene-1_1-st-v1',
      'CHOGene-2_0-st-v1',
      'FelGene-1_0-st-v1',
      'CHOGene-2_1-st-v1',
      'FelGene-1_1-st-v1',
      
     ],
    },



    IMPORT_ILLUMINA_WG_ARRAYS => {
    
     IIDREGEXP => '^>(\S+):(\S+).*$',

     IFIELDORDER => {
        -name       => 1,
        -array_chip => 0,
        -array      => 0,
      },

     ARRAY_PARAMS => {

      'Default' => {
        -vendor => 'ILLUMINA',
        -format  => 'EXPRESSION',
        -type    => 'OLIGO',
        -class   => 'ILLUMINA_WG',

        -is_probeset_array      => 0,
        -is_linked_array        => 0,
        -has_sense_interrogation  => 0,

      },
      
      'RatRef-12' => {
        -name => 'RatRef-12_V1',
        -vendor => 'ILLUMINA',
        -format  => 'EXPRESSION',
        -type    => 'OLIGO',
        -class   => 'ILLUMINA_WG',

        -is_probeset_array      => 0,
        -is_linked_array        => 0,
        -has_sense_interrogation  => 0,

      },
     },

     INPUT_FORMAT => 'FASTA',

     ARRAYS_WITH_DEFAULT_PARAMS => [
      'MouseWG_6_V1',
      'MouseWG_6_V2',
      'MouseRef-8_V2',
      'MouseRef-8',
      'HumanWG_6_V1',
      'HumanHT-12_V3',
      'HumanHT-12_V4',
      'HumanRef-8_V3',
      'HumanWG_6_V2',
      'HumanWG_6_V3',
     ],
    },



    IMPORT_ILLUMINA_INFINIUM_ARRAYS => {
    
     IIDREGEXP => '^>(\S+):(\S+).*$',
     
     IFIELDORDER => {
        -name       => 1,
        -array_chip => 0,
        -array      => 0,
      },
    
    ARRAY_PARAMS => {
      'Default' => {
        -vendor  => 'ILLUMINA',
        -format  => 'METHYLATION',
        -type    => 'OLIGO',
        -class   => 'ILLUMINA_INFINIUM',

        -is_probeset_array      => 0,
        -is_linked_array        => 0,
        -has_sense_interrogation  => 0,

      },
    },
    INPUT_FORMAT => 'FASTA',

     ARRAYS_WITH_DEFAULT_PARAMS => [
      'HumanMethylation27',
      'HumanMethylation450',
    ],
  },



    # CODELINK

    IMPORT_CODELINK_ARRAYS => {
    
     IIDREGEXP => '^>(\S+):(\S+).*$',

     IFIELDORDER => {
        -name       => 1,
        -array_chip => 0,
        -array      => 0,
      },

     ARRAY_PARAMS => {
      'Default' => {
        -vendor => 'CODELINK',
        -format  => 'EXPRESSION',
        -type    => 'OLIGO',
        -class   => 'CODELINK',

        -is_probeset_array      => 0,
        -is_linked_array        => 1,
        -has_sense_interrogation  => 0,
      },
     },

     INPUT_FORMAT => 'FASTA',

     ARRAYS_WITH_DEFAULT_PARAMS => [
        'CODELINK',
     ],
    },



    # AGILENT

    IMPORT_AGILENT_ARRAYS => {

     IIDREGEXP => '^>(\S+?):(\S+)\s*(.*)$',
     #IIDREGEXP => '^>(\S+):(.+)', #EG HACK

     IFIELDORDER => {
        -name       => 1,
        -array_chip => 0,
        -array      => 0,
        -description => 2,
      },

     ARRAY_PARAMS => {
      'Default' => {
        -vendor => 'AGILENT',
        -format  => 'EXPRESSION',
        -type    => 'OLIGO',
        -class   => 'AGILENT',
        
        -is_probeset_array      => 0,
        -is_linked_array        => 1,
        -has_sense_interrogation  => 0,

      },

      'SurePrint_GPL16709_4x44k' => {
       -name => 'SurePrint_GPL16709_4x44k',
       -vendor => 'AGILENT',
       -format  => 'EXPRESSION',
       -type    => 'OLIGO',
       -class   => 'AGILENT',

        -is_probeset_array      => 0,
        -is_linked_array        => 1,
        -has_sense_interrogation  => 0,

       skip_config => {skip_reps =>1, skip_non_unique_names=>1},
      },

      'SurePrint_GPL7083_4x44k' =>
      {
       -name => 'SurePrint_GPL7083_4x44k',
       -vendor => 'AGILENT',
       -format  => 'EXPRESSION',
       -type    => 'OLIGO',
       -class   => 'AGILENT',

        -is_probeset_array      => 0,
        -is_linked_array        => 1,
        -has_sense_interrogation  => 0,

       skip_config => {skip_reps =>1, skip_non_unique_names=>1},
      },

      'CGH_44b' => {
        -name => 'CGH_44b',
        -vendor => 'AGILENT',
        -format  => 'CGH',
        -type    => 'OLIGO',
        -class   => 'AGILENT',

        -is_probeset_array      => 0,
        -is_linked_array        => 1,
        -has_sense_interrogation  => 0,

      },
     },

     INPUT_FORMAT => 'FASTA',

     ARRAYS_WITH_DEFAULT_PARAMS => [
        'G2518A',
        'G2519F', # This is a generic name!
        'A-MEXP-2203',
        'WholeGenome_4x44k_v1',
        'WholeGenome_4x44k_v2',
        'SurePrint_G3_GE_8x60k',
        'SurePrint_G3_GE_8x60k_v2',
        'WholeGenome_4x44k_v3',
        '012795',
        '015061',
        '020186',
        'GPL13394',
        'GPL14144',
        'GPL8304',
        'GPL14143',
        'GPL8303',
        'GPL14142',
        'GPL13914',
        'GPL19516',
        'GPL14146',
        'GPL14145',
        'GPL14372',
        'AgilentTiling',
        'WholeGenome',
        'OligoArray_012795',
        'AGILENT_059389_Custom_Chicken_GE_8X60k',
	'Arraystar',
	'WholeGenome_4x44k',
	'CHO2agl44v1',
	'037725_HamArrayV',
	'Agilent_8x15K',
        'GPL10157',
        'GPL10158',
        'GPL17465',
        'GPL19384',
        'GPL6848',
        'GPL18606',
        'GPL23389',
        'GPL10248',
        'GPL15189',
        'GPL17689',
        'GPL15204',
        'GPL23307',
        'GPL20908',
        'GPL20928',
        'GPL15217',
        'GPL15190',
        'GPL13795',
        'GPL13723',
        'GPL13522',
        'GPL16032',
        'GPL19564',
        'GPL16776',
        'GPL18520',
        'GPL19666',
        'GPL26966',
        ## Plant arrays
        'G2519F-015058', # Rice
        'G2519F-015059', # Arabidopsis
        'G2519F-015241', # Rice
        'G2519F-016047', # Corn
        'G2519F-016772', # Soybean
        'G2519F-021113', # Tobacco
        'G2519F-021169', # Arabidopsis
        'G2519F-021623', # Barley
        'G2519F-022270', # Tomato
        'G2519F-022297', # Wheat
        'G2519F-022520', # Brassica
        'G2519F-022523', # Cotton
        'G2519F-022524', # Medicago
        'G4136A-011839', # Arabidopsis
        'G4136B-013324', # Arabidopsis
        'G4138A-012106', # Rice
        'G4142A-012600', # Arabidopsis
        'GPL14629',
        'GPL14664',
        'GPL15450',
        'GPL15747',
        'GPL15799',
        'GPL17670',
        'GPL17686',
        'GPL20686',
        'GPL20834',
        'GPL20900',
        'GPL21244',
        'GPL21361',
        'GPL21860',
        'GPL22083',
        'GPL23036',
        'GPL7244',
        'GPL7301',
        'GPL7302',
        'GPL7735',
        'GPL7801',
        'GPL9060',
        'GPL9074',

      ],
    },



    IMPORT_WUSTL_ARRAYS => {
       IIDREGEXP => '^>(\S+):(\S+)',

       IFIELDORDER => {
         -name       => 1,
         -array_chip => 0,
         -array      => 0,
       },

       ARRAY_PARAMS => {
         'Default' => {
           -vendor => 'WUSTL',
           -format  => 'EXPRESSION',
           -type    => 'OLIGO',
           -class   => 'WUSTL',

            -is_probeset_array      => 0,
            -is_linked_array        => 1,
            -has_sense_interrogation  => 0,

         },
       },

       INPUT_FORMAT => 'FASTA',

       ARRAYS_WITH_DEFAULT_PARAMS => [
        'WUSTL-C_elegans',
       ],
    },



    IMPORT_SLRI_ARRAYS => {
       IIDREGEXP => '^>(\S+):(\S+)',

       IFIELDORDER => {
         -name       => 1,
         -array_chip => 0,
         -array      => 0,
       },

       ARRAY_PARAMS => {
         'Default' => {
           -vendor => 'SLRI',
           -format  => 'EXPRESSION',
           -type    => 'OLIGO',
           -class   => 'SLRI',

            -is_probeset_array      => 0,
            -is_linked_array        => 1,
            -has_sense_interrogation  => 0,

         },
       },

       INPUT_FORMAT => 'FASTA',

       ARRAYS_WITH_DEFAULT_PARAMS => [
        'GPL3518',
       ],
     },



    IMPORT_UCSF_ARRAYS =>
     {
       IIDREGEXP => '^>(\S+):(\S+)',#Need to add desc field here

       IFIELDORDER => {
         -name       => 1,
         -array_chip => 0,
         -array      => 0,
       },

       ARRAY_PARAMS => {
         'Default' => {
           -vendor => 'UCSF',
           -format  => 'EXPRESSION',
           -type    => 'OLIGO',
           -class   => 'UCSF',

            -is_probeset_array      => 0,
            -is_linked_array        => 1,
            -has_sense_interrogation  => 0,


         },
       },

       INPUT_FORMAT => 'FASTA',

       ARRAYS_WITH_DEFAULT_PARAMS => [
        'GPL9450',
       ],
     },



    IMPORT_NIMBLEGEN_MODENCODE_ARRAYS =>
     {
       IIDREGEXP => '^>(\S+):(\S+)',

       IFIELDORDER => {
          -name       => 1,
          -array_chip => 0,
          -array      => 0,
       },

       ARRAY_PARAMS => {
         'Default' => {
           -vendor => 'NIMBLEGEN',
           -format  => 'EXPRESSION',
           -type    => 'OLIGO',
           -class   => 'NIMBLEGEN_MODENCODE',

            -is_probeset_array      => 0,
            -is_linked_array        => 1,
            -has_sense_interrogation  => 0,

         },
       },

       INPUT_FORMAT => 'FASTA',

       ARRAYS_WITH_DEFAULT_PARAMS => [
        'GPL9450',
        'Nimblegen_modencode',
        'GPL8673',
       ],
     },



    # NIMBLEGEN

    IMPORT_NIMBLEGEN_ARRAYS =>
    {

       IIDREGEXP => '^>probe:(\S+?):(\w+)\-.*[TranscriptCluster|ProbeSet]ID=(\S+);',

       IFIELDORDER => {
          -name       => 1,
          -array_chip => 0,
          -array      => 0,
          -probe_set   => 2,
       },


       ARRAY_PARAMS => {
          'Default' => {
            -vendor => 'NIMBLEGEN',
            -format => 'EXPRESSION',
            -type    => 'OLIGO',
            -class   => 'NIMBLEGEN',

            -is_probeset_array       => 1,
            -is_linked_array         => 1,
            -has_sense_interrogation => 0,
          },
       },

       INPUT_FORMAT => 'FASTA',

       ARRAYS_WITH_DEFAULT_PARAMS => [
         'NimbleGen_13K',
         'GPL13762',
         'GPL21301',
         'GPL14994',
         'GPL14562',
       ],
    },

    # NIMBLEGEN_ZEBRAFISH

    IMPORT_NIMBLEGEN_ZEBRAFISH_ARRAYS =>
    {

       IIDREGEXP => '^>probe:(\S+?):.*ProbeID=(\S+);.*[TranscriptCluster|ProbeSet]ID=(\S+);',

       IFIELDORDER => {
          -name       => 1,
          -array_chip => 0,
          -array      => 0,
          -probe_set   => 2,
       },


       ARRAY_PARAMS => {
          'Default' => {
            -vendor => 'NIMBLEGEN',
            -format => 'EXPRESSION',
            -type    => 'OLIGO',
            -class   => 'NIMBLEGEN',

            -is_probeset_array       => 1,
            -is_linked_array         => 1,
            -has_sense_interrogation => 0,
          },
       },

       INPUT_FORMAT => 'FASTA',

       ARRAYS_WITH_DEFAULT_PARAMS => [
         'GPL10076',
         'GPL10392',
         'GPL13318',
         'GPL13784',
         'GPL14375',
         'GPL14607',
         'GPL17210',
         'GPL21560',
         'GPL22527',
         'GPL7338',
       ],
    },




    # PHALANX

    # Human
    # ftp://ftp.phalanxbiotech.com/pub/probe_sequences/hoa

    # Mouse
    # ftp://ftp.phalanxbiotech.com/pub/probe_sequences/moa

    IMPORT_PHALANX_ARRAYS => {
     IIDREGEXP => '^>(\S+):(\S+).*$',

     IFIELDORDER => {
        -name       => 1,
        -array_chip => 0,
        -array      => 0,
      },

     ARRAY_PARAMS => {
        'Default' => {
          -vendor => 'PHALANX',
          -format  => 'EXPRESSION',
          -type    => 'OLIGO',
          -class   => 'PHALANX',

          -is_probeset_array      => 0,
          -is_linked_array        => 1,
          -has_sense_interrogation  => 0,

        },
     },

     INPUT_FORMAT => 'FASTA',

     ARRAYS_WITH_DEFAULT_PARAMS => [
      'OneArray',
     ],

    },



    # LEIDEN

    IMPORT_LEIDEN_ARRAYS => {
    
     IIDREGEXP => '^>(\S+):(\S+)\s*(.*)$',

    IFIELDORDER => {
      -name        => 1,
      -array_chip  => 0,
      -array       => 0,
      -description => 2,
    },

    ARRAY_PARAMS => {
      'Default' => {
        -vendor => 'LEIDEN',
        -format  => 'EXPRESSION',
        -type    => 'OLIGO',
        -class   => 'LEIDEN',

        -is_probeset_array      => 0,
        -is_linked_array        => 1,
        -has_sense_interrogation  => 0,

      },
    },

     INPUT_FORMAT => 'FASTA',

     ARRAYS_WITH_DEFAULT_PARAMS => [
      'LEIDEN2',
      'LEIDEN3',
     ],
    },



    # STEMPLE

    IMPORT_STEMPLE_LAB_SANGER_ARRAYS => {
    
      IIDREGEXP => '^>(\S+):(\S+)\s*(.*)$',

      IFIELDORDER => {
        -name       => 1,
        -array_chip => 0,
        -array      => 0,
        -description   => 2,
      },

      ARRAY_PARAMS => {
        'Default' => {
          -vendor => 'STEMPLE_LAB_SANGER',
          -format  => 'EXPRESSION',
          -type    => 'OLIGO',
          -class   => 'STEMPLE_LAB_SANGER',
          
          -is_probeset_array      => 0,
          -is_linked_array        => 1,
          -has_sense_interrogation  => 0,

        },
      },

     INPUT_FORMAT => 'FASTA',

     ARRAYS_WITH_DEFAULT_PARAMS => [
      'MattArray1',
      'MattArray2',
     ],
    },



    IMPORT_CATMA_ARRAYS => {
    
     IIDREGEXP => '^>(\S+):(\S+)',

      IFIELDORDER => {
        -name       => 1,
        -array_chip => 0,
        -array      => 0,
      },

      ARRAY_PARAMS => {
        'Default' => {
          -vendor => 'CATMA',
          -format  => 'EXPRESSION',
          -type    => 'OLIGO',
          -class   => 'CATMA',
        },
      },

     INPUT_FORMAT => 'FASTA',

     ARRAYS_WITH_DEFAULT_PARAMS => [
      'CATMA',
      'CATMA_v5',
     ],
   },



   IMPORT_NSF_ARRAYS => {
     IIDREGEXP => '^>(\S+):(\S+)',

     IFIELDORDER => {
        -name       => 1,
        -array_chip => 0,
        -array      => 0,
      },

     ARRAY_PARAMS => {
        'Default' => {
          -vendor  => 'NSF',
          -format  => 'EXPRESSION',
          -type    => 'OLIGO',
          -class   => 'NSF',
        },
      },

     INPUT_FORMAT => 'FASTA',

     ARRAYS_WITH_DEFAULT_PARAMS => [
      'BGIYale',
      'NSF20K',
      'NSF45K',
     ],
    },
   }
);



sub import {
  my ($callpack) = caller(0);   # Name of the calling package
  my $pack = shift;             # Need to move package off @_

  # Get list of variables supplied, or else everything
  my @vars = @_ ? @_ : keys( %Config );
  return unless @vars;

  # Predeclare global variables in calling package
  eval "package $callpack; use vars qw("
    . join(' ', map { '$'.$_ } @vars) . ")";
  die $@ if $@;


  foreach (@vars) {
    if ( defined $Config{$_} ) {
      no strict 'refs';
	    # Exporter does a similar job to the following
	    # statement, but for function names, not
	    # scalar variables:
	    *{"${callpack}::$_"} = \$Config{ $_ };
    } else {
	    die "Error: Config: $_ not known\n";
    }
  }
}

1;
