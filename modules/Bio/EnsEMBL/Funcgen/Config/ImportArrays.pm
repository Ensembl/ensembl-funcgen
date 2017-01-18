# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016] EMBL-European Bioinformatics Institute
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
It is entirely dependant on the arrays.env environment which can be used
to set up and run the pipeline in an easy and interactive way. This contains
all possible configurations which will then be set dynamically by the RunnableDB
for each instance using the input_id as a key into a separate ImportArrays.conf
file, listed here as ARRAY_FORMAT_FILE.


The layout of the configuration is a set of hashes,
each one keyed by logic name. There is also a DEFAULT hash,
which is used as the default for all logic names (this
was the configuration pattern stolen from Exonerate2Genes,
although in this case it's very unlikely you will need to have
different configs by logic name).

=head1 CONTACT

=cut


#package Bio::EnsEMBL::Analysis::Config::ImportArrays;
package Bio::EnsEMBL::Funcgen::Config::ImportArrays;

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
#      DNADB => {
#       -dbname          => $ENV{'DNADB_NAME'},
#       -host            => $ENV{'DNADB_HOST'},
#       -port            => $ENV{'DNADB_PORT'},
#       -user            => $ENV{'DNADB_USER'},
#       -pass            => $ENV{'DNADB_PASS'},
#       -species         => $ENV{'SPECIES'},
#       -multispecies_db => $ENV{'DNADB_MULTISPECIES_DB'},
#       -species_id      => $ENV{'DNADB_SPECIES_ID'}
#      },
#      OUTPUT_DIR           => $ENV{'WORK_DIR'},

     # Regular expression for parsing file headers
     #
     IIDREGEXP =>  '^>probe:(\S+):(\S+):(\S+:\S+;).*$',

     IFIELDORDER  => {},
     ARRAY_PARAMS => {},
     ARRAYS_WITH_DEFAULT_PARAMS => [],
    },
    IMPORT_AFFY_UTR_ARRAYS => {

      IIDREGEXP => '^>probe:(\S+):(\S+):(\S+:\S+;).*$',

      IFIELDORDER => {
        -name       => 2,
        -array_chip => 0,
        -array      => 0,
        -probe_set  => 1
      },

      ARRAY_PARAMS => {

      Default =>  {
        #-name    => 'Porcine',
        -vendor  => 'AFFY',
        -format  => 'EXPRESSION',
        -type    => 'OLIGO',
        -class   => 'AFFY_UTR',
      },

      platypus_exon => {
        #-name    => 'platypus_exon',
        -vendor  => 'CUSTOM',
        -format  => 'EXPRESSION',
        -type    => 'OLIGO',
        -class   => 'CUSTOM',
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
        'ATH1-121501',
        'Barley1',
        'Rice',
        'Poplar',
        'Vitis_Vinifera',
      ],

    },

    IMPORT_AFFY_ST_ARRAYS => {

     IIDREGEXP => '^>probe:(\S+?):([0-9]+).*[TranscriptCluster|ProbeSet]ID=(\S+);',

     #Can't use ProbeID=([0-9]+) as control probes only have there ProbeID in the concat'd full name string
     #Hence the match will fail.


     IFIELDORDER => {
                     -name       => 1,
                     -array_chip => 0,
                     -array      => 0,
                     -probe_set   => 2,
                    },

     ARRAY_PARAMS =>
     {
       'Default' => {
                     -name => 'RaEx-1_0-st-v1',
                     -vendor => 'AFFY',
                     #-setsize => undef,
                     -format  => 'EXPRESSION',
                     -type    => 'OLIGO',
                     -class   => 'AFFY_ST',
                    },

       'HTA-2_0' => {
                     -name => 'HTA-2_0',
                     -vendor => 'AFFY',
                     -format  => 'EXPRESSION',
                     -type    => 'OLIGO',
                     -description => 'Human Transcriptome Array 2.0',
                     -class   => 'AFFY_ST',
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
     ],
    },

    IMPORT_ILLUMINA_WG_ARRAYS =>
    {
     IIDREGEXP => '^>(\S+):(\S+).*$',

     IFIELDORDER => {
                     -name       => 1,
                     -array_chip => 0,
                     -array      => 0,
                     #-probe_set   => 2,#This could be annotation
                    },

     ARRAY_PARAMS =>
     {

      'Default' => {
                         -name => 'MouseWG_6_V1',
                         -vendor => 'ILLUMINA',
                         #-setsize => undef,
                         -format  => 'EXPRESSION',
                         -type    => 'OLIGO',
                         -class   => 'ILLUMINA_WG',
                        },
     },

     INPUT_FORMAT => 'FASTA',

     ARRAYS_WITH_DEFAULT_PARAMS => [
      'MouseWG_6_V1',
      'MouseWG_6_V2',
      'MouseRef-8_V2',
      'HumanWG_6_V1',
      'HumanHT-12_V3',
      'HumanHT-12_V4',
      'HumanRef-8_V3',
      'HumanWG_6_V2',
      'HumanWG_6_V3',
      'RatRef-12_V1',
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
        #-name    => 'HumanMethylation27',
        -vendor  => 'ILLUMINA',
        -format  => 'METHYLATION',
        -type    => 'OLIGO',
        -class   => 'ILLUMINA_INFINIUM',
      },
    },
    INPUT_FORMAT => 'FASTA',

     ARRAYS_WITH_DEFAULT_PARAMS => [
	    'HumanMethylation27',
	    'HumanMethylation450',
    ],
  },

    #CODELINK

    IMPORT_CODELINK_ARRAYS =>
    {
     IIDREGEXP => '^>(\S+):(\S+).*$',

     IFIELDORDER => {
                     -name       => 1,
                     -array_chip => 0,
                     -array      => 0,

                     #-probe_set   => 2,#This could be annotation
                    },

     ARRAY_PARAMS => {

                      'Default' => {
                                    -name => 'CODELINK',
                                    -vendor => 'CODELINK',
                                    #-setsize => undef,
                                    -format  => 'EXPRESSION',
                                    -type    => 'OLIGO',
                             	     -class   => 'CODELINK',
                                   },
     },

     INPUT_FORMAT => 'FASTA',

     ARRAYS_WITH_DEFAULT_PARAMS => [
	    'CODELINK',
     ],
    },

    #AGILENT

    IMPORT_AGILENT_ARRAYS =>
    {

     IIDREGEXP => '^>(\S+):(\S+)\s*(.*)$',
     #IIDREGEXP => '^>(\S+):(.+)', #EG HACK

     IFIELDORDER => {
                     -name       => 1,
                     -array_chip => 0,
                     -array      => 0,
                     -description => 2,
                     #-probe_set   => 2,#This could be annotation
                    },

     ARRAY_PARAMS =>
     {
      'Default' => {
                   -name => 'G2518A',
                   -vendor => 'AGILENT',
                   #-setsize => undef,
                   -format  => 'EXPRESSION',
                   -type    => 'OLIGO',

                   -class   => 'AGILENT',
                  },

      #Rabbit only

      #This naming was erronoes and non-species specific
      #'SurePrint_G2519F_4x44k' => {
      #                             -name => 'SurePrint_G2519F_4x44k',
      #                             -vendor => 'AGILENT',
      #                            #-setsize => undef,
      #                            -format  => 'EXPRESSION',
      #                            -type    => 'OLIGO',
      #
      #                            -class   => 'AGILENT',
      #                            },

      'SurePrint_GPL16709_4x44k' =>
      {
       -name => 'SurePrint_GPL16709_4x44k',
       -vendor => 'AGILENT',
       #-setsize => undef,
       -format  => 'EXPRESSION',
       -type    => 'OLIGO',
       #-description => '',
       -class   => 'AGILENT',
       skip_config => {skip_reps =>1, skip_non_unique_names=>1},
      },

      'SurePrint_GPL7083_4x44k' =>
      {
       -name => 'SurePrint_GPL7083_4x44k',
       -vendor => 'AGILENT',
       #-setsize => undef,
       -format  => 'EXPRESSION',
       -type    => 'OLIGO',
       #-description => '',
       -class   => 'AGILENT',
       skip_config => {skip_reps =>1, skip_non_unique_names=>1},
      },

      'CGH_44b' => {
                    -name => 'CGH_44b',
                    -vendor => 'AGILENT',
                    #-setsize => undef,
                    -format  => 'CGH',
                    -type    => 'OLIGO',
                    -class   => 'AGILENT',
                   },
     },

     INPUT_FORMAT => 'FASTA',

     ARRAYS_WITH_DEFAULT_PARAMS => [
	    'G2518A',
	    'G2519F',
	    'A-MEXP-2203',
	    'G2519F-015059',
	    'G2519F-021169',
	    'G2519F-015241',
	    'G4138A-012106',
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
	    'GPL14146',
	    'GPL14145',
	    'GPL14372',
      ],
    },

    IMPORT_WUSTL_ARRAYS =>
    {
       IIDREGEXP => '^>(\S+):(\S+)',#Need to add desc field here

       IFIELDORDER => {
         -name       => 1,
         -array_chip => 0,
         -array      => 0,
       },

       ARRAY_PARAMS => {

         'Default' => {
           -name => 'WUSTL-C_elegans',
           -vendor => 'WUSTL',
           #-setsize => undef,
           -format  => 'EXPRESSION',
           -type    => 'OLIGO',

           -class   => 'WUSTL',
         },
       },

       INPUT_FORMAT => 'FASTA',

       ARRAYS_WITH_DEFAULT_PARAMS => [
        'WUSTL-C_elegans',
       ],
    },

    IMPORT_SLRI_ARRAYS =>
     {
       IIDREGEXP => '^>(\S+):(\S+)',#Need to add desc field here

       IFIELDORDER => {
         -name       => 1,
         -array_chip => 0,
         -array      => 0,
       },

       ARRAY_PARAMS => {

         'Default' => {
           -name => 'GPL3518',
           -vendor => 'SLRI',
           #-setsize => undef,
           -format  => 'EXPRESSION',
           -type    => 'OLIGO',

           -class   => 'SLRI',
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
           -name => 'GPL9450',
           -vendor => 'UCSF',
           #-setsize => undef,
           -format  => 'EXPRESSION',
           -type    => 'OLIGO',

           -class   => 'UCSF',
         },
       },

       INPUT_FORMAT => 'FASTA',

       ARRAYS_WITH_DEFAULT_PARAMS => [
        'GPL9450',
       ],
     },


    #PHALANX
    #Human
    #ftp://ftp.phalanxbiotech.com/pub/probe_sequences/hoa
    #Mouse
    #ftp://ftp.phalanxbiotech.com/pub/probe_sequences/moa
    IMPORT_PHALANX_ARRAYS =>
    {
     IIDREGEXP => '^>(\S+):(\S+).*$',

     IFIELDORDER => {
                     -name       => 1,
                     -array_chip => 0,
                     -array      => 0,
                     #-probe_set   => 2,#This could be annotation
                    },

     ARRAY_PARAMS => {
                      'Default' => {
                                     -name => 'OneArray',
                                     -vendor => 'PHALANX',
                                     #-setsize => undef,
                                     -format  => 'EXPRESSION',
                                     -type    => 'OLIGO',

                                     -class   => 'PHALANX',
                                    },
                     },

     INPUT_FORMAT => 'FASTA',

     ARRAYS_WITH_DEFAULT_PARAMS => [
      'OneArray',
     ],

    },

    #LEIDEN

    IMPORT_LEIDEN_ARRAYS =>
    {
     IIDREGEXP => '^>(\S+):(\S+)\s*(.*)$',

     IFIELDORDER => {
                     -name        => 1,
                     -array_chip  => 0,
                     -array       => 0,
                     -description => 2,
                     #-probe_set   => 2,#This could be annotation
                    },

     ARRAY_PARAMS => {
                      #Danio
                      'Default' => {
                                    -name => 'LEIDEN2',
                                    -vendor => 'LEIDEN',
                                    #-setsize => undef,
                                    -format  => 'EXPRESSION',
                                    -type    => 'OLIGO',

                                    -class   => 'LEIDEN',
                                   },
                     },

     INPUT_FORMAT => 'FASTA',

     ARRAYS_WITH_DEFAULT_PARAMS => [
      'LEIDEN2',
      'LEIDEN3',
     ],
    },

    #STEMPLE

    IMPORT_STEMPLE_LAB_SANGER_ARRAYS =>
    {
     IIDREGEXP => '^>(\S+):(\S+)\s*(.*)$', #Need to add desc field here

     IFIELDORDER => {
                     -name       => 1,
                     -array_chip => 0,
                     -array      => 0,
                     -description   => 2,
                    },

     ARRAY_PARAMS => {
                      #Danio
                      'Default' => {
                                       -name => 'MattArray1',
                                       -vendor => 'STEMPLE_LAB_SANGER',
                                       #-setsize => undef,
                                       -format  => 'EXPRESSION',
                                       -type    => 'OLIGO',

                                       -class   => 'STEMPLE_LAB_SANGER',
                                      },
                     },

     INPUT_FORMAT => 'FASTA',

     ARRAYS_WITH_DEFAULT_PARAMS => [
      'MattArray1',
      'MattArray2',
     ],
    },

    IMPORT_CATMA_ARRAYS =>
    {
     IIDREGEXP => '^>(\S+):(\S+)', #Need to add desc field here

     IFIELDORDER => {
                     -name       => 1,
                     -array_chip => 0,
                     -array      => 0,
                    },

     ARRAY_PARAMS => {
                      'Default' => {
                                  -name => 'CATMA',
                                  -vendor => 'CATMA',
                                  #-setsize => undef,
                                  -format  => 'EXPRESSION',
                                  -type    => 'OLIGO',

                                  -class   => 'CATMA',
                                 },
                     },
     INPUT_FORMAT => 'FASTA',

     ARRAYS_WITH_DEFAULT_PARAMS => [
      'CATMA',
     ],
   },

   IMPORT_NSF_ARRAYS =>
   {
     IIDREGEXP => '^>(\S+):(\S+)', #Need to add desc field here

     IFIELDORDER => {
                     -name       => 1,
                     -array_chip => 0,
                     -array      => 0,
                    },

     ARRAY_PARAMS => {
                      'Default' => {
                                    -name => 'BGIYale',
                                    -vendor => 'NSF',
                                    #-setsize => undef,
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
