#!/usr/bin/env perl
=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2017] EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either 405ress or implied.
See the License for the specific language governing permissions and
limitations under the License.

=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

  Questions may also be sent to the Ensembl help desk at
  <http://www.ensembl.org/Help/Contact>.

=head1 trim_peaks_to_seq_region_boundaries.pl

=head1 SYNOPSIS

=head1 DESCRIPTION

perl scripts/sequencing/describe_computation.pl \
  --registry /homes/mnuhn/work_dir_ersa/lib/ensembl-funcgen/registry.pm \
  --species homo_sapiens \
  --peak_calling HeLa_S3_Max_SWEmbl_default_ENCODE

perl scripts/sequencing/describe_computation.pl \
  --registry /homes/mnuhn/work_dir_ersa/lib/ensembl-funcgen/registry.pm \
  --species homo_sapiens \
  --peak_calling Monocytes_CD14___PB__Roadmap_H3K27me3_ccat_histone_Roadmap_Epigenomics

peak_callings="
  GM12878_Jund_SWEmbl_default_ENCODE
  GM12878_PolIII_SWEmbl_default_ENCODE
  H1ESC_H4K5ac_SWEmbl_default_Roadmap_Epigenomics
  H1ESC_H3K18ac_SWEmbl_default_Roadmap_Epigenomics
  GM18526_NFKB_SWEmbl_default_ENCODE
  GM12892_NFKB_SWEmbl_default_ENCODE
  GM19099_NFKB_SWEmbl_default_ENCODE
  GM12878_H3K79me2_ccat_histone_ENCODE
  GM12878_Tr4_SWEmbl_default_ENCODE
  H1ESC_H3K36me3_ccat_histone_ENCODE
  H1ESC_H4K20me1_ccat_histone_Roadmap_Epigenomics
  H1ESC_H3K27me3_ccat_histone_ENCODE
  Thymus_H3K4me1_ccat_histone_Roadmap_Epigenomics
  H1ESC_H4K20me1_ccat_histone_ENCODE
  Monocytes_CD14___PB__Roadmap_H3K9me3_ccat_histone_Roadmap_Epigenomics
  Monocytes_CD14___PB__Roadmap_H3K27me3_ccat_histone_Roadmap_Epigenomics
  H1ESC_H3K4me3_SWEmbl_default_ENCODE
  Monocytes_CD14___PB__Roadmap_H3K36me3_ccat_histone_Roadmap_Epigenomics
  HeLa_S3_Cfos_SWEmbl_default_ENCODE
  Monocytes_CD14___PB__Roadmap_H3K4me1_ccat_histone_Roadmap_Epigenomics
  IMR90_H3K4me3_SWEmbl_default_Roadmap_Epigenomics
  H1ESC_H3K27ac_SWEmbl_default_ENCODE
  Monocytes_CD14___PB__Roadmap_H3K4me3_SWEmbl_default_Roadmap_Epigenomic
  HeLa_S3_Cmyc_SWEmbl_default_ENCODE
  HeLa_S3_E2F6_SWEmbl_default_ENCODE
  GM18951_NFKB_SWEmbl_default_ENCODE
  H1ESC_H3K27ac_SWEmbl_default_Roadmap_Epigenomics
  H1ESC_H3K4me2_SWEmbl_default_ENCODE
  H1ESC_H3K9ac_SWEmbl_default_ENCODE
  GM18505_NFKB_SWEmbl_default_ENCODE
  Right_Atrium_H3K4me1_ccat_histone_Roadmap_Epigenomics
  H1ESC_H3K4me2_SWEmbl_default_Roadmap_Epigenomics
  GM15510_NFKB_SWEmbl_default_ENCODE
  IMR90_H3K4me1_ccat_histone_Roadmap_Epigenomics
  IMR90_H3K4ac_SWEmbl_default_Roadmap_Epigenomics
  Thymus_H3K27ac_SWEmbl_default_Roadmap_Epigenomics
  H1ESC_H3K9ac_SWEmbl_default_Roadmap_Epigenomics
  IMR90_H3K79me2_ccat_histone_Roadmap_Epigenomics
  H1ESC_H2BK15ac_SWEmbl_default_Roadmap_Epigenomics
  IMR90_H3K27me3_ccat_histone_Roadmap_Epigenomics
  IMR90_H3K23ac_SWEmbl_default_Roadmap_Epigenomics
  H1ESC_H3K4me1_ccat_histone_Roadmap_Epigenomics
  IMR90_H3K36me3_ccat_histone_Roadmap_Epigenomics
  IMR90_H3K27ac_SWEmbl_default_Roadmap_Epigenomics
"

for peak_calling in $peak_callings
do
  perl scripts/sequencing/describe_computation.pl \
    --registry /homes/mnuhn/work_dir_ersa/lib/ensembl-funcgen/registry.pm \
    --species homo_sapiens \
    --peak_calling $peak_calling
done 2>/dev/null | less


=cut

use strict;
use Data::Dumper;
use Getopt::Long;
use Bio::EnsEMBL::DBSQL::DBConnection;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Utils::Logger;

my $registry;
my $species;
my $peak_calling_name;

my %config_hash = (
  'registry'      => \$registry,
  'species'       => \$species,
  'peak_calling'  => \$peak_calling_name,
);

my $result = GetOptions(
  \%config_hash,
  'registry=s',
  'species=s',
  'peak_calling=s',
);

Bio::EnsEMBL::Registry->load_all($registry);

my $logger = Bio::EnsEMBL::Utils::Logger->new();
$logger->init_log;

my $peak_calling_adaptor = Bio::EnsEMBL::Registry->get_adaptor($species, 'funcgen', 'peakcalling');

$logger->info("Fetching $peak_calling_name\n");

my $peak_calling = $peak_calling_adaptor->fetch_by_name($peak_calling_name);

if (! defined $peak_calling) {
  $logger->error("Can't find peak calling $peak_calling_name!");
  $logger->error(Dumper($peak_calling_adaptor));
  
  die;
}

use Template;
my $tt = Template->new;

use Bio::EnsEMBL::Funcgen::Template::PeakCallingDescription qw( 
  PEAK_CALLING_TXT_TEMPLATE
);

my $description_template = PEAK_CALLING_TXT_TEMPLATE;

use Number::Format qw( :all );
use File::Spec;

$tt->process(
  \$description_template, 
  {
    peak_calling  => $peak_calling,
    
    canonpath => sub {
      my $path = shift;
      return File::Spec->canonpath($path)
    },
    
    bool_to_yes_no => sub {
      my $boolean = shift;
      if ($boolean) {
        return 'yes'
      }
      return 'no'
    },
    
    round_percent => sub {
      my $number = shift;
      return sprintf("%.2f", $number) . '%';
    },
    default_round => sub {
      my $number = shift;
      return sprintf("%.2f", $number);
    },
    scientific_notation => sub {
      my $number = shift;
      return sprintf("%.2e", $number);
    },
    format_number => sub {
      my $number = shift;
      return format_number($number);
    }
  }
)
    || die $tt->error;


$logger->finish_log;


