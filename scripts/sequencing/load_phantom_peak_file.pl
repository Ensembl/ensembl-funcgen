#!/usr/bin/env perl
=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2021] EMBL-European Bioinformatics Institute

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

=head1 load_phantom_peak_file.pl

=head1 SYNOPSIS

=head1 DESCRIPTION

=cut

use strict;
use Data::Dumper;
use Getopt::Long;
use Bio::EnsEMBL::Utils::Logger;

my $result_file;
my $alignment_name;
my $registry;
my $species;
my $failed;

my %config_hash = (
  "result_file"     => \$result_file,
  "alignment_name"  => \$alignment_name,
  "registry"        => \$registry,
  "species"         => \$species,
  "failed"          => \$failed,
);

my $result = GetOptions(
  \%config_hash,
  'alignment_name=s',
  'result_file=s',
  'registry=s',
  'species=s',
  'failed=s',
);

my $error_message = undef;

if (! $result_file) {
  die("The result_file parameter was not specified!");
}
if (! -e $result_file) {
  $error_message = "The result_file ($result_file) specified on the command line does not exist!";
  $failed = 1;
}
if (! $alignment_name) {
  die("The alignment_name parameter was not specified!");
}

my $logic_name = 'phantom peak quality tools';
my @phantom_peak_analysis_details = (
        -logic_name      => $logic_name,
        -program         => 'run_spp.R',
        -parameters      => '',
        -description     => 'Computes enrichment and quality measures and fragment lengths for ChIP-seq/DNase-seq/FAIRE-seq/MNase-seq data',
        -display_label   => 'phantom peak quality tools',
        -displayable     => undef,
);

my $logger = Bio::EnsEMBL::Utils::Logger->new();
$logger->init_log;

use Bio::EnsEMBL::Registry;
Bio::EnsEMBL::Registry->load_all($registry, 1, 1, 0, 1);

my $phantom_peak_adaptor = Bio::EnsEMBL::Registry->get_adaptor($species, 'funcgen', 'PhantomPeak');

my $analysis_adaptor = Bio::EnsEMBL::Registry->get_adaptor($species, 'funcgen', 'analysis');
my $analysis = $analysis_adaptor->fetch_by_logic_name($logic_name);

my $alignment_adaptor = Bio::EnsEMBL::Registry->get_adaptor($species, 'funcgen', 'alignment');
my $alignment = $alignment_adaptor->fetch_by_name($alignment_name);
my $alignment_id = $alignment->dbID;

if (! $analysis ) {
      $logger->info("No analysis with logic name $logic_name found. Creating one.");
      $analysis = Bio::EnsEMBL::Analysis->new(@phantom_peak_analysis_details);
      $analysis_adaptor->store($analysis);
}
my $analysis_id = $analysis->dbID;

if ($failed) {

  use Bio::EnsEMBL::Funcgen::PhantomPeak;
  my $phantom_peak = Bio::EnsEMBL::Funcgen::PhantomPeak->new(
      -alignment_id  => $alignment_id,
      -run_failed    => 1,
      -error_message => undef,
      -analysis_id   => $analysis_id,
  );
  $phantom_peak_adaptor->store($phantom_peak);
  exit;
}


open IN, $result_file;

while (my $current_line = <IN>) {
  chomp $current_line;

  my @f = split "\t", $current_line;

  (
    my $filename,
    my $numReads,
    my $estFragLenTriple,
    my $corr_estFragLenTriple,
    my $phantomPeak,
    my $corr_phantomPeak,
    my $argmin_corr,
    my $min_corr,
    my $NSC,
    my $RSC,
    my $QualityTag,
  ) = @f;
  
  (
    my $estFragLen,
    my $estFragLen2,
    my $estFragLen3,
  ) = split ',', $estFragLenTriple;
  
  if ($estFragLen2 eq '') {
    $estFragLen2 = undef;
  }
  if ($estFragLen3 eq '') {
    $estFragLen3 = undef;
  }
  
  (
    my $corr_estFragLen,
    my $corr_estFragLen2,
    my $corr_estFragLen3,
  ) = split ',', $corr_estFragLenTriple;
  
  if ($corr_estFragLen2 eq '') {
    $corr_estFragLen2 = undef;
  }
  if ($corr_estFragLen3 eq '') {
    $corr_estFragLen3 = undef;
  }

  # This can happen, if the data is very bad. In that case the quality tag 
  # will indicate that with
  #
  if ($RSC eq 'Inf') {
    $RSC = undef;
    die unless ($QualityTag == 2);
  }

  use Bio::EnsEMBL::Funcgen::PhantomPeak;
  my $phantom_peak = Bio::EnsEMBL::Funcgen::PhantomPeak->new(
      -analysis_id         => $analysis_id,
      -alignment_id        => $alignment_id,
      -num_reads           => $numReads,
      -est_frag_len        => $estFragLen,
      -est_frag_len_2      => $estFragLen2,
      -est_frag_len_3      => $estFragLen3,
      -corr_est_frag_len   => $corr_estFragLen,
      -corr_est_frag_len_2 => $corr_estFragLen2,
      -corr_est_frag_len_3 => $corr_estFragLen3,
      -phantom_peak        => $phantomPeak,
      -corr_phantom_peak   => $corr_phantomPeak,
      -argmin_corr         => $argmin_corr,
      -min_corr            => $min_corr,
      -nsc                 => $NSC,
      -rsc                 => $RSC,
      -quality_tag         => $QualityTag,
      -run_failed          => 0,
      -error_message       => undef,

  );
  eval {
    $phantom_peak_adaptor->store($phantom_peak);
  };
  if ($@) {
    my $error_message = $@;
    my $already_exists = $error_message =~ /alignment_id_unique/;
    
    if (!$already_exists) {
      die($error_message);
    }
  }


}

