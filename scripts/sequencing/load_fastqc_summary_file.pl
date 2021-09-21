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

=head1 load_fastqc_summary_file.pl

=head1 SYNOPSIS

load_fastqc_summary_file.pl       --read_file_id #read_file_id#         --summary_file #fastqc_summary_file#  --work_dir #tempdir# --registry #reg_conf# --species #species#

=head1 DESCRIPTION

=cut

use strict;
use Data::Dumper;
use Getopt::Long;

my $summary_file;
my $read_file_id;
my $work_dir;
my $registry;
my $species;

my %config_hash = (
  "summary_file" => \$summary_file,
  "read_file_id" => \$read_file_id,
  "registry"     => \$registry,
  "work_dir"     => \$work_dir,
  "species"      => \$species,
);

my $result = GetOptions(
  \%config_hash,
  'read_file_id=s',
  'summary_file=s',
  'work_dir=s',
  'species=s',
  'registry=s',
);

die unless(-e $summary_file);
die unless($read_file_id);
die("No registry!") unless($registry);

use Bio::EnsEMBL::Registry;
Bio::EnsEMBL::Registry->load_all($registry, 1, 1, 0, 1);

open IN, $summary_file;

my %parsed;

while (my $current_line = <IN>) {
  chomp $current_line;
  my @f = split "\t", $current_line;
  $parsed{$f[1]} = $f[0];
}

close(IN);

# $VAR1 = {
#           'Sequence Duplication Levels' => 'PASS',
#           'Per base N content' => 'PASS',
#           'Per sequence GC content' => 'PASS',
#           'Overrepresented sequences' => 'WARN',
#           'Per base sequence content' => 'PASS',
#           'Kmer Content' => 'FAIL',
#           'Per sequence quality scores' => 'WARN',
#           'Sequence Length Distribution' => 'PASS',
#           'Per tile sequence quality' => 'FAIL',
#           'Basic Statistics' => 'PASS',
#           'Per base sequence quality' => 'WARN',
#           'Adapter Content' => 'PASS'
#         };

use Bio::EnsEMBL::Funcgen::FastQC;
my $fastqc = Bio::EnsEMBL::Funcgen::FastQC->new(
    -read_file_id                     => $read_file_id,
    -basic_statistics                 => $parsed{'Basic Statistics'},
    -per_base_sequence_quality        => $parsed{'Per base sequence quality'},
    -per_tile_sequence_quality        => $parsed{'Per tile sequence quality'},
    -per_sequence_quality_scores      => $parsed{'Per sequence quality scores'},
    -per_base_sequence_content        => $parsed{'Per base sequence content'},
    -per_sequence_gc_content          => $parsed{'Per sequence GC content'},
    -per_base_n_content               => $parsed{'Per base N content'},
    -sequence_length_distribution     => $parsed{'Sequence Length Distribution'},
    -sequence_duplication_levels      => $parsed{'Sequence Duplication Levels'},
    -overrepresented_sequences        => $parsed{'Overrepresented sequences'},
    -adapter_content                  => $parsed{'Adapter Content'},
    -kmer_content                     => $parsed{'Kmer Content'},
);

my $fastqc_adaptor = Bio::EnsEMBL::Registry->get_adaptor($species, 'funcgen', 'fastqc');

my $existing_entry = $fastqc_adaptor->fetch_by_read_file_id($read_file_id);

if (defined $existing_entry) {
  warn("There already is an entry for this read file.");
  warn("Deleting existing enty");
  $fastqc_adaptor->_delete($existing_entry);
}

$fastqc_adaptor->store($fastqc);





