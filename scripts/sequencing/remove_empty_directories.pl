#!/usr/bin/env perl
=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2022] EMBL-European Bioinformatics Institute

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

perl scripts/sequencing/remove_empty_directories.pl \
  --registry /homes/mnuhn/work_dir_ersa/lib/ensembl-funcgen/registry.pm \
  --species homo_sapiens \
  --dry_run 1 \
  --data_root_dir /hps/nobackup/production/ensembl/mnuhn/chip_seq_analysis/dbfiles

=cut

use strict;
use Data::Dumper;
use Getopt::Long;
use Bio::EnsEMBL::DBSQL::DBConnection;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Utils::Logger;

use Bio::EnsEMBL::Funcgen::Utils::GoodUtils qw( create_species_assembly_path );

my $registry;
my $species;
my $data_root_dir;
my $dry_run;

my %config_hash = (
  'registry'      => \$registry,
  'species'       => \$species,
  'data_root_dir' => \$data_root_dir,
  'dry_run'       => \$dry_run,
);

my $result = GetOptions(
  \%config_hash,
  'registry=s',
  'species=s',
  'data_root_dir=s',
  'dry_run=s',
);

Bio::EnsEMBL::Registry->load_all($registry);
my $species_assembly_path = create_species_assembly_path($species);

use Cwd 'abs_path';
my $path_to_files = abs_path($data_root_dir . '/' . $species_assembly_path);

my $logger = Bio::EnsEMBL::Utils::Logger->new();
$logger->init_log;

use File::Find;

finddepth(
  sub {
    my $directory = $File::Find::name;
    
    if (! -d $directory) {
      return;
    }
    
    if (folder_is_empty($directory)) {
      $logger->info("Empty: $directory\n");
      rmdir($directory);
    }
    return;
  },  $path_to_files
);

$logger->finish_log;

sub folder_is_empty {
    my $dirname = shift;
    opendir(my $dh, $dirname) or die "Not a directory";
    my @stuff = grep { $_ ne "." && $_ ne ".." } readdir($dh);
    
    my $num_files = scalar(@stuff);
    #print " $num_files ";
    return $num_files == 0;
}



