#!/usr/bin/env perl
=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2025] EMBL-European Bioinformatics Institute

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

perl scripts/sequencing/presort_peak_table.pl \
  --registry /homes/mnuhn/work_dir_regbuild_testrun/lib/ensembl-funcgen/registry.with_previous_version.human_regbuild_testdb12.pm \
  --species homo_sapiens

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
);

my $result = GetOptions(
  \%config_hash,
  'registry=s',
  'species=s',
);

Bio::EnsEMBL::Registry->load_all($registry);

my $logger = Bio::EnsEMBL::Utils::Logger->new();
$logger->init_log;

my $funcgen_adaptor = Bio::EnsEMBL::Registry->get_DBAdaptor($species, 'funcgen');

my $dbc = $funcgen_adaptor->dbc;

my $sth = $dbc->prepare("drop table if exists peak_sorted;");
$sth->execute;
$sth->finish;

# $sth = $dbc->prepare("create table peak_sorted as select * from peak limit 10;");
# $sth->execute;
# $sth->finish;

$sth = $dbc->prepare("show create table peak;");
$sth->execute;

my $x = $sth->fetchrow_hashref;
my $sql = $x->{'Create Table'};

print Dumper($sql);

$sth = $dbc->prepare("rename table peak to peak_backup;");
$sth->execute;

$sth = $dbc->prepare($sql);
$sth->execute;

$sth = $dbc->prepare("
  insert into 
    peak
  select 
    * 
  from 
    peak_backup 
  order by 
    seq_region_id, 
    seq_region_start, 
    peak_calling_id
");
$sth->execute;

# drop table peak;
# rename table peak_backup to peak;

$logger->finish_log;

