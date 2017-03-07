#!/usr/bin/env perl

use strict;
use Data::Dumper;
use Bio::EnsEMBL::Registry;
use Getopt::Long;

=head1

create_bedtools_genome_file.pl \
  --registry /homes/mnuhn/work_dir_probemapping/lib/ensembl-funcgen/registry.pm \
  --species  homo_sapiens \
  --bedtools_genome_file bedtools_genome_file

=cut

my $bedtools_genome_file;
my $registry;
my $species;

GetOptions (
   'registry=s'             => \$registry,
   'species=s'              => \$species,
   'bedtools_genome_file=s' => \$bedtools_genome_file,
);

use Bio::EnsEMBL::Utils::Logger;
my $logger = Bio::EnsEMBL::Utils::Logger->new();
$logger->init_log;

Bio::EnsEMBL::Registry->load_all($registry);
my $funcgen_db_adaptor = Bio::EnsEMBL::Registry->get_DBAdaptor($species, 'core');
my $funcgen_dbc = $funcgen_db_adaptor->dbc;

my $mysql_base_cmd = 'mysql'
  . ' -N'
  . ' --host='     . $funcgen_dbc->host 
  . ' --port='     . $funcgen_dbc->port 
  . ' --user='     . $funcgen_dbc->username 
  . ' --password=' . $funcgen_dbc->password 
  . ' '            . $funcgen_dbc->dbname
  . ' -e '
;

my $sql = q(select seq_region.name, seq_region.length from seq_region join seq_region_attrib using (seq_region_id) join attrib_type using (attrib_type_id) where code=\"toplevel\" group by seq_region.name, seq_region.length order by seq_region.name);

my $cmd = $mysql_base_cmd . "'" . $sql . "'";

my $piped_command = qq(bash -o pipefail -c "$cmd | sort -k1,1 > $bedtools_genome_file");

$logger->info("Running:\n");
$logger->info("$piped_command\n");

my $exit_code = system($piped_command);

if ($exit_code != 0) {
  $logger->error("Failure when running command\n$cmd\n");
}

$logger->info("Done.\n");
$logger->finish_log;

