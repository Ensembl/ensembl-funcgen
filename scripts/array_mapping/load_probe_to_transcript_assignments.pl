#!/usr/bin/env perl

use strict;
use Data::Dumper;
use Bio::EnsEMBL::Registry;
use Getopt::Long;

=head1

load_probeset_to_transcript_assignments.pl \
  --registry /homes/mnuhn/work_dir_probemapping/lib/ensembl-funcgen/registry.pm \
  --species  homo_sapiens \
  --array_name foobar \
  --probe_transcript_assignments_file /nfs/nobackup/ensembl/mnuhn/array_mapping/temp/homo_sapiens/probeset_to_transcript_file.pl

=cut


my $probe_transcript_assignments_file;
my $species;
my $registry;
my $array_name;

GetOptions (
   'registry=s'    => \$registry,
   'species=s'     => \$species,
   'array_name=s'  => \$array_name,
   'probe_transcript_assignments_file=s' => \$probe_transcript_assignments_file,
);

use Bio::EnsEMBL::Utils::Logger;
my $logger = Bio::EnsEMBL::Utils::Logger->new();
$logger->init_log;

Bio::EnsEMBL::Registry->load_all($registry);
my $funcgen_db_adaptor = Bio::EnsEMBL::Registry->get_DBAdaptor($species, 'funcgen');
my $funcgen_dbc = $funcgen_db_adaptor->dbc;

my $mysql_base_cmd = 'mysql'
  . ' --host='     . $funcgen_dbc->host 
  . ' --port='     . $funcgen_dbc->port 
  . ' --user='     . $funcgen_dbc->username 
  . ' --password=' . $funcgen_dbc->password 
  . ' '            . $funcgen_dbc->dbname
  . ' -e '
;

my @load_command = map { $mysql_base_cmd . "'" . $_ . "'" } (
  'load data local infile "' . $probe_transcript_assignments_file . '" into table probe_transcript (probe_id, stable_id, description);',
);

foreach my $current_load_command (@load_command) {

  $logger->info("Running:\n");
  $logger->info("$current_load_command\n");
  
  my $exit_code = system($current_load_command);
  
  if ($exit_code != 0) {
    $logger->error("Failure when running command\n$current_load_command\n");
  }
}

$logger->info("Done.\n");
$logger->finish_log;
