#!/usr/bin/env perl

use strict;
use Data::Dumper;
use Bio::EnsEMBL::Registry;
use Getopt::Long;

=head1

mysql $(mysql-ens-reg-prod-1-ensadmin details mysql) mnuhn_alnnew_homo_sapiens_funcgen_86_38 < sql/patch_87_88_c.sql
mysql $(mysql-ens-reg-prod-1-ensadmin details mysql) mnuhn_alnnew_homo_sapiens_funcgen_86_38 < sql/patch_87_88_d.sql
mysql $(mysql-ens-reg-prod-1-ensadmin details mysql) mnuhn_alnnew_homo_sapiens_funcgen_86_38 < sql/patch_87_88_e.sql
mysql $(mysql-ens-reg-prod-1-ensadmin details mysql) mnuhn_alnnew_homo_sapiens_funcgen_86_38 < sql/patch_87_88_f.sql
mysql $(mysql-ens-reg-prod-1-ensadmin details mysql) mnuhn_alnnew_homo_sapiens_funcgen_86_38 < sql/patch_87_88_g.sql
mysql $(mysql-ens-reg-prod-1-ensadmin details mysql) mnuhn_alnnew_homo_sapiens_funcgen_86_38 < sql/patch_87_88_h.sql
mysql $(mysql-ens-reg-prod-1-ensadmin details mysql) mnuhn_alnnew_homo_sapiens_funcgen_86_38 < sql/patch_87_88_i.sql

mysql $(mysql-ens-reg-prod-1-ensadmin details mysql) mnuhn_alnnew_homo_sapiens_funcgen_86_38 -e 'load data local infile "/nfs/nobackup/ensembl/mnuhn/array_mapping/temp/probeset_transcript_assignments.pl" into table probeset_transcript (probeset_id, stable_id);'

mysql $(mysql-ens-reg-prod-1-ensadmin details mysql) mnuhn_alnnew_homo_sapiens_funcgen_86_38 -e 'load data local infile "/nfs/nobackup/ensembl/mnuhn/array_mapping/temp/probe_transcript_assignments.pl" into table probe_transcript (probe_id, stable_id);'

load_probe_to_transcript_assignments.pl \
  --registry /homes/mnuhn/work_dir_probemapping/lib/ensembl-funcgen/registry.pm \
  --species  homo_sapiens \
  --probeset_transcript_assignments /nfs/nobackup/ensembl/mnuhn/array_mapping/temp/probeset_transcript_assignments.pl \
  --probeset_transcript_rejections  /nfs/nobackup/ensembl/mnuhn/array_mapping/temp/probeset_transcript_rejections.pl \
  --probe_transcript_assignments    /nfs/nobackup/ensembl/mnuhn/array_mapping/temp/probe_transcript_assignments.pl \

=cut


my $probeset_transcript_assignments_file;
my $probeset_transcript_rejections_file;
my $probe_transcript_assignments_file;
my $registry;
my $species;

GetOptions (
   'registry=s' => \$registry,
   'species=s'  => \$species,
   'probeset_transcript_assignments_file=s' => \$probeset_transcript_assignments_file,
   'probeset_transcript_rejections_file=s'  => \$probeset_transcript_rejections_file,
   'probe_transcript_assignments_file=s'    => \$probe_transcript_assignments_file,
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
  'truncate probe_transcript;',
  'truncate probeset_transcript;',
  'load data local infile "' . $probeset_transcript_assignments_file . '" into table probeset_transcript (probeset_id, stable_id);',
  'load data local infile "' . $probe_transcript_assignments_file    . '" into table probe_transcript    (probe_id,    stable_id, description);',
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

