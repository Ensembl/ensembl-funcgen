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

=head1 load_argenrich_file.pl

=head1 SYNOPSIS

load_argenrich_qc_file.pl \
    --argenrich_file        /hps/nobackup/production/ensembl/mnuhn/chip_seq_analysis/temp_dir/qc_chance/GM18526_NFKB_ChIP-Seq_ENCODE86/argenrich_outfile.txt \
    --user   ensadmin    \
    --pass   ensembl    \
    --port   4545    \
    --host   mysql-ens-reg-prod-2.ebi.ac.uk    \
    --dbname mnuhn_testdb2_homo_sapiens_funcgen_91_38    \
    --work_dir /hps/nobackup/production/ensembl/mnuhn/chip_seq_analysis/temp_dir/qc_chance/GM18526_NFKB_ChIP-Seq_ENCODE86 \
    --experiment_name GM18526_NFKB_ChIP-Seq_ENCODE86

=head1 DESCRIPTION

=cut

use strict;
use Data::Dumper;
use Getopt::Long;
use Bio::EnsEMBL::DBSQL::DBConnection;
use Bio::EnsEMBL::Utils::Logger;
use Bio::EnsEMBL::DBSQL::DBAdaptor;

my $argenrich_file;
my $dry_run;
my $user;
my $pass;
my $host;
my $port;
my $dbname;
my $work_dir;
my $signal_alignment_name;
my $control_alignment_name;
my $failed;
my $species;

my %config_hash = (
  'argenrich_file'  => \$argenrich_file,
  'dry_run'         => \$dry_run,
  'user'            => \$user,
  'pass'            => \$pass,
  'port'            => \$port,
  'host'            => \$host,
  'dbname'          => \$dbname,
  'work_dir'        => \$work_dir,
  'signal'          => \$signal_alignment_name,
  'control'         => \$control_alignment_name,
  'failed'          => \$failed,
  'species'         => \$species,
);

my $result = GetOptions(
  \%config_hash,
  'argenrich_file=s',
  'dry_run',
  'user=s',
  'pass=s',
  'port=s',
  'host=s',
  'dbname=s',
  'work_dir=s',
  'signal=s',
  'control=s',
  'failed=s',
  'species=s',
);


if (! -e $argenrich_file && !$failed) {
  die;
}

my $logger = Bio::EnsEMBL::Utils::Logger->new();
$logger->init_log;

my @tracking_db_connection_details = (
  -user     => $user,
  -pass     => $pass,
  -host     => $host,
  -port     => $port,
  -dbname   => $dbname,
);
my $dbc = Bio::EnsEMBL::DBSQL::DBConnection->new(@tracking_db_connection_details);

use Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor;
my $dba = Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor->new(
    -dbconn => $dbc,
);

my $alignment_adaptor = $dba->get_adaptor('Alignment');

my $signal_alignment  = $alignment_adaptor->fetch_by_name($signal_alignment_name);
my $control_alignment = $alignment_adaptor->fetch_by_name($control_alignment_name);
  
if (! defined $signal_alignment) {
    die("Can't fetch signal alignment!");
}
if (! defined $control_alignment) {
    die("Can't fetch control alignment!");
}

my $signal_alignment_id  = $signal_alignment->dbID;
my $control_alignment_id = $control_alignment->dbID;

my $analysis_id = create_analysis_if_not_exists($dbc);

if ($failed) {
  my $chance = Bio::EnsEMBL::Funcgen::Chance->new(
    -signal_alignment_id   => $signal_alignment_id,
    -control_alignment_id  => $control_alignment_id,
    -analysis_id           => $analysis_id,
    -run_failed            => 1,
    -error_message         => undef,
  );

  my $chance_adaptor = Bio::EnsEMBL::Registry->get_adaptor($species, 'funcgen', 'chance');

  eval {
    $chance_adaptor->store($chance);
  };
  exit;
}

open IN, $argenrich_file;

my %key_value_pairs;
LINE: while (my $current_line = <IN>) {
  chomp $current_line;
  #print " - $current_line\n";
  my @f = split '=', $current_line;
  next LINE unless(@f == 2);
  
  $f[1] = 'null' if ($f[1] eq 'NA');
  
  $key_value_pairs{$f[0]} = $f[1];
}
use Hash::Util qw( lock_hash );
lock_hash(%key_value_pairs);

use Bio::EnsEMBL::Funcgen::Chance;

my $chance = Bio::EnsEMBL::Funcgen::Chance->new(
  -signal_alignment_id   => $signal_alignment_id,
  -control_alignment_id  => $control_alignment_id,
  -analysis_id           => $analysis_id,
  -p                     => $key_value_pairs{'p'}, 
  '-q'                   => $key_value_pairs{'q'}, 
  -divergence            => $key_value_pairs{'divergence'}, 
  -z_score               => $key_value_pairs{'z_score'}, 
  -percent_genome_enriched => $key_value_pairs{'percent_genome_enriched'}, 
  -input_scaling_factor    => $key_value_pairs{'input_scaling_factor'}, 
  -differential_percentage_enrichment            => $key_value_pairs{'differential_percentage_enrichment'},
  -control_enrichment_stronger_than_chip_at_bin  => $key_value_pairs{'Control enrichment stronger than ChIP at bin'},
  -first_nonzero_bin_at                          => $key_value_pairs{'Zero-enriched IP, maximum difference at bin'},
  -pcr_amplification_bias_in_Input_coverage_of_1_percent_of_genome => $key_value_pairs{'PCR amplification bias in Input, coverage of 1% of genome'},
  -run_failed     => 0,
  -error_message => undef,
);

my $chance_adaptor = Bio::EnsEMBL::Registry->get_adaptor($species, 'funcgen', 'chance');

eval {
  $chance_adaptor->store($chance);
};

if ($@) {
  my $error_message = $@;
  my $already_exists = $error_message =~ /signal_control_alignment_unique/;
  
  if (!$already_exists) {
    die($error_message);
  }
}

# my $sql = qq(insert into chance (
#       signal_alignment_id, control_alignment_id, analysis_id, p, q, divergence, z_score, percent_genome_enriched, input_scaling_factor, differential_percentage_enrichment,
#       control_enrichment_stronger_than_chip_at_bin,
#       first_nonzero_bin_at,
#       pcr_amplification_bias_in_Input_coverage_of_1_percent_of_genome, path
#     ) values (
#     $signal_alignment_id,
#     $control_alignment_id,
#     $analysis_id,
#     $key_value_pairs{'p'}, 
#     $key_value_pairs{'q'}, 
#     $key_value_pairs{'divergence'}, 
#     $key_value_pairs{'z_score'}, 
#     $key_value_pairs{'percent_genome_enriched'}, 
#     $key_value_pairs{'input_scaling_factor'}, 
#     $key_value_pairs{'differential_percentage_enrichment'},
#     $key_value_pairs{'Control enrichment stronger than ChIP at bin'},
#     $key_value_pairs{'Zero-enriched IP, maximum difference at bin'},
#     $key_value_pairs{'PCR amplification bias in Input, coverage of 1% of genome'},
#     '$work_dir'
#   )
# );
# 
# my $sql_processor;
# if ($dry_run) {
#   $sql_processor = sub {
#     my $sql = shift;
#     $logger->info("Dry run: " . $sql . "\n");
#   };
# } else {
#   $sql_processor = sub {
#     my $sql = shift;
#     $logger->info("Running: " . $sql . "\n");
#      $dbc->do($sql);
#   };
# }
# 
# eval {
#   $sql_processor->($sql);
# };
# 
# if ($@) {
#   my $error_message = $@;
#   my $already_exists = $error_message =~ /signal_control_alignment_unique/;
#   
#   if (!$already_exists) {
#     die($error_message);
#   }
# }

$logger->finish_log;
exit;

sub create_analysis_if_not_exists {

  my $dbc = shift;
  my $logic_name = 'QC Chance';
  my @phantom_peak_analysis_details = (
	  -logic_name      => $logic_name,
	  -program         => 'argenrich.R',
	  -parameters      => '',
	  -description     => '',
	  -display_label   => 'Chance',
	  -displayable     => undef,
  );
  
  my $dba = Bio::EnsEMBL::DBSQL::DBAdaptor->new(
    -DBCONN => $dbc,  
  );
  my $analysis_adaptor = $dba->get_AnalysisAdaptor();
  my $analysis = $analysis_adaptor->fetch_by_logic_name($logic_name);

  if (! $analysis && ! $dry_run) {
	$logger->info("No analysis with logic name $logic_name found. Creating one.");
	$analysis = Bio::EnsEMBL::Analysis->new(@phantom_peak_analysis_details);
	$analysis_adaptor->store($analysis);
  }
  my $analysis_id = $analysis->dbID;
  return $analysis_id;
}
