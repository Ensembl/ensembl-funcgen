#!/usr/bin/perl


use strict;
use warnings;
use feature qw(say);
use Config::Tiny;


use Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Utils::SqlHelper;
use Bio::EnsEMBL::Funcgen::Utils::EFGUtils  qw(dump_data);
use Bio::EnsEMBL::Registry;

main();

sub main {

  my $cfg = Config::Tiny->new;
  $cfg = Config::Tiny->read('db_test.ini') or die "Could not find cfg";
  my $db = connect_adaptor($cfg);

  my $rfa   = $db->get_RegulatoryFeatureAdaptor;
  my $epga  = $db->get_EpigenomeAdaptor;
  my $rf    = $rfa->fetch_by_stable_id('ENSR00001381357');
  say $rf->feature_type->name;
  my $hash = {};
  for my $ra (@{$rf->regulatory_activity}) {
    $hash->{activity}->{$ra->activity}++;
     $hash->{epigenomes}->{$ra->epigenome->name}++;
  }
  $hash->{epigenome_count} = $rf->epigenome_count;
  $hash->{feature_type_name} = $rf->feature_type->name;
  say dump_data($hash,1,1);
  #say dump_data($rf->regulatory_activity->[0],1,1);
}
sub connect_adaptor {
  my ($cfg) = @_;

  my $db = Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor->new (
    -user       => $cfg->{efg_db}->{user},
    -pass       => $cfg->{efg_db}->{pass},
    -host       => $cfg->{efg_db}->{host},
    -port       => $cfg->{efg_db}->{port},
    -dbname     => $cfg->{efg_db}->{dbname},
    -dnadb_user => $cfg->{dna_db}->{user},
    -dnadb_name => $cfg->{dna_db}->{dbname},
    -dnadb_host => $cfg->{dna_db}->{host},
    );

  $db->dbc->do("SET sql_mode='traditional'");
  say $db->dbc->host;
  say $db->dbc->dbname;

  return $db;
}

sub connect_registry {
 my $hash = {};
  my $registry = 'Bio::EnsEMBL::Registry';

  $registry->load_registry_from_db(
  -host => 'ensembldb.ensembl.org',
  -user => 'anonymous',
#  -port   => 3337,

  );

  $hash->{slice} = $registry->get_adaptor('human', 'core', 'slice');
  $hash->{slice}->dbc->disconnect_if_idle;

  return $hash;
}
