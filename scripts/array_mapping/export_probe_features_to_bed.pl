#!/usr/bin/env perl

use strict;
use Data::Dumper;
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor;
use Getopt::Long;

=head1

date
export_probe_features_to_bed.pl \
  --registry /homes/mnuhn/work_dir_probemapping/lib/ensembl-funcgen/registry.pm \
  --species  homo_sapiens \
  --file     /nfs/nobackup/ensembl/mnuhn/probe2transcript/probe_features_nogroup.test.bed \
date

=cut

my $registry;
my $species;
my $file;

GetOptions (
   'registry=s' => \$registry,
   'species=s'  => \$species,
   'file=s'     => \$file,
);

use Bio::EnsEMBL::Utils::Logger;
my $logger = Bio::EnsEMBL::Utils::Logger->new();
$logger->init_log;

Bio::EnsEMBL::Registry->load_all($registry);
my $funcgen_db_adaptor = Bio::EnsEMBL::Registry->get_DBAdaptor($species, 'funcgen');

my $sql = '
    select
      seq_region.name, 
      seq_region_start, 
      seq_region_end, 
      seq_region_strand, 
      cigar_line, 
      mismatches, 
      probe_feature_id, 
      probe_id, 
      probe.
      name, 
      probe_set_id, 
      probe_set.name,
      array.name,
      array.vendor,
      array.class,
      is_probeset_array,
      is_linked_array,
      has_sense_interrogation
    from
      probe_feature join probe using(probe_id)
      left join probe_set using(probe_set_id)
      join array_chip using(array_chip_id)
      join array using(array_id)
      join seq_region using(seq_region_id)
  ';

# Prevent memory issues from buffering
$funcgen_db_adaptor->dbc->db_handle->{mysql_use_result} = 1;

my $helper = $funcgen_db_adaptor->dbc->sql_helper;
open my $fh, '>', $file;
  
$logger->info("Running $sql\n");

my $last_line_printed;

$helper->execute_no_return(
  -SQL      => $sql,
  -CALLBACK => sub {
    my $row  = shift;
    
    my $line = join "\t", @$row;
    
    if ($line ne $last_line_printed) {
      $fh->print( $line );
      $fh->print( "\n"  );
      $last_line_printed = $line;
    }
  },
);

$logger->info("Done.\n");
$logger->finish_log;

