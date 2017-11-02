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
my $core_db_adaptor    = Bio::EnsEMBL::Registry->get_DBAdaptor($species, 'core');

my $slice_adaptor = $core_db_adaptor->get_SliceAdaptor;
# Include non-reference seq regions, important for mouse.
my $slices = $slice_adaptor->fetch_all('toplevel', undef, 1);

my %seq_region_id_to_name_lookup;

use Hash::Util qw( lock_keys );

for my $slice (@$slices) {
    $seq_region_id_to_name_lookup{$slice->get_seq_region_id} = $slice->seq_region_name;
}

lock_keys(%seq_region_id_to_name_lookup);

my $sql = '
    select
        seq_region_id, 
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
        probe_feature 
        join probe using(probe_id)
        left join probe_set using(probe_set_id)
        join array_chip on(array_chip.array_chip_id=probe.array_chip_id)
        join array using(array_id)
  ';

# Prevent memory issues from buffering
$funcgen_db_adaptor->dbc->db_handle->{mysql_use_result} = 1;

my $helper = $funcgen_db_adaptor->dbc->sql_helper;
open my $fh, '>', $file || die("Can't open $file for writing!");
  
$logger->info("Running $sql\n");

my $last_line_printed;

$helper->execute_no_return(
  -SQL      => $sql,
  -CALLBACK => sub {
    my $row  = shift;
    
    $row->[0] = $seq_region_id_to_name_lookup{$row->[0]};
    
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

