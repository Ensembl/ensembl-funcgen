#!/usr/bin/env perl

use strict;
use Data::Dumper;
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor;
use Getopt::Long;

=head1

calculate_arrays_per_object.pl \
  --registry /homes/mnuhn/work_dir_probemapping/lib/ensembl-funcgen/registry.pm \
  --species homo_sapiens \
  --arrays_per_object_file /nfs/nobackup/ensembl/mnuhn/probe2transcript/arrays_per_object.pl \
  --probeset_sizes_file    /nfs/nobackup/ensembl/mnuhn/probe2transcript/probeset_sizes.pl    \
  --object_names_file      /nfs/nobackup/ensembl/mnuhn/probe2transcript/object_names.pl      \

=cut

my $debug;

my $registry;
my $species;
my $arrays_per_object_file;
my $probeset_sizes_file;
my $object_names_file;

GetOptions (
   'registry=s'                => \$registry,
   'species=s'                 => \$species,
   'arrays_per_object_file=s'  => \$arrays_per_object_file,
   'probeset_sizes_file=s'     => \$probeset_sizes_file,
   'object_names_file=s'       => \$object_names_file,
);

use Bio::EnsEMBL::Utils::Logger;
my $logger = Bio::EnsEMBL::Utils::Logger->new();
$logger->init_log;

Bio::EnsEMBL::Registry->load_all($registry);
my $funcgen_db_adaptor = Bio::EnsEMBL::Registry->get_DBAdaptor($species, 'funcgen');
my $array_adaptor      = Bio::EnsEMBL::Registry->get_adaptor  ($species, 'funcgen', 'array');

# use Bio::EnsEMBL::Funcgen::Config::ArrayFormatConfig;
# my $array_format_config = Bio::EnsEMBL::Funcgen::Config::ArrayFormatConfig->new;

my $arrays_per_object_fh;
my $probeset_sizes_fh;
my $object_names_fh;

open $arrays_per_object_fh, '>' , $arrays_per_object_file;
open $probeset_sizes_fh,    '>' , $probeset_sizes_file;
open $object_names_fh,      '>' , $object_names_file;

my $arrays = $array_adaptor->fetch_all;

foreach my $current_array (@$arrays) {

  my $array_name  = $current_array->name;
  my $array_class = $current_array->class;

  $logger->info("Processing " . $array_class . "\t" . $array_name . "\n" );
#   my $current_array_format_config = $array_format_config->for_array_class($array_class);
  
  (
    my $arrays_per_object, 
    my $probeset_sizes, 
    my $object_names,
  ) = cache_arrays_per_object(
    $funcgen_db_adaptor, 
#     $current_array_format_config,
    $current_array->is_probeset_array,
    $array_name,
  );
  
  $arrays_per_object_fh->print(Dumper({
    array_name        => $array_name,
    array_class       => $array_class,
    arrays_per_object => $arrays_per_object,
  }));
  
  $probeset_sizes_fh->print(Dumper({
    array_name     => $array_name,
    array_class    => $array_class,
    probeset_sizes => $probeset_sizes,
  }));
  
  $object_names_fh->print(Dumper({
    array_name     => $array_name,
    array_class    => $array_class,
    object_names   => $object_names
  }));
  
}

$arrays_per_object_fh ->close;
$probeset_sizes_fh    ->close;
$object_names_fh      ->close;

$logger->finish_log;

=head2 cache_arrays_per_object

 Description: Stores the list of objects in each array
 Arg1: Bio::EnsEMBL::DBSQL::DBAdaptor object
 Arg2: Hash ref containing
 - array_names: array ref
 - array_config: hash ref with vendor configs
 Returntype: arrayref with three hashrefs:
 - List of array names for each object_id
 - Integer for each object_id
 - String name for each object_id

=cut

sub cache_arrays_per_object {
  my ($probe_db, $is_probeset_array, $array_name) = @_;

  my $sql;
  if($is_probeset_array) {
    # Find all probesets belonging to the chosen microarrays
    $sql = '
      select
        probe_set_id,
        probe_set.name,
        array.name, 
        count(probe.probe_id)
      from
        probe
        join probe_set  using(probe_set_id)
        join array_chip on(array_chip.array_chip_id=probe.array_chip_id)
        join array      using(array_id)
      where
        array.name=?
      group by 
        probe.probe_set_id, 
        array.name
    ';
  } else {
    # Find all probes belonging to the chosen microarrays
    # Sometimes the same probe has different names across arrays, hence we concatenate those
    $sql = '
      select
        probe.probe_id, 
        group_concat(probe.name separator "#"), 
        array.name, 
        count(probe.probe_id)
      from
        probe
        join array_chip using (array_chip_id)
        join array using (array_id)
      where
        array.name = ?
      group by
        probe.probe_id, 
        array.name
    ';
  }

  my $sth = $probe_db->dbc->prepare($sql);
  my ($object_id, $array, $probeset_size, $object_name);
  
  $sth->bind_param(1, $array_name);
  $sth->execute();
  $sth->bind_columns(\$object_id, \$object_name, \$array, \$probeset_size);

  my %arrays_per_object = ();
  my %probeset_sizes = ();
  my %object_names = ();
  while($sth->fetch()) {
    $arrays_per_object{$object_id} ||= [];
    push @{$arrays_per_object{$object_id}}, $array;
    $object_names{$object_id} ||= [];
    if ($is_probeset_array) {
      push @{$object_names{$object_id}}, $object_name;
      if ( defined $probeset_sizes{$object_id} && $probeset_size !=  $probeset_sizes{$object_id}) {
        warn("Found probeset(dbID=$object_id) with differing size between arrays:\ti".join(", ", @{$arrays_per_object{$object_id}})." and $array($probeset_size)\n");
      }
      if ($probeset_sizes{$object_id} < $probeset_size) {
        $probeset_sizes{$object_id} = $probeset_size;
      }
    } else {
      push @{$object_names{$object_id}}, (split/#/, $object_name);
    }
  }
  $sth->finish();

  return (
    \%arrays_per_object, 
    \%probeset_sizes, 
    \%object_names,
  );
}

