#!/usr/bin/env perl

use strict;
use JSON;
use Bio::EnsEMBL::Registry;
use File::Spec;
use Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Hive::DBSQL::DBConnection;
use Getopt::Long;

# export_qc_mapped_reads.pl --registry /nfs/users/nfs_m/mn1/work_dir_ftp/lib/ensembl-funcgen/registry.pm --species homo_sapiens --output_file testun/qc_mapped_reads.json 
# export_qc_mapped_reads.pl --registry /nfs/users/nfs_m/mn1/work_dir_ftp/lib/ensembl-funcgen/registry.pm --species homo_sapiens | less

my $registry;
my $species;
my $output_file;

GetOptions (
   'registry=s'    => \$registry,
   'species=s'     => \$species,
   'output_file=s' => \$output_file,
);

Bio::EnsEMBL::Registry->load_all($registry);

use Bio::EnsEMBL::Utils::Logger;
my $logger = Bio::EnsEMBL::Utils::Logger->new();

my $funcgen_db_adaptor = Bio::EnsEMBL::Registry->get_DBAdaptor($species, 'funcgen');

my $json = JSON->new->utf8;
$json->pretty(1);
$json->canonical(1);

my $db_connection = $funcgen_db_adaptor->dbc;

use Bio::EnsEMBL::Funcgen::Ftp::FetchQCRelatedData;
my $fetchQCRelatedData = Bio::EnsEMBL::Funcgen::Ftp::FetchQCRelatedData->new(
  -db_connection => $db_connection,
);

my $sth = $db_connection->prepare(
  qq(
    select 
      result_set.name as alignment_name,
      epigenome.display_label as epigenome,
      bam_file.path bam_file,
      bigwig_file.path bigwig_file,
      total.qc_passed_reads num_total, 
      experiment.is_control,
      mapped.qc_passed_reads num_mapped, 
      mapped.qc_passed_reads/total.qc_passed_reads proportion_of_reads_mapped
    from 
      result_set_qc_flagstats total join result_set_qc_flagstats mapped on (total.category="in total" and mapped.category="mapped" and total.result_set_id=mapped.result_set_id)
      join result_set on (total.result_set_id=result_set.result_set_id) 
      join epigenome using (epigenome_id)
      join experiment on (result_set.experiment_id=experiment.experiment_id)
      left join dbfile_registry bam_file    on (bam_file.table_name='result_set'    and bam_file.table_id=result_set.result_set_id    and bam_file.file_type='BAM')
      left join dbfile_registry bigwig_file on (bigwig_file.table_name='result_set' and bigwig_file.table_id=result_set.result_set_id and bigwig_file.file_type='BIGWIG')
  )
);

$sth->execute;

my $output_fh;
if ($output_file) {
  $logger->info("The features will be written to " . $output_file ."\n");

  use File::Basename;
  my $ftp_dir = dirname($output_file);

  use File::Path qw(make_path);
  make_path($ftp_dir);

  use IO::File;
  $output_fh = IO::File->new(">$output_file");
} else {
  $output_fh = *STDOUT;
}

my $is_first = 1;

$output_fh->print("[\n");

while (my $hash_ref = $sth->fetchrow_hashref) {

  $hash_ref->{sequence_files} = $fetchQCRelatedData->fetch_input_subset_data_for_result_set($hash_ref->{alignment_name});
  $hash_ref->{epigenome}      = $fetchQCRelatedData->fetch_xrefs_for_epigenome($hash_ref->{epigenome});

  translate($hash_ref);

  if ($is_first) {
    $is_first = undef;
  } else {
    $output_fh->print(",\n");
  }
  $output_fh->print($json->encode($hash_ref));
}

$output_fh->print("\n]");

sub translate {
  my $hash_ref = shift;
  
  $hash_ref->{alignment} = {
    name        => $hash_ref->{alignment_name},
    bam_file    => $hash_ref->{bam_file},
    bigwig_file => $hash_ref->{bigwig_file},
    is_control  => $hash_ref->{is_control},
  };
  delete $hash_ref->{alignment_name};
  delete $hash_ref->{bam_file};
  delete $hash_ref->{bigwig_file};
  delete $hash_ref->{is_control};
}
