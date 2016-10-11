#!/usr/bin/env perl

use strict;
use JSON;
use Bio::EnsEMBL::Registry;
use File::Spec;
use Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Hive::DBSQL::DBConnection;
use Getopt::Long;

# export_qc_proportion_of_reads_in_peaks.pl --registry /nfs/users/nfs_m/mn1/work_dir_ftp/lib/ensembl-funcgen/registry.pm --species homo_sapiens --output_file testun/qc_proportion_of_reads_in_peaks.json 
# export_qc_proportion_of_reads_in_peaks.pl --registry /nfs/users/nfs_m/mn1/work_dir_ftp/lib/ensembl-funcgen/registry.pm --species homo_sapiens | less
 
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
      analysis_description.display_label as analysis,
      feature_type.name as feature_type,
      epigenome.display_label as epigenome,
      feature_set_qc_prop_reads_in_peaks.prop_reads_in_peaks as proportion_of_reads_in_peaks,
      feature_set_qc_prop_reads_in_peaks.total_reads as total_number_of_reads,
      result_set.name alignment_name,
      experiment.is_control,
      bam_file.path bam_file,
      bigwig_file.path bigwig_file
    from 
      feature_set_qc_prop_reads_in_peaks 
      join feature_set using (feature_set_id)
      join analysis_description on (analysis_description.analysis_id=feature_set.analysis_id)
      join feature_type using (feature_type_id)
      join epigenome using (epigenome_id)
      join data_set using (feature_set_id)
      join supporting_set using (data_set_id)
      join result_set on (supporting_set.supporting_set_id=result_set.result_set_id and supporting_set.type="result")
      join experiment on (result_set.experiment_id=experiment.experiment_id)
      left join dbfile_registry bam_file    on (bam_file.table_name='result_set'    and bam_file.table_id=result_set.result_set_id    and bam_file.file_type='BAM')
      left join dbfile_registry bigwig_file on (bigwig_file.table_name='result_set' and bigwig_file.table_id=result_set.result_set_id and bigwig_file.file_type='BIGWIG')
    order by result_set.name
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

  $hash_ref->{epigenome} = $fetchQCRelatedData->fetch_xrefs_for_epigenome($hash_ref->{epigenome});

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


  delete $hash_ref->{path};
  delete $hash_ref->{analysis_id};
  delete $hash_ref->{feature_set_id};
  delete $hash_ref->{feature_set_qc_prop_reads_in_peaks_id};
}
