#!/usr/bin/env perl

use strict;
use JSON;
use Bio::EnsEMBL::Registry;
use File::Spec;
use Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Hive::DBSQL::DBConnection;
use Getopt::Long;

# export_qc_chance.pl --registry /nfs/users/nfs_m/mn1/work_dir_ftp/lib/ensembl-funcgen/registry.pm --species homo_sapiens --output_file testun/qc_phantom_peaks.json 
# export_qc_chance.pl --registry /nfs/users/nfs_m/mn1/work_dir_ftp/lib/ensembl-funcgen/registry.pm --species homo_sapiens | less

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

my $sth = $db_connection->prepare(
  qq(
    select 
      signal_bam_file.path signal_bam_file, 
      signal_bigwig_file.path signal_bigwig_file, 
      alignment.name signal_alignment_name,
      control_bam_file.path control_bam_file,
      control_bigwig_file.path control_bigwig_file,
      control_alignment.name control_alignment_name,
      analysis_description.display_label analysis,
      epigenome.display_label as epigenome,
      chance.* 
    from 
      chance 
      join alignment on (chance.signal_alignment_id=alignment.alignment_id) 
      left join data_file signal_bam_file    on (signal_bam_file.data_file_id    = alignment.bam_file_id)
      left join data_file signal_bigwig_file on (signal_bigwig_file.data_file_id = alignment.bigwig_file_id)
      join analysis_description on (chance.analysis_id=analysis_description.analysis_id)
      join experiment the_signal using (experiment_id)
      join epigenome using (epigenome_id)
      left join experiment control on (the_signal.control_id=control.experiment_id)
      join alignment control_alignment on (control.experiment_id=control_alignment.experiment_id)
      left join data_file control_bam_file    on ( control_bam_file.data_file_id = control_alignment.bam_file_id )
      left join data_file control_bigwig_file on ( control_bam_file.data_file_id = control_alignment.bigwig_file_id )
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

$output_fh->print("[\n");

my $is_first = 1;

while (my $hash_ref = $sth->fetchrow_hashref) {

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
  
  $hash_ref->{control_alignment} = {
    name        => $hash_ref->{control_alignment_name},
    bam_file    => $hash_ref->{control_bam_file},
    bigwig_file => $hash_ref->{control_bigwig_file},
  };
  delete $hash_ref->{control_alignment_name};
  delete $hash_ref->{control_bam_file};
  delete $hash_ref->{control_bigwig_file};

  $hash_ref->{signal_alignment} = {
    name        => $hash_ref->{signal_alignment_name},
    bam_file    => $hash_ref->{signal_bam_file},
    bigwig_file => $hash_ref->{signal_bigwig_file},
  };
  delete $hash_ref->{signal_alignment_name};
  delete $hash_ref->{signal_bam_file};
  delete $hash_ref->{signal_bigwig_file};


  delete $hash_ref->{path};
  delete $hash_ref->{analysis_id};
  delete $hash_ref->{name};
  delete $hash_ref->{chance_id};
}
