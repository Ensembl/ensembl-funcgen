#!/usr/bin/env perl

use strict;
use Data::Dumper;
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor;
use Getopt::Long;

  use Bio::EnsEMBL::Funcgen::Probe;
  use Bio::EnsEMBL::Funcgen::Array;


=head1

  compute_read_file_statistics.pl \
    -registry /homes/mnuhn/work_dir_ersa/lib/ensembl-funcgen/registry_new_mouse_encode_data.pm \
    -species mus_musculus \
    -read_file_id 1223

=cut

my $registry;
my $species;
my $read_file_id;

GetOptions (
   'registry=s'     => \$registry,
   'species=s'      => \$species,
   'read_file_id=s' => \$read_file_id,
);

Bio::EnsEMBL::Registry->load_all($registry);

my $read_file_adaptor = Bio::EnsEMBL::Registry->get_adaptor($species, 'funcgen', 'ReadFile');
my $read_file = $read_file_adaptor->fetch_by_dbID($read_file_id);

my $fastq_file = $read_file->file;

if (! -e $fastq_file) {
  die("Fastq file for read file with id ${read_file_id}: $fastq_file");
}

my $file_size = -s $fastq_file;

use Bio::EnsEMBL::Funcgen::Utils::Fastq::Parser;
my $parser = Bio::EnsEMBL::Funcgen::Utils::Fastq::Parser->new;

my $cmd = "zcat $fastq_file |";

open my $fh, $cmd;

my $number_of_reads = 0;
my $average_length  = 0;

while (
    (my $record = $parser->parse_next_record($fh))
#     && 
#     ($number_of_reads < 20)
  ) {

  my $sequence = $record->[1];
  
  $number_of_reads++;
#   print "$sequence\n";
  
  my $sequence_length = length($sequence);
  
  $average_length = $average_length - $average_length  / $number_of_reads;
  $average_length = $average_length + $sequence_length / $number_of_reads;
  
#   print "$sequence\n";
#   print "$sequence_length\n";
#   print "$average_length\n";
  
}

$read_file->number_of_reads($number_of_reads);
$read_file->file_size($file_size);
$read_file->read_length($average_length);

$read_file_adaptor->update($read_file);

print <<REPORT

File: $fastq_file

number of reads = $number_of_reads
file size       = $file_size
average length  = $average_length

REPORT
;


