#!/usr/bin/env perl


use strict;
use warnings;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Funcgen::Utils::EFGUtils qw( run_system_cmd run_backtick_cmd );

#fastaexplode -f GCA_000001405.15_GRCh38_no_alt_analysis_set.fna -d tmp_ncbi

my ($user, $fasta_dir, $host, $dbname) = @ARGV;

my $db = Bio::EnsEMBL::DBSQL::DBAdaptor->new
 (-host    => $host,
  -dbname  => $dbname,
  -user    => $user);
  
my $slice_a = $db->get_SliceAdaptor;

#This should default to GRCh38
my @slices = @{$slice_a->fetch_all('toplevel', undef, undef, 1)};
#no inc non ref, but inc dups for correct slice names.

#Iterate through the slices, identify the repspective ncbi file
#replace the header & validate the length

#This is the correct header format
#>GL000200.1 dna:supercontig supercontig:GRCh37:GL000200.1:1:187035:1 supercontig GL000200.1
#which is from the sequence_dump.pl script

#$slice->seq_region_name().' '.$slice->moltype().':'.$slice->coord_system_name().' '.$slice->name();
opendir(DIR, $fasta_dir) or die $!;
my @fastas = readdir(DIR);

foreach my $slice(@slices){
  my $name = $slice->seq_region_name;  
  print "Validating slice $name\n";
  
  #Matching criterion are not exhaustive, but appear to work in this case
  $name = 'M' if $name eq 'MT';  
  my @files = grep(/[_r]$name\./, @fastas);
  
  #v1 files are caught due to the v being transposed to a . handy
  #Validating slice KI270750.1
  #Found fastas:
  #chrUn_KI270750v1.fa
  
  if(! @files){
    @files = grep(/$name/, @fastas);      
  }
  
  print "Found fastas:\n\t".join("\n\t", @files)."\n";
   
  if(scalar(@files) != 1){
    die('Failed to identify unique fasta file');  
  }
  
  my $file = $files[0];
  my $header = '>'.$slice->seq_region_name().' '.$slice->moltype().':'.$slice->coord_system_name().' '.$slice->name();
  
  my $h_length = length($header) -1 ; #length counts \n, so we need to -1
  my $s_length = $slice->length;
  
  my $out_file = "$fasta_dir/../ncbi_converted/$file";
  my $cmd = "sed 's/>.*/$header/' $fasta_dir/$file > $out_file";
  run_system_cmd($cmd); 
  
  #wc count also counts \n's so deduct the number of lines -1 for the last record
  $cmd = "wc -ml $out_file";
  my ($num_lines, $f_length) = split(' ', run_backtick_cmd($cmd));
  
  #warn "$f_length $num_lines";
  $f_length -= $num_lines + 1;
  my $r_length = ($h_length + $s_length);
  
  if($f_length != $r_length){
    die("Length mismatch:\t File $f_length\tExpected record length $r_length");  
  }    
}


1;