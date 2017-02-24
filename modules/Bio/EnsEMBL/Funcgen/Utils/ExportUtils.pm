package Bio::EnsEMBL::Funcgen::Utils::ExportUtils;

use warnings;
use strict;

use base qw( Exporter );
use vars qw( @EXPORT_OK );

@EXPORT_OK = qw(
  assert_source_files_exist
  assert_destination_file_names_uniqe
);

sub assert_source_files_exist {

  my @file_list = @_;
  my @not_existing_files;
  foreach my $current_file (@file_list) {
    if (!-e $current_file) {
      push @not_existing_files, $current_file;
    }
  }
  return @not_existing_files;
}

sub assert_destination_file_names_uniqe {

  my $source_file_to_destination_file_map = shift;
  my $destination_file_to_more_than_one_source_file = find_duplicates($source_file_to_destination_file_map);
  
  my $duplicates_found = keys %$destination_file_to_more_than_one_source_file;
  
  if ($duplicates_found) {
    use Data::Dumper;
    use Carp;
    confess(
      "There is a problem with naming the new files. The following have more "
      . "than one source file that would map to it: " 
      . Dumper(\$destination_file_to_more_than_one_source_file)
    );
  }
  return
}

sub find_duplicates {

  my $source_file_to_destination_file_map = shift;
  
  # Map from destination file to the source files that map to it. If
  # there is more than one, we have a problem.
  #
  my %destination_file_to_source_file_list;
  
  my %bam_file_to_from = reverse %$source_file_to_destination_file_map;
  my @destination_file_list = values %bam_file_to_from;
  
  for my $current_destination_file (@destination_file_list) {
    if (!exists $destination_file_to_source_file_list{$current_destination_file}) {
      $destination_file_to_source_file_list{$current_destination_file} = [];
    }
    my $source_file = $bam_file_to_from{$current_destination_file};
    push @{$destination_file_to_source_file_list{$current_destination_file}}, $source_file;
  }
  
  # Like %destination_file_to_source_file_list, but only with those 
  # destination files to which more than one source file maps. These are the
  # ones that will cause problems.
  #
  my %destination_file_to_more_than_one_source_file;
  foreach my $destination_file (keys %destination_file_to_source_file_list) {
  
    my $number_of_source_files_mapped_to_this_file = @{$destination_file_to_source_file_list{$destination_file}};
    
    if ($number_of_source_files_mapped_to_this_file > 1) {
      $destination_file_to_more_than_one_source_file{$destination_file} = $destination_file_to_source_file_list{$destination_file};
    }
  }
  return \%destination_file_to_more_than_one_source_file;
}

1;
