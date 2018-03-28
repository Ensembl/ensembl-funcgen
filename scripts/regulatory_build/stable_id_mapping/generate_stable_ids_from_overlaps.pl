use strict;
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor;
use Data::Dumper;
use Getopt::Long;

my %options;

GetOptions (
  \%options,
  "all_overlaps=s",
  "new_regulatory_features=s",
  "outfile=s",
  "stable_id_prefix=s",
);

my $all_overlaps            = $options{all_overlaps};
my $new_regulatory_features = $options{new_regulatory_features};
my $outfile                 = $options{outfile};
my $stable_id_prefix        = $options{stable_id_prefix};

open my $all_overlaps_fh, '<', $all_overlaps;

my $overlaps;
my $max_seen_stable_id_number;

($overlaps, $max_seen_stable_id_number) = filter_for_compatible_overlaps($all_overlaps_fh);

(
  my $stable_id_hash,
  my $num_mapped_stable_ids,
  my $num_new_stable_ids
  
)= generate_stable_id_assignments(
  $stable_id_prefix, 
  $overlaps, 
  $max_seen_stable_id_number
);

print "Writing stable id assignments to $outfile\n";

open my $out_fh, ">" . $outfile;
foreach my $regulatory_feature_id (keys %$stable_id_hash) {

  my $stable_id = $stable_id_hash->{$regulatory_feature_id};
  $out_fh->print(
    join "\t", $regulatory_feature_id, $stable_id
  );
  $out_fh->print("\n");
}
$out_fh->close();

print "$num_mapped_stable_ids stable ids were mapped.\n";
print "$num_new_stable_ids stable ids were newly assigned.\n";

sub generate_stable_id_assignments {

  my $stable_id_prefix          = shift;
  my $overlaps                  = shift;
  my $max_seen_stable_id_number = shift;

  # Go through overlaps in order of increasing overlap length. This means that you should always
  # overwrite an overlap with a later one.
  # 
  # The overlap length is in column 13, so the index is 12.
  sub cmp_overlaps {
    return $a->[12] <=> $b->[12];
  }
  
  my %old_pref = ();
  my %new_pref = ();
  foreach my $entry (sort cmp_overlaps @{$overlaps}) {
    $old_pref{$entry->[4]} = $entry->[10];
    $new_pref{$entry->[10]} = $entry->[4];
  }

  my %stable_id_hash = ();
  my $next_free_id = $max_seen_stable_id_number + 1;
  open my $in, "<", $new_regulatory_features;

  my $num_mapped_stable_ids = 0;
  my $num_new_stable_ids    = 0;

  while (my $line = <$in>) {
    chomp $line;
    my @items = split "\t", $line;
    
    # This is the name of the feature in the bed file. e.g.: tfbs_2, ctcf_1, 
    # dnase_2, dnase_3, tfbs_3, tfbs_4, proximal_1, proximal_2, tss_1, 
    # tfbs_6, tfbs_7, ctcf_2, ctcf_3, tfbs_8 etc.
    #
    my $new_stable_id_number = $items[4];

    my $stable_id = undef;

    if (exists $new_pref{$new_stable_id_number} && $old_pref{$new_pref{$new_stable_id_number}} eq $new_stable_id_number) {
      $stable_id = $new_pref{$new_stable_id_number};
      $num_mapped_stable_ids++;
    } else {
      $stable_id = $next_free_id;
      $next_free_id += 1;
      $num_new_stable_ids++;
    }

    # Creating stable id string, composed of prefix + 11 digit integer, front padded with 0s
    # 
    my $regulatory_feature_id = $items[5];
    $stable_id_hash{$regulatory_feature_id} = $stable_id_prefix . sprintf("%011d", $stable_id);

  }
  close $in;
  return \%stable_id_hash, $num_mapped_stable_ids, $num_new_stable_ids;
}

sub filter_for_compatible_overlaps {

  my $bedfile_fh = shift;

  my $max_seen_stable_id_number = 0;
  my @overlaps;

  while (my $current_bed_file_line = <$bedfile_fh>) {

    chomp $current_bed_file_line;
    my @bed_file_fields = split "\t", $current_bed_file_line;
    
    # Feature types must be identical for transferring the stable id.
    #
    my $regulatory_feature_type_old             = $bed_file_fields[3];
    my $regulatory_feature_old_stable_id_number = $bed_file_fields[4];
    my $regulatory_feature_type_new             = $bed_file_fields[9];
    
    if ($regulatory_feature_type_old eq $regulatory_feature_type_new) {
      print "Compatible     $current_bed_file_line\n";
      push @overlaps, \@bed_file_fields;
    } else {
      print "Not compatible $current_bed_file_line\n";
    }

    # Look for maximum ID number among the old features
    if ($max_seen_stable_id_number < $regulatory_feature_old_stable_id_number) {
      $max_seen_stable_id_number = $regulatory_feature_old_stable_id_number;
    }
  }
  return (  \@overlaps, $max_seen_stable_id_number );
}
