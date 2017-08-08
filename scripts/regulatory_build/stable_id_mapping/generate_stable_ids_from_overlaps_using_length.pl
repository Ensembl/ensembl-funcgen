use strict;
use Data::Dumper;
use Getopt::Long;

my $all_overlaps;
my $source_regulatory_features;
my $target_regulatory_features;
my $outfile;
my $stable_id_prefix;

GetOptions (
'all_overlaps=s'                => \$all_overlaps,
'source_regulatory_features=s'  => \$source_regulatory_features,
'target_regulatory_features=s'  => \$target_regulatory_features,
'outfile=s'                     => \$outfile,
'stable_id_prefix=s'            => \$stable_id_prefix,
);

## Identify the maximum stable id seen in the "old" regulatory feature dataset
open my $source_regulatory_features_fh, '<', $source_regulatory_features;
my $max_seen_stable_id_number = find_max_seen_stable_id( $source_regulatory_features_fh );
$source_regulatory_features_fh->close();

## Retain overlaps between the same regulatory feature types and
## for each one of those only the one with the longest overlap
open my $all_overlaps_fh, '<', $all_overlaps;
my $stable_id_mappings_ref_hash = filter_overlaps( $all_overlaps_fh );
$all_overlaps_fh->close();

open my $target_regulatory_features_fh, "<", $target_regulatory_features;

(
my $stable_id_hash,
my $num_mapped_stable_ids,
my $num_new_stable_ids
)= generate_stable_id_assignments(
$target_regulatory_features_fh,
$stable_id_prefix,
$stable_id_mappings_ref_hash,
$max_seen_stable_id_number
);

$target_regulatory_features_fh->close();
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

sub find_max_seen_stable_id {
  
  my $bedfile_fh = shift;
  
  my $max_seen_stable_id_number;
  
  while (my $current_bed_file_line = <$bedfile_fh>) {
    # 18	76429380	76430144	Open chromatin	00000105157	535842
    chomp $current_bed_file_line;
    my ($chr, $start, $end, $regulatory_feature_type, $source_regulatory_feature_stable_id, $source_regulatory_feature_db_id) = split "\t", $current_bed_file_line;
    
    # Look for maximum ID number among the old features
    if ($max_seen_stable_id_number < $source_regulatory_feature_stable_id) {
      $max_seen_stable_id_number = $source_regulatory_feature_stable_id;
    }
  }
  
  print "The maximum stable id seen is: $max_seen_stable_id_number\n";
  return( $max_seen_stable_id_number );
}

sub filter_overlaps {
  
  my $all_overlaps_fh = shift;
  
  my $same_feature_type_overlaps = retain_same_feature_type_overlaps($all_overlaps_fh);
 
  my %old_preference = ();
  my %new_preference = ();

  my $no_overlaps_retained=0;
  my $no_overlaps_dismissed=0;
  
  # The 13th column in th BED file refers to the length of the overlap as reported by bedtools intersect.
  # which is the 12th index in the array
  foreach my $overlap ( sort { $b->[12] <=> $a->[12] } @{$same_feature_type_overlaps} ) {
    my $source_regulatory_feature_stable_id  = $overlap->[4];
    my $target_regulatory_feature_db_id = $overlap->[11];
    
    if( ! (exists $old_preference{$source_regulatory_feature_stable_id} || exists $new_preference{$target_regulatory_feature_db_id} ) ) {
      $old_preference{$source_regulatory_feature_stable_id} = $target_regulatory_feature_db_id;
      $new_preference{$target_regulatory_feature_db_id} = $source_regulatory_feature_stable_id;
      
      $no_overlaps_retained += 1;
      
    } else {
      
      $no_overlaps_dismissed += 1;
      #print "Either the $old_regulatory_feature_stable_id or $new_regulatory_feature_type_db_id have already been used\n";
    }
  }
  
  print "$no_overlaps_retained overlap were retained.\n";
  print "$no_overlaps_dismissed overlaps were discarded.\n";
  
  return ( \%new_preference );
}

sub retain_same_feature_type_overlaps {
  my $bedfile_fh = shift;
  
  my @same_feature_type_overlaps;

  while (my $current_bed_file_line = <$bedfile_fh>) {
    
    chomp $current_bed_file_line;
    
    #1	13201	13800	Enhancer	00000000001	623456	1	13371	13724	Open chromatin	00000341931	160079	353
    my ($source_chr, $source_start, $source_end, $source_regulatory_feature_type, $source_regulatory_feature_stable_id, $source_regulatory_feature_db_id, $target_chr, $target_start, $target_end, $target_regulatory_feature_type, $target_regulatory_feature_stable_id, $target_regulatory_feature_db_id, $overlap_length)  = split "\t", $current_bed_file_line;
    
    # Feature types must be identical for transferring the stable id.
    if ($source_regulatory_feature_type eq $target_regulatory_feature_type) {
      #print "Compatible $current_bed_file_line\n";
      my @bed_file_line = ( $source_chr, $source_start, $source_end, $source_regulatory_feature_type, $source_regulatory_feature_stable_id, $source_regulatory_feature_db_id, $target_chr, $target_start, $target_end, $target_regulatory_feature_type, $target_regulatory_feature_stable_id, $target_regulatory_feature_db_id, $overlap_length );   
      push @same_feature_type_overlaps, \@bed_file_line;
    }
  }
  
  return( \@same_feature_type_overlaps );
}

sub generate_stable_id_assignments {
  
  my $target_regulatory_features_fh  = shift;
  my $stable_id_prefix               = shift;
  my $stable_id_mappings_ref_hash    = shift;
  my $max_seen_stable_id_number      = shift;
  
  my %stable_id_mappings_hash = %$stable_id_mappings_ref_hash;
  
  my @no_stable_id_mappings = keys %$stable_id_mappings_ref_hash;
  print "Recovered ". @no_stable_id_mappings." mappings!!\n";
  
  my %stable_id_hash = ();
  my $next_free_id = $max_seen_stable_id_number + 1;
  
  my $num_mapped_stable_ids = 0;
  my $num_new_stable_ids    = 0;
  
  while (my $line = <$target_regulatory_features_fh>) {
    chomp $line;
    #15	102118695	102119230	TF binding site	00000368862	1
    my ($target_chr, $target_start, $target_end, $target_regulatory_feature_type, $target_regulatory_feature_stable_id, $target_regulatory_feature_db_id) = split "\t", $line;
    
    my $updated_stable_id = undef;
    
    if ( exists $stable_id_mappings_hash{$target_regulatory_feature_db_id} ) {
      $updated_stable_id = $stable_id_mappings_hash{$target_regulatory_feature_db_id};
      $num_mapped_stable_ids++;
    } else {
      $updated_stable_id = $next_free_id;
      $next_free_id += 1;
      $num_new_stable_ids++;
    }
    
    # Creating stable id string, composed of prefix + 11 digit integer, front padded with 0s
    $stable_id_hash{$target_regulatory_feature_db_id} = $stable_id_prefix . sprintf("%011d", $updated_stable_id);
    
  }
  return \%stable_id_hash, $num_mapped_stable_ids, $num_new_stable_ids;
}





