use strict;
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor;
use Data::Dumper;
use Getopt::Long;

=head1

perl scripts/regulatory_build/stable_id_mapping/map_coordinates_to_assembly.pl \
  --url 'mysql://ensro@ens-genomics2:3306/homo_sapiens_core_85_38?group=core&species=homo_sapiens' \
  --species homo_sapiens \
  --bed_file_before_mapping /lustre/scratch109/ensembl/funcgen/mn1/stable_id_mapping/homo_sapiens_GrCh38_new.bed \
  --bed_file_after_mapping /lustre/scratch109/ensembl/funcgen/mn1/stable_id_mapping/homo_sapiens_GrCh38_mapped_to_GrCh37.bed \
  --assembly_from GRCh38 \
  --assembly_to GRCh37

=cut

my $bed_file_before_mapping;
my $bed_file_after_mapping;

my $assembly_from;
my $assembly_to;

my $species;

# The url to the core database that know how to map between the two assemblies.
my $url;

my %options;
GetOptions (
   'url=s'                     => \$url,
   'species=s'                 => \$species,
   'bed_file_before_mapping=s' => \$bed_file_before_mapping,
   'bed_file_after_mapping=s'  => \$bed_file_after_mapping,
   'assembly_from=s'           => \$assembly_from,
   'assembly_to=s'             => \$assembly_to,
);

Bio::EnsEMBL::Registry->load_registry_from_url($url, 1);

my $core_slice_adaptor = Bio::EnsEMBL::Registry->get_adaptor($species, 'core', 'slice');

open my $bed_file_before_fh, "<$bed_file_before_mapping";
open my $bed_file_after_fh,  ">$bed_file_after_mapping";

LINE: while (my $current_bed_file_line = <$bed_file_before_fh>) {

  chomp $current_bed_file_line;
  my @bed_file_fields = split "\t", $current_bed_file_line;

  my $seq_region_name = $bed_file_fields[0];
  my $start           = $bed_file_fields[1];
  my $end             = $bed_file_fields[2];
  
  my @fetch_by_region_parameters = (
    'chromosome',
    $seq_region_name,
    $start,
    $end,
    1,
    $assembly_from,
  );
  my $seq_region_from = $core_slice_adaptor->fetch_by_region(@fetch_by_region_parameters);
  
  if (! defined $seq_region_from) {
    print "Couldn't find ${seq_region_name}:${start}-${end} on chromosome level, skipping.";
    next LINE;
  }
  
  my $projected_slices = $seq_region_from->project(
    'chromosome', 
    $assembly_to
  );
  if (@$projected_slices==0) {
    print("Can't map ${seq_region_name}:${start}-${end}, skipping\n");
    next LINE;
  }
  if (@$projected_slices>1) {
    print("Got ". @$projected_slices ." slices, skipping\n");
    next LINE;
#     my @all_slices = reverse sort { $a->length <=> $b->length } map { $_->to_Slice } @$projected_slices;
#     
#     foreach my $current_projected_slice (@all_slices) {
#       print " - " . $current_projected_slice->length . "\n";
#     }
  }
  my $projected_slice = $projected_slices->[0]->to_Slice;

#   print
#     join " ",
#       $seq_region_name,
#       $start,
#       $end;
# 
#   print " -> ";
# 
#   print 
#     join " ",
#       $projected_slice->coord_system->name,
#       $projected_slice->coord_system->version,
#       $projected_slice->seq_region_name,
#       $projected_slice->start,
#       $projected_slice->end,
#       $projected_slice->strand,
#   ;
#   print "\n";
  
  $bed_file_fields[0] = $projected_slice->seq_region_name;
  $bed_file_fields[1] = $projected_slice->start;
  $bed_file_fields[2] = $projected_slice->end;
  
  my $new_bed_line = join "\t", @bed_file_fields;
  
  $bed_file_after_fh->print($new_bed_line);
  $bed_file_after_fh->print("\n");
}

$bed_file_before_fh->close();
$bed_file_after_fh->close;


