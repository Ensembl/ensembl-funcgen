#!/usr/bin/env perl

=head1 LICENSE


Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2018] EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

  Questions may also be sent to the Ensembl help desk at
  <http://www.ensembl.org/Help/Contact>.

=cut

# Project the VISTA enhancer genomic coordinates from one assembly to another prior to loading them on the funcgen db 
#
# Example usage:
# project_external_features_vista.pl \
#       -registry registry.pm \
#       -species=mus_musculus \
#       -old_assembly="NCBIM37" \
#       -new_assembly="GRCm38"
#       -in=/nfs/production/panda/ensembl/funcgen/source_files/old/warehouse/ensembl07/source_files/external_features/VISTA/VISTA_4Apr17_header \
#       -out=VISTA_4Apr17_header_mouse_89_38 \

use strict;
use warnings;
use feature qw(say);

use Data::Dumper;
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Getopt::Long;


my ($registry, $species, $old_assembly, $new_assembly, $in, $out);

GetOptions
  (
   "registry=s"      => \$registry,
   "species=s"       => \$species,
   "old_assembly=s", => \$old_assembly,
   "new_assembly=s", => \$new_assembly,
   "in=s",           => \$in,
   "out=s",          => \$out,
  );

Bio::EnsEMBL::Registry->load_all($registry);

my $core_db_adaptor = Bio::EnsEMBL::Registry->get_DBAdaptor($species, 'core');
my $core_dbc = $core_db_adaptor->dbc;

my $vista_species_name;
if ($species eq 'homo_sapiens') {
  $vista_species_name = 'Human';
}
elsif ($species eq 'mus_musculus') {
  $vista_species_name = 'Mouse';
}

my $slice_adaptor  = $core_db_adaptor->get_SliceAdaptor; 
open(my $fh_out,'>',$out) or die "Can not open/access '".$out."'\n$!";

# The first lines of the FASTA file downloaded look like:
# >Human|chr16:86430087-86430726 | element 1 | positive  | neural tube[12/12] | hindbrain (rhombencephalon)[12/12] | limb[3/12] | cranial nerve[8/12]

open(my $fh, '<', $in) or die "Can not open/access '".$in."'\n$!";
  while(my $line = <$fh>){
    next if($line !~ /^>$vista_species_name/);
    chomp($line);
    my @line = split('\|',$line);
    $line[1] =~ /chr([^:]+):(\d+)-(\d+)/;
    # $1=chromosome; $2=start; $3=end

    my $slice = $slice_adaptor->fetch_by_region('chromosome', $1, $2, $3, 1, $old_assembly);
    
    my $tmp   = $slice->project('chromosome', $new_assembly);
    
    if(scalar @$tmp != 1){
      say "Found [".scalar @$tmp."] projections: $line";
      next;
    }
    
    my $ps    = shift(@{$tmp});
    next if(! defined $ps);
    my $slice_37 = $ps->to_Slice;
    my @slice_37 = split(':',$slice_37->name);
    
   # chromosome:GRCh38:X:124216977:124218737:1
    my $new_coords  = 'chr'.$slice_37[2].':'.$slice_37[3].'-'.$slice_37[4];
    my $new_line = ">$vista_species_name|$new_coords |";
    shift @line ;
    shift @line ;
    $new_line .= join('|', @line);
    say $fh_out $new_line;
  }
close($fh);


