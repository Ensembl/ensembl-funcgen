#!/usr/bin/env perl

=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2024] EMBL-European Bioinformatics Institute

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
  developers list at <ensembl-dev@ebi.ac.uk>.

  Questions may also be sent to the Ensembl help desk at
  <helpdesk@ensembl.org>.

=head1 NAME

  A script to compute the number of bases covered by regulatory features both in total and by feature type

  perl scripts/regulatory_build/compute_enhancer_coverage_bedtools.pl --registry /homes/mnuhn/work_dir_regbuild_script/lib/ensembl-funcgen/registry.pm --species homo_sapiens --tempdir foobar
  
  perl scripts/regulatory_build/compute_enhancer_coverage_bedtools.pl --registry /homes/mnuhn/work_dir_regbuild_script/lib/ensembl-funcgen/registry.pm --species mus_musculus --tempdir foobar

=cut

use strict;
use Getopt::Long;
use Bio::EnsEMBL::Registry;
use Data::Dumper;
use Bio::EnsEMBL::DBSQL::SliceAdaptor;
use Bio::EnsEMBL::Slice;
use Bio::EnsEMBL::Utils::Logger;

my %options;
GetOptions (
    \%options,
    "species|s=s",
    "registry|r=s",
    "tempdir|t=s",
 );

use Hash::Util qw( lock_keys );
lock_keys( %options );

my $species  = $options{'species'};
my $registry = $options{'registry'};
my $tempdir  = $options{'tempdir'};

my $skip_dumping_if_files_exist_already = 0;

Bio::EnsEMBL::Registry->load_all($registry);

my $core_adaptor    = Bio::EnsEMBL::Registry->get_DBAdaptor($species, 'core');
my $funcgen_adaptor = Bio::EnsEMBL::Registry->get_DBAdaptor($species, 'funcgen');

my $logger = Bio::EnsEMBL::Utils::Logger->new();
$logger->init_log;

my $feature_set_adaptor        = $funcgen_adaptor->get_FeatureSetAdaptor;
my $external_feature_adaptor   = $funcgen_adaptor->get_ExternalFeatureAdaptor;
my $regulatory_build_adaptor   = $funcgen_adaptor->get_RegulatoryBuildAdaptor;
my $regulatory_feature_adaptor = $funcgen_adaptor->get_RegulatoryFeatureAdaptor;
my $regulatory_build           = $regulatory_build_adaptor->fetch_current_regulatory_build;

use File::Path qw(make_path remove_tree);
make_path($tempdir);

my $regulatory_features_file = $tempdir . "/regulatory_features.bed";

my $have_to_dump_regulatory_features = 
  ! -e $regulatory_features_file || (!$skip_dumping_if_files_exist_already && -e $regulatory_features_file);

if ($have_to_dump_regulatory_features) {

  $logger->info("Exporting regulatory features to $regulatory_features_file ...");
  
  export_regulatory_features($regulatory_features_file);
  sort_bed_file             ($regulatory_features_file);
  
  $logger->info(" done.\n");
}

if (!$have_to_dump_regulatory_features) {
  $logger->info("Skipping dumping of regulatory features.");
}

my $regulatory_build_statistic_adaptor = $funcgen_adaptor->get_RegulatoryBuildStatisticAdaptor;

$logger->info("Removing enhancer statistics ...");

my $dbc = $funcgen_adaptor->dbc;

$dbc->do('delete from regulatory_build_statistic where statistic like "total_enhancers_checked_vista"');
$dbc->do('delete from regulatory_build_statistic where statistic like "num_enhancers_overlapping_vista"');
$dbc->do('delete from regulatory_build_statistic where statistic like "total_enhancers_checked_fantom"');
$dbc->do('delete from regulatory_build_statistic where statistic like "num_enhancers_overlapping_fantom"');

$logger->info(" done.\n");

use constant {
  FANTOM_FEATURE_SET_NAME => 'FANTOM',
  VISTA_FEATURE_SET_NAME  => 'VISTA enhancer set',
};


my @experimentally_verified_enhancer_feature_set_names = (
  FANTOM_FEATURE_SET_NAME,
  VISTA_FEATURE_SET_NAME,
);

my %feature_set_to_file_name = (
  &FANTOM_FEATURE_SET_NAME => 'fantom.bed',
  &VISTA_FEATURE_SET_NAME  => 'vista_enhancers.bed',
);

use Hash::Util qw( lock_keys );
lock_keys(%feature_set_to_file_name);

my @experimentally_verified_enhancer_feature_sets
  = grep { 
    defined $_
  } map {
    $feature_set_adaptor->fetch_by_name($_);
  } @experimentally_verified_enhancer_feature_set_names;

my $slice_adaptor = $core_adaptor->get_SliceAdaptor;
my $karyotype_slice = $slice_adaptor->fetch_all_karyotype;

use Bio::EnsEMBL::Funcgen::Utils::RegulatoryBuildStatUtils qw ( 
  REGULATORY_FEATURE_TYPES
  ENHANCER
);

foreach my $feature_set (@experimentally_verified_enhancer_feature_sets) {
    $logger->info("feature_set: " . $feature_set->name . "\n");
}

$logger->info("Checking enhancers\n");

foreach my $experimentally_verified_enhancer_feature_set (@experimentally_verified_enhancer_feature_sets) {

    my $current_feature_set = $experimentally_verified_enhancer_feature_set->name;

    my $total_enhancers_checked       = 0;
    my $num_enhancers_overlapping     = 0;
    my $num_enhancers_not_overlapping = 0;

    my $file_name = $tempdir . '/' . $feature_set_to_file_name{$current_feature_set};

    $logger->info("Writing $current_feature_set to $file_name\n");

    open 
        my $features_fh, 
        '>', 
        $file_name 
        || die("Can't open $file_name!");

    $logger->info("Iterating over karyotype slices\n");
    
    my $number_in_total = 0;

    foreach my $current_karyotype_slice (@$karyotype_slice) {

        my $seq_region_name = $current_karyotype_slice->seq_region_name;
        my $length          = $current_karyotype_slice->length;
        
        $logger->info("    Exporting $current_feature_set features on chromosome $seq_region_name");
        
        my $feature_iterator = $external_feature_adaptor
        ->fetch_Iterator_by_Slice_FeatureSets(
            $current_karyotype_slice, 
            [ $experimentally_verified_enhancer_feature_set ]
        );
        
        my $number_on_slice = 0;
        
        FEATURE:
        while ($feature_iterator->has_next) {
        
            # This may or may not actually be an enhancer
            my $current_enhancer_feature = $feature_iterator->next;
            
            # Check it really is an enhancer, sadly enhancers are hodgepodged together
            # with other stuff in the same feature set.
            #
            if ($current_enhancer_feature->feature_type->class ne ENHANCER) {
                next FEATURE;
            }
            # When checking FANTOM, only use either permissive or robust otherwise there will be duplicate enhancers.
            if (
              $current_feature_set eq FANTOM_FEATURE_SET_NAME 
              #&& $current_enhancer_feature->feature_type->name ne "FANTOM permissive enhancer"
              && $current_enhancer_feature->feature_type->name ne "FANTOM robust enhancer"
            ) {
                next FEATURE;
            }
            $number_on_slice++;
            $number_in_total++;
            
            my $bed_line = 
                join 
                    "\t",
                    $seq_region_name, 
                    $current_enhancer_feature->start, 
                    $current_enhancer_feature->end;

            $features_fh->print($bed_line);
            $features_fh->print("\n");
        }
        $logger->info(" got: $number_on_slice\n");
    }
    $features_fh->close;
    $logger->info("    Got a total of $number_in_total $current_feature_set features.\n");
    $logger->info("    Sorting $file_name\n");
    
    sort_bed_file($file_name);
    
    #my $count_overlaps_cmd = "bedtools intersect -a $regulatory_features_file -b $file_name | wc -l";
    
    # Some features may overlap more than one regulatory feature
    #my $count_overlaps_cmd = "bedtools intersect -loj -a $file_name -b $regulatory_features_file | cut -f 6 | grep -v '\\\-1' | sort | uniq | wc -l";
    
    my $count_overlaps_cmd = "bedtools intersect -u -a $file_name -b $regulatory_features_file | sort | uniq | wc -l";
    
    
    $logger->info("    Computing and counting overlaps using this command: $count_overlaps_cmd\n");
    
    my $num_overlaps = `$count_overlaps_cmd`;
    chomp($num_overlaps);
    
    use Scalar::Util qw( looks_like_number );
    
    if (! looks_like_number($num_overlaps)) {
        die("num_overlaps $num_overlaps doesn't look like a number!");
    }
    
    $logger->info("    $current_feature_set features have $num_overlaps overlaps ($number_in_total total) with regulatory features.\n");
    
    use Bio::EnsEMBL::Funcgen::RegulatoryBuildStatistic;
    if ($current_feature_set eq VISTA_FEATURE_SET_NAME) {
        $regulatory_build_statistic_adaptor->store(
            Bio::EnsEMBL::Funcgen::RegulatoryBuildStatistic->new(
                -statistic           => "total_enhancers_checked_vista",
                -value               => $number_in_total,
                -regulatory_build_id => $regulatory_build->dbID,
            )
        );
        $regulatory_build_statistic_adaptor->store(
            Bio::EnsEMBL::Funcgen::RegulatoryBuildStatistic->new(
                -statistic           => "num_enhancers_overlapping_vista",
                -value               => $num_overlaps,
                -regulatory_build_id => $regulatory_build->dbID,
            )
        );
    }
    if ($current_feature_set eq FANTOM_FEATURE_SET_NAME) {
        $regulatory_build_statistic_adaptor->store(
            Bio::EnsEMBL::Funcgen::RegulatoryBuildStatistic->new(
                -statistic           => "total_enhancers_checked_fantom",
                -value               => $number_in_total,
                -regulatory_build_id => $regulatory_build->dbID,
            )
        );
        $regulatory_build_statistic_adaptor->store(
            Bio::EnsEMBL::Funcgen::RegulatoryBuildStatistic->new(
                -statistic           => "num_enhancers_overlapping_fantom",
                -value               => $num_overlaps,
                -regulatory_build_id => $regulatory_build->dbID,
            )
        );
    }
}

$logger->info("Done checking enhancers.\n");

$logger->finish_log;
exit(0);

sub export_regulatory_features {

    my $regulatory_features_file = shift;

    my $max = 20;
    my $i = 0;
    my $stop_after_maximum_reached = 0;
    
    my $iterator = $regulatory_feature_adaptor->fetch_Iterator_by_RegulatoryBuild($regulatory_build);

    open 
        my $regulatory_features_fh, 
        '>', 
        $regulatory_features_file 
        || die("Can't open $regulatory_features_file!");

    REGULATORY_FEATURE:
    while ($iterator->has_next) {

        my $current_regulatory_feature = $iterator->next;
        
        my $bed_line = join "\t", (
            $current_regulatory_feature->slice->seq_region_name,
            $current_regulatory_feature->bound_seq_region_start,
            $current_regulatory_feature->bound_seq_region_end,
        );
        
        $regulatory_features_fh->print($bed_line);
        $regulatory_features_fh->print("\n");
        
        $i++;
        if ($stop_after_maximum_reached && $i == $max) {
            last REGULATORY_FEATURE
        }
    }
    $regulatory_features_fh->close;
}

sub sort_bed_file {
    my $bed_file = shift;
    system("bedSort $bed_file $bed_file");
}
