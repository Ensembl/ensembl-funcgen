#!/usr/bin/env perl

use strict;
use Bio::EnsEMBL::Registry;
use Data::Dumper;
use Carp;
use Bio::EnsEMBL::Utils::Logger;

use strict;
use Getopt::Long;

my $species;
my $registry;
my $dry_run;
my $partition_by_experimental_group;

use Bio::EnsEMBL::Funcgen::PeakCallingPlan::Constants qw ( :all );

=head1

classify_epigenome_to_segmentation_run.pl \
    --species  homo_sapiens \
    --registry /homes/mnuhn/work_dir_regbuild_testrun/lib/ensembl-funcgen/registry.with_previous_version.human_regbuild_testdb6.pm \
    --dry_run 1 \
    --partition_by_experimental_group 0

=cut

GetOptions (
   'species=s'    => \$species,
   'registry=s'   => \$registry,
   'dry_run=s'    => \$dry_run,
   'partition_by_experimental_group=s' => \$partition_by_experimental_group,
);

my $logger = Bio::EnsEMBL::Utils::Logger->new();
$logger->init_log;

$logger->info("registry   = " . $registry   . "\n");
$logger->info("species    = " . $species    . "\n");
$logger->info("dry_run    = " . $dry_run    . "\n");
$logger->info(
    "partition_by_experimental_group = " 
    . $partition_by_experimental_group . "\n"
);

use Bio::EnsEMBL::Registry;
Bio::EnsEMBL::Registry->load_all($registry);

my $funcgen_dba = Bio::EnsEMBL::Registry->get_DBAdaptor($species, 'funcgen');
my $epigenome_adaptor          = $funcgen_dba->get_EpigenomeAdaptor;
my $experiment_adaptor         = $funcgen_dba->get_ExperimentAdaptor;
my $experimental_group_adaptor = $funcgen_dba->get_ExperimentalGroupAdaptor;
my $peak_calling_adaptor       = $funcgen_dba->get_PeakCallingAdaptor;
my $feature_type_adaptor       = $funcgen_dba->get_FeatureTypeAdaptor;

my $all_epigenomes = $epigenome_adaptor->fetch_all;

=head2

drop table if exists segmentation_cell_tables;
CREATE TABLE if not exists segmentation_cell_tables (
  superclass       varchar(255) NOT NULL,
  class            varchar(255) NOT NULL,
  epigenome        varchar(255) NOT NULL,
  feature_type     varchar(255) NOT NULL,
  signal_bam_path  varchar(255) NOT NULL,
  control_bam_path varchar(255) NOT NULL
) ENGINE=InnoDB DEFAULT CHARSET=latin1;

=cut

my $superclass = 'default';

my $segmentation_feature_types = [
    "CTCF",
    "DNase1",
    "H3K4me1", 
    "H3K4me2", 
    "H3K4me3", 
    "H3K9ac", 
    "H3K9me3",
    "H3K27ac", 
    "H3K27me3", 
    "H3K36me3", 
];

my $segmentation_feature_type_objects 
    = [ 
        map { 
            $feature_type_adaptor->fetch_by_name($_) 
        } @$segmentation_feature_types 
    ];

my $experimental_groups = $experimental_group_adaptor->fetch_all;

foreach my $experimental_group (@$experimental_groups) {

    $logger->info("Classifying all epigenomes from " . $experimental_group->name . "\n");

    my %class_to_epigenome_name = classify_epigenomes(
        $all_epigenomes, 
        $experimental_group
    );

    use Data::Dumper;
    print Dumper(\%class_to_epigenome_name);
    
    $logger->info("\tCreating segmentation cell table\n");

    create_segmentation_cell_table(
        \%class_to_epigenome_name,
        generate_sql_callback({ 
            dry_run => $dry_run 
        }),
        $experimental_group
    );
}

$logger->finish_log;

sub generate_sql_callback {

    my $param = shift;
    
    my $dry_run = $param->{dry_run};

    my $callback = sub {
    
        my $param = shift;
        
        my $class                      = $param->{class};
        my $epigenome                  = $param->{epigenome};
        my $segmentation_feature_type  = $param->{segmentation_feature_type};
        my $signal_bam_data_file       = $param->{signal_bam_data_file};
        my $control_bam_data_file      = $param->{control_bam_data_file};
        my $experimental_group         = $param->{experimental_group};

        my $epigenome_production_name      = $epigenome->production_name;
        my $segmentation_feature_type_name = $segmentation_feature_type->name;
        my $signal_bam_data_file_path      = $signal_bam_data_file->path;

        my $control_bam_data_file_path;
        if (defined $control_bam_data_file) {
            $control_bam_data_file_path = $control_bam_data_file->path;
        }
        
        my $superclass = 'default';
        
        if ($partition_by_experimental_group) {
            $superclass = $experimental_group->production_name;
        }
        
        my @values_to_insert = (
            $superclass,
            $class,
            $epigenome_production_name,
            $segmentation_feature_type_name,
            $signal_bam_data_file_path,
            $control_bam_data_file_path
        );
        
        my $sql_cmd = "insert into segmentation_cell_tables values ("
            . (
                join ", ", map { '"' . $_ . '"' } @values_to_insert
            )
            . ");";
        if ($dry_run) {
        
            $logger->info(
                "\t" . (join "\t", @values_to_insert) . "\n"
            );
            
        } else {
            $logger->info("\t$sql_cmd\n");
            $funcgen_dba->dbc->do($sql_cmd);
        }
   };
   return $callback;
}

sub create_segmentation_cell_table {

    my $class_to_epigenome_name        = shift;
    my $build_row_callback             = shift;
    my $restrict_to_experimental_group = shift;

    CLASS:
    foreach my $current_class (sort keys %$class_to_epigenome_name) {
        
        if ($current_class eq SEGMENTATION_CLASS_NONE) {
            next CLASS;
        }

        #print "\t$current_class\n";
        
        my $epigenomes_in_class = $class_to_epigenome_name->{$current_class};
        
        EPIGENOME:
        foreach my $epigenome_production_name (@$epigenomes_in_class) {
        
            if (! defined $epigenome_production_name) {
                print "No epigenome for segmentation class $current_class\n";
                next EPIGENOME;
            }
        
            my $epigenome = $epigenome_adaptor->fetch_by_production_name($epigenome_production_name);
            
            if (! defined $epigenome) {
                die("Couldn't find epigenome with name $epigenome_production_name!");
            }
            
            SEGMENTATION_FEATURE_TYPE:
            foreach my $segmentation_feature_type (@$segmentation_feature_type_objects) {
            
                my $peak_callings = $peak_calling_adaptor
                    ->fetch_all_by_Epigenome_FeatureType(
                        $epigenome, 
                        $segmentation_feature_type
                    );
                if (@$peak_callings == 0) {
                    #print "No peak calling for ".$epigenome -> name.", ".$segmentation_feature_type->name."\n";
                    next SEGMENTATION_FEATURE_TYPE;
                }
                
                PEAK_CALLING:
                foreach my $peak_calling (@$peak_callings) {
                    my $signal_alignment  = $peak_calling->fetch_signal_Alignment;
                    my $control_alignment = $peak_calling->fetch_control_Alignment;
                    
                    my $control_bam_data_file;
                    
                    if (defined $control_alignment) {
                        $control_bam_data_file = $control_alignment->fetch_bam_DataFile;
                    }
                    
                    my $signal_bam_data_file  = $signal_alignment->fetch_bam_DataFile;
                    
                    my $experiment = $peak_calling->fetch_Experiment;
                    my $experimental_group = $experiment->get_ExperimentalGroup;
                    
                    # Include, if the experiment matches the one that is 
                    # being restricted to or if there is no such restriction.
                    #
                    my $include_this = 
                        (
                            (defined $restrict_to_experimental_group) 
                            && ($experimental_group->name eq $restrict_to_experimental_group->name)
                        )
                        || ( ! defined $restrict_to_experimental_group );
                    
                    if (! $include_this) {
                        next PEAK_CALLING;
                    }
                    
                    $build_row_callback->({
                        class                      => $current_class,
                        epigenome                  => $epigenome,
                        segmentation_feature_type  => $segmentation_feature_type,
                        signal_bam_data_file       => $signal_bam_data_file,
                        control_bam_data_file      => $control_bam_data_file,
                        experimental_group         => $experimental_group,
                    });
                }
            }
        }
    }
}

sub classify_epigenomes {

    my $epigenomes = shift;
    my $experimental_group = shift;

    my %class_to_epigenome_name;
    
    EPIGENOME:
    foreach my $current_epigenome (@$epigenomes) {

        my $class = classify_epigenome($current_epigenome, $experimental_group);

        if ($class eq SEGMENTATION_CLASS_NONE) {
            next EPIGENOME;
        }
        if (! exists $class_to_epigenome_name{$class}) {
            $class_to_epigenome_name{$class} = [];
        }
        push @{$class_to_epigenome_name{$class}}, $current_epigenome->production_name;
    }
    return %class_to_epigenome_name;
}

sub classify_epigenome {

    my $epigenome          = shift;
    my $experimental_group = shift;
    
    my $feature_types = fetch_feature_types_for_epigenome($epigenome, $experimental_group);
    
    my $has_minimum_required_feature_types_for_segmentation 
        = has_minimum_required_feature_types_for_segmentation($feature_types);

    my $has_CTCF                = has_CTCF               ($feature_types);
    my $has_H3K27ac             = has_H3K27ac            ($feature_types);
    my $has_H3K4me3_and_H3K4me1 = has_H3K4me3_and_H3K4me1($feature_types);
    
    # H3K27ac can be missing as long as H3K4me3 and H3K4me1 are available.
    
    if (! $has_minimum_required_feature_types_for_segmentation) {
        return SEGMENTATION_CLASS_NONE;
    }
    if (! $has_H3K4me3_and_H3K4me1) {
        return SEGMENTATION_CLASS_NONE;
    }

    if (  $has_H3K27ac &&   $has_CTCF) { return SEGMENTATION_CLASS_CTCF               }
    if (  $has_H3K27ac && ! $has_CTCF) { return SEGMENTATION_CLASS_NO_CTCF            }
    if (! $has_H3K27ac &&   $has_CTCF) { return SEGMENTATION_CLASS_CTCF_NO_H3K27AC    }
    if (! $has_H3K27ac && ! $has_CTCF) { return SEGMENTATION_CLASS_NO_CTCF_NO_H3K27AC }
    
    return SEGMENTATION_CLASS_NONE;
}

sub fetch_feature_types_for_epigenome {

    my $epigenome = shift;
    my $experimental_group_name = shift;
    
    my $experiments;
    
    if ($experimental_group_name) {
        $experiments = $experiment_adaptor->_fetch_all_by_Epigenome_ExperimentalGroup(
            $epigenome, 
            $experimental_group_name
        );
     } else {
        $experiments = $experiment_adaptor->fetch_all_by_Epigenome($epigenome);
     }
    
    
    my @feature_types = map { $_->feature_type } @$experiments;
    return \@feature_types;
}

sub has_CTCF {

    my $feature_types = shift;

    return has_required_feature_types({
        required => [ 'CTCF' ],
        present  => $feature_types,
    });
}

sub has_H3K27ac {

    my $feature_types = shift;

    return has_required_feature_types({
        required => [ 'H3K27ac' ],
        present  => $feature_types,
    });
}

sub has_H3K4me3_and_H3K4me1 {

    my $feature_types = shift;

    return has_required_feature_types({
        required => [ 'H3K4me3', 'H3K4me1' ],
        present  => $feature_types,
    });
}

sub has_minimum_required_feature_types_for_segmentation {

    my $feature_types = shift;
    
    my @minimum_required_feature_types_for_segmentation = qw(
        H3K4me1 H3K4me3 H3K27me3 H3K36me3
    );
    
    return has_required_feature_types({
        required => \@minimum_required_feature_types_for_segmentation,
        present  => $feature_types,
    });
}

sub has_required_feature_types {

    my $param = shift;
    
    my $required_feature_types = $param->{required};
    my $present_feature_types  = $param->{present};
    
    my $has_all_required_feature_types = 1;
    
    REQUIRED_FEATURE_TYPE:
    foreach my $required_feature_type (@$required_feature_types) {
    
        my $is_present = grep { 
            $_->name eq $required_feature_type 
        } @$present_feature_types;
        
        if (! $is_present) {
            $has_all_required_feature_types = 0;
            last REQUIRED_FEATURE_TYPE;
        }
        next REQUIRED_FEATURE_TYPE;
    }
    return $has_all_required_feature_types;
}

