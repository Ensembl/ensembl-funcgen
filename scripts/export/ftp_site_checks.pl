#!/usr/bin/env perl

use strict;
use Bio::EnsEMBL::Registry;
use File::Spec;
use Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Hive::DBSQL::DBConnection;
use Getopt::Long;
use Data::Dumper;

=head2

ftp_site_checks.pl \
  --registry /homes/mnuhn/work_dir_regbuild_testrun/lib/ensembl-funcgen/registry.with_previous_version.testdb12.pm \
  --species mus_musculus

ftp_site_checks.pl \
  --registry /homes/mnuhn/work_dir_regbuild_testrun/lib/ensembl-funcgen/registry.ftp_grch37.pm \
  --species homo_sapiens \
  --ftp_dir /nfs/production/panda/ensembl/funcgen/mnuhn/human_grch37_fixed_ftp_regulatory_features

=cut


my $registry;
my $species;
my $ftp_dir;

GetOptions (
   'registry=s' => \$registry,
   'species=s'  => \$species,
   'ftp_dir=s'  => \$ftp_dir,
);

Bio::EnsEMBL::Registry->load_all($registry);

use Bio::EnsEMBL::Utils::Logger;
my $logger = Bio::EnsEMBL::Utils::Logger->new();

my $funcgen_adaptor = Bio::EnsEMBL::Registry->get_DBAdaptor( $species, 'funcgen');
my $dbc = $funcgen_adaptor->dbc;

check_regulatory_activity_directories_exist($dbc);

check_regulatory_activity_files_exist($dbc);

check_regulatory_feature_numbers({
    dbc     => $dbc,
    ftp_dir => $ftp_dir,
    species => $species,
});

check_file_content_consistent({
    ftp_dir => $ftp_dir,
    species => $species,
    funcgen_adaptor => $funcgen_adaptor,
});

check_activity_summaries_consistent({
    dbc     => $dbc,
    ftp_dir => $ftp_dir,
    species => $species,
});

$logger->finish_log;

exit;

sub create_epigenome_cache {

    my $funcgen_adaptor = shift;
    my $epigenome_adaptor = $funcgen_adaptor->get_EpigenomeAdaptor;
    my $all_epigenomes = $epigenome_adaptor->fetch_all;
    
    my $epigenome_cache = { 
        map  { $_->production_name => $_ } 
        grep { $_->production_name } 
        @$all_epigenomes 
    };
    return $epigenome_cache;
}

sub check_file_content_consistent {

    my $param = shift;
    
    my $funcgen_adaptor = $param->{funcgen_adaptor};
    my $ftp_dir = $param->{ftp_dir};
    my $species = $param->{species};
    
    my $dbc = $funcgen_adaptor->dbc;
    
    my $epigenome_cache = create_epigenome_cache($funcgen_adaptor);
    
    my $regulatory_activity_file_names = get_expected_regulatory_activity_file_names({
        dbc     => $dbc,
        ftp_dir => $ftp_dir,
    });
    
    my $regulatory_feature_adaptor = $funcgen_adaptor->get_RegulatoryFeatureAdaptor;
    
    foreach my $regulatory_activity_file_name (@$regulatory_activity_file_names) {
    
        my $file_name = $regulatory_activity_file_name->{file_name};
        my $epigenome_production_name
            = $regulatory_activity_file_name->{epigenome_production_name};
        
        my $epigenome = $epigenome_cache->{$epigenome_production_name};
        
        if (! defined $epigenome) {
            die;
        }
        
        open my $fh, "zcat $file_name |" or die($file_name);
        
        while (my $current_line = <$fh>) {
            chomp($current_line);
            
            my @f = split "\t", $current_line;
            
            my $attributes = $f[8];
            
            # Attributes look like this:
            #
            # activity=NA;bound_end=102119230;bound_start=102118695;description=TF binding site na in Lung;epigenome=Lung;feature_type=TF binding site;regulatory_feature_stable_id=ENSR00000368862
            #
            my %attribute_hash = map { split /=/ } split ';', $attributes;
            
            if (! exists $attribute_hash{regulatory_feature_stable_id}) {
                die("Couldn't find stable id in: $attributes");
            }
            my $stable_id = $attribute_hash{regulatory_feature_stable_id};
            
            my $regulatory_feature_from_file = {
                seq_region_name => $f[0],
                start           => $f[3],
                end             => $f[4],
                %attribute_hash
            };
            
            my $regulatory_feature_from_db = $regulatory_feature_adaptor->fetch_by_stable_id($regulatory_feature_from_file->{regulatory_feature_stable_id});
            
            my $start_ok        = $regulatory_feature_from_file->{start}           == $regulatory_feature_from_db->{start};
            my $end_ok          = $regulatory_feature_from_file->{end}             == $regulatory_feature_from_db->{end};
            my $feature_type_ok = $regulatory_feature_from_file->{feature_type}    eq $regulatory_feature_from_db->feature_type->name;
            my $seq_region_ok   = $regulatory_feature_from_file->{seq_region_name} eq $regulatory_feature_from_db->slice->seq_region_name;
            
            my $regulatory_activity_from_db   = $regulatory_feature_from_db->regulatory_activity_for_epigenome($epigenome)->activity;
            my $regulatory_activity_from_file = $regulatory_feature_from_file->{activity};
            
            my $activity_ok = $regulatory_feature_from_file->{activity} eq $regulatory_activity_from_db;
            
            if (! $seq_region_ok) {
                die("Seq regions don't agree in $file_name! File: ".$regulatory_feature_from_file->{seq_region_name}." DB:" . $regulatory_feature_from_db->slice->name);
            }
            if (! $activity_ok) {
                die("Regulatory activities don't agree: File: $regulatory_activity_from_file DB: $regulatory_activity_from_db!");
            }
            if (! $feature_type_ok) {
                die("Feature types don't agree in $file_name File: ".$regulatory_feature_from_file->{feature_type}." DB: ".$regulatory_feature_from_db->feature_type->name."!");
            }
            if (! ($start_ok && $end_ok)) {
                die("Coordinates don't agree in $file_name!");
            }
        }
        
        close($fh);
        $logger->info("Ok - File content for $file_name\n");
    }
    return;
}

sub fetch_activity_summary_from_db {

    my $param = shift;
    
    my $dbc = $param->{dbc};
    my $epigenome_production_name = $param->{epigenome_production_name};

    my $sql = "
        select 
            activity, 
            count(regulatory_activity_id) as number
        from 
            regulatory_activity 
            join epigenome using (epigenome_id) 
        where 
            epigenome.production_name = '$epigenome_production_name'
        group by 
            activity
    ";

    my $sth = $dbc->prepare($sql);
    $sth->execute;
    #my $result = $sth->fetchall_hashref('activity');
    my $result = $sth->fetchall_arrayref();
    # $VAR1 = [
    #           [
    #             'INACTIVE',
    #             107014
    #           ],
    #           [
    #             'REPRESSED',
    #             1621
    #           ],
    #           [
    #             'POISED',
    #             13030
    #           ],
    #           [
    #             'ACTIVE',
    #             25195
    #           ],
    #           [
    #             'NA',
    #             133167
    #           ]
    #         ];

    # Hashify
    my $activity_summary_from_db = { map { ($_->[0] => $_->[1])  } @$result };
    # $VAR1 = {
    #           'POISED' => 13030,
    #           'INACTIVE' => 107014,
    #           'ACTIVE' => 25195,
    #           'REPRESSED' => 1621,
    #           'NA' => 133167
    #         };

    return $activity_summary_from_db;
}

sub check_activity_summaries_consistent {

    my $param = shift;
    
    my $dbc     = $param->{dbc};
    my $ftp_dir = $param->{ftp_dir};
    my $species = $param->{species};
    
    my $regulatory_activity_file_names = get_expected_regulatory_activity_file_names({
        dbc     => $dbc,
        ftp_dir => $ftp_dir,
    });
    
    my @valid_activities = ('INACTIVE', 'REPRESSED', 'POISED', 'ACTIVE', 'NA');
    
    foreach my $regulatory_activity_file_name (@$regulatory_activity_file_names) {
    
        my $file_name = $regulatory_activity_file_name->{file_name};
        my $epigenome_production_name
            = $regulatory_activity_file_name->{epigenome_production_name};
        
        my $activities_from_file = {};
        
        foreach my $activity (@valid_activities) {
            my $number_of_na_regulatory_features_from_file
                = fetch_number_of_regulatory_features_with_activity_from_file(
                    $file_name,
                    $activity
                );
            $activities_from_file->{$activity} = $number_of_na_regulatory_features_from_file;
        }

        my $activity_summary_from_db = fetch_activity_summary_from_db({
            dbc     => $dbc,
            epigenome_production_name => $epigenome_production_name,
        });
        
        foreach my $activity (@valid_activities) {
        
            my $counts_match = $activities_from_file->{$activity} == $activity_summary_from_db->{$activity};
            
            if ($counts_match) {
                $logger->info("Ok - numbers identical in $epigenome_production_name for $activity ".$activities_from_file->{$activity}."\n");
            } else {
                $logger->error("Not ok - numbers not identical in $epigenome_production_name for $activity\n");
            }
        }
    }
    return;
}

sub fetch_number_of_regulatory_features_with_activity_from_file {

    my $file_name = shift;
    my $activity  = shift;

    my $cmd = "zcat $file_name | grep 'activity=$activity;' | wc -l";
    #die($cmd);
    my $number_of_regulatory_features = `$cmd`;
    chomp $number_of_regulatory_features;
    return 0 + $number_of_regulatory_features;
}

sub check_regulatory_feature_numbers {

    my $param = shift;
    
    my $dbc     = $param->{dbc};
    my $ftp_dir = $param->{ftp_dir};
    my $species = $param->{species};
    
    my $file = locate_regulatory_build_file($ftp_dir, $species);

    my $number_of_regulatory_features_from_db   = fetch_number_of_regulatory_features_from_db($dbc);
    my $number_of_regulatory_features_from_file = fetch_number_of_regulatory_features_from_file($file);

    my $numbers_identical   = $number_of_regulatory_features_from_db == $number_of_regulatory_features_from_file;
    my $number_greater_zero = $number_of_regulatory_features_from_db > 0;

    if ($numbers_identical) {
        $logger->info("Ok - numbers identical\n");
    } else {
        $logger->error("Not ok - numbers not identical\n");
    }

    if ($number_greater_zero) {
        $logger->info("Ok - have regulatory features\n");
    } else {
        $logger->error("Not ok - no regulatory features\n");
    }
    return;
}
sub check_regulatory_activity_files_exist {

    my $dbc = shift;
    my $epigenome_production_names = fetch_epigenome_production_names($dbc);
    
    my $expected_regulatory_activity_file_names = get_expected_regulatory_activity_file_names({
        dbc     => $dbc,
        ftp_dir => $ftp_dir,
    });
    
    foreach my $expected_file (@$expected_regulatory_activity_file_names) {

        my $file_name                 = $expected_file->{file_name};
        my $pattern                   = $expected_file->{pattern};
        my $epigenome_production_name = $expected_file->{epigenome_production_name};
    
        my $file_exists = -e $file_name;
        
        if ($file_exists) {
            $logger->info("Ok - file exists: $file_name\n");
        } else {
            $logger->error("Not ok - file for $epigenome_production_name ($file_name) not found using pattern $pattern\n");
        }
    }
    return;
}

sub check_regulatory_activity_directories_exist {

    my $dbc = shift;
    my $epigenome_production_names = fetch_epigenome_production_names($dbc);
    
    my $regulatory_activity_directory_names = get_expected_regulatory_activity_directory_names({
        dbc     => $dbc,
        ftp_dir => $ftp_dir,
    });
    
    foreach my $expected_directory (@$regulatory_activity_directory_names) {
    
        my $dir_name = $expected_directory->{dir_name};
        my $epigenome_production_name = $expected_directory->{epigenome_production_name};
        
        my $directory_exists = -d $dir_name;
        if ($directory_exists) {
            $logger->info("Ok - directory exists: $dir_name\n");
        } else {
            $logger->error("Not ok - directory for " . $epigenome_production_name . " missing: $dir_name\n");
        }
    }
    return;
}

sub get_expected_regulatory_activity_file_names {

    my $param = shift;
    
    my $dbc     = $param->{dbc};
    my $ftp_dir = $param->{ftp_dir};

    my $epigenome_production_names = fetch_epigenome_production_names($dbc);
    my $regulatory_activity_directory_names 
        = get_expected_regulatory_activity_directory_names($param);

    my $expected_regulatory_activity_directory_names = [
        map
            {
                $ftp_dir . '/RegulatoryFeatureActivity/' . $_->{dir_name}
            } @$regulatory_activity_directory_names
    ];
    
    my @regulatory_activity_file_names;
    foreach my $directory (@$regulatory_activity_directory_names) {
    
        my $dir_name = $directory->{dir_name};
        my $epigenome_production_name = $directory->{epigenome_production_name};
    
        # E.g.:
        # homo_sapiens.GRCh37.Fetal_Adrenal_Gland.Regulatory_Build.regulatory_activity.20180503.gff.gz
        my $pattern 
            = $dir_name
                . '/' . $species 
                . '.*.'
                . $epigenome_production_name
                . '.Regulatory_Build.regulatory_activity.*.gff.gz'
        ;
        my @expected_files = glob $pattern;
        my $expected_file  = $expected_files[0];
        push 
            @regulatory_activity_file_names, 
            {
                file_name                 => $expected_file,
                pattern                   => $pattern,
                epigenome_production_name => $epigenome_production_name,
            }
        ;
    }
    return \@regulatory_activity_file_names;
}

sub get_expected_regulatory_activity_directory_names {

    my $param = shift;
    
    my $dbc     = $param->{dbc};
    my $ftp_dir = $param->{ftp_dir};

    my $epigenome_production_names = fetch_epigenome_production_names($dbc);

    my $expected_regulatory_activity_directory_names = [
        map
            {
                {
                    dir_name 
                        => build_regulatory_activity_directory_name({
                            dbc                       => $dbc,
                            ftp_dir                   => $ftp_dir,
                            epigenome_production_name => $_,
                        }),
                    epigenome_production_name => $_
                }
            } @$epigenome_production_names
    ];
    return $expected_regulatory_activity_directory_names;
}

sub build_regulatory_activity_directory_name {

    my $param = shift;
    
    my $dbc     = $param->{dbc};
    my $ftp_dir = $param->{ftp_dir};
    my $epigenome_production_name = $param->{epigenome_production_name};

    my $regulatory_activity_directory_name
        = $ftp_dir 
            . '/RegulatoryFeatureActivity/' 
            . $epigenome_production_name
    ;
    return $regulatory_activity_directory_name;
}

sub locate_regulatory_build_file {

    my $ftp_dir = shift;
    my $species = shift;
    
    # E.g.: /nfs/production/panda/ensembl/funcgen/mnuhn/human_grch37_fixed_ftp_regulatory_features/homo_sapiens.GRCh37.Regulatory_Build.regulatory_features.20180503.gff.gz
    #
    my $file = glob $ftp_dir . '/' . $species . '*.Regulatory_Build.regulatory_features.*.gff.gz';

    if (! -e $file) {
        die;
    }
    return $file;
}

sub fetch_epigenome_production_names {

    my $dbc = shift;

    my $sql = 'select production_name from epigenome join regulatory_build_epigenome using (epigenome_id)';

    my $sth = $dbc->prepare($sql);
    $sth->execute;
    my $result = $sth->fetchall_arrayref;
    my $list_of_names = [ map { $_->[0] } @$result ];
    
    return $list_of_names;
}

sub fetch_number_of_regulatory_features_from_db {

    my $dbc = shift;

    my $sql = 'select count(*) num from regulatory_feature';

    my $sth = $dbc->prepare($sql);
    $sth->execute;
    my $result = $sth->fetchall_arrayref;
    my $number_of_regulatory_features = $result->[0]->[0];
    
    return $number_of_regulatory_features;
}

sub fetch_number_of_regulatory_features_from_file {

    my $file_name = shift;

    my $number_of_regulatory_features = `zcat $file_name | wc -l`;
    chomp $number_of_regulatory_features;
    return 0 + $number_of_regulatory_features;
}
