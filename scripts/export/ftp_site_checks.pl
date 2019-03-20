#!/usr/bin/env perl

use strict;
use Bio::EnsEMBL::Registry;
use File::Spec;
use Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Hive::DBSQL::DBConnection;
use Getopt::Long;
use Data::Dumper;
use Carp;

use Bio::EnsEMBL::Funcgen::Utils::FtpUtils qw (
  check_file_content_consistent 
  get_all_regulatory_activity_file_infos_from_ftp
  fetch_epigenome_production_names
  get_expected_regulatory_activity_directory_names
);

=head2

ftp_site_checks.pl \
  --registry /homes/mnuhn/work_dir_regbuild_testrun/lib/ensembl-funcgen/registry.with_previous_version.human_regbuild_testdb22.pm \
  --species homo_sapiens \
  --ftp_dir /hps/nobackup/production/sds-flicek-bp/blueprint_fastq_files/mnuhn_regbuild_pipeline/rb_grch38_testdb22/ftp \
  --check check_activity_summaries_consistent

=cut

my $registry;
my $species;
my $ftp_dir;
my $check;
my $check_file;

GetOptions (
   'registry=s'   => \$registry,
   'species=s'    => \$species,
   'ftp_dir=s'    => \$ftp_dir,
   'check=s'      => \$check,
   'check_file=s' => \$check_file,
);

Bio::EnsEMBL::Registry->load_all($registry);

use Bio::EnsEMBL::Utils::Logger;
my $logger = Bio::EnsEMBL::Utils::Logger->new;

# Doesn't work when passed in constructor
$logger->{'loglevel'} = 4;

$logger->init_log;

$logger->info('registry   = ' . $registry   . "\n");
$logger->info('species    = ' . $species    . "\n");
$logger->info('ftp_dir    = ' . $ftp_dir    . "\n");
$logger->info('check      = ' . $check      . "\n");
$logger->info('check_file = ' . $check_file . "\n");

$logger->info('loglevel    = ' . $logger->{'loglevel'} . "\n");

my $funcgen_adaptor = Bio::EnsEMBL::Registry->get_DBAdaptor( $species, 'funcgen');
my $dbc = $funcgen_adaptor->dbc;

if (! defined $check || $check eq 'check_regulatory_activity_directories_exist') {

  $logger->info("Running check_regulatory_activity_directories_exist\n");
  check_regulatory_activity_directories_exist($dbc);
}

if (! defined $check || $check eq 'check_regulatory_activity_files_exist') {

  $logger->info("Running check_regulatory_activity_files_exist\n");
  check_regulatory_activity_files_exist($dbc);
}

if (! defined $check || $check eq 'check_regulatory_feature_numbers') {

  $logger->info("Running check_regulatory_feature_numbers\n");
  check_regulatory_feature_numbers({
      dbc     => $dbc,
      ftp_dir => $ftp_dir,
      species => $species,
  });
}

if (! defined $check || $check eq 'check_file_content_consistent') {

  my $regulatory_activity_file_name_info = get_all_regulatory_activity_file_infos_from_ftp({
      dbc     => $dbc,
      ftp_dir => $ftp_dir,
      species => $species,
  });
  
  my $check_these = $regulatory_activity_file_name_info;
  
  if (defined $check_file) {
    
    my @found_files = grep {
      $_->{file_name} eq $check_file
    } @$regulatory_activity_file_name_info;
    
    if (@found_files == 0) {
      confess("Can't find $check_file!");
    }
    
    if (@found_files > 1) {
      confess("Found more than one file matching $check_file!");
    }
    
    $check_these = [ $found_files[0] ];
  }
  
  #die(Dumper($regulatory_activity_file_name_info));
  #die(Dumper($check_these));

  $logger->info("Running check_file_content_consistent\n");
  check_file_content_consistent({
      ftp_dir => $ftp_dir,
      species => $species,
      funcgen_adaptor => $funcgen_adaptor,
      check_these => $check_these
  });
}

if (! defined $check || $check eq 'check_activity_summaries_consistent') {

  $logger->info("Running check_activity_summaries_consistent\n");
  check_activity_summaries_consistent({
      dbc     => $dbc,
      ftp_dir => $ftp_dir,
      species => $species,
  });
}

$logger->finish_log;

exit;

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
        species => $species,
    });
    
    my @valid_activities = ('INACTIVE', 'REPRESSED', 'POISED', 'ACTIVE', 'NA');
    
    foreach my $regulatory_activity_file_name (@$regulatory_activity_file_names) {
    
        $logger->debug("     Checking " . $regulatory_activity_file_name->{file_name} . "\n");
    
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
        
        my $passes = 1;
        
        foreach my $activity (@valid_activities) {
        
            $logger->debug("          file: " . $activities_from_file->{$activity} . " db: " . $activity_summary_from_db->{$activity} . "\n");
        
            my $counts_match = $activities_from_file->{$activity} == $activity_summary_from_db->{$activity};
            
            if ($counts_match) {
                $logger->info("Ok - numbers identical in $epigenome_production_name for $activity ".$activities_from_file->{$activity}."\n");
            } else {
                $logger->warning("Not ok - numbers not identical in $epigenome_production_name for $activity\n");
                $passes = 0;
            }
        }
        if (! $passes) {
          $logger->error("Test has failed, see errors above.\n");
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
        species => $species,
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
        
        if (@expected_files > 1) {
          $logger->error("Found more than one file in directory! " . Dumper(\@expected_files));
        }
        
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

sub locate_regulatory_build_file {

    my $ftp_dir = shift;
    my $species = shift;
    
    # E.g.: /nfs/production/panda/ensembl/funcgen/mnuhn/human_grch37_fixed_ftp_regulatory_features/homo_sapiens.GRCh37.Regulatory_Build.regulatory_features.20180503.gff.gz
    #
    my $file = glob $ftp_dir . '/' . $species . '*.Regulatory_Build.regulatory_features.*.gff.gz';
    
    if (! -e $file) {
        die "Can't find regulatory features file in $ftp_dir";
    }
    return $file;
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
