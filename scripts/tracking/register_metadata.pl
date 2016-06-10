#!/usr/bin/env perl

# Copyright [1999-2016] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

use strict;
use warnings;

# use diagnostics;
use autodie;
use feature qw(say);

use Cwd 'abs_path';

# use Data::Dumper;
use Config::Tiny;
use Bio::EnsEMBL::Utils::Logger;
use Bio::EnsEMBL::Utils::Exception qw(throw);
use Bio::EnsEMBL::OntologyXref;
use Getopt::Long;
use File::Basename;

# use DateTime;

#TODO Dry run implementation
#TODO Confirm object creation
#TODO POD for every subroutine
#TODO Usage subroutine
#TODO Logger
#TODO Register controls first
#TODO Input Subset exists error - throw exception
#TODO Experiment exists error - throw exception

#TODO Input subset tracking table: registration date?
#TODO EFO IDs
#TODO Partial registration warning

main();

sub main {
    my ( $csv, $config, $help, $dry );

    GetOptions(
        "i=s" => \$csv,
        "c=s" => \$config,
        "h"   => \$help,
        "n"   => \$dry,
    );

    usage() if $help;
    usage() unless $csv;
    usage() unless $config;

    my $cfg = Config::Tiny->new;
    $cfg = Config::Tiny->read($config);

    my $logger = Bio::EnsEMBL::Utils::Logger->new(

        # -LOGAUTO     => 1,
        # -LOGAUTOBASE => $cfg->{log}->{logautobase},
        # -LOGLEVEL    => $cfg->{log}->{loglevel},
    );
    $logger->init_log();

    $logger->info( "CSV file used: " . abs_path($csv) . "\n",       0, 1 );
    $logger->info( "Config file used: " . abs_path($config) . "\n", 0, 1 );

    $logger->info( 'Connecting to ' . $cfg->{dna_db}->{dbname} . '... ',
        0, 0 );

    my $dba                  = connect_to_trackingDB($cfg);
    my $tracking_db_adaptors = get_trackingDB_adaptors($cfg);
    my $db_entry_adaptor     = Bio::EnsEMBL::DBSQL::DBEntryAdaptor->new($dba);
    $logger->info( "done\n", 0, 1 );

    $logger->info( 'Reading metadata from ' . abs_path($csv) . '... ', 0, 0 );
    my $data = get_data_from_csv( abs_path($csv), $logger );
    my $entries = keys %{$data};
    $logger->info( "done\n",                0, 1 );
    $logger->info( $entries . "imported\n", 1, 1 );

    verify_analysis_featureType_expGroup( $data, $tracking_db_adaptors,
        $logger );

    for my $accession ( keys %{$data} ) {
        $logger->info( "Registering $accession\n", 0, 1 );

        my $entry = $data->$accession;
        my $objects;

        fetch_analysis( $entry, $tracking_db_adaptors, $objects );

        fetch_epigenome( $entry, $tracking_db_adaptors, $objects );

        fetch_ontology_xref( $entry, $db_entry_adaptor, $objects );

        fetch_feature_type( $entry, $tracking_db_adaptors, $objects );

        fetch_exp_group( $entry, $tracking_db_adaptors, $objects );

        fetch_experiment( $entry, $tracking_db_adaptors, $objects );

        fetch_input_subset();

        $logger->info( "Successful Registration\n", 0, 1 );

    }

}

sub connect_to_trackingDB {
    my ($cfg) = @_;

    my $dba = Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor->new(
        -user       => $cfg->{efg_db}->{user},
        -pass       => $cfg->{efg_db}->{pass},
        -host       => $cfg->{efg_db}->{host},
        -port       => $cfg->{efg_db}->{port},
        -dbname     => $cfg->{efg_db}->{dbname},
        -dnadb_name => $cfg->{dna_db}->{dbname},
    );
    $dba->dbc->do("SET sql_mode='traditional'");

    return \$dba;
}

sub get_trackingDB_adaptors {
    my ($cfg) = @_;
    my %tracking_db_adaptors;

    # Tracking DB hidden from user, hence no get_TrackingAdaptor method.
    # TrackingAdaptor->new() does not YET accept DBAdaptor object

    my $tracking_adaptor = Bio::EnsEMBL::Funcgen::DBSQL::TrackingAdaptor->new(
        -user       => $cfg->{efg_db}->{user},
        -pass       => $cfg->{efg_db}->{pass},
        -host       => $cfg->{efg_db}->{host},
        -port       => $cfg->{efg_db}->{port},
        -dbname     => $cfg->{efg_db}->{dbname},
        -species    => $cfg->{general}->{species},
        -dnadb_user => $cfg->{dna_db}->{user},
        -dnadb_pass => $cfg->{dna_db}->{pass},
        -dnadb_host => $cfg->{dna_db}->{host},
        -dnadb_port => $cfg->{dna_db}->{port},
        -dnadb_name => $cfg->{dna_db}->{dbname},
    );

    my $dba = $tracking_adaptor->db;

    $tracking_db_adaptors{ep}  = $dba->get_EpigenomeAdaptor();
    $tracking_db_adaptors{ft}  = $dba->get_FeatureTypeAdaptor();
    $tracking_db_adaptors{an}  = $dba->get_AnalysisAdaptor();
    $tracking_db_adaptors{eg}  = $dba->get_ExperimentalGroupAdaptor();
    $tracking_db_adaptors{ex}  = $dba->get_ExperimentAdaptor();
    $tracking_db_adaptors{iss} = $dba->get_InputSubsetAdaptor();
    $tracking_db_adaptors{rs}  = $dba->get_ResultSetAdaptor();
    $tracking_db_adaptors{rf}  = $dba->get_RegulatoryFeatureAdaptor();
    $tracking_db_adaptors{fs}  = $dba->get_FeatureSetAdaptor();
    $tracking_db_adaptors{ds}  = $dba->get_DataSetAdaptor();
    $tracking_db_adaptors{af}  = $dba->get_AnnotatedFeatureAdaptor();

    return 1;
}

sub get_data_from_csv {
    my ( $csv, $logger ) = @_;
    my %data;

    open my $csv_fh, '<', $csv;

    while ( readline $csv_fh ) {
        my ($accession, $epigenome,           $feature_type,
            $br,        $tr,                  $is_control,
            $md5,       $local_url,           $analysis,
            $exp_group, $ontology_accessions, $controlled_by,
        ) = split /\t/;

        $data{$accession} ||= {};

        $data{$accession}->{accession}           = $accession;
        $data{$accession}->{epigenome_name}      = $epigenome;
        $data{$accession}->{feature_type}        = $feature_type;
        $data{$accession}->{br}                  = $br;
        $data{$accession}->{tr}                  = $tr;
        $data{$accession}->{is_control}          = $is_control;
        $data{$accession}->{md5}                 = $md5;
        $data{$accession}->{local_url}           = $local_url;
        $data{$accession}->{analysis}            = $analysis;
        $data{$accession}->{exp_group}           = $exp_group;
        $data{$accession}->{ontology_accessions} = $ontology_accessions;
        $data{$accession}->{controlled_by}       = $controlled_by;

    }

    close $csv_fh;

    return \%data;
}

sub verify_analysis_featureType_expGroup {
    my ( $data, $tracking_db_adaptors, $logger ) = @_;
    my $abort = 0;
    my %to_register;

    for my $entry ( values %{$data} ) {

        my $analysis = $tracking_db_adaptors->{an}
            ->fetch_by_logic_name( $entry->{analysis} );

        my $feature_type = $tracking_db_adaptors->{ft}
            ->fetch_by_name( $entry->{feature_type} );

        my $exp_group = $tracking_db_adaptors->{eg}
            ->fetch_by_name( $entry->{exp_group} );

        if ( !defined $analysis ) {

            $to_register{Analysis} ||= [];
            push @{ $to_register{Analysis} }, $entry->{analysis};
            $abort = 1;

        }

        if ( !defined $feature_type ) {

            $to_register{Feature_Type} ||= [];
            push @{ $to_register{Feature_Type} }, $entry->{feature_type};
            $abort = 1;

        }

        if ( !defined $exp_group ) {

            $to_register{Experimental_Group} ||= [];
            push @{ $to_register{Experimental_Group} }, $entry->{exp_group};
            $abort = 1;

        }
    }

    if ($abort) {

        for my $object ( keys %to_register ) {
            for my $missing ( @{ $to_register{$object} } ) {
                $logger->warn( 'Register ' . $object . ': ' . $missing . "\n",
                    0, 1 );
            }
        }

        $logger->error( 'Aborting registration' . "\n", 0, 1 );
        exit 0;
    }

    return 1;
}

sub fetch_analysis {
    my ( $entry, $tracking_db_adaptors, $objects ) = @_;

    my $analysis = $tracking_db_adaptors->{an}
        ->fetch_by_logic_name( $entry->{analysis} );

    if ( !defined $analysis ) {
        throw "Register Analysis: '" . $entry->{analysis} . "'";
    }

    $objects->{analysis} = $analysis;

    return 1;
}

sub fetch_epigenome {
    my ( $entry, $tracking_db_adaptors, $objects ) = @_;

    my $epigenome_name = $entry->{epigenome_name};
    my $epigenome
        = $tracking_db_adaptors->{ep}->fetch_by_name($epigenome_name);

    if ( !defined $epigenome ) {
        $epigenome = Bio::EnsEMBL::Funcgen::Epigenome->new(
            -name          => $epigenome_name,
            -display_label => $epigenome_name,

            # -description   => '',

            # -gender        => $data->{sex},
            # -tissue => 'blood',
        );
        $tracking_db_adaptors->{ep}->store($epigenome);
    }

    $objects->{epigenome} = $epigenome;

    return 1;
}

sub fetch_feature_type {
    my ( $entry, $tracking_db_adaptors, $objects ) = @_;

    my $ft_name = $entry->{feature_type};

    my $feature_type = $tracking_db_adaptors->{ft}->fetch_by_name($ft_name);

    if ( !defined $feature_type ) {
        throw "Create new FeatureType for " . $ft_name;
    }

    $objects->{feature_type} = $feature_type;

    return 1;
}

sub fetch_exp_group {
    my ( $entry, $tracking_db_adaptors, $objects ) = @_;

    my $exp_group
        = $tracking_db_adaptors->{eg}->fetch_by_name( $entry->{exp_group} );

    if ( !defined $exp_group ) {
        throw "Create an entry for the experimental group: "
            . $entry->{exp_group};
    }

    $objects->{exp_group} = $exp_group;

    return 1;
}

sub fetch_experiment {
    my ( $entry, $tracking_db_adaptors, $objects, $cfg ) = @_;

    my $exp_name
        = $entry->{epigenome_name} . '_'
        . $entry->{feature_type} . '_'
        . $cfg->{general}->{study};

    my $experiment = $tracking_db_adaptors->{ex}->fetch_by_name($exp_name);

    if ( !defined $experiment ) {

        my $control_id;
        unless ( $entry->{is_control} ) {
            $control_id = $tracking_db_adaptors->{iss}
                ->fetch_by_name( $entry->{controlled_by} )->dbID();
        }

        $experiment = Bio::EnsEMBL::Funcgen::Experiment->new(
            -NAME               => $exp_name,
            -EPIGENOME          => $objects->{epigenome},
            -FEATURE_TYPE       => $objects->{feature_type},
            -EXPERIMENTAL_GROUP => $objects->{experimental_group},
            -IS_CONTROL         => $entry->{is_control},
            -CONTROL_ID         => $control_id,
        );

        $tracking_db_adaptors->{ex}->store($experiment);

  # $tracking_db_adaptors->{tr}->store_tracking_info( $experiment, $tr_info );

    }

    $objects->{experiment} = $experiment;

    return 1;
}

sub fetch_input_subset {
    my ( $entry, $tracking_db_adaptors, $objects ) = @_;

    my $iss
        = $tracking_db_adaptors->{iss}->fetch_by_name( $entry->{accession} );

    if ( !defined $iss ) {

        $iss = Bio::EnsEMBL::Funcgen::InputSubset->new(
            -name                 => $entry->{accession},
            -analysis             => $objects->{analysis},
            -epigenome            => $objects->{epigenome},
            -experiment           => $objects->{experiment},
            -feature_type         => $objects->{feature_type},
            -is_control           => $entry->{is_control},
            -biological_replicate => $entry->{br},
            -technical_replicate  => $entry->{tr},
        );
        $tracking_db_adaptors->{iss}->store($iss);

        my $tr_info->{info} = {

            # availability_date => '2015-10-13',
            # download_url      => $data->{download_url},
            # download_date     => $data->{download_date},
            local_url => $entry->{local_url},
            md5sum    => $entry->{md5sum},

            # notes             => $data->{experiment_name},
        };
        $tracking_db_adaptors->{tr}->store_tracking_info( $iss, $tr_info );

    }

    push( @{ $objects->{input_subsets} }, $iss );

    return 1;
}

sub fetch_ontology_xref {
    my ( $entry, $db_entry_adaptor, $objects ) = @_;

    my @ontology_accessions = split /;/, $entry->{ontology_accessions};
    my %valid_linkage_annotations
        = ( 'SAMPLE' => 1, 'TISSUE' => 1, 'CONDITION' => 1 );

    for my $ontology_accession (@ontology_accessions) {

        my ( $linkage_annotation, $primary_id ) = split /-/,
            $ontology_accession;
        my ($dbname) = split /:/, $primary_id;

        if ( !$valid_linkage_annotations{$linkage_annotation} ) {
            throw
                "Invalid linkage_annotation, please use 'SAMPLE', 'TISSUE' or 'CONDITION'";
        }

        my $ontology_xref = Bio::EnsEMBL::OntologyXref->new(
            -primary_id         => $primary_id,
            -dbname             => $dbname,
            -linkage_annotation => $linkage_annotation,

            # -display_id         => "",
            # -description        => "",
        );

        my $ignore_release = 1;

        my $epigenome_id = 1;

        $db_entry_adaptor->store( $ontology_xref, $epigenome_id, "epigenome",
            $ignore_release );
    }

}

sub usage {
    say "Usage: ";
    exit 0;
}
