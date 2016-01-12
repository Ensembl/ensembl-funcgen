#!/usr/bin/perl

use strict;
use warnings;
use diagnostics;
use autodie;
use feature qw(say);
use Data::Dumper qw( Dumper );
use Config::Tiny;
use DateTime;

use Bio::EnsEMBL::Funcgen::Utils::EFGUtils qw( dump_data get_date);
use Bio::EnsEMBL::Utils::Exception qw(throw);
use Bio::EnsEMBL::Funcgen::DBSQL::TrackingAdaptor;
use Bio::EnsEMBL::Funcgen::CellType;
use Bio::EnsEMBL::Funcgen::FeatureType;
use Bio::EnsEMBL::Funcgen::Experiment;
use Bio::EnsEMBL::Funcgen::ExperimentalGroup;
use Bio::EnsEMBL::Funcgen::InputSubset;

main();

sub main {
    my $config = $ENV{'CTTV020_DIR'} . '/register.conf' or die;

    my $cfg = Config::Tiny->new;
    $cfg = Config::Tiny->read($config);

    _connect_to_trackingDB($cfg);
    _get_trackingDB_adaptors($cfg);

    my $data;
    open my $in, '<', $ARGV[0];

    my $key;

    while (<$in>) {
        chomp;

        my @columns  = split /\t/;
        $key++;

        # my $filename_no_extension = substr $filename, 0, 9;

        (   $data->{$key}->{filename},
            $data->{$key}->{cell_type_name},
            $data->{$key}->{feature_type},

            $data->{$key}->{replicate},

            $data->{$key}->{md5sum}
        ) = @columns;

        $data->{$key}->{analysis} = $cfg->{generic}->{analysis};
        $data->{$key}->{study_id} = $cfg->{generic}->{study_id};

        $data->{$key}->{download_url}
            = 'Multiple files merged. Check notes on table input_subset_tracking';

        # = 'iRODS /seq/'
        # . $data->{$key}->{run_id} . '/'
        # . $data->{$key}->{run_id} . '_'
        # . $data->{$key}->{lane} . '#'
        # . $data->{$key}->{tag_index} . '.cram';

        $data->{$key}->{local_url}
            = $ENV{'WAREHOUSE_DIR'} . '/' . $data->{$key}->{filename};

        $data->{$key}->{download_date}
            = DateTime->now( time_zone => "local" )->ymd('-');
    }

    close $in;

    my $objects = {};

    foreach my $key ( keys %{$data} ) {
        fetch_analysis( $cfg, $data->{$key}, $objects );

        # say $data->{$key}->{cell_type_name};
        create_experiment_name( $data->{$key} );

        # say $data->{$key}->{experiment_name};
        my $exp_name       = $data->{$key}->{experiment_name};
        my $exp_group_name = $data->{$key}->{experimental_group_name}
            = 'CTTV020';

        fetch_cell_type( $cfg, $data->{$key}, $objects );

        fetch_feature_type( $cfg, $data->{$key}, $objects );
        fetch_exp_group( $cfg, $data->{$key}, $objects );
        fetch_experiment( $cfg, $data->{$key}, $objects, $exp_name );

        fetch_input_subset( $cfg, $data->{$key}, $objects, $data->{$key}->{filename} );

        # foreach my $key ( sort keys %{$objects} ) {
        #     say "$key: " . ref( $objects->{$key} );
        # }
    }

}

###############################################################################
#                          create_and_store_experimental_group
###############################################################################

sub fetch_exp_group {
    my ( $cfg, $data, $objects ) = @_;

    my $eg = $cfg->{tr_a}->{eg}
        ->fetch_by_name( $data->{experimental_group_name} );
    if ( !defined $eg ) {

# throw "Create an entry for the experimental group: " . $data->{experimental_group_name};
        $eg = Bio::EnsEMBL::Funcgen::ExperimentalGroup->new(
            -NAME       => $data->{experimental_group_name},
            -URL        => 'http://www.targetvalidation.org/',
            -is_project => 1,
            -description =>
                'Centre for Therapeutic Target Validation, Epigenome profiles of GSK cells'
        );

        $cfg->{tr_a}->{eg}->store($eg);

        # my $tr_info->{info} = {

      #     # experiment_id     => $experiment->dbID,
      #     notes => 'Experimental_Group:' . $data->{experimental_group_name},
      # };
      # $cfg->{tr_a}->{eg}->store_tracking_info( $eg, $tr_info );

    }
    $objects->{experimental_group} = $eg;
}

###############################################################################
#                          create_and_store_input_subset
###############################################################################
sub fetch_input_subset {
    my ( $cfg, $data, $objects, $file ) = @_;

    # my $iss = $cfg->{tr_a}->{iss}->fetch_by_name($file, $data->{experiment_name});
    my $iss;
    if ( !defined $iss ) {
        my $is_control = ( $data->{feature_type} eq 'WCE' ) ? 1 : 0;

        $iss = Bio::EnsEMBL::Funcgen::InputSubset->new(
            -analysis     => $objects->{analysis},
            -cell_type    => $objects->{cell_type},
            -experiment   => $objects->{experiment},
            -feature_type => $objects->{feature_type},
            -is_control   => $is_control,
            -name         => $file,
            -replicate    => $data->{replicate},
        );
        $cfg->{tr_a}->{iss}->store($iss);

        my $tr_info->{info} = {
            availability_date => '2015-10-13',
            download_url      => $data->{download_url},
            download_date     => $data->{download_date},
            local_url         => $data->{local_url},
            md5sum            => $data->{md5sum},
            notes             => $data->{experiment_name},
        };

        if ( $file eq 'ctcf_input.fastq.gz' ) {
            $tr_info->{info}->{notes}
                .= ' This file is a merge of: '
                . ' 16723_2#28.cram, 16723_3#28.cram, 16723_4#28.cram, 16723_5#28.cram, 16723_6#28.cram, 16723_2#25.cram, 16723_3#25.cram, 16723_4#25.cram, 16723_5#25.cram, 16723_6#25.cram, 15688_6#8.cram, 15688_7#8.cram, 15688_8#8.cram, 15688_6#10.cram, 15688_7#10.cram, 15688_8#10.cram, 16723_2#31.cram, 16723_3#31.cram, 16723_4#31.cram, 16723_5#31.cram, 16723_6#31.cram'
                . ' They can be retrieved from iRODS /seq/16723 & /seq/15688' ;
        }

        if ( $file eq 'hist_input.fastq.gz' ) {
            $tr_info->{info}->{notes}
                .= ' This file is a merge of: '
                . '16493_8#35.cram, 16557_5#35.cram, 16557_6#35.cram, 16557_7#35.cram, 16557_8#35.cram, 16493_8#24.cram, 16557_5#24.cram, 16557_6#24.cram, 16557_7#24.cram, 16557_8#24.cram, 15688_6#3.cram, 15688_7#3.cram, 15688_8#3.cram, 15688_6#6.cram, 15688_7#6.cram, 15688_8#6.cram, 16493_8#36.cram, 16557_5#36.cram, 16557_6#36.cram, 16557_7#36.cram, 16557_8#36.cram'
                . ' They can be retrieved from iRODS /seq/16493 & /seq/16557 & /seq/15688';
        }

        $cfg->{tr_a}->{tr}->store_tracking_info( $iss, $tr_info );

    }

    push( @{ $objects->{input_subsets} }, $iss );

    return;
}

###############################################################################
#                          create_and_store_experiment
###############################################################################
# Could cache experiment, but again, readability vs efficiency

sub fetch_experiment {
    my ( $cfg, $data, $objects, $exp_name ) = @_;

    my $experiment = $cfg->{tr_a}->{ex}->fetch_by_name($exp_name);
    if ( !defined $experiment ) {

        my $exp_name = $data->{experiment_name};
        $experiment = Bio::EnsEMBL::Funcgen::Experiment->new(
            -NAME               => $data->{experiment_name},
            -CELL_TYPE          => $objects->{cell_type},
            -FEATURE_TYPE       => $objects->{feature_type},
            -ARCHIVE_ID         => $data->{experiment_id},
            -EXPERIMENTAL_GROUP => $objects->{experimental_group},
            -DESCRIPTION        => $data->{description},
        );

        $cfg->{tr_a}->{ex}->store($experiment);

        my $tr_info->{info} = {

            # experiment_id     => $experiment->dbID,
            notes => 'Experiment:' . $data->{experiment_name},
        };
        $cfg->{tr_a}->{tr}->store_tracking_info( $experiment, $tr_info );

    }
    $objects->{experiment} = $experiment;
}

###############################################################################
#                           fetch_feature_type
###############################################################################
sub fetch_feature_type {
    my ( $cfg, $data, $objects, $json ) = @_;

    my $ft_name = $data->{feature_type};
 
    my $feature_type = $cfg->{tr_a}->{ft}->fetch_by_name($ft_name);

    if ( !defined $feature_type ) {
        throw "Create new FeatureType for " . $ft_name;

        # if ( $json->{target}->{investigated_as}->[0] eq
        #     'transcription factor' )
        # {

        # $feature_type = Bio::EnsEMBL::Funcgen::FeatureType->new(
        #     -name         => $ft_name,
        #     -class        => 'Transcription Factor',
        #     -description  => $ft_name . ' Transcription Factor Binding',
        #     -so_accession => 'SO:0000235',
        #     -so_name      => 'TF_binding_site',
        # );
        # $cfg->{tr_a}->{ft}->store($feature_type);
        # }

        # else {
        #     my $id  = $data->{experiment_id};
        #     my $url = $data->{experiment_url};
        #     throw "Register FeatureType: '$ft_name' URL: $id $url";
        # }
    }

    $objects->{feature_type} = $feature_type;

    return;
}

###############################################################################
#                           fetch_cell_feature_type
###############################################################################
sub fetch_cell_type {
    my ( $cfg, $data, $objects ) = @_;

    my $ct_name   = $data->{cell_type_name};
    my $cell_type = $cfg->{tr_a}->{ct}->fetch_by_name($ct_name);

    if ( !defined $cell_type ) {
        $cell_type = Bio::EnsEMBL::Funcgen::CellType->new(
            -name          => $data->{cell_type_name},
            -display_label => $data->{cell_type_name},
            -description   => '',
            -gender        => $data->{sex},
            -tissue        => 'blood',
        );
        $cfg->{tr_a}->{ct}->store($cell_type);
    }

    $objects->{cell_type} = $cell_type;

    return;
}

###############################################################################
#                          create_experiment_name
###############################################################################

sub create_experiment_name {
    my ($data) = @_;

    # say dump_data($data,1,1);die;

    $data->{experiment_name}
        = $data->{cell_type_name} . '_'
        . $data->{feature_type} . '_'
        . $data->{study_id};

    return;

    # my $name;
    # $name .= '_'         . $data->{feature_type};
    # $name .= '_ENCODE_'  . $data->{lab};
    # $name .= '_BR'       . $file->{bio_replicate};
    # $data->{experiment_name} = $data->{cell_type_name} . $name;
    # return $name;

}

###############################################################################
#                          create_cell_type_name
###############################################################################

sub create_cell_type_name {
    my ( $cfg, $data, $file ) = @_;

    my $name;
    $name .= $data->{cell_type};

    if ( $data->{br} ) {
        $name .= ':' . $data->{br};
    }

    $data->{cell_type_name} = $name;
    return;

}

###############################################################################
#                              fetch_analysis
###############################################################################
sub fetch_analysis {
    my ( $cfg, $data, $objects ) = @_;

    my $anal = $cfg->{tr_a}->{an}->fetch_by_logic_name( $data->{analysis} );

    if ( !defined $anal ) {
        throw "Register Analysis: '" . $data->{analysis} . "'";
    }

    $objects->{analysis} = $anal;

    return;
}

###############################################################################
#                            _connect_to_trackingDB
###############################################################################
sub _connect_to_trackingDB {
    my ($cfg) = @_;

    # say dump_data($cfg->{efg_db},1,1);
    my $db_a = Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor->new(
        -user       => $cfg->{efg_db}->{user},
        -pass       => $cfg->{efg_db}->{pass},
        -host       => $cfg->{efg_db}->{host},
        -port       => $cfg->{efg_db}->{port},
        -dbname     => $cfg->{efg_db}->{dbname},
        -dnadb_name => $cfg->{dna_db}->{dbname},
    );
    $db_a->dbc->do("SET sql_mode='traditional'");
    say "\nConnected to trDB: " . $cfg->{efg_db}->{dbname} . "\n";

    return ( $cfg->{dba_tracking} = $db_a );
}

###############################################################################
#                            _get_trackingDB_adaptors
###############################################################################
sub _get_trackingDB_adaptors {
    my ($cfg) = @_;

    # Tracking DB hidden from user, hence no get_TrackingAdaptor method.
    # TrackingAdaptor->new() does not YET accept DBAdaptor object

    $cfg->{tr_a}->{tr} = Bio::EnsEMBL::Funcgen::DBSQL::TrackingAdaptor->new(
        -user       => $cfg->{efg_db}->{user},
        -pass       => $cfg->{efg_db}->{pass},
        -host       => $cfg->{efg_db}->{host},
        -port       => $cfg->{efg_db}->{port},
        -dbname     => $cfg->{efg_db}->{dbname},
        -species    => $cfg->{generic}->{species},
        -dnadb_user => $cfg->{dna_db}->{user},
        -dnadb_pass => $cfg->{dna_db}->{pass},
        -dnadb_host => $cfg->{dna_db}->{host},
        -dnadb_port => $cfg->{dna_db}->{port},
        -dnadb_name => $cfg->{dna_db}->{dbname},
    );

    my $db_a = $cfg->{tr_a}->{tr}->db;

    $cfg->{tr_a}->{ct} = $db_a->get_CellTypeAdaptor();
    $cfg->{tr_a}->{ft} = $db_a->get_FeatureTypeAdaptor();
    $cfg->{tr_a}->{an} = $db_a->get_AnalysisAdaptor();

    $cfg->{tr_a}->{eg} = $db_a->get_ExperimentalGroupAdaptor();

    $cfg->{tr_a}->{ex}  = $db_a->get_ExperimentAdaptor();
    $cfg->{tr_a}->{iss} = $db_a->get_InputSubsetAdaptor();

    $cfg->{tr_a}->{rs} = $db_a->get_ResultSetAdaptor();
    $cfg->{tr_a}->{rf} = $db_a->get_RegulatoryFeatureAdaptor();

    $cfg->{tr_a}->{fs} = $db_a->get_FeatureSetAdaptor();
    $cfg->{tr_a}->{ds} = $db_a->get_DataSetAdaptor();
    $cfg->{tr_a}->{af} = $db_a->get_AnnotatedFeatureAdaptor();
}
