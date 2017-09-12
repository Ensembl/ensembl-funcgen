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
use autodie;
use feature qw(say);

use Cwd 'abs_path';
use Config::Tiny;
use Bio::EnsEMBL::Utils::Logger;
use Bio::EnsEMBL::Utils::Exception qw(throw);
use Bio::EnsEMBL::OntologyXref;
use Bio::EnsEMBL::DBEntry;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor;
use Bio::EnsEMBL::DBSQL::DBEntryAdaptor;
use Getopt::Long;
use File::Basename;
use List::MoreUtils qw(uniq);

# use DateTime;

#TODO Dry run implementation
#TODO Rollback function
#TODO POD for every subroutine
#TODO Logger
#TODO PerlCritic
#TODO check for external db availability in verify_basic_objects
#TODO healthchecks, ie.invalid br/tr values, invalid local/download url
#TODO use state variable for $control_db_ids

main();

sub main {
    my ( $csv, $config, $help, $dry );

    # ----------------------------
    # read command line parameters
    # ----------------------------
    GetOptions(
        'i=s'  => \$csv,
        'c=s'  => \$config,
        'h'    => \$help,
        'help' => \$help,
        'n'    => \$dry,      #not implemented
    );

    # ------------------------------------------------------
    # display usage and exit if anything critical is missing
    # ------------------------------------------------------
    if ( $help || !$csv || !$config ) {
        usage();
        exit;
    }

    # -----------------
    # initialize logger
    # -----------------
    my $logger = Bio::EnsEMBL::Utils::Logger->new();
    $logger->init_log();

    # -------------------------------------
    # check that config and csv files exist
    # -------------------------------------
    if ( !-e $config ) {
        $logger->error(
            'Config file ' . abs_path($config) . ' doesn\'t exist!',
            0, 1 );
    }
    if ( !-e $csv ) {
        $logger->error(
            'Input csv file ' . abs_path($csv) . ' doesn\'t exist!',
            0, 1 );
    }

    $logger->info( 'CSV file used: ' . abs_path($csv) . "\n",       0, 1 );
    $logger->info( 'Config file used: ' . abs_path($config) . "\n", 0, 1 );

    # --------------------
    # read the config file
    # --------------------
    my $cfg = Config::Tiny->read($config);

    # ----------------
    # get current date
    # ----------------
    my ( $sec, $min, $hour, $mday, $mon, $year, $wday, $yday, $isdst )
        = localtime(time);
    $cfg->{date} = ( 1900 + $year ) . '-' . ++$mon . '-' . $mday;

    # -------------------
    # read input csv file
    # -------------------
    $logger->info( 'Reading metadata from ' . abs_path($csv) . '... ', 0, 0 );
    my ( $control_data, $signal_data )
        = get_data_from_csv( abs_path($csv), $logger );
    $logger->info( "done\n", 0, 1 );

    my $control_entries = keys %{$control_data};
    my $signal_entries  = keys %{$signal_data};
    $logger->info( $control_entries . " control entries found\n", 1, 1 );
    $logger->info( $signal_entries . " signal entries found\n",   1, 1 );

    # ------------------------------------------------------------
    # connect to funcgen tracking db, fetch all necessary adaptors
    # ------------------------------------------------------------
    $logger->info( 'Connecting to ' . $cfg->{efg_db}->{dbname} . '... ',
        0, 0 );
    my $adaptors = fetch_adaptors($cfg);
    $logger->info( "done\n", 0, 1 );

    # -----------------------------------------------
    # verify that fundamental objects exist in the db
    # -----------------------------------------------
    verify_basic_objects_in_db( $control_data, $signal_data, $adaptors,
        $logger );

    # ------------
    # registration
    # ------------
    my $control_db_ids = {};

    # control files have to be registered first, because their db ids are used
    # in the registration of signal files
    for my $entry ( values %{$control_data}, values %{$signal_data} ) {
        register( $logger, $entry, $adaptors, $cfg, $control_db_ids );
    }

    return 1;
}

sub get_data_from_csv {
    my ( $csv, $logger ) = @_;
    my ( %control_data, %signal_data );

    open my $csv_fh, '<', $csv;

    while ( readline $csv_fh ) {
        chomp;

        if (/^accession/i) {
            next;    # ignore input file header
        }

        my ($accession,             $epigenome_name,     $feature_type_name,
            $br,                    $tr,
            $paired_end_tag,        $multiple,
            $gender,
            $md5,                   $local_url,          $analysis_name,
            $exp_group_name,        $ontology_xref_accs, $xref_accs,
            $epigenome_description, $controlled_by,      $download_url,
            $info
        ) = split /\t/;

        my $entry = {};
        $entry->{accession}             = $accession;
        $entry->{epigenome_name}        = $epigenome_name;
        $entry->{feature_type_name}     = $feature_type_name;
        $entry->{br}                    = $br;
        $entry->{tr}                    = $tr;
        $entry->{paired_end_tag}        = $paired_end_tag;
        $entry->{multiple}              = $multiple;
        $entry->{md5}                   = $md5;
        $entry->{gender}                = $gender;
        $entry->{local_url}             = $local_url;
        $entry->{analysis_name}         = $analysis_name;
        $entry->{exp_group_name}        = $exp_group_name;
        $entry->{ontology_xref_accs}    = $ontology_xref_accs;
        $entry->{xref_accs}             = $xref_accs;
        $entry->{epigenome_description} = $epigenome_description;
        $entry->{controlled_by}         = $controlled_by;
        $entry->{download_url}          = $download_url;
        $entry->{info}                  = $info;

        $entry = verify_entry_metadata( $entry, $logger );

        if ( $control_data{$accession} || $signal_data{$accession} ) {
            $logger->error( 'Accession ' . $accession . ' is NOT unique!',
                0, 0 );
        }

        if ( $feature_type_name eq 'WCE' ) {
            $entry->{is_control} = 1;
            $control_data{$accession} = $entry;
        }
        else {
            $entry->{is_control} = 0;
            $signal_data{$accession} = $entry;
        }

    }

    close $csv_fh;

    return ( \%control_data, \%signal_data );
}

sub verify_entry_metadata {
    my ( $entry, $logger ) = @_;

    my @mandatory = qw(
        accession         epigenome_name
        feature_type_name br
        tr                md5
        local_url         analysis_name
        exp_group_name
    );

    for my $man (@mandatory) {
        if ( !$entry->{$man} ) {
            $logger->error(
                'There is no ' . $man . ' value for ' . $entry->{accession},
                0, 0 );
        }
    }

    my @optional = qw(
        gender        ontology_xref_accs
        xref_accs     epigenome_description
        controlled_by download_url
        info
    );

    for my $opt (@optional) {
        if ( $entry->{$opt} eq '-' ) {
            $entry->{$opt} = undef;
        }
    }

    return $entry;
}

# avoid using a register.conf file, use the existing mechanism
#
# use registry instead
# list and expose TrackingAdaptor in DBAdaptor
sub fetch_adaptors {
    my ($cfg) = @_;
    my %adaptors;

    my $dba = Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor->new(
        -user    => $cfg->{efg_db}->{user},
        -pass    => $cfg->{efg_db}->{pass},
        -host    => $cfg->{efg_db}->{host},
        -port    => $cfg->{efg_db}->{port},
        -dbname  => $cfg->{efg_db}->{dbname},
        -species => $cfg->{general}->{species},
    );

    $adaptors{epigenome}    = $dba->get_EpigenomeAdaptor();
    $adaptors{feature_type} = $dba->get_FeatureTypeAdaptor();
    $adaptors{analysis}     = $dba->get_AnalysisAdaptor();
    $adaptors{exp_group}    = $dba->get_ExperimentalGroupAdaptor();
    $adaptors{experiment}   = $dba->get_ExperimentAdaptor();
    $adaptors{read_file_experimental_configuration}    = $dba->get_ReadFileExperimentalConfigurationAdaptor();
    $adaptors{read_file}    = $dba->get_ReadFileAdaptor();
    

    $adaptors{db}       = $dba;
    $adaptors{db_entry} = Bio::EnsEMBL::DBSQL::DBEntryAdaptor->new($dba);

    return \%adaptors;
}

sub verify_basic_objects_in_db {
    my ( $control_data, $signal_data, $adaptors, $logger ) = @_;
    my $abort = 0;
    my %to_register;

    for my $entry ( values %{$control_data}, values %{$signal_data} ) {

        my $analysis = $adaptors->{analysis}
            ->fetch_by_logic_name( $entry->{analysis_name} );

        my $feature_types = $adaptors->{feature_type}
            ->fetch_all_by_name( $entry->{feature_type_name} );

        my $exp_group = $adaptors->{exp_group}
            ->fetch_by_name( $entry->{exp_group_name} );

        if ( !defined $analysis ) {
            $to_register{Analysis} //= [];
            push @{ $to_register{Analysis} }, $entry->{analysis_name};
            $abort = 1;
        }

        if ( !defined $feature_types->[0] ) {
            $to_register{Feature_Type} //= [];
            push @{ $to_register{Feature_Type} }, $entry->{feature_type_name};
            $abort = 1;
        }

        if ( !defined $exp_group ) {
            $to_register{Experimental_Group} //= [];
            push @{ $to_register{Experimental_Group} },
                $entry->{exp_group_name};
            $abort = 1;
        }
    }

    if ($abort) {
        for my $object ( keys %to_register ) {
            for my $missing ( uniq @{ $to_register{$object} } ) {
                $logger->warning(
                    'Register ' . $object . ': ' . $missing . "\n",
                    0, 1 );
            }
        }

        $logger->error( 'Aborting registration' . "\n", 0, 1 );
    }

    return 1;
}

sub register {
    my ( $logger, $entry, $adaptors, $cfg, $control_db_ids ) = @_;

    $logger->info( 'Registering ' . $entry->{accession} . "\n", 0, 1 );

    my $analysis = $adaptors->{analysis}
        ->fetch_by_logic_name( $entry->{analysis_name} );

    my $epigenome
        = $adaptors->{epigenome}->fetch_by_name( $entry->{epigenome_name} );

    if ( !$epigenome ) {
        $epigenome = store_epigenome( $entry, $adaptors );
    }

    if ( $entry->{ontology_xref_accs} ) {
        store_ontology_xref( $entry, $adaptors, $epigenome );
    }

    if ( $entry->{xref_accs} ) {
        store_db_xref( $entry, $adaptors, $epigenome );
    }

    my $feature_type = fetch_feature_type( $entry, $adaptors, $analysis );
    
    if (! defined $feature_type) {
      die("Couldn't fetch feature type with name " . $entry->{feature_type_name} . "!");
    }
    
    my $exp_group
        = $adaptors->{exp_group}->fetch_by_name( $entry->{exp_group_name} );

    my $experiment_name;
    if ( $entry->{is_control} ) {
        $experiment_name
            = create_control_experiment_name( $entry, $cfg, $epigenome,
            $adaptors );
    }
    else {
        $experiment_name
            = create_signal_experiment_name( $entry, $cfg, $epigenome,
            $adaptors );
    }

    my $experiment = $adaptors->{experiment}->fetch_by_name($experiment_name);

    # don't use control_db_ids, fetch directly from db for every signal file
    if ( !$experiment ) {
        $experiment = store_experiment(
            $entry,     $experiment_name, $adaptors, $control_db_ids,
            $epigenome, $feature_type,    $exp_group
        );
    }

    # avoid passing too many parameters to subroutines
    store_read_file( $logger, $entry, $adaptors, $cfg, $analysis,
        $epigenome, $experiment, $feature_type );

    $logger->info( "Successful Registration\n", 1, 1 );

    return 1;
}

sub store_epigenome {
    my ( $entry, $adaptors ) = @_;

    my $production_name
        = create_epigenome_production_name( $entry->{epigenome_name} );

    my $epigenome = Bio::EnsEMBL::Funcgen::Epigenome->new(
        -name            => $entry->{epigenome_name},
        -display_label   => $entry->{epigenome_name},
        -description     => $entry->{epi_description},
        -production_name => $production_name,
        -gender          => $entry->{gender}
    );

    $adaptors->{epigenome}->store($epigenome);

    return $epigenome;
}

sub fetch_feature_type {
    my ( $entry, $adaptors, $analysis ) = @_;

    my $ft_name = $entry->{feature_type_name};

    my $feature_type;

    if ( $entry->{is_control} ) {
        $feature_type = $adaptors->{feature_type}
            ->fetch_by_name( $ft_name, 'DNA', $analysis );
    }
    else {
        $feature_type = $adaptors->{feature_type}->fetch_by_name($ft_name);
    }

    return $feature_type;
}

sub create_control_experiment_name {
    my ( $entry, $cfg, $epigenome, $adaptors ) = @_;

    my $experiment_name;
    my $number = 1;

    do {
        $experiment_name
            = $epigenome->{production_name} . '_'
            . $entry->{feature_type_name} . '_'
            . $entry->{analysis_name} . '_no'
            . $number . '_'
            . $entry->{exp_group_name};

        if ( $cfg->{general}->{release} ) {
            $experiment_name .= $cfg->{general}->{release};
        }

        $experiment_name =~ s/\s//g;
        $number++;
    } while ( $adaptors->{experiment}->fetch_by_name($experiment_name) );

    return $experiment_name;
}

sub create_signal_experiment_name {
    my ( $entry, $cfg, $epigenome, $adaptors ) = @_;

    my $experiment_name;

    $experiment_name
        = $epigenome->{production_name} . '_'
        . $entry->{feature_type_name} . '_'
        . $entry->{analysis_name} . '_'
        . $entry->{exp_group_name};

    if ( $cfg->{general}->{release} ) {
        $experiment_name .= $cfg->{general}->{release};
    }

    $experiment_name =~ s/\s//g;

    return $experiment_name;
}

sub store_experiment {
    my ($entry,     $experiment_name, $adaptors, $control_db_ids,
        $epigenome, $feature_type,    $exp_group
    ) = @_;

    my $control_experiment;
    if ( !$entry->{is_control} && $entry->{controlled_by} ) {

        my $control_db_id = $control_db_ids->{ $entry->{controlled_by} };

        $control_experiment
            = $adaptors->{experiment}->fetch_by_dbID($control_db_id);
    }
    
    my $experiment = Bio::EnsEMBL::Funcgen::Experiment->new(
        -NAME               => $experiment_name,
        -EPIGENOME          => $epigenome,
        -FEATURE_TYPE       => $feature_type,
        -EXPERIMENTAL_GROUP => $exp_group,
        -IS_CONTROL         => $entry->{is_control},
        -CONTROL            => $control_experiment,
    );

    $adaptors->{experiment}->store($experiment);

    if ( $entry->{is_control} ) {
        $control_db_ids->{ $entry->{accession} } = $experiment->dbID();
    }

    return $experiment;
}

sub store_read_file {
    my ( $logger, $entry, $adaptors, $cfg, $analysis, $epigenome, $experiment,
        $feature_type )
        = @_;

    my $read_file = $adaptors->{read_file}
        ->fetch_by_name( $entry->{accession}, $experiment );

    if ($read_file) {
        $logger->warning(
            'A read file entry for accession '
                . $entry->{accession}
                . ' with experiment name '
                . $experiment->{name}
                . ' already exists in DB! ',
            0, 1
        );
    }
    
    my $paired_end_tag = $entry->{paired_end_tag};
    if ($paired_end_tag eq '-') {
      $paired_end_tag = undef;
    }
    
    use Bio::EnsEMBL::Funcgen::ReadFile;
    use Bio::EnsEMBL::Funcgen::ReadFileExperimentalConfiguration;
    
    my $read_file_experimental_configuration = Bio::EnsEMBL::Funcgen::ReadFileExperimentalConfiguration->new(
    
      -read_file => Bio::EnsEMBL::Funcgen::ReadFile->new(
          -name           => $entry->{accession},
          -analysis       => $analysis,
          -is_paired_end  => undef,
          -paired_with    => undef,
          -file_size      => undef,
          -read_length    => undef,
          -md5sum         => $entry->{md5},
          -file           => $entry->{local_url},
          -notes          => $entry->{info},
      ),
      -experiment            => $experiment,
      -biological_replicate  => $entry->{br},
      -technical_replicate   => $entry->{tr},
      -paired_end_tag        => $paired_end_tag,
      -multiple              => $entry->{multiple},
    );

    $adaptors->{read_file_experimental_configuration}->store($read_file_experimental_configuration);

    return 1;
}

sub store_ontology_xref {
    my ( $entry, $adaptors, $epigenome_obj ) = @_;

    my @ontology_accessions = split /;/, $entry->{ontology_xref_accs};
    my %valid_linkage_annotations
        = ( 'SAMPLE' => 1, 'TISSUE' => 1, 'CONDITION' => 1 );

    for my $ontology_accession (@ontology_accessions) {

        my ( $linkage_annotation, $primary_id ) = split /-/,
            $ontology_accession;
        my ($dbname) = split /:/, $primary_id;

        if ( !$valid_linkage_annotations{$linkage_annotation} ) {
            throw
                'Invalid linkage_annotation, please use \'SAMPLE\', \'TISSUE\' or \'CONDITION\'';
        }

        my $ontology_xref = Bio::EnsEMBL::OntologyXref->new(
            -primary_id         => $primary_id,
            -dbname             => $dbname,
            -linkage_annotation => $linkage_annotation,
        );

        my $ignore_release = 1;

        my $epigenome_id = $epigenome_obj->dbID();

        $adaptors->{db_entry}
            ->store( $ontology_xref, $epigenome_id, 'epigenome',
            $ignore_release );
    }

    return 1;
}

sub store_db_xref {
    my ( $entry, $adaptors, $epigenome_obj ) = @_;

    my @xref_accessions = split /;/, $entry->{xref_accs};

    for my $xref_acc (@xref_accessions) {
        my ( $primary_id, $dbname ) = split /-/, $xref_acc;

        my $xref = Bio::EnsEMBL::DBEntry->new(
            -primary_id => $primary_id,
            -dbname     => $dbname,
            -display_id => undef, # if not set to undef, an empty string will be stored in the display_label column
            -info_text  => undef # if not set to undef, an empty string will be stored in the info_text column
        );

        my $ignore_release = 1;

        my $epigenome_id = $epigenome_obj->dbID();

        $adaptors->{db_entry}
            ->store( $xref, $epigenome_id, 'epigenome', $ignore_release );
    }

    return 1;
}

sub create_epigenome_production_name {
    my ($epigenome_name) = shift;

    $epigenome_name =~ s/\:/_/g;
    $epigenome_name =~ s/\+//g;
    $epigenome_name =~ s/\(//g;
    $epigenome_name =~ s/\)//g;
    $epigenome_name =~ s/\-/_/g;
    $epigenome_name =~ s/\./_/g;
    $epigenome_name =~ s/\//_/g;
    $epigenome_name =~ s/ /_/g;

    return $epigenome_name;
}

sub usage {
    my $usage = << 'END_USAGE';

Usage: register_metadata.pl -i <input_file> -c <config_file>

Options:
-i input_file:  this is the tab delimited text file that contains the metadata
-c config_file: this is the configuration file that contains the database connection details
-n dry_run:     run the script without storing any data into the database, NOT YET IMPLEMENTED
-h help:        shows this help message
 
END_USAGE

    say $usage;

    return 1;
}
