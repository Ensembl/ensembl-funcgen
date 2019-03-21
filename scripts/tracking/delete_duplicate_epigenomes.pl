#!/usr/bin/env perl

# Copyright [1999-2018] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
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
use Getopt::Long qw(GetOptions);
use Config::Tiny;
use Bio::EnsEMBL::Utils::Logger;
use Bio::EnsEMBL::Utils::Exception qw(throw);
use DBI;
use Cwd 'abs_path';
use Try::Tiny;



main();

sub main {

    my ($cfgFile, $help);

    GetOptions(
                'c=s' => \$cfgFile,
                'h=s' => \$help,);

    # -----------------
    # check parameters
    # -----------------           
    if ( $help || !$cfgFile) {
        usage();
        exit;
    }

    # -----------------
    # initialize logger
    # -----------------
    my $logger = Bio::EnsEMBL::Utils::Logger->new();
    $logger->init_log();

    # ---------------------------
    # check if config file exists
    # ---------------------------
    if (not(-f $cfgFile)){
        die "$cfgFile file can not be found\n";
    }

    # ---------------------------
    # load config file
    # ---------------------------
    my $cfg = Config::Tiny->read($cfgFile);

    $logger->info( 'Config file used: ' . abs_path($cfg) . "\n", 0, 1 );

    # ---------------------------
    # connect to DB
    # ---------------------------
    my $connection = DBI->connect("DBI:mysql:".$cfg->{efg_db}->{dbname}.":".$cfg->{efg_db}->{host}.":".$cfg->{efg_db}->{port}, $cfg->{efg_db}->{user}, $cfg->{efg_db}->{pass}, {RaiseError => 1});
    
    $logger->info( 'Connecting to DB: ' . $cfg->{efg_db}->{dbname} . "\n", 0, 1 );

    # ------------------------------------
    # get duplicate epigenomes using ihec
    # ------------------------------------
    $logger->info( 'Getting duplicate ihecs...'. "\n", 0, 1 );
    my @duplicated_ihecs = get_duplicate_ihecs($connection);
    my $num_duplicates = scalar(@duplicated_ihecs);
    $logger->info( 'Found '. $num_duplicates . ' duplicated epigenomes using ihec'. "\n", 0, 1 );

    # ------------------------------------------
    # remove old duplicate epigenomes using ihec
    # ------------------------------------------
    $logger->info( 'Removing old duplicate epigenomes' . "\n", 0, 1 );

    foreach my $ihec (@duplicated_ihecs) {

        #get old epigenome
        my %old_epigenome = get_older_duplicate_epigenome($ihec, 'EpiRR', $connection);

        if (!%old_epigenome){
            next;
        }

        #get epigenome_id from the new epigenome 
        my $new_epigenome = get_newer_duplicate_epigenome($ihec, 'EpiRR', $connection);

        #remove old epigenome
        remove_epigenome($old_epigenome{'epigenome_id'}, $connection);

        #transfer old names to the new epigenome
        transfer_old_names(\%old_epigenome, $new_epigenome, $connection);

    }

    # -----------------------------------------------
    # get duplicate epigenomes using ENCODE accession
    # --------------------------------------------------
    $logger->info( 'Getting duplicate by Encode accession...'. "\n", 0, 1 );
    my @duplicated_encode = get_duplicate_encode($connection);
    my $num_duplicates_encode = scalar(@duplicated_encode);
    $logger->info( 'Found '. $num_duplicates_encode . ' duplicated epigenomes using ENCODE accession'. "\n", 0, 1 );

    # ----------------------------------------------------
    # remove old duplicate epigenomes by ENCODE accession
    # -----------------------------------------------------
    $logger->info( 'Removing old duplicate epigenomes' . "\n", 0, 1 );

    foreach my $encode_accession (@duplicated_encode) {

        #get old epigenome
        my %old_epigenome = get_older_duplicate_epigenome($encode_accession, 'ENCODE', $connection);

        if (!%old_epigenome){
            next;
        }

        #get epigenome_id from the new epigenome 
        my $new_epigenome = get_newer_duplicate_epigenome($encode_accession, 'ENCODE', $connection);

        #remove old epigenome
        remove_epigenome($old_epigenome{'epigenome_id'}, $connection);

        #transfer old names to the new epigenome
        transfer_old_names(\%old_epigenome, $new_epigenome, $connection);

    }




    # ------------------------------------------
    # remove in cascade
    # ------------------------------------------
    $logger->info( 'Removing in cascade' . "\n", 0, 1 );

    #remove experiment
    remove_experiment($connection);

    #remove experimental group
    remove_experimental_group($connection);

    #remove read_file_experimental_configuration
    remove_read_file_experimental_configuration($connection);

    #remove read_file
    remove_read_file($connection);

    #remove object_xref Epigenome
    remove_object_xref_epigenome($connection);

    #remove object_xref ReadFile
    remove_object_xref_read_file($connection);

    #remove xref
    remove_xref($connection);

    $connection->disconnect();
    $logger->info( 'Process ended' . "\n", 0, 1 );
    return 1;

}

sub get_duplicate_encode{
    my $connection = shift;
    
    my $sql = "select dbprimary_acc
                from epigenome
                left join object_xref on (epigenome.epigenome_id = object_xref.ensembl_id and ensembl_object_type = 'epigenome')
                left join xref using (xref_id)
                join external_db using (external_db_id)
                where db_name = 'ENCODE'
                group by dbprimary_acc
                having count(1)> 1";
    my $sth = $connection->prepare($sql);
    $sth->execute();

    my @duplicated_encode;
    my @row;
    while (@row = $sth->fetchrow_array) {
        push (@duplicated_encode, $row[0]);
    }

    return @duplicated_encode;

}

sub transfer_old_names{
    my $old_epigenome = shift;
    my $new_epigenome = shift;
    my $connection = shift;
    my $sql = 'UPDATE epigenome
            SET epigenome.name = ? , epigenome.short_name =?, epigenome.description = ?, epigenome.production_name = ?
            WHERE epigenome.epigenome_id = ? ';

    my $sth = $connection->prepare($sql);
    $sth->bind_param( 1, $old_epigenome->{'name'} );
    $sth->bind_param( 2, $old_epigenome->{'short_name'} );
    $sth->bind_param( 3, $old_epigenome->{'description'} );
    $sth->bind_param( 4, $old_epigenome->{'production_name'} );
    $sth->bind_param( 5, $new_epigenome );
    $sth->execute();

    return;
}

sub remove_xref{
    my $connection = shift;

    my $sql = "Delete xref from xref
            LEFT JOIN object_xref ON object_xref.xref_id = xref.xref_id
            where object_xref.xref_id is null ";

    my $sth = $connection->prepare($sql);
    $sth->execute();

    return;

}

sub remove_object_xref_epigenome{

    my $connection = shift;

    my $sql = "Delete object_xref from object_xref
            LEFT JOIN epigenome ON epigenome.epigenome_id = object_xref.ensembl_id
            where ensembl_object_type = 'Epigenome' and epigenome.epigenome_id is null ";

    my $sth = $connection->prepare($sql);
    $sth->execute();

    return;

}

sub remove_object_xref_read_file{

    my $connection = shift;

    my $sql = "Delete object_xref from object_xref
            LEFT JOIN read_file ON read_file.read_file_id = object_xref.ensembl_id
            where ensembl_object_type = 'ReadFile' and read_file.read_file_id is null ";

    my $sth = $connection->prepare($sql);
    $sth->execute();

    return;

}



sub remove_read_file{
    my $connection = shift;

    my $sql = "Delete read_file from read_file
            LEFT JOIN read_file_experimental_configuration ON read_file_experimental_configuration.read_file_id = read_file.read_file_id
            where read_file_experimental_configuration.read_file_id is null ";

    my $sth = $connection->prepare($sql);
    $sth->execute();

    return;

}

sub remove_read_file_experimental_configuration{
    my $connection = shift;

    my $sql = "Delete read_file_experimental_configuration from read_file_experimental_configuration
            LEFT JOIN experiment ON experiment.experiment_id = read_file_experimental_configuration.experiment_id
            where experiment.experiment_id is null ";

    my $sth = $connection->prepare($sql);
    $sth->execute();

    return;

}

sub remove_experimental_group{
    my $connection = shift;

    my $sql = "Delete experimental_group from experimental_group
            LEFT JOIN experiment ON experiment.experimental_group_id=experimental_group.experimental_group_id
            where experiment.experimental_group_id is null ";

    my $sth = $connection->prepare($sql);
    $sth->execute();

    return;

}

sub remove_experiment{
    my $connection = shift;

    my $sql = "Delete experiment from experiment
            left join epigenome on epigenome.epigenome_id = experiment.epigenome_id
            where epigenome.epigenome_id is null ";

    my $sth = $connection->prepare($sql);
    $sth->execute();

    return;

}


sub remove_epigenome{
    my $epigenome = shift;
    my $connection = shift;

    my $sql = "delete from epigenome where epigenome_id = ? ";

    my $sth = $connection->prepare($sql);
    $sth->bind_param( 1, $epigenome );
    $sth->execute();

    return;

}

sub get_older_duplicate_epigenome{
    my $dbprimary_acc = shift;
    my $db_name = shift;
    my $connection = shift;

    my $sql = "select epigenome_id, epigenome.name, epigenome.short_name, epigenome.description, epigenome.production_name
            from epigenome
            left join object_xref on (epigenome.epigenome_id = object_xref.ensembl_id and ensembl_object_type = 'epigenome')
            left join xref using (xref_id)
            join external_db using (external_db_id)
            where dbprimary_acc = ? and external_db.db_name = ?
            order by epigenome_id ASC";

    my $sth = $connection->prepare($sql);
    $sth->bind_param( 1, $dbprimary_acc );
    $sth->bind_param( 2, $db_name );
    $sth->execute();
    my @row = $sth->fetchrow_array;
    my %epi_values;
    #if (!($row[1] =~ m/^new_/g)) {
        print "Found!!! ".$row[1]."\n";
        $epi_values{'epigenome_id'}=$row[0];
        $epi_values{'name'}=$row[1];
        $epi_values{'short_name'}=$row[2];
        $epi_values{'description'}=$row[3];
        $epi_values{'production_name'}=$row[4];
    #}
    return %epi_values;
}

sub get_newer_duplicate_epigenome{
    my $dbprimary_acc = shift;
    my $db_name = shift;
    my $connection = shift;

    my $sql = "select epigenome_id
            from epigenome
            left join object_xref on (epigenome.epigenome_id = object_xref.ensembl_id and ensembl_object_type = 'epigenome')
            left join xref using (xref_id)
            join external_db using (external_db_id)
            where dbprimary_acc = ? and external_db.db_name = ?
            order by epigenome_id DESC";

    my $sth = $connection->prepare($sql);
    $sth->bind_param( 1, $dbprimary_acc );
    $sth->bind_param( 2, $db_name );
    $sth->execute();
    my @row = $sth->fetchrow_array;
    return $row[0];
}



sub get_duplicate_ihecs{
    my $connection = shift;
    
    my $sql = "select dbprimary_acc
                from epigenome
                left join object_xref on (epigenome.epigenome_id = object_xref.ensembl_id and ensembl_object_type = 'epigenome')
                left join xref using (xref_id)
                join external_db using (external_db_id)
                where db_name = 'EPIRR'
                group by dbprimary_acc
                having count(1)> 1";
    my $sth = $connection->prepare($sql);
    $sth->execute();

    my @duplicated_ihecs;
    my @row;
    while (@row = $sth->fetchrow_array) {
        push (@duplicated_ihecs, $row[0]);
    }

    return @duplicated_ihecs

}

sub usage {
    my $usage = << 'END_USAGE';

Usage: delete_duplicate_epigenomes.pl -c <config_file>

Options:
-c config_file:			configuration file that contains the database connection details
-h help:				help message
 
END_USAGE

    say $usage;

    return 1;
}