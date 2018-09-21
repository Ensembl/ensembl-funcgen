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
use REST::Client;
use JSON;
use Data::Dumper;



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
    $logger->info( 'Connecting to DB: ' . $cfg->{efg_db}->{dbname} . "\n", 0, 1 );

    my $connection = DBI->connect("DBI:mysql:".$cfg->{efg_db}->{dbname}.":".$cfg->{efg_db}->{host}.":".$cfg->{efg_db}->{port}, $cfg->{efg_db}->{user}, $cfg->{efg_db}->{pass}, {RaiseError => 1});

    # ---------------------------
    # get the new epigenomes
    # ---------------------------
    $logger->info( 'Getting New epigenomes ... ' .  "\n", 0, 1 );

    my @new_epigenomes = get_the_new_epigenomes($connection);

    # -------------------------------------------
    # get ENCODE reference epigegenome accessions
    # --------------------------------------------
    $logger->info( 'Getting reference epigegenome accessions ... ' .  "\n", 0, 1 );

    get_encode_accessions(\@new_epigenomes);

    # -------------------------------------------
    # create new names
    # --------------------------------------------
    $logger->info( 'Creating new Epigenome names ... ' .  "\n", 0, 1 );

    get_new_names($connection, \@new_epigenomes);

    # -------------------------------------------
    # Update names in DB
    # --------------------------------------------
    $logger->info( 'Updating new Epigenome names ... ' .  "\n", 0, 1 );

    update_names($connection, \@new_epigenomes);

    # -------------------------------------------
    # END
    # --------------------------------------------
    $logger->info( 'FIN' .  "\n", 0, 1 );


}

sub update_names{
    my $connection = shift;
    my $new_names = shift;

    foreach my $epigenome (@{$new_names}){

        #print $epigenome->{old_name}."\t".$epigenome->{new_name}."\n";


        my $sql = 'UPDATE epigenome
            SET epigenome.name = ? , epigenome.display_label =?, epigenome.description = ?, epigenome.production_name = ?
            WHERE epigenome.epigenome_id = ? ';

        my $sth = $connection->prepare($sql);
        $sth->bind_param( 1, $epigenome->{new_name} );
        $sth->bind_param( 2, $epigenome->{label} );
        $sth->bind_param( 3, $epigenome->{description} );
        $sth->bind_param( 4, $epigenome->{new_name} );
        $sth->bind_param( 5, $epigenome->{id} );
        $sth->execute();

    }

    return;

}

sub get_new_names{
    my $connection = shift;
    my $encode_accessions = shift;

    my %done_names;

    foreach my $accession_data (@{$encode_accessions}){
        
        my $accession = $accession_data->{accession};

        my $url = 'https://www.encodeproject.org/reference-epigenomes/'.$accession.'/?frame=embedded&format=json';
        my $headers = {Accept => 'application/json'};
        my $rEClient = REST::Client->new();
        $rEClient->GET($url, $headers);
        my $rEResponse = decode_json($rEClient->responseContent());

        my $name = $rEResponse->{'biosample_term_name'}[0];

        my $name_without_special_chars = get_name_without_special_chars($name);

        my $name_exists = check_name_already_exists($connection, $name_without_special_chars);

        my $num = (split '_', $accession_data->{old_name})[-1];

        if ($name_exists == 1 || $num > 1){
            $name_without_special_chars = $name_without_special_chars.'_'.$accession;
            $name = $name.'_'.$accession;
        }

        $accession_data->{new_name} = $name_without_special_chars;
        $accession_data->{label} = $name;
        $accession_data->{description} = $rEResponse->{'description'};

        if (exists $done_names{$name_without_special_chars}){
            if ($done_names{$name_without_special_chars} == 0){
                $accession_data->{new_name} = $accession_data->{new_name}.'_'.$accession;
                $accession_data->{label} =$accession_data->{label} .'_'.$accession;
                $done_names{$name_without_special_chars} = 1;
            }else{
                my $numeration = $done_names{$name_without_special_chars};
                $accession_data->{new_name} = $accession_data->{new_name}.'_'.$accession.'_'.$numeration;
                $accession_data->{label} =$accession_data->{label} .'_'.$accession.'_'.$numeration;
                $done_names{$name_without_special_chars} = $numeration + 1;
            }
        }else{
            $done_names{$name_without_special_chars} = 0;
        }

        

    }





    return;

}

sub check_name_already_exists{
    my $connection = shift;
    my $name = shift;
    
    my $sql = "select epigenome.name from epigenome where epigenome.name = ?";
    my $sth = $connection->prepare($sql);
    $sth->bind_param( 1, $name);
    $sth->execute();
    
    my $found = 0;
    while ($sth->fetch()){
        $found = 1;
        last;
    }

    return $found;

}

sub get_name_without_special_chars {
    my ($epigenome_name) = shift;

    $epigenome_name =~ s/\:/_/g;
    $epigenome_name =~ s/\+//g;
    $epigenome_name =~ s/\(//g;
    $epigenome_name =~ s/\)//g;
    $epigenome_name =~ s/\-/_/g;
    $epigenome_name =~ s/\./_/g;
    $epigenome_name =~ s/\//_/g;
    $epigenome_name =~ s/\,/_/g;
    $epigenome_name =~ s/ /_/g;

    return $epigenome_name;
}

sub get_encode_accessions{
    my $new_epigenomes = shift;
    
    foreach my $epigenome (@{$new_epigenomes}){
        my $epigenome_name = $epigenome->{old_name};
        my $encode_accession = (split '_', $epigenome_name)[-2];
        my %enc_accessions;
        $epigenome->{accession}=$encode_accession;;
    }

    return;

}

sub get_the_new_epigenomes{
    my $connection = shift;
    
    my $sql = "select epigenome_id, epigenome.name from epigenome where epigenome.name like 'new_%'";
    my $sth = $connection->prepare($sql);
    $sth->execute();

    my @epigenomes;
    my @row;
    while (@row = $sth->fetchrow_array) {
        my %new_epigenomes;
        $new_epigenomes{id}=@row[0];
        $new_epigenomes{old_name}=@row[1];
        push (@epigenomes, \%new_epigenomes);
    }

    return @epigenomes;

}


sub usage {
    my $usage = << 'END_USAGE';

Usage: delete_duplicate_epigenomes.pl -c <config_file>

Options:
-c config_file:         configuration file that contains the database connection details
-h help:                help message
 
END_USAGE

    say $usage;

    return 1;
}