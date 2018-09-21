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
use DBI;
use Getopt::Long qw(GetOptions);
use Config::Tiny;


    my $cfgFile;
    my $new_epigenomes_start_id;

    GetOptions(
                'c=s' => \$cfgFile,
                'i=s' => \$new_epigenomes_start_id);

    # -----------------
    # check parameters
    # -----------------           
    if (!$cfgFile) {
 	   exit;
    }


my $cfg = Config::Tiny->read($cfgFile);

my $connection = DBI->connect("DBI:mysql:".$cfg->{efg_db}->{dbname}.":".$cfg->{efg_db}->{host}.":".$cfg->{efg_db}->{port}, $cfg->{efg_db}->{user}, $cfg->{efg_db}->{pass}, {RaiseError => 1});

my @no_ihec_epigenomes = get_no_ihec_epigenomes($connection, $new_epigenomes_start_id);

foreach my $epi_no_ihec (@no_ihec_epigenomes){
	my %old_epigenome = %{$epi_no_ihec};
	my $old_epigenome_id = (keys %old_epigenome)[0];
	my $old_epigenome_name = $old_epigenome{$old_epigenome_id};

	my @corr_new_epi = get_relative_epigenomes($connection, $old_epigenome_name, $new_epigenomes_start_id);
	foreach my $rel_epi (@corr_new_epi){
		my %new_epigenome = %{$rel_epi};
		my $new_epigenome_id = (keys %new_epigenome)[0];
		my $new_epigenome_name = $new_epigenome{$new_epigenome_id};

		print $old_epigenome_id."\t".$old_epigenome_name."\t".$new_epigenome_id."\t".$new_epigenome_name."\n";
	}



}




sub get_no_ihec_epigenomes{
    my $connection = shift;
    my $start_id = shift;
    
    my $sql = "select 
   				distinct epigenome.epigenome_id, epigenome.name,
   				experimental_group.name 
			   from 
   				epigenome 
   			   	join experiment using (epigenome_id) 
               	join feature_type using (feature_type_id) 
   			   	join experimental_group using (experimental_group_id) 
   			   	left join object_xref on (epigenome_id = ensembl_id and ensembl_object_type = 'Epigenome') 
   			   	left join xref using (xref_id) 
   			   	left join external_db on (external_db.external_db_id = xref.external_db_id and db_name = 'EpiRR') 
			   where class in (
 				'Histone', 
 				'Open Chromatin',
 				'Transcription Factor',
 				'Polymerase'
				) 
				and external_db.external_db_id is null
				and epigenome_id <= ?";
    my $sth = $connection->prepare($sql);
    $sth->bind_param( 1, $start_id );
    $sth->execute();

    my @epi_no_ihec;
    my @row;
    while (@row = $sth->fetchrow_array) {
    	my %old_epigenome;
    	$old_epigenome{$row[0]} = $row[1];
        push (@epi_no_ihec, \%old_epigenome);
    }

    return @epi_no_ihec;

}

sub get_relative_epigenomes{
	my $connection = shift;
	my $epi = shift;
	my $start_id = shift;

	my $sql = "select epigenome.epigenome_id, epigenome.name 
			   from epigenome 
			   where epigenome_id > ".$start_id." 
			   and (epigenome.name like '%".$epi."%' or epigenome.description like '%".$epi."%')";

	my $sth = $connection->prepare($sql);
    
    $sth->execute();

    my @epi_no_ihec;
    my @row;
    while (@row = $sth->fetchrow_array) {
    	my %new_epigenome;
    	$new_epigenome{$row[0]} = $row[1];
        push (@epi_no_ihec, \%new_epigenome);
    }

    return @epi_no_ihec;
}


