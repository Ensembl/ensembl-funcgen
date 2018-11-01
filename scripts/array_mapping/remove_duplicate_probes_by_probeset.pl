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
use Bio::EnsEMBL::Registry;
use Getopt::Long qw(GetOptions);
use DBI;

my $registry_file;
my $species;

GetOptions (
   'registry=s'           => \$registry_file,
   'species=s'            => \$species,
);

#------------------------
# get connection the DB
#------------------------
print "Connecting to the database ...\n";

Bio::EnsEMBL::Registry->load_all($registry_file);
my $db_adapt = Bio::EnsEMBL::Registry->get_DBAdaptor($species, 'funcgen');
my $connection = $db_adapt->dbc();

#---------------------------------------------
# get arrays with duplicate probes by probeset
#---------------------------------------------
print "Searching for arrays with duplicate probes by probeset ...\n";

my $sql = "select distinct array_chip.array_chip_id
			from array join array_chip using (array_id) join probe using (array_chip_id) 
			where array.is_probeset_array = true 
			group by array.name, probe.name, probe.probe_set_id 
			having count(distinct probe.probe_id)>1";

my $sth = $connection->prepare($sql);
$sth->execute();
my @row;
my @array_chip_id;
while (@row = $sth->fetchrow_array) {
    push (@array_chip_id, $row[0]);
}

my $array_count = scalar @array_chip_id;
print "Number of arrays found with duplicate probes by probeset: " . $array_count . "\n";

if ($array_count > 0){
	#---------------------------------------------
	# get probes with duplicates by probeset
	#---------------------------------------------
	print "Searching for duplicate probes by probeset in the arrays found ...\n";

	$sql = "select probe.name, probe.probe_set_id, probe.probe_seq_id
		from probe
		where array_chip_id in(" . join (',', map { '?' } @array_chip_id) . ")";

	$sql .=	" group by probe.name, probe.probe_set_id,  probe.probe_seq_id
		having count(distinct probe.probe_id)>1";

	print "Query: ".$sql."\n";

	$sth = $connection->prepare($sql);
	my $count = 1;
	foreach my $array (@array_chip_id){
		$sth->bind_param( $count, $array);
		$count ++;
	}

	$sth->execute();

	my @duplicate_probes;
	while (@row = $sth->fetchrow_array) {
		my %duplicate_probe_data;
		$duplicate_probe_data{'name'}=$row[0];
		$duplicate_probe_data{'probe_set_id'}=$row[1];
		$duplicate_probe_data{'probe_seq_id'}=$row[2];
		push (@duplicate_probes, \%duplicate_probe_data);
	}

	my $probe_count = scalar @duplicate_probes;
	print "Number of duplicate probes by probeset: " . $probe_count . "\n";

	if ($probe_count > 0){
		#---------------------------------------------
		# select probes to be deleted
		#---------------------------------------------
		print "Searching for probes to be deleted ...\n";

		my @probes_to_delete;
		foreach my $duplicate (@duplicate_probes){
			$sql = "select probe_id
					from probe
					where probe.name = ? and probe_set_id = ? and probe_seq_id = ?";

			$sth = $connection->prepare($sql);
		    $sth->bind_param( 1, $duplicate->{'name'} );
		    $sth->bind_param( 2, $duplicate->{'probe_set_id'} );
		    $sth->bind_param( 3, $duplicate->{'probe_seq_id'} );
		    $sth->execute();

		    #skip first row
		    @row = $sth->fetchrow_array;

		    while (@row = $sth->fetchrow_array) {
		    	push (@probes_to_delete, $row[0]);
			}
		}

		#---------------------------------------------
		# delete duplicate probes
		#---------------------------------------------
		print "Deleting probes ...\n";

		my $deleted_counter = 0;
		foreach my $to_delete (@probes_to_delete){
			$sql = "delete from probe where probe_id = ?";
			$sth = $connection->prepare($sql);
			$sth->bind_param( 1, $to_delete );
			$sth->execute();

			$deleted_counter ++;
		}

		print "Number of deleted probes: " . $deleted_counter . "\n";

	}

}


























