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

GetOptions('c=s' => \$cfgFile);

# -----------------
# check parameters
# -----------------           
if (!$cfgFile) {
	   exit;
}

my $cfg = Config::Tiny->read($cfgFile);

my $connection = DBI->connect("DBI:mysql:".$cfg->{efg_db}->{dbname}.":".$cfg->{efg_db}->{host}.":".$cfg->{efg_db}->{port}, $cfg->{efg_db}->{user}, $cfg->{efg_db}->{pass}, {RaiseError => 1});

my $release = $cfg->{general}->{release};

# -----------------------------------
# get experiments with 'new_' prefix 
# -----------------------------------
my $sql = "select experiment_id, feature_type_id, epigenome_id, experimental_group_id, is_control
		   from experiment 
		   where experiment.name like 'new_%'";

print "get experiments with 'new_' prefix ...\n";

my $sth = $connection->prepare($sql);     
$sth->execute();
my @row;
my @experiment_rows;
while (@row = $sth->fetchrow_array) {
	my %experiment_data;
	$experiment_data{id} = @row[0];
	$experiment_data{is_control}=@row[4];
	$experiment_data{new_name} = get_experiment_new_name(@row[1], @row[2], @row[3], @row[4], $release, $connection, \@experiment_rows);

	push @experiment_rows, \%experiment_data;
}


# -----------------------------------
# update experiments with the new name
# -----------------------------------

foreach my $update_data (@experiment_rows){

	$sql = "update experiment 
		set experiment.name = ?
		where experiment_id = ?";
	
	$sth = $connection->prepare($sql);
	$sth->bind_param( 1, $update_data->{new_name});
	$sth->bind_param( 2, $update_data->{id});

	$sth->execute();

}

# -----------------------------------
# END
# -----------------------------------

print "FIN\n";

sub get_experiment_new_name {

	my $feature_type_id = shift;
	my $epigenome_id = shift;
	my $experimental_group_id = shift;
	my $is_control = shift;
	my $release = shift;
	my $connection = shift;
	my $experiment_rows = shift;

	my $new_name;

	my $analysis_name = get_analysis_name($feature_type_id, $connection);

	if ($is_control == 1){
		$new_name = create_control_experiment_name($feature_type_id, $epigenome_id, $experimental_group_id, $analysis_name, $release, $connection, $experiment_rows);
	}
	if ($is_control == 0){
		$new_name = create_signal_experiment_name($feature_type_id, $epigenome_id, $experimental_group_id, $analysis_name, $release, $connection);
	}

	return $new_name;

}

sub get_analysis_name {
	my $feature_id = shift;
	my $connection = shift;

	my $sql = "select analysis.logic_name
		   		from feature_type join analysis using (analysis_id)
		   		where feature_type_id = ?";

	my $sth = $connection->prepare($sql);
	$sth->bind_param( 1, $feature_id);     
	$sth->execute();
	my @row = $sth->fetchrow_array;

	return @row[0];

}

sub get_feature_name {
	my $feature_id = shift;
	my $connection = shift;

	my $sql = "select feature_type.name
		   		from feature_type 
		   		where feature_type_id = ?";


	my $sth = $connection->prepare($sql);
	$sth->bind_param( 1, $feature_id);     
	$sth->execute();
	my @row = $sth->fetchrow_array;

	return @row[0];

}

sub get_epigenome_name {
	my $epigenome_id = shift;
	my $connection = shift;

	my $sql = "select epigenome.name
		   		from epigenome 
		   		where epigenome_id = ?";


	my $sth = $connection->prepare($sql);
	$sth->bind_param( 1, $epigenome_id);     
	$sth->execute();
	my @row = $sth->fetchrow_array;

	return @row[0];

}

sub get_experimental_group_name {
	my $experimental_group_id = shift;
	my $connection = shift;

	my $sql = "select experimental_group.name
		   		from experimental_group 
		   		where experimental_group_id = ?";


	my $sth = $connection->prepare($sql);
	$sth->bind_param( 1, $experimental_group_id);
	$sth->execute();
	my @row = $sth->fetchrow_array;

	return @row[0];

}

sub check_experiment_name{
	my $experiment_name = shift;
	my $connetion = shift;
	my $experiment_rows = shift;

	my $exists = 0;

	my $sql = "select experiment_id
		   		from experiment 
		   		where UPPER(TRIM(experiment.name)) = ?";


	my $sth = $connection->prepare($sql);
	$sth->bind_param( 1, uc ($experiment_name));     
	$sth->execute();
	if($sth->rows) {

  		$exists = 1;
	}

	if ($exists == 0){
		foreach my $exp_row (@{$experiment_rows}){
			if (uc ($experiment_name) eq uc($exp_row->{new_name})){
				$exists = 1;
				last;
			}
		}
	}

	return $exists;
}

sub create_control_experiment_name {
    my ( $feature_type_id, $epigenome_id, $experimental_group_id, $analysis_name, $release, $connection, $experiment_rows) = @_;

    my $number = 1;

    my $epigenome_name = get_epigenome_name($epigenome_id, $connection);
    my $feature_name = get_feature_name($feature_type_id, $connection);
    my $exp_group_name = get_experimental_group_name($experimental_group_id, $connection);

    my $experiment_name;
    my $experiment_exists =0;

    do {

        $experiment_name = $epigenome_name . '_' . $feature_name . '_' . $analysis_name .'_no' . $number . '_'.$exp_group_name;
 
        if ( $release ) {
            $experiment_name .= '_'.$release;
        }

        $experiment_name =~ s/\s//g;
        $number++;

        $experiment_exists = check_experiment_name($experiment_name, $connection, $experiment_rows);

    } while ( $experiment_exists == 1 );

    return $experiment_name;
}

sub create_signal_experiment_name {
    my ( $feature_type_id, $epigenome_id, $experimental_group_id, $analysis_name, $release, $connection ) = @_;

    my $epigenome_name = get_epigenome_name($epigenome_id, $connection);
    my $feature_name = get_feature_name($feature_type_id, $connection);
    my $exp_group_name = get_experimental_group_name($experimental_group_id, $connection);

    my $experiment_name = $epigenome_name . '_' . $feature_name . '_' . $analysis_name .'_'. $exp_group_name;

    if ( $release ) {
        $experiment_name .= '_'.$release;
    }

    $experiment_name =~ s/\s//g;

    return $experiment_name;
}



