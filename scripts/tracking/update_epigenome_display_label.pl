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
if (!$cfgFile || !$new_epigenomes_start_id) {
	   exit;
}

my $cfg = Config::Tiny->read($cfgFile);

my $connection = DBI->connect("DBI:mysql:".$cfg->{efg_db}->{dbname}.":".$cfg->{efg_db}->{host}.":".$cfg->{efg_db}->{port}, $cfg->{efg_db}->{user}, $cfg->{efg_db}->{pass}, {RaiseError => 1});

my $sql = "select epigenome.epigenome_id, epigenome.name
		   from epigenome 
		   where epigenome_id >= ?";


my $sth = $connection->prepare($sql);
$sth->bind_param( 1, $new_epigenomes_start_id);
$sth->execute();

my @row;
while (@row = $sth->fetchrow_array) {
	my $new_display_label = buid_new_display_label(@row[1]);
	#print $new_display_label."\n";
	update_display_label($connection, @row[0], $new_display_label);
}

# -------------------------------------------
# END
# --------------------------------------------
print 'FIN' .  "\n";


sub buid_new_display_label{
	my $epigenome_name = shift;
	$epigenome_name =~ s/_/ /g;
	$epigenome_name =~ s/  / /g;
	return $epigenome_name;
}

sub update_display_label{
	my $connection = shift;
	my $epigenome_id = shift;
	my $display_label = shift;

	my $sql = "Update epigenome
				set epigenome.display_label = ?
				where epigenome.epigenome_id = ?";

	my $sth = $connection->prepare($sql);
	$sth->bind_param( 1, $display_label );
	$sth->bind_param( 2, $epigenome_id );
	$sth->execute();

	return;
}