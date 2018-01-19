#!/usr/bin/env perl

=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2018] EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

  Questions may also be sent to the Ensembl help desk at
  <http://www.ensembl.org/Help/Contact>.

=head1 NAME

 merge_and_index_collections.pl
  
=head1 SYNOPSIS

 merge_and_index_collections.pl [options]

=head1 OPTIONS

 Mandatory params:

  --dbname          
  --dbuser          
  --dbpass          
  --dbport            
  --dbhost            
  --data_dir        The data dir of the input seq_region .col files
  --result_set_name The name of the corresponding ResultSet

 Optional params:
  --force_overwrite Forces over-writing of existing output files
 

=head1 DESCRIPTION

B<This program> merges seq_region specific BLOB collection files into one indexed file.

=cut

# To do
# 1 Add target dir to copy to and updte dbfile_registry_dir?
# 2 Add some states to the DB?
# 3 Add pod2usage

use strict;
use warnings;
use Getopt::Long;
use Pod::Usage; 
use Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor;

my ($dbuser, $dbpass, $dbname, $dbport, $dbhost, $dnadb_user, $dnadb_host, $dnadb_port, $dnadb_name);
my ($force, $data_dir, $rset_name, $dnadb_pass);
my @tmp_args = @ARGV;

GetOptions 
 (
  'data_dir=s'   => \$data_dir,
  'result_set=s' => \$rset_name,
  
  #DB connection
  'dbname=s' => \$dbname,
  'dbpass=s' => \$dbpass,
  'dbport=i' => \$dbport,
  'dbhost=s' => \$dbhost,
  'dbuser=s' => \$dbuser,
			
  #DNADB connection
  'dnadb_user=s' => \$dnadb_user,
  'dnadb_host=s' => \$dnadb_host,
  'dnadb_port=i' => \$dnadb_port,
  'dnadb_name=s' => \$dnadb_name,
  'dnadb_pass=s' => \$dnadb_pass,

  #Run modes
  'force_overwrite' => \$force,
  'help'            => sub { pod2usage(-verbose => 2); },
  #todo add slices/skip_slices support here
) or pod2usage(
               -exitval => 1,
               -message => "Params are:\t@tmp_args"
              );

if ( !( $dbuser && $dbname && $dbhost) ) {
   die('You must provide mysql -dbuser, -dbhost and -dbname arguments');
}
#Leave all other mandatory arg checking to API


#Set up DBAdaptor 
my $db = Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor->new
 (
  -host       => $dbhost,
  -dbname     => $dbname,
  -user       => $dbuser,
  -pass       => $dbpass,
  -port       => $dbport,
  -dnadb_host => $dnadb_host,
  -dnadb_user => $dnadb_user,
  -dnadb_port => $dnadb_port,
  -dnadb_name => $dnadb_name,
  -dnadb_pass => $dnadb_pass,
);

#test connection
#todo move these to DBAdaptor::new with no_test_connection flag?
$db->dbc->db_handle;
$db->dnadb->dbc->db_handle;

my @slices    = @{$$db->dnadb->get_SliceAdaptor->fetch_all('toplevel', 
                                                            undef, 
                                                            undef, 
                                                            1)};#inc dups
my ($rset)       = @{$db->get_ResultSetAdaptorr->fetch_all_by_name($rset_name)};

print STDOUT "Merging slice collections for ResultSet:\t$rset_name\n";
eval { $db->get_ResultFeatureAdaptor->merge_and_index_slice_collections($rset, 
                                                                        \@slices, 
                                                                        $data_dir, 
                                                                        $force); };
                                                                                                                                                
if($@){
  pod2usage(
            -exitval => 1,
            -message => "$@\nFailed to merge_and_index_slice_collections",
           );
}                                                                        

print STDOUT "Finished.\n";


#Should add IMPORTED states here, this is done the CollectionWriter as we need access to the Importer

1;
