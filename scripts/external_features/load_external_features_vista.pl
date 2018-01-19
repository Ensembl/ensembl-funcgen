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

=cut

# Load VISTA enhancers to a funcgen database

# Example usage:
# load_external_features_vista.pl \
#       -registry   <registry.pm file> \
#       -type      Vista \
#       -species   mus_musculus \
#       -clobber \
#       -tee \
#       -logfile  <log file> \
#       <VISTA enhancers file>

use strict;
use warnings;
use feature qw(say);

use Data::Dumper;
use Getopt::Long;
use Pod::Usage;
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor;

#my ($host, $user, $dnadb_host, $dnadb_user, $pass, $port, $dbname, $dnadb_name, $dnadb_port);
my ($registry, $species, $type, $clobber, $old_assembly, $new_assembly, $archive);
#my ($, $dnadb_pass, $release);

#$main Helper params
#Set here to avoid used once error
$main::_tee     = 0;
$main::_log_file   = undef;
$main::_debug_file = undef;

my @tmp_args = @ARGV;

GetOptions
  (
   "registry=s"    => \$registry,
   "species=s"     => \$species,
   "type=s",       => \$type,
   "clobber",      => \$clobber,  #koThis is not behaving like a boolean flag??
   "tee",          => \$main::_tee,
   "logfile=s"     => \$main::_log_file,
   "debugfile=s"   => \$main::_debug_file,
   "archive=s",    => \$archive,
   "old_assembly=s", => \$old_assembly,
   "new_assembly=s", => \$new_assembly,
   "help",         => \&usage,
  );

#   "host=s",       => \$host,
#   "dnadb_host=s"  => \$dnadb_host,
#   "user=s",       => \$user,
#   "dnadb_user=s", => \$dnadb_user,
#   "species=s",    => \$species,
#   "pass=s",       => \$pass,
#   "port=i",       => \$port,
#   "dbname=s",     => \$dbname,
#   "dnadb_name=s", => \$dnadb_name,
#   "dnadb_port=s"  => \$dnadb_port,
#   "dnadb_pass=s"  => \$dnadb_pass,
#   "type=s",       => \$type,
#   "tee",          => \$main::_tee,
#   "logfile=s"     => \$main::_log_file,
#   "debugfile=s"   => \$main::_debug_file,
#   "archive=s",    => \$archive,
#   "clobber",      => \$clobber,#This is not behaving like a boolean flag??
#   "help",         => \&usage,
#   "old_assembly=s", => \$old_assembly,
#   "new_assembly=s", => \$new_assembly,
#  )
#  or pod2usage( -exitval => 1,
#                -message => "Params are:\t@tmp_args"
#              );

print "load_external_features.pl @tmp_args\n";
my @files = @ARGV;

#convert usage to pod2usage and move usage docs to pod at top

#usage() if (!$host || !$user || !$dbname || !$type);

Bio::EnsEMBL::Registry->load_all($registry);

my $funcgen_db_adaptor = Bio::EnsEMBL::Registry->get_DBAdaptor($species, 'funcgen');
my $funcgen_dbc = $funcgen_db_adaptor->dbc;

if(! eval "{require Bio::EnsEMBL::Funcgen::Parsers::$type; 1}"){
  die("Did not find a parser module: Bio::EnsEMBL::Funcgen::Parsers::$type\n$@");
}

my $parser = "Bio::EnsEMBL::Funcgen::Parsers::$type"->new
  (
   -db            => $funcgen_db_adaptor,
   -clobber       => $clobber,
   -archive       => $archive,
   -no_disconnect => 1, #Should never need to disconnect with these imports
  );

$parser->parse_and_load(\@files, $old_assembly, $new_assembly);


sub usage {
  print <<EOF; exit(0);

Usage: perl $0 <options>

  -host    Database host to connect to

  -port    Database port to connect to

  -dbname  Database name to use

  -user    Database username

  -pass    Password for user

  -type    Type of regulatory features to upload; must be present in analysis table and
           have a correspoinding parser e.g. vista, cisred, miranda.

  -file    File to load data from; must be parseable by parser corresponding to type

  -clobber Delete all features from type feature_set(s) before loading

  -archive Specify a version suffix to move the old feature_set to e.g. 47 will move feature_set to feature_set_name_v47

  -old_assembly, -new_assembly If *both* these options are specified, co-ordinates will be
           projected from the old assembly to the new assembly before being stored. This
           relies on a suitable mapping between the two assemblies having already been
           carried out.

  -help    This message

EOF

}

