#!/usr/bin/perl

=head1 NAME

generate_DAS_config.pl

=head1 SYNOPSIS

generate_DAS_config.pl -h ensembldb.ensembl.org -p 3306 -u anonymous -d homo_sapiens_funcgen_45_36g -H localhost -P 9000

=head1 DESCRIPTION

Produces ProServer DAS server configuration file based on a given
eFG database.

=head1 LICENCE

This code is distributed under an Apache style licence. Please see
http://www.ensembl.org/info/about/code_licence.html for details.

=head1 AUTHOR

Stefan Graf <graef@ebi.ac.uk>, Ensembl Functional Genomics

=head1 CONTACT

Please post comments/questions to the Ensembl development list
<ensembl-dev@ebi.ac.uk>

=cut

use strict;
use warnings;
use Data::Dumper;

use Bio::EnsEMBL::Utils::Exception qw(verbose throw warning info);
use Bio::EnsEMBL::Analysis::Tools::Logger qw(logger_verbosity logger_info);

my $utils_verbosity = 'WARNING';
my $logger_verbosity = 'OFF';
verbose($utils_verbosity);
logger_verbosity($logger_verbosity);

use Getopt::Std;
my %opts;
getopts('h:p:u:w:P:d:H:', \%opts);

my $dbhost = $opts{h} or throw("Need to specify dbhost via option -h!");
my $dbport = $opts{p} or throw("Need to specify dbport via option -p!");
my $dbuser = $opts{u} or throw("Need to specify dbuser via option -u!");
my $dbpass = $opts{w} || "";
my $dbname = $opts{d} or throw("Need to specify dbname via option -d!");

my $dashost = $opts{H} or throw("Need to specify DAS server hostname via option -H!");
my $dasport = $opts{P} or throw("Need to specify DAS server hostname via option -P!");

$| = 1;
$| = 1;

use Bio::EnsEMBL::Registry;
Bio::EnsEMBL::Registry->load_registry_from_db
    (
	 -host => 'ensembldb.ensembl.org',
	 -user => 'anonymous',
     #-verbose => "1"
     );
my $cdb = Bio::EnsEMBL::Registry->get_DBAdaptor('human', 'core');

my $db = Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor->new
    (
     -host   => $dbhost,
     -user   => $dbuser,
     -dbname => $dbname,
     -pass   => $dbpass,
     -port   => $dbport,
     -dnadb  => $cdb
     );

use Bio::EnsEMBL::Funcgen::FeatureSet;
my $fsa = $db->get_FeatureSetAdaptor();
my $fsets = $fsa->fetch_all();

#print Dumper $fsets;
my @fset_names = map $_->name(), @{$fsets};
#print Dumper @fset_names;

&print_header();
foreach my $fset (sort @fset_names) {

    my $name = $fset;
    $name =~ s/::.+// if ($fset =~ m/^Overlap/);
    &print_config($name, $fset);

}

sub print_header() 
{
    print <<EOH;
[general]
prefork=1
maxclients=1
port=$dasport
hostname=$dashost
;response_hostname=das.example.com
;response_port=80
;response_protocol=https
;response_baseuri=/frontend
;oraclehome=/usr/local/oracle
;ensemblhome=/usr/local/ensembl
pidfile=/nfs/acari/graef/src/ensembl-functgenomics/DAS/log/efg-das.pid
logfile=/nfs/acari/graef/src/ensembl-functgenomics/DAS/log/efg-das.log
 
EOH
}

sub print_config() 
{
    my ($name, $fset) = @_;

    print <<EOC;
[$name]
state             = on
adaptor           = efg_feature_set
transport         = dbi
host              = $dbhost
port              = $dbport
dbname            = $dbname
username          = $dbuser
description       = [Homo sapiens] eFG feature
source            = SOURCE
type              = TYPE
category          = CATEGORY
feature_set       = $fset

EOC

}


1;
