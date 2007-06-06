#!/software/bin/perl

=head1 NAME

generate_DAS_config.pl

=head1 SYNOPSIS

this script will ...

=head1 DESCRIPTION

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

use Getopt::Std;
my %opts;
getopts('h:p:u:P:d:', \%opts);

my $dbhost = $opts{h} || 'ens-genomics1';
my $dbport = $opts{p} || 3306;
my $dbuser = $opts{u} || 'ensadmin';
my $dbpass = $opts{P} || 'ensembl';
my $dbname = $opts{d} || 'sg_homo_sapiens_funcgen_45_36g';

$| = 1;

use Bio::EnsEMBL::Utils::Exception qw(verbose throw warning info);
use Bio::EnsEMBL::Analysis::Tools::Logger qw(logger_verbosity logger_info);

my $utils_verbosity = 'WARNING';
my $logger_verbosity = 'OFF';
verbose($utils_verbosity);
logger_verbosity($logger_verbosity);

use Bio::EnsEMBL::Registry;
Bio::EnsEMBL::Registry->load_registry_from_db
    (
     -host => 'ens-livemirror',
     -user => 'ensro',
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
port=9876
hostname=bc-9-1-03.internal.sanger.ac.uk
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
port              = 3306
dbname            = $dbname
username          = ensro
description       = [Homo sapiens] eFG feature
source            = SOURCE
type              = TYPE
category          = CATEGORY
feature_set       = $fset

EOC

}


1;
