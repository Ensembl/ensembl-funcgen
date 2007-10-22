#!/usr/bin/perl

=head1 NAME

generate_DAS_config.pl

=head1 SYNOPSIS

./generate_DAS_config.pl \
    -h dbhost \
    -p 3306 \
    -u ensro \
    -d homo_sapiens_funcgen_47_36i \
    -s homo_sapiens \
    -v 46_36h \
    -H DAShost \
    -P 9000

=head1 DESCRIPTION

Produces ProServer DAS server configuration file based on a given
eFG database and a list of links to automatically attach tracks.

=head1 LICENCE

This code is distributed under an Apache style licence. Please see
http://www.ensembl.org/info/about/code_licence.html for details.

=head1 AUTHOR

Stefan Graf, Ensembl Functional Genomics

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
getopts('h:p:u:w:d:s:v:H:P:lg:c:', \%opts);

my $dbhost = $opts{h} or throw("Need to specify dbhost via option -h!");
my $dbport = $opts{p} or throw("Need to specify dbport via option -p!");
my $dbuser = $opts{u} or throw("Need to specify dbuser via option -u!");
my $dbpass = $opts{w} || "";
my $dbname = $opts{d} or throw("Need to specify dbname via option -d!");
my $species = $opts{s} or throw("Need to specify spieces via option -s!");
my $data_version = $opts{v} or throw("Need to specify data version via option -v!");

my $dashost = $opts{H} or throw("Need to specify DAS server hostname via option -H!");
my $dasport = $opts{P} or throw("Need to specify DAS server hostname via option -P!");

my $location;
$location = 'gene=STAT1';
$location = 'gene='.$opts{g} if ($opts{g});
$location = 'c='.$opts{c} if ($opts{c});

$| = 1;

my %set = (
           ctcf_ren_BR1_TR1_ =>    { color => 'contigblue2', name => 'IMR90_CTCF' },
           Nessie_NG_STD_2_ctcf_ren_BR1 => { color => 'contigblue2', name => 'IMR90_CTCF_Nessie'},
           ctcf_ren_IMPORT_TileMap => { color => 'contigblue2', name => 'ctcf_ren_IMPORT_TileMap'},
           ctcf_ren_IMPORT_Chipotle => { color => 'contigblue2', name => 'ctcf_ren_IMPORT_Chipotle'},

           GM06990_DNASE_IMPORT => { color => 'contigblue2', name => 'GM06990_DNASE' },
           CD4_DNASE_IMPORT =>     { color => 'contigblue2', name => 'CD4_DNASE' },
           
           CD4_CTCF => { color => 'contigblue1', name => 'CD4_CTCF', plot => 'hist' },

           CD4_H2AZ => { color => 'contigblue1', name => 'CD4_H2AZ', plot => 'hist' },

           CD4_H2BK5me1 => { color => 'contigblue1', name => 'CD4_H2BK5me1', plot => 'hist' },

           CD4_H3K27me1 => { color => 'contigblue1', name => 'CD4_H3K27me1', plot => 'hist' },
           CD4_H3K27me2 => { color => 'contigblue1', name => 'CD4_H3K27me2', plot => 'hist' },
           CD4_H3K27me3 => { color => 'contigblue1', name => 'CD4_H3K27me3', plot => 'hist' },

           CD4_H3K36me1 => { color => 'contigblue1', name => 'CD4_H3K36me1', plot => 'hist' },
           CD4_H3K36me3 => { color => 'contigblue1', name => 'CD4_H3K36me3', plot => 'hist' },

           CD4_H3K4me1 => { color => 'contigblue1', name => 'CD4_H3K4me1', plot => 'hist' },
           CD4_H3K4me2 => { color => 'contigblue1', name => 'CD4_H3K4me2', plot => 'hist' },
           CD4_H3K4me3 => { color => 'contigblue1', name => 'CD4_H3K4me3', plot => 'hist' },

           CD4_H3K79me1 => { color => 'contigblue1', name => 'CD4_H3K79me1', plot => 'hist' },
           CD4_H3K79me2 => { color => 'contigblue1', name => 'CD4_H3K79me2', plot => 'hist' },
           CD4_H3K79me3 => { color => 'contigblue1', name => 'CD4_H3K79me3', plot => 'hist' },

           CD4_H3K9me1 => { color => 'contigblue1', name => 'CD4_H3K9me1', plot => 'hist' },
           CD4_H3K9me2 => { color => 'contigblue1', name => 'CD4_H3K9me2', plot => 'hist' },
           CD4_H3K9me3 => { color => 'contigblue1', name => 'CD4_H3K9me3', plot => 'hist' },

           CD4_H3R2me1 => { color => 'contigblue1', name => 'CD4_H3R2me1', plot => 'hist' },
           CD4_H3R2me2 => { color => 'contigblue1', name => 'CD4_H3R2me2', plot => 'hist' },

           CD4_H4K20me1 => { color => 'contigblue1', name => 'CD4_H4K20me1', plot => 'hist' },
           CD4_H4K20me3 => { color => 'contigblue1', name => 'CD4_H4K20me3', plot => 'hist' },

           CD4_H4R3me2 => { color => 'contigblue1', name => 'CD4_H4R3me2', plot => 'hist' },

           CD4_PolII => { color => 'contigblue1', name => 'CD4_PolII', plot => 'hist' },

           RegulatoryFeatures =>   { color => 'red3', name => 'RegulatoryFeatures' }
);


use Bio::EnsEMBL::Registry;
Bio::EnsEMBL::Registry->load_registry_from_db
    (
     -host => 'ensembldb.ensembl.org',
     -user => 'anonymous',
     #-verbose => "1"
     );
my $cdb = Bio::EnsEMBL::Registry->get_DBAdaptor("$species", 'core');

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

my $dsa = $db->get_DataSetAdaptor();
my $fsa = $db->get_FeatureSetAdaptor();


open(HTML, "> efg_das_$species.html")
    or throw "Can't open file";
    
&print_html_header();


&print_header();
my %seen = ();
foreach my $dset (sort {$a->name cmp $b->name} @{$dsa->fetch_all()}) {

    #print $dset->name, ' ', $dset->supporting_set_type, "\n";

    if ($dset->supporting_set_type() eq 'result') {
        
        &print_result_config($dset);
        
    }   
     
    &print_feature_config($dset);
        
}


&print_html_footer();

sub print_result_config() 
{
    my ($dset) = @_;

    my @rsets = sort {$a->name cmp $b->name} @{$dset->get_supporting_sets()};
    #map {print Dumper $_->name } @rsets;
    my $name = $rsets[0]->name;
    
    # don't duplicate entries
    return if exists $seen{$name};
    $seen{$name}='';

    $name =~ s/\s+/_/g;
    my $dsn = $name.'_'.$dset->supporting_set_type;
    my $dset_id = $dset->dbID;
    my $type = $dset->supporting_set_type;
    my $cell_type = $dset->cell_type->name;
    my $desc = "[$species - $cell_type] $dsn ($type set)";

    print <<EOC;
[$dsn]
state             = on
adaptor           = efg_result_set
transport         = efg_transport
host              = $dbhost
port              = $dbport
dbname            = $dbname
species           = $species
data_version      = $data_version
username          = $dbuser
description       = $desc
source            = SOURCE
type              = $type
cell_type         = $cell_type
category          = CATEGORY
data_set_id       = $dset_id

EOC


    if ($opts{l} && $set{$dsn}) {
        print HTML '<p><a href="'.
            'http://www.ensembl.org/'.$species.
            '/contigview?'.$location.';'.
            'add_das_source=(name='.$set{$dsn}{name}.'+url=http://'.$dashost.':'.$dasport.
            '/das+dsn='.$dsn.
            '+type=ensembl_location_chromosome'.
            '+color='.$set{$dsn}{color}.'+strand=r+labelflag=n'.
            '+group=n+score=s+fg_data=o+active=1)">'.$desc.'</a></p>'."\n";


    } else {

        print HTML '<p><a href="'.
            'http://www.ensembl.org/'.$species.
            '/contigview?'.$location.';'.
            'add_das_source=(name='.$dsn.'+url=http://'.$dashost.':'.$dasport.
            '/das+dsn='.$dsn.
            '+type=ensembl_location_chromosome'.
            '+color=red+strand=r+labelflag=n'.
            '+group=n+score=s+fg_data=o+active=1)">'.$desc.'</a></p>'."\n";


    }



}

sub print_feature_config() 
{
    my ($dset) = @_;

    my $type = $dset->supporting_set_type();
    my $fset = $dset->product_FeatureSet();

    my $cell_type = $fset->cell_type ? $fset->cell_type->name : '';

    my $name = $fset->name;
    $name =~ s/\s+/_/g;
    my $dsn = $name;
    my $fset_id = $fset->dbID;
    my $fset_type = $fset->type;
    my $desc = "[$species - $cell_type] $dsn ($fset_type feature set)";

    print <<EOC;
[$dsn]
state             = on
adaptor           = efg_feature_set
transport         = dbi
host              = $dbhost
port              = $dbport
dbname            = $dbname
species           = $species
data_version      = $data_version
username          = $dbuser
description       = $desc
source            = SOURCE
type              = $type
feature_set_type  = $fset_type
category          = CATEGORY
feature_set_id    = $fset_id

EOC

    if ($opts{l} && exists $set{$fset->name}) {
        
        my $plot ='';
        if (exists $set{$fset->name}{plot}){
            if ($set{$fset->name}{plot} eq 'cgrad') {
                $plot = '+score=c+fg_grades=100+fg_data=n+fg_max=100';
            } elsif ($set{$fset->name}{plot} eq 'hist') {
                $plot = '+score=h';
            } elsif ($set{$fset->name}{plot} eq 'tiling') {
                $plot = '+score=t';
            }
        }
        
        print HTML '<p><a href="'.
            'http://www.ensembl.org/'.$species.
            '/contigview?'.$location.';'.
            'add_das_source=(name='.$set{$fset->name}{name}.'+url=http://'.$dashost.':'.$dasport.
            '/das+dsn='.$dsn.
            '+type=ensembl_location_chromosome'.
            $plot.
            '+color='.$set{$fset->name}{color}.'+strand=r+labelflag=n'.
            '+group=n+active=1)">'.$desc.'</a></p>'."\n";

    } else {
    
        print HTML '<p><a href="'.
            'http://www.ensembl.org/'.$species.
            '/contigview?'.$location.';'.
            'add_das_source=(name='.$dsn.'+url=http://'.$dashost.':'.$dasport.
            '/das+dsn='.$dsn.
            '+type=ensembl_location_chromosome'.
            '+color=blue+strand=r+labelflag=n'.
            '+group=n+active=1)">'.$desc.'</a></p>'."\n";
        
        
    }

}

sub print_header() 
{
    print <<EOH;
[general]
prefork=2
maxclients=4
port=$dasport
hostname=$dashost
;response_hostname=das.example.com
;response_port=80
;response_protocol=https
;response_baseuri=/frontend
;oraclehome=/usr/local/oracle
;ensemblhome=/usr/local/ensembl
pidfile=/nfs/acari/graef/src/ensembl-functgenomics/DAS/log/efg-das_${species}.pid
logfile=/nfs/acari/graef/src/ensembl-functgenomics/DAS/log/efg-das_${species}.log
 
EOH
}

sub print_html_header() 
{

    print HTML <<EOHTML;
<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">
<html xmlns="http://www.w3.org/1999/xhtml" xml:lang="en-gb"  lang="en-gb">
<head><title>Ensembl Functional Genomics DAS links</title></head>
<body>
<h3>DAS data sources for</h3>
<h2>$dbname on $dashost:$dasport</h2>
EOHTML

print HTML '<p>selected location: <a href="'.
'http://www.ensembl.org/'.$species.
'/contigview?'.$location.';">'.$location.'</a></p>'."\n";

}

sub print_html_footer() 
{

    print HTML '</body></html>';

}


1;
