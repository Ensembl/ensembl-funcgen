#!/sw/arch/bin/perl
##!/software/bin/perl

=head1 NAME

run_import_type.pl

=head1 SYNOPSIS

run_import_type.pl -pass ensembl

=head1 DESCRIPTION

This script imports feature and cell types that have been defined 
in the types hash. For cell types the class needs to be an empty string.

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
use Getopt::Long;

### Edit this hash according to your needs #####################################
#
# Note: for cell types no class is needed, use an empty string instead
my %types = 
    (#          type   description                        class
     #'HeLa' => ['Cell','Human Epithelial Carcinoma Cells',''],
     #'GM06990' => ['Cell','Human B-Lymphocyte Cells',''],
     #'U2OS' => ['Cell','Human Bone Osteosarcoma Epithelial Cells',''],
     'CD4' => ['Cell','Human CD4 T-Cells',''],
     #'IMR90' => ['Cell','Human Fetal Lung Fibroblast',''],
     #'HL-60' => ['Cell','Human Promyelotic Leukemia Cells',''],
     #'HepG2' => ['Cell','Human Hepatocellular Liver Carcinoma Cells',''],
     #'Lymphoblastoid' => ['Cell','Human Lymphoblastoid Cells',''],

#     'MEF' => ['Cell','Mouse Embryonic Fibroblasts (MEF)',''],
#     'ES' => ['Cell','Mouse Embryonic Stem Cells',''],
#     'ESHyb' => ['Cell','Mouse Embryonic Stem Cells (hybrid)',''],
#     'NPC' => ['Cell','Mouse Neural Progenitor Cells (NPC)',''],
#
#     'MCV6' => ['Cell','Mouse B Lymphocyte Cells, partially reprogrammed (MCV6)',''],
#     'MCV8' => ['Cell','Mouse B Lymphocyte Cells, partially reprogrammed (MCV8)',''],
#     'MCV81' => ['Cell','Mouse B Lymphocyte Cells, partially reprogrammed (MCV8.1)',''],
#     'BIV1' => ['Cell','Mouse B Lymphocyte Cells, partially reprogrammed (BIV1)',''],
#
#     'Brain' => ['Cell','Mouse Whole Brain Tissue',''],
#
#     'E14.undiff' => ['Cell','Mouse Undifferentiated ES Cells (male)',''],
#     'E14.diff' => ['Cell','Mouse Differentiated ES Cells (4 days at RA treated, male)',''],
#     'E14.eb' => ['Cell','Mouse Differentiated ES Cells (10 days grown to Embryoid Bodies, male)',''],
#
#     'LF2.undiff' => ['Cell','Mouse Undifferentiated ES Cells (female)',''],
#     'LF2.diff' => ['Cell','Mouse Differentiated ES Cells (4 days at RA treated, female)',''],
#     'LF2.eb' => ['Cell','Mouse Differentiated ES Cells (10 days grown to Embryoid Bodies, female)',''],
#
#     'XT67E1.undiff' => ['Cell','Mouse Undifferentiated ES Cells with deleted Xist allele (female)',''],
#     'XT67E1.diff' => ['Cell','Mouse Differentiated ES Cells with deleted Xist allele (4 days at RA treated, female)',''],
#     'XT67E1.eb' => ['Cell','Mouse Differentiated ES Cells with deleted Xist allele (10 days grown to Embryoid Bodies, female)',''],
#
#     #'CTCF' => ['Feature','CCCTC-Binding Factor','Insulator'],
     'DNase1' => ['Feature','DNase1 Hypersensitive Site','DNA'],
#     #'H2AK5ac' => ['Feature','Histone 2A Lysine 5 Acetylation','Histone'],
#     #'H2AK9ac' => ['Feature','Histone 2A Lysine 9 Acetylation','Histone'],
#     #'H2AZ' => ['Feature','Histone 2A variant Z','Histone'],
#     #'H2BK120ac' => ['Feature','Histone 2B Lysine 120 Acetylation','Histone'],
#     #'H2BK12ac' => ['Feature','Histone 2B Lysine 12 Acetylation','Histone'],
#     #'H2BK20ac' => ['Feature','Histone 2B Lysine 20 Acetylation','Histone'],
#     #'H2BK5ac' => ['Feature','Histone 2B Lysine 5 Acetylation','Histone'],
#     #'H2BK5me1' => ['Feature','Histone 2B Lysine 5 Mono-Methylation','Histone'],
#     'H3' => ['Feature','Histone 3','Histone'],
#     'H3ac' => ['Feature','Histone 3 Acetylation','Histone'],
#     'H3K14ac' => ['Feature','Histone 3 Lysine 14 Acetylation','Histone'],
#     'H3K18ac' => ['Feature','Histone 3 Lysine 18 Acetylation','Histone'],
#     'H3K23ac' => ['Feature','Histone 3 Lysine 23 Acetylation','Histone'],
#     'H3K27ac' => ['Feature','Histone 3 Lysine 27 Acetylation','Histone'],
#     'H3K27me1' => ['Feature','Histone 3 Lysine 27 Mono-Methylation','Histone'],
#     'H3K27me2' => ['Feature','Histone 3 Lysine 27 Di-Methylation','Histone'],
#     'H3K27me3' => ['Feature','Histone 3 Lysine 27 Tri-Methylation','Histone'],
#     'H3K36ac' => ['Feature','Histone 3 Lysine 36 Acetylation','Histone'],
#     'H3K36me1' => ['Feature','Histone 3 Lysine 36 Mono-Methylation','Histone'],
     'H3K36me3' => ['Feature','Histone 3 Lysine 36 Tri-Methylation','Histone'],
#     'H3K4ac' => ['Feature','Histone 3 Lysine 4 Acetylation','Histone'],
#     'H3K4me1' => ['Feature','Histone 3 Lysine 4 Mono-Methylation','Histone'],
#     'H3K4me2' => ['Feature','Histone 3 Lysine 4 Di-Methylation','Histone'],
     'H3K4me3' => ['Feature','Histone 3 Lysine 4 Tri-Methylation','Histone'],
#     'H3K79me1' => ['Feature','Histone 3 Lysine 79 Mono-Methylation','Histone'],
#     'H3K79me2' => ['Feature','Histone 3 Lysine 79 Di-Methylation','Histone'],
#     'H3K79me3' => ['Feature','Histone 3 Lysine 79 Tri-Methylation','Histone'],
#     'H3K9ac' => ['Feature','Histone 3 Lysine 9 Acetylation','Histone'],
#     'H3K9me1' => ['Feature','Histone 3 Lysine 9 Mono-Methylation','Histone'],
#     'H3K9me2' => ['Feature','Histone 3 Lysine 9 Di-Methylation','Histone'],
#     'H3K9me3' => ['Feature','Histone 3 Lysine 9 Tri-Methylation','Histone'],
#     'H3R2me1' => ['Feature','Histone 3 Arginine 2 Mono-Methylation','Histone'],
#     'H3R2me2' => ['Feature','Histone 3 Arginine 2 Di-Methylation','Histone'],
#     'H4ac' => ['Feature','Histone 4 Acetylation','Histone'],
#     'H4K12ac' => ['Feature','Histone 4 Lysine 12 Acetylation','Histone'],
#     'H4K16ac' => ['Feature','Histone 4 Lysine 16 Acetylation','Histone'],
#     'H4K20me1' => ['Feature','Histone 4 Lysine 20 Mono-Methylation','Histone'],
#     'H4K20me3' => ['Feature','Histone 4 Lysine 20 Tri-Methylation','Histone'],
#     'H4K5ac' => ['Feature','Histone 4 Lysine 5 Acetylation','Histone'],
#     'H4K8ac' => ['Feature','Histone 4 Lysine 8 Acetylation','Histone'],
#     'H4K91ac' => ['Feature','Histone 4 Lysine 91 Acetylation','Histone'],
#     'H4R3me2' => ['Feature','Histone 4 Arginine 3 Di-Methylation','Histone'],
#     'PolII' => ['Feature','RNA Polymerase II','Polymerase'],
#     'WCE' => ['Feature','Mouse Whole Cell Extract','DNA'],
     );

################################################################################

my ($pass,$port,$host,$user,$dbname,$species,$help,$man,
    $version,$assembly,$data_version,$outdir);

$host = $ENV{EFG_HOST};
$port = $ENV{EFG_PORT};
$user = $ENV{EFG_WRITE_USER};
$species = $ENV{SPECIES};

$version=$ENV{VERSION};
$assembly=$ENV{ASSEMBLY};
$data_version = $ENV{DATA_VERSION};
$dbname = $ENV{EFG_DBNAME};

GetOptions (
            "pass|p=s"         => \$pass,
            "port=s"           => \$port,
            "host|h=s"         => \$host,
            "user|u=s"         => \$user,
            "dbname|d=s"       => \$dbname,
            "species=s"        => \$species,
            "data_version|v=s" => \$data_version,
            "help|?"           => \$help,
            "man|m"            => \$man,
            "outdir|o=s"     => \$outdir
            );

### defaults ###
$port = 3306 if !$port;
$species = 'homo_sapiens' if !$species;

### check options ###

throw("Must specify mandatory database hostname (-host).\n") if ! defined $host;
throw("Must specify mandatory database username. (-user)\n") if ! defined $user;
throw("Must specify mandatory database password (-pass).\n") if ! defined $pass;
throw("Must specify mandatory database name (-dbname).\n") if ! defined $dbname;
throw("Must specify mandatory database data version, like 47_36i (-data_version).\n") 
    if !$data_version;

$| = 1;

use Bio::EnsEMBL::Utils::Exception qw(throw);
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Funcgen::Utils::EFGUtils qw(open_file);
use Bio::EnsEMBL::Funcgen::RegulatoryFeature;

# use ensembldb as we may want to use an old version

my $cdb = Bio::EnsEMBL::DBSQL::DBAdaptor->new
    (
     -host => $ENV{CORE_HOST},
     -port => $ENV{CORE_PORT},
     -user => $ENV{CORE_USER},
     #-host => 'ens-staging',
     #-port => 3306,
     #-user => 'ensro',
     -dbname => $ENV{CORE_DBNAME},
     -species => $species,
     );

my $db = Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor->new
    (
     -host   => $host,
     -user   => $user,
     -dbname => $dbname,
     -species => $species,
     -pass   => $pass,
     -port   => $port,
     -dnadb  => $cdb
     );

#print Dumper $db;


foreach my $ft (keys %types) {
   
    my $cmd = 
        "$ENV{EFG_SRC}/scripts/import/import_type.pl ". 
        "-user $user  ".
        "-host $host  ".
        "-port $port ".
        "-dbname $dbname ".
        "-type $types{$ft}[0]Type  ".
        "-name $ft ".
        "-description '$types{$ft}[1]' ".
        "-class  '$types{$ft}[2]' ".
        "-pass $pass";
    #print Dumper $cmd;
    system($cmd) == 0
        or throw("Can't execute command $cmd");
}

1;
