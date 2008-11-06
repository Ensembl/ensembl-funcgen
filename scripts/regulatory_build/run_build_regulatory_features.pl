#!/software/bin/perl
###!/usr/bin/perl

=head1 NAME

run_build_regulatory_features.pl -- wrapper script to run the
build_regulatory_features.pl to create "Ensembl Regulatory Build"

=head1 SYNOPSIS

run_build_regulatory_features.pl -host host -user user -pass password 


Options:

  Mandatory
    -host|h          Host for eFG DB
    -user|u          User for eFG DB
    -pass|p          Password for eFG DB
    -dbname|d        Name of eFG DB
    -data_version|v  Version of data in eFG DB (e.g. 51_36m)
    -outdir|o        Name of outputut directory

    -cdbhost
    -cdbport
    -cdbuser
    -cdbpass

  Optional
    -jobs
    -help|?
    -man|m

=head1 DESCRIPTION


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
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor;

use Getopt::Long;

my ($host,$port,$user,$pass,$dbname,$species,$data_version,
    $cdbhost,$cdbport,$cdbuser,$cdbpass,$jobs,
    $help,$man);

$|=1;

# database parameters for eFG DB
$host = $ENV{EFG_HOST};
$port = $ENV{EFG_PORT};
$user = $ENV{EFG_WRITE_USER};
$dbname = $ENV{EFG_DBNAME};
$species = $ENV{SPECIES};
$data_version = $ENV{DATA_VERSION};

# default core database 
$cdbhost = 'ens-staging';
$cdbport = 3306;
$cdbuser = 'ensro';
#$cdbhost = 'ensembldb.ensembl.org';
#$cdbuser = 'anonymous';
#$cdbport = 3306;

# build specififc paramters
my $outdir = '/lustre/work1/ensembl/graef/RegBuild/v52';

my @focus = qw 
    (CD4_CTCF_e52 CD4_DNASE_IMPORT GM06990_DNASE_IMPORT 
     Nessie_NG_STD_2_ctcf_ren_BR1);

my @attributes = qw
    (
     CD4_H2AZ_e52 CD4_H2BK5me1_e52 CD4_H3K27me1_e52 CD4_H3K27me2_e52 
     CD4_H3K27me3_e52 CD4_H3K36me1_e52 CD4_H3K36me3_e52 CD4_H3K4me1_e52 
     CD4_H3K4me2_e52 CD4_H3K4me3_e52 CD4_H3K79me1_e52 CD4_H3K79me2_e52 
     CD4_H3K79me3_e52 CD4_H3K9me1_e52 CD4_H3K9me2_e52 CD4_H3K9me3_e52 
     CD4_H3R2me1_e52 CD4_H3R2me2_e52 CD4_H4K20me1_e52 CD4_H4K20me3_e52 
     CD4_H4R3me2_e52 CD4_PolII_e52 
     Wiggle_H3K27me3 Wiggle_H3K36me3 Wiggle_H3K4me3 Wiggle_H3K79me3 
     Wiggle_H3K9me3 Wiggle_H4K20me3 
     CD4_H2AK5ac_e52 CD4_H2AK9ac_e52 CD4_H2BK120ac_e52 CD4_H2BK12ac_e52 
     CD4_H2BK20ac_e52 CD4_H2BK5ac_e52 CD4_H3K14ac_e52 CD4_H3K18ac_e52 
     CD4_H3K23ac_e52 CD4_H3K27ac_e52 CD4_H3K36ac_e52 CD4_H3K4ac_e52 
     CD4_H3K9ac_e52 CD4_H4K12ac_e52 CD4_H4K16ac_e52 CD4_H4K5ac_e52 
     CD4_H4K8ac_e52 CD4_H4K91ac_e52
     );

#     CD4_H2AZ_e51 CD4_H2BK5me1_e51 CD4_H3K27me1_e51 CD4_H3K27me2_e51 
#     CD4_H3K27me3_e51 CD4_H3K36me1_e51 CD4_H3K36me3_e51 CD4_H3K4me1_e51 
#     CD4_H3K4me2_e51 CD4_H3K4me3_e51 CD4_H3K79me1_e51 CD4_H3K79me2_e51 
#     CD4_H3K79me3_e51 CD4_H3K9me1_e51 CD4_H3K9me2_e51 CD4_H3K9me3_e51 
#     CD4_H3R2me1_e51 CD4_H3R2me2_e51 CD4_H4K20me1_e51 CD4_H4K20me3_e51 
#     CD4_H4R3me2_e51 CD4_PolII_e51 
#     Wiggle_H3K27me3 Wiggle_H3K36me3 Wiggle_H3K4me3 Wiggle_H3K79me3 
#     Wiggle_H3K9me3 Wiggle_H4K20me3 
#     CD4_H2AK5ac_e51 CD4_H2AK9ac_e51 CD4_H2BK120ac_e51 CD4_H2BK12ac_e51 
#     CD4_H2BK20ac_e51 CD4_H2BK5ac_e51 CD4_H3K14ac_e51 CD4_H3K18ac_e51 
#     CD4_H3K23ac_e51 CD4_H3K27ac_e51 CD4_H3K36ac_e51 CD4_H3K4ac_e51 
#     CD4_H3K9ac_e51 CD4_H4K12ac_e51 CD4_H4K16ac_e51 CD4_H4K5ac_e51 
#     CD4_H4K8ac_e51 CD4_H4K91ac_e51



GetOptions (
            "host|h=s"       => \$host,
            "port=s"         => \$port,
            "user|u=s"       => \$user,
            "pass|p=s"       => \$pass,
            "dbname|d=s"     => \$dbname,
            "species=s"      => \$species,
            "data_version|v=s" => \$data_version,
            "cdbhost=s"      => \$cdbhost,
            "cdbport=s"      => \$cdbport,
            "cdbuser=s"      => \$cdbuser,
            "cdbpass=s"      => \$cdbpass,
            "jobs=s"         => \$jobs,
            "help|?"         => \$help,
            "man|m"          => \$man,
            );

### check options ###

throw("Must specify mandatory database hostname (-host).\n") if ! defined $host;
throw("Must specify mandatory database username. (-user)\n") if ! defined $user;
throw("Must specify mandatory database password (-pass).\n") if ! defined $pass;
throw("Must specify mandatory database name (-dbname).\n") if ! defined $dbname;
throw("Must specify mandatory database data version, like 47_36i (-data_version).\n") 
     if !$data_version;

my $cdb = Bio::EnsEMBL::DBSQL::DBAdaptor->new
    (
     -host => $cdbhost,
     -port => $cdbport,
     -user => $cdbuser,
     -dbname => $species.'_core_'.$data_version,
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
     -dnadb  => $cdb,
     );

# determine the number of toplevel slices we can analyze 
# (excluding mitochondria sequences (MT)!!!)
# this bit needs to be moved t the pipeline helper that it also can be called 
# by the actual run script

my @slices = ();
my $sa = $db->get_SliceAdaptor();
#foreach my $s (sort {$a->name cmp $b->name} @{$sa->fetch_all('toplevel', undef, 1)}) {
foreach my $s (sort {$a->name cmp $b->name} @{$sa->fetch_all('toplevel')}) {
        
    next if ($s->seq_region_name =~ m/^MT/);
    push @slices, $s;
    
}

#map {print Dumper($_->name, $_->get_seq_region_id);} @slices;

warn("There are ".scalar @slices." slices available to analyse.");

$jobs = '1-'.scalar(@slices) unless (defined $jobs);


# parse both focus and attribute sets and check that they exist
my (%focus_fsets, %attrib_fsets);
my $fsa = $db->get_FeatureSetAdaptor();

map { my $fset = $fsa->fetch_by_name($_);
      throw("Focus set $_ does not exist in the DB") 
          if (! defined $fset); 
      $focus_fsets{$fset->dbID} = $fset; 
  } @focus;
#print Dumper %focus_fsets;

map { 
    my $fset = $fsa->fetch_by_name($_);
    throw("Attribute set $_ does not exist in the DB") 
        if (! defined $fset); 
    $attrib_fsets{$fset->dbID()} = $fset; 
} @attributes;
#print Dumper %attrib_fsets;

# make sure that attribute sets also contain focus sets (Do we really need this?)
map { $attrib_fsets{$_} = $focus_fsets{$_} } keys %focus_fsets;

### build and import regbuild strings by feature_set_id and feature_type_id

my ($sql);

my $sth = $db->dbc->prepare("select * from meta where meta_key='regbuild.feature_set_ids'");
$sth->execute();
my $noe = scalar @{$sth->fetchall_arrayref};

if ($noe == 0) {

    $sql = 'insert into meta (meta_key, meta_value) values ("regbuild.feature_set_ids", "'
        .join(',', map {$_->dbID} sort {$a->name cmp $b->name} values %attrib_fsets).'");';
    #print "$sql\n";
    eval {
        $db->dbc->do($sql);
    };
    throw("Couldn't store regbuild.feature_set_ids in meta table.") if ($@);

} else {

    throw("regbuild.feature_set_ids already exists in meta table. You need to move it\n"
          ."  away before you can run the regulatory build.");
    
}


$sth = $db->dbc->prepare("select * from meta where meta_key='regbuild.feature_type_ids'");
$sth->execute();
$noe = scalar @{$sth->fetchall_arrayref};

if ($noe == 0) {

    $sql = 'insert into meta (meta_key, meta_value) values ("regbuild.feature_type_ids", "'
        .join(',', map {$_->feature_type->dbID} sort {$a->name cmp $b->name} values %attrib_fsets).'");';
    #print "$sql\n";
    eval {
        $db->dbc->do($sql);
    };
    throw("Couldn't store regbuild.feature_type_ids in meta table.") if ($@);

} else {

    throw("regbuild.feature_set_ids already exists in meta table. You need to move it\n"
          ."  away before you can run the regulatory build.");
    
}


# submit jobs to the farm

(my $bsubhost = $host) =~ s/-/_/g;
my $bsub = "bsub -oo '$outdir/RegulatoryBuild_%I.log' -J 'RegulatoryBuild_[".$jobs."]' ".
    "-R 'select[my".$bsubhost."<80] rusage[my".$bsubhost."=10:duration=10]'";
my $cmd = "$ENV{EFG_SRC}/scripts/regulatory_build/build_regulatory_features.pl"
    .' -host '.$host
    .' -port '.$port
    .' -user '.$user
    .' -pass '.$pass
    .' -dbname '.$dbname
    .' -data_version '.$data_version
    .' -focus '.join(',', @focus)
    .' -attrib '.join(',', @attributes)
    .' -outdir '.$outdir
    .' -dump_regulatory_features '
    .' -clobber '
    .' -stats '
    .' -write_features'
    ;

print $bsub, "    ", $cmd, "\n";

system("$bsub $cmd") &&
    throw ("Can't submit job array to farm");

__END__

#!/bin/sh

if [ $# -lt 4 ]; then

	echo "Usage: $0 <host> <port> <user> <password>"
	exit;

fi

HOST=$1
shift
PORT=$1
shift
USER=$1
shift
PASS=$1
shift

DATA_VERSION='51_36m'
DBNAME='sg_homo_sapiens_funcgen_'${DATA_VERSION}

OUTDIR='/lustre/work1/ensembl/graef/RegBuild/v51'

FOCUS=\
CD4_CTCF_e51,\
CD4_DNASE_IMPORT,\
GM06990_DNASE_IMPORT,\
Nessie_NG_STD_2_ctcf_ren_BR1

TARGET=\
CD4_H2AZ_e51,\
CD4_H2BK5me1_e51,\
CD4_H3K27me1_e51,\
CD4_H3K27me2_e51,\
CD4_H3K27me3_e51,\
CD4_H3K36me1_e51,\
CD4_H3K36me3_e51,\
CD4_H3K4me1_e51,\
CD4_H3K4me2_e51,\
CD4_H3K4me3_e51,\
CD4_H3K79me1_e51,\
CD4_H3K79me2_e51,\
CD4_H3K79me3_e51,\
CD4_H3K9me1_e51,\
CD4_H3K9me2_e51,\
CD4_H3K9me3_e51,\
CD4_H3R2me1_e51,\
CD4_H3R2me2_e51,\
CD4_H4K20me1_e51,\
CD4_H4K20me3_e51,\
CD4_H4R3me2_e51,\
CD4_PolII_e51,\
Wiggle_H3K27me3,\
Wiggle_H3K36me3,\
Wiggle_H3K4me3,\
Wiggle_H3K79me3,\
Wiggle_H3K9me3,\
Wiggle_H4K20me3,\
CD4_H2AK5ac_e51,\
CD4_H2AK9ac_e51,\
CD4_H2BK120ac_e51,\
CD4_H2BK12ac_e51,\
CD4_H2BK20ac_e51,\
CD4_H2BK5ac_e51,\
CD4_H3K14ac_e51,\
CD4_H3K18ac_e51,\
CD4_H3K23ac_e51,\
CD4_H3K27ac_e51,\
CD4_H3K36ac_e51,\
CD4_H3K4ac_e51,\
CD4_H3K9ac_e51,\
CD4_H4K12ac_e51,\
CD4_H4K16ac_e51,\
CD4_H4K5ac_e51,\
CD4_H4K8ac_e51,\
CD4_H4K91ac_e51

case "$MODE" in 
    "dump")


        # first dump annotated features from database 
        # for further processing

        bsub -o "$OUTDIR/RegulatoryBuild_afDump_$FOCUS.log" -J RegBuild_afDump \
        $EFG_SRC/scripts/regulatory_build/build_regulatory_features.pl \
            -host $HOST \
            -port $PORT \
            -user $USER \
            -pass $PASS \
            -dbname $DBNAME \
            -data_version $DATA_VERSION \
            -outdir $OUTDIR \
            -dump \
            -focus $FOCUS \
            -target $TARGET \
            $@
        ;;
    
    * )

        # build regulatory features and dump to file when dumping 
        # features has been sucessfully completed

        #            -w 'done("RegBuild_afDump")' \

        #bsub -o "$OUTDIR/RegulatoryBuild_$FOCUS.log" -J RegBuild_[1-25] \
        $EFG_SRC/scripts/regulatory_build/build_regulatory_features.pl \
            -host $HOST \
            -port $PORT \
            -user $USER \
            -pass $PASS \
            -dbname $DBNAME \
            -data_version $DATA_VERSION \
            -focus $FOCUS \
            -target $TARGET \
            -outdir $OUTDIR \
            -dump_features \
            -clobber \
            -stats \
            -write_features \
            $@

           # -gene_signature \
           # -seq_name 10 \
       ;;
esac
          
