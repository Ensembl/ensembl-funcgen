#!/usr/bin/env perl

=head1 NAME

get_windows.pl -- 

=head1 SYNOPSIS

Generates slices (sequence windows) based on given window size (--win_size)
and optional shift (-shift; default: --window_size).

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
use Getopt::Long;

my ($pass,$port,$host,$user,$dbname,$species,$help,$man, $data_version,
	$debug, $file, $bin_size, $frag_length);

GetOptions (
            'pass|p:s'         => \$pass,
            'port:i'           => \$port,
            'host|h=s'         => \$host,
            'user|u=s'         => \$user,
            'dbname|d=s'       => \$dbname,
            'data_version|v=s' => \$data_version,
            'species=s'        => \$species,
            'help|?'           => \$help,
            'man|m'            => \$man,
            'debug'            => \$debug,
            'file=s'             => \$file,
            'bin_size=i'       => \$bin_size,
            'frag_length=i'    => \$frag_length,
            );

### defaults ###
#if (!$port) {
#	$port = 3306;
#	warn("No port specified, using default '$port'.")
#}
#if (!$species) {
#	$species = 'homo_sapiens';
#	warn("No species specified, using default '$species'.")
#}

### check options ###
throw("Must specify mandatory bin size (-bin_size).\n") 
    if ! defined $bin_size;
throw("Must specify mandatory fragment length (-frag_length).\n") 
    if ! defined $frag_length;

#throw("Must specify mandatory database name (-dbname).\n") if ! defined $dbname;
#throw("Must specify mandatory database data version, like 47_36i (-data_version).\n") 
#    if !$data_version;

$| = 1;

use Bio::EnsEMBL::Utils::Exception qw(verbose throw warning info);
#use Bio::EnsEMBL::DBSQL::DBAdaptor;
#use Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor;
#use Bio::EnsEMBL::Funcgen::Utils::EFGUtils qw(open_file);

#my $cdb = Bio::EnsEMBL::DBSQL::DBAdaptor->new
#    (
#     #-host => 'ensembldb.ensembl.org',
#     #-port => 3306,
#     #-user => 'anonymous',
#	 ### local forwards
#     -host => '127.0.0.1',
#     #-port => 33060, # ensembldb.ensembl.org (local forward)
#     -port => 33064, # ens-livemirror (local forward)
#     -user => 'ensro',
#     -dbname => $species.'_core_'.$data_version,
#     -species => $species,
#     );
#print Dumper $cdb;

#my $db = Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor->new
#    (
#     -host   => $host,
#     -user   => $user,
#     -dbname => $dbname,
#     -species => $species,
#     -pass   => $pass,
#     -port   => $port,
#     -dnadb  => $cdb
#     );
#print Dumper $db;

## Ensembl core and FG adaptors
#my $fsa = $db->get_FeatureSetAdaptor();
#my $dsa = $db->get_DataSetAdaptor();
#my $fta = $db->get_FeatureTypeAdaptor();
#my $afa = $db->get_AnnotatedFeatureAdaptor();
#my $rfa = $db->get_RegulatoryFeatureAdaptor();
#my $aa = $db->get_AnalysisAdaptor();

#my $sa = $cdb->get_SliceAdaptor();
#my $ga = $cdb->get_GeneAdaptor();
#my $ta = $cdb->get_TranscriptAdaptor();


my (@bin, $start_bin, $start_bin_start, $end_bin, $end_bin_start,
    $seq, $read_start, $read_end, $read_length, $ori, $read_extend);

# Get infile with features to project
if ($ENV{LSB_JOBINDEX}) {
    $file = $ARGV[$ENV{LSB_JOBINDEX}-1];
}
print $file, "\n";

throw("Must specify mandatory file name (-file) or".
      " run by wrapper script 'run_project_features'.\n")
    unless $file;

#open(FILE, $file)
open(FILE, "gzip -dc $file |")
    or throw ("Can't open file $file");

my $binsize = sprintf("%03d", $bin_size);
(my $out = $file) =~ s,_reads\.bed,_profile_${binsize}.bed,;

open(OUT, "| gzip -c > $out")
    or throw ("Can't open out file $out");

while (<FILE>) {

    chomp;
    my @col = split("\t");
    #print Dumper @col;

    if (defined $seq && $seq ne $col[0]) {

        &write_bins();

        @bin = ();

    }

    $seq = $col[0];
    $read_start = $col[1];
    $read_end = $col[2];
    $read_length = $read_end-$read_start+1;
    $ori = $col[5];

    throw("read is longer ($read_length) than specified fragment length ($frag_length)")
        if ($frag_length<$read_length);

    $read_extend = $frag_length-$read_length;
    #printf "%s\t%d\t%d\t%d\t%s\n", 
  	#$seq, $read_start, $read_end, $read_length, $ori;

    # extend reads to given fragment length
    
    if ($ori eq '+') {
        #print "forward\n";
        $read_end+=$read_extend;
    } else {
        #print "reverse\n";
        $read_start-=$read_extend;
        $read_start=1 if ( $read_start < 1 );
    }

    # update read length
    $read_length = $read_end-$read_start+1;
    #printf "%s\t%d\t%d\t%d\t%s\n", 
  	#$seq, $read_start, $read_end, $read_length, $ori;

    # determine bins that are covered by the read and add 
    # coverage to bin score
    
    $start_bin = sprintf("%d", ($read_start-1)/$bin_size);
    #start pos of start bin
    $start_bin_start = ($start_bin*$bin_size)+1;

    $end_bin = sprintf("%d", ($read_end-1)/$bin_size);
    #start pos of end bin
    $end_bin_start = ($end_bin*$bin_size)+1;

    #printf "%s\t%d\t%d\t%d\t%s\t%d\t%d\t%d\t%d\n", 
    #$seq, $read_start, $read_end, $read_length,	$ori,
    #$start_bin, $start_bin_start, 
    #$end_bin, $end_bin_start;

    if ($start_bin == $end_bin) { ## read start and end together in one bin

        $bin[$start_bin] += $read_length/$bin_size;
    
    } else {

        #printf "%8d\t%8d\n", $start_bin, $end_bin;

        $bin[$start_bin] += (($start_bin_start+$bin_size)-$read_start)/$bin_size;
        #print "(($start_bin_start+$bin_size)-$read_start)/$bin_size = "
        #    .(($start_bin_start+$bin_size)-$read_start)/$bin_size."\n";
        
        for (my $i=$start_bin+1; $i<$end_bin; $i++) {

            #print $i, "\n";
            $bin[$i]++;

        }

        $bin[$end_bin] += ($read_end-$end_bin_start+1)/$bin_size;
        #print "($read_end-$end_bin_start+1)/$bin_size = "
        #    .($read_end-$end_bin_start+1)/$bin_size."\n";

    }

}

&write_bins;

close FILE;
close OUT;

sub write_bins () {

    my ($bin_start, $bin_end);

    for (my $i=0; $i<=$#bin; $i++) {

        $bin_start = $i*$bin_size+1;
        $bin_end = $bin_start+$bin_size-1;
        
        if (defined $bin[$i]) {
            printf OUT "%s\t%d\t%d\t.\t%.1f\n", 
            $seq, $bin_start, $bin_end, $bin[$i];
        }

    }

    return 1;

}

1;
