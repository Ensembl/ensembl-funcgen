#!/usr/bin/env perl

=head1 NAME

build_profile.pl

=head1 SYNOPSIS

build_profile.pl --bin_size X --frag_size Y --files /file/path1 /file/path2 ...

=head1 DESCRIPTION

This script ...

=head1 LICENSE

  Copyright (c) 1999-2009 The European Bioinformatics Institute and
  Genome Research Limited.  All rights reserved.

  This software is distributed under a modified Apache license.
  For license details, please see

    http://www.ensembl.org/info/about/code_licence.html

=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <ensembl-dev@ebi.ac.uk>.

  Questions may also be sent to the Ensembl help desk at
  <helpdesk@ensembl.org>.


=cut

use strict;
use warnings;
use Pod::Usage;
use Getopt::Long;
use Bio::EnsEMBL::Utils::Exception qw(verbose throw warning info);

my ($file, $bin_size, $frag_length, @files);
#$dbname
my $params_msg = "Params are:\t@ARGV";

GetOptions (
			'files=s{,}'       => \@files,
            'bin_size=i'       => \$bin_size,
            'frag_length=i'    => \$frag_length,
			'help|?'           => sub { pos2usage(-exitval => 0, -message => $params_msg);},
			'man|m'            => sub { pos2usage(-exitval => 0, -message => $params_msg, verbose => 2);},
		   ) or pod2usage ( -exitval => 1,
							-message => $params_msg
						  );

if (@ARGV){
  pod2usage( -exitval =>1,
			 -message => "You have specified incomplete options. $params_msg");
}


### check options ###
throw("Must specify mandatory bin size (-bin_size).\n") if ! defined $bin_size;
throw("Must specify mandatory fragment length (-frag_length).\n") if ! defined $frag_length;


my (@bin, $start_bin, $start_bin_start, $end_bin, $end_bin_start,
    $seq, $read_start, $read_end, $read_length, $ori, $read_extend);

#Is this used?
# Get infile with features to project
if ($ENV{LSB_JOBINDEX}) {
    $file = $files[$ENV{LSB_JOBINDEX}-1];
}
else{

  if(scalar(@files) > 1){
	throw('You have specified more than one file, maybe you want to submit this to the farm using run_build_profile.sh|BuildBedProfile');
  }

  $file = $files[0];
}

if(! -e $file){
  throw("File does not exist:\t$file\nMust provide at least one file path to build a profile");
}


print "Building profile for:\t$file\n";

open(FILE, "gzip -dc $file |")
    or throw ("Can't open file $file");

my $binsize = sprintf("%03d", $bin_size);
(my $out = $file) =~ s,_reads\.bed,_profile_${binsize}.bed,;

open(OUT, "| gzip -c > $out")
    or throw ("Can't open out file $out");

while (<FILE>) {
    chomp;
    my @col = split("\t");
  
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
