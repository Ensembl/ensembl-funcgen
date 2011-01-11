#!/usr/bin/env perl

=head1 LICENSE

  Copyright (c) 1999-2011 The European Bioinformatics Institute and
  Genome Research Limited.  All rights reserved.

  This software is distributed under a modified Apache license.
  For license details, please see

    http://www.ensembl.org/info/about/code_licence.html

=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <ensembl-dev@ebi.ac.uk>.

  Questions may also be sent to the Ensembl help desk at
  <helpdesk@ensembl.org>.

=head1 NAME

sam2bed.pl

=head1 SYNOPSIS

 sam2bed.pl [ -on_farm -man -help ] -files  (file.sam[.gz])+

=head1 DESCRIPTION

This script loads converts sam(inc gzipped) format files to a bed format files.
Unampped reads are filtered as appropriate and the MAPQ score is used as the bed score.

=cut


use warnings;
use strict;
use Pod::Usage;
use Getopt::Long;
use Bio::EnsEMBL::Funcgen::Utils::EFGUtils qw(get_file_format is_gzipped open_file);


# To do
# Utilise Bio::DB::Sam
# Add filters on quality etc. See IDs code
# Add GetOpt::long support:
#     -No zip output?
#     -submit to farm

#It works using bwa output. Not sure if it is fully SAM compliant!
#This does not use the Binary format!

my ($on_farm, @files);

my @tmp_args = @ARGV;

GetOptions (
			"on_farm"    => \$on_farm,
			"files=s{,}" => \@files,
			'man'        => sub { pod2usage(-exitval => 0, -verbose => 2); },
			'help'       => sub { pod2usage(-exitval => 0, -message => "Params are:\t@tmp_args"); }
		   ) or pod2usage( -exitval => 1,
						   -message => "Params are:\t@tmp_args"
						 );


print "sam2bed.pl @tmp_args\n" if $on_farm;

#Not implemented option functions yet

if (! @files){
  pod2usage( -exitval => 1,
			 -message => "You must supply some sam filenames to convert");
}


if((scalar(@files) > 1) &&
   ! $on_farm){
  die("You must supply the -on_farm flag if you want to run on more than one file");
  
}

foreach my $file(@files){
  
  if($on_farm){
	

	my $bsub_cmd = "bsub -q normal -J sam2bed:$file -o ${file}.sam2bed.out -e ${file}.sam2bed.err ".$ENV{'EFG_SRC'}."/scripts/miscellaneous/sam2bed.pl -files $file";
	system($bsub_cmd);
	next;
  }

 
  if(&get_file_format($file) ne 'sam'){
	 pod2usage( -exitval => 1,
				-message => "Not a sam format file:\t$file");
  }
  
 
  my $gz = (&is_gzipped($file)) ? '.gz' : '';
  my $file_operator = '';
  $file_operator = "gzip -dc %s |"  if $gz;
  
  my $infile = open_file($file, $file_operator);
  my $outfile = $file;
  my $regex = ($gz) ? '\.sam\.gz' : '\.sam$';
  $outfile =~ s/${regex}/\.bed/;
  $outfile .= '.gz';
  $outfile = open_file($outfile, '| gzip -c > %s');

  print "Converting file to bed format:\t$file\n";
  my @cache;

  while(<$infile>){
	next if(/^@/);#Skip header
	chomp;
	#(Query,flag,ref,pos,mapq,cigar,mrnm,mpos,isize,seq,quak,opt)
	#my ($name, $flag, $slice_name, $pos, $mapq, $cigar, $mrnm, $mpos, $isize, $read, $quak, $opt) = split("\t");
	my ($name, $flag, $slice_name, $pos, $mapq, undef, undef, undef, undef, $read) = split("\t");
	next if $flag & 4;#Unmapped read. 
	
	#Query strand is in the 5th(index == 4) bit...  I'm assuming the reference strand never changes
	#i.e. 2*4 = 16 bitwise 16 & 16 == 16 else 0
	my $strand = ($flag & 16) ? '-' : '+';
	my (undef, undef ,$seq_region_name) = split(":", $slice_name);
	
	#Can we put something better than 100 in score position?
	#Do we really want the names in the bed file?
	#SEQ_REGION_NAME, START, END, FEATURE_NAME, SCORE STRAND
	push @cache, join("\t", ($seq_region_name, $pos, ($pos +length($read) -1), $name, $mapq, $strand));
	#print $outfile join("\t", ($seq_region_name, $pos, ($pos +length($read) -1), $name, $mapq, $strand))."\n";

	if(scalar(@cache) == 1000){
	  print $outfile join("\n", @cache)."\n";
	  @cache = ();
	}
  }

  print $outfile join("\n", @cache)."\n";

  close $infile;
  close $outfile;
}
1;
