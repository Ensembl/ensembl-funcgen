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

ensembl-efg create_probe_fasta.pl
  
=head1 SYNOPSIS

create_probe_fasta.pl [options]

Options:

Mandatory
  -instance|i      Instance name
  -format|f        Data format
  -group|g         Group name


=head1 OPTIONS

=over 8

=item B<-instance|i>

Mandatory:  Instance name for the data set, this is the directory where the native data files are located

=item B<-format|f>

Mandatory:  The format of the data files e.g. nimblegen

=over 8

=item B<-group|g>

Mandatory:  The name of the experimental group

=over 8

=item B<-data_root>

The root data dir containing native data and pipeline data, default = $ENV{'EFG_DATA'}

=item B<-debug>

Turns on and defines the verbosity of debugging output, 1-3, default = 0 = off

=over 8

=item B<-log_file|l>

Defines the log file, default = "${instance}.log"

=item B<-help>

Print a brief help message and exits.

=item B<-man>

Prints the manual page and exits.

=back

=head1 DESCRIPTION

B<This program> will read use the ensembl-efg API to parse array data failes and create fasta files for the remapping pipeline.

=cut

BEGIN{
	if(! defined $ENV{'EFG_DATA'}){
		if(-f "~/src/ensembl-efg/scripts/.efg"){
			system (". ~/src/ensembl-efg/scripts/.efg");
		}else{
			die ("This script requires the .efg file available from ensembl-efg\n".
				 "Please source it before running this script\n");
		}
	}
}
	
use Getopt::Long;
use Pod::Usage;
use Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Funcgen::Utils::Helper;


GetOptions (
	 	    "pass|p=s"     => \$pass,
			"port|l=s"     => \$port,
			"host|h=s"     => \$host,
			"user|u=s"     => \$user,
			"dbname=s"     => \$dbname,
			#"species|s=s"  => \$species,
			#"data_version|d=s" => \$data_version,
			"array_name=s" => \$array_name,
			"input_dir=s"  => \$input_dir,
			"output_dir=s" => \$output_dir,
			#should have MAGE flag here? or would this be format?
			"help|?"       => \$help,
			"man|m"        => \$man,
			"verbose=s"    => \$verbose,
		   );





if(! defined $input_dir || ! defined $array_name){
  die("You must supply the following args -input_dir -array_name\n");
}

print "Input Dir:\t$input_dir\nDesign Name:\t$design_name\nOutput Dir:\t$output_dir\n";


#check vars here

=pod

#SLURP PROBE POSITIONS
my $file = $input_dir."/DesignFiles/${design_name}.pos";
open(IN, $file) || die ("Cannot open file:\t$file");
my @probe_pos;
map (s/\r*\n//, @probe_pos = <IN>);
close(IN);



#REGION POSITIONS
$file = $input_dir."/DesignFiles/${design_name}.ngd";
open(IN, $file) || die ("Cannot open file:\t$file");
my ($line, $seq_id, $build, $chr, $loc, $start, $stop,  %regions);
#Need to add build id mappins in Array/VendorDefs.pm

while ($line = <IN>){
	next if $. == 1;#can we just ignore this? doing the test each time will slow it down
	$line =~ s/\r*\n//;

	#What about strand?

	($seq_id, $build, $chr, $loc) = split/\|/, $line;

	$seq_id =~ s/\s+[0-9]*$//;
	$chr =~ s/chr//;
	($start, $stop) = split/-/, $loc;


	#Do we need seq_id check here for validity?
	#overkill?
	if(exists $regions{$seq_id}){
		croak("Duplicate regions\n");
	}else{

		$chr = 23 if ($chr eq "X");
		$chr = 24 if ($chr eq "Y");
		

		$regions{$seq_id} = {
							 start => $start,
							 stop  => $stop,
							 chrom => $chr,
							 build => $build,
							};
	}

}

close(IN);


=cut


#Problem with NR probes in an array set(either by name or seq, not hashed on seq at present)
#We could simply use a hash here to see whether we've seen the probe ID before
#Or we can print all NR probes to file, sort and then remove replicates by looping and reprinting/skipping.
#Go for hash at present as work with >14million probes with less than 2GB memory



#Validate NDF files

my $array = $


$file = $input_dir."/${design_name}.ndf";#He

open(IN, $file) || die ("Cannot open file:\t$file");
$file = $output_dir."/${design_name}_probe.fasta";
open(FASTA, ">$file") || die ("Cannot open file:\t$file");

my ($xref_id, $seq, $probe_id, @features);
my $length = 50;##HARDCODED!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

while($line = <IN>){
	next if $. == 1;#or should we have format check here for files
	$line =~ s/\r*\n//;
	$loc = "";
	
	#PROBE_DESIGN_ID	CONTAINER	DESIGN_NOTE	SELECTION_CRITERIA	SEQ_ID	PROBE_SEQUENCE	MISMATCH	MATCH_INDEX	FEATURE_ID	ROW_NUM	COL_NUM	PROBE_CLASS	PROBE_ID	POSITION	DESIGN_ID	X	Y
	#Shall we use x and y here to creete a temp cel file for checking in R?
	#(undef, undef, undef, undef, $xref_id, $seq, undef, undef, 
	# undef, undef, undef, undef, $probe_id) = split/\t/, $line;
	
	my @tmp = split/\t/, $line;

	#print "@tmp\n";

	#next if $tmp[11] ne "experimental";
	next if $tmp[1] !~ /BLOCK1/;
	#5 & 12

	###PROBE FEATURES
	#Put checks in here for build?
	#grep for features, need to handle more than one mapping

#	@features = grep (/\s+$probe_id\s+/, @probe_pos);

#	croak("Multiple probe_features: feature code needs altering") if (scalar(@features > 1));

	#Need to handle controls/randoms here
	#won't have features but will have results!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	#The format of the pos file looks like it should have all the data required, but 
	# chromsome is missing, first undef :(

#	if(@features){
		
#		foreach $feature(@features){
#			#$feature =~ s/\r*\n//;
#			
#			#SEQ_ID	CHROMOSOME	PROBE_ID	POSITION	COUNT
#			($seq_id, undef, undef, $start, undef) = split/\t/, $feature;
#
#			if(exists $regions{$seq_id}){			
#				$loc .= $regions{$seq_id}{'chrom'}.":${start}-".($start + $length).";";
#			}
#			else{ croak("No regions defined for $seq_id"); }
#		}
#	}

	#else{#CONTROL/RANDOM
	#	#Enter placeholder features to enable result entries
	#	print PROBE_FEATURE "\t0\t0\t0\t0\t0\t${pid}\t0\t0\tNA\n";
	##	push @{$feature_map{$probe_id}}, $fid;	
	#	$fid++;
	#}
	

	

	
	#filter controls/randoms?  Or would it be sensible to see where they map
	#Do we need to wrap seq here?
	#print FASTA ">${probe_id}\t$xref_id\t$loc\n$seq\n";

	print FASTA ">probe:${design_name}:$tmp[12]\n$tmp[5]\n";

	$loc = "";
}

print "Processed ".($.- 1)." probes\n";

close(IN);
close(FASTA);
