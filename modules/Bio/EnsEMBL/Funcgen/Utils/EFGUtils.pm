
=head1 NAME

Bio::EnsEMBL::Funcgen::Utils::EFGUtils

=head1 DESCRIPTION

This module collates a variety of miscellaneous methods


=head1 SYNOPSIS

  BEGIN
  {
    unshift(@INC,"/path/of/local/src/modules");
  }

  use Utils;

  &Utils::send_mail($to_address, $title, $message);


=head2 FILES


=head2 NOTES



=head2 AUTHOR(S)

Nathan Johnson njohnson@ebi.ac.uk

=cut

###############################################################################

package Bio::EnsEMBL::Funcgen::Utils::EFGUtils;
require Exporter;
@ISA = qw(Exporter);
@EXPORT_OK = qw(get_date species_name get_month_number species_chr_num open_file);

use strict;
use Time::Local;
use FileHandle;
use Carp;


sub get_date{
	my ($format, $file) = @_;

	my ($time, $sec, $min, $hour, $mday, $mon, $year, $wday, $yday, $isdst);	

	($sec, $min, $hour, $mday, $mon, $year, $wday, $yday, $isdst) = (defined $file) ? 
	  localtime((stat($file))[9]) : localtime();

	
	if((! defined $format && ! defined $file) || $format eq "date"){
		$time = ($year+1900)."-".$mday."-".($mon+1);	
	}
	elsif($format eq "time"){#not working!
		$time = "${hour}:${min}:${sec}";
	}
	else{#add mysql formats here, datetime etc...
		croak("get_date does not handle format:\t$format");
	}

	return $time;
}


#migrate this data to defs file!!??
#must contain all E! species and any other species which are used in local DB extractions
#NEED TO ADD FLY!!

sub species_name{
  my($species) = @_;
  my %species_names = (
		       "HOMO_SAPIENS", "human",
		       "MUS_MUSCULUS", "mouse",
		       "RATTUS_NORVEGICUS", "rat",
		       "CANIS_FAMILIARIS", "dog",
		       "PAN_TROGOLODYTES", "chimp",
		       "GALLUS_GALLUS", "chicken",
		       "SACCHAROMYCES_CEREVISIAE", "yeast",
		       "HUMAN",  "HOMO_SAPIENS",
		       "MOUSE", "MUS_MUSCULUS",
		       "RAT","RATTUS_NORVEGICUS",
		       "DOG", "CANIS_FAMILIARIS",
		       "CHIMP", "PAN_TROGOLODYTES",
		       "CHICKEN", "GALLUS_GALLUS",
		       "YEAST", "SACCHAROMYCES_CEREVISIAE",
		      );

  return $species_names{uc($species)};
}

sub get_month_number{
  my($mon) = @_;
  my %month_nos =(
		  "jan", "01",
		  "feb", "02",
		  "mar", "03",
		  "apr", "04",
		  "may", "05",
		  "jun", "06",
		  "jul", "07",
		  "aug", "08",
		  "sep", "09",
		  "oct", "10",
		  "nov", "11",
		  "dec", "12",
		 );
  return $month_nos{lc($mon)};
}


sub species_chr_num{
	my ($species, $val) = @_;

	($species = lc($species)) =~ s/\_/ /;

	my %species_chrs = (
						'homo sapiens' => {(
											x  => 23,
											y  => 24,
											mt => 25, 
										   )},
						
						'mus musculus' => {(
											x  => 20,
											y  => 21,
											mt => 22,
										   )},
						
						'rattus norvegicus' =>  {(
												  x  => 21,
												  y  => 22,
												  mt => 23,
												 )},
						

					   );

	die("species not defined in chromosome hash") if(! exists $species_chrs{$species});

	return (exists $species_chrs{$species}{lc($val)}) ? $species_chrs{$species}{lc($val)} : $val;


}

sub open_file{
	my ($operator, $file) = @_;
	
	my $fh = new FileHandle "$operator $file";
	
	if(! defined $fh){
		croak("Failed to open $operator $file");
	}

	return $fh;
}


1;
