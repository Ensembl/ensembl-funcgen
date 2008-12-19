
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
@EXPORT_OK = qw(get_date species_name get_month_number species_chr_num open_file median mean run_system_cmd backup_file);

use Bio::EnsEMBL::Utils::Exception qw( throw );
use strict;
use Time::Local;
use FileHandle;
use Carp;


sub get_date{
	my ($format, $file) = @_;

	my ($time, $sec, $min, $hour, $mday, $mon, $year, $wday, $yday, $isdst);	


	throw("File does not exist or is not a regular file:\t$file") if ! -f $file;


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

	($species = lc($species)) =~ s/ /_/;

	my %species_chrs = (
						homo_sapiens => {(
										  'x' => 23,
										  'y' => 24,
										  'mt' => 25, 
										 )},
						
						mus_musculus => {(
										  'x'  => 20,
										  'y'  => 21,
										  'mt' => 22,
										   )},
						
						rattus_norvegicus =>  {(
												'x'  => 21,
												'y'  => 22,
												'mt' => 23,
											   )},
					   );

	die("species not defined in chromosome hash") if(! exists $species_chrs{$species});

	return (exists $species_chrs{$species}{lc($val)}) ? $species_chrs{$species}{lc($val)} : $val;
}

#Sort should always be done in the caller if required

sub median{
  my $scores = shift;

  return undef if (! @$scores);

  my ($median);
  my $count = scalar(@$scores);
  my $index = $count-1;
  #need to deal with lines with no results!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  #deal with one score fastest
  return  $scores->[0] if ($count == 1);
  
  #taken from Statistics::Descriptive
  #remeber we're dealing with size starting with 1 but indices starting at 0
  
  if ($count % 2) { #odd number of scores
    $median = $scores->[($index+1)/2];
  }
  else { #even, get mean of flanks
    $median = ($scores->[($index)/2] + $scores->[($index/2)+1] ) / 2;
  }


  return $median;
}


sub mean{
  my $scores = shift;
  
  my $total = 0;

  map $total+= $_, @$scores;
  my $mean = $total/(scalar(@$scores));

  return $mean;

}

sub open_file{
  my ($file, $operator) = @_;
	
  $operator ||= '<';

  my $fh = new FileHandle "$operator $file";
	
  if(! defined $fh){
	croak("Failed to open $operator $file");
  }

  return $fh;
}



################################################################################

=head2 run_system_cmd

 Description : Method to control the execution of the standard system() command

 ReturnType  : none

 Example     : $Helper->debug(2,"dir=$dir file=$file");

 Exceptions  : throws exception if system command returns none zero

=cut

################################################################################


#Move most of this to EFGUtils.pm
#Maintain wrapper here with throws, only warn in EFGUtils

sub run_system_cmd{
  my ($command, $no_exit) = @_;

  my $redirect = '';

  #$self->debug(3, "system($command)");
  
  # decide where the command line output should be redirected

  #This should account for redirects

  #if ($self->{_debug_level} >= 3){

  #  if (defined $self->{_debug_file}){
  #    $redirect = " >>".$self->{_debug_file}." 2>&1";
  #  }
  #  else{
  #    $redirect = "";
  #  }
  #}
  #else{
    #$redirect = " > /dev/null 2>&1";
  #}

  # execute the passed system command
  my $status = system("$command $redirect");
  my $exit_code = $status >> 8; 
 
  if ($status == -1) {	
	warn "Failed to execute: $!\n";
  }    
  elsif ($status & 127) {
	warn sprintf("Child died with signal %d, %s coredump\nError:\t$!",($status & 127),($status & 128) ? 'with' : 'without');
  }    
  elsif($status != 0) {	
	warn sprintf("Child exited with value %d\nError:\t$!\n", $exit_code); #get the true exit code
  }
 
  if ($exit_code != 0){
		  
    if (! $no_exit){
      throw("System command failed:\t$command\n");
    }
    else{
      warn("System command returned non-zero exit code:\t$command\n");
    }
  }
  
  #reverse boolean logic for perl...can't do this anymore due to tab2mage successful non-zero exit codes :/

  return $exit_code;
}


sub backup_file{
  my $file_path = shift;

  throw("Must define a file path to backup") if(! $file_path);

  if (-f $file_path) {
    #$self->log("Backing up:\t$file_path");
    system ("mv ${file_path} ${file_path}.".`date '+%T'`) == 0 || return 0;
  }

  return 1;

}


1;
