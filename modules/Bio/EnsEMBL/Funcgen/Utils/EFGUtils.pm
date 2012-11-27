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

Bio::EnsEMBL::Funcgen::Utils::EFGUtils

=head1 DESCRIPTION

This module collates a variety of miscellaneous methods.


=head1 SYNOPSIS

  BEGIN
  {
    unshift(@INC,"/path/of/local/src/modules");
  }

  use Utils;

  &Utils::send_mail($to_address, $title, $message);

=cut


# No API/Object based methods in here

###############################################################################

package Bio::EnsEMBL::Funcgen::Utils::EFGUtils;
require Exporter;
@ISA = qw(Exporter);
@EXPORT_OK = qw(get_date species_name get_month_number species_chr_num 
				open_file median mean run_system_cmd backup_file 
				is_gzipped is_sam is_bed get_file_format strip_param_args 
				generate_slices_from_names strip_param_flags
				get_current_regulatory_input_names add_external_db);

use Bio::EnsEMBL::Utils::Exception qw( throw );
use File::Path qw (mkpath);
use File::Basename qw (dirname);
use strict;
use Time::Local;
use FileHandle;
use Carp;

sub get_date{
	my ($format, $file) = @_;

	my ($time, $sec, $min, $hour, $mday, $mon, $year, $wday, $yday, $isdst);	


	throw("File does not exist or is not a regular file:\t$file") if $file && ! -f $file;


	($sec, $min, $hour, $mday, $mon, $year, $wday, $yday, $isdst) = (defined $file) ? 
	  localtime((stat($file))[9]) : localtime();

	#print "	($sec, $min, $hour, $mday, $mon, $year, $wday, $yday, $isdst)\n";
	
	if((! defined $format && ! defined $file) || $format eq "date"){
		$time = ($year+1900)."-".$mday."-".($mon+1);	
	}
	elsif($format eq "time"){#not working!
		$time = "${hour}:${min}:${sec}";
	}
	elsif($format eq "timedate"){#
	  $time = localtime();
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
  my ($scores, $sort) = shift;

  return undef if (! @$scores);


  my ($median);
  my $count = scalar(@$scores);
  my $index = $count-1;
  #need to deal with lines with no results!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  #deal with one score fastest
  return  $scores->[0] if ($count == 1);

	 
  if($sort){
	#This is going to sort the reference here, so will affect
	#The array in the caller
	#We need to deref to avoid this
  }
  
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

#Should really extend this to detect previous file?
#Or do in caller?



sub open_file{
  my ($file, $operator, $file_permissions) = @_;
  
  $operator ||= '<';

  if ($operator !~ /%/) {
    $operator = "$operator $file";
  } else {
  	#We have some piping to do
  	$operator = sprintf($operator, $file);
  }

  #Get dir here and create if not exists
  my $dir = dirname($file);  

  my $mkpath_opts = {verbose => 1};
  $mkpath_opts->{mode} = $file_permissions if defined $file_permissions;

  mkpath($dir, $mkpath_opts) if(! -d $dir);
  my $fh = new FileHandle "$operator";


  #This does not catch incorrectly defined named pipes

  if (! defined $fh) {
    croak("Failed to open $operator");
  }


  #Have to chmod here as umask will over-ride permissions passed to FileHandle
  if (defined $file_permissions) {

    #Catch non-numeric here as chmod still returns true
    if ($file_permissions =~ /[^0-9]/) {
      croak("Failed to change $file permissions using:\t$file_permissions");
    }

    #chmod requires a literal octal number e.g. 0775 not '0775'
    #should catch numbers as strings here, but perl makes this very hard to test
    #Can't even system this as if we build the cmd line with an octal it will be converted to a decimal
    #These is still no way of testing for non-octal number or string
    #eval/sprintf will also not fail if there are non-octal digits i.e. 1999
	
    #eval will treat octal number and string as true octal number
    #else will pass non-octal string/number which we can't catch
    chmod(eval($file_permissions), $file);
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
 
  #We're not catchign error message here!

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


sub get_file_format{
  my $file = shift;

  my $format = &is_bed($file);

  if(! $format){
	$format =  &is_sam($file);

	#Add more testes here
  }
  
  
  return $format;
}

sub is_gzipped {
  my ($file, $fail_if_compressed) = @_;

  throw ("File does not exist:\t$file") if ! -e $file;

  open(FILE, "file -L $file |")
	or throw("Can't execute command 'file' on '$file'");
  my $file_info = <FILE>;
  close FILE;

  my $gzip = 0;

  if($file_info =~ /compressed data/){

	if($file_info =~ /gzip/){
	  $gzip = 1;
	}
	else{
	  throw("File is compressed but not with gzip, please unzip or gzip:\t$file_info");
	}
  }

  return $gzip;
}

#Change these to also return the gz status

sub is_sam{
  my $file = shift;

  warn "Only checking file suffix for is_sam";
  #Could check for header here altho this is not mandatory!
  #Can we use web format guessing code?

  my $gz = (&is_gzipped($file, 1)) ? '.gz' : '';

  return ($file =~ /.sam${gz}/) ? 'sam' : 0;
}

#need is bam here too!

sub is_bed {
  my $file = shift;

  #Use open_file here!
  
  if(&is_gzipped($file, 1)){
	
    open(FILE, "zcat $file 2>&1 |") or throw("Can't open file via zcat:\t$file");
  } 
  else{
    open(FILE, $file) or throw("Can't open file:\t$file");
  }

  my @line;
  #$verbose =1;


  while (<FILE>) {
	chomp;
	@line = split("\t", $_);
	last;
  }
  close FILE;
  
  if (scalar @line < 6) {
	warn "Infile '$file' does not have 6 or more columns. We expect bed format:\t".
	  "CHROM START END NAME SCORE STRAND.\n";
	return 0;
	#} elsif ($line[0] !~ m/^((chr)?[MTXYNT_\d]+)$/) {
	#    warn ("1st column must contain name of seq_region (e.g. chr1 or 1) in '$file'");
	#    return 0;
	#Commented this out for now due to HSCHR_RANDOM seqs
	#How does the webcode handle this?
  } 
  elsif ($line[1] !~ m/^\d+$/ && $line[2] =~ m/^\d+$/) {
	warn "2nd and 3rd column must contain start and end respectively in '$file'\n";
	return 0;
  } 
  elsif ($line[5] !~ m/^[+-\.]$/) {
	warn "6th column must define strand (either +, - or .) in '$file'\n";
	return 0;
  }

  return 'bed';
}


#These subs are useful for implementing
#a farm mode in a run script, where a script can
#submit itself to the farm as slice based jobs

#strip cmd line params and associated arguments from a list 
#should not be used to remove flag options i.e. no following args
#as this may cause removal of any following @ARGV;
#Can this be used on flattened args hash?

sub strip_param_args{
  my ($args, @strip_params) = @_;

  my $param_name;
  my $seen_opt = 0;

  foreach my $i(0..$#{$args}){

	if($args->[$i] =~ /^[-]+/){
	  $seen_opt = 0;#Reset seen opt if we seen a new one
	  
	  ($param_name = $args->[$i]) =~ s/^[-]+//;

	  if(grep/^${param_name}$/, @strip_params){
		$seen_opt = 1;
	  }
	}

	#$args->[$i] = '' if $args->[$i] =~ /^[-]+farm/;#Only remove current flag
	#$seen_opt = 1 if $args->[$i] =~ /^[-]+skip_slices/;
	#$seen_opt = 1 if $args->[$i] =~ /^[-]+slice/;#Don't have full param name incase we have just specified -slice
	
	$args->[$i] = '' if $seen_opt;#Remove option and args following option
  }

  return $args;
}


sub strip_param_flags{
  my ($args, @strip_params) = @_;

  my @args = @$args;

  foreach my $flag(@strip_params){
	@args = grep(!/[-]+${flag}$/, @args);
  }

  return \@args;
}

#Generates slices from names or optionally alll default top level nonref
#slice ref args can be array ref (inc empty) or undef

sub generate_slices_from_names{
  my ($slice_adaptor, $slice_names, $skip_slices, $level, $non_ref, $inc_dups, $assembly) = @_;
  my (@slices, $slice, $sr_name, $have_slice_names, $have_skip_slices);

  #Test if $assembly is old?
 
  #Validate array ref skip/slice name args
  
  if(defined $slice_names){
    if(ref($slice_names) ne 'ARRAY'){
      throw('Slice names argument must be an ARRAYREF');
    }
    
    $have_slice_names = 1 if @$slice_names;
  }

  if(defined $skip_slices){
    if(ref($skip_slices) ne 'ARRAY'){
      throw('Skip slices argument must be an ARRAYREF');
    }

    $have_skip_slices = 1 if @$skip_slices;
  }

  #Generate slices

  if($have_slice_names){

    foreach my $name(@$slice_names){
      $slice = $slice_adaptor->fetch_by_region(undef, $name, undef, undef, undef, $assembly);
      
      #WHy is this failing for hap regions?
      
      if(! $slice){
        
        #Need to eval this as it will break with incorrect formating
        
        eval { $slice = $slice_adaptor->fetch_by_name($name) };
        
        if(! $slice){
          throw("Could not fetch slice by region or name:\t".$name);
        }
      }
      
      $sr_name = $slice->seq_region_name;
      
      next if(grep/^${sr_name}$/, @$skip_slices);
      push @slices, $slice;
    }
  }
  elsif($level){
   
    #Can't guarantee what level will have an assembly version
    #Can only throw if we have set toplevel and assembly
    #as this will alway default to the current assembly
    #if($assembly && 
    #   ($level eq 'toplevel') ){
    #  throw('You cannot specify an assembly version with the toplevel coordinate system');
    #}
    #Let core API handle this?
    

    my @tmp_slices = @{$slice_adaptor->fetch_all($level, $assembly, $non_ref, $inc_dups)};
    
    if($have_skip_slices){
      
      foreach $slice(@tmp_slices){
        $sr_name = $slice->seq_region_name;
        push @slices, $slice if ! grep/^${sr_name}$/, @$skip_slices;
      }
    }
    else{
      @slices = @tmp_slices;
    }
  }
  else{
    throw('You must either pass an arrayref of slice names or specify the toplevel flag');
  }


  if(! @slices){
    throw("You have specified slice_names and skip_slices paramters which have generated no slices.\nslice_names:\t".join(' ',@$slice_names)."\nskip_slices:\t".join(' ', @$skip_slices));
  }

  return \@slices;
}


# Tracking DB methods
# Move to DBAdaptor? Can we add this as a separate package in the same module?

sub get_current_regulatory_input_names{
  my ($tdb, $efg_db, $focus) = @_;

  warn "Move get_current_regulatory_input_names to the TrackingAdaptor";

  #Validate is production?
  my $sql;

  

  if($focus){
    $focus = 'Focus';
    $sql   = 'SELECT efgdb_set_name from dataset where is_focus=true and is_current=true and species="'.$efg_db->species.'"';
  }
  else{
    $focus = 'Non-focus';
    #0 rather than false so we don't get NULLs
    $sql = 'SELECT efgdb_set_name from dataset where is_focus=0 and is_current=true and species="'.$efg_db->species.'"';
  }
 
  

  #Currently efgdb_set_name can either be data_set or feature_set name!
  #Need to standardise this

  my @prd_names = @{$tdb->db_handle->selectcol_arrayref($sql)};
  my @names;
  my @failed_sets;

  foreach my $prd_name(@prd_names){

    $sql = "SELECT name from feature_set where name like '${prd_name}%'";
    my @tmp_names =  @{$efg_db->dbc->db_handle->selectcol_arrayref($sql)};
    
    #This is causing problems with multiple feature sets with differing analyses
    
    #Do this via InputSets(using query extension?) instead of using like?
    
    #This is very hacky right now to get it to work
    #Need to standardise and review tracking db data.

    if(scalar(@tmp_names) > 1){
      
      $sql = "SELECT name from feature_set where name ='${prd_name}_ccat_histone'";
      @tmp_names =  @{$efg_db->dbc->db_handle->selectcol_arrayref($sql)};
      
      if(scalar(@tmp_names) == 1){
        push @names, $tmp_names[0];
      }else{
        push @failed_sets, $prd_name;
      }
    }
    elsif(scalar(@tmp_names) == 0){
      push @failed_sets, $prd_name;
    }
    else{
      push @names, $tmp_names[0];
    }
  }

  if(@failed_sets){
    throw("Failed to find unique $focus FeatureSets for production dataset names:\n\t".
          join("\n\t", @failed_sets)."\n");
  }

  return @names;
}

#Handy function to add an external_db entry when needed
sub add_external_db{
  my ($efg_db, $db_name,$db_release,$db_display_name) = @_;
  my $sql = "select external_db_id from external_db where db_name='$db_name' and db_release='$db_release'";
  my ($db_id) =  $efg_db->dbc->db_handle->selectrow_array($sql);
  if($db_id){ 
    warn "External DB $db_name $db_release already exists in db with db_id $db_id\n";
  } else {
    #TODO check if it there was a failure
    $efg_db->dbc->do("insert into external_db (db_name, db_release, status, dbprimary_acc_linkable, priority, db_display_name, type) values('$db_name', '$db_release', 'KNOWNXREF', '1', '5', '$db_display_name', 'MISC')");
  }

}

1;
