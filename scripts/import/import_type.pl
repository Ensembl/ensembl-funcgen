#!/usr/bin/env perl

=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2018] EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

  Questions may also be sent to the Ensembl help desk at
  <http://www.ensembl.org/Help/Contact>.

=head1 NAME

import_type.pl
  
=head1 SYNOPSIS

import_type.pl [options]

The script will import a new CellType, FeatureType or Analysis.


=head1 OPTIONS

=over 8

=item B<-name|n>

Mandatory:  

=item B<-type|t>

Mandatory:  The type you are trying to import e.g. FeatureType or CellType


=item B<-species|s>

Species name for the array.

=item B<-help>

Print a brief help message and exits.

=item B<-man>

Prints the manual page and exits.

=back

=head1 DESCRIPTION

B<This program> imports FeatureType or CellType entries either from paramters or from a flat file

=cut

#To do
#1 POD Docs
#2 Helper logs
#3 Remove/Implement analysis?



BEGIN{
  if (! defined $ENV{'EFG_DATA'}) {
	if (-f "~/src/ensembl-funcgen/scripts/.efg") {
	  system (". ~/src/ensembl-funcgen/scripts/.efg");
	} else {
	  die ("This script requires the .efg file available from ensembl-funcgen\n".
		   "Please source it before running this script\n");
	}
  }
}
	

#use Bio::EnsEMBL::Root; #Only used for rearrange see pdocs
#Roll own Root object to handle debug levels, logging, dumps etc.

### MODULES ###
use Getopt::Long;
#use Carp;#For dev only? cluck not exported by default Remove this and implement in Helper
use Pod::Usage;
#POSIX? File stuff
use File::Path;
#use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Utils::Exception qw( throw warning );
use Bio::EnsEMBL::Funcgen::Utils::EFGUtils qw (open_file run_system_cmd backup_file);
use Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Funcgen::FeatureType;
use Bio::EnsEMBL::Funcgen::CellType;
use Bio::EnsEMBL::Analysis;
use Data::Dumper;
use strict;

$| = 1;							#autoflush
my ($pass, $dbname, $array_name, $line, $label, $dnadb_user, $dnadb_port,$dnadb_name);
my ($clobber, $type, $desc, $file, $class, $logic_name, $name, $dnadb_host);
my ($anal_db, $db_version, $db_file, $program, $program_version, $program_file);
my ($gff_source, $gff_feature, $module, $module_version, $parameters, $created);
my ($displayable, $web_data, $species, $gender);

#Need to change these to match EFG_USER EFG_HOST EFG_PORT
#And then test
my $user = "ensadmin";
my $host = 'ens-genomics1';
my $port = '3306';


#should also build cache and generate nr file?
#this depends on id/name field refering to unique seq
#same name can't refer to more than one seq
my @tmp_args = @ARGV;

GetOptions (
			#general params
			"file|f=s"        => \$file,
			"pass|p=s"        => \$pass,
			"port=s"          => \$port,
			"host|h=s"        => \$host,
			"dnadb_host=s"    => \$dnadb_host,
			"dnadb_user=s"    => \$dnadb_user,
			"dnadb_port=s"    => \$dnadb_port,
			"dnadb_name=s"    => \$dnadb_name,
			"user|u=s"        => \$user,
			"dbname|d=s"      => \$dbname,
			"species=s"       => \$species,
			"help|?"          => sub { pos2usage(-exitval => 0, -message => "Params are:\t@tmp_args"); },
			"man|m"           => sub { pos2usage(-exitval => 0, -verbose => 2, -message => "Params are:\t@tmp_args"); },
			"type|t=s"        => \$type,
			'clobber'         => \$clobber,#update old entries?
			#Cell/Feature params
			"class=s"         => \$class,#FeatureType only
			"display_label=s" => \$label,
			"name=s"          => \$name,
			"description=s"   => \$desc,
			"gender=s"        => \$gender,
			#analysis opts
			"logic_name=s"    => \$logic_name,
			"db=s"            => \$anal_db,
			"db_version=s"    => \$db_version,
			"db_file=s"       => \$db_file,
			"program=s"       => \$program,
			"program_version=s" => \$program_version,
			"program_file=s"    => \$program_file,
			"gff_source=s"      => \$gff_source,
			"gff_feature=s"     => \$gff_feature,
			"module=s"          => \$module,
			"module_version=s"  => \$module_version,
			"parameters=s"      => \$parameters,
			"created=s"         => \$created,
			"displayable=s"     => \$displayable,
			"web_data=s"        => \$web_data,
		   ) or pod2usage(
						 -exitval => 1,
						 -message => "Params are:\t@tmp_args");


#This should work for any object so long as we set u the config correctly

my %type_config = (
				   'FeatureType' => {(
									  class            => 'Bio::EnsEMBL::Funcgen::FeatureType',
									  fetch_method     => 'fetch_by_name',
									  fetch_arg       => '-name',
									  mandatory_params => {(
															-name        => $name,
														   )},
									  optional_params  => {(
															-class       => $class,
															-description => $desc,
															
														   )},
									 )},

				   'CellType' => {(
								   class            => 'Bio::EnsEMBL::Funcgen::CellType',
								   fetch_method     => 'fetch_by_name',
								   fetch_arg       => '-name',
								   mandatory_params => {(
														 -name          => $name,,
														)},
								   optional_params  => {(
														 -gender        => $gender,
														 -display_label => $label,
														 -description   => $desc,
														)},

								  )},
				   
				   'Analysis' => {(
								   class            => 'Bio::EnsEMBL::Analysis',
								   fetch_method => 'fetch_by_logic_name',
								   fetch_arg   => '-logic_name',
								   
								   #DB
								   #DB_VERSION
								   #DB_FILE
								   #PROGRAM
								   #PROGRAM_VERSION
								   #PROGRAM_FILE
								   #GFF_SOURCE
								   #GFF_FEATURE
								   #MODULE
								   #MODULE_VERSION
								   #PARAMETERS
								   #CREATED
								   #LOGIC_NAME
								   #DESCRIPTION
								   #DISPLAY_LABEL
								   #DISPLAYABLE
								   #WEB_DATA


								   #this is assumed mandatory params as they are not forced in Analysis->new

								   #select a.logic_name, a.db, a.db_version, a.db_file, a.program, a.program_version, a.program_file, a.parameters, a.module, a.module_version, a.created, a.gff_source, a.gff_feature, ad.description, ad.display_label, ad.displayable, ad.web_data from analysis a left join analysis_description ad on a.analysis_id=ad.analysis_id where logic_name not like "%Probe%Align" and logic_name !='probe2transcript';


								   mandatory_params => {(
														 -logic_name => $logic_name,
														)},
								   optional_params => {(
														-db => $anal_db,
														-db_version => $db_version,
														-db_file => $db_file,
														-program => $program,
														-program_version => $program_version,
														-program_file => $program_file,
														-gff_source => $gff_source,
														-gff_feature => $gff_feature,
														-module => $module,
														-module_version => $module_version,
														-parameters => $parameters,
														-created => $created,
														-description => $desc, #DESCRIPTION
														-display_label => $label,#DISPLAY_LABEL
														-displayable => $displayable,
														-web_data => $web_data,
													   )},
								  )},
				  
				   
				  );

#generic mandatory params
if(!(exists $type_config{$type} && $dbname && $pass && $species)){
  throw("Mandatory parameters not met -dbname($dbname) -pass($pass) -species($species) or $type config is not yet accomodated");
}

if($file && ($name || $logic_name)){
  throw('Cannot specify -file and -name|logic_name, use one or other');
}


#now do type specific checking

my $db = Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor->new(
													  -dbname  => $dbname,
													  -port    => $port,
													  -pass    => $pass,
													  -host    => $host,
													  -user    => $user,
													  -dnadb_host => $dnadb_host,
													  -dnadb_port => $dnadb_port,
													  -dnadb_user => $dnadb_user,
													  -dnadb_name => $dnadb_name,
													  -species => $species,
													 );
#test db connections
$db->dbc->db_handle;
$db->dnadb->dbc->db_handle;

my $info;

if(defined $file){
  $info = "s from file:\t$file";
}
elsif($name){
  $info = ":\t$name";
}
elsif($logic_name){
  $info = ":\t$logic_name";
}

print "Importing ${type}${info}\n";

my $fetch_method = $type_config{$type}->{'fetch_method'};
my $obj_class = $type_config{$type}->{'class'};
my $method = 'get_'.$type.'Adaptor';
my $adaptor = $db->$method();
my ($field, @values);

if($file){

  #parse file here
  #parse headers to match to params and call relevant sub


  my @fields = keys %{$type_config{$type}{mandatory_params}};
  push @fields, keys %{$type_config{$type}{optional_params}};

  map $_=~ s/^-//, @fields;

  #Could set header hash here? using helper?
  #Or should this be in EFGUtils?

  #would need to clean hash values here
  my $in = open_file($file);

  $line = <$in>;
  chomp $line;
  my @header = split /\t/, $line;

  #mysql -hens-genomics1 -uensro -e "select name, class, description from feature_type where class in('Histone', 'Regulatory Feature', 'Open Chromatin', 'Insulator')" homo_sapiens_funcgen_55_37

  my $hposns = set_header_hash(\@header, \@fields);


  while ($line = <$in>){
	next if $line =~ /^#/;
	next if $line =~ /^$/;

	chomp $line;

	@values = split /\t/, $line;

	#warn "values are ".join(' x ', @values);

	
	#This will clean all the old values
	foreach my $param(keys %{$type_config{$type}{mandatory_params}}){
	  ($field = $param) =~ s/^-//;
	  $type_config{$type}{mandatory_params}{$param} = ($values[$hposns->{$field}] eq 'NULL') ? undef : $values[$hposns->{$field}];
	}
	

	foreach my $param(keys %{$type_config{$type}{optional_params}}){
	  ($field = $param) =~ s/^-//;


	  if(defined $values[$hposns->{$field}]){

		$type_config{$type}{optional_params}{$param} = ($values[$hposns->{$field}] eq 'NULL') ? undef : $values[$hposns->{$field}];
	  }
	  else{

		if($type ne 'Analysis'){
		  warn "WARNING:\t$field is not defined for ".$type_config{$type}{mandatory_params}{-name}."\n";
		}
		else{
		  warn "WANRING:\t$field is not defined for ".$type_config{$type}{mandatory_params}{-logic_name}."\n";
		}
	  }
	}
		
	&import_type;
  }
}
else{
  #Values already set
  &import_type;
}


sub import_type{

  #check mandatorys here

  foreach my $man_param(keys %{$type_config{$type}{'mandatory_params'}}){
	throw ("$man_param not defined") if ! defined $type_config{$type}->{'mandatory_params'}->{$man_param};
  }


  #test if already present
  my $obj = $adaptor->$fetch_method($type_config{$type}{mandatory_params}->{$type_config{$type}->{'fetch_arg'}});
  
  if(defined $obj){
	warn("Found pre-existing $type object:\t".$type_config{$type}{mandatory_params}->{$type_config{$type}->{'fetch_arg'}}.
		 "\nClobber/Update not yet implementing, skipping import\n");
	
  }
  else{
	
	#if(exists $type_config{$type}->{'mandatory_params'}{'-logic_name'} &&
	 #  ($type_config{$type}->{'mandatory_params'}{'-logic_name'} eq 'Nessie_NG_STD_2')){

	
	 # warn Data::Dumper::Dumper($type_config{$type}->{'mandatory_params'});
	 # warn Data::Dumper::Dumper($type_config{$type}->{'optional_params'});
	  

	#}
	$obj = new $obj_class(%{$type_config{$type}->{'mandatory_params'}}, 
						   %{$type_config{$type}->{'optional_params'}});

	print "Storing $type ".$type_config{$type}{mandatory_params}->{$type_config{$type}->{'fetch_arg'}}."\n";

	$adaptor->store($obj);
  }

  return;
}


#Should use helper for this and add logging too

sub set_header_hash{
  my ($header_ref, $fields) = @_;
	
  my %hpos;

  for my $x(0..$#{$header_ref}){
	$hpos{$header_ref->[$x]} = $x;
  }	


  if($fields){

    foreach my $field(@$fields){
	  
      if(! exists $hpos{$field}){
		throw("Header does not contain mandatory field:\t${field}");
      }
    }
  }
  
  return \%hpos;
}
