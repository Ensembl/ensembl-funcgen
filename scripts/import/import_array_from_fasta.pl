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

import_array_from_fasta.pl
  
=head1 SYNOPSIS

import_array_from_fasta.pl [options] file

This is a short cut import with little validation. Expects a pre-process non-redundant fasta file, 
otherwise any probe mapping done based on the output of this script will generate duplicate probes 
and features.

Options:

Mandatory


Optional


=head1 OPTIONS

=over 8

=item B<-name|n>

Mandatory:  Instance name for the data set, this is the directory where the native data files are located

=item B<-format|f>

Mandatory:  The format of the data files e.g. nimblegen

=over 8

=item B<-group|g>

Mandatory:  The name of the experimental group

=over 8

=item B<-data_root>

The root data dir containing native data and pipeline data, default = $ENV{'EFG_DATA'}

=over 8

=item B<-fasta>

Flag to turn on dumping of all probe_features in fasta format for the remapping pipeline

=item B<-norm>

Normalisation method, deafult is the Bioconductor vsn package which performs generalised log ratio transformations

=item B<-species|s>

Species name for the array.

=item B<-debug>

Turns on and defines the verbosity of debugging output, 1-3, default = 0 = off

=over 8

=item B<-help>

Print a brief help message and exits.

=item B<-man>

Prints the manual page and exits.

=back

=head1 DESCRIPTION

B<This program> takes a input redundant probe name fasta file and generates an NR probe dbID fasta file.

=cut


#add @INC stuff here, or leave to .bashrc/.efg?

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
use Bio::EnsEMBL::Funcgen::Array;
use Bio::EnsEMBL::Funcgen::ArrayChip;
use Bio::EnsEMBL::Funcgen::Probe;
use strict;

$| = 1;							#autoflush
my ($pass, $dbname, $help, $man, $array_name, $line);
my ($format, $clobber, $vendor, $desc, $file);
#my $reg = "Bio::EnsEMBL::Registry";
my $data_dir = $ENV{'EFG_DATA'};
my $user = "ensadmin";
my $host = 'ens-genomics1';
my $port = '3306';
my $out_dir = '.';


#this should import using the API
#taking array name vendor args to populate the appropriate array/arary_chip records
#or parse them from the info line?
#currently just generates and imports flat file

#should also build cache and generate nr file?
#this depends on id/name field refering to unique seq
#same name can't refer to more than one seq

GetOptions (
			"file|f=s"       => \$file,
			"pass|p=s"     => \$pass,
			"port=s"     => \$port,
			"host|h=s"     => \$host,
			"user|u=s"     => \$user,
			"dbname|d=s"   => \$dbname,
			"help|?"       => \$help,
			"man|m"        => \$man,
			"array|a=s"    => \$array_name,
			"vendor|v=s" => \$vendor,
			'clobber'          => \$clobber,
			"format=s"       => \$format,
			"description=s" => \$desc,
			"outdir|o=s" => \$out_dir,
		   );


pod2usage(1) if $help;
pod2usage(-exitstatus => 0, -verbose => 2) if $man;


#do mandatory params and checking here


if(!($array_name && $dbname && $vendor && $pass)){
  throw('Mandatory parameters not met, more here please');
}


my $db = Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor->new(
													  -dbname => $dbname,
													  -port   => $port,
													  -pass   => $pass,
													  -host   => $host,
													  -user   => $user,
													 );




my $array_a = $db->get_ArrayAdaptor();
my $array_chip_a = $db->get_ArrayChipAdaptor();
my $probe_a = $db->get_ProbeAdaptor();

#we cant use the same methodology too store the array/arraychips as they are not in the 
#original format i.e. there will be no status control on array chip level

my $array = $array_a->fetch_by_name_vendor($array_name, $vendor);

if($array){
  
  if(! $clobber){
	
	throw("Array already exists, specify -clobber only if you are aboslutely sure you want to'");
  }else{
	warn "clobber not yer implementd";
  }
}

$array = Bio::EnsEMBL::Funcgen::Array->new(
										   -NAME        => $array_name,
										   -FORMAT      => uc($format),
										   -VENDOR      => uc($vendor),
										   -TYPE        => 'OLIGO',
										   -DESCRIPTION => $desc,
										  );


($array) = @{$array_a->store($array)};

my $array_chip = Bio::EnsEMBL::Funcgen::ArrayChip->new(
													   -ARRAY_ID    => $array->dbID(),
													   -NAME        => $array_name,
													   -DESIGN_ID   => $array_name,
													  );

($array_chip) = @{$array_chip_a->store($array_chip)};

#don't really need to do this?
$array->add_ArrayChip($array_chip);



#set up the nr fasta file
my $nr_fasta = open_file($out_dir."/${array_name}_nr.fasta", '>');

#open in file
my $in = open_file($file);

#do the biz
my ($aname, $pid);
my $ac_id = $array_chip->dbID();

while($line = <$in>){
  chomp $line;

  if($line =~ /^>/){
	(undef, $aname, $pid) = split/\:/, $line;
  }
  else{#found seq line
	
	my $probe = Bio::EnsEMBL::Funcgen::Probe->new(
											   -NAME          => $aname.':'.$pid,
											   -LENGTH        => length($line),
											   -ARRAY         => $array,
											   -ARRAY_CHIP_ID => $ac_id,
											   -CLASS         => 'EXPERIMENTAL',
											  );

	($probe) = @{$probe_a->store($probe)};


	print $nr_fasta '>'.$probe->dbID()."\n${line}\n";

  }
}

close($in);
close($nr_fasta);

