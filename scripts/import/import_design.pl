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

ensembl-efg import_design.pl
  
=head1 SYNOPSIS

import_design.pl [options]

Options:

Mandatory
  -name|n          Instance name
  -format|f        Data format
  -group|g         Group name

Optional
  -pass|p          The password for the target DB, if not defined in GroupDefs.pm
  -data_root       The root data dir
  -dbname          Defines the eFG dbname if it is not standard
  -ssh             Forces use of tcp protocol for using ssh forwarded ports from remote machine(remove and just specify 127.0.0.1 as host?)
  -array_set       Flag to treat all chip designs as part of same array
  -array_name      Name of the set of array chips
  -result_set      Name to give the raw/normalised result set.
  -fasta           Fasta dump flag
  -norm|n          Normalisation method (default=vsn)
  -species|s       Species name any standard e! species alias(will be reset to dbname/latin name e.g. "homo_sapiens")
  -location        Physical location of experimental group
  -contact         Contact details for experimental group
  -debug           Debug level (1-3)
  -log_file        Defines the log file
  -help            Brief help message
  -man             Full documentation

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

=item B<-log_file|l>

Defines the log file, default = "${instance}.log"

=item B<-help>

Print a brief help message and exits.

=item B<-man>

Prints the manual page and exits.

=back

=head1 DESCRIPTION

B<This program> take several options, including an definitions file to parse and import array data into the ensembl-efg DB

=cut


#add @INC stuff here, or leave to .bashrc/.efg?

BEGIN{
	if(! defined $ENV{'EFG_DATA'}){
		if(-f "~/src/ensembl-funcgen/scripts/.efg"){
			system (". ~/src/ensembl-funcgen/scripts/.efg");
		}else{
			die ("This script requires the .efg file available from ensembl-funcgen\n".
				 "Please source it before running this script\n");
		}
	}
}
	

### MODULES ###
use Getopt::Long;
#use Carp;#For dev only? cluck not exported by default Remove this and implement in Helper
use Pod::Usage;
#POSIX? File stuff
use File::Path;
use Bio::EnsEMBL::Funcgen::Importer;
use Bio::EnsEMBL::Analysis;
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Utils::Exception qw( throw warning );

use strict;

$| = 1;#autoflush
my ($loc, $contact, $group, $pass, $dbname, $ssh, $notes_file, $format, $port, $host, $fasta, $data_dir, $name, $assembly);
my ($data_version, $help, $logname, $man, $species, $array_set, $array_file, $array_name, $vendor, $recover, $mage_tab, $user);
my $reg = "Bio::EnsEMBL::Registry";
my $output_dir = ".";


#Definitely need some sort of Defs modules for each array?

$main::_debug_level = 0;
$main::_tee = 0;

GetOptions (
			"format|f=s"   => \$format,
			"vendor|v=s"   => \$vendor,
			"name|n=s"     => \$name,
			"assembly|a=s" => \$assembly,
	    "pass|p=s"     => \$pass,
	    "port|l=s"     => \$port,
	    "host|h=s"     => \$host,
	    "user|u=s"     => \$user,
	    "ssh"          => \$ssh,
	    "dbname=s"     => \$dbname,
	    "group|g=s"    => \$group,
	    "species|s=s"      => \$species, 
	    "data_version|d=s" => \$data_version,
	    "array_set"    => \$array_set,
	    "array_name=s" => \$array_name,
			
			#only required if no standard dir format available
	    "notes_file=s" => \$notes_file,
	    "array_file=s" => \$array_file,

	    "mage_tab=s"   => \$mage_tab,#what is this for..doc please
	    "debug=i"    => \$main::_debug_level,
	    "output_dir=s" => \$output_dir,
	    "fasta"        => \$fasta,
	    "recover|r"    => \$recover,
	    "location=s"   => \$loc,
	    "contact=s"    => \$contact,
	    "tee"          => \$main::_tee,
	    "log_file=s"   => \$main::_log_file,
	    "debug_file=s" => \$main::_debug_file,
	    #debug level
	    #should have MAGE flag here? or would this be format?
	    "help|?"       => \$help,
	    "man|m"        => \$man,
	   );


pod2usage(1) if $help;
pod2usage(-exitstatus => 0, -verbose => 2) if $man;


#MANDATORY PARAMS HERE!!!!!!!!

#Build and validate all these in Experiment::new? We only need these for importing/analysing???
#need to test this before making
warn "change the implementation of Helper to set up log files etc only after we have initialised/tested Importer";

mkpath $output_dir;#log/debug files fail in Helper without this
#system("mkdir -p $output_dir -m 0755");#not working?


chmod 0755, $output_dir;


# pass as args?
if(defined $array_name){
  $logname = $array_name;
}else{
  ($logname = $array_file) =~ s/.*\///;
}

$main::_log_file = "${output_dir}/${logname}.log" if(! defined $main::_log_file);
$main::_debug_file = "${output_dir}./${logname}.dbg" if(! defined $main::_debug_file);

### SET UP IMPORTER (FUNCGENDB/DNADB) ###

my $Imp = Bio::EnsEMBL::Funcgen::Importer->new
  (
   -format      => $format,
   -vendor      => $vendor,
   -group       => $group,
   -pass        => $pass,
   -host        => $host,
   -user        => $user,
   -port        => $port,
   -ssh         => $ssh,
   -dbname      => $dbname,
   -array_set   => $array_set,
   -array_name  => $array_name,
   -assembly    => $assembly,
   -name        => $name,
   -mage_tab    => $mage_tab,#remove?
   #-data_version => $data_version,
   -data_root   => $data_dir,
   -output_dir  => $output_dir,
   -recover     => $recover,
   -dump_fasta  => $fasta,
   -species     => $species,
   -location    => $loc,
   -contact     => $contact,
  );





my $anal_a = $Imp->db->get_AnalysisAdaptor();
my $u_anal = Bio::EnsEMBL::Analysis->new(
				       -logic_name      => 'UScore',
				       -db              => 'NULL',
				       -db_version      => 'NULL',
				       -db_file         => 'NULL',
				       -program         => 'NULL',
				       -program_version => 'NULL',
				       -program_file    => 'NULL',
				       -gff_source      => 'NULL',
				       -gff_feature     => 'NULL',
				       -module          => 'NULL',
				       -module_version  => 'NULL',
				       -parameters      => 'NULL',
				       -created         => 'NULL',
				       -description     => 'Uniqueness score',
				       -display_label   => 'UScore',
				       -displayable     => 1,
				      );


my $mas_anal = Bio::EnsEMBL::Analysis->new(
					   -logic_name      => 'MASCycles',
					   -db              => 'NULL',
					   -db_version      => 'NULL',
					   -db_file         => 'NULL',
					   -program         => 'NULL',
					   -program_version => 'NULL',
					   -program_file    => 'NULL',
					   -gff_source      => 'NULL',
					   -gff_feature     => 'NULL',
					   -module          => 'NULL',
					   -module_version  => 'NULL',
					   -parameters      => 'NULL',
					   -created         => 'NULL',
					   -description     => 'Maskless array photolithographic cycles',
					   -display_label   => 'MASCycles',
					   -displayable     => 1,
					     );


my $tm_anal = Bio::EnsEMBL::Analysis->new(
					  -logic_name      => 'NimblegenTM',
					  -db              => 'NULL',
					  -db_version      => 'NULL',
					  -db_file         => 'NULL',
					  -program         => 'NULL',
					  -program_version => 'NULL',
					  -program_file    => 'NULL',
					  -gff_source      => 'NULL',
					  -gff_feature     => 'NULL',
					  -module          => 'NULL',
					  -module_version  => 'NULL',
					  -parameters      => 'NULL',
					  -created         => 'NULL',
					  -description     => 'Melting Temperature',
					  -display_label   => 'NimblegenTM',
					  -displayable     => 1,
					 );


my $tile_anal = Bio::EnsEMBL::Analysis->new(
					    -logic_name      => 'TileMap',
					    -db              => 'NULL',
					    -db_version      => 'NULL',
					    -db_file         => 'NULL',
					    -program         => 'NULL',
					    -program_version => 'NULL',
					    -program_file    => 'NULL',
					    -gff_source      => 'NULL',
					    -gff_feature     => 'NULL',
					    -module          => 'NULL',
					    -module_version  => 'NULL',
					    -parameters      => 'NULL',
					    -created         => 'NULL',
					    -description     => 'TileMap',
					    -display_label   => 'TileMap',
					    -displayable     => 1,
					 );
#Maniatis from: Bolton & McCarthy PNAS 84:1390 ?
#this checks if already stored
$anal_a->store($u_anal);
$anal_a->store($mas_anal);
$anal_a->store($tm_anal);
$anal_a->store($tile_anal);

#Validate, parse and import all array design data
$Imp->init_array_import;
$Imp->read_array_data($notes_file);
$Imp->read_probe_data($array_file);




1;
