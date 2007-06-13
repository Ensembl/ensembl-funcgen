#!/software/bin/perl -w


####!/opt/local/bin/perl -w


=head1 NAME

import_type.pl
  
=head1 SYNOPSIS

import_array_from_fasta.pl [options]

The script will import a new CellType, FeatureType or Analysis.


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
	if (-f "~/src/ensembl-functgenomics/scripts/.efg") {
	  system (". ~/src/ensembl-functgenomics/scripts/.efg");
	} else {
	  die ("This script requires the .efg file available from ensembl-functgenomics\n".
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
my ($pass, $dbname, $help, $man, $array_name, $line, $label);
my ($clobber, $type, $desc, $file, $class, $logic_name, $name);
my $user = "ensadmin";
my $host = 'ens-genomics1';
my $port = '3306';



#this should import using the API
#taking array name vendor args to populate the appropriate array/arary_chip records
#or parse them from the info line?
#currently just generates and imports flat file

#should also build cache and generate nr file?
#this depends on id/name field refering to unique seq
#same name can't refer to more than one seq

GetOptions (
			"file|f=s"        => \$file,
			"pass|p=s"        => \$pass,
			"port=s"          => \$port,
			"host|h=s"        => \$host,
			"user|u=s"        => \$user,
			"dbname|d=s"      => \$dbname,
			"help|?"          => \$help,
			"man|m"           => \$man,
			"type|t=s"        => \$type,
			"class=s"         => \$class,
			"display_label=s" => \$label,
			"name=s"          => \$name,
			"logic_name=s"    => \$logic_name,
			'clobber'         => \$clobber,#update old entries?
			"description=s"   => \$desc,
		   );


pod2usage(1) if $help;
pod2usage(-exitstatus => 0, -verbose => 2) if $man;

#This should work for any object so long as we set u the config correctly

my %type_config = (
				   'FeatureType' => {(
									  fetch_method     => 'fetch_by_name',
									  fetch_args       => [$name],
									  mandatory_params => {(
															-name        => $name,
														   )},
									  optional_params  => {(
															-class       => $class,
															-description => $desc,
															
														   )},
									 )},

				   'CellType' => {(
								   fetch_method     => 'fetch_by_name',
								   fetch_args       => [$name],
								   mandatory_params => {(
														 -name          => $name,
														)},
								   optional_params  => {(
														 -display_label => $label,
														 -description   => $desc,
														)},

								  )},
				   
				   'Analysis' => {(
								   fetch_method => 'fetch_by_logic_name',
								   fetch_args   => [$logic_name],
								   mandatory_params => {()},
								   optional_params => {()},
								  )},
				  
				   
				  );

#generic mandatory params
if(!(exists $type_config{$type} && $dbname && $pass)){
  throw('Mandatory parameters not met, more here please');
}


#now do type specific checking

my $db = Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor->new(
													  -dbname => $dbname,
													  -port   => $port,
													  -pass   => $pass,
													  -host   => $host,
													  -user   => $user,
													 );

my ($adaptor, $method, $obj_class);
my $fetch_method = $type_config{$type}->{'fetch_method'};

if($file){

  #parse file here
  #parse headers to match to params and call relevant sub


  throw('File import not yet implemented');


  #would need to clean hash values here
  my $in = open_file($file);
  my ($line);

  while ($line = <$in>){
	next if $line =~ /^#/;

	chomp $line;

	if($line =~ /^>/){#found new header
	  
	  #set type and adaptor
	  #validate header vs. mandatory keys

	}else{

	  #clean param values here
	  #call import_type
	  

	}

  }



}else{
  $obj_class = 'Bio::EnsEMBL::Funcgen::'.$type;
  $method = 'get_'.$type.'Adaptor';
  $adaptor = $db->$method();
  &import_type();
}


sub import_type{

  #check mandatorys here

  foreach my $man_param(keys %{$type_config{$type}->{'mandatory_params'}}){
	throw ("$man_param not defined") if ! defined $type_config{$type}->{'mandatory_params'}->{$man_param};
  }


  #test if already present
  my $obj = $adaptor->$fetch_method(@{$type_config{$type}->{'fetch_args'}});

  if(defined $obj){
	warn("Found pre-existing $type object:\t".join("\t", @{$type_config{$type}->{'fetch_args'}}).
		 "\nClobber/Update not yet implementing, skipping import\n");
	
  }
  else{
	$obj = new $obj_class(%{$type_config{$type}->{'mandatory_params'}}, 
						   %{$type_config{$type}->{'optional_params'}});

	$adaptor->store($obj);
  }

  return;
}
