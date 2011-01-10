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

generate_DAS_config.pl

=head1 SYNOPSIS

 GenerateDASConfig|generate_DAS_config.pl [options]

 e.g

 Using the efg environment function:

 GenerateDASConfig -dbhost $DB_HOST -dbport $DB_PORT -dbuser $DB_RO_USER -dbname homo_sapiens_funcgen_47_36i -species homo_sapiens -das_config $EFG_SRC/config/DAS -das_name efg -das_host DAShost -das_port  9000

=head1 DESCRIPTION

This script writes DAS configuration for all DAS_DISPLAYABLE sets or Hydra source tables from a given DB. Most of the parameter use default values which can be set in the efg environment.

=head1 OPTIONS

 Mandatory: If running without -only_headers (default).
    DB parameters
    --dbhost  DB host name, default is $DB_HOST
    --dbport  DB port, default is $DB_PORT
    --dbuser  DB user name, default is $DB_RO_USER
    --dbname  Name of DB, default is $DB_NAME
    --species Latin name of species as used in DB name, e.g. homo_sapiens. Default is $SPECIES

 Mandatory: If running without -no_header (default).
    DAS parameters
    --das_config  default is $EFG_DAS_CONFIG
    --das_name    Das instance name
    --das_host    default is $EFG_DAS_HOST
    --das_port    default is $EFG_DAS_PORT
    
 Optional:

    DB parameters
    --dbpass      default is $DB_PASS

    DAS paramters
    --hydra_instance  Unique name of source hydra instance, this is to differentiate between different host DBs.
                      Default is 'eFG_SOURCETYPE:DBNAME@DBHOST:DBPORT
                      e.g eFG_features:homo_sapiens_funcgen_56_37a@mydbhost:3306
                      WARNING: This will be visible in the URL, so set this if you are concerned about privacy.
                      e.g. human_sources_1, human_sources_2 etc.
                      NOTE: StartDASServer will not catch any duplication of hydra_instance names between host
                      So any duplication will result in sources being overwritten.
    --severroot       ProServer root code directory (default = $SRC/Bio-Das-ProServer)
    --maintainer      email of DAS server admin
    --maxclients      Maximum DAS clients, default is 20
    --region          Region for feature capability test (default = 17:35640000,35650000)
    --assemblies      Names of coordinate system e.g. GRCh37 (default = default CoordSystem in dnadb)

    DNA DB parameters, default is to use ensembldb.ensembl.org
    --dnadb_host  Core DB host name, default is $DNADB_HOST
    --dnadb_port  Core DB port, default is $DNADB_PORT
    --dnadb_user  Core DB user name, default is $DNADB_USER
    --dnadb_name  Name of core DB, default is $DNADB_NAME
    --dnadb_pass  Core DB password, default is $DNADB_PASS

    Run modes
    --no_headers   Only prints source config
    --only_headers Only prints DAS server config
    --help         Prints this documentation and exits
    --not_hydra    Generates standard sources, default is dynamic Hydra sources
    --source_types List of Hydra source types to generate, default is: feature, result, bed
    

    Individual Set handling, for use with --not_hydra
    --set_type        Set type e.g. result, feature, regualtory
    --set_name        Name of set
    #-set_colour      Colour of track for given set e.g. contigblue1, contigblue2, red3
    #-set_plot        Plot type for given set e.g. hist or tiling

    HTML/Feature link params, used to generate page with attachement links
    --link_gene       Name of gene e.g. STAT1
    or
    --link_region     Loci e.g. 2:191541121-191588181

    #Add more here for default colours?

 Other:
   --help Prints a short help message and exits
   --man  Prints the man page and exits

=cut


#TO DO
#This needs to use the DAS DISPLAYABLE status to auto set up for a given list of hosts
#Or use default available config host file to poll the existing DB
#Or add a given dbhost using either the 
#Then we cat the individual host files with the general section to make the complete ini file
#What about losing sources that have been set up explicitly without DAS_DISPLAYABLE?
#This script should set DAS_DISPLAYABLE!
#Add option to clear away all previous source config
#We need func to list sources given host, cell/feature type, experiment name?
#Validate species vs registry.
#Implement location link! Currently hardcoded to some human loci


use strict;
use warnings;
use Pod::Usage;
use Getopt::Long;

#use Bio::EnsEMBL::Utils::Exception qw(verbose throw warning info);
#use Bio::EnsEMBL::Analysis::Tools::Logger qw(logger_verbosity logger_info);
#Use Helper instead of Logger?

use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor;

my ($set_type, $dnadb, $set_name, $das_name, @source_types, $hydra_instance);
my (@adaptor_names, $no_headers, $headers_only, $link_region, $cs_version, @cs_versions);
my ($set_colour, $set_plot, $set_display_name, $location);


my $maintainer = 'your.email@ddress.co.uk';
my $coords = '';

#ENV DEFAULTS

#These should be set in efg.env
my $serverroot  = $ENV{'EFG_DAS_HOME'};
my $das_config  = $ENV{'EFG_DAS_CONFIG'};
my $das_port    = $ENV{'EFG_DAS_PORT'};
my $das_host    = $ENV{'EFG_DAS_HOST'};

#These are not that useful as we will only be doing this once per DB for hydra
#So not much point in have a specific env file?

#Change these to use EFG_READ_USER?
#Grr should these have EFG prefix?
my $dbname      = $ENV{'DB_NAME'};
my $dbhost      = $ENV{'DB_HOST'};
my $dbuser      = $ENV{'DB_RO_USER'};
my $dbport      = $ENV{'DB_PORT'};
my $dbpass      = $ENV{'DB_PASS'};
my $species     = $ENV{'SPECIES'};
my $dnadb_name  = $ENV{'DNADB_NAME'};
my $dnadb_host  = $ENV{'DNADB_HOST'};
my $dnadb_user  = $ENV{'DNADB_USER'};
my $dnadb_port  = $ENV{'DNADB_PORT'};
my $dnadb_pass  = $ENV{'DNADB_PASS'};



my $styleshome = $serverroot.'/stylesheets';
#coordshome?

my %types = (
			 'feature' => {(
								adaptor   => 'ensembl_funcgen_set',
								hydra     => 'ensembl_funcgen',  
								transport => 'ensembl_funcgen',
								#Is this relevant for feature_sets?
								#As always use DAS_DISPLAYABLE feature sets?
								basename  => '',
							  )},
			 'result'  => {(
								adaptor   => 'ensembl_funcgen_set',
								hydra     => 'ensembl_funcgen',  
								transport => 'ensembl_funcgen',
								basename  => '',
							  )},
			 
			 'bed'        => {(
							   adaptor   => 'ensembl_funcgen_reads',
							   hydra     => 'ensembl_funcgen',  
							   transport => 'ensembl_funcgen',
							   basename  => "basename\t= bed\n",
							   #This is used to fetch tables like "$basename"
							   #So has to be mysql compliant wildcards?
							   #Actually, no.
							   #Hydra dbi does not handle wild cards well 
							   #and corrupts table name
							   #handle in SourceAdaptor/efg_reads.pm
							 )},
			 
			 #should we separate reads and profile?
			);

my $ini_password= '';
my $max_clients = 20;
#my $das_name = 'efg';
my $link_gene = 'STAT1';
my $features_region = '17:35640000,35650000';
my $not_hydra = 0;
my $pod_params = "Params are:\t".join(' ', @ARGV);


my %plots = (
			 hist   => '+score=h',
			 cgrad  => '+score=c+fg_grades=100+fg_data=n+fg_max=100',
			 tiling => '+score=t',
			);

my %default_colours = (
					   result  => 'red',
					   feature => 'blue',
					  );

my %display_params = (
				  result  => '+score=s+fg_data=o+',
				  feature => '',
				 );

#Need to write this to config file and use this to generate das xsl
#We should always use config unless we have a value which isn't defined in config
#Then use default or individual set params if the set name matches
#Then write back to config file so we don't have to do this again
#Hence can sey up each set individually without over writing other set config
#This does mean regenerating the file each time, but thi sis little over head.

my %set_config = (
           ctcf_ren_BR1_TR1_ =>    { color => 'contigblue2', name => 'IMR90_CTCF' },
           Nessie_NG_STD_2_ctcf_ren_BR1 => { color => 'contigblue2', name => 'IMR90_CTCF_Nessie'},
           ctcf_ren_IMPORT_TileMap => { color => 'contigblue2', name => 'ctcf_ren_IMPORT_TileMap'},
           ctcf_ren_IMPORT_Chipotle => { color => 'contigblue2', name => 'ctcf_ren_IMPORT_Chipotle'},

           GM06990_DNASE_IMPORT => { color => 'contigblue2', name => 'GM06990_DNASE' },
           CD4_DNASE_IMPORT =>     { color => 'contigblue2', name => 'CD4_DNASE' },
           
           CD4_CTCF => { color => 'contigblue1', name => 'CD4_CTCF', plot => 'hist' },

           CD4_H2AZ => { color => 'contigblue1', name => 'CD4_H2AZ', plot => 'hist' },

           CD4_H2BK5me1 => { color => 'contigblue1', name => 'CD4_H2BK5me1', plot => 'hist' },

           CD4_H3K27me1 => { color => 'contigblue1', name => 'CD4_H3K27me1', plot => 'hist' },
           CD4_H3K27me2 => { color => 'contigblue1', name => 'CD4_H3K27me2', plot => 'hist' },
           CD4_H3K27me3 => { color => 'contigblue1', name => 'CD4_H3K27me3', plot => 'hist' },

           CD4_H3K36me1 => { color => 'contigblue1', name => 'CD4_H3K36me1', plot => 'hist' },
           CD4_H3K36me3 => { color => 'contigblue1', name => 'CD4_H3K36me3', plot => 'hist' },

           CD4_H3K4me1 => { color => 'contigblue1', name => 'CD4_H3K4me1', plot => 'hist' },
           CD4_H3K4me2 => { color => 'contigblue1', name => 'CD4_H3K4me2', plot => 'hist' },
           CD4_H3K4me3 => { color => 'contigblue1', name => 'CD4_H3K4me3', plot => 'hist' },

           CD4_H3K79me1 => { color => 'contigblue1', name => 'CD4_H3K79me1', plot => 'hist' },
           CD4_H3K79me2 => { color => 'contigblue1', name => 'CD4_H3K79me2', plot => 'hist' },
           CD4_H3K79me3 => { color => 'contigblue1', name => 'CD4_H3K79me3', plot => 'hist' },

           CD4_H3K9me1 => { color => 'contigblue1', name => 'CD4_H3K9me1', plot => 'hist' },
           CD4_H3K9me2 => { color => 'contigblue1', name => 'CD4_H3K9me2', plot => 'hist' },
           CD4_H3K9me3 => { color => 'contigblue1', name => 'CD4_H3K9me3', plot => 'hist' },

           CD4_H3R2me1 => { color => 'contigblue1', name => 'CD4_H3R2me1', plot => 'hist' },
           CD4_H3R2me2 => { color => 'contigblue1', name => 'CD4_H3R2me2', plot => 'hist' },

           CD4_H4K20me1 => { color => 'contigblue1', name => 'CD4_H4K20me1', plot => 'hist' },
           CD4_H4K20me3 => { color => 'contigblue1', name => 'CD4_H4K20me3', plot => 'hist' },

           CD4_H4R3me2 => { color => 'contigblue1', name => 'CD4_H4R3me2', plot => 'hist' },

           CD4_PolII => { color => 'contigblue1', name => 'CD4_PolII', plot => 'hist' },

           RegulatoryFeatures =>   { color => 'red3', name => 'RegulatoryFeatures' }
);

#What other options do we want here
#restart server?
#rehash config? Always do this?

#Shall we let the env handle cating the individual db config and the header?
#The cat'd file needs to be named after the $das_instance


GetOptions(
   	   #Run modes
		   'no_headers'   => \$no_headers,
		   'only_headers' => \$headers_only,

		   #Source DB params
		   'dbhost=s' => \$dbhost,
		   'dbuser=s' => \$dbuser,
		   'dbpass=s' => \$dbpass,
		   'dbport=i' => \$dbport,
		   'dbname=s' => \$dbname,
		   'species=s'=> \$species,
		   'dnadb_host=s' => \$dnadb_host,
		   'dnadb_user=s' => \$dnadb_user,
		   'dnadb_pass=s' => \$dnadb_pass,
		   'dnadb_port=i' => \$dnadb_port,
		   'dnadb_name=s' => \$dnadb_name,
		   'hydra_instance=s' => \$hydra_instance,
		   
		   #DAS Server params
		   'das_host=s'   => \$das_host,
		   'das_config=s' => \$das_config,
		   'das_port=i'   => \$das_port,#$ENV{'EFG_DAS_PORT'} 9876?
		   'maxclients=i' => \$max_clients,
		   'das_name=s'   => \$das_name,
		   'styleshome=s' => \$styleshome,
		   'maintainer=s' => \$maintainer,
		   'severroot=s'  => \$serverroot,

		 		   
		   #Coord system config
		   'features_region=s' => \$features_region,
		   'assemblies=s{,}'   => \@cs_versions,

		   #'default_colour=s' => \$default_colour,

		   #Individual set
		   #Change this to have result_set feature_set
		   #Then move checking from shell func to here
		   #Or just leave as we're moving to hydra
		   'set_type=s'   => \$set_type,
		   'set_name=s'   => \$set_name,
		   'set_colour=s' => \$set_colour,#Not yet implemented
		   'set_plot=s'   => \$set_plot,#Not yet implemented
		   'set_display_name=s' => \$set_display_name,#Not yet implemented
		   #'set_link_gene=s'
		   #'set_link_region=s

		   #Hydra options
		   'not_hydra'         => \$not_hydra,
		   'source_types=s{,}' => \@source_types,#This pushes, so set default after

		   #HTML/Feature link params
		   #'local_port'    => \$local_port,#Will this work, isn't the data integrated on the server side?
		   'link_gene=s' => \$link_gene,#Not yet implemented
		   'link_region=s' => \$link_region,#Not yet implemented

		   #Helper params
		   #'tee'                    => \$main::_tee,#Not yet implemented
		   #'filename'               => \$main::_log_file,#Not yet implemented
		   
		   #add a reduced log to minimize memory usage?
           'help'                   => sub { pod2usage(-exitval => 0, -message => $pod_params); },
		   'man'                    => sub { pod2usage(-exitval => 0, -message => $pod_params, verbose => 2); },
		  ) or pod2usage(
						 -exitval => 1,
						 -message => $pod_params
						);

$| = 1;

#Set some more defaults
@source_types = keys(%types) if ! @source_types;


if(! defined $main::_log_file){
  $main::_no_log = 1;
  $main::_tee    = 1;
}

# Some header only mandatory params
if(! ($das_config && -d $das_config)){

  pod2usage(
			-exitval => 1,
			-message => "Could not find -das_config directory:\t$das_config",
		   )
}

# Mandatory DB params
if(! $headers_only && 
   ((! ($dbhost && $dbuser && $dbname && $das_host && $species)) ||
	($set_name && ! $set_type) ||
	($set_type && ! $set_name))){
  pod2usage(
			-exitval => 1,
			-message => "You have omited some mandatory parameters required for generating DAS source configurations\n$pod_params",
		   )
}

#Validate/Over-ride das_host


my $uname = `uname`;
chomp $uname;

my $hostname = ($uname eq 'Darwin') ? `hostname` : `hostname -f`;
chomp $hostname;

if($hostname ne $das_host){
  warn "WARNING:\t You have specified the das_host $das_host but appear to be running from $hostname
WARNING:\tMaybe you want to reset change the das_host or run this from a different host?\n";
}


#Check we are not trying to turn on standard sources with hydra mode

if(! $not_hydra){

  if($set_type || $set_name){
	die('Cannot yet specify individual set(-set_type|name) config with hydra source');
	#We would have to do this by matching the source name to the table names somehow?
	#Or maybe specifying a config string which could be split?
	#Would still have to match to the table names, so probably best to set up separate sources
  }

  foreach my $source_type(@source_types){
	
	if(! exists $types{$source_type}){
	  die("Hydra source type $source_type is not valid. Please use one of the following:\t".
		  join("\t", keys(%types)));
	}
  }
}
#else{#Standard only validation?
#
#}

#Do some set validation here
# plot?

if ($link_region && $link_gene){
  die('Must specify only one link e.g -link_region 2:191541121-191588181 or -link_gene STAT1');
}

#Could do with validating the region here
$location = "gene=${link_gene}" if $link_gene;
$location = "r=${link_region}"  if $link_region;
 

if(! $no_headers){
  
  if(! $das_name){
	die('You must provide an -instance name to write das configuration header files');
  }

  my $prefork = int($max_clients/2);
  my $das_instance="${das_name}.${das_host}.${das_port}";
  $maintainer = "maintainer   = $maintainer\n";

  my $header_file = $das_config."/${das_instance}.config.header";
  print ":: Generating DAS ini header:\t$header_file\n";
  open (OUT, ">$header_file") || die("Cannot open header file:\t$header_file");

  #Do we need both prefork and maxclients?
  #Is prefork deprecated?
  
  #hard coded coords here, but need to change to slice defined by $location

  #Is this absolutely necessary?
  #This appears only to be required if we want the default styles/coords home
  #Why does this not default to this anyway?
  if(! -d $serverroot){
	warn "serverroot directory does not exist:\t$serverroot\n".
	  'Default styles/coordhome will not be found?\n';
  }
  else{
	$serverroot = "serverroot   = $serverroot";
  }


  if(! -d $styleshome){
	die('Could not find -styleshome directory:\t$styleshome');
  }

  

  #Do we need to add ensemblhome or bioperlhome to this?
  #Only if we launch with these in the PERL5LIB!
  #Do we need to add coordshome here?

  print OUT "[general]
prefork      = ${prefork}
maxclients   = ${max_clients}
port         = ${das_port}
hostname     = ${das_host}
styleshome   = ${styleshome}
pidfile      = ${das_config}/${das_instance}.pid
logfile      = ${das_config}/${das_instance}.log
${maintainer}${serverroot}\n\n";

#;response_hostname=das.example.com
#;response_port=80
#;response_protocol=https
#;response_baseuri=/frontend
#;oraclehome=/usr/local/oracle
#;ensemblhome=/usr/local/ensembl
  close(OUT);

  $header_file = $das_config."/${das_instance}.html.header";
  print ":: Generating HTML header:\t$header_file\n";
  open (OUT, ">$header_file") || die("Cannot open header file:\t$header_file");

  
  print OUT <<EOHTML;
<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">
<html xmlns="http://www.w3.org/1999/xhtml" xml:lang="en-gb"  lang="en-gb">
<head><title>Ensembl Functional Genomics DAS links</title></head>
<body>
<h3>DAS data sources available on $das_host:$das_port</h3>
EOHTML

  close(OUT);
}


  my ($sets, $set_class, $name, $desc, $display_name, $colour, $plot);
my ($sources_file, $html_file, $type);



#Potential to want to use two different dnadbs for the same efg db if we have 
#More than one coord system
#This would require building two config files for each cs and dnadb
#Naming these differently so we don't overwrite
#Leave this for now as we can do it by just editing manually the config
#And it is most likely that both CSs will be in the same DB if it is Mouse, Rat or Human

if(! $headers_only){
 
  $ini_password = "\npassword          = $dbpass" if $dbpass;

  # Set DBs & Adaptors
  if($dnadb_name){
	$dnadb_pass ||= $dbpass;#used in config

	$dnadb = Bio::EnsEMBL::DBSQL::DBAdaptor->new
	  (
	   -host    => $dnadb_host  || $dbhost,
	   -user    => $dnadb_user  || $dbuser,
	   -dbname  => $dnadb_name,
	   -pass    => $dnadb_pass,
	   -port    => $dnadb_port  || 3306,
	   -species => $species,
	   -group   => 'core',
	);
	$dnadb->dbc->db_handle;	#Test connection
  }

  my $db = Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor->new(
     -host    => $dbhost,
     -user    => $dbuser,
     -dbname  => $dbname,
     -pass    => $dbpass,
     -port    => $dbport,
     -dnadb   => $dnadb,
	 -species => $species,
	 -group   => 'funcgen',
     );

  $db->dbc->db_handle;#Test connection
  

  #Set DAS species name and default assembly
  my $das_species = ucfirst($species);
  $das_species =~ s/_/ /;
  
  
  if(@cs_versions){#validate
	
	foreach my $cs_version(@cs_versions){

	  my $found_cs = 0;
	
	  foreach my $cs (@{$db->dnadb->get_CoordSystemAdaptor->fetch_all_by_name('chromosome')}){
		$found_cs = 1 if $cs_version eq $cs->version;
	  }
	
	  if(! $found_cs){
		die("Could not find the CoordSystem $cs_version in the dnadb $dnadb_name");
	  }
	}
  }
  else{#get default
	
	foreach my $cs (@{$db->dnadb->get_CoordSystemAdaptor->fetch_all_by_name('chromosome')}){
	  
	  @cs_versions = ($cs->version) if $cs->is_default;
	}
  }
  
  map $_ =~ s/([0-9]+)/_$1/, @cs_versions;
  my $coords = join(",Chromosome,$das_species -> $features_region; ", @cs_versions).
		",Chromosome,$das_species -> $features_region;";

  if(! $not_hydra){

	warn "WARNING:\tCannot generate source attachment links for hydra sources\n";


	$sources_file = $das_config."/${das_name}.${dbhost}.${dbport}.${dbname}.hydra.sources";
	print ":: Generating DAS @source_types ini sources:\t$sources_file\n";
	open (OUT, ">$sources_file") || die("Cannot open sources file:\t$sources_file");

	foreach my $type(@source_types){
	  
	  #To allow multiple source DBs
	  #We need the name to define the source DB and type uniquely
	  #define types separately so we can turn them off in the config without hacking the code
	  #Do we need to turn them off without dropping them from the DB?
	  #Move the table away from the DB OR implement some status checking?
	  #This info will however be visible to user!!!

	
	  #Do we need separate hydra/adaptor classes for bed/feature_set/result_set?

	  #basename is tables hydra dbi looks for in DB using mysql like "$basename%"	

	  #Do we want to be able to change this prefix?
	  #Or can we drop it all together?
	  
	  #foreach my $cs_version(@cs_versions){
		
	  my $source_name = $hydra_instance || $dbname.'@'.$dbhost.':'.$dbport;
	  $source_name = 'eFG_'.$type."s:".$source_name;

	  print OUT "\n[${source_name}]
state           = on\n".
  $types{$type}{basename}.
	"adaptor         = ".$types{$type}->{'adaptor'}."
hydra           = ".$types{$type}->{'hydra'}."
transport       = ".$types{$type}->{'transport'}."
set_type        = $type
host            = $dbhost
port            = $dbport
species         = $species
username        = ${dbuser}${ini_password}
dbname          = $dbname
dnadb_host      = ".$db->dnadb->dbc->host."
dnadb_port      = ".$db->dnadb->dbc->port."
dnadb_user      = ".$db->dnadb->dbc->username."
dnadb_pass      = ".$db->dnadb->dbc->password."
dnadb_name      = ".$db->dnadb->dbc->dbname."
autodisconnect  = no
coordinates     = $coords
\n\n";

#;skip_registry = 1
#;category      = sequencing
#;method        = Solexa 1G	
#;basename      = 
	  
	  
	}
	close(OUT);

  }
  else{#Old standard source adaptor setup


	my %adaptors = (
					result     => $db->get_ResultSetAdaptor,
					feature    => $db->get_FeatureSetAdaptor,
					#experimental => $db->get_ExperimentalSetAdaptor,
					#We need this to point to one file or a blob in the db
					#Need to extend the schema/API to handle this?

					

				   );
	
	
	#Generate sources file for each type

	$sources_file = $das_config."/${das_name}.${dbhost}.${dbport}.${dbname}.standard.sources";
	print ":: Generating DAS ini sources:\t$sources_file\n";
	open (OUT, ">$sources_file") || die("Cannot open sources file:\t$sources_file");
	  
	  $html_file = $das_config."/${das_name}.${dbhost}.${dbport}.${dbname}.standard.html";
	  print ":: Generating DAS source links:\t$html_file\n";
	  open (HTML, ">$html_file") || die("Cannot open html file:\t$html_file");


	print HTML "<h2>Served from $dbname</h2>
<p>selected location: <a href='http://www.ensembl.org/${species}/contigview?${location}'>${location}</a></p>";
	
	foreach my $aname(keys %adaptors){
	  $plot = '';
	  $set_class = ucfirst($aname).'Set';
	  
	  #$sources_file = $das_config."/${das_name}.${dbhost}.${dbport}.${dbname}.${set_class}.sources";
	  #print ":: Generating DAS $set_class ini sources:\t$sources_file\n";
	  #open (OUT, ">$sources_file") || die("Cannot open sources file:\t$sources_file");
	  
	  #$html_file = $das_config."/${das_name}.${dbhost}.${dbport}.${dbname}.${set_class}.html";
	  #print ":: Generating DAS $set_class source links:\t$html_file\n";
	  #open (HTML, ">$html_file") || die("Cannot open html file:\t$html_file");
	  
	  $sets = $adaptors{$aname}->fetch_all;
	  
	  foreach my $set (sort @{$sets}) {
		$name = $set->name;
		
		#Skip if this is not the -set_name exists in set_config or does not has DAS_DISPLAYABLE status
		if ( (! (($set_name && ($set_name eq $name)) || (exists $set_config{$name}))) &&
			 ! $set->adaptor->has_stored_status('DAS_DISPLAYABLE', $set)){
		  next;
		}
		
		if (! $set->adaptor->has_stored_status('DAS_DISPLAYABLE', $set)){
		  #Do this here or in the env?
		  
		  #Need to eval this
		  eval { $set->adaptor->store_status('DAS_DISPLAYABLE', $set); };
		  
		  warn "Failed to store DAS_DISPLAYABLE status for $set_class:\t$name\n" if $@;
		}
		
	  
		
		#Need to handle individual set override here
		
		$desc = $set->display_label;
		$display_name = (exists $set_config{$name}{name})   ? $set_config{$name}{name}   : $name;
		$colour       = (exists $set_config{$name}{colour}) ? $set_config{$name}{colour} : $default_colours{$set_class};
		#$type = ($set_class eq 'result_set') ? 'result' : $set->type;
		
		#my $desc = "[$species - $cell_type] $dsn ($fset_type feature set)";
		#my $desc = "[$species - $cell_type] $dsn ($type set)";

		print  ":  Configuring $set_class\t$name\n";
		
		#Also build config array here to print back to file.



		#Do we need feature_query?
		#feature_query= field0 = %segment and field2 >= %start and field1 <= %end
		#coordinates  = NCBI_36,Chromosome,Homo sapiens -> X:1,2000000

		#Can we config style sheets from here?
		#TYPE id= $type?
		#annotated and DNA Methlyation are being set in ensembl_feature_set
		
		#Also need to add config here such that we can build das_xsl necessary for ensembl to auto configure tracks
		#Speak to James/Andy about necessary elements for species, coordinate system
		
		
   		#Need to test each set for DISPLAYABLE_CSVERSION

		#foreach my $cs_version(@cs_versions){

		  print OUT "\n[${display_name}]
state             = on
adaptor        = ".$types{$aname}->{'adaptor'}."
transport      = ".$types{$aname}->{'transport'}."
host              = $dbhost
port              = $dbport
dbname            = $dbname
species           = $species
set_type          = $aname
username          = ${dbuser}${ini_password}\n".#description       = $desc #Now handled by adaptor
#remove set_name an just use set_id and dsn?
"set_name          = $name
set_id            = ".$set->dbID."
dnadb_host      = ".$db->dnadb->dbc->host."
dnadb_port      = ".$db->dnadb->dbc->port."
dnadb_user      = ".$db->dnadb->dbc->username."
dnadb_pass      = ".$db->dnadb->dbc->password."
dnadb_name      = ".$db->dnadb->dbc->dbname."
coordinates     = $coords
";
		
#source            = SOURCE
#category          = CATEGORY
		}
		#feature_query     = $location
		#What is CATEGORY, SOURCE and TYPE? TYPE was set_type for feature_set
		#Are these das style sheet config options?
		

		#TYPE = result, annotated, external, regulatory?
		#Do we even need this as we can access the FeatureAdaptor from the Set it's self?
		#Can we omit or put in header if hardcoded?
		

		#set_name is included so we can manually edit the source name
		#But have an easy readable reference to the actual set
		#Use the dbID for fetching in the adaptor
		#feature_set       = $name\n";
		
		#These links are now used to attach the sources with a given format
		#We should really some more config and let the adaptor pass this format info as xss
		
		#result_sets
		
		
		if (exists $set_config{$name}{plot}){
		  
		  if(! exists $plots{$set_config{$name}{plot}}){
			die("You have specified an invalid plot type for $set_class $name:\t".$set_config{$name}{plot});
		  }
		  $plot = $plots{$set_config{$name}{plot}};
		}
		




		#This is currently not working!
	  #print HTML '<p><a href="'.
	#	'http://www.ensembl.org/'.$species.
	#	  '/contigview?'.$location.';'.
	#		'add_das_source=(name='.$display_name.'+url=http://'.$das_host.':'.$das_port.
	#		  '/das+dsn='.$display_name.
	#			'+type=ensembl_location_chromosome'.
	#			  $plot.
	#				'+color='.$colour.'+strand=r+labelflag=n'.
	#				  '+group=n+'.$display_params{$set_class}.'active=1)">'.
	#					$desc.'</a></p>'."\n";
		
		
		#No need to set_config for individual set if defined
		
		
	#}
	}
	close(OUT);
	close(HTML);
  }

  
  #Now need to rewrite config file if we have a new set_name defined.
  #Do we need to backup old file first? Maybe just to tmp.filename, 
  #rather than one for each process, so we have just one backup.

}






#DAS_HOST This will need checking against run host
#DAS_HOME

#Will this host be exposed? Do we have to set up a forward for this?
#We need rto define a config dir in the environment
#This will contain individual config for each instance which will be cat'd together to get the DAS
#server to point at different mysql instances?

#We need to pass some extra info here

#Do we need $sepcies in here?
#If not then we can remove it unless we do not specify a $dnadb_name



1;
