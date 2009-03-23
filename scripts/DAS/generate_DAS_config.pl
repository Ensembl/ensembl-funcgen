#!/software/bin/perl

=head1 NAME

generate_DAS_config.pl

=head1 SYNOPSIS

generate_DAS_config.pl \
    -dbhost dbhost \
    -dbport 3306 \
    -dbuser ensro \
    -dbname homo_sapiens_funcgen_47_36i \
    -species homo_sapiens \
    -das_config $EFG_SRC/config/DAS\
    -das_name efg \
    -das_host DAShost \
    -das_port  9000

=head1 DESCRIPTION

This script writes DAS configuration for all DAS_DISPLAYABLE sets from a given DB.

=head1 OPTIONS

 Mandatory: If running without -only_header (default).
    DB parameters
    -dbhost  DB host name
    -dbport  DB port
    -dbuser  DB user name
    -dbname  Name of DB
    -species Latin name of species as used in DB name

 Mandatory: If running without -no_header (default).
    DAS parameters
    -das_config
    -das_name 
    -das_host 
    -das_port 

 Optional:
    -dbpass      DB password
    -maxclients  Maximum DAS clients

    DNA DB parameters, default is to use ensembldb.
    -dnadb_host  DB host name
    -dnadb_port  Core DB port
    -dnadb_user  Core DB user name
    -dnadb_name  Name of core DB
    -dnadb_pass  Core DB password
    
    Run modes
    -no_header   Only prints source config
	-header_only Only prints DAS server config
    -help        Prints this documentation and exits

    Individual Set handling
    -set_type        Set type e.g. result, feature, regualtory
	-set_name        Name of set
    #-set_colour      Colour of track for given set e.g. contigblue1, contigblue2, red3
    #-set_plot        Plot type for given set e.g. hist or tiling

	HTML link params, used to generate page with attachement links
	#-link_gene    Name og gene
	or
 	#-link_loci    Loci e.g. 19:1243:1567??? Is this correct?

    #Add more here for default colours?


=head1 LICENCE

This code is distributed under an Apache style licence. Please see
http://www.ensembl.org/info/about/code_licence.html for details.

=head1 CONTACT

Please post comments/questions to the Ensembl development list
<ensembl-dev@ebi.ac.uk>

=cut

#This needs to use the DAS DISPLAYABLE status to auto set up for a given list of hosts
#Or use default available config host file to poll the existing DB
#Or add a given dbhost using either the 
#Then we cat the individual host files with the general section to make the complete ini file
#What about losing sources that have been set up explicitly without DAS_DISPLAYABLE?
#This script should set DAS_DISPLAYABLE!

#We need func to list sources given host, cell/feature type, experiment name?


use strict;
use warnings;
use Data::Dumper;
use Pod::Usage;
use Getopt::Long;
#use Bio::EnsEMBL::Registry;
#use Bio::EnsEMBL::Utils::Exception qw(verbose throw warning info);
#use Bio::EnsEMBL::Analysis::Tools::Logger qw(logger_verbosity logger_info);
#Use Helper instead of Logger?
#use Bio::EnsEMBL::Funcgen::FeatureSet;

use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor;

my ($dbhost, $dbport, $dbuser, $dbpass, $dbname, $das_port, $das_host, $species);
my ($das_config, $set_type);
my ($dnadb_host, $dnadb_port, $dnadb_user, $dnadb_pass, $dnadb_name, $dnadb, $set_name);
my (@adaptor_names, $no_headers, $headers_only, $link_gene, $link_loci);
my ($set_colour, $set_plot, $set_display_name);
my $ini_password= '';
my $max_clients = 5;
my $das_name = 'efg';
my $location = 'gene=STAT1';
#$location = 'gene='.$opts{g} if ($opts{g});
#$location = 'c='.$opts{c} if ($opts{c});

my $pod_params = "Params are:\t".join(' ', @ARGV);

#Need to print this to a config file which can be read/updated by this script or edited by hand
#for each instance

my %plots = (
			 hist   => '+score=h',
			 cgrad  => '+score=c+fg_grades=100+fg_data=n+fg_max=100',
			 tiling => '+score=t',
			);

my %default_colours = (
					   ResultSet  => 'red',
					   FeatureSet => 'blue',
					  );

my %display_params = (
				  ResultSet  => '+score=s+fg_data=o+',
				  FeatureSet => '',
				 );

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
		   'headers_only' => \$headers_only,

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
		   
		   #DAS Server params
		   'das_host=s'   => \$das_host,
		   'das_config=s'   => \$das_config,
		   'das_port=i'   => \$das_port,#$ENV{'EFG_DAS_PORT'} 9876?
		   'maxclients=i' => \$max_clients,
		   'das_name=s'   => \$das_name,

		   
		   #'default_colour=s' => \$default_colour,

		   #Individual set
		   'set_type=s'   => \$set_type,
		   'set_name=s'   => \$set_name,
		   'set_colour=s' => \$set_colour,#Not yet implemented
		   'set_plot=s'   => \$set_plot,#Not yet implemented
		   'set_display_name=s' => \$set_display_name,#Not yet implemented

		   #HTML link params
		   'link_gene=s' => \$link_gene,#Not yet implemented
		   'link_loci=s' => \$link_loci,#Not yet implemented

		   #Helper params
		   #'tee'                    => \$main::_tee,#Not yet implemented
		   #'filename'               => \$main::_log_file,#Not yet implemented
		   
		   #add a reduced log to minimize memory usage?
           'help'                   => sub { pos2usage(-exitval => 0, -message => $pod_params); },
		  ) or pod2usage(
						 -exitval => 1,
						 -message => $pod_params
						);

$| = 1;


if(! defined $main::_log_file){
  $main::_no_log = 1;
  $main::_tee    = 1;
}

#Do we need these? Or just use Helper?
#my $utils_verbosity = 'WARNING';
#my $logger_verbosity = 'OFF';
#verbose($utils_verbosity);
#logger_verbosity($logger_verbosity);


#This needs to be done in env, as we never restart from this script!
#my $local_host = `hostname --long`;


# Some header only mandatory params
if(! ($das_config && -d $das_config)){
  pod2usage(
			-exitval => 1,
			-message => $pod_params,
		   )
}

#Do some set validation here plot?

if(! $no_headers){
  my $prefork = int($max_clients/2);
  my $das_instance="${das_name}.${das_host}:${das_port}";

  my $header_file = $das_config."/${das_instance}.config.header";
  print ":: Generating DAS ini header:\t$header_file\n";
  open (OUT, ">$header_file") || die("Cannot open header file:\t$header_file");

  #Do we need both prefork and maxclients?
  #Is prefork deprecated?

  print OUT "[general]
prefork=${prefork}
maxclients=${max_clients}
port=${das_port}
hostname=${das_host}
pidfile=${das_config}/${das_name}.pid
logfile=${das_config}/${das_name}.log";

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
<h3>DAS data sources for</h3>
<h2>$dbname on $das_host:$das_port</h2>
<p>selected location: <a href='http://www.ensembl.org/${species}/contigview?${location};">${location}</a></p>
EOHTML

  close(OUT);
}


  my ($sets, $set_class, $name, $desc, $display_name, $colour, $plot);
my ($sources_file, $html_file, $type);


if(! $headers_only){

  $ini_password = "\npassword = $dbpass\n" if $dbpass;

  # Mandatory DB params
  if((! ($dbhost && $dbuser && $dbname && $das_host && $species)) ||
	 ($set_name && ! $set_type) ||
	 ($set_type && ! $set_name)){
	pod2usage(
			  -exitval => 1,
			  -message => $pod_params,
			)
  }

  # Set DBs & Adaptors
  if($dnadb_name){
	
	$dnadb = Bio::EnsEMBL::DBSQL::DBAdaptor->new
	  (
	   -host    => $dnadb_host  || $dbhost,
	   -user    => $dnadb_user  || $dbuser,
	   -dbname  => $dnadb_name,
	   -pass    => $dnadb_pass  || $dbpass,
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

  my %adaptors = (
				  result  => $db->get_ResultSetAdaptor,
				  feature => $db->get_FeatureSetAdaptor,
				  #experimental => $db->get_ExperimentalSetAdaptor,
				  #We need this to point to one file or a blob in the db
				  #Need to extend the schema/API to handle this?
				 );


 
  #Generate sources file for each type
 

  foreach my $aname(keys %adaptors){
	$plot = '';
	$set_class = ucfirst($aname).'Set';
	$sources_file = $das_config."/${das_name}.${dbhost}.${dbport}.${dbname}.${set_class}.sources";
	print ":: Generating DAS $set_class ini sources:\t$sources_file\n";
	open (OUT, ">$sources_file") || die("Cannot open sources file:\t$sources_file");

	$html_file = $das_config."/${das_name}.${dbhost}.${dbport}.${dbname}.${set_class}.html";
	print ":: Generating DAS $set_class source links:\t$html_file\n";
	open (HTML, ">$html_file") || die("Cannot open html file:\t$html_file");

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
	  $type = ($set_class eq 'ResultSet') ? 'result' : $set->type;

	  #my $desc = "[$species - $cell_type] $dsn ($fset_type feature set)";
	  #my $desc = "[$species - $cell_type] $dsn ($type set)";

	  print  ":  Configuring $set_class\t$name\n";

	  #Also build config array here to print back to file.


	  print OUT "\n[${display_name}]
state             = on
adaptor           = $set_class
transport         = ensembl_funcgen
host              = $dbhost
port              = $dbport
dbname            = $dbname
source            = SOURCE
type              = $type
category          = CATEGORY
username          = ${dbuser}${ini_password}
description       = [ $species ] $desc
set_name          = $name
set_id            = ".$set->dbID."\n";


#What is CATEGORY, SOURCE and TYPE? TYPE was set_type for feature_set



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

	  
	  print HTML '<p><a href="'.
		'http://www.ensembl.org/'.$species.
		  '/contigview?'.$location.';'.
			'add_das_source=(name='.$display_name.'+url=http://'.$das_host.':'.$das_port.
			  '/das+dsn='.$display_name.
				'+type=ensembl_location_chromosome'.
				  $plot.
					'+color='.$colour.'+strand=r+labelflag=n'.
					  '+group=n+'.$display_params{$set_class}.'active=1)">'.
						$desc.'</a></p>'."\n";


	  #No need to set_config for individual set if defined


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
