#!/usr/bin/env perl

=head1 LICENSE

  Copyright (c) 1999-2012 The European Bioinformatics Institute and
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

create_ftp_dir_links.pl -- Creates ftp dir structure and links for nfs files

=head1 SYNOPSIS

create_ftp_dir_links.pl -dbhost $DB_HOST -dbuser $DB_USER -dbname $DB_NAME [ -dbport $DB_PORT -dbpass $DB_PASS -ftp_root $FTP_DIR ]

=head1 DESCRIPTION

Creates ftp dir structure and links for nfs files i.e. collection files. Queries the DB
to identify which files are current and creates links to relevant nfs assembly specififc
directories.

NOTE: update_DB_for_release.pl and healtchecks should be run before this script.


=head1 PARAMETERS

=over

=head4 B<Mandatory>

=item B<-dbhost>       Database host server

=item B<-dbname>       Name of the database

=back

=head4 B<Optional>

=over

=item B<-dbuser>       Username of the database. Default = ensro.

=item B<-dbpass>       Password for the database user

=item B<-dbport>       Port of the MySQL server. Default = 3306.

=item B<-dnadb_host>   Core database host server

=item B<-dnadb_user>   Core database username

=item B<-dnadbdb_pass> Core database user password

=item B<-dnadb_port>   Core database port

=item B<-dnadb_name>   Core database name

=item B<-ftp_root>     Root directory for ftp staging dir. Default = /lustre/scratch103/ensembl/funcgen/output/ftp

=item B<-help>         Prints a short help message

=item B<-man>          Prints a short help man page

=back

=cut

#To do

# Add support for FeatureSets? Would we ever want flat files for these?


use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Registry;
use Cwd;

my @tmp_args = @ARGV;
my ($host, $pass, $dbname, $dnadbhost, $dnadbport, $dnadbuser, $dnadbname, $dnadbpass);

#Default values
my $force    = 0;
my $port     = 3306;
my $user     = 'ensro';
my $ftp_root = '/lustre/scratch103/ensembl/funcgen/output/ftp';

#get command line options

my @set_types = ('result');#, 'feature');

print "create_ftp_dir_links.pl @tmp_args\n";

GetOptions (
			'dnadb_host=s'       => \$dnadbhost,
			'dnadb_user=s'       => \$dnadbuser,
			'dnadb_port=i'       => \$dnadbport,
			'dnadb_pass=s'       => \$dnadbpass,
			'dnadb_name=s'       => \$dnadbname,
			'dbhost=s'           => \$host,
			'dbuser=s'           => \$user,
			'dbport=i'           => \$port,
			'dbpass=s'           => \$pass,
			'dbname=s'           => \$dbname,
			'ftp_root=s'         => \$ftp_root,
			'force'              => \$force,
			"help|?"             => sub { pos2usage(-exitval => 0, -verbose => 1, -message => "Params are:\t@tmp_args"); },
			"man|m"              => sub { pos2usage(-exitval => 0, -verbose => 2, -message => "Params are:\t@tmp_args"); },

		   )  or pod2usage( -exitval => 1 ); #Catch unknown opts


if(!$host || !$user || !$dbname )  {  print "Missing connection parameters (use -h for help)\n"; exit 1; }

#Create database connections


 
my $efgdb = Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor->new
  (
   -host    => $host,
   -port    => $port,
   -user    => $user,
   -dbname  => $dbname,
#   -species => $species, #Not strictly necessary
   -pass    => $pass,
   
   -dnadb_host   => $dnadbhost,
   -dnadb_port   => $dnadbport,
   -dnadb_user   => $dnadbuser,
   -dnadb_name => $dnadbname,
   -dnadb_pass   => $dnadbpass,
  );


#Test connections
$efgdb->dbc->db_handle;
$efgdb->dnadb->dbc->db_handle;
my $reg = 'Bio::EnsEMBL::Registry';


#Test schema_version/species/assembly and create release dir
my $metac          = $efgdb->get_MetaContainer;
my $schema_version = $metac->get_schema_version;
my $species        = $metac->single_value_by_key('species.production_name');

if(! $schema_version ||
   ($dbname !~ /_${schema_version}_/)){
  die("Meta schema.version is not present in meta($schema_version) or does not match the dbname($dbname)");
}

if(! $species ||
   ($dbname !~ /^${species}_/)){
  die("Meta species.production_name  is not present in meta($species) or does not match the dbname($dbname)");
}

my @assm_vers = @{$efgdb->dbc->db_handle->selectall_arrayref('SELECT distinct(version) from coord_system where is_current=1 and version is not NULL AND version !=""')};
@assm_vers = @{$assm_vers[0]};


if($#assm_vers != 0){
  die("Cannot find unique current coord_system.version(@assm_vers). Have you run update_DB_for_release.pl?")
}

my $assm = $assm_vers[0];
my $rel_dir = $ftp_root."/release-${schema_version}/data_files/${species}/${assm}";


if(! -d $rel_dir){
  system("mkdir -p $rel_dir") == 0 or die("Failed to create FTP release dir:\t$rel_dir");
}

#Now check we have top level data_file dir link to nfs

my $df_dir_link = $ftp_root."/data_files/${species}/${assm}";

if(! -l $df_dir_link){ #Let's create it
  print "Top level data_files dir link does not exist:\t$df_dir_link\n";
  

  #Should grab dbfile_data_root here
  my $dbfile_data_root    = $efgdb->get_MetaContainer->single_value_by_key('dbfile.data_root');
 
  if(! ($dbfile_data_root &&
		-d $dbfile_data_root) ){
	die("Unable to find dbfile.data_root:\t$dbfile_data_root");
  }
  else{
	my $cmd = "ln -s $dbfile_data_root $df_dir_link";
	print "Creating top level data_files dir link:\t$cmd\n";
	system($cmd) == 0 or die("Failed to create top level data_file dir link");
  }
}




#Loop through set types and create links
my $dir = getcwd;
my ($feature_class, $source_dir, $target_dir);



#Log errors and caryr on or die straight away?
#Use Helper?

foreach my $set_type(@set_types){

  #warn "getting $set_type setadaptor";
  #my $set_adaptor = $reg->get_adaptor($species, 'funcgen', $set_type.'setadaptor');

  my $method      = 'get_'.ucfirst($set_type).'SetAdaptor';
  my $set_adaptor = $efgdb->$method;
  

  my $num_dirs = 0;

  foreach my $set (@{$set_adaptor->fetch_all}){

	#Check if is DISPLAYABLE and have data files

	#This is ResultSet specific after here
	($source_dir = $set->dbfile_data_dir) =~ s'/+'/'g;
	#Let's tidy up the source dir as we have numerous multiple slashes ///
	
	if( ! ($set->is_displayable &&
		   $source_dir)){
	  next;
	}
	 
	#We don't actually want to link directly to the data files, but relatively to
	#../../data_files/SPECIES/ASSEMBLY
	#This should also contain links, but to /nfs/ensnfs-dev/staging
	#This should already be in place but need to add check in here somewhere

	



	$feature_class = ($set->set_type eq 'result') ? 'result_feature' : $set->feature_class.'_feature';
  	my $fclass_dir = $rel_dir."/${feature_class}";
	
	if(! -d $fclass_dir){
	  system("mkdir -p $fclass_dir") == 0 or die("Failed to create FTP release dir:\t${fclass_dir}");
	}
	
	chdir($fclass_dir) || die("Failed to move to FTP release dir:\t${fclass_dir}");;


	($target_dir = $source_dir) =~ s'/$'';
	($target_dir = $target_dir) =~ s'/.*/'';
	
	#Now we need to redefine source_dir as link to data_files relative to target dir
	$source_dir = "../../../../../data_files/${species}/${assm}/${feature_class}/${target_dir}";
	

	#SANITY CHECKING
	
	if(! -d $source_dir){ #Is is a directory?
	  #This will be a directory as the link is at the data_file level
	  die("Source dbfile_data_dir does not exist for ".$set->name."\t".$source_dir);
	}
	else{                 #Does is contain any data?
	  opendir(DirHandle, $source_dir) || die("Failed to open source dir:\t$source_dir");

	  #Need to catch error here

	  my $num_files = 0;

	  while (readdir(DirHandle)) {
		$num_files++;
		#Could check for expected suffixes or non-empty files
		#check for name match
		#But we are verging on a HC here.
	  }
	  
	  closedir(DirHandle);

	  if(! $num_files){
		die("Found 0 files in source directory:\t${source_dir}");
	  }
	}	


	if(-e $target_dir &&
	   (! -l $target_dir) ){
	  die("$target_dir already exists but is not a link to $source_dir");
	}

	my $cmd = "ln -sf $source_dir $target_dir";
	system("$cmd") == 0 or die("Failed to create link using:\t${cmd}");
	#print "$cmd\n";
	$num_dirs ++;
  }

  print "Linked $num_dirs $set_type directories\n";
}
	
	
chdir($dir);


print "Now print command to send to webteam which will follow links and rsync to the ftp dir\n";
