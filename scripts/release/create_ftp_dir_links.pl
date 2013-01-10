#!/usr/bin/env perl

=head1 LICENSE

  Copyright (c) 1999-2013 The European Bioinformatics Institute and
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

=item B<-nfs_root>     Root directory for nfs dir. Default = /nfs/ensnfs-dev/staging

=item B<-link_type>    Link types to create e.g. nfs or ftp, default generates both.

=item B<-list_source_only> Doesn't create anything, simply lists the source dirs which would be used

=item B<-help>         Prints a short help message

=item B<-man>          Prints a short help man page

=back

=cut

#To do

# 1 Add support for FeatureSets? Would we ever want flat files for these?

# 2 Add list files, for use by webteam as they need access to true files patch for copying to the mirrors
# we are moving towards only ever linking from a release specific dir.

# 3 Make files read only

# 4


# General notes

# master files should never change between releases without being renamed!
# Hence,  should never need to overwrite in target, does rsync do this anyway i.e. copy identical files
# can we use checksums to avoid this?
# rsync should never delete from target if not in source(should never have to anyway as we should always have all the files on master staging?)
# we should really move toward only release specific files on staging, and back up master somewhere else.



# ultimately live and archive should run off /nfs/ensnfs-master|live
# this should probably contain a simlinked dir for release specific files for convinience (altho could remove this as not needed if we are copying from staging)

# ftp assm copy should be mindful that this simply links to the nfs assm copy, which may in future only have current files (if we archive old files elsewhere)
# that would mean release dir and master dir would effectively be the same!

# what are issues here wrt updating safety
# where is redundancy/back ups?
# snapshots and old file backup is this enough?


#issues around redundancy of feature_class subdirs between groups
#i.e. we don't have regulation or core dirs, just rnaseq, result_feature e.g.

#/nfs/ensnfs-dev/staging/homo_sapiens/GRCh37/rnaseq
#/nfs/ensnfs-dev/staging/homo_sapiens/GRCh37/result_feature

#should probably be

#/nfs/ensnfs-dev/staging/homo_sapiens/GRCh37/core/rnaseq
#/nfs/ensnfs-dev/staging/homo_sapiens/GRCh37/regulation|funcgen/result_feature

#funcgen|regulation naming issues wrt matching db name versus visibility on ftp site?
#should never reach ftp site?


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
my ($list_source_only, $nfs_assm_dir, $nfs_root, $ftp_root);
my $ftp_mnt = '../';

#Default values
#my $force      = 0;
my $port       = 3306;
my $user       = 'ensro';
my $link_type   = 'ftp'; #We don't want to create nfs by default.
my %link_config = (nfs => {}, ftp=>{});#empty config for validating link types
my @link_types = keys %link_config;


#get command line options
#should this be feature types?
my @set_types = ('result');#, 'dna_methylation');
my @feature_classes = ('result_feature', 'dna_methylation_feature');

print "create_ftp_dir_links.pl @tmp_args\n";

GetOptions 
  (
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
   'nfs_root=s'         => \$nfs_root,
   'link_type=s'        => \$link_type,
   'list_source_only'   => \$list_source_only,
   'old_ftp_mnt'        => \$ftp_mnt,
   #'force'              => \$force,
   "help|?"             => sub { pod2usage(-exitval => 0, -verbose => 1); },
   "man|m"              => sub { pod2usage(-exitval => 0, -verbose => 3); },
   
  )  or pod2usage( -exitval => 1 ); #Catch unknown opts, culprit will be printed by GetOptions


$ftp_mnt = '' if $ftp_mnt eq '1';


if($link_type){

  if(! exists $link_config{$link_type}){
    die("You have specified an invalid -link_type($link_type), must be one of:\t@link_types");
  }
  else{
    @link_types = ($link_type);
  }
}


if( ! (defined $nfs_root &&
       -d $nfs_root) ){
  die('You have not specified a valid -nfs_root');
}

if( ! (defined $ftp_root &&
       -d $ftp_root) ){
  die('You have not specified a valid -ftp_root');
}



if(!$host || !$user || !$dbname )  {  die("Missing connection parameters (use -h for help)"); }

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

if(! defined $schema_version ||
   ($dbname !~ /_${schema_version}_/)){
  die("Meta schema.version is not present in meta($schema_version) or does not match the dbname($dbname)");
}

if(! defined $species ||
   ($dbname !~ /(.*_)*${species}_/)){
  die("Meta species.production_name is not present in meta, or does not match the dbname($dbname !~ /^[.*_]${species}_/)");
}

my @assm_vers = @{$efgdb->dbc->db_handle->selectall_arrayref('SELECT distinct(version) from coord_system where is_current=1 and version is not NULL AND version !=""')};
@assm_vers = @{$assm_vers[0]};


if($#assm_vers != 0){
  die("Cannot find unique current coord_system.version(@assm_vers). Have you run update_DB_for_release.pl?")
}

my $assm = $assm_vers[0];

my $dbfile_data_root = $nfs_root.'/'.$species.'/'.$assm;

# CONFIG
#To enable copying of links around, we use source paths relative to the target directory (rather than full paths)

$link_config{nfs} =
  { 
   root     => $nfs_root,
   release_target_path  => $nfs_root."/release_${schema_version}/${species}/${assm}",
   
   #assuming we have cd'd into a subdir of the release_target_path
   #e.g. $nfs_root."/release-${schema_version}/${species}/${assm}/feature_class
   relative_source_path => "../../../../${species}/${assm}",
   
  };

   
$link_config{ftp} =
  {
   root     => $ftp_root,#do we even need this?
   release_target_path  => $ftp_root."/release-${schema_version}/data_files/${species}/${assm}",
   
   #assuming we have cd'd into a subdir of the release_target_path
   #e.g. $ftp_root."/release-${schema_version}/data_files/${species}/${assm}/feature_class
   relative_source_path => "../../../../../data_files/${species}/${assm}",
  };
                                      


#Looping is a bit higgeldy-piggeldy here as we want to dealing with all links for a 
#given data set at the same time
#hence we have to interate over @link_types a few times.

my %num_dirs;

foreach my $link_type(@link_types){
  $num_dirs{$link_type} = 0;

  #Check we have top level data_file dir link to nfs
  #Create top level data_file dir which simply links to nfs
  if($link_type eq 'ftp'){

    #this is wrong below here
    #nee to look into $dbfile_data_root to get feature types?
    #no this needs to come from config or DB
    #as we may have feature types from other DBs in there.


    foreach my $fclass(@feature_classes){

      my $fclass_dir_link = $ftp_root."/data_files/${species}/${assm}/${fclass}";
      my $fclass_nfs_dir  = $dbfile_data_root."/${fclass}";
      
      if(! -l $fclass_dir_link){ #Let's create it
        print "Top level data_files dir link does not exist:\t$fclass_dir_link\n";
      
        #We build dbfile_data_root above, as it has been removed from meta
        
        if(! -d $dbfile_data_root){
          die("dbfile.data_root is not a valid directory:\t$dbfile_data_root");
        }
        else{
          my $cmd = "ln -s $fclass_nfs_dir $fclass_dir_link";
          print "Creating top level data_files dir link:\t$cmd\n";
          system($cmd) == 0 or die("Failed to create ${fclass}_feature data_file dir link");
        }
      }
    }
  }
  
  my $rel_dir = $link_config{$link_type}->{release_target_path};

  if(! -d $rel_dir){
    system("mkdir -p $rel_dir") == 0 or die("Failed to create $link_type release dir:\t$rel_dir");
  }

  #remove old links here?
}



#Loop through set types and create links
#my $dir = getcwd;
my ($feature_class, $ftp_source_dir, $source_dir, $target_dir, @true_paths);
  
foreach my $set_type (@set_types) {
    
  #warn "getting $set_type setadaptor";
  #my $set_adaptor = $reg->get_adaptor($species, 'funcgen', $set_type.'setadaptor');    
  my $method      = 'get_'.ucfirst($set_type).'SetAdaptor';
  my $set_adaptor = $efgdb->$method;
    
  foreach my $set (@{$set_adaptor->fetch_all}) {
      
    #Check if is DISPLAYABLE and have data files
      
    #This is ResultSet specific after here
    ($source_dir = $set->dbfile_data_dir) =~ s'/+'/'g;
    #Let's tidy up the source dir as we have numerous multiple slashes ///
	
    if ( ! ($set->is_displayable &&
            $source_dir)) {
      next;
    }
	 
 
    #get true source here to list
    #will always be on nfs

    #do we even need dbfile_registry?
    #subdir is easily auto-generated
    #from feature_type and set name
    #keep for now for flexibility

    $feature_class = $set->feature_class.'_feature';
    #Strip leading /, to avoid empty string as first element
    ($target_dir = $source_dir) =~ s'^/'';
    $target_dir = (split '/', $target_dir)[1];
    #It should always be the 2nd element, even if we have a full file path
    #if this ever does not match the source dir name, 
    #then we need add the $target_dir back into the ln cmd and
    #implement an overwrite mode which removes the existing link

    push @true_paths, $dbfile_data_root."/${feature_class}/${target_dir}";


    if (! $list_source_only) {

      foreach my $link_type (@link_types) {
        my $rel_dir = $link_config{$link_type}->{release_target_path};
        my $fclass_dir = $rel_dir."/${feature_class}";
        
        if (! -d $fclass_dir) {
          system("mkdir -p $fclass_dir") == 0 or die("Failed to create $link_type release dir:\t${fclass_dir}");
        }
	
        chdir($fclass_dir) || die("Failed to move to $link_type release dir:\t${fclass_dir}");;

	
        #Now we need to redefine source_dir as link to data_files relative to target dir
        $source_dir =  $link_config{$link_type}->{relative_source_path}."/${feature_class}/${target_dir}";
      
   
        #SANITY CHECKING
      
        #This is now dependant on link_type?
        #should do this once really 
     

        if (! -d $source_dir) { #Is is a directory?
          #This will be a directory as the link is at the data_file level
          die("Source dbfile_data_dir does not exist for ".$set->name."\t".$source_dir."\nFrom $fclass_dir");
        } else {                #Does is contain any data?
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
        
          if (! $num_files) {
            die("Found 0 files in source directory:\t${source_dir}");
          }
        }	
      
      
        if (-e $target_dir){
          
          if(! -l $target_dir){
            #this is not quite true as we are not testing the path
            die("$target_dir already exists but is not a link to $source_dir");
          }
          else{ #Assume this is okay
            #count skipped links? or fail and implement -overwrite mode?

            #overwrite mode would remove existing link
            # as if source and target names match then new link would be made within the target
            #with the name of the source dir e.g. 'target/source',  not simply as 'target'
            next;
          }
        }

        if($link_type eq 'ftp'){
          #alter source dir here for old ftp mnt
          #after we have tested the source files
          #rather than before in the config
          $source_dir = $ftp_mnt.$source_dir;
        }
     
        my $cmd = "ln -sf $source_dir "; #$target_dir";
        system("$cmd") == 0 or die("Failed to create link using:\t${cmd}\nFrom $fclass_dir");

        #print "$cmd\n";
        $num_dirs{$link_type} ++;
      }
    }
  }

  foreach my $link_type (@link_types){
    print "Linked ".$num_dirs{$link_type}." $link_type $set_type directories\n";
  }

}
	


print "Source file paths:\n\t".join("\n\t", @true_paths)."\n";

	
#chdir($dir);


print "Need to add print of rsync command to send to webteam which will follow links and rsync to the ftp dir?\n";
