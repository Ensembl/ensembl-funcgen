#!/usr/bin/env perl

=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2017] EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either 405ress or implied.
See the License for the specific language governing permissions and
limitations under the License.

=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

  Questions may also be sent to the Ensembl help desk at
  <http://www.ensembl.org/Help/Contact>.

=head1 NAME

downloadad_inputset_data.pl

=head1 SYNOPSIS

download_inputset_data.pl -e experiment_name -s species -c cell_type -f feature_type  [-list_only]

=head1 DESCRIPTION

Downloads data for input sets registered in the data tracking database.

=head1 OPTIONS

=over

=item B<-help>

Gives this help message

=item B<-experiment_name|e>

Experiment name. (e.g. Mikkelsen07_PMID17603471)

=item B<-species|s>

Species for the dataset (e.g. homo_sapiens)

=item B<-cell_type|c>

Cell Type of the dataset (e.g. K562)

=item B<-f>

Feature Type for the dataset (e.g. Dnase1)

=item B<-list_only>

When specified, it does not download, but just lists the available datasets, data sources and their download status.

=item B<-force_download>

When specified, previous files will be replaced (e.g. in case a previous download was not successfull).

=item B<-allow_incomplete>

Downloads a dataset if even if not all urls for that dataset are available for download.
Otherwise only downloads if all urls are available (default behaviour).
Usually used in combination with -date and -override

=item B<-date>

When specified, it downloads all sets with availabiity date below the given date (defaults to current date)
Format: YYYY-MM-DD

=item B<-output_folder>

When specified, it downloads all data to this folder
By default uses the 'repository' meta entry in the tracking db

=item B<-ignore_date>

When specified, it ignores the availability date when downloading the set

=item B<-no_exit>

Skip over errors when downloading data

=item B<-tracking_dbhost>

Host where the data tracking database is (defaults to ens-genomics1/2 for mouse/human respectively)

=item B<-tracking_dbuser>

User of the data tracking database (defaults to ensadmin)

=item B<-tracking_dbpass>

Password for the data tracking database user

=item B<-tracking_dbbport>

Port of the host where the data tracking database

=item B<-tracking_dbname>

Name of the tracking database (defaults to "tracking_SPECIES_NAME_funcgen_SCHEMA_BUILD")

=back


=head1 SEE ALSO

ensembl-funcgen/scripts/tracking/register_new_inputset.pl




=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

  Questions may also be sent to the Ensembl help desk at
  <http://www.ensembl.org/Help/Contact>.



=cut


# TO DO
# We have add schema_build requirement here unless we re-instate dbname requirement
# Just changte this back to dbname reqiurement. won't this introduce a risk of using
# an old tracking DB?
# This needs to use InputSet/Subset and pass these to the TrackingAdaptor to update the DB.
# Update register_new_inputset first? No need as we already have input_sets loaded.
# Change this to list_only and exit by default?

use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use Bio::EnsEMBL::Funcgen::DBSQL::TrackingAdaptor;
use File::Path qw(make_path);

#use LWP::Simple;


my ($help, $experiment, $species, $cell_type, $feature_type, $list_only);
my ($allow_incomplete, $repository, $no_exit);
my $force_download = 0;
my $ignore_date = 0;
my ($tdb_name, $tdb_host, $tdb_pass, $tdb_port);
my ($dnadb_name, $dnadb_host, $dnadb_pass, $dnadb_port, $dnadb_user, $dnadb_assembly);
my $tdb_user = 'ensadmin';
my $date = '';

my @tmp_args = @ARGV;

my %species_hosts = (
					homo_sapiens => 'ens-genomics2',
					mus_musculus => 'ens-genomics1',	
				   );

#get command line options

print "download_dataset.pl @ARGV\n";

GetOptions 
  (
   'species|s=s'         => \$species,
   'experiment_name|e=s' => \$experiment,
   'cell_type|c=s'       => \$cell_type,
   'feature_type|f=s'    => \$feature_type,
   'list_only'           => \$list_only,
   'force_download'      => \$force_download,
   'ignore_date'         => \$ignore_date,
   'allow_incomplete'    => \$allow_incomplete,
   'date=s'              => \$date,
   'output_folder=s'     => \$repository,
   'no_exit'             => \$no_exit,

   #schema_build!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   #optional
   'tracking_dbhost=s'   => \$tdb_host,
   'tracking_dbport=s'   => \$tdb_port,
   'tracking_dbuser=s'   => \$tdb_user,
   'tracking_dbname=s'   => \$tdb_name,
   'tracking_dbpass=s'   => \$tdb_pass,
   'dnadb_host=s'        => \$dnadb_host, 
   'dnadb_user=s'        => \$dnadb_user,
   'dnadb_port=i'        => \$dnadb_port,
   'dnadb_pass=s'        => \$dnadb_pass,
   'dnadb_assembly=i'    => \$dnadb_assembly,
   'dnadb_name=s'        => \$dnadb_name,
   
   "help|h"              => \$help,
  )  or pod2usage( -exitval => 1 ); #Catch unknown opts

pod2usage(1) if ($help);


#Validate params

if($date && $ignore_date){
  die("You have chosen mutually exclusive params:\t-date(${date}) and -ignore_date");
}
elsif($ignore_date){
  $date = 'IGNORE';
}

=pod

if(! $species ||
   ($tdb_host && $tdb_name) ){
  
  die("You must either specify a:\n\t-species to use tracking DB defaults\nOR".
	  "\n\t-tracking_dbhost and -tracking_dbname");
}
elsif( ($species && $tdb_name) &&
	   ($tdb_name !~ /${species}/) ){
  die("You have specified mismatched -species ($species) and -tracking_dbname ($tdb_name)");
}
elsif($species && 
	  ! $tdb_host ){
  #Need to check we have default host for this species

  if(! exists $species_hosts{$species}){
	die("No defaults available for $species, please specify -tracking_dbhost or -species as one of:\t".
		join(', ', keys %species_hosts));
  }
}

=cut

# Connect to Tracking DB
# We don't want this to be available via $efgdb->get_TrackingAdaptor
# 
my $tracking_adaptor = Bio::EnsEMBL::Funcgen::DBSQL::TrackingAdaptor->new
  (
   -user     => $tdb_user,
   -pass     => $tdb_pass,
   -dbname   => $tdb_name,
   -host     => $tdb_host,
   -species  => $species,
   -dnadb_host     => $dnadb_host, 
   -dnadb_user     => $dnadb_user,
   -dnadb_port     => $dnadb_port,
   -dnadb_pass     => $dnadb_pass,
   -dnadb_assembly => $dnadb_assembly,
   -dnadb_name     => $dnadb_name,
  );


$species ||= $tracking_adaptor->db->species;

$repository ||= $tracking_adaptor->repository_path($repository);
#$repository .= '/'.$species;

if(! -d $repository){
  die("Repository does not seem to exist\t$repository");
}


#This needs to change to get InputSets or InputSubsets
#Change this to use standard set name format?
#Grab experiment, feature_type, cell_type first and pass to InputSetAdaptor method
#How are we going to filter downloaded is NULL
#This would require query extension from the input(sub)set adaptor to the tracking tables. NO!
#Filter in tracking adaptor
#DONE Add in param_hash processing to InputSet::fetch_all adaptor to enable complex queries with no mandatory args
#DONE Need to add support for new input_set/subset fields replicate (and control?)
#will defo need separate InputSubsetAdaptor once/if we decouple InputSet & InputSubset
##Do we still use the optional subset function? in InputSet::add_new_subset?
#look at analysis hack in InputSet::new
#Need to handle 255 replicate values(segmentation can be changed to 0?)
#DONE Update InputSubset description
#Can remove from??? Have we patched this?
#WTF? We have replicate, is_control duplciated between input_subset and input_subset_tracking???!!!
#This is becuase we have merged rep in input_subset, can get rid of this once we migrate true files to input_subset
#can't have ds.downloaded as status until this is done, as we can't use status methods on input_subset_tracking!
#md5?

#Validate experiment, ftype, ctype and build contraints hash

my %constraints;
my $efg_db        = $tracking_adaptor->db;
my $exp_adaptor   = $efg_db->get_ExperimentAdaptor;
my $ftype_adaptor = $efg_db->get_FeatureTypeAdaptor;
my $ctype_adaptor = $efg_db->get_CellTypeAdaptor;
#my $inp_set_adaptor = $efg_db->get_InputSetAdaptor;
my $inp_sset_adaptor = $efg_db->get_InputSubsetAdaptor;

#Could probably do this dynamically by changing 
#keys to class  name e.g. FeatureType rather than feature_type
#Then grab adaptors directly and call appropriate generic methods
#i.e. BaseAdaptor::define_contraints_by_names
#issue around plural names and generating adaptor names
#e.g. CellType vs CellTypes
#Maintain as singular and simply push single value onto array if not array ref.
#Think about this a little more
 
if($experiment){
  $constraints{experiments} = [$exp_adaptor->fetch_by_name($experiment)];
  
  if(! defined $constraints{experiments}->[0]){
	die("Could not fetch Experiment:\t$experiment");
  }
}
  
if($feature_type){
  $constraints{feature_types} = [$ftype_adaptor->fetch_by_name($feature_type)];
  
  if(! defined $constraints{feature_types}->[0]){
	die("Could not fetch FeatureType:\t$feature_type");
  }
}
  
if($cell_type){
  $constraints{cell_types} = [$ctype_adaptor->fetch_by_name($cell_type)];
  
  if(! defined $constraints{cell_types}->[0]){
	die("Could not fetch CellType:\t$cell_type");
  }
}

#my @input_sets = @{$inp_set_adaptor->fetch_all({constraints => \%constraints})};
my @input_ssets = @{$inp_sset_adaptor->fetch_all({constraints => \%constraints})};
my ($exp, $ctype, $ftype, $inp_set_id, $iss_id, $url, @inp_ssets, @iss_tracking_rows);
my ($rep, $target_dir, $target_file, $exp_dir);
my $action = ($list_only) ? 'Listing' : 'Downloading';

#my %inp_set_info;
my $tracking_info;
#foreach my $inp_set(@input_sets){
foreach my $inp_sset(@input_ssets){
  #$inp_set_id = $inp_set->dbID;
  #$exp        = $inp_set->get_Experiment;
  #$ctype      = $inp_set->cell_type->name;
  #$ftype      = $inp_set->feature_type->name;
  #$exp        = $inp_sset->experiment->archive_id;
  #$ctype      = $inp_sset->cell_type->name;
  #$ftype      = $inp_sset->feature_type->name;
  #$exp_dir    = $exp."/".$ctype."_".$ftype;
  
  #@iss_tracking_rows = @{$tracking_adaptor->fetch_InputSubset_tracking_info($inp_set, $force_download, $date)};
  $tracking_info = $tracking_adaptor->fetch_tracking_info($inp_sset, $force_download, $date);



  #Would it be better to set the trackign data in the InputSubset object?

  #Check to see if we have all the iss's for this input_set
  #else skip as we don't want to run with a subset of the files
  #This might be because the DOWNLOADED status was set incorrectly, 
  #or maybe some of the files have different availablility dates?
  
  #@inp_ssets = @{$inp_set->get_InputSubsets};
  @inp_ssets = ($inp_sset);

  #if(scalar(@inp_ssets) != scalar(@iss_tracking_rows)){
	#warn 'For InputSet '.$inp_set->name.' there is a mismatch between the number of registered InputSubsets('.
	#  scalar(@inp_ssets).') and the number of tracking info rows fetch('.
#		scalar(@iss_tracking_rows).") using -date(${date}) and -force_download(${force_download})\n".
#		  'Please set above params or correct database';
  if(! $tracking_info){
      warn 'Failed to get tracking information for InputSubSet '.$inp_sset->name.
       "\nUsing -date(${date}) and -force_download(${force_download})\n".
      'Please set above params or correct database';
	next;
  }

  $iss_id = $inp_sset->dbID;
  $url    = $tracking_info->{download_url};
  #$rep    = $inp_sset->replicate;

  #foreach my $iss_info_ref(@iss_tracking_rows){
	#($iss_id, $url, $rep) = @$iss_info_ref;#ignoring is_control and not_pooled
	#$url_replicate{$url} = $rep;
	#$expers_download{$inp_set_id}{$url} = $iss_id;

	#my ($dir, @pieces, $file);

	#foreach my $url (sort keys %{$expers_download{$ds_id}}){	  

	$target_dir = $repository."/".$inp_sset->experiment->archive_id;#.$rep;
	
  if((! $list_only) && ! -d $target_dir){

    if(! make_path($target_dir)){ 
      warn "Could not create directory:\t$target_dir";
      next; 
    }
  }
	
	($target_file = $url) =~ s'.*/'';
	$target_file = $target_dir.'/'.$target_file;

	#@pieces = split("/",$url);
	#$file = $dir."/". $pieces[-1];
	
	print "\t${action}: ".$target_file."\n";
	  
	if(! $list_only){
      #Save to a file... if not already saved (or if we want to override a file already saved)
      #Check insted in the data tracking database to see if file finished ok or not..
      #if(is_success(getstore($url,$file))){
      #Another option: lwp-download
      #Use SRA and structure...

      #BEFORE DOWNLOADING NEED TO RESET STATUS... IN CASE SOMETHING GOES WRONG!
	  $tracking_adaptor->set_download_status_by_input_subset_id($iss_id, 'to_null');
	  

    #warn "Running wget -nv $url -O $target_file";

    #TODO Use Perl-based download system... hopefully integrated with SRA...
    if(system("wget -nv $url -O $target_file")==0){
      #TODO Check md5sum...?
  		# WRITE LOG...
      $tracking_adaptor->set_download_status_by_input_subset_id($iss_id, undef, $target_file);		
    } 
	  else {
      my $error = "Problems downloading source:\t$url\n\t\t\t\tto:\t${target_file}";
		
		if($no_exit){
		  warn $error;
		}
		else{
		  die($error);
		}
      }
    }
  #}
}


1;
