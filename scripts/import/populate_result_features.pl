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

=head1 DESCRIPTION

This script utilises the ResultFeature collector code to generate pre-computed windows from a given ResultSet and store them in the result_feature table.


=head1 USAGE

populate_result_features.pl [options]

Mandatory

Optional

=head2 SYNOPSIS


=cut

# InputSet/result_feature population
# Now we have two potantial ways of importing InputSets
# 1 Direct to feature_sets using the Importer.pm
# 2 Via this script if we add support for generating InputSets
# We need to incorporate/mirror some of this code into the Importer
# So we can take advantage of the Bed parser. Can we disentangle this from 
# the rest of the import parser code , so we can use it in isolation? (flat file output)
# 3 Add rollback level?


use strict;
use Data::Dumper;
use Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Funcgen::Utils::EFGUtils qw (strip_param_args strip_param_flags generate_slices_from_names);
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Pod::Usage;
use Getopt::Long;


my ($db_name, $db_host, $db_port, $db_user, $db_pass, $farm, $species);
my ($dnadb_name, $dnadb_host, $dnadb_port, $dnadb_user, $dnadb_pass, $force);
my ($rset_name, $old_assm, $new_assm, $dnadb, @slices, @skip_slices);#, $delete
my ($output_dir, $skip_zero_window, @window_sizes);

my @tmp_args = @ARGV;

print "populate_result_features.pl @tmp_args\n";


GetOptions(
		   #Mandatory
		   'db_host=s'        => \$db_host,
           'db_user=s'        => \$db_user,
           'db_port=i'        => \$db_port,
           'db_pass=s'        => \$db_pass,
           'db_name=s'        => \$db_name,
		   'result_set=s'     => \$rset_name,
		   'farm'             => \$farm,
		   'force'            => \$force,#force continue, stop jobs hanging on farm if RESULT_FEATURE_SET
           #'transcript_species_id=i' => \$transcript_species_id, #See probe2transcript.pl for how to implement multi species support
           #'transcript_multi_species' => \$transcript_multi_species,

		   #Optional
		   'dnadb_host=s'     => \$dnadb_host,
           'dnadb_user=s'     => \$dnadb_user,
		   'dnadb_port=i'     => \$dnadb_port,
           'dnadb_pass=s'     => \$dnadb_pass,
		   'dnadb_name=s'     => \$dnadb_name,
		   'output_dir=s'     => \$output_dir,
		   'old_assm=s'       => \$old_assm,
		   'new_assm=s'       => \$new_assm, 
		   'skip_slices=s{,}' => \@skip_slices,
		   'slices=s{,}'      => \@slices,
		   'skip_zero_window' => \$skip_zero_window,
		   'window_sizes=s{,}' => \@window_sizes,
		   #Do we need window_sizes here for projected windows?
	 	   'species=s'      => \$species,#Not really needed


		   #'delete'                 => \$delete,#How are we going to handle this?
		   #'debug'                  => \$debug,

		   #Helper params
		   #'tee'                    => \$main::_tee,#Should always be 1? as we run this on the farm and specify an -o and -e
		   #'filename'               => #\$main::_log_file,#Specified in the run script
		   
		   #add a reduced log to minimize memory usage?
           'help'                   => sub { pos2usage(-exitval => 0, -message => "Params are:\t@tmp_args"); }
		  ) or pod2usage(
						 -exitval => 1,
						 -message => "Params are:\t@tmp_args"
						);


if(@ARGV){
  die("You have trailing arguments which are unrecognised:\t@ARGV\n");
}

if(! $farm && ! @slices){
  die("To run with all slices you should specify the -farm option, or to run locally specify one -slice.\n");
}

if($new_assm && $skip_zero_window){
  die('You cannot specify both -new_assm and -skip_zero_window');
}

if(@window_sizes){
  warn 'You have specified custom window sizes, if these do no match those in Bio::EnsEMBL::Collector::ResultFeautre then this may lead to missing data, should be used for testing only!';
}

if($dnadb_name){
  $dnadb = Bio::EnsEMBL::DBSQL::DBAdaptor->new(
											   -dbname  => $dnadb_name,
											   -host    => $dnadb_host,
											   -port    => $dnadb_port,
											   -user    => $dnadb_user,
											   -pass    => $dnadb_pass,
											   -species => $species,
											   -type => 'core',
											  );
}


#This will add the default chromosome CS, but not any other levels
my $efg_db = Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor->new(
														  -dbname  => $db_name,
														  -host    => $db_host,
														  -port    => $db_port,
														  -user    => $db_user,
														  -pass    => $db_pass,
														  -species => $species,
														  -type    => 'funcgen',
														  -dnadb   => $dnadb,
														 );

#test connections(altho this should already be done in the run script)
$efg_db->dbc->db_handle;
$dnadb = $efg_db->dnadb;
$dnadb->dbc->db_handle;

my $rset_adaptor  = $efg_db->get_ResultSetAdaptor;
my $rfeat_adaptor = $efg_db->get_ResultFeatureAdaptor;
my $slice_adaptor = $efg_db->get_SliceAdaptor;
@slices = @{&generate_slices_from_names($slice_adaptor, \@slices, \@skip_slices, 'toplevel', 1, 1, $old_assm)};#non ref, inc dups
#non_dup here will load full chr collections for all hap/par regions!
#Need to implement fetch_normalised slice_projections in collection feature adaptors?
#This will fail for old assembly as we can have toplevel

 
#grab result_set
#should provide some default options which do this automatically for every rset which is displayable, but doesn't have RESULT_FEATURE_SET status.
my ($rset) = @{$rset_adaptor->fetch_all_by_name($rset_name)};
die("Could not fetch ResultSet using name:\t$rset_name") if ! $rset;


#Sanity check that we have a RESULT_FEATURE_SET
#If we are skipping_zero_level
#Do this here instead of in Collector to avoid submitting doomed jobs
#Move all this to the Collector?


if($skip_zero_window && ! $rset->has_status('RESULT_FEATURE_SET')){
  die("You are trying to -skip_zero_window for a ResultSet which is not yet a RESULT_FEATURE_SET.\n".
	  "It is likely you are trying to generate bins from assembly projected ResultFeatures.\n".
	  "Please check your 0 window size projected ResultFeatures and set the RESULT_FEATURE_SET status.\n".
	  "Or if you really want to omit the natural resolution, please do this in \@Bio::EnsEMBL::Utils::Collector::window_sizes\n");
  #Could do it here, but it is probably more safe like this
  #Altho new_assm automatically uses 0 window_size even if it is not set.
}
elsif($rset->has_status('RESULT_FEATURE_SET') &&
	  ! $skip_zero_window){
  die("You are trying to load ResultSet which is already a RESULT_FEATURE_SET, maybe you want to rollback_ResultFeatures or -skip_zero_window if you are trying to regenerate windows from a projected ResultSet");
  
}
elsif($skip_zero_window && $rset->has_status('RESULT_FEATURE_SET')){
  print "\nYou have specified -skip_zero_window for a RESULT_FEATURE_SET. This may duplicate data for window_sizes != 0\n";


  if(! $force){

	my $msg = "Continue?[y|n]"; 
	my $response;

	while (! defined $response){
	  print $msg;
	  $response = <STDIN>;
	  chomp $response;
	  exit if $response eq 'n';
	  $response = undef if $response ne 'y';
	}
	
	push @tmp_args, '-force';
	#$force = 1;
  }
}

#Set job parameters
if($farm){
  $output_dir ||= $ENV{'EFG_DATA'}."/output/${db_name}/result_features/${rset_name}";

  if(! -d $output_dir){
	system("mkdir -p $output_dir") == 0 || die("cannot make output directory:\t$output_dir\n$?");
  }

  @tmp_args = @{&strip_param_args(\@tmp_args, ('slices', 'skip_slices'))};
  @tmp_args = @{&strip_param_flags(\@tmp_args, ('farm'))}; 
}





SLICE: foreach my $slice(@slices){
  my $sr_name = $slice->seq_region_name;
    
  foreach my $done_slice(@skip_slices){
  
	if($sr_name eq $done_slice){
	  warn "Skipping slice:\t".$slice->name;
	  next SLICE;
	}
	
  }

  my $cmd = "perl $ENV{EFG_SRC}/scripts/import/populate_result_features.pl @tmp_args -slices ".$slice->name;
  
  if($farm){
    #hugemem -R "select[mem>20000] rusage[mem=20000] -M 20000000"

	my $bsub = "bsub -q long -J populate_result_features:${rset_name}:${sr_name} -o $output_dir/result_features.${sr_name}.out -e $output_dir/result_features.${sr_name}.err $cmd";
	
	print "Submitting $bsub\n";
	
	system($bsub) == 0 || die("Failed to submit $bsub\n$?");

	warn "Need to wait here and add RESULT_FEATURE_SET status if no errors found";
	#This is a pipeline thing, no?

  }
  else{

	my %wsize_conf = (@window_sizes) ? (-WINDOW_SIZES => \@window_sizes) : ();

	#Expose more config through the options?
	my %config = (-NEW_ASSEMBLY => $new_assm, 
				  -SKIP_ZERO_WINDOW => $skip_zero_window,
				  %wsize_conf
				 );

	foreach my $slice(@slices){
	  $rfeat_adaptor->store_window_bins_by_Slice_ResultSet($slice, $rset, %config);
	}
  }
}

if($farm){

  #Just get current assm for this warn
  if(! $new_assm){
	$new_assm = $slice_adaptor->fetch_by_region('chromosome', 1)->coord_system->version;
  }

  warn "You must manually check all jobs have run correctly before setting RESULT_FEATURE_SET, DAS_DISPLAYABLE, IMPORTED_${new_assm} states for ".$rset->name."\n" if ! $rset->has_status('RESULT_FEATURE_SET');

  if($new_assm){
	warn "To generate the remaning window sizes your must now resubmit the jobs with -skip_zero_window\n";
  }
}
