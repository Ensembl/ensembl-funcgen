
=head1 DESCRIPTION

This script utilises the ResultFeature collector code to generate pre-computed windows from a given ResultSet and store them in the result_feature table.


=head1 USAGE

populate_result_features.pl [options]

Mandatory

Optional

=head2 SYNOPSIS

=head1 TO DO

Write run_populate_result_features.sh
Handle pre-release dnadbs
POD

=head1 CONTACT

Post comments or questions to the Ensembl development list: ensembl-dev@ebi.ac.uk

=head1 CVS

 $Log: not supported by cvs2svn $

=cut



use strict;
use Data::Dumper;
use Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Pod::Usage;
use Getopt::Long;


my ($db_name, $db_host, $db_port, $db_user, $db_pass, $farm, $species);
my ($dnadb_name, $dnadb_host, $dnadb_port, $dnadb_user, $dnadb_pass, $force);
my ($slice_name, $rset_name, $old_assm, $new_assm, $dnadb, @slices, @skip_slices);#, $delete
my ($output_dir, $skip_zero_window, @window_sizes);

my @tmp_args = @ARGV;

GetOptions(
		   #Mandatory
		   'db_host=s'        => \$db_host,
           'db_user=s'        => \$db_user,
           'db_port=i'        => \$db_port,
           'db_pass=s'        => \$db_pass,
           'db_name=s'        => \$db_name,
		   'result_set=s'     => \$rset_name,
		   'slice_name=s'     => \$slice_name,
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
		   #'probe_species_id=i' => \$probe_species_id,
		   #'probe_multi_species' => \$probe_multi_species,
		   'old_assm=s'       => \$old_assm,
		   'new_assm=s'       => \$new_assm, 
		   'skip_slices=s{,}' => \@skip_slices,
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

if(! $farm && ! $slice_name){
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
														  -type => 'funcgen',
														  -dnadb   => $dnadb,
														 );

#test connections(altho this should already be done in the run script)
$efg_db->dbc->db_handle;
$dnadb = $efg_db->dnadb;
$dnadb->dbc->db_handle;

my $rset_adaptor  = $efg_db->get_ResultSetAdaptor;
my $rfeat_adaptor = $efg_db->get_ResultFeatureAdaptor;
my $slice_adaptor = $efg_db->get_SliceAdaptor;

#generate slice/s
if($slice_name){
  @slices = ($slice_adaptor->fetch_by_name($slice_name));
  die("Could not generate valid slice from slice name:\t$slice_name") if ! defined $slices[0];#This should already be done in the runner
}
else{
  #We are not handling top level here!
  #This is okay as unassembled seq_regions are the same between assemblies
  #And should inherit the old features
  #This may not be quite true if we have different version of supercontigs
  #but these *should* be named differently in the DB
  #Will not catch contig/supercontig > chromosome assembly projection
  #But the assembly projection does not yet support this anyway
  @slices = @{$slice_adaptor->fetch_all('chromosome', $old_assm)};
}

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
  $output_dir ||= $ENV{'EFG_DATA'}.'/result_features/'.$db_name;

  if(! -d $output_dir){
	system("mkdir -p $output_dir") == 0 || die("cannot make output directory:\t$output_dir\n$?");
  }

  #We basically want to strip out farm, slice_name and skip_slices
  #This would be easier with something hash based
  my $seen_opt = 0;

  foreach my $i(0..$#tmp_args){
	$seen_opt = 0 if $tmp_args[$i] =~ /^[-]+/;#Reset seen opt if we seen a new one

	$tmp_args[$i] = '' if $tmp_args[$i] =~ /^[-]+farm/;#Only remove current flag
	$seen_opt = 1 if $tmp_args[$i] =~ /^[-]+skip_slices/;
	$seen_opt = 1 if $tmp_args[$i] =~ /^[-]+slice/;#Don't have full param name incase we have just specified -slice
	
	$tmp_args[$i] = '' if $seen_opt;#Remove args following option
  }

}


SLICE: foreach my $slice(@slices){
  my $sr_name = $slice->seq_region_name;
    
  foreach my $done_slice(@skip_slices){
  
	if($sr_name eq $done_slice){
	  warn "Skipping slice:\t".$slice->name;
	  next SLICE;
	}
	
  }

  my $cmd = "perl $ENV{EFG_SRC}/scripts/import/populate_result_features.pl @tmp_args -slice_name ".$slice->name;
  
  if($farm){
    #hugemem -R "select[mem>20000] rusage[mem=20000] -M 20000000"

	my $bsub = "bsub -q long -o $output_dir/result_features.${sr_name}.out -e $output_dir/result_features.${sr_name}.err $cmd";
	
	print "Submitting $bsub\n";
	
	system($bsub) == 0 || die("Failed to submit $bsub\n$?");

	warn "Need to wait here and add RESULT_FEATURE_SET status if no errors found";
	#This is a pipeline thing, no?
	warn "Need to kick of 2nd stage of project bin generation";
	#We need to manually check output and 

  }
  else{

	#Expose more config through the options?
	my %config = (-NEW_ASSEMBLY => $new_assm, 
				  -SKIP_ZERO_WINDOW => $skip_zero_window,
				  -WINDOW_SIZES => \@window_sizes
				 );

	foreach my $slice(@slices){
	  $rfeat_adaptor->store_window_bins_by_Slice_ResultSet($slice, $rset, %config);
	}
  }
}

if($farm){
  warn "You must manually check all jobs have run correctly before setting RESULT_FEATURE_SET status for ".$rset->name."\n";

  if($new_assm){
	warn "To generate the remaning window sizes your must now resubmit the jobs with -skip_zero_window\n";
  }
}
