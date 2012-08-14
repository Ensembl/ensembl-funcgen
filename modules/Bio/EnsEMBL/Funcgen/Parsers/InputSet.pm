#
# EnsEMBL module for Bio::EnsEMBL::Funcgen::Parsers::InputSet
#

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

Bio::EnsEMBL::Funcgen::Parsers::InputSet

=head1 SYNOPSIS

  use vars qw(@ISA);
  @ISA = qw(Bio::EnsEMBL::Funcgen::Parsers::InputSet);

  
=head1 DESCRIPTION

This is a base class to support simple file format parsers. For simple imports the vendor is
set to the parser type i.e. the file format.  The generic read_and_import_simple_data assumes
a one line per feature format, other format need there own read_and_import_format_data method, 
which will need defining in the result_data config element. Features are stored either as 
ResultFeature collections or AnnotatedFeatures dependan ton the input feature class.

=cut

# To do
# Add Parsers for BAM/SAM
# Rename to InputSet
# Handle mysqlimport for large data sets e.g. reads
# Incorporate collection code
# Implement matrix storage

package Bio::EnsEMBL::Funcgen::Parsers::InputSet;

use Bio::EnsEMBL::Funcgen::AnnotatedFeature;
use Bio::EnsEMBL::Funcgen::SegmentationFeature;
use Bio::EnsEMBL::Utils::Exception qw( throw warning deprecate );
use Bio::EnsEMBL::Funcgen::Utils::EFGUtils qw(species_chr_num open_file is_gzipped);
use Bio::EnsEMBL::Utils::Argument qw( rearrange );
#use Bio::EnsEMBL::Funcgen::Utils::Helper;
use strict;

#config stuff, move to BaseImporter?
use Bio::EnsEMBL::Analysis;
use Bio::EnsEMBL::Funcgen::FeatureType;


use base qw(Bio::EnsEMBL::Funcgen::Parsers::BaseImporter); #@ISA change to parent with perl 5.10

#use vars qw(@ISA);
#@ISA = qw(Bio::EnsEMBL::Funcgen::Utils::Helper);

my %valid_types = (
				   result       => undef,
				   annotated    => undef,
				   segmentation => undef,
				  );


=head2 new

  Example    : my $self = $class->SUPER::new(@_);
  Description: Constructor method for Bed class
  Returntype : Bio::EnsEMBL::Funcgen::Parsers::Simple
  Exceptions : throws if caller is not Importer
  Caller     : Bio::EnsEMBL::Funcgen::Parsers:Simple
  Status     : at risk

=cut


sub new{
  my $caller = shift;
  
  my $class = ref($caller) || $caller;
  my $self  = $class->SUPER::new(@_);

  ($self->{'input_set_name'}, 
   $self->{'input_feature_class'}, 
   $self->{total_features}, 
   $self->{force},                  #Is this generic enough to go in Importer? used by store_window_bins_by_Slice_Parser
   $self->{dbfile_data_root},       #only appropriate for result input_feature_class
   $self->{merged_replicates},
  ) = rearrange
    (
     ['input_set_name', 'input_feature_class', 'total_features', 'force',
      'dbfile_data_root', 'merged_replicates'],
     @_);


  #Could potentially take fields params directly to define a custom format
  #Take direct field mappings, plus special fields which needs parsing differently
  #i.e. default is tab delimited, and GFF would define Attrs field as compound field and provide special parsing and field mapping
  

  throw("This is a skeleton class for Bio::EnsEMBL::Importer, should not be used directly") 
	if(! $self->isa("Bio::EnsEMBL::Funcgen::Importer"));
  
  $self->{'config'} =  
	{(
	  #can we omit these?
      array_data   => [],#['experiment'],
      probe_data   => [],#["probe"],
      norm_method => undef,
	  #protocols => {()},
	  'results_data' => ["and_import"], 
     )};

  #set up feature params
  $self->{'_feature_params'} = {};
  $self->{'_dbentry_params'} = [];
   
  return $self;
}


sub output_file{
  my ($self, $output_file) = @_;

  $self->{'output_file'} = $output_file if $output_file;
  return $self->{'output_file'};
}

sub input_file{
  my ($self, $input_file) = @_;

  $self->{'input_file'} = $input_file if $input_file;
  return $self->{'input_file'};
}



=head2 set_config

  Example    : my $self->set_config;
  Description: Sets attribute dependent config
  Returntype : None
  Exceptions : None
  Caller     : Bio::EnsEMBL::Funcgen::Importer
  Status     : at risk

=cut


sub set_config{
  my $self = shift;

  #Move all this to new when we fix the inheritance in Importer

  #We could set input_set_name to experiment name
  #But we would have to make warning in define_and_validate_sets mention -input_set_name

  throw('Must provide an -input_set name for a '.uc($self->vendor).' import') if ! defined $self->input_set_name();

  #Mandatory checks
  if(! defined $self->feature_analysis){
	throw('Must define a -feature_analysis parameter for '.uc($self->vendor).' imports');
  }

  #some convenience methods
  $self->{dbentry_adaptor}   = $self->db->get_DBEntryAdaptor;
  $self->{input_set_adaptor} = $self->db->get_InputSetAdaptor;
  
  my $fclass_name = $self->{input_set_adaptor}->build_feature_class_name($self->input_feature_class);
  $self->{'feature_class'} = 'Bio::EnsEMBL::Funcgen::'.$fclass_name;

  #We need to undef norm method as it has been set to the env var
  $self->{config}{norm_method} = undef;

  #dirs are not set in config to enable generic get_dir method access
  $self->{dbfile_data_root} ||= $self->get_dir('output');#Required for Collector

  my $adaptor_method         = 'get_'.$fclass_name.'Adaptor';

  #cannot validate here using 'can' as it doesn't appear to work for autoloaded methods :(
  eval { $self->{feature_adaptor} =  $self->db->$adaptor_method; };

  if($@){
    throw("Failed to $adaptor_method, ".$self->input_feature_class.' is not a valid feature class');
  }


  #Validate slices
  $self->slices($self->{'slices'}) if defined $self->{'slices'};

  #Move to new when we sort out inheritance
  $self->validate_and_store_config([$self->name]);
  #Could use input_set_name here?
  #This was to support >1 input set per experiment (name)

  #This current breaks for no config imports
  #i.e. standard Bed import e.g. result_feature collections
  #segmentation imports use Bed and config
  #allow no config imports in BaseImporter?
  #or ultimately set the params as part of the user_config?

  return;
}


#Some of this can probably move to define_and_validate_sets?

sub define_sets{
  my ($self) = @_;

  my $eset = $self->db->get_InputSetAdaptor->fetch_by_name($self->input_set_name);
  
  if(! defined $eset){

    #temporary fix until we implement input_subset_file replicate/control config
    #and merge validite_files in here
    #do not currently have non merged replicate import
    #need to define this param in all the run scripts
    my $replicate = 0 if $self->merged_replicates;

    $eset = Bio::EnsEMBL::Funcgen::InputSet->new
      (
       -name          => $self->input_set_name(),
       -experiment    => $self->experiment(),
       -feature_type  => $self->feature_type(),
       -cell_type     => $self->cell_type(),
       -vendor        => $self->vendor(),
       -format        => $self->format(),
       -analysis      => $self->feature_analysis,
       -feature_class => $self->input_feature_class,
       -replicate     => $replicate,
      );


    
    ($eset)  = @{$self->db->get_InputSetAdaptor->store($eset)};
  }

  #Use define_and_validate with fetch/append as we may have a pre-existing set
  #This now needs to handle ResultSets based on InputSets

  my $dset = $self->define_and_validate_sets
	(
	 -dbadaptor    => $self->db,
	 -name         => $self->input_set_name,#or name?
	 -feature_type => $self->feature_type,
	 -cell_type    => $self->cell_type,
	 -analysis     => $self->feature_analysis,
	 -feature_class=> $self->input_feature_class, 
	 -description  => $self->feature_set_description,
	 #-append          => 1,#Omit append to ensure we only have this eset
	 -recovery     => $self->recovery,
	 -supporting_sets => [$eset],
	 -slices        => $self->slices,
	 #Can't set rollback here, as we don't know until we've validated the files
	 #Can't validate the files until we have the sets.
	 #So we're doing this manually in validate_files
	);

  #We are now using IMPORTED to define wheather a FeatureSet has been imported succesfully
  #However we already have IMPORTED on the InputSubSet
  #We should add it to FeatureSet to remain consistent.
  #See Helper::define_and_validate_sets for more notes on
  #potential problems with FeatureSet IMPORTED status
 

  #define_and_validate_sets should also replicate ResultSets?
  #Questionable, mapped reads are never normalised across replicates
  #There are generally used as input for peak calling individually.
  #So files in this instance are expected to be separate parts of the same replicate
  #e.g. different chromosomes
  #Force one input file?
  #What if we want to link several assays(feature/cell_types) to the same experiment?

  $self->{'_data_set'} = $dset;

  

  #set output_set and new_data as attrs, or return?
  #change define_sets to define_sets_and_validate_input_files
  #eset is the supporting set for this import
  my $output_set;


  if ($self->input_feature_class eq 'result' ||
      $self->input_feature_class eq 'dna_methylation') {
    
    $output_set = $dset->get_supporting_sets->[0];
    $eset       = $output_set->get_InputSets->[0];

    $self->result_set($output_set); #required for ResultFeature Collector and Bed Parser
  } 
  else {                      #annotated/segmentation
    $output_set = $dset->product_FeatureSet;
    $eset       = $dset->get_supporting_sets->[0];
  }

  return ($dset, $output_set, $eset);
}




#we have rollback functionality incorporated here
#why is this not integrated into define_sets?
#they are inherantly linked and for an InputSet import
#we always want to validate files.
#No way of defining replicate/control config just yet
#This should be changed to support both tracking DB and direct import.
#to correctly define the sets, we need access to the files and in the
#rep/control config
#We can do this implicitly by simply letting the order define the rep number
#-result_files rep1 rep2 rep2
#-control_files rep1 rep2
#-merged_replicates would over-ride the above and trigger a check for >1 file
#let's implement -merged_replicates for now and come back to this order based config
#rename -result_files to -input_subset_files?
#Will this support the collection import?
#Need to be mindful about the true input files and what pipeline will actually use
#as it may be some merged intermediate i.e. a sam/bed file rather than a fastq


#should we also store the dbfile_registry data here?
#this maybe done before generating the actual files


sub validate_files{
  my ($self, $prepare) = @_;

  #Get file

  if (! @{$self->result_files()}) {
    my $list = "ls ".$self->get_dir('input').'/'.$self->input_set_name().'*.'; #.lc($self->vendor);#could use vendor here? Actually need suffix attr
    my @rfiles = `$list`;
    $self->result_files(\@rfiles);
  }

  #We don't yet support multiple files
  if (scalar(@{$self->result_files()}) >1) {
    warn('Found more than one '.$self->vendor." file:\n".
         join("\n", @{$self->result_files()})."\nThe InputSet parser does not yet handle multiple input files(e.g. replicates).".
         "  We need to resolve how we are going handle replicates with random cluster IDs");
    #do we even need to?
  }
  
  #Here were are tracking the import of individual files by adding them as InputSubSets
  #Recovery would never know what to delete
  #So would need to delete all, Hence no point in setting status?
  #We do not rollback IMPORTED data here.  This is done via separate scripts
  #To reduce the rick of accidentally deleting/overwriting data by leaving a stry -rollback
  #flag in the run script

  ### VALIDATE FILES ###
  #We need validate all the files first, so the import doesn't fall over half way through
  #Or if we come across a rollback halfway through
  my (%new_data, $eset);
  my $dset = $self->data_set;
   
  if ((scalar(@{$self->slices}) > 1) &&
      ! $prepare) {
    throw('Validate files does not yet support multiple Slice rollback');
  }
  
  #This all assumes that there is only ever 1 InputSet
  
  if ($self->input_feature_class eq 'result' ||
      $self->input_feature_class eq 'dna_methylation' ) {
    $eset = $dset->get_supporting_sets->[0]->get_InputSets->[0];
  } else {                      #annotated/segmentation
    $eset =  $dset->get_supporting_sets->[0];
  }
  
  
  #IMPORTED status here may prevent
  #futher slice based imports
  #so we have wait to set this until we know all the slices 
  #are loaded, unless we store slice based IMPORTED states
  #We currently get around this be never settign IMPORTED for slice based jobs
  #and always rolling back by slice before import

  #This loop supports multiple files
  my (@rollback_sets, %file_paths);
  my $auto_rollback = ($self->rollback) ? 0 : 1;

  
  #set the max seen replicate value
  my $replicate = 0;

  foreach my $iss( @{$eset->get_InputSubsets} ){
    
    if($iss->replicate >  $replicate){
      $replicate = $iss->replicate;
    }
  }


  foreach my $filepath ( @{$self->result_files} ) {
    my ($filename, $sub_set);
    chomp $filepath;
   	($filename = $filepath) =~ s/.*\///;
    $file_paths{$filename} = $filepath;			 
    $filename =~ s/^prepared\.// if $self->prepared; #reset filename to that originally used to create the Inputsubsets

    $self->log('Validating '.$self->vendor." file:\t$filename");
    throw("Cannot find ".$self->vendor." file:\t$filepath") if(! -e $filepath); #Can deal with links
		
    if ( $sub_set = $eset->get_subset_by_name($filename) ) {
      #IMPORTED status here is just for the file
      #Any changes to analysis or coord_system should result in different InputSubset(file)
      #Will only ever be imported into one Feature|ResultSet

      #Currently conflating recover_unimported and rollback
      #as they serve the same purpose until we implement InputSubset level recovery

	
      if ( $sub_set->has_status('IMPORTED') ) {
        $new_data{$filepath} = 0;
        $self->log("Found previously IMPORTED InputSubset:\t$filename");
      } else {
        $self->log("Found existing InputSubset without IMPORTED status:\t$filename");
        push @rollback_sets, $sub_set;
      }
    } else {
      $replicate ++;

      $self->log("Found new InputSubset:\t${filename}");
      throw("Should not have found new 'prepared' file:\t$filename") if $self->prepared;	  
      $new_data{$filepath} = 1;

      $replicate = 0 if $self->merged_replicates;
      #need to add is_control support here too.
      #tidy this up, force full creation here and remove convenience code from InputSet?
      #where else is it used?

      $sub_set = Bio::EnsEMBL::Funcgen::InputSubset->new
        (
         -name       => $filename,
         -replicate  => $replicate,
         -input_set  => $eset,
         # should also handle these here
         -is_control  => 0,#for result files, need to handle control files separately
         #to be defined
         #-display_url => $display_url,
         #-archive_id  => $archive_id,
        );

      $self->input_set_adaptor->store_InputSubsets([$sub_set]);
    }
  }


  #Does -recover allow a single extra new file to be added to an existing InputSet?
 
  if (@rollback_sets &&         #recoverable sets i.e. exists but not IMPORTED
      ( (! $self->recovery) && (! $self->rollback) ) ) {
    throw("Found partially imported InputSubsets:\n\t".join("\n\t", (map $_->name, @rollback_sets))."\n".
          "You must specify -recover or -rollback to perform a full rollback");

    if ($self->recovery) {
      #Change these to logger->warn
      $self->log("WARNING::\tCannot yet rollback for just an InputSubset, rolling back entire set? Unless slices defined");
      $self->log("WARNING::\tThis may be deleting previously imported data which you are not re-importing..list?!!!\n");
    }
  }
  
  
  if ($self->rollback) {
    #Check we have all existing InputSubsets files before we do full rollback
    #Can probably remove this if we support InputSubset(file/slice) level rollback
    $self->log('Rolling back all InputSubsets');
    @rollback_sets = @{$eset->get_InputSubsets};

    foreach my $isset (@rollback_sets) {
	  
      if (! exists $file_paths{$isset->name}) {
        throw("You are attempting a multiple InputSubset rollback without specifying an existing InputSubset:\t".$isset->name.
              "\nAborting rollback as data will be lost. Please specifying all existing InputSubset file names");
      }
    }
  }


  foreach my $esset (@rollback_sets) {
    #This needs to be mapped to the specified filepaths
    my $fp_key = $esset->name;
    $fp_key = 'prepared.'.$fp_key if $self->prepared;

    $new_data{$file_paths{$fp_key}} = 1;
    $self->log("Revoking states for InputSubset:\t\t\t".$esset->name);
    $eset->adaptor->revoke_states($esset);

    if (! $prepare) {
      #This was to avoid redundant rollback in prepare step
      $self->log("Rolling back InputSubset:\t\t\t\t".$esset->name);
	  
      if ($self->input_feature_class eq 'result' ||
          $self->input_feature_class eq 'dna_methylation') {
        #Can we do this by slice for parallelisation?
        #This will only ever be a single ResultSet due to Helper::define_and_validate_sets
        #flags are rollback_results and force(as this won't be a direct input to the product feature set)
        $self->rollback_ResultSet($self->data_set->get_supporting_sets->[0], 1, $self->slices->[0], 1);
        #Do no have rollback_InputSet here as we may have parallel Slice based imports running
      } else {                  #annotated/segmentation
        $self->rollback_FeatureSet($self->data_set->product_FeatureSet, undef, $self->slices->[0]);
        $self->rollback_InputSet($eset);
        last;
      }			
    }
  }

  return \%new_data;
}




sub set_feature_separator{
  my ($self, $separator) = @_;

  #How do we test if something undefined was passed?
  #Rather than nothing passed at all?
  #Can't do this as this is the accessor
  #Need to split method
 
  throw('Must provide a valid feature separator') if ( (! defined $separator) || ($separator eq '') ); 

  $self->{'_feature_separator'} = $separator;

}

### SIMPLE ACCESSORS ###
# Some of these can be called for each record
# Trim the access time as much as possible

#This is to support input_set/input_subset.replicate
#Not appropriate for array experiment imports

sub merged_replicates{   return $_[0]->{merged_replicates}; }

sub input_feature_class{ return $_[0]->{input_feature_class}; }

sub input_set_name{      return $_[0]->{input_set_name}; }  #Set in new

sub feature_adaptor{     return $_[0]->{feature_adaptor}; }

sub dbentry_adaptor{     return $_[0]->{dbentry_adaptor}; }

sub input_set_adaptor{   return $_[0]->{input_set_adaptor}; }

sub set{                 return $_[0]->{set}; } #Feature or Result, set in define_sets

sub data_set{            return $_[0]->{_data_set}; }

sub feature_separator{   return $_[0]->{_feature_separator}; }

sub feature_params{      return $_[0]->{_feature_params}; }

sub dbentry_params{      return $_[0]->{_dbentry_params}; }

sub input_gzipped{       return $_[0]->{input_gzipped}; }


sub input_file_operator{
  my ($self, $op) = @_;
  #Should be set in format parser
  $self->{'input_file_operator'} = $op if defined $op;

  return $self->{'input_file_operator'};
}

# prepare boolean, simply stores the sets and preprocesses the file
# so we don't get each batch job trying to sort etc


#Still need to implement prepare in other Parsers!!

sub read_and_import_data{
  my ($self, $prepare) = @_;

  my $action = ($prepare) ? 'preparing' : 'importing';

  $self->log("Reading and $action ".$self->vendor()." data");
  my ($filename, $fh, $f_out, %feature_params, @lines);

  if($prepare && ! $self->isa('Bio::EnsEMBL::Funcgen::Parsers::Bed')){
	throw('prepare mode is only currently implemented for the Bed parser');
  }
  
 
  #Test for conflicting run modes
  if($prepare &&
	 ($self->batch_job || $self->prepared)){
	#prepare should be called once by the runner, not in each batch_job
	#don't prepare if already prepared
	throw('You cannot run read_and_import_data in prepare mode with a -batch_job or -prepared job');
  }

  my ($dset, $output_set, $eset) = $self->define_sets;
  
  #We also need to account for bsub'd slice based import

  
  #If we can do these the other way around we can get define_sets to rollback the FeatureSet
  #Cyclical dependency for the sets :|
  my $new_data = $self->validate_files($prepare);
  my $seen_new_data = 0;
 
  
  ### READ AND IMPORT FILES ###
  foreach my $filepath(@{$self->result_files()}) {
	chomp $filepath;
	
	($filename = $filepath) =~ s/.*\///;
	$self->input_file($filepath); #This is only used by Collector::ResultFeature::reinitialise_input method
	
	if($new_data->{$filepath} ){	#This will currently autovivify!
	  $seen_new_data = 1;
	  $self->{'input_gzipped'} = &is_gzipped($filepath);
	  
	  $filepath = $self->pre_process_file($filepath, $prepare) if $self->can('pre_process_file');
	  
	  $self->log_header(ucfirst($action).' '.$self->vendor." file:\t".$filepath);

	  #We need to be able to optional open pipe to gzip | sort here
	  #i.e. define open command
	  $fh = open_file($filepath, $self->input_file_operator);

	  #This my become way too large for some reads files
	  #Currently no problems
	  #This is not working as we are sorting the file!
	  #$self->parse_header($fh) if $self->can('parse_header');
	  
	  #For result features some times we want to run 
	  #locally and just sort without dumping
	  #i.e if we are not a batch job
	  #as there is no need to dump if it is a single process
	  

	  #Should this be prepared?
    
	  if((($self->input_feature_class eq 'result' || $self->input_feature_class eq 'dna_methylation' ) 
        && ! $prepare)){
      
      #Use the ResultFeature Collector here
      #Omiting the 0 wsize
      #How are we going to omit 0 wsize when doing the fetch?
      #simply check table name in ResultSet?
      
      #Should we do this for multiple chrs?
      #or fail here
      # we need to pass self
      #for access to get_Features_by_Slice
		#which should be in the specific parser e.g Bed
      
		#Will this not clash with standard ResultFeature::get_ResultFeatures_by_Slice?
		#Could really do with separating the pure file parsers from the importer code, so these can be reused
		#by other code. Then simply use Bed import parser for specific import functions and as wrapper to 
		#Bed file parser
		#So should really have
		#Parsers::File::Bed
		#and
		#Parsers::Import::Bed
		#This can probably wait until we update BioPerl and just grab the Bed parser from there?

		my $slices = $self->slices;

		#Should this be caught in new?
		if(! @$slices){
		  throw("You must define a slice to generate ResultFeature Collections from InputSet:\t".$eset->name);
		}
		

		if(scalar(@$slices) > 1){
		  throw("InputSet parser does not yet support multi-Slice import for ResultFeature collections\n"
				."Please submit these to the farm as single slice jobs");
		}

		#restrict to just 1 slice as we don't yet support disk seeking
		#if the slices are not in the same order as they appear in the file
		#also we want to parallelise this
		
		#Set as attr for parse_Features_by_Slice in format sepcific Parsers

	
		$self->file_handle(open_file($filepath, $self->input_file_operator));
		
		
		foreach my $slice(@$slices){
		  $self->feature_adaptor->store_window_bins_by_Slice_Parser($slice, $self, 
																	(
																	 #Force needs reimplementing here?
																	 -force            => $self->{force},
																	 -dbfile_data_root => $self->{dbfile_data_root},
																	));
		}  

		warn "Need to update InputSubset status to IMPORTED after all slices have been loaded";
		#Do we even need to set RESULT_FEATURE_SET for input_set ResultFeatures?

		

		warn "Closing $filename\nDisregard the following 'Broken pipe' warning";

		#Closing the read end of a pipe before the process writing to it at the other end 
		#is done writing results in the writer receiving a SIGPIPE. If the other end can't 
		#handle that, be sure to read all the data before closing the pipe.
		#This suggests the gzip pipe has not finished reading, but it *should* be at the end of the file?
		#$SIG{PIPE} = 'IGNORE'; #Catch with subref and warn instead?
		#Or maybe an eval will catch the STDERR better?
		#sub catch_SIGPIPE{
		#  my $sig = shift @_;
		#  print " Caught SIGPIPE: $sig $1 \n";
		#  return;
		#  
		#}
		#$SIG{PIPE} = \&catch_SIGPIPE;
		#my $retval =  eval { no warnings 'all'; $fh->close };
		#if($@){
		#  warn "after eval with error $@\nretval is $retval";
		#}
		#Neither of these catch gzip: stdout: Broken pipe
		
		#IO::UnCompress::Gunzip?


		$fh->close;
	  }
	  else{
		  
	  
		#Revoke FeatureSet IMPORTED state here incase we fail halfway through
		$output_set->adaptor->revoke_status('IMPORTED', $output_set) if ($output_set->has_status('IMPORTED') && (! $prepare));

		#What about IMPORTED_"CSVERSION"
		#This may leave us with an incomplete import which still has
		#an IMPORTED_CSVERSION state
		#We need to depend on IMPORTED for completeness of set
		#DAS currently only uses IMPORTED_CSVERSION
		#This is okayish but we also need to write HCs for any sets 
		#which do not have IMPORTED state!
		my ($line, @outlines, $out_fh);

		
		if($prepare && ! $self->batch_job){
		  #Assume we want gzipped output
		  #filename is actull based on input, so may not have gz in file name
		  $out_fh = open_file($self->output_file, "| gzip -c > %s");
		}
		

		while(defined ($line=<$fh>)){
		  #Generic line processing
		  #Move these to parse_line?
		  $line =~ s/\r*\n//o;
		  next if $line =~ /^\#/;	
		  next if $line =~ /^$/;

		  #This has now been simplified to process_line method
		  #my @fields = split/\t/o, $line;
		  #start building parameters hash
		  #foreach my $field_index(@{$self->get_field_indices}){
		  #  my $field = $self->get_field_by_index($field_index);
		  #  $feature_params = ($field =~ /^-/) ? $fields[$field_index] : $self->$field($fields[$field_index]);
		  #  }	


		  #We also need to enable different parse line methods if we have different file
		  #e.g. cisRED
		  #Code refs?

		  
		  if($self->parse_line($line, $prepare)){
			$self->count('total parsed lines');

			#Cache or print to sorted file
			if($prepare && ! $self->batch_job){

			  if(scalar(@outlines) >1000){
				print $out_fh join("\n", @outlines)."\n";
				@outlines = ();
			  }
			  else{
				push @outlines, $line;
			  }
			}
		  }
		}
	 
		close($fh);

		#Print last of sorted file
		if($prepare && ! $self->batch_job){
		  print $out_fh join("\n", @outlines)."\n";
		  close($out_fh);
		  @outlines = ();
		}

		if(! $prepare){
		  #Now we need to deal with anything left in the read cache
		  $self->process_params if $self->can('process_params');
	  
		  #To speed things up we may need to also do file based import here with WRITE lock?
		  #mysqlimport will write lock the table by default?
		 		  
		  #reset filename to that originally used to create the Inputsubsets
		  $filename =~ s/^prepared\.// if $self->prepared;

		  my $sub_set = $eset->get_subset_by_name($filename);
		  $sub_set->adaptor->store_status('IMPORTED', $sub_set) if ! $self->batch_job;
		}
	  }

	  
	  if($prepare){
		$self->log("Finished preparing import from:\t$filepath");
	  }
	  else{
		#Need to tweak this for slice based import
		$self->log('Finished importing '.$self->counts('features').' '.
				   $output_set->name." features from:\t$filepath");
		
	  }
	
	  
	  #This currently fails here if the uncaught file sort was not successful

	  foreach my $key (keys %{$self->counts}){
		$self->log("Count $key:\t".$self->counts($key));
	  }	  
	}
  }

  #Here we should set IMPORTED on the FeatureSet
  #We could also log the first dbID of each feature in a subset to facilitate subset rollback
  #in feature table
  #this would be sketchy at best
  #delete from annotated_feature where annotated_feature_id >= $first_subset_feature_id and feature_set_id=$feature_set_id
  #This may already have IMPORTED status as we don't revoke the status whilst
  #updating to protect the feature set due to lack of supportingset tracking
  #see Helper::defined_and_validate_sets for more notes.
  #Is there any point in setting it if we don't revoke it?
  #To allow consistent status handling across sets. Just need to be aware of fset status caveat.
  #Also currently happens with ResultFeatures loaded by slice jobs, as this may already be set by a parallel job

  if(! $prepare){
	$output_set->adaptor->set_imported_states_by_Set($output_set) if $seen_new_data && ! $self->batch_job;
	$self->log("No new data, skipping result parse") if ! grep /^1$/o, values %{$new_data};
	$self->log("Finished parsing and importing results");
  }
    
  return;
}
  

#Should be called from format parser e.g. BED, GFF, eQTL etc
#Why don't we pass feature_params and dbentry_params directly?

sub load_feature_and_xrefs{
  my $self = shift;
		  
  #warn "Loading ".($self->{_counts}{features}+1).' '.$self->feature_params->{-FEATURE_TYPE}->name."\n";
  #This now only fails once on the first run and then
  #Need to count based on feature_type?

  #new rather than new fast here as we want to validate the import
  my $feature = $self->{feature_class}->new(%{$self->feature_params});
  ($feature) = @{$self->feature_adaptor->store($feature)};
  $self->count('features');

  #Add count based on FeatureType, should be ftype name and analysis to reflect unique ftype key?



  ##This needs to be handled in caller as we are validating loci?
  #if($self->ucsc_coords){
  #	$start += 1;
  #	$end   += 1;
  #  }
  
  #This needs to be put in a separate sub and called by the caller
  #if(!  $self->cache_slice($chr)){
  #  warn "Skipping AnnotatedFeature import, cound non standard chromosome: $chr";
  #}
  #else{
  #grab seq if dump fasta and available		  
  #my $seq;
  #if(exists $self->feature_params->{'sequence'}){
  #	  $seq = $self->feature_params->{'sequence'};
  #	  delete $self->feature_params->{'sequence'};
  #	}
  #	else{
  #	  $self->log('No fasta sequence available for '.$self->feature_params->display_label);
  #	}
  #  }
	  
  #dump fasta here
  #if ($self->dump_fasta){
  #	$self->{'_fasta'} .= $self->generate_fasta_header($feature)."\n$seq\n";
  #  }
  
  #Store the xrefs
 
  foreach my $dbentry_hash(@{$self->{'_dbentry_params'}}){
	my $ftype = $dbentry_hash->{feature_type};
	delete $dbentry_hash->{feature_type};

	my $dbentry = Bio::EnsEMBL::DBEntry->new(%{$dbentry_hash});
	$self->dbentry_adaptor->store($dbentry, $feature->dbID, $ftype, 1);#1 is ignore release flag
	#count here? no count in caller
  }

  
  #Clean data cache
  $self->{'_feature_params'} = {};
  $self->{'_dbentry_params'} = [];

  return $feature;
}

#This should really be handled in Bio::EnsEMBL::Feature?
#Move to Helper?

sub set_strand{
  my ($self, $strand) = @_;

  my $ens_strand = 0;

  my %strand_vals = (
					 '1'  => 1,
					 '0'  => 0,
					 '-1' => -1,
					 '+'  => 1,
					 '-'  => -1,
					 '.'  => 0,
					);

  if($strand){
	
	if(exists $strand_vals{$strand}){
	  $ens_strand = $strand_vals{$strand};
	}
	else{
	  throw("Could not identify strand value for $strand");
	}
  }

  return $ens_strand;
}

sub total_features{
  my ($self, $total) = @_;

  $self->{'total_features'} = $total if defined $total;
  return $self->{'total_features'};
}

#Currently only required for Bed::parse_Features_by_Slice

#filehandle


sub file_handle{
  my ($self, $fh) = @_;

  $self->{'file_handle'} = $fh if defined $fh;
  return $self->{'file_handle'};
}

sub result_set{
  my ($self, $rset) = @_;

  #already tested/created by self
  
  $self->{'result_set'} = $rset if $rset;
  return $self->{'result_set'};
}


1;
