#
# EnsEMBL module for Bio::EnsEMBL::Funcgen::Parsers::InputSet
#

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
use Bio::EnsEMBL::Utils::Exception qw( throw warning deprecate );
use Bio::EnsEMBL::Funcgen::Utils::EFGUtils qw(species_chr_num open_file is_gzipped);
use Bio::EnsEMBL::Utils::Argument qw( rearrange );
use Bio::EnsEMBL::Funcgen::Utils::Helper;
use strict;


use vars qw(@ISA);
@ISA = qw(Bio::EnsEMBL::Funcgen::Utils::Helper);

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
  
  #Defaults
  #$self->{dbfile_data_root} = $self->get_dir('output');

  ($self->{'input_set_name'}, 
   $self->{'input_feature_class'}, 
   $self->{'slices'}, 
   $self->{total_features}, 
   $self->{force},#Is this generic enough to go in Importer? Similar to recover?
   $self->{dbfile_data_root}, #only appropriate for result input_feature_class
  ) = rearrange(['input_set_name', 'input_feature_class', 'slices', 'total_features', 'force', 'dbfile_data_root'], @_);


  

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
  
  $self->{'counts'}   = {};
  $self->{'slices'}   = [];
  $self->{'seq_region_names'} = [];#Used for slice based import


  return $self;
}

#These can be called for each record
#So we want to trim the access time as much as possible
#Should be slightly faster than using the reg
#Over lots of records will make some difference
#Only accessor as set in config for speed

sub input_feature_class{#annotated or result
  return $_[0]->{'input_feature_class'};
}

sub annotated_feature_adaptor{
  return $_[0]->{'annotated_feature_adaptor'};
}

sub result_feature_adaptor{
  return $_[0]->{'result_feature_adaptor'};
}


sub dbentry_adaptor{
  return $_[0]->{'dbentry_adaptor'};
}

sub input_set_adaptor{
  return $_[0]->{'input_set_adaptor'};
}

#This is either feature or result set
#Should be set in define sets
sub set{
  return $_[0]->{'set'};
}


sub slice_adaptor{
 my $self = shift;
  return $self->{'slice_adaptor'};
}


=head2 input_set_name
  
  Example    : my $input_set_name = $imp->input_set_name();
  Description: Getter for InputSet name
  Arg [1]    : optional - InputSet name
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : at risk

=cut

sub input_set_name{
  #Set in new
  return $_[0]->{'input_set_name'};
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

  if($self->input_feature_class ne 'result' && $self->input_feature_class ne 'annotated'){
	throw('You must define a valid set_type (result or annotated) to import using '.ref($self));
	#This will current print Importer but will reveal the correct parser when inheritance is fixed
  }


  #We need to undef norm method as it has been set to the env var
  $self->{'config'}{'norm_method'} = undef;

  #dirs are not set in config to enable generic get_dir method access
  $self->{dbfile_data_root} ||= $self->get_dir('output');#Required for Collector

  #some convenience methods
  $self->{'annotated_feature_adaptor'} = $self->db->get_AnnotatedFeatureAdaptor if $self->input_feature_class eq 'annotated';
  $self->{'result_feature_adaptor'} = $self->db->get_ResultFeatureAdaptor if $self->input_feature_class eq 'result';#Is this the right adaptor?
  $self->{'dbentry_adaptor'}           = $self->db->get_DBEntryAdaptor;
  $self->{'input_set_adaptor'}  = $self->db->get_InputSetAdaptor();
  $self->{'slice_adaptor'} = $self->db->dnadb->get_SliceAdaptor;

  $self->slices($self->{'slices'}) if defined $self->{'slices'};

  return;
}



sub define_sets{
  my ($self) = @_;

  my $eset = $self->db->get_InputSetAdaptor->fetch_by_name($self->input_set_name);
  
  if(! defined $eset){
	$eset = Bio::EnsEMBL::Funcgen::InputSet->new
	  (
	   -name         => $self->input_set_name(),
	   -experiment   => $self->experiment(),
	   -feature_type => $self->feature_type(),
	   -cell_type    => $self->cell_type(),
	   -vendor       => $self->vendor(),
	   -format       => $self->format(),
	   -analysis     => $self->feature_analysis,
	   -feature_class => $self->input_feature_class,
	  );
	($eset)  = @{$self->db->get_InputSetAdaptor->store($eset)};

  }

  #Use define_and_validate with fetch/append as we may have a pre-existing set
  #This now needs to handle ResultSets based on InputSets


  my $dset = $self->define_and_validate_sets
	(
	 -dbadaptor    => $self->db,
	 -name         => $self->input_set_name,
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
 
  return $self->{'_data_set'};

}

sub data_set{
  my $self = shift;
  return $self->{'_data_set'};
}

sub validate_files{
  my ($self, $prepare) = @_;

  #Get file
  if (! @{$self->result_files()}) {
	my $list = "ls ".$self->get_dir('input').'/'.$self->input_set_name().'*.';#.lc($self->vendor);#could use vendor here? ACtually need suffix attr
	my @rfiles = `$list`;
	$self->result_files(\@rfiles);
  }

  #We don't yet support multiple files
  if(scalar(@{$self->result_files()}) >1){
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
  my $recover_unimported = 0;
  my $dset = $self->data_set;
   
  if((scalar(@{$self->slices}) > 1) &&
	 ! $prepare){
	throw('Validate files does not yet support multiple Slice rollback');
  }

  #This all assumes that there is only ever 1 InputSet

  if ($self->input_feature_class eq 'annotated'){
	$eset =  $dset->get_supporting_sets->[0]; 
  }
  elsif($self->input_feature_class eq 'result'){
	$eset = $dset->get_supporting_sets->[0]->get_InputSets->[0];
  }

  #IMPORTED status here may prevent
  #futher slice based imports
  #so we have wait to set this until we know all the slices 
  #are loaded, unless we store slice based IMPORTED states
  #We currently get around this be never settign IMPORTED for slice based jobs
  #and always rolling back by slice before import
  


  foreach my $filepath( @{$self->result_files} ) {
	chomp $filepath;

 	my $filename;
	($filename = $filepath) =~ s/.*\///;
	my $sub_set;
	$self->log('Validating '.$self->vendor." file:\t$filename");
	throw("Cannot find ".$self->vendor." file:\t$filepath") if(! -e $filepath);#Can deal with links
	
	#reset filename to that originally used to create the Inputsubsets
	$filename =~ s/^prepared\.// if $self->prepared;
	
	if( $sub_set = $eset->get_subset_by_name($filename) ){
	  #IMPORTED status here is just for the file
	  #Any changes to analysis or cs would result in different Inputsubset/file
	  #so we don't need to account for that.
	  #Also will only ever be imported into one Feature|ResultSet

	  if($recover_unimported){
		$new_data{$filepath} = 1;
		next;
	  }

	  if( $sub_set->has_status('IMPORTED') ){
		$new_data{$filepath} = 0;
		$self->log("Found previously IMPORTED InputSubset:\t$filename");
	  } 
	  else{
		$self->log("Found partially IMPORTED InputSubset:\t$filename");
		$recover_unimported = 1;
		$new_data{$filepath} = 1;
		
		if(! $prepare){

		  #InputSet may be peaks or reads so how are we going to rollback?
		  #Need to pass parameter to Importer for feature/set type
		  #Given an InputSet in isolation there is no way of determining where it's
		  #features are stored. Do we need to add set_type to input_set?

		  if ( $self->recovery && $recover_unimported ) {
			$self->log("Rolling back results for InputSubset:\t".$filename);
			#Change these to logger->warn
			$self->log("WARNING::\tCannot yet rollback for just an InputSubset, rolling back entire set? Unless slices defined");
			$self->log("WARNING::\tThis may be deleting previously imported data which you are not re-importing..list?!!!\n");
			
			if($self->input_feature_class eq 'annotated'){
			  $self->rollback_FeatureSet($self->data_set->product_FeatureSet, undef, $self->slices->[0]);
			  $self->rollback_InputSet($eset);
			  last;
			}			
			elsif($self->input_feature_class eq 'result'){
			  #Can we do this by slice for parallelisation?
			  #This will only ever be a single ResultSet due to Helper::define_and_validate_sets
			  #flags are rollback_results and force(as this won't be a direct input to the product feature set)
			  $self->rollback_ResultSet($self->data_set->get_supporting_sets->[0], 1, $self->slices->[0], 1);
			}
		  }
		  elsif( $recover_unimported ){
			throw("Found partially imported InputSubSet:\t$filepath\n".
				  "You must specify -recover  to perform a full roll back for this InputSet:\t".$eset->name);
		  }
		}
	  }
	}
	else{

	  #This should never happed for prepared data
	  throw("Should not have found new 'prepared' file:\t$filename") if $self->prepared;
	  
	  $self->log("Found new InputSubset:\t${filename}");
	  $new_data{$filepath} = 1;
	  $sub_set = $eset->add_new_subset($filename);
	  $self->input_set_adaptor->store_InputSubsets([$sub_set]);
	}
  }

  #Set all the new if we have rolled back due to a recovery.
  if ($recover_unimported){

	foreach my $esset(@{$eset->get_InputSubsets}){
	  $new_data{$esset->name} = 1; 
	  $eset->adaptor->revoke_states($esset);
	}
  }

  return \%new_data;
}


#Separate setter and getter for speed;

sub set_feature_separator{
  my ($self, $separator) = @_;

  #How do we test if something undefined was passed?
  #Rather than nothing passed at all?
  #Can't do this as this is the accessor
  #Need to split method
 
  throw('Must provide a valid feature separator') if ( (! defined $separator) || ($separator eq '') ); 

  $self->{'_feature_separator'} = $separator;

}

sub feature_separator{
  my $self = shift;
  return $self->{'_feature_separator'};
}

#getter only
sub feature_params{
  my $self = shift;
  return $self->{'_feature_params'};
}

#getter only
sub dbentry_params{
  my $self = shift;
  return $self->{'_dbentry_params'};
}


sub counts{
  my ($self, $count_type) = @_;

  if($count_type){
	$self->{'_counts'}{$count_type} ||=0;
	return 	$self->{'_counts'}{$count_type};
  }
 
  return $self->{'_counts'}
}


#Move this to Importer?

sub slices{
  my ($self, $slices) = @_;

  if(defined $slices){

	if (ref($slices) ne 'ARRAY'){
	  throw("-slices parameter must be an ARRAYREF of Bio::EnsEMBL::Slices (i.e. not $slices)");
	}

	foreach my $slice(@$slices){
	  
	  if(! ($slice && ref($slice) && $slice->isa('Bio::EnsEMBL::Slice'))){
		throw("-slices parameter must be Bio::EnsEMBL::Slices (i.e. not $slice)");
	  }
	  
	  #Removed cache_slice from here as this was
	  #preventing us from identifying the seq_region in an input file

	  my $full_slice = $self->slice_adaptor->fetch_by_name($slice->name);

	  if(($slice->start != 1) ||
		 ($slice->end != $full_slice->end)){
		throw("InputSet Parser does not yet accomodate partial Slice based import i.e. slice start > 1 or slice end < slice length:\t".$slice->name);
		
	  }

	  push @{$self->{seq_region_names}}, $slice->seq_region_name;
	}
	$self->{'slices'} = $slices;
  }

  return $_[0]->{slices};
}


sub count{
  my ($self, $count_type) = @_;

  $self->{'_counts'}{$count_type} ||=0;
  $self->{'_counts'}{$count_type}++;
  return;
}


sub input_file_operator{
  my ($self, $op) = @_;
  #Should be set in format parser
  $self->{'input_file_operator'} = $op if defined $op;

  return $self->{'input_file_operator'};
}

sub input_gzipped{
  return  $_[0]->{'input_gzipped'};
}

# prepare boolean, simply stores the sets and preprocesses the file
# so we don't get each batch job trying to sort etc


#Still need to implement prepare in other Parsers!!

sub read_and_import_data{
  my ($self, $prepare) = @_;

  my $action = ($prepare) ? 'preparing' : 'importing';

  $self->log("Reading and $action ".$self->vendor()." data");
  my ($eset, $filename, $output_set, $fh, $f_out, %feature_params, @lines);

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

  my $dset = $self->define_sets;
  
  #We also need to account for bsub'd slice based import
  #seq alignments loaded into a ResultSet
  #Cannot have 0 window for ChIP Seq alignments
  #As this would mean storing all the individual reads
  #Hence we need to remap to a new assm before we import!
    
  if ($self->input_feature_class eq 'annotated'){
	$output_set = $dset->product_FeatureSet;
	$eset =  $dset->get_supporting_sets->[0]; 
  }
  elsif($self->input_feature_class eq 'result'){
	$output_set = $dset->get_supporting_sets->[0];
	$eset = $output_set->get_InputSets->[0];
	$self->result_set($output_set);#required for ResultFeature Collector and Bed Parser
  }
  
  #If we can do these the other way araound we can get define_sets to rollback the FeatureSet
  #Cyclical dependency for the sets :|
  my $new_data = $self->validate_files($prepare);
  my $seen_new_data = 0;
 
  
  ### READ AND IMPORT FILES ###
  foreach my $filepath(@{$self->result_files()}) {
	chomp $filepath;
	($filename = $filepath) =~ s/.*\///;

	$self->input_file($filepath);
	#This is only used by Collector::ResultFeature::reinitialise_input method


	#We're checking for recover here, as we have to reload all if just one has been screwed up.
	
	if( $new_data->{$filepath} ){
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

	  if((($self->input_feature_class eq 'result') && ! $prepare)){
		#(($self->input_feature_class eq 'result') && (! $self->batch_job))){   #Local run on just 1 chr
		#

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
		  $self->result_feature_adaptor->store_window_bins_by_Slice_Parser($slice, $self, 
																		   (
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

		#This slurp may need to go if data gets to large
		  
	  
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
		  
		  $self->log('Finished importing '.$self->counts('features').' '.
					 $output_set->name." features from:\t$filepath");
		  
		  
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
	$self->log("No new data, skipping result parse") if ! grep /1/,values %{$new_data};
	$self->log("Finished parsing and importing results");
  }
    
  return;
}
  

#Should be called from format parser e.g. BED, GFF, eQTL etc

sub load_feature_and_xrefs{
  my $self = shift;

  my $seq;

  #grab seq if dump fasta and available
  #if($self->dump_fasta){
	
  #if(exists $self->feature_params->{'sequence'}){
#	  $seq = $self->feature_params->{'sequence'};
#	  delete $self->feature_params->{'sequence'};
#	}
#	else{
#	  $self->log('No fasta sequence available for '.$self->feature_params->display_label);
#	}
#  }
		  
		  
  my $feature = Bio::EnsEMBL::Funcgen::AnnotatedFeature->new(%{$self->feature_params});
  ($feature) = @{$self->annotated_feature_adaptor->store($feature)};
  $self->count('features');


  ##This needs to be handled in caller as we are validating loci
  #if($self->ucsc_coords){
  #	$start += 1;
  #	$end   += 1;
  #  }
  
  #This needs to be put in a separate sub and called by the caller
  #if(!  $self->cache_slice($chr)){
  #  warn "Skipping AnnotatedFeature import, cound non standard chromosome: $chr";
  #}
  #else{
		  
		  
  #dump fasta here
  #if ($self->dump_fasta){
  #	$self->{'_fasta'} .= $self->generate_fasta_header($feature)."\n$seq\n";
  #  }
  
  #Now we need to store the xrefs
  
  #my $edb_ref = $self->db->dbc->db_handle->selectrow_arrayref('select db_name from external_db where db_name="ensembl_variation_Variation"');

  #if(! defined $edb_ref){
#	throw('Could not find external_db ensembl_variation_Variation');
#  }

 # my ($edbname) = @{$edb_ref};


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


  if(! defined $strand){
	return 0;
  }
  elsif($strand eq '+'){
	return 1;
  }
  elsif($strand eq '-'){
	return -1;
  } 

  return 0;

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
