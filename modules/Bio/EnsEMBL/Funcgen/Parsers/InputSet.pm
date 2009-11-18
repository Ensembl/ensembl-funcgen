#
# EnsEMBL module for Bio::EnsEMBL::Funcgen::Parsers::InputSet
#

=head1 NAME

Bio::EnsEMBL::Funcgen::Parsers::InputSet

=head1 SYNOPSIS

  use vars qw(@ISA);
  @ISA = qw(Bio::EnsEMBL::Funcgen::Parsers::InputSet);

  
=head1 DESCRIPTION

This is a base class to support simple file format parsers. For simple imports the vendor is
set to the parser type i.e. the file format.  The generic read_and_import_simple_data assumes
a one line per feature format, other format need there own read_and_import_format_data method, 
which will need defining in the result_data config element.


=head1 AUTHOR

This module was created by Nathan Johnson.

=head1 CONTACT

Post questions to the EnsEMBL development list ensembl-dev@ebi.ac.uk

=head1 METHODS

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
use Bio::EnsEMBL::Funcgen::Utils::EFGUtils qw(species_chr_num open_file);
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
  
  ($self->{'input_set_name'}) = rearrange(['input_set_name'], @_);

  #No rollback flag yet to void losing old data which we are not reimporting


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
  $self->{'set_feature_type'} = 'annotated'; 

  return $self;
}

#These can be called for each record
#So we want to trim the access time as much as possible
#Should be slightly faster than using the reg
#Over lots of records will make some difference
#Only accessor as set in config for speed

sub set_feature_type{
  return $_[0]->{'set_feature_type'};
}

sub annotated_feature_adaptor{
  return $_[0]->{'annotated_feature_adaptor'};
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


  #We could set input_set_name to experiment name
  #But we would have to make warning in define_and_validate_sets mention -input_set_name

  throw('Must provide an -input_set name for a '.uc($self->vendor).' import') if ! defined $self->input_set_name();

  #Mandatory checks
  if(! defined $self->feature_analysis){
	throw('Must define a -feature_analysis parameter for '.uc($self->vendor).' imports');
  }


  #We need to undef norm method as it has been set to the env var
  $self->{'config'}{'norm_method'} = undef;

  #dir are not set in config to enable generic get_dir method access


  #some convenience methods
  $self->{'annotated_feature_adaptor'} = $self->db->get_AnnotatedFeatureAdaptor if $self->set_feature_type eq 'annotated';
  $self->{'result_feature_adaptor'} = $self->db->get_ResultFeatureAdaptor if $self->set_feature_type eq 'result';#Is this the right adaptor?
  $self->{'dbentry_adaptor'}           = $self->db->get_DBEntryAdaptor;
  $self->{'input_set_adaptor'}  = $self->db->get_InputSetAdaptor();
  $self->{'slice_adaptor'} = $self->db->dnadb->get_SliceAdaptor;
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
	 -type         => $self->set_feature_type, 
	 -description  => $self->feature_set_description,
	 #-append          => 1,#Omit append to ensure we only have this eset
	 -recovery     => $self->recovery,
	 -supporting_sets => [$eset],
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
  my $self = shift;


  #Get file
  if (! @{$self->result_files()}) {
	my $list = "ls ".$self->get_dir('input').'/'.$self->input_set_name().'*.';#.lc($self->vendor);#could use vendor here? ACtually need suffix attr
	my @rfiles = `$list`;
	$self->result_files(\@rfiles);
  }
  
  if (scalar(@{$self->result_files()}) >1) {
	warn('Found more than one '.$self->vendor." file:\n".
		 join("\n", @{$self->result_files()})."\nThe InputSet parser does not yet handle replicates.".
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
  my (%new_data);
  my $recover_unimported = 0;
  my ($eset) = @{$self->data_set->get_supporting_sets};
  
  foreach my $filepath( @{$self->result_files} ) {
	chomp $filepath;

 	my $filename;
	($filename = $filepath) =~ s/.*\///;
	my $sub_set;
	$self->log('Validating '.$self->vendor." file:\t$filename");
	throw("Cannot find ".$self->vendor." file:\t$filename") if(! -e $filepath);#Can deal with links
	

	if( $sub_set = $eset->get_subset_by_name($filename) ){
	  
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
		
		if ( $self->recovery && $recover_unimported ) {
		  $self->log("Rolling back results for InputSubset:\t".$filename);
		  #Change these to logger->warn
		  $self->log("WARNING::\tCannot yet rollback for just an InputSubset, rolling back entire set");
		  $self->log("WARNING::\tThis may be deleting previously imported data which you are not re-importing..list?!!!\n");
		  $self->rollback_FeatureSet($self->data_set->product_FeatureSet);
		  $self->rollback_InputSet($eset);
		  last;
		}
		elsif( $recover_unimported ){
		  throw("Found partially imported InputSubSet:\t$filepath\n".
				"You must specify -recover  to perform a full roll back for this InputSet:\t".$eset->name);
		}
	  }
	}
	else{
	  $self->log("Found new InputSubset:\t${filename}");
	  $new_data{$filepath} = 1;
	  $sub_set = $eset->add_new_subset($filename);
	  $self->input_set_adaptor->store_InputSubsets([$sub_set]);
	}
  }

  #Set all the new if we have rolled back due to a recovery.
  if ($recover_unimported){

	foreach my $esset(@{$eset->get_subsets}){
	  #map $new_data{$_} = 1, keys %new_data if $recover_unimported;
	  $new_data{$esset->name} = 1; 
	  $eset->adaptor->revoke_states($esset);
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

sub read_and_import_data{
  my $self = shift;
    
  $self->log("Reading and importing ".$self->vendor()." data");
  my ($filename, $fh, $f_out, %feature_params, @lines);
  my $dset   = $self->define_sets;
  my $fset   = $dset->product_FeatureSet;

  #If we can do these the other way araound we can get define_sets to rollback the FeatureSet
  #Cyclical dependency for the sets :|
  my ($eset) = @{$dset->get_supporting_sets};   
  my $new_data = $self->validate_files;
  my $seen_new_data = 0;

  ### READ AND IMPORT FILES ###
  foreach my $filepath(@{$self->result_files()}) {
	chomp $filepath;
	($filename = $filepath) =~ s/.*\///;

	#We're checking for recover here, as we have to reload all if just one has been screwed up.
	
	if( $new_data->{$filepath} ){
	  $seen_new_data = 1;

	  #Do standard gzip test first
	  my $compressed_data =  `file -L $filepath` or die "Can't execute 'file -L $filepath'";
	  $self->{'input_gzipped'} = 1 if $compressed_data =~ /gzip/;
	  if($compressed_data =~ /compressed/ && ! $self->input_gzipped){
		throw("Bio::Ensembl::Funcgen::Parsers::ExternalSet only handles gzip compressed files, please uncompress $filepath manually before rerunning");
	  }
	  

	  $filepath = $self->pre_process_file($filepath) if $self->can('pre_process_file');

	  $self->log_header('Reading '.$self->vendor." file:\t".$filepath);

	  #We need to be able to optional open pipe to gzip | sort here
	  #i.e. define open command
	  $fh = open_file($filepath, $self->input_file_operator);

	  #This my become way too large for some reads files
	  #Currently no problems
	  #This is not working as we are sorting the file!
	  #$self->parse_header($fh) if $self->can('parse_header');

	  my @lines = <$fh>;
	  close($fh);
	  
	  #Revoke FeatureSet IMPORTED state here incase we fail halfway through
	  $fset->adaptor->revoke_status('IMPORTED', $fset) if $fset->has_status('IMPORTED');
	  #What about IMPORTED_"CSVERSION"
	  #This may leave us with an incomplete import which still has
	  #an IMPORTED_CSVERSION state
	  #We need to depend on IMPORTED for completeness of set
	  #DAS currently only uses IMPORTED_CSVERSION
	  #This is okayish but we also need to write HCs for any sets 
	  #which do not have IMPORTED state!

	 	  
	  foreach my $line (@lines) {
		#Generic line processing
		#Move these to parse_line?
		$line =~ s/\r*\n//o;
		next if $line =~ /^\#/;	
		next if $line =~ /^$/;

		$self->count('Total lines');

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

		$self->parse_line($line);		

		#Here we depedant on whether we are parsing reads and using collection
		#We need to pass a coderef to the collection somehow, so this can read the data
		#either from a query or from a file
		#

	  }
	 
	  #Now we need to deal with anything left in the read cache
	  $self->process_params if $self->can('process_params');
	  
	  #To speed things up we may need to also do file based import here with WRITE lock?
	  #mysqlimport will write lock the table by default?
	 
	  
	  $self->log('Finished importing '.$self->counts('features').' '.
				 $fset->name." features from:\t$filepath");
	  
	  foreach my $key (%{$self->counts}){
		$self->log("Count $key:\t".$self->counts($key));
	  }

	  my $sub_set = $eset->get_subset_by_name($filename);
	  $sub_set->adaptor->store_status('IMPORTED', $sub_set);
	  
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

  $self->set_imported_states_by_Set($fset) if $seen_new_data;

  $self->log("No new data, skipping result parse") if ! grep /1/,values %{$new_data};
  $self->log("Finished parsing and importing results");  
    
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

1;
