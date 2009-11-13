#
# EnsEMBL module for Bio::EnsEMBL::Funcgen::Parsers::Bed
#

=head1 NAME

Bio::EnsEMBL::Funcgen::Parsers::Bed

=head1 SYNOPSIS

  my $parser_type = "Bio::EnsEMBL::Funcgen::Parsers::Bed";
  push @INC, $parser_type;
  my $imp = $class->SUPER::new(@_);


=head1 DESCRIPTION

This is a definitions class which should not be instatiated directly, it 
normally set by the Importer as the parent class.  Bed contains meta 
data and methods specific to data in bed format, to aid 
parsing and importing of experimental data.

=head1 CONTACT

Post questions to the EnsEMBL development list ensembl-dev@ebi.ac.uk

=head1 METHODS

=cut

package Bio::EnsEMBL::Funcgen::Parsers::Bed;

use Bio::EnsEMBL::Funcgen::Parsers::InputSet;
use Bio::EnsEMBL::Funcgen::FeatureSet;
use Bio::EnsEMBL::Funcgen::AnnotatedFeature;
use Bio::EnsEMBL::Utils::Exception qw( throw warning deprecate );
use Bio::EnsEMBL::Funcgen::Utils::EFGUtils qw(open_file);
use Bio::EnsEMBL::Utils::Argument qw( rearrange );
use Bio::EnsEMBL::Funcgen::Utils::Helper;
use File::Basename;
use strict;

use Data::Dumper;

use vars qw(@ISA);
@ISA = qw(Bio::EnsEMBL::Funcgen::Parsers::InputSet);

#To do
# Extend this so we can import features to annotated_feature
# and profiles/reads to result_feature
# This will replace individual tables for current bed das sources
# Altho' idea of separate tables for each set is not necessarily a bad one
# We would have to have extra table or fields details registration
# Easier admin i.e. drop tables
# What about partitioning? Irrelevant if we are moving to matrix
# Not easy to patch dynamically named tables? Can't assign values to user varaible from query!
# Would need to be a stored procedure
# Could have separate matrix files for result_features as we probably wouldn't want 
# to do any cross set querying at this level.
# Could we also use this in the RunnableDBs?
# Simply separate those methods required by both normal bed annotated_feature import
# and RunnableDB based annotated_feature import
# Importing into result_feature requires a result_set which assumes a chip experiment
# We need to alter result_set such that it can be optionally associated with ExperimentalSets
# Or should we just go for a different table? read_feature?
# NOTE: PARITION by key may run into problems if different sets have different windows
# Would also need to modify ResultFeatureAdaptor?

=head2 new

  Arg[0]     : hash containing optional attributes:
                          -bed_reads => 0|1, #Set input as read alignments rather than peak calls (default is 0)
  Example    : my $self = $class->SUPER::new(@_);
  Description: Constructor method for Bed class
  Returntype : Bio::EnsEMBL::Funcgen::Parsers::Bed
  Exceptions : throws if Experiment name not defined or if caller is not Importer
  Caller     : Bio::EnsEMBL::Funcgen::Importer
  Status     : at risk

=cut

sub new{
  my $caller = shift;
  
  my $class = ref($caller) || $caller;

  my $self = $class->SUPER::new(@_);
  

  #Bed specific params
  my ($bed_reads) = rearrange(['BED_READS'], @_);

  #This should be annotated or result which can then be used in define_and_validate_sets in Helper

  $self->{'set_feature_type'} = 'result' if $bed_reads;

  throw("This is a skeleton class for Bio::EnsEMBL::Importer, should not be used directly") 
	if(! $self->isa("Bio::EnsEMBL::Funcgen::Importer"));
  
  $self->{'config'} =  
	{(
      #order of these method arrays is important!
	  #remove method arrays and execute serially?
      array_data   => [],#['experiment'],
      probe_data   => [],#["probe"],
      results_data => ["and_import"],
      norm_method => undef,

	  #Need to make these definable?
	  #have protocolfile arg and just parse tab2mage protocol section format
	  #SEE NIMBLEGEN FOR EXAMPLE
	  protocols => {()},
     )};
	

  return $self;
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

  $self->SUPER::set_config;

  #Probably need to add MAQ/BOWTIE/BWA as analyses

  #dir are not set in config to enable generic get_dir method access
  #Need to define output dir if we are processing reads
  


  return;
}


sub pre_process_file{
  my ($self, $filepath) = @_;

  #We really just need to set file_operator here dependant on compression

  #If we are loading the reads then we want to gunzip anyway?
  #Not if we are using ensembl style table (i.e. does not match bed format)
  #Therefore need to write another output file?
  #Or do we just incorporate this as 0 size window in result_feature?


  #Will these %'s be interpolated here and fail in openfile?

  if($self->input_gzipped){
	$self->input_file_operator("gzip -dc %s | sort -k 1,3 |");
  }
  else{
	$self->input_file_operator("sort -k 1,3 %s |");
  }

  #We also need to set optional output filehandle here
  #So we can mysqlimport instead of using API to import
  #This should be put in output_dir
  #So we can point to an input file anywhere, and have our outputs in our own directory
  #We should also create links to input files in input_dir!
  #Let's just have one big file so there is no risk of having competing processes between imports


  
  if(! defined  $self->output_file && $self->set_feature_type eq 'result'){
	#output_file is currently only required for read mysqlimport loading
	#This could also be used by SAM/BAM so put in InputSet

	#Or do we need a file for each discrete set?
	#Each import consititutes one discrete data set
	#So this is okay for replicates in a result set
	#But not for different cell/feature types
	#Use one file
	#This is imposed case in INputSet::validate_files
	my ($name) = fileparse($filepath);

	$name =~ s/\.gz// if $self->input_gzipped;
	#We can't pipe this to gzip as we need to mysqlimport it, which does not support gzip?
	#gzip -dc file.sql | mysql would work but would be much slower due ti individual inserts
   	$self->output_file($self->get_dir('output')."/result_feature.${name}.txt");
  }


  #return $new_path;
  return $filepath;
}


sub parse_line{
  my ($self, $line) = @_;

  return if $line !~ /^chr/io;

  my ($chr, $start, $end, $name, $score, $strand, @other_fields) = split/\t/o, $line;				  
  #my ($chr, $start, $end, $score) = split/\t/o, $line;#Mikkelson hack	
  #Validate variables types here beofre we get a nasty error from bind_param?

  #Any more valid BED fields here?
  # thickStart - The starting position at which the feature is drawn thickly (for example, the start codon in gene displays).
  # thickEnd - The ending position at which the feature is drawn thickly (for example, the stop codon in gene displays).
  # itemRgb - An RGB value of the form R,G,B (e.g. 255,0,0). If the track line itemRgb attribute is set to "On", this RBG value will determine the display color of the data contained in this BED line. NOTE: It is recommended that a simple color scheme (eight colors or less) be used with this attribute to avoid overwhelming the color resources of the Genome Browser and your Internet browser.
  # blockCount - The number of blocks (exons) in the BED line.
  # blockSizes - A comma-separated list of the block sizes. The number of items in this list should correspond to blockCount.
  # blockStarts - A comma-separated list of block starts. All of the blockStart positions should be calculated relative to chromStart. The number of items in this list should correspond to blockCount. 
  
  $strand = $self->set_strand($strand);			   

  if($self->ucsc_coords){
	$start +=1;
  }
		   
  if(!  $self->cache_slice($chr)){
	warn "Skipping AnnotatedFeature import, cound non standard chromosome: $chr";
  }
  else{
	
	#This is generic count handled in InputSet 
	$self->count('features');
		 
	#This is currently only loading peaks?
	#We need to plug the collection code in here for reads?
	#Also need to optionally have as result set or feature sets dependant on bed_type = reads | peaks
	
	#Need to convert all of this to use InputSet parser methods
	#This is single line format, so we don't need to use $self->feature_params or record_separator

	my $feature = Bio::EnsEMBL::Funcgen::AnnotatedFeature->new
	  (
	   -START         => $start,
	   -END           => $end,
	   -STRAND        => $strand,
	   -SLICE         => $self->cache_slice($chr),
	   -ANALYSIS      => $self->data_set->product_FeatureSet->analysis,
	   -SCORE         => $score,
	   -DISPLAY_LABEL => $name,
	   -FEATURE_SET   => $self->data_set->product_FeatureSet,
	  );

					 
	$self->annotated_feature_adaptor->store($feature);
			 
	#Unlikely we will ever need this
	#dump fasta here
	#if ($self->dump_fasta()){
	#  #We need to prefix this with the dbID?
	#  #And put the fasta in a host:port:dbname specfific dir?
	#  $fasta .= '>'.$name."\n".$self->cache_slice($chr)->sub_Slice($start, $end, 1)->seq()."\n";
	#}
  }
  
  return;
}

1;
