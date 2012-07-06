
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

Bio::EnsEMBL::Funcgen::RegulatoryFeature

=head1 SYNOPSIS

 use Bio::EnsEMBL::Registry;
 use Bio::EnsEMBL::Funcgen::RegulatoryFeature;
 my $reg = Bio::EnsEMBL::Registry->load_adaptors_from_db
     (
      -host    => 'ensembldb.ensembl.org',
      -user    => 'anonymous'
     );

 my $regfeat_adaptor = $reg->get_adaptor($species, 'funcgen', 'RegulatoryFeature');


 ### Creating/storing a RegulatoryFeature Set ###
 my $feature = Bio::EnsEMBL::Funcgen::RegulatoryFeature->new
     (
      -SLICE         => $chr_1_slice,
      -START         => 1_000_000,
	  -END           => 1_000_024,
      -STRAND        => 0,
      -DISPLAY_LABEL => $text,
      -FEATURE_SET   => $fset,
      -FEATURE_TYPE  => $reg_ftype,
      -ATTRIBUTE_CACHE => \%attr_cache,
     );

 my ($stored_feat) = @{$regfeat_adaptor->store([$feature])};


 ### Fetching some RegualtoryFeatures
 my @regfeats = @{$regfeat_adaptor->fetch_all_by_Slice_FeatureSets($slice, \@fsets)};


 ### Print the bound and core loci
 print join(' - ', ($reg_feat->bound_start,
                    $reg_feat->start,
                    $reg_feat->end,
                    $reg_feat->bound_end)."\n";


 ### Getting some supporting evidence for a RegualtoryFeatures
 my @reg_attrs = @{$reg_feat->regulatory_attributes('annotated')};


=head1 DESCRIPTION

A RegulatoryFeature object represents the output of the Ensembl RegulatoryBuild:
    http://www.ensembl.org/info/docs/funcgen/regulatory_build.html

It is comprises many possible histone, transcription factor, polymerase and open
chromatin features, which have been combined to provide a summary view and
classification of the regulatory status at a given loci.


=head1 SEE ALSO

Bio::EnsEMBL:Funcgen::DBSQL::RegulatoryFeatureAdaptor
Bio::EnsEMBL::Funcgen::SetFeature

=cut


package Bio::EnsEMBL::Funcgen::RegulatoryFeature;

use Bio::EnsEMBL::Utils::Argument qw( rearrange );
use Bio::EnsEMBL::Utils::Exception qw( throw );
use strict;
use warnings;

use base qw(Bio::EnsEMBL::Funcgen::SetFeature); #@ISA


=head2 new

  Arg [-SLICE]             : Bio::EnsEMBL::Slice - The slice on which this feature is located.
  Arg [-START]             : int - The start coordinate of this feature relative to the start of the slice
                             it is sitting on. Coordinates start at 1 and are inclusive.
  Arg [-END]               : int -The end coordinate of this feature relative to the start of the slice
                    	     it is sitting on. Coordinates start at 1 and are inclusive.
  Arg [-FEATURE_SET]       : Bio::EnsEMBL::Funcgen::FeatureSet - Regulatory Feature set
  Arg [-FEATURE_TYPE]      : Bio::EnsEMBL::Funcgen::FeatureType - Regulatory Feature sub type
  Arg [-BINARY_STRING]     : (optional) string - Regulatory Build binary string
  Arg [-STABLE_ID]         : (optional) string - Stable ID for this RegualtoryFeature e.g. ENSR00000000001
  Arg [-DISPLAY_LABEL]     : (optional) string - Display label for this feature
  Arg [-ATTRIBUTE_CACHE]   : (optional) HASHREF of feature class dbID|Object lists
  Arg [-PROJECTED]         : (optional) boolean - Flag to specify whether this feature has been projected or not
  Arg [-dbID]              : (optional) int - Internal database ID.
  Arg [-ADAPTOR]           : (optional) Bio::EnsEMBL::DBSQL::BaseAdaptor - Database adaptor.

  Example    : my $feature = Bio::EnsEMBL::Funcgen::RegulatoryFeature->new(
										                                  -SLICE         => $chr_1_slice,
									                                      -START         => 1000000,
									                                      -END           => 1000024,
                                                                          -DISPLAY_LABEL => $text,
									                                      -FEATURE_SET   => $fset,
                                                                          -FEATURE_TYPE  => $reg_ftype,
                                                                          -ATTRIBUTE_CACHE => \%attr_cache,
                                                                         );


  Description: Constructor for RegulatoryFeature objects.
  Returntype : Bio::EnsEMBL::Funcgen::RegulatoryFeature
  Exceptions : None
  Caller     : General
  Status     : Stable

=cut

sub new {
  my $caller = shift;
  my $class = ref($caller) || $caller;
  my $self = $class->SUPER::new(@_);

  my ($stable_id, $attr_cache, $bin_string, $projected)
    = rearrange(['STABLE_ID', 'ATTRIBUTE_CACHE', 'BINARY_STRING', 'PROJECTED'], @_);

  #None of these are mandatory at creation
  #under different use cases
  $self->{binary_string} = $bin_string if defined $bin_string;
  $self->{stable_id}     = $stable_id  if defined $stable_id;
  $self->{projected}     = $projected  if defined $projected;
  $self->attribute_cache($attr_cache)  if $attr_cache;

  return $self;
}


=head2 display_label

  Example    : my $label = $feature->display_label();
  Description: Getter for the display label of this feature.
  Returntype : String
  Exceptions : None
  Caller     : General
  Status     : Stable

=cut

sub display_label {
  my $self = shift;

  if(! defined $self->{'display_label'}){
	$self->{'display_label'}  = $self->feature_type->name.' Regulatory Feature';

	if( defined $self->cell_type ){
	  $self->{'display_label'} .= ' - '.$self->cell_type->name;
	}
  }

  return  $self->{display_label};
}

=head2 display_id

  Example    : print $feature->display_id();
  Description: This method returns a string that is considered to be
               the 'display' identifier. In this case the stable Id is
               preferred
  Returntype : String
  Exceptions : none
  Caller     : web drawing code, Region Report tool
  Status     : Stable

=cut

sub display_id {  return $_[0]->{stable_id}; }


=head2 binary_string

  Arg [1]    : optional string - binary string from regulatory build
  Example    : my $bin_string = $feature->binary_string();
  Description: Getter and setter for the binary_string for this feature.
  Returntype : String
  Exceptions : None
  Caller     : General
  Status     : At Risk - May change to BLOB, remove setter functionality

=cut

sub binary_string{
  my ($self, $bin_string)  = @_;

  if (defined $bin_string){
	#added v67
	warn "RegualtoryFeature::binary_string setter functionality is being removed\n";
	$self->{binary_string} = $bin_string;
  }

  return $self->{binary_string};
}


=head2 stable_id

  Arg [1]    : (optional) string - stable_id e.g ENSR00000000001
  Example    : my $stable_id = $feature->stable_id();
  Description: Getter and setter for the stable_id attribute for this feature. 
  Returntype : string
  Exceptions : None
  Caller     : General
  Status     : At risk - setter functionality to be removed

=cut

sub stable_id {
  my $self = shift;

  if (@_){
    #added v67
    warn "RegualtoryFeature::stable_id setter functionality is being removed\n";
    $self->{stable_id} = shift;
  }

  return $self->{stable_id};
}


=head2 regulatory_attributes

  Arg [1]    : String (optional) - Class of feature e.g. annotated or motif
  Example    : print "Regulatory Attributes:\n\t".join("\n\t", (map $_->feature_type->name, @{$feature->regulatory_attributes()}))."\n";
  Description: Getter for the regulatory_attributes for this feature.
  Returntype : ARRAYREF
  Exceptions : Throws if feature class not valid
  Caller     : General
  Status     : At Risk

=cut


sub regulatory_attributes{
  my ($self, $feature_class) = @_;

  #Incorporating the MFs like this does cause some redundancy in the DB
  #But will speed up display of the RegFeat image including the MFs
  #Redefine the cache to have class keys e.g. TFBS, OpenChromatin, Histone Mods
  #Can't do this as we need the table key to be able to fetch the features
  #Really need something to be able to draw the image first, then create the zmenu details later.

  my %adaptors = (
				  'annotated' => $self->adaptor->db->get_AnnotatedFeatureAdaptor,
				  'motif'     => $self->adaptor->db->get_MotifFeatureAdaptor,
				  #external
				 );

  my @fclasses;

  if(defined $feature_class){

	if(exists $adaptors{lc($feature_class)}){
	  @fclasses = (lc($feature_class));
	}
	else{
	  throw("The feature class you specified is not valid:\t$feature_class\n".
			"Please use one of:\t".join(', ', keys %adaptors));
	}
  }
  else{
	@fclasses = keys %adaptors;
  }

  foreach my $fclass(@fclasses){
	#Now structured as hash to facilitate faster has_attribute method
	#Very little difference to array based cache

	my @attr_dbIDs = keys %{$self->{'attribute_cache'}{$fclass}};

	
	if(scalar(@attr_dbIDs) > 0){
	  
	  if( ! ( ref($self->{'regulatory_attributes'}{$fclass}->[0])  &&
			  ref($self->{'regulatory_attributes'}{$fclass}->[0])->isa('Bio::EnsEMBL::Feature') )){
		
		$adaptors{$fclass}->force_reslice(1);#So we don't lose attrs which aren't on the slice
		$self->{'regulatory_attributes'}{$fclass} = $adaptors{$fclass}->fetch_all_by_dbID_list(\@attr_dbIDs, $self->slice);

		#Having problems here if we are trying to project between Y PAR and X
		#Current dest_slice mapping code simply changes the start end values assuming the slice is correct
		#currently no test for seq_region name match
		

		#foreach my $attr(@{	$self->{'regulatory_attributes'}{$fclass}}){
		#  warn "$attr ".$attr->dbID." ".$attr->feature_Slice->name."\n";
		#}


		$adaptors{$fclass}->force_reslice(0);

		#Problems here with attrs not being returning when they do not lie on dest slice
		#i.e. core projected to cell line, but dest slice only over laps a region of the core which
		#actually has no attrs.
		#either use the feature_Slice and reslice everthing to the dest slice
		#or skip test in attr obj_frm_sth?
		#

		#This method transfers to the query slice, do not use fetch_by_dbID
		#It also should use _final_clause
		#This is currently only specified in the MotifFeatureAdaptor
		#as these are required to be sorted to relate to the structure string

		#but we are stll storing in has where order is not preserved!!
		#so this will not match order of underlying strcture!

		#separate so we can have ordered array returned
		#do we need redundant caches?
		#defo need db id cache for 'has' methods
		
		#foreach my $attr(@{$fclass_attrs}){
		#  $self->{'regulatory_attributes'}{$fclass}{$attr->dbID} = $attr;
		#}
	  }
	}
	else{
	  $self->{'regulatory_attributes'}{$fclass} = [];
	}
  }

  return [ map { @{$self->{'regulatory_attributes'}{$_}} } @fclasses ];
}

=head2 has_attribute

  Arg [1]    : Attribute Feature dbID
  Arg [2]    : Attribute Feature class e.g. motif or annotated
  Example    : if($regf->has_attribute($af->dbID, 'annotated'){ #do something here }
  Description: Identifies whether this RegualtoryFeature has a given attribute
  Returntype : Boolean
  Exceptions : Throws if args are not defined
  Caller     : General
  Status     : Stable

=cut


sub has_attribute{
  my ($self, $dbID, $fclass) = @_;

  throw('Must provide a dbID and a Feature class argument') if ! $dbID && $fclass;

  return exists ${$self->attribute_cache}{$fclass}{$dbID};
}

=head2 get_focus_attributes

  Arg [1]    : None
  Example    : my @focus_attrs = @{$regf->get_focus_attributes};
  Description: Getter for the focus features of this RegualtoryFeature, used to defined the core region
  Returntype : ARRAYREF
  Exceptions : None
  Caller     : General
  Status     : Stable

=cut

sub get_focus_attributes{
  my $self = shift;

  if(! exists $self->{'focus_attributes'} ||
	 ! @{$self->{'focus_attributes'}}){
	$self->_sort_attributes;
  }


  return $self->{'focus_attributes'};
}


=head2 get_nonfocus_attributes

  Arg [1]    : None
  Example    : my @non_focus_attrs = @{$regf->get_nonfocus_attributes};
  Description: Getter for the non-focus features of this RegulatoryFeature, used to defined 
               the non core region i.e. the whiskers.
  Returntype : ARRAYREF
  Exceptions : None
  Caller     : General
  Status     : Stable

=cut

sub get_nonfocus_attributes{
  my $self = shift;

  #Test focus here as we may not have any nonfocus
  #But focus will show that we have sorted already
  if(! exists $self->{'focus_attributes'} ||
	 ! @{$self->{'focus_attributes'}}){
	$self->_sort_attributes;
  }

  return $self->{'nonfocus_attributes'};
}

#Add pod here

sub _sort_attributes{
  my $self = shift;

  $self->{'focus_attributes'} = [];
  $self->{'nonfocus_attributes'} = [];

  foreach my $attrf(@{$self->regulatory_attributes}){

	if($attrf->isa('Bio::EnsEMBL::Funcgen::MotifFeature') ||
	   $attrf->feature_set->is_focus_set){
	  push @{$self->{'focus_attributes'}}, $attrf;
	}
	else{
	  push @{$self->{'nonfocus_attributes'}}, $attrf;
	}
  }

  return;
}


=head2 attribute_cache

  Arg [1]    : optional - HASHREF of attribute table keys with values as either a list of attribute 
               feature dbIDs or objects. If passing object, any MotifFeature objects should be in position
               order with respect to the slice.
  Example    : $feature->attribute_cache(\%attribute_feature_info);
  Description: Setter for the regulatory_attribute cache for this feature. This is a short cut method used by the 
               regulatory build and the webcode to avoid unnecessary fetching and enable enable lazy loading 
  Returntype : HASHREF
  Exceptions : Throws if trying to overwrite existing cache
  Caller     : RegulatoryFeatureAdaptor.pm and build_regulatory_features.pl
  Status     : At Risk

=cut


sub attribute_cache{
  my ($self, $attr_hash) = @_;

#  if(! defined $attr_hash){
#	$self->regulatory_attributes; #Fetch the attrs?
#
#
#	#Do we need to do this now we have separated the caches?
#
#  }

  if(defined $attr_hash){

	foreach my $fclass(keys %{$attr_hash}){

	  if(exists $self->{'attribute_cache'}{$fclass}){
		throw("You are trying to overwrite a pre-existing regulatory attribute cache entry for feature class:\t$fclass");
	  }
	  else{
		$self->{'attribute_cache'}{$fclass} = $attr_hash->{$fclass};
	  }
	}
  }

  return $self->{'attribute_cache'} || {};
}


=head2 bound_start

  Example    : my $bound_start = $feature->bound_start();
  Description: Getter for the bound_start attribute for this feature.
               Gives the 5' most start value of the underlying attribute
               features.
  Returntype : string
  Exceptions : None
  Caller     : General
  Status     : Stable

=cut

sub bound_start {
  my $self = shift;
  $self->get_underlying_structure if ! defined $self->{'bound_start'};

  return $self->{'bound_start'};
}


=head2 bound_end

  Example    : my $bound_end = $feature->bound_start();
  Description: Getter for the bound_end attribute for this feature.
               Gives the 3' most end value of the underlying attribute
               features.
  Returntype : string
  Exceptions : None
  Caller     : General
  Status     : Stable

=cut

sub bound_end {
  my $self = shift;
  $self->get_underlying_structure if ! defined $self->{'bound_end'};

  return $self->{'bound_end'};
}


=head2 bound_seq_region_start

  Example    : my $bound_sr_start = $feature->bound_seq_region_start;
  Description: Getter for the seq_region bound_start attribute for this feature.
               Gives the 5' most start value of the underlying attribute
               features.
  Returntype : string
  Exceptions : None
  Caller     : General
  Status     : Stable

=cut

sub bound_seq_region_start {
  my $self = shift;

  if(! defined $self->{bound_seq_region_start}){
	
	if($self->slice->strand == 1){
	  $self->{bound_seq_region_start} = $self->slice->start + $self->bound_start - 1;
	}
	else{ #strand = -1
	  $self->{bound_seq_region_start} = $self->slice->end - $self->bound_end + 1;
	}
  }

  return $self->{bound_seq_region_start};
}


=head2 bound_seq_region_end

  Example    : my $bound_sr_end = $feature->bound_seq_region_end;
  Description: Getter for the seq_region bound_end attribute for this feature.
               Gives the 3' most end value of the underlying attribute
               features.
  Returntype : string
  Exceptions : None
  Caller     : General
  Status     : Stable

=cut


sub bound_seq_region_end {
  my $self = shift;

  if(! defined $self->{bound_seq_region_end}){
	
	if($self->slice->strand == 1){
	  $self->{bound_seq_region_end} = $self->slice->start + $self->bound_end - 1;
	}
	else{ #strand = -1
	  $self->{bound_seq_region_end} = $self->slice->end - $self->bound_start + 1;
	}
  }

  return $self->{bound_seq_region_end};
}






=head2 is_projected

  Arg [1]    : optional - boolean
  Example    : if($regf->is_projected){ #do something different here }
  Description: Getter/Setter for the projected attribute.
  Returntype : boolean
  Exceptions : None
  Caller     : General
  Status     : At risk - remove setter functionality

=cut

sub is_projected {
  my $self = shift;
  
  if(@_){
	#added v67
	warn "RegulatoryFeature::is_projected setter functionality is being removed\n";
	$self->{'projected'} = shift;
  }
  
  return $self->{'projected'};
}


=head2 get_underlying_structure

  Example    : $self->get_underlying_structure() if(! exists $self->{'bound_end'});
  Description: Getter for the bound_end attribute for this feature.
               Gives the 3' most end value of the underlying attribute
               features.
  Returntype : string
  Exceptions : None
  Caller     : General
  Status     : At Risk

=cut

#This should really be precomputed and stored in the DB to avoid the MF attr fetch
#Need to be aware of projecting here, as these will expire if we project after this method is called

sub get_underlying_structure{
  my $self = shift;

  if(! defined $self->{underlying_structure}){

	my @attrs = @{$self->regulatory_attributes()};

	if(! @attrs){
	  throw('No underlying regulatory_attribute features to get_underlying_structure for dbID '.$self->dbID);
	  #This should never happen even with a projection build
	}
	else{


	  #We only need to set the bounds when storing on full slice/seq_region values
	  #else they should be fetched from the DB

	  if(! defined $self->{'bound_start'}){

		my (@start_ends);

		foreach my $attr(@attrs){
		  push @start_ends, ($attr->start, $attr->end);
		}

		#Accounts for core region, where data may be absent on this cell type
		push @start_ends, ($self->start, $self->end);

		@start_ends = sort { $a <=> $b } @start_ends;

		$self->{'bound_end'} = pop @start_ends;
		$self->{'bound_start'} = shift @start_ends;

		#Need to account for projection build here
		#i.e. attr extremeties may not extend past core start/end

		if($self->is_projected){
		  $self->{'bound_end'}   = $self->end   if $self->end   > $self->{'bound_end'};
		  $self->{'bound_start'} = $self->start if $self->start < $self->{'bound_start'};
		}
	  }

	  #Now deal with MotifFeature loci
	  my @mf_loci;

	  foreach my $mf(@{$self->regulatory_attributes('motif')}){
		push @mf_loci, ($mf->start, $mf->end);
	  }

	  $self->{underlying_structure} = [$self->{'bound_start'}, $self->start, @mf_loci, $self->end, $self->{'bound_end'}];
	}
  }

  return $self->{underlying_structure};
}


=head2 is_unique_to_FeatureSets

  Arg[1]     : optional - ARRAYREF of regualtory Bio::EnsEMBL::Funcgen::FeatureSet objects
                          Default is FeatureSet of given RegulatoryFeature, else need to be 
                          defined explicitly.
  Arg[2]     : optional - HASHREF Params hash: 
                                    {
                                     include_projected => 0|1, # Boolean, include 'projected' features
                                    }
  Example    : if($reg_feat->is_unique_to_FeatureSets($fsets)}{  
                   #then do some analysis here
               }
  Description: Identifies whether this RegualtoryFeature is unique to a set of FeatureSets.
  Returntype : boolean
  Exceptions : Throw is arguments not stored or valid.
  Caller     : General
  Status     : At risk

=cut

#Probably want to add in an FeatureType constraint here
#e.g. so we can compare active vs inactive or poised promoters

#omit include_multi doesn't make sense here

sub is_unique_to_FeatureSets{
  my ($self, $fsets, $params_hash) = @_;

  $fsets ||= [$self->feature_set];
  my @fset_ids;
  
  
  #define to avoid deref fails below.
  $params_hash ||= {};
  if(ref($params_hash) ne 'HASH'){
	throw("The params hash argument must be a valid HASHREF:\t".ref($params_hash));
  }


  foreach my $fset(@$fsets){
	#assume we have an adaptor set
	$self->adaptor->db->is_stored_and_valid('Bio::EnsEMBL::Funcgen::FeatureSet', $fset);
	
	if($fset->feature_class ne 'regulatory'){
	  throw('Found non-regulatory FeatureSet');
	}

	push @fset_ids, $fset->dbID;
  }
  
  my $stable_id;
  ($stable_id = $self->stable_id) =~ s/^[A-Z0]+//;

  
  my @other_rf_ids = @{$self->adaptor->_fetch_other_dbIDs_by_stable_feature_set_ids
						 ($stable_id, 
							 \@fset_ids, 
						  { include_projected => $params_hash->{include_projected}} )};
  
  return (@other_rf_ids) ? 0 : 1;
}



=head2 get_other_RegulatoryFeatures

  Arg[1]     : optional - ARRAYREF of regualtory Bio::EnsEMBL::Funcgen::FeatureSet objects
                          Default is FeatureSet of given RegulatoryFeature, else need to be 
                          defined explicitly.
  Arg[2]     : optional - HASHREF Params hash: 
                                    {
                                     include_projected => 0|1, # Boolean, include 'projected' features
                                     include_multicell => 0|1, # Boolean, include MultiCell features
                                    }
  Example    : my @other_fsets = @{$reg_feat->get_other_FeatureSets($fsets)};
  Description: Gets other RegualtoryFeatures (linked via the stable ID) which are present in the 
               specified list of FeatureSets.
  Returntype : ARRAYREF of Bio::EnsEMBL::Funcgen::RegulatoryFeature objects
  Exceptions : Throw is arguments not stored or valid.
  Caller     : General
  Status     : At risk

=cut

sub get_other_RegulatoryFeatures{
  my ($self, $fsets, $params_hash) = @_;
  
  #define to avoid deref fails below.
  $params_hash ||= {};
  if(ref($params_hash) ne 'HASH'){
	throw("The params hash argument must be a valid HASHREF:\t".ref($params_hash));
  }

  $fsets ||= [$self->feature_set];
  my @fset_ids;
  
  foreach my $fset(@$fsets){
	#assume we have an adaptor set
	$self->adaptor->db->is_stored_and_valid('Bio::EnsEMBL::Funcgen::FeatureSet', $fset);
	
	if($fset->feature_class ne 'regulatory'){
	  throw('Found non-regulatory FeatureSet');
	}

	push @fset_ids, $fset->dbID;
  }
    
  my $stable_id;
  ($stable_id = $self->stable_id) =~ s/^[A-Z0]+//;

  my @other_fsets_ids = @{$self->adaptor->_fetch_other_dbIDs_by_stable_feature_set_ids
							($stable_id, \@fset_ids, 
							 {
							  include_projected => $params_hash->{include_projected},
							  include_multicell => $params_hash->{include_multicell},
							 })};

  return $self->adaptor->fetch_all_by_dbID_list(\@other_fsets_ids);
}



1;

__END__
