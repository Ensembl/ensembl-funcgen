#
# Ensembl module for Bio::EnsEMBL::DBSQL::Funcgen::AnnotatedFeatureAdaptor
#
# You may distribute this module under the same terms as Perl itself

=head1 NAME

Bio::EnsEMBL::DBSQL::Funcgen::AnnotatedFeatureAdaptor - A database adaptor for fetching and
storing AnnotatedFeature objects.

=head1 SYNOPSIS

my $afa = $db->get_AnnotatedFeatureAdaptor();

my $features = $afa->fetch_all_by_Slice($slice);
$features = $afa->fetch_all_by_Slice_Target($slice, $target);

=head1 DESCRIPTION

The AnnotatedFeatureAdaptor is a database adaptor for storing and retrieving
AnnotatedFeature objects.

=head1 AUTHOR

This module was created by Nathan Johnson.

This module is part of the Ensembl project: http://www.ensembl.org/

=head1 CONTACT

Post comments or questions to the Ensembl development list: ensembl-dev@ebi.ac.uk

=head1 METHODS

=cut

use strict;
use warnings;

package Bio::EnsEMBL::Funcgen::DBSQL::AnnotatedFeatureAdaptor;

use Bio::EnsEMBL::Utils::Exception qw( throw warning );
use Bio::EnsEMBL::Funcgen::AnnotatedFeature;
use Bio::EnsEMBL::Funcgen::DBSQL::BaseFeatureAdaptor;

use vars qw(@ISA);
@ISA = qw(Bio::EnsEMBL::Funcgen::DBSQL::BaseFeatureAdaptor);# Bio::EnsEMBL::Funcgen::DBSQL::BaseAdaptor);
#why is baseadaptor req'd? can we remove this?


#fetch_all_by_Slice_Experiment
#fetch_all_by_Slice_Analysis handled by BaseFeatureAdaptor fetch_all_by_Slice_score, which also takes logic name
#make wrapper method for useability?

#re-write the methods to take names aswell as objects? 


=head2 fetch_all_by_Slice_FeatureType

  Arg [1]    : Bio::EnsEMBL::Slice
  Arg [2]    : Bio::EnsEMBL::FeatureType
  Arg [3]    : (optional) string - analysis logic name
  Example    : my $slice = $sa->fetch_by_region('chromosome', '1');
               my $features = $ofa->fetch_by_Slice_FeatureType($slice, $ft);
  Description: Retrieves a list of features on a given slice, specific for a given FeatureType.
  Returntype : Listref of Bio::EnsEMBL::AnnotatedFeature objects
  Exceptions : Throws if no FeatureType object provided
  Caller     : General
  Status     : At Risk

=cut

sub fetch_all_by_Slice_FeatureType {
  my ($self, $slice, $type, $logic_name) = @_;
	

  throw('Need type as parameter') if ! $type->isa("Bio::EnsEMBL::Funcgen::FeatureType");
  my $ft_id = $type->dbID();
  
  my $constraint = qq( af.feature_set_id = fs.feature_set_id AND fs.feature_type_id = '$ft_id');

  $constraint = $self->_logic_name_to_constraint($constraint, $logic_name);

  return $self->SUPER::fetch_all_by_Slice_constraint($slice, $constraint);
}


=head2 fetch_all_by_Slice_FeatureSet

  Arg [1]    : Bio::EnsEMBL::Slice
  Arg [2]    : Bio::EnsEMBL::FeatureType
  Arg [3]    : (optional) string - analysis logic name
  Example    : my $slice = $sa->fetch_by_region('chromosome', '1');
               my $features = $ofa->fetch_by_Slice_FeatureType($slice, $ft);
  Description: Retrieves a list of features on a given slice, specific for a given FeatureType.
  Returntype : Listref of Bio::EnsEMBL::AnnotatedFeature objects
  Exceptions : Throws if no FeatureType object provided
  Caller     : General
  Status     : At Risk

=cut

sub fetch_all_by_Slice_FeatureSet {
  my ($self, $slice, $fset) = @_;
	

  throw('Need type as parameter') if ! $fset->isa("Bio::EnsEMBL::Funcgen::FeatureSet");
  my $fs_id = $fset->dbID();
  
  my $constraint = qq( af.feature_set_id =$fs_id );

  #could have individual logic_names for each annotated feature here?
  #$constraint = $self->_logic_name_to_constraint($constraint, $logic_name);

  return $self->SUPER::fetch_all_by_Slice_constraint($slice, $constraint);
}

=head2 fetch_all_by_Slice_FeatureSets

  Arg [1]    : Bio::EnsEMBL::Slice
  Arg [2]    : listref Bio::EnsEMBL::FeatureSet objects
  Example    : my $slice = $sa->fetch_by_region('chromosome', '1');
               my $features = $ofa->fetch_by_Slice_FeatureType($slice, $ft);
  Description: Retrieves a list of features on a given slice, specific for a given FeatureType.
  Returntype : Listref of Bio::EnsEMBL::AnnotatedFeature objects
  Exceptions : Throws if no FeatureType object provided
  Caller     : General
  Status     : At Risk

=cut

sub fetch_all_by_Slice_FeatureSets {
  my ($self, $slice, $fsets) = @_;
	
  my @fs_ids;
  foreach my $fset (@$fsets) {
	  throw('Not a FeatureSet object') 
		  if ! (ref($fset) && $fset->isa("Bio::EnsEMBL::Funcgen::FeatureSet"));
	  push (@fs_ids, $fset->dbID());
  }

  my $fs_ids = join(',', @fs_ids);
  my $constraint = qq( af.feature_set_id IN ($fs_ids) );

  #could have individual logic_names for each annotated feature here?
  #$constraint = $self->_logic_name_to_constraint($constraint, $logic_name);

  return $self->SUPER::fetch_all_by_Slice_constraint($slice, $constraint);
}



# Redefine BaseFeatureAdaptor method as analysis now abstracted to feature_set
# Given a logic name and an existing constraint this will
# add an analysis table constraint to the feature.  Note that if no
# analysis_id exists in the columns of the primary table then no
# constraint is added at all
# DO WE HAVE TO CALL THIS EXPLICITLY in this adaptor, or will BaseAdaptor use redefined method?





sub _logic_name_to_constraint {
  my $self = shift;
  my $constraint = shift;
  my $logic_name = shift;

  

  return $constraint if (!$logic_name);


  #make sure that an analysis_id exists in the primary table
  #my ($prim_tab) = $self->_tables();
  #my $prim_synonym = $prim_tab->[1];

  #my $found_analysis=0;
  #foreach my $col ($self->_columns) {
  #  my ($syn,$col_name) = split(/\./,$col);
  #  next if($syn ne $prim_synonym);
  #  if($col_name eq 'analysis_id') {
  #    $found_analysis = 1;
  #    last;
  #  }
  #}

  #if(!$found_analysis) {
  #  warning("This feature is not associated with an analysis.\n" .
  #          "Ignoring logic_name argument = [$logic_name].\n");
  #  return $constraint;
  #}

  my $aa = $self->db->get_AnalysisAdaptor();
  my $an = $aa->fetch_by_logic_name($logic_name);

  if(! $an) {
	  #warn or throw?
	  warn("No analysis associated with logic_name $logic_name");
	  return undef;
  }

  my $an_id = $an->dbID();

  $constraint .= ' AND' if($constraint);
  $constraint .= " feature_set.analysis_id = $an_id";
  return $constraint;
}



=head2 fetch_all_by_Slice_Experiment

  Arg [1]    : Bio::EnsEMBL::Slice
  Arg [2]    : Bio::EnsEMBL::Funcgen::Experiment
  Arg [3]    : (optional) string - logic name
  Example    : my $slice = $sa->fetch_by_region('chromosome', '1');
               my $features = $ofa->fetch_by_Slice_experiment_id($slice, $exp);
  Description: Retrieves a list of features on a given slice, specific for a given experiment.
  Returntype : Listref of Bio::EnsEMBL::AnnotatedFeature objects
  Exceptions : Throws if no Experiment object defined
  Caller     : General
  Status     : At Risk

=cut

sub fetch_all_by_Slice_Experiment {
  my ($self, $slice, $exp, $logic_name) = @_;
	

  if (! ($exp && $exp->isa("Bio::EnsEMBL::Funcgen::Experiment"))){
    throw("Need to pass a valid Bio::EnsEMBL::Funcgen::Experiment");
  }
  

  throw("Need to write DataSet and ResultSet to track back to experiment");
  
  
  my $constraint = qq( af.feature_set_id = fs.feature_set_id AND ep.experiment_id='$exp->dbID()' );
  
  $constraint = $self->_logic_name_to_constraint($constraint, $logic_name);
  
  return $self->SUPER::fetch_all_by_Slice_constraint($slice, $constraint, $logic_name);
}


=head2 _tables

  Args       : None
  Example    : None
  Description: PROTECTED implementation of superclass abstract method.
               Returns the names and aliases of the tables to use for queries.
  Returntype : List of listrefs of strings
  Exceptions : None
  Caller     : Internal
  Status     : At Risk

=cut

sub _tables {
  my $self = shift;
	
  return (
		  [ 'annotated_feature', 'af' ],
		  [ 'feature_set', 'fs'],#this is required for fetching on analysis or cell_type or feature_type
		 );
			#target?
}

=head2 _columns

  Args       : None
  Example    : None
  Description: PROTECTED implementation of superclass abstract method.
               Returns a list of columns to use for queries.
  Returntype : List of strings
  Exceptions : None
  Caller     : Internal
  Status     : At Risk

=cut

sub _columns {
  my $self = shift;
  
  return qw(
			af.annotated_feature_id  af.seq_region_id
			af.seq_region_start      af.seq_region_end
			af.seq_region_strand     af.coord_system_id
			af.feature_set_id        af.display_label
			af.score
	   );
}


#change this to experiment_prediction query?
#do we need this?

=head2 _default_where_clause

  Args       : None
  Example    : None
  Description: PROTECTED implementation of superclass abstract method.
               Returns an additional table joining constraint to use for
			   queries.
  Returntype : List of strings
  Exceptions : None
  Caller     : Internal
  Status     : At Risk

=cut

sub _default_where_clause {
  my $self = shift;
	
  #return 'af.feature_set_id = fs.feature_set_id';

  return;
}

=head2 _final_clause

  Args       : None
  Example    : None
  Description: PROTECTED implementation of superclass abstract method.
               Returns an ORDER BY clause. Sorting by oligo_feature_id would be
			   enough to eliminate duplicates, but sorting by location might
			   make fetching features on a slice faster.
  Returntype : String
  Exceptions : None
  Caller     : generic_fetch
  Status     : At Risk

=cut


sub _final_clause {
  return ' ORDER BY af.seq_region_id, af.seq_region_start, af.annotated_feature_id';
}


=head2 _objs_from_sth

  Arg [1]    : DBI statement handle object
  Example    : None
  Description: PROTECTED implementation of superclass abstract method.
               Creates AnnotatedFeature objects from an executed DBI statement
			   handle.
  Returntype : Listref of Bio::EnsEMBL::AnnotatedFeature objects
  Exceptions : None
  Caller     : Internal
  Status     : At Risk

=cut

sub _objs_from_sth {
	my ($self, $sth, $mapper, $dest_slice) = @_;


	#For EFG this has to use a dest_slice from core/dnaDB whether specified or not.
	#So if it not defined then we need to generate one derived from the species_name and schema_build of the feature we're retrieving.

	# This code is ugly because caching is used to improve speed

	#my $sa = $self->db->get_SliceAdaptor();

	
	my ($sa, $old_cs_id);
	$sa = $dest_slice->adaptor->db->get_SliceAdaptor() if($dest_slice);#don't really need this if we're using DNADBSliceAdaptor?

	#Some of this in now probably overkill as we'll always be using the DNADB as the slice DB
	#Hence it should always be on the same coord system
	#my $aa = $self->db->get_AnalysisAdaptor();
	#my $ct_adaptor = $self->db->get_CellTypeAdaptor();
	my $fset_adaptor = $self->db->get_FeatureSetAdaptor();
	my @features;
	my (%fset_hash, %slice_hash, %sr_name_hash, %sr_cs_hash);

	my (
	    $annotated_feature_id,  $seq_region_id,
	    $seq_region_start,      $seq_region_end,
	    $seq_region_strand,     $cs_id,
	    $fset_id,               $display_label,
	    $score
	);

	$sth->bind_columns(
					   \$annotated_feature_id,  \$seq_region_id,
					   \$seq_region_start,  \$seq_region_end,
					   \$seq_region_strand, \$cs_id,
					   \$fset_id,        \$display_label,
					   \$score
					  );


	#This needs doing properly!!
	#my $epsth = $self->prepare("SELECT experiment_id 
    #                                FROM experiment_prediction 
    #                                WHERE annotated_feature_id = ?");

	my $asm_cs;
	my $cmp_cs;
	my $asm_cs_name;
	my $asm_cs_vers;
	my $cmp_cs_name;
	my $cmp_cs_vers;
	if ($mapper) {
		$asm_cs      = $mapper->assembled_CoordSystem();
		$cmp_cs      = $mapper->component_CoordSystem();
		$asm_cs_name = $asm_cs->name();
		$asm_cs_vers = $asm_cs->version();
		$cmp_cs_name = $cmp_cs->name();
		$cmp_cs_vers = $cmp_cs->version();
	}

	my $dest_slice_start;
	my $dest_slice_end;
	my $dest_slice_strand;
	my $dest_slice_length;
	my $dest_slice_sr_name;
	if ($dest_slice) {
		$dest_slice_start   = $dest_slice->start();
		$dest_slice_end     = $dest_slice->end();
		$dest_slice_strand  = $dest_slice->strand();
		$dest_slice_length  = $dest_slice->length();
		$dest_slice_sr_name = $dest_slice->seq_region_name();
	}


	my $last_feature_id = -1;
	FEATURE: while ( $sth->fetch() ) {

	    #Need to build a slice adaptor cache here?
	    #Would only ever want to do this if we enable mapping between assemblies??
	    #Or if we supported the mapping between cs systems for a given schema_build, which would have to be handled by the core api
	    
		#get core seq_region_id
		$seq_region_id = $self->get_core_seq_region_id($seq_region_id);
		

	    
	    if($old_cs_id && ($old_cs_id != $cs_id)){
	      throw("More than one coord_system for feature query, need to implement SliceAdaptor hash?");
	    }
	    
	    $old_cs_id = $cs_id;
	    
	    
	    #Need to make sure we are restricting calls to Experiment and channel(i.e. the same coord_system_id)
	    
	    $sa ||= $self->db->get_SliceAdaptor($cs_id);
	    
	    
	    
	    # This assumes that features come out sorted by ID
	    next if ($last_feature_id == $annotated_feature_id);
	    $last_feature_id = $annotated_feature_id;
	    
	    # Get the analysis object
		#$analysis_hash{$analysis_id} = $aa->fetch_by_dbID($analysis_id) if(! exists $analysis_hash{$analysis_id});
		
		#possiblity of circular reference here?
		#Get the FeatureSet object
		$fset_hash{$fset_id} = $fset_adaptor->fetch_by_dbID($fset_id) if(! exists $fset_hash{$fset_id});

   
	    # Get the slice object
	    my $slice = $slice_hash{'ID:'.$seq_region_id};
	    
	    if (!$slice) {
	      $slice                              = $sa->fetch_by_seq_region_id($seq_region_id);
	      $slice_hash{'ID:'.$seq_region_id} = $slice;
	      $sr_name_hash{$seq_region_id}     = $slice->seq_region_name();
	      $sr_cs_hash{$seq_region_id}       = $slice->coord_system();
	    }
	    
	    my $sr_name = $sr_name_hash{$seq_region_id};
	    my $sr_cs   = $sr_cs_hash{$seq_region_id};
	    
	    # Remap the feature coordinates to another coord system if a mapper was provided
	    if ($mapper) {
	      
	      throw("Not yet implmented mapper, check equals are Funcgen calls too!");
	      
	      ($sr_name, $seq_region_start, $seq_region_end, $seq_region_strand)
			= $mapper->fastmap($sr_name, $seq_region_start, $seq_region_end, $seq_region_strand, $sr_cs);
	      
	      # Skip features that map to gaps or coord system boundaries
	      next FEATURE if !defined $sr_name;
	      
	      # Get a slice in the coord system we just mapped to
	      if ( $asm_cs == $sr_cs || ( $cmp_cs != $sr_cs && $asm_cs->equals($sr_cs) ) ) {
		$slice = $slice_hash{"NAME:$sr_name:$cmp_cs_name:$cmp_cs_vers"}
		  ||= $sa->fetch_by_region($cmp_cs_name, $sr_name, undef, undef, undef, $cmp_cs_vers);
	      } else {
		$slice = $slice_hash{"NAME:$sr_name:$asm_cs_name:$asm_cs_vers"}
		  ||= $sa->fetch_by_region($asm_cs_name, $sr_name, undef, undef, undef, $asm_cs_vers);
	      }
	    }
	    
	    # If a destination slice was provided convert the coords
	    # If the destination slice starts at 1 and is forward strand, nothing needs doing
	    if ($dest_slice) {
	      unless ($dest_slice_start == 1 && $dest_slice_strand == 1) {
		if ($dest_slice_strand == 1) {
		  $seq_region_start = $seq_region_start - $dest_slice_start + 1;
		  $seq_region_end   = $seq_region_end   - $dest_slice_start + 1;
		} else {
		  my $tmp_seq_region_start = $seq_region_start;
		  $seq_region_start        = $dest_slice_end - $seq_region_end       + 1;
		  $seq_region_end          = $dest_slice_end - $tmp_seq_region_start + 1;
		  $seq_region_strand      *= -1;
		}
	      }
	      
	      # Throw away features off the end of the requested slice
	      next FEATURE if $seq_region_end < 1 || $seq_region_start > $dest_slice_length
		|| ( $dest_slice_sr_name ne $sr_name );
	      
	      $slice = $dest_slice;
	    }
	    

	    #my @exp_ids = map $_ = "@$_", @{$self->dbc->db_handle->selectall_arrayref($sql)};
	    #$epsth->bind_param(1, $annotated_feature_id,    SQL_INTEGER);
	    #$epsth->execute();
	    #my @exp_ids = map $_ = "@$_", @{$epsth->fetchall_arrayref()};


		#RegulatoryFeature hack
		my ($reg_type, $stable_id, $reg_attrs, $vector);
		
		#look at analysis to avoid set name issues
		if($fset_hash{$fset_id}->analysis->logic_name() eq 'RegulatoryRegion'){

		  ($stable_id, $vector) = split/:/, $display_label;
		  
		  $display_label = 'Regulatory Feature';
		  
		  #RegulatoryFeature hack



		  #We don't consider the non-epi feature bits as these are only used to
		  #cluster and build the patterns, not to assign a classification
		  #as this would prevent us from finding novel regions

		  my %reg_class_regexs = (
								  #'1....(10|01).'  => 'Gene end associated', 
								  #'1...1...'        => 'Promoter associated',#orig
								  '1...1.....'        => 'Promoter associated',
								  '1.0.001...' => 'Non-gene associated',
								  '11..01....' => 'Gene associated',
								 );



		  #omit TSS and TES from here?
		  my @reg_feature_attrs = ('DNase1', 'CTCF', 'H4K20me3', 'H3K27me3', 
								   'H3K36me3', 'H3K4me3', 'H3K79me3', 'H3K9me3', 'TSS Proximal', 'TES Proximal'); 
		  my @vector = split//, $vector;

		  foreach my $i(0..7){#$#vector){
			push @$reg_attrs, $reg_feature_attrs[$i] if $vector[$i];
		  }


		  foreach my $regex(keys %reg_class_regexs){

			if($vector =~ /$regex/){

			  #warn "$vector matches ".$reg_class_regexs{$regex}."\t$regex\n";

			  throw('Found non-mutually exclusive regexs') if $reg_type;
			  $reg_type = $reg_class_regexs{$regex};
			}
			
		  }

		  $reg_type ||= 'Unclassified';
		}

	
	    push @features, $self->_new_fast( {
										   'start'          => $seq_region_start,
										   'end'            => $seq_region_end,
										   'strand'         => $seq_region_strand,
										   'slice'          => $slice,
										   'analysis'       => $fset_hash{$fset_id}->analysis(),
										   'adaptor'        => $self,
										   'dbID'           => $annotated_feature_id,
										   'score'          => $score,
										   'display_label'  => $display_label,
										   'feature_set'    => $fset_hash{$fset_id},
										   'regulatory_type'=> $reg_type,
										   'regulatory_attributes' => $reg_attrs,
										   'stable_id'      => $stable_id,
										  } );
	}
	
	return \@features;
}

=head2 _new_fast

  Args       : Hashref to be passed to AnnotatedFeature->new_fast()
  Example    : None
  Description: Construct a AnnotatedFeature object using quick and dirty new_fast.
  Returntype : Bio::EnsEMBL::Funcgen::AnnotatedFeature
  Exceptions : None
  Caller     : _objs_from_sth
  Status     : At Risk

=cut

sub _new_fast {
	my $self = shift;
	
	my $hash_ref = shift;
	return Bio::EnsEMBL::Funcgen::AnnotatedFeature->new_fast($hash_ref);
}

=head2 store

  Args       : List of Bio::EnsEMBL::Funcgen::AnnotatedFeature objects
  Example    : $ofa->store(@features);
  Description: Stores given AnnotatedFeature objects in the database. Should only
               be called once per feature because no checks are made for
			   duplicates. Sets dbID and adaptor on the objects that it stores.
  Returntype : Listref of stored AnnotatedFeatures
  Exceptions : Throws if a list of AnnotatedFeature objects is not provided or if
               the Analysis, CellType and FeatureType objects are not attached or stored
  Caller     : General
  Status     : At Risk

=cut

sub store{
	my ($self, @pfs) = @_;
	
	if (scalar(@pfs) == 0) {
		throw('Must call store with a list of AnnotatedFeature objects');
	}
	
	my $sth = $self->prepare("
		INSERT INTO annotated_feature (
			seq_region_id,   seq_region_start,
			seq_region_end,  seq_region_strand,
                        coord_system_id, feature_set_id,
			display_label,   score
		) VALUES (?, ?, ?, ?, ?, ?, ?, ?)
	");
	
	#my $epsth = $self->prepare("INSERT INTO experiment_prediction (
	#experiment_id, annotated_feature_id)
    #                          VALUES (?, ?)");
	
	my $db = $self->db();
	#my $analysis_adaptor = $db->get_AnalysisAdaptor();
	
  FEATURE: foreach my $pf (@pfs) {
		
		if( !ref $pf || !$pf->isa('Bio::EnsEMBL::Funcgen::AnnotatedFeature') ) {
			throw('Feature must be an AnnotatedFeature object');
		}
		
		if ( $pf->is_stored($db) ) {
			#does not accomodate adding Feature to >1 feature_set
			warning('AnnotatedFeature [' . $pf->dbID() . '] is already stored in the database');
			next FEATURE;
		}
		
		
		
		#Have to do this for Analysis separately due to inheritance, removed defined as constrained in Feature->new
		#Redundancy with Analysis in FeatureSet
		if ( ! $pf->analysis->is_stored($db)) {
			throw('A stored Bio::EnsEMBL::Analysis must be attached to the AnnotatedFeature objects to be stored.');
		}
		
		if (! $pf->feature_set->is_stored($db)) {
			throw('A stored Bio::EnsEMBL::Funcgen::FeatureSet must be attached to the AnnotatedFeature objects to be stored.');
		}

		#sanity check analysis matches feature_set analysis
		if($pf->analysis->dbID() != $pf->feature_set->analysis->dbID()){
			throw("AnnotatedFeature analysis(".$pf->analysis->logic_name().") does not match FeatureSet analysis(".$pf->feature_set->analysis->logic_name().")\n".
				  "Cannot store mixed analysis sets");
		}
		#Complex analysis to be stored as one in analysis table, or have feature_set_prediciton link table?
		#Or only have single analysis feature which can contribute to multi analysis "regulons"
		#Or can we have multiple entries in feature_set with the same id but different analyses?
		#This would still not be specific for each feature, nor would the annotated_feature analysis_id
		#reflect all the combined analyses.  Maybe just the one which contributed most?

		# Store the analysis if it has not been stored yet
		#$analysis_adaptor->store( $pf->analysis()) if ( !$pf->analysis->is_stored($db) );
		#could this potentially store the same on multiple times?

		my $seq_region_id;
		($pf, $seq_region_id) = $self->_pre_store($pf);
		
		$sth->bind_param(1, $seq_region_id,             SQL_INTEGER);
		$sth->bind_param(2, $pf->start(),               SQL_INTEGER);
		$sth->bind_param(3, $pf->end(),                 SQL_INTEGER);
		$sth->bind_param(4, $pf->strand(),              SQL_TINYINT);
		$sth->bind_param(5, $pf->coord_system_id(),     SQL_INTEGER);
		$sth->bind_param(6, $pf->feature_set->dbID(),   SQL_INTEGER);
		$sth->bind_param(7, $pf->display_label(),       SQL_VARCHAR);
		$sth->bind_param(8, $pf->score(),                SQL_DOUBLE);
		
		$sth->execute();
		$pf->dbID( $sth->{'mysql_insertid'} );

		#foreach my $exp_id(@{$pf->experiment_ids()}){
		#  $epsth->bind_param(1, $exp_id);
		#  $epsth->bind_param(2, $original->dbID());
		#  $epsth->execute();
		#}
		
		$pf->adaptor($self);
	}

  return \@pfs;
}

#re-written for non-standard feature table

=head2 fetch_all_by_logic_name 

  Arg [1] : string $logic_name the logic name of the type of features to obtain 
  Example : $fs = $a->fetch_all_by_logic_name('foobar'); 
  Description: Returns a listref of features created from the database. only 
               features with an analysis of type $logic_name will be returned. 
  Returntype : listref of Bio::EnsEMBL::SeqFeatures 
  Exceptions : thrown if $logic_name 
  Caller : General 
  Status : Stable 

=cut 

sub fetch_all_by_logic_name { 
  my $self = shift; 
  my $logic_name = shift || throw( "Need a logic_name" ); 

  my $constraint;
  my $an = $self->db->get_AnalysisAdaptor->fetch_by_logic_name($logic_name); 
  $constraint = ' af.feature_set_id=fs.feature_set_id and fs.analysis_id = '.$an->dbID() if($an);

  return (defined $constraint) ? $self->generic_fetch($constraint) : undef;
} 




=head2 list_dbIDs

  Args       : None
  Example    : my @feature_ids = @{$ofa->list_dbIDs()};
  Description: Gets an array of internal IDs for all OligoFeature objects in
               the current database.
  Returntype : List of ints
  Exceptions : None
  Caller     : ?
  Status     : Medium Risk

=cut

sub list_dbIDs {
	my $self = shift;
	
	return $self->_list_dbIDs('annotated_feature');
}







1;
