#
# Ensembl module for Bio::EnsEMBL::Funcgen::FeatureType
#

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

=head1 NAME

Bio::EnsEMBL::Funcgen::FeatureType - A module to represent a FeatureType. i.e. the target of an experiment.

=head1 SYNOPSIS

use Bio::EnsEMBL::Funcgen::FeatureType;


=head1 DESCRIPTION

This is a simple class to represent information about a FeatureType, containing the name i.e Brno nomenclature or other controlled/validated name relevant to the class (HISTONE, PROMOTER etc), and description. This module is part of the Ensembl project: http://www.ensembl.org/

=cut

#To do
# add coding_transcript/gene methods.  Store as xrefs or custom feature_type_coding table? (miRanda etc)
#

package Bio::EnsEMBL::Funcgen::FeatureType;

use strict;
use warnings;
use Bio::EnsEMBL::Utils::Argument  qw( rearrange ) ;
use Bio::EnsEMBL::Utils::Exception qw( throw );

use base qw(Bio::EnsEMBL::Funcgen::Storable);


=head2 new

  Arg [-name]         : String - name of FeatureType
  Arg [-class]        : String - class of FeatureType
  Arg [-description]  : String - descriptiom of FeatureType
  Arg [-analysis]     : optional Bio::EnsEMBL::Analysis used to generate FeatureType
  Arg [-so_accession] : optional String - Sequence ontology accession
  Arg [-so_name]      : optional String - Sequence ontology name

  Example    : my $ft = Bio::EnsEMBL::Funcgen::FeatureType->new
                           (
                            -name  => "H3K9Me",
                            -class => "HISTONE",
                            -description => "Generalised methylation of Histone 3 Lysine 9",
                            -analysis => $analysis,
                            -so_name  => $so_name,
                            -so_accession => $so_accession
                           );
  Description: Constructor method for FeatureType class
  Returntype : Bio::EnsEMBL::Funcgen::FeatureType
  Exceptions : Throws if name or class not defined
               Throws if analysis is defined but not valid
  Caller     : General
  Status     : Stable

=cut

sub new {
  my $caller = shift;
  my $obj_class = ref($caller) || $caller;
  my $self = $obj_class->SUPER::new(@_);

  my ($name, $desc, $class, $analysis, $so_acc, $so_name) =
    rearrange(['NAME', 'DESCRIPTION', 'CLASS', 'ANALYSIS', 'SO_ACCESSION', 'SO_NAME'], @_);

  throw("Must supply a FeatureType name\n") if ! defined $name;
  throw("Must supply a FeatureType class\n") if ! defined $class;

  #Direct assignments here prevent set arg test in getter only method
  $self->{name}         = $name;
  $self->{class}        = $class;
  $self->{description}  = $desc    if defined $desc;
  $self->{so_name}      = $so_name if defined $so_name;
  $self->{so_accession} = $so_acc  if defined $so_acc;

  if($analysis){

    if(ref($analysis) ne 'Bio::EnsEMBL::Analysis'){
      throw('Optional Analysis parameter must be a valid Bio::EnsEMBL::Analysis');
      #is_stored checks done in other fetch and store methods
    }

    $self->{analysis} = $analysis;
  }

  return $self;
}


=head2 _creates_broad_peaks

  Returns:
    true, if the given feature_type creates_broad_peaks, 
    false otherwise.

=cut
sub _creates_broad_peaks {

  my $self = shift;
  my $feature_type_name = $self->name;

  my @broad_peak_feature_type_names = qw(
    H3K36me3
    H3K27me3
    H2AK5ac 
    H2BK12ac
    H3K14ac 
    H3K23me2
    H3K4me1 
    H3K79me1
    H3K79me2
    H3K9me1 
    H3K9me3 
    H4K20me1
    H4K8ac
  );
  
  my $is_broad_peak_feature_type 
    = grep { 
      $_ eq $feature_type_name 
    } @broad_peak_feature_type_names;
  
  if ($is_broad_peak_feature_type) {
    return 1;
  }
  return;
}

=head2 name

  Example    : my $name = $ft->name;
  Description: Getter of name attribute for FeatureType objects
  Returntype : String
  Exceptions : None
  Caller     : General
  Status     : Stable

=cut

sub name { return $_[0]->{name}; }

=head2 description

  Example    : my $desc = $ft->description;
  Description: Getter of description attribute for FeatureType objects.
  Returntype : String
  Exceptions : None
  Caller     : General
  Status     : Stable

=cut

sub description {  return $_[0]->{description}; }


=head2 class

  Example    : my $ft_class = $ft->class;
  Description: Getter of class attribute for FeatureType objects.
  Returntype : String
  Exceptions : None
  Caller     : General
  Status     : Stable

=cut

sub class{ return $_[0]->{class}; }

=head2 so_accession

  Example    : my $ft_class = $ft->class;
  Description: Getter of sequence ontoloy accession for FeatureType objects.
  Returntype : String
  Exceptions : None
  Caller     : General
  Status     : At Risk

=cut

sub so_accession{  return $_[0]->{so_accession}; }

=head2 so_name

  Example    : my $so_name = $ft->so_name;
  Description: Getter of sequence ontology name  for FeatureType objects.
  Returntype : String
  Exceptions : None
  Caller     : General
  Status     : At Risk

=cut

sub so_name{   return $_[0]->{so_name}; }

=head2 analysis

  Example    : my $ft_anal = $ft->analysis;
  Description: Getter of the Analysis for FeatureType objects.
  Returntype : Bio::EnsEMBL::Analysis
  Exceptions : None
  Caller     : General
  Status     : At risk

=cut

sub analysis{ return $_[0]->{analysis}; }


=head2 evidence_type_label

  Example    : my $track_label = $fsets[0]->feature_type->evidence_type_label.' MultiCell';
  Description: Getter for short evidence type label used in track label and field headers etc.
  Returntype : string
  Exceptions : None
  Caller     : Web code
  Status     : At risk

=cut

sub evidence_type_label{
  my $self = shift;

  if(! exists $self->{evidence_type_label}){
    $self->{evidence_type_label} =
      $self->adaptor->get_regulatory_evidence_info($self->regulatory_evidence_type)->{label};
  }

  return $self->{evidence_type_label};
}


=head2 evidence_type_name

  Example    : my $name = $fsets[0]->feature_type->evidence_type_name;
  Description: Getter for evidence type name used in browser.
  Returntype : string
  Exceptions : None
  Caller     : Web code
  Status     : At risk

=cut

sub evidence_type_name{
  my $self = shift;

  if(! exists $self->{evidence_type_name}){
    $self->{evidence_type_name} =
      $self->adaptor->get_regulatory_evidence_info($self->regulatory_evidence_type)->{name};
  }

  return $self->{evidence_type_name};
}


=head2 evidence_type_long_name

  Example    : my $long_name = $fsets[0]->feature_type->evidence_type_long_name;
  Description: Getter for evidence type name used in browser.
  Returntype : string
  Exceptions : None
  Caller     : Web code
  Status     : At risk

=cut

sub evidence_type_long_name{
  my $self = shift;

  if(! exists $self->{evidence_type_long_name}){
    $self->{evidence_type_long_name} =
      $self->adaptor->get_regulatory_evidence_info($self->regulatory_evidence_type)->{long_name};
  }

  return $self->{evidence_type_long_name};
}


=head2 is_core_evidence

  Example    : if($ftype->is_core_evidence){#this is a TFBS or DNase}
  Description: Returns true if this FeatureType is used to defin the core region
               of RegulatoryFeatures
  Returntype : Boolean
  Exceptions : None
  Caller     : Regulatory build
  Status     : At risk

=cut

sub is_core_evidence{
  my $self = $_[0];

  if(! defined $self->{is_core_evidence}){

    if($self->regulatory_evidence_type eq 'core'){
      $self->{is_core_evidence} = 1;
    }
    else{
      $self->{is_core_evidence} = 0;
    }
  }

  return $self->{is_core_evidence};
}

=head2 regulatory_evidence_type

  Example    : my $re_type = $ftype->regulatory_evidence_type
  Description: Returns the regulatory evidence type i.e. core or non_core
  Returntype : String
  Exceptions : None
  Caller     : General
  Status     : At risk

=cut

sub regulatory_evidence_type{
  my $self = $_[0];

  if(! defined $self->{regulatory_evidence_type}){
    $self->{regulatory_evidence_type} = $self->adaptor->get_regulatory_evidence_type($self->class);
  }

  return $self->{regulatory_evidence_type};
}

=head2 get_all_coding_gene_stable_ids

  Example    : my @gene_dbentries = @{ $storable->get_all_coding_gene_stable_ids };
  Description: Retrieves Ensembl Gene stable IDs (xrefs) for this FeatureType.
               This does _not_ include the corresponding translations
               DBEntries (see get_all_DBLinks).

               This method will attempt to lazy-load DBEntries from a
               database if an adaptor is available and no DBEntries are present
               on the transcript (i.e. they have not already been added or
               loaded).
  Returntype : Listref of strings (stable ids)
  Exceptions : none
  Caller     : general
  Status     : at risk

=cut

#Need to add optional Transcript/Gene param so we can filter
#Filter here or would be better to restrict in sql query ni DBEntryAdaptor?

sub get_all_coding_gene_stable_ids {
    my $self = shift;

    #We wouldn't need this if we made the xref schema multi species
    #my $species = $self->adaptor->db->species;
    my $species
        = Bio::EnsEMBL::Registry->get_alias( $self->adaptor->db->species );

    if ( !$species ) {
        throw(
            'You must specify a DBAdaptor -species to retrieve DBEntries based on the external_db.db_name'
        );
    }

    #safety in case we get Homo sapiens
    ( $species = lc($species) ) =~ s/ /_/;

    my $DBEntries = $self->_get_all_DBEntries( $species . '_core_Gene' );
    my @stable_ids;

    for my $DBEntry (@$DBEntries){
      push @stable_ids, $DBEntry->primary_id();
    }

    return \@stable_ids;
}

=head2 _get_all_DBEntries

  Arg[1]     : string - External DB name e.g. ensembl_core_Gene
  Arg[2]     : string - External DB type
  Example    : my @dbentries = @{ $set_feature->get_all_DBEntries };
  Description: Retrieves DBEntries (xrefs) for this SetFeature.
               This does _not_ include the corresponding translations
               DBEntries (see get_all_DBLinks).

               This method will attempt to lazy-load DBEntries from a
               database if an adaptor is available and no DBEntries are present
               on the SetFeature (i.e. they have not already been added or
               loaded).
  Returntype : Listref of Bio::EnsEMBL::DBEntry objects
  Exceptions : none
  Caller     : general, get_all_DBLinks
  Status     : Stable - at risk move to storable

=cut


#We could add 3rd arg here which would be xref(info_)type e.g. Gene/Transcript etc.
#Move info_type to ox.linkage_type to sit along side linkage_annotated


sub _get_all_DBEntries {
  my $self = shift;
  my $ex_db_exp = shift;
  my $ex_db_type = shift;

  my $cache_name = "dbentries";

  if(defined($ex_db_exp)){
    $cache_name .= $ex_db_exp;
  }
  if(defined($ex_db_type)){
    $cache_name .= $ex_db_type;
  }

  #Need to add tests for valid objects for xrefs

  # if not cached, retrieve all of the xrefs for this gene

  #This is not using the caching optimally
  #It seems for naive(ex_db_exp,ex_db_type) queries we create a naive cache
  #This means that further more specific queries will make another query and not use the cache


  if( (! defined $self->{$cache_name}) && $self->adaptor() ){

  my @tables = $self->adaptor->_tables;
  @tables = split/_/, $tables[0]->[0];#split annotated_feature
  my $object_type = join('', (map ucfirst($_), @tables));#change to AnnotatedFeature

    $self->{$cache_name} =
      $self->adaptor->db->get_DBEntryAdaptor->_fetch_by_object_type($self->dbID(), $object_type, $ex_db_exp, $ex_db_type);
  }
  elsif( ! defined $self->{$cache_name} ){
  throw('You must have set and adaptor to be able to get_all_DBEntries');
  }


  $self->{$cache_name} ||= [];

  return $self->{$cache_name};
}





=head2 compare

  Arg[1]     : Bio::EnsEMBL::Funcgen::FeatureType
               The analysis to compare to
  Example    : none
  Description: returns 1 if this FeatureType is the same
               returns 0 if there is a mistmatch apart from the dbID/DB/Adaptor
  Returntype : Boolean
  Exceptions : Throws if arg is not valid
  Caller     : General
  Status     : At risk

=cut

#simplified version of Analysis:compare
#move to storable and take method args

sub compare{
  my ($self, $ftype) = @_;

  if(ref($ftype) ne 'Bio::EnsEMBL::Funcgen::FeatureType'){
    throw('You must pass a valid Bio::EnsEMBL::Funcgen::FeatureType to compare');
  }

  my $same = 1;

  foreach my $methodName ( 'name', 'class', 'so_accession','so_name','description'){

    if( defined $self->$methodName() && ! $ftype->can($methodName )) {
      $same = 0;
	  last;
    }
    if( defined $self->$methodName() && ! defined $ftype->$methodName() ) {
      $same = 0;
	  last;
    }

    if( defined($ftype->$methodName()) && defined($self->$methodName()) &&
        ( $ftype->$methodName() ne $ftype->$methodName() )) {
      $same = 0;
	  last;
    }
  }


  #This would be in a wrapper method
  if($self->analysis && $ftype->analysis){

	if($self->analysis->compare($ftype->analysis)){
	  #analysis compare returns the opposite of what you expect
	  $same = 0;
	}
  }
  elsif( ! ((! $self->analysis) && (! $ftype->analysis)) ){#Only one has analysis
    $same = 0;
  }

  return $same;
}

=head2 reset_relational_attributes

  Arg[1]     : UNDEF - placeholder for compliance with similar methods in
               other modules
  Arg[2]     : Flag to avoid reseting of adaptor and dbID

  Description: Resets all the relational attributes of a given FeaturSet
               Useful when creating a cloned object for migration beween DBs
  Returntype : None
  Exceptions : Throws if any of the parameters are not defined or invalid.
  Caller     : Migration code
  Status     : At risk

=cut

sub reset_relational_attributes{
  my ($self, $params_hash, $no_db_reset) = @_;

  # Undef  dbID and adaptor by default
  if(! $no_db_reset){
    $self->{adaptor} = undef;
    $self->{dbID}    = undef;
  }

  # For future use. Analysis is optional, as it only is used at
  # Segmentation
  if(defined $self->{analysis}){
    my $msg = 'Linked Analysis found. This is usually the case with ';
    $msg   .= 'Segmenation data, for which reset has not been implemented';
    throw($msg);

    my ($analysis) = rearrange(['ANALYSIS',], @$params_hash);

    if(defined $analysis){
      if(! (ref($analysis) eq 'Bio::EnsEMBL::Analysis') ){
        throw('You must pass a valid Bio::EnsEMBL::Analysis');
      }
    }

    $self->{analysis}     = $analysis;
  }


  return;
}

=head2 compare_to

  Args[1]    : Bio::EnsEMBL::Funcgen::Storable (mandatory)
  Args[2]    : Boolean - Optional 'shallow' - no object methods compared
  Args[3]    : Arrayref - Optional list of FeatureType method names each
               returning a Scalar or an Array or Arrayref of Scalars.
               Defaults to: name class description so_accession so_name
  Args[4]    : Arrayref - Optional list of FeatureType method names each
               returning a Storable or an Array or Arrayref of Storables.
               Defaults to: analysis
  Example    : my %shallow_diffs = %{$rset->compare_to($other_rset, 1)};
  Description: Compare this FeatureType to another based on the defined scalar
               and storable methods.
  Returntype : Hashref of key attribute/method name keys and values which differ.
               Keys will always be the method which has been compared.
               Values can either be a error string, a hashref of diffs from a
               nested object, or an arrayref of error strings or hashrefs where
               a particular method returns more than one object.
  Exceptions : None
  Caller     : Import/migration pipeline
  Status     : At Risk

=cut

sub compare_to {
  my ($self, $obj, $shallow, $scl_methods, $obj_methods) = @_;

  $scl_methods ||= [qw(name class description so_accession so_name)];
  $obj_methods ||= [qw(analysis)];

  return $self->SUPER::compare_to($obj, $shallow, $scl_methods,
                                  $obj_methods);
}

1;

