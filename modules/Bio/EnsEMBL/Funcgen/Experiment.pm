=head1 LICENSE

Copyright [1999-2013] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute

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
  developers list at <ensembl-dev@ebi.ac.uk>.

  Questions may also be sent to the Ensembl help desk at
  <helpdesk@ensembl.org>.


=head1 NAME

Bio::EnsEMBL::Funcgen::Experiment

=head1 SYNOPSIS

use Bio::EnsEMBL::Funcgen::Experiment;

my $exp = Bio::EnsEMBL::Funcgen::Experiment->new
               (
				-ADAPTOR             => $self,
				-NAME                => $name,
				-EXPERIMENTAL_GROUP  => $experimental_group,
				-DATE                => $date,
				-PRIMARY_DESIGN_TYPE => 'binding_site_indentification',
				-DESCRIPTION         => $description,
				-ARCHIVE_ID          => $archive_id,
               );

my $db_adaptor = Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor->new(...);
my $exp_adaptor = $db_adaptor->get_ExperimentAdaptor();
my $exp = $exp_adaptor->fetch_by_name($exp_name)

=head1 DESCRIPTION

The Experiment class represents an instance of an experiment i.e. a discrete set of data

=cut


################################################################################

package Bio::EnsEMBL::Funcgen::Experiment;

use warnings;
use strict;

use Bio::EnsEMBL::Utils::Argument  qw( rearrange );
use Bio::EnsEMBL::Utils::Exception qw( throw  );
use base qw( Bio::EnsEMBL::Funcgen::Storable );


=head2 new

  Arg [-NAME]                : String - experiment name
  Arg [-EXPERIMENTAL_GROUP]  : Bio::EnsEMBL::Funcgen ExperimentalGroup associated with this experiment
  Arg [-DATE]                : String - Date of the experiment (YYYY-MM-DD)
  Arg [-PRIMARY_DESIGN_TYPE] : String - MGED term for the primary design of teh experiment e.g. binding_site_identification
  Arg [-DESCRIPTION]         : String

  Example    : my $array = Bio::EnsEMBL::Funcgen::Experiment->new
                              (
							   -NAME                => $name,
							   -EXPERIMENTAL_GROUP  => $group,
							   -DATE                => $date,
							   -PRIMARY_DESIGN_TYPE => $p_design_type,
							   -DESCRIPTION         => $description,
					          );

  Description: Creates a new Bio::EnsEMBL::Funcgen::Experiment object.
  Returntype : Bio::EnsEMBL::Funcgen::Experiment
  Exceptions : Throws if name not defined
               Throws if ExperimentalGroup not valid
  Caller     : General
  Status     : Medium Risk

=cut

sub new {
	my $caller = shift;
	my $class  = ref($caller) || $caller;
	my $self   = $class->SUPER::new(@_);

	my ($name, $group, $date, $p_dtype, $desc, $archive_id, $data_url, $xml_id, $xml)
		= rearrange( ['NAME', 'EXPERIMENTAL_GROUP', 'DATE', 'PRIMARY_DESIGN_TYPE',
					  'DESCRIPTION','ARCHIVE_ID', 'DATA_URL', 'MAGE_XML', 'MAGE_XML_ID'], @_ );


    #Added in v68
    #Remove in v69
    if($data_url || $archive_id){
      throw('The -data_url and -archive_id parameters have been moved to the InputSubSet class');
    }


	#Mandatory attr checks

	if(ref($group) ne 'Bio::EnsEMBL::Funcgen::ExperimentalGroup'){
	  throw("Must pass a valid stored Bio::EnsEMBL::Funcgen::ExperimentalGroup object");
	}

	if(! defined $name){
	  throw('You must provide a name parameter');
	}

	#test date format here?


	#Direct assignment here so we avoid setter test in methods
	$self->{name}                = $name;
	$self->{group}               = $group;
	$self->{date}                = $date       if defined $date;
	$self->{primary_design_type} = $p_dtype    if defined $p_dtype; #MGED term for primary design type
	$self->{description}         = $desc       if defined $desc;

	#Maintain setter funcs here as these are populated after initialisation
	$self->mage_xml_id($xml_id) if defined $xml_id;
	$self->mage_xml($xml)       if defined $xml;

	return $self;
}


=head2 name

  Example     : my $exp_name = $exp->name;
  Description : Getter for the experiment name
  Returntype  : String
  Exceptions  : None
  Caller      : General
  Status      : Stable

=cut

sub name{
  return shift->{name};
}


=head2 experimental_group

  Example     : my $exp_group_name = $exp->experimental_group()->name();
  Description : Getter for the experimental group
  Returntype  : Bio::EnsEMBL::Funcgen::ExperimentalGroup
  Exceptions  : None
  Caller      : General
  Status      : At risk - to be deprecated

=cut

#change in webcode before deprecating

sub experimental_group{
  return shift->{group};
}


=head2 get_ExperimentalGroup

  Example     : my $exp_group_name = $exp->experimental_group()->name();
  Description : Getter for the experimental group
  Returntype  : Bio::EnsEMBL::Funcgen::ExperimentalGroup
  Exceptions  : None
  Caller      : General
  Status      : Stable

=cut

sub get_ExperimentalGroup{ return shift->{group}; }


=head2 date

  Example     : my $exp_date = $exp->date;
  Description : Getter for the date
  Returntype  : String
  Exceptions  : None
  Caller      : General
  Status      : Stable

=cut

sub date{
  return shift->{date};
}


=head2 description

  Example     : my $exp_desc = $exp->description
  Description : Getter for the experiment description
  Returntype  : String
  Exceptions  : None
  Caller      : General
  Status      : At risk - Not used, was stable until v64

=cut

sub description{
  return shift->{description};
}


=head2 primary_design_type

  Example     : my $pdt = $exp->primary_design_type;
  Description : Getter for the primary design type
  Returntype  : String - MGED term
  Exceptions  : None
  Caller      : General
  Status      : At risk

=cut

sub primary_design_type{
  return shift->{primary_design_type};
}


=head2 mage_xml

  Arg [1]     : string(optional) - MAGE XML
  Example     : my $xml = $exp->mage_xml();
  Description : Getter/Setter for the mage_xml attribute
  Returntype  : String
  Exceptions  : None
  Caller      : General
  Status      : at risk

=cut

sub mage_xml{
  my $self          = shift;	
  $self->{mage_xml} = shift if @_;

  if(! exists $self->{mage_xml} && $self->mage_xml_id()){
	$self->{mage_xml} = $self->adaptor->fetch_mage_xml_by_Experiment($self);
  }

  return (exists $self->{'mage_xml'}) ? $self->{'mage_xml'} : undef;
}


=head2 mage_xml_id

  Arg [1]     : int (optional) - mage_xml_id
  Example     : $exp->group_db_id('1');
  Description : Getter/Setter for the mage_xml attribute
  Returntype  : String
  Exceptions  : None
  Caller      : General
  Status      : at risk

=cut

sub mage_xml_id{
  my $self             = shift;	
  $self->{mage_xml_id} = shift if @_;
  return $self->{mage_xml_id};
}


=head2 get_InputSubsets

  Example     : my @issets = @{$exp->getInputSubsets()}
  Description : Retrieves all InputSubsets associated with this experiment.
  Returntype  : Arrayref of Bio::EnsEMBL::Funcgen::InputSubsets
  Exceptions  : None
  Caller      : General
  Status      : At risk

=cut

sub get_InputSubsets{
  my $self = shift;
 
  if(! exists $self->{'input_subsets'}){
    $self->{'input_subsets'} = {};
  
    foreach my $isset(@{$self->adaptor->db->get_InputSubsetAdaptor->
                          fetch_all_by_Experiments([$self])}){
      $self->{'input_subsets'}->{$isset->dbID} = $isset;
    }
  }

  return [values %{$self->{'input_subsets'}}];
}


=head2 get_ExperimentalChips

  Example     : my $exp_chips = @{$exp->get_ExperimentalChips()}
  Description : Retrieves all ExperiemntalChips
  Returntype  : Listref of ExperimentalChips
  Exceptions  : None
  Caller      : General
  Status      : At risk

=cut

sub get_ExperimentalChips{
  my $self = shift;

  if(! exists $self->{experimental_chips}){
    $self->{experimental_chips} = {};
  
    foreach my $echip(@{$self->adaptor->db->get_ExperimentalChipAdaptor->
                          fetch_all_by_experiment_dbID($self->dbID())}){
      $self->{experimental_chips}->{$echip->unique_id()} = $echip;
    }
  }

  return [values %{$self->{experimental_chips}}];
}

=head2 add_ExperimentalChip

  Example     : $exp_chip = $exp->add_ExperimentalChip($exp_chip)
  Description : Adds and stores an ExperiemntalChip for this Experiment
  Returntype  : Bio::EnsEMBL::Funcgen::ExperimentalChip
  Exceptions  : Throws is not passed a valid stored Bio::EnsENBML::Funcgen::ExperimentalChip
  Caller      : General
  Status      : At risk

=cut

sub add_ExperimentalChip{
  my $self  = shift; 
  my $echip = shift;
 

 throw("Must pass a valid stored Bio::EnsEMBL::Funcgen::ExperimentalChip object")
    if(! $echip->isa("Bio::EnsEMBL::Funcgen::ExperimentalChip") || ! $echip->dbID());

  if(! exists $self->{'experimental_chips'}){
    $self->get_ExperimentalChips();
    $self->{'experimental_chips'}->{$echip->unique_id()} = $echip;
    #do this here without checking to avoid probelm of retrieving first stored chip
  }elsif(exists  $self->{'experimental_chips'}->{$echip->unique_id()}){
    warn("You cannot add the same ExperimentalChip(".$echip->unique_id().")to an Experiment more than once, check your code");
  }else{
    $self->{'experimental_chips'}->{$echip->unique_id()} = $echip;
  }

  return;
}

=head2 get_ExperimentalChip_by_unique_id

  Example     : $exp_chip = $exp->add_ExperimentalChip($exp_chip)
  Description : Adds and stores an ExperiemntalChip for this Experiment
  Returntype  : Bio::EnsEMBL::Funcgen::ExperimentalChip
  Exceptions  : Throws if no uid supplied
  Caller      : General
  Status      : At risk

=cut

sub get_ExperimentalChip_by_unique_id{
  my $self = shift;
  my $uid  = shift;
  
  my ($echip);

  throw("Must supply a ExperimentalChip unque_id") if(! defined $uid);

  $self->{'experimental_chips'} || $self->get_ExperimentalChips();

  if(exists $self->{'experimental_chips'}->{$uid}){
    $echip = $self->{'experimental_chips'}->{$uid};
  }
  #should we warn here if not exists?

  return $echip;
}


=head2 get_ExperimentalChip_unique_ids

  Example     : foreach my $uid(@{$self->experiment->get_ExperimentalChip_unique_ids()}){ ... }
  Description : retrieves all ExperimentalChip unique_ids
  Returntype  : ListRef
  Exceptions  : None
  Caller      : General
  Status      : At risk

=cut

sub get_ExperimentalChip_unique_ids{
  my $self = shift;

  $self->{'experimental_chips'} || $self->get_ExperimentalChips();

  return [keys %{ $self->{'experimental_chips'}}];
}




=head2 reset_relational_attributes

  Arg[1] : Hashref containing the following mandatory parameters:
            -experimental_group => Bio::EnsEMBL::ExperimentalGroup

  Description: Resets all the relational attributes of a given Experiment
               Useful when creating a cloned object for migration beween DBs
  Returntype : None
  Exceptions : Throws if any of the parameters are not defined or invalid.
  Caller     : Migration code
  Status     : At risk

=cut

sub reset_relational_attributes{
  my ($self, $params_hash, $no_db_reset) = @_;

  my ($experimental_group) = rearrange(['EXPERIMENTAL_GROUP'], %$params_hash);

  #is_stored (in corresponding db) checks will be done in store method

  if(! (defined $experimental_group &&
        ref($experimental_group) eq 'Bio::EnsEMBL::Funcgen::ExperimentalGroup') ){
    my $msg = 'You must pass a valid Bio::EnsEMBL::Funcgen::ExperimentalGroup ';
    $msg .= 'not ' . ref($experimental_group);
    throw($msg);
  }

  $self->{group}    = $experimental_group;

  #Undef the dbID and adaptor by default
  if(! $no_db_reset){
    $self->{adaptor} = undef;
    $self->{dbID}    = undef;
  }

  return;
}


=head2 compare_to

Args[1]    : Bio::EnsEMBL::Funcgen::Storable (mandatory)
Args[2]    : Boolean - Optional 'shallow' - no object methods compared
Args[3]    : Arrayref - Optional list of Experiment method names each
             returning a Scalar or an Array or Arrayref of Scalars.
             Defaults to: name date primary_design_type description mage_xml_id
Args[4]    : Arrayref - Optional list of Experiment method names each
             returning a Storable or an Array or Arrayref of Storables.
             Defaults to: experimental_group
Example    : my %shallow_diffs = %{$rset->compare_to($other_rset, 1)};
Description: Compare this Experiment to another based on the defined scalar
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

  $scl_methods ||= [qw(name date primary_design_type description mage_xml_id)];
  $obj_methods ||= [qw(experimental_group)];

  return $self->SUPER::compare_to($obj, $shallow, $scl_methods,
                                  $obj_methods);
}

1;

