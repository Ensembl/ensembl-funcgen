#
# Ensembl module for Bio::EnsEMBL::Funcgen::Epigenome
#


=head1 LICENSE

Copyright [1999-2016] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute

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

Bio::EnsEMBL::Funcgen::Epigenome - A module to represent a Epigenome.

=head1 SYNOPSIS

use Bio::EnsEMBL::Funcgen::Epigenome;

#Fetch from adaptor
my $epigenome = $epigenome_adaptor->fetch_by_name($epigenome_name);

#Create from new
my $epigenome = Bio::EnsEMBL::Funcgen::Epigenome->new
  (
   -name          => 'H1-ESC',
   -display_label => 'H1-ESC',
   -description   => 'Human Embryonic Stem Cell',
   -efo_id        => 'efo:EFO_0003042',
   -tissue        => 'embryonic stem cell',
  );

print $epigenome->name.' is a '.$epigenome->description."\n";

#H1-ESC is a Human Embryonic Stem Cell


=head1 DESCRIPTION

This is a simple class to represent information about an epigenome.  This may represent
harvested cells, a cell line or a more generic tissue type.

=cut

package Bio::EnsEMBL::Funcgen::Epigenome;

use strict;
use warnings;

use Bio::EnsEMBL::Utils::Argument  qw( rearrange ) ;
use Bio::EnsEMBL::Utils::Exception qw( throw deprecate );

use parent qw(Bio::EnsEMBL::Funcgen::Storable);

my %valid_genders = (male   => 1,
                     female => 1,
                     hermaphrodite => 1,
                     mixed => 1,
                     unknown => 1 
                    );

=head2 new

  Arg [1]    : String - name of Epigenome
  Arg [2]    : String (optional) - display label of Epigenome. Defaults to name
  Arg [3]    : String (optional) - description of Epigenome
  Arg [4]    : String (optional) - gender e.g. male, female or hermaphrodite
  Arg [5]    : String (optional) - Experimental Factor Ontology ID e.g. efo:EFO_0002869

  Example              : my $epigenome = Bio::EnsEMBL::Funcgen::Epigenome->new
                                    (
                                     -name          => "U2OS",
                                     -display_label => "U20S",
                                     -description   => "Human Bone Osteosarcoma Epithelial Cells",
                                     -gender        => 'female',
                                     -efo_id        => 'efo:EFO_0002869',
                                    );

  Description: Constructor method for Epigenome class
  Returntype : Bio::EnsEMBL::Funcgen::Epigenome
  Exceptions : Throws if name not defined or gender is invalid
  Caller     : General
  Status     : Stable

=cut

#-type/class => "TISSUE", enum? Mandatory.
#remove display label?

sub new {
  my $caller = shift;
  my $class = ref($caller) || $caller;
  my $self = $class->SUPER::new(@_);

  my ($name, $dlabel, $desc, $gender, $ontology_accession, $tissue, $production_name) = rearrange
    (['NAME', 'DISPLAY_LABEL', 'DESCRIPTION','GENDER', 'ONTOLOGY_ACCESSION', 'TISSUE', 'PRODUCTION_NAME'], @_);

  throw("Must supply an Epigenome name") if ! defined $name;

  if (defined $gender) {

    if ( ! exists $valid_genders{lc($gender)} ) { #enum will not force this so validate here
      throw("Gender $gender not valid, must be one of:\t".join(' ', keys %valid_genders));
    }

    $self->{gender} = $gender;
  }

  #Set explicitly to enable faster getter only methods
$self->{name}               = $name;
$self->{display_label}      = $dlabel || $name;
$self->{description}        = $desc if defined $desc;
$self->{ontology_accession} = $ontology_accession if defined $ontology_accession;
$self->{tissue}             = $tissue if defined $tissue;
$self->{production_name}    = $production_name;

  

  return $self;
}

=head2 name

  Example    : my $name = $epigenome->name;
  Description: Getter of name attribute for Epigenome objects
  Returntype : String
  Exceptions : None
  Caller     : General
  Status     : Stable

=cut

sub name { return $_[0]->{name}; }

sub production_name { return $_[0]->{production_name}; }

=head2 gender

  Example    : my $gender = $epigenome->gender();
  Description: Getter for the gender attribute for Epigenome objects
  Returntype : String
  Exceptions : None
  Caller     : General
  Status     : Stable

=cut

sub gender {  return $_[0]->{gender}; }


=head2 description

  Example    : my $desc = $epigenome->description();
  Description: Getter of description attribute for Epigenome objects
  Returntype : String
  Exceptions : None
  Caller     : General
  Status     : Stable

=cut

sub description {  return $_[0]->{description}; }


=head2 display_label

  Example    : my $display_label = $epigenome->display_label();
  Description: Getter of display_label attribute for Epigenome objects.
  Returntype : String
  Exceptions : None
  Caller     : General
  Status     : At Risk

=cut

sub display_label {  return $_[0]->{display_label}; }


=head2 efo_id

  Example    : my $efo_id = $epigenome->efo_id;
  Description: Getter of the Experimental Factor Ontology ID
  Returntype : String
  Exceptions : None
  Caller     : General
  Status     : At Risk

=cut

sub efo_id {
    deprecate(
        "Bio::EnsEMBL::Funcgen::Epigenome::efo_id() has been deprecated and will be removed in Ensembl release 89."
            . " Please use Bio::EnsEMBL::Funcgen::Epigenome::ontology_accession() instead."
    );
    return $_[0]->{ontology_accession};
}

=head2 efo_db_entry

  Example    : $epigenome->efo_accession
  Description: Returns the accession in the Experimental Factor Ontology (EFO) for this epigenome.
  Returntype : String
  Exceptions : 
  Caller     : 
  Status     : 

=cut

sub efo_accession {
  my $self = shift;
  my $efo_db_entry = $self->efo_db_entry;
  if (! defined $efo_db_entry) {
    return undef;
  }
  return $efo_db_entry->primary_id
}

=head2 efo_db_entry

  Example    : print $epigenome->efo_db_entry->primary_id;
  Description: Returns the DBEntry of the external reference to the Experimental Factor Ontology (EFO).
  Returntype : Bio::EnsEMBL::Funcgen::DBEntry
  Exceptions : 
  Caller     : 
  Status     : 

=cut

sub efo_db_entry {
  my $self = shift;
  return $self->_unique_db_entry('EFO');
}

=head2 epirr_accession

  Example    : print $epigenome->epirr_accession;
  Description: Returns the EpiRR accession for this epigenome.
  Returntype : String
  Exceptions : 
  Caller     : 
  Status     : 

=cut

sub epirr_accession {
  my $self = shift;
  my $epirr_db_entry = $self->epirr_db_entry;
  if (! defined $epirr_db_entry) {
    return undef;
  }
  return $epirr_db_entry->primary_id
}

=head2 epirr_db_entry

  Example    : print $epigenome->epirr_db_entry->primary_id;
  Description: Returns the DBEntry of the external reference to EpiRR.
  Returntype : Bio::EnsEMBL::Funcgen::DBEntry
  Exceptions : 
  Caller     : 
  Status     : 

=cut

sub epirr_db_entry {
  my $self = shift;
  return $self->_unique_db_entry('EpiRR')
}

sub _unique_db_entry {

  my $self = shift;
  my $external_db_name = shift;

  my $dbentry_adaptor = $self->adaptor->db->get_DBEntryAdaptor;

  my $efo_db_entry = $dbentry_adaptor->_fetch_by_object_type(
    $self->dbID,
    'epigenome',
    $external_db_name
  );
  
  if (
    (ref $efo_db_entry eq 'ARRAY')
    && (@$efo_db_entry == 1)
  ) {
    return $efo_db_entry->[0];
  }
  if (
    (ref $efo_db_entry eq 'ARRAY')
    && (@$efo_db_entry == 0)
  ) {
#     warn("No $external_db_name id defined for " . $self->name . "!\n");
    return undef;
  }
  use Data::Dumper;
  throw("Unexpected return value for $external_db_name id!\n" . Dumper($efo_db_entry));
}

=head2 ontology_accession

  Example    : my $ontology_accession = $epigenome->ontology_accession;
  Description: Getter of the ontology accession
  Returntype : String
  Exceptions : None
  Caller     : General
  Status     : At Risk

=cut

sub ontology_accession{  return $_[0]->{ontology_accession}; }


=head2 tissue

  Example    : my $tissue = $epigenome->tissue;
  Description: Getter of the tissue attribute for a given epigenome
  Returntype : String
  Exceptions : None
  Caller     : General
  Status     : At Risk

=cut

sub tissue{   return $_[0]->{tissue}; }

=head2 reset_relational_attributes

  Arg[1]     : Hashref containing linked objects - KEPT for compliance
               with similar methods in other modules
  Arg[2]     : Flag to avoid reseting of adaptor and dbID

  Description: Undefs Adaptor and dbID
               Useful when creating a cloned object for migration beween DBs
  Returntype : None
  Exceptions : None
  Caller     : Migration code
  Status     : At risk

=cut

sub reset_relational_attributes{
  my ($self, $params_hash, $no_db_reset) = shift;

  # undef the dbID and adaptor
  if(! $no_db_reset){
    $self->{adaptor} = undef;
    $self->{dbID}    = undef;
  }

  return;
}

=head2 compare_to

  Args[1]    : Bio::EnsEMBL::Funcgen::Storable (mandatory)
  Args[2]    : Boolean - Optional 'shallow' - no object methods compared
  Args[3]    : Arrayref - Optional list of Epigenome method names each
               returning a Scalar or an Array or Arrayref of Scalars.
               Defaults to: name table_name feature_class get_all_states
  Args[4]    : Arrayref - Optional list of Epigenome method names each
               returning a Storable or an Array or Arrayref of Storables.
               Defaults to: feature_type epigenome analysis get_support
  Example    : my %shallow_diffs = %{$rset->compare_to($other_rset, 1)};
  Description: Compare this Epigenome to another based on the defined scalar
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

  $scl_methods ||= [qw(name display_label description gender efo_id tissue)];

  return $self->SUPER::compare_to($obj, $shallow, $scl_methods,
                                  $obj_methods);
}


1;

