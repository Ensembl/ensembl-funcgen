#
# Ensembl module for Bio::EnsEMBL::Funcgen::InputSubset
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

Bio::EnsEMBL::InputSubset - A module to represent InputSubset object.
 

=head1 SYNOPSIS

use Bio::EnsEMBL::Funcgen::InputSubset;

my $input_subset = Bio::EnsEMBL::Funcgen::InputSubset->new
                    (
                     -DBID        => $dbID,
                     -ADAPTOR     => $self,
                     -NAME        => $name,
                     -INPUT_SET   => $iset,
                     -archive_id  => $archive_id,
                     -display_url => $display_url,
                     -replicate   => $iss_rep,
                     -is_control  => $is_control,
                    );



=head1 DESCRIPTION

An InputSubset object represents an individual distinct input within a given InputSet. This 
normally translates to single file or replicate. There is no dedicated InputSubsetAdaptor,
store and fetch functionality is embedded within the InputSetAdaptor.

=cut

use strict;
use warnings;

package Bio::EnsEMBL::Funcgen::InputSubset;

use Bio::EnsEMBL::Utils::Argument qw( rearrange );
use Bio::EnsEMBL::Utils::Exception qw( throw );
use Bio::EnsEMBL::Funcgen::Storable;

use vars qw(@ISA);
@ISA = qw(Bio::EnsEMBL::Funcgen::Storable);


=head2 new

  Example    : my $eset = Bio::EnsEMBL::Funcgen::InputSubset->new
                            (
                             -DBID        => $dbID,
                             -ADAPTOR     => $self,
                             -NAME        => $name,
                             -INPUT_SET   => $iset,
                             -archive_id  => $archive_id,
                             -display_url => $display_url,
                             -replicate   => $iss_rep,
                             -is_control  => $is_control,
                            );


  Description: Constructor for InputSubset objects.
  Returntype : Bio::EnsEMBL::Funcgen::InputSubset
  Exceptions : Throws if no name defined
  Caller     : InputSetAdaptor
  Status     : At risk

=cut

sub new {
  my $caller = shift;
	
  my $class = ref($caller) || $caller;
	
  my $self = $class->SUPER::new(@_);
	
  #do we need to add $fg_ids to this?  Currently maintaining one feature_group focus.(combi exps?)
  my ($name, $eset, $archive_id,
	 $display_url, $rep, $is_control)
    = rearrange(['NAME', 'INPUT_SET', 'ARCHIVE_ID', 
				 'DISPLAY_URL', 'REPLICATE', 'IS_CONTROL'], @_);
  
  
  throw('Must provide a name argument') if ! defined $name;
  

  #can't do is_stored_and_valid here as we don't adaptor

  if(! (ref($eset) && 
        $eset->isa('Bio::EnsEMBL::Funcgen::InputSet')
        && $eset->dbID())){
    throw('Must provide a valid stored Bio::EnsEMBL::Funcgen::InputSet argument');
  }
  

  $self->{name}        = $name;
  $self->{input_set}   = $eset;
  $self->{archive_id}  = $archive_id;
  $self->{display_url} = $display_url;
  $self->{replicate}   = $rep;
  $self->{is_control}  = $is_control;

  $eset->_add_new_subset($self);

  return $self;
}


=head2 name

  Example    : my $name = $exp_sset->name();
  Description: Getter for the name of this InputSubset.
  Returntype : String
  Exceptions : None
  Caller     : General
  Status     : At Risk

=cut

sub name { return $_[0]->{name}; }


=head2 input_set

  Example    : my $input_set = $input_sset->input_set;
  Description: Getter for the input_set attribute of this InputSubset.
  Returntype : Bio::EnsEMBL::Funcgen::InputSet
  Exceptions : None
  Caller     : General
  Status     : At Risk

=cut

sub input_set {    return $_[0]->{input_set}; }


=head2 archive_id

  Example    : my $archive_id = $inp_sset->archive_id;
  Description: Getter for the archive of this InputSubset.
  Returntype : String
  Exceptions : None
  Caller     : General
  Status     : At Risk

=cut

sub archive_id { return $_[0]->{archive_id}; }


=head2 display_url

  Example    : my $url = $inp_sset->displau_url;
  Description: Getter for the display_url of this InputSubset.
  Returntype : String
  Exceptions : None
  Caller     : General
  Status     : At Risk

=cut

sub display_url{ return $_[0]->{display_url}; }


=head2 replicate

  Example    : my $rep = $inp_sset->replicate;
  Description: Getter for the replicate attribute of this InputSubset.
  Returntype : Integer
  Exceptions : None
  Caller     : General
  Status     : At Risk

=cut

sub replicate { return $_[0]->{replicate}; }


=head2 is_control

  Example    : if($input_sset->is_control){ # Do some control specific stuff here }
  Description: Getter for the is_control attribute of this InputSubset.
  Returntype : Boolean
  Exceptions : None
  Caller     : General
  Status     : At Risk

=cut

sub is_control { return $_[0]->{is_control}; }

1;

