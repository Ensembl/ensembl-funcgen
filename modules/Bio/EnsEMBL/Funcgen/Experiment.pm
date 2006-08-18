
=head1 NAME

Bio::EnsEMBL::Funcgen::Experiment
  
=head1 SYNOPSIS

use Bio::EnsEMBL::Funcgen::OligoArray;

my $array = Bio::EnsEMBL::Funcgen::Experiment->new(
												   -ADAPTOR             => $self,
												   -NAME                => $name,
												   -GROUP               => $group,
												   -DATE                => $date,
												   -PRIMARY_DESIGN_TYPE => $p_design_type,
												   -DESCRIPTION         => $description,
                                                   );

my $db_adaptor = Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor->new(...);
my $exp_adaptor = $db_adaptor->get_ExperimentAdaptor();
my $exp = $exp_adaptor->fetch_by_name($exp_name)

=head1 DESCRIPTION

An Experiment object represents an experiment instance . The data
are stored in the experiment, egroup, target, design_type and 
experimental_variable tables.


=head1 AUTHOR

This module was created by Nathan Johnson.

This module is part of the Ensembl project: http://www.ensembl.org/

=head1 CONTACT

Post comments or questions to the Ensembl development list: ensembl-dev@ebi.ac.uk

=head1 METHODS

=cut


################################################################################

package Bio::EnsEMBL::Funcgen::Experiment;

use warnings;
use strict;

use Bio::EnsEMBL::Utils::Argument qw( rearrange );
use Bio::EnsEMBL::Utils::Exception qw( throw warning );
use Bio::EnsEMBL::Storable;


use vars qw(@ISA);

@ISA = qw(Bio::EnsEMBL::Storable);



=head2 new

  Arg [-NAME]: string - the name of this experiment
  Arg [-]: string - the type of thi

  Example    : my $array = Bio::EnsEMBL::Funcgen::Experiment->new(
												                  -NAME                => $name,
												                  -GROUP               => $group,
												                  -DATE                => $date,
												                  -PRIMARY_DESIGN_TYPE => $p_design_type,
												                  -DESCRIPTION         => $description,
                                                 				 );

#experimental_variables?
#design_type
#target(s)?

  Description: Creates a new Bio::EnsEMBL::Funcgen::Experiment object.
  Returntype : Bio::EnsEMBL::Funcgen::Experiment
  Exceptions : None
  Caller     : General
  Status     : Medium Risk

=cut

sub new {
	my $caller = shift;

	my $class = ref($caller) || $caller;

	my $self = $class->SUPER::new(@_);

	my ($name, $group_id, $group, $date, $p_dtype, $desc)
		= rearrange( ['NAME', 'GROUP_ID', 'GROUP', 'DATE', 'PRIMARY_DESIGN_TYPE', 'DESCRIPTION'], @_ );
	
	$self->name($name)          if defined $name;
	$self->group_id($group_id)  if defined $group_id;
	$self->group($group)        if defined $group;
	$self->date($date)          if defined $date;
	$self->primary_design_type($p_dtype)    if defined $p_dtype;
	$self->description($desc)   if defined $desc;

	return $self;
}


### GENERIC ACCESSOR METHODS ###

sub name{
	my ($self) = shift;	

	$self->{'name'} = shift if(@_);

	return $self->{'name'};
}

sub group_id{
	my ($self) = shift;	

	$self->{'grup_id'} = shift if(@_);

	return $self->{'group_id'};
}



sub group{
	my ($self) = shift;	

	if(@_){
		$self->{'group'} = shift;
	}elsif(! exists $self->{'group'}){
		$self->get_group_details();	
	}
		

	return $self->{'group'};
}

sub location{
	my ($self) = shift;	

	if(@_){
		$self->{'location'} = shift;
	}elsif(! exists $self->{'location'}){
		$self->get_group_details();		
	}
  
	return $self->{'location'};
}

sub get_group_details{
  my $self = shift;
  
  #do this by dbid
  throw("Not yet implemented");


}


sub contact{
	my ($self) = shift;	

	if(@_){
		$self->{'contact'} = shift;
	}elsif(! exists $self->{'contact'}){
		$self->get_group_details();
	}

	return $self->{'group'};
}


sub date{
	my $self = shift;

	$self->{'date'} = shift if(@_);

	return $self->{'date'};
}

sub description{
	my $self = shift;

	$self->{'description'} = shift if(@_);

	return $self->{'description'};
}


sub primary_design_type{
	my ($self) = shift;
	
	$self->{'primary_design_type'} = shift if(@_);

	return $self->{'primary_design_type'};
}


#methods?
#lazy load design_types and exp_variables
#target?  Is this a one to one?



1;

