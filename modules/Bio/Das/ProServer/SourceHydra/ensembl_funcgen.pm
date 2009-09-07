=pod

=head1 NAME

Bio::Das::ProServer::SourceHydra::ensembl_funcgen

=head1 DESCRIPTION

Dynamic source adaptor broker, serves hydra source adaptors from 
Result/FeatureSets with status 'DAS_DISPLAYABLE'


=head1 LICENSE

  Copyright (c) 1999-2009 The European Bioinformatics Institute and
  Genome Research Limited.  All rights reserved.

  This software is distributed under a modified Apache license.
  For license details, please see

    http://www.ensembl.org/info/about/code_licence.html

=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <ensembl-dev@ebi.ac.uk>.

  Questions may also be sent to the Ensembl help desk at
  <helpdesk@ensembl.org>.


=cut


#To do
# 1 Incorporate result sets
# 2 Change timeout?
# 3 More debug stuff
# 4 Docs

package Bio::Das::ProServer::SourceHydra::ensembl_funcgen;
use strict;
use warnings;
use base qw(Bio::Das::ProServer::SourceHydra);
use English qw(-no_match_vars);
use Data::Dumper;
use Carp;
use Readonly;

our $VERSION       = do { my ($v) = (q$Revision: 1.2 $ =~ /\d+/mxg); $v; };
Readonly::Scalar our $CACHE_TIMEOUT => 30;

#########
# the purpose of this module:
#
sub sources {
	#warn "<<< hydra::sources >>>";
  my ($self)    = @_;
  my $hydraname = $self->{'dsn'};
  my $sql       = $self->config->{'query'};
  my $now       = time;

  #########
  # flush the table cache *at most* once every $CACHE_TIMEOUT
  # This may need signal triggering to have immediate support
  #
  if($now > ($self->{'_sourcecache_timestamp'} || 0)+$CACHE_TIMEOUT) {
    $self->{'debug'} and carp qq(Flushing table-cache for $hydraname);
    delete $self->{'_sources'};
    $self->{'_sourcecache_timestamp'} = $now;
  }
  
  # Use the configured query to find the names of the sources
  if(!exists $self->{'_sources'}) {
    $self->{'_sources'} = [];
    eval {
      $self->{'debug'} and carp('Fetching sources using '.ucfirst($self->config->{'set_type'}).'SetAdaptor');
	  my $adaptor_method = 'get_'.ucfirst($self->config->{'set_type'}).'SetAdaptor';
	  my $adaptor = $self->transport->adaptor->$adaptor_method;

	  #map is slow here?
	  @{$self->{'_sources'}} = map {
		$hydraname.'_'.$_->name;
	
	  } @{$adaptor->fetch_all('DAS_DISPLAYABLE')};
      $self->{'debug'} and carp qq(@{[scalar @{$self->{'_sources'}}]} sources found);
      1;

    } or do {
      carp "Error scanning database: $EVAL_ERROR";
      delete $self->{'_sources'};
    };
  }

  return @{$self->{'_sources'} || []};
}

1;

