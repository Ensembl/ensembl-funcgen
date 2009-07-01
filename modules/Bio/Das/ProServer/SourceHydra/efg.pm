#########
# Author:  Stefan Graf
# Maintainer: sg550@cam.ac.uk
# Created: 2009-05-11
# Last Modified: 2009-05-11
# Builds DAS features from eFG databases
#
# DBI-driven sourceadaptor broker
#
package Bio::Das::ProServer::SourceHydra::efg;
use strict;
use warnings;
use base qw(Bio::Das::ProServer::SourceHydra);
use English qw(-no_match_vars);
use Data::Dumper;
use Carp;
use Readonly;

our $VERSION       = do { my ($v) = (q$Revision: 1.1 $ =~ /\d+/mxg); $v; };
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
      $self->{'debug'} and carp qq(Fetching sources using feature_set adaptor);
	  @{$self->{'_sources'}} = map {
		  $hydraname.'_'.$_->name;
		  #$_->name if $_->name =~ m/^HMM.*A488_00/
	  } @{$self->transport->feature_set_adaptor->fetch_all()};
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

=head1 NAME

Bio::Das::ProServer::SourceHydra::efg

=head1 AUTHOR

Stefan Graf <graef@ebi.ac.uk or sg550@cam.ac.uk>.

=head1 LICENSE AND COPYRIGHT

Copyright (c) 2009 EMBL-EBI

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself.  See DISCLAIMER.txt for
disclaimers of warranty.

=cut
