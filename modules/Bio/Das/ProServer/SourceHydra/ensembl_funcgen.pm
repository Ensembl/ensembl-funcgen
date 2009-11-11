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
# 5 Integrate ensembl_funcgen_reads, so we can set bed_ tables as DAS_DISPLAYABLE


package Bio::Das::ProServer::SourceHydra::ensembl_funcgen;
use strict;
use warnings;
use base qw(Bio::Das::ProServer::SourceHydra);
use English qw(-no_match_vars);
use Data::Dumper;
use Carp;
use Readonly;

our $VERSION       = do { my ($v) = (q$Revision: 1.4 $ =~ /\d+/mxg); $v; };
Readonly::Scalar our $CACHE_TIMEOUT => 30;


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

	  #We need to get coord system here from source name and use this to bring back sets 
	  #with IMPORTED_ASSEMBLY information
	  #Need to do this for every coords version!

	  #map is slow here?
	#  @{$self->{'_sources'}} = map {
		#We need to add CS version here to differentiate between the same set on different assemblies
		#Let's do this only to those which are available on two versions?
		#Or should we add this by default?
	  #How will transport know which assembly to use?
		#Can we see in Ensembl that a source is on a different assembly?
	#	$hydraname.'_'.$_->name;
	
	  #} @{$adaptor->fetch_all('DAS_DISPLAYABLE')};


	  #This should now only ever return one which should match the hydra source name
	  #Can't have multiple coordinates entries unless all sources are present on both?
	  #Would also need to support assembly in features query! 
	  my @cs_versions = @{$self->transport->get_coord_system_versions};

	  foreach my $set(@{$adaptor->fetch_all('DAS_DISPLAYABLE')}){

		foreach my $cs_version(@cs_versions){

		  if($set->has_status("IMPORTED_${cs_version}")){
			#Have to add the version to the end of the source name?
			push @{$self->{'_sources'}}, $hydraname.'_'.$set->name.':'.$cs_version;
		  }
		}
	  }

      $self->{'debug'} and carp qq(@{[scalar @{$self->{'_sources'}}]} sources found);
	  1;#To make the eval return true

    } or do {
      carp "Error scanning database: $EVAL_ERROR";
      delete $self->{'_sources'};
    };
  }

  return @{$self->{'_sources'} || []};
}

1;

