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

our $VERSION       = do { my ($v) = (q$Revision: 1.5 $ =~ /\d+/mxg); $v; };
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
  

  #Add support for basename eq bed and nothing else



  # Use the configured query to find the names of the sources
  if(!exists $self->{'_sources'}) {
    $self->{'_sources'} = [];
    eval {

	  my $basename = $self->config->{'basename'};	
	  ####CAN'T ACTUALLY HAVE MULTIPLE COORDS FOR THE SAME HYDRA
	  ####AS THE COORD SYSTEM IS AT HYDRA/TRANSPORT LEVEL
	  ####NOT THE SOURCE LEVEL, SO THIS GET'S RESET FOR SUBSEQUENT
	  ####SOURCES
	  my @cs_versions = @{$self->transport->get_coord_system_versions};

	  if(scalar(@cs_versions) > 1){
		carp("Hydra sources can currently only handle one coordinate system, please define separate config for each assembly:\t".join("\t", @cs_versions))
	  }

	  my $cs_version = $cs_versions[0];

	  if ($basename){

		if($basename eq 'bed'){
		  my $l = length $basename;
		  $self->{'debug'} and carp qq(Fetching tables like $basename%);
		  my $dbh =  $self->transport->adaptor->dbc->db_handle;
		  my $sth = $dbh->prepare(qq(SHOW TABLES LIKE "$basename%"));
		  $sth->execute();
	
		  
		  foreach my $table_ref(@{$sth->fetchall_arrayref()}){
			my ($table) = @$table_ref;
			#We now use the status table for custom bed tables
			#Which ties their use to the efg schema
			#Check for DAS_DISPLAYABLE here?
			#This is going to change when we import as result_features
			
			#Now check for IMPORTED_CSVERSION status
			#Do we need to bind here to protext against injection from $basename and coordinates?
			#Use default table_id of 1
			my $sql = "select sn.name from status s, status_name sn where s.table_name='$table' and s.table_id=1 and s.status_name_id=sn.status_name_id and sn.name='IMPORTED_${cs_version}'";
			my ($status_name) = $dbh->selectrow_array($sql);
			push @{$self->{'_sources'}}, $hydraname.($table =~ /^.{$l}(.*)$/mx)[0] if $status_name;
		  }
		}
		else{
		  #throw here as we don't support any other basenames
		  carp "The ensembl_funcgen hydra only supports a basename of bed, omit basename from config for set support";
		}
	  }
	  else{
		
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

	
	
		foreach my $set(@{$adaptor->fetch_all('DAS_DISPLAYABLE')}){
		 
		  foreach my $cs_version(@cs_versions){
			
			if($set->has_status("IMPORTED_${cs_version}")){
			  #Have to add the version to the end of the source name?
			  push @{$self->{'_sources'}}, $hydraname.'_'.$set->name.':'.$cs_version;
			}
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

