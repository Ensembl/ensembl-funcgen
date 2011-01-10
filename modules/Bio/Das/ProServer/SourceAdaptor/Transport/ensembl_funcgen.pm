=head1 LICENSE

  Copyright (c) 1999-2011 The European Bioinformatics Institute and
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

ensembl_funcgen.pm

=head1 DESCRIPTION

Transport module to retrieve Annotated/ResultFeatures from an efg database

=cut

#To do
# 1 More debug stuff
# 2 Docs


package Bio::Das::ProServer::SourceAdaptor::Transport::ensembl_funcgen;

use strict;
use Carp;
use Data::Dumper;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Registry;
use Bio::Das::ProServer::SourceAdaptor::Transport::generic;

use vars qw(@ISA);
@ISA = qw(Bio::Das::ProServer::SourceAdaptor::Transport::generic);



#No init as data_version now either explicilty set as part of dnadb config
#Or autoselect in the DBAdaptor
#Add some debug stuff?

sub adaptor {
  my $self = shift;


    if(! defined $self->{'_adaptor'}) {
		  
	  #Can get problems here with Registry persisting after transport has gone out of scope 
	  #and been destroyed.
	  #This happens for every source after the 1st DB, as the DB comparison in ConfigRegistry is
	  #always against the first stored

	  #Ensembl API was never designed to run like this!

	  #Let's test the registry version to see if it matches what we want and reinstate it?
	  #Otherwise we are going to have to deal with incremented species names in the DBAdaptor
	  #No, this should be done in ConfigRegistry

	  #do we need to parse coordinates here to set dnadb_assembly?
	  #We should really be able to run this from the same dnadb so long as it supports it
	  #This mean we need to set a version to pass to SliceAdaptor->fetch_by_region
	  #We also need to test wether it is available in the dnadb that is defined or detected

	  #Should we make this optional?  Just rely on default dnadb chosen by efg?
	  my $dnadb = Bio::EnsEMBL::DBSQL::DBAdaptor->new
		( 
		 -host   => $self->config->{'dnadb_host'},
		 -dbname => $self->config->{'dnadb_name'}, 
		 -port   => $self->config->{'dnadb_port'},
		 -user   => $self->config->{'dnadb_user'},
		 -pass   => $self->config->{'dnadb_pass'},
		 -group  => 'core',
		);

	  #Now check we have the required CS version in the dnadb?
	  #This is already done in general_DAS_config.pl
	  
	  $self->{'_adaptor'} ||= Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor->new
		(
		 -host    => $self->config->{'host'}     || 'localhost',
		 -user    => $self->config->{'username'} || 'ensro',
		 -dbname  => $self->config->{'dbname'}   || carp('No dbname defined for '.$self->{'dsn'}.' in your das_source.ini file'),
		 -species => $self->config->{'species'}  || carp('No species defined for '.$self->{'dsn'}.' in your das_source.ini file'),
		 -pass    => $self->config->{'password'},
		 -port    => $self->config->{'port'},
		 -dnadb   => $dnadb,
		 #-group   => 'funcgen',#Now set as default
		 #The following are never set by GenerateConfig
		 #So ensembldb would be used, which is fine for just fetching slice info
		 #Will need to set this if we don't have it on ensembldb!

		 #-dnadb_assembly => $self->config->{'dnadb_assembly'},
		 #-dnadb_host => $self->config->{'dnadb_host'},
		 #-dnadb_port => $self->config->{'dnadb_post'},
		 #-dnadb_user => $self->config->{'dnadb_user'},
		 #-dnadb_pass => $self->config->{'dnadb_pass'},
		 #Need to add dnadbname here!
		);

	  #print Dumper $self->{'_adaptor'} if ($self->{'debug'});

	  #Test the connection
	  $self->{'_adaptor'}->dbc->db_handle;

	  #Set disconnect when inactive to 1 to enable autoreconnect
	  $self->{'_adaptor'}->dbc->disconnect_when_inactive(1);
	  $self->{'_adaptor'}->dnadb->dbc->disconnect_when_inactive(1);
    }
  
  return $self->{'_adaptor'};	
}


#To support dbi hydra
sub dbh{
  my $self = shift;
  return $self->adaptor->dbc->db_handle;
}

#Do we even really need this?
#Ensembl API replaces what should be in here

sub chromosome_by_region {
  my ($self, $chr, $start, $end, $strand) = @_;

  #$cs_version should be set once for this transport based on the coordinates config or the set name?

  #This should really be region by name, and use the toplevel?
  #However, we can't use top level for old assemblies
  
  return $self->adaptor->dnadb->get_SliceAdaptor->fetch_by_region('chromosome', $chr, $start, $end, $strand, $self->coord_system_version);
}


sub get_coord_system_versions{
  my ($self) = @_;

  if(! defined $self->{'_coord_versions'}){
	$self->{'_coord_versions'} = [];
	my $coords = $self->config->{'coordinates'};
	$coords =~ s/;$//;
	my @cs_versions = split /;\s*/, $coords;
	map { $_ =~ s/\s*,.*//} @cs_versions;
	map { $_ =~ s/_//} @cs_versions;
	@{$self->{'_coord_versions'}} = @cs_versions;
  }

  return $self->{'_coord_versions'};
}

sub disconnect {
  my $self = shift;
  return unless (exists $self->{'_adaptor'});

  $self->{'_adaptor'}->dbc->db_handle->disconnect();
  delete $self->{'_adaptor'};
  $self->{'debug'} and print STDERR "$self PERFORMED DBA DISCONNECT\n";
}



sub fetch_set{
  my ($self, $hydra_set_name) = @_;


  #Don't know about methods hydra or dsn here, this is done in SourceAdaptor init
  #If we have a set_name, we know this is a hydra set

  my ($errmsg, $set);

  if($hydra_set_name){
	$errmsg = ' name '.$hydra_set_name;
  }
  else{
	$errmsg = ' dbID '.$self->config->{'set_id'};
  }

  #Fetch the set
  if($self->config->{'set_type'} eq 'feature'){

	if($hydra_set_name){
		$set = $self->adaptor->get_FeatureSetAdaptor->fetch_by_name($hydra_set_name);
	  }
	else{
	  #we have a set_id in the ini
	  $set = $self->adaptor->get_FeatureSetAdaptor->fetch_by_dbID($self->config->{'set_id'});
	}
  }
  elsif($self->config->{'set_type'} eq 'result'){

	if($hydra_set_name){
	  #Do we need to add cell type, ftype and analysis to config to enable unique set fetch?
	  #Would also have to add these to the dsn and title to avoid overwriting
	  #Just assume we are dealing with uniquely named sets for now
	  my @sets = @{$self->adaptor->get_ResultSetAdaptor->fetch_all_by_name($hydra_set_name)};
	  
	  if(scalar(@sets) >1){
		  die "Found non-unique ResultSet name ".$self->{'set_name'};
		}
	  
	  $set = $sets[0];
	}
	else{
	  #we have a set_id in the ini
	  $set = $self->adaptor->get_ResultSetAdaptor->fetch_by_dbID($self->config->{'set_id'});
	}
  }
  else{
	#Can we die here instead or will that prevent other sources from being loaded?
	#This should be eval'd in the caller
	die 'Bio::Das::ProServer::SourceAdaptor::ensembl_set does not support set_type '.$self->config->{'set_type'};
  }

  #We need to test feature_set here
  if (! defined  $set){
	die 'Cannot find '.ucfirst($self->config->{'set_type'})."Set with $errmsg";
  }
  
  return $set;
}

#Do we need to rewrite this to use the LAST_DAS_UPDATE status or meta entry?
#Rather than checking the whole DB? This should be a lot faster and may avoid timeouts

sub last_modified {
  my ($self) = @_;
  my $dbc = $self->adaptor()->dbc();
  my $sth = $dbc->prepare(q(SHOW TABLE STATUS));
  $sth->execute();
  my $server_text = [sort { $b cmp $a } ## no critic
                     keys %{ $sth->fetchall_hashref('Update_time') }
                    ]->[0]; # server local time
  $sth->finish();
  $sth = $dbc->prepare(q(SELECT UNIX_TIMESTAMP(?) as 'unix'));
  $sth->execute($server_text); # sec since epoch
  my $server_unix = $sth->fetchrow_arrayref()->[0];
  $sth->finish();
  return $server_unix;
}

sub version {
  my ($self) = @_;
  return Bio::EnsEMBL::Registry->software_version();
}



1;

__END__

#DO NOT NEED THIS
#And in fact seems to cause problems with mysterious disconnections on and autreconnect handle
#sub DESTROY {
#  my $self = shift;

#  $self->disconnect();
#}




#These are from the standard ensembl transport
#We could resurrect these(or inherit from standard ensembl transport)
#Currently not being used


sub chromosomes {
  my ($self, $species, $group) = @_;
  return $self->slice_adaptor($species, $group)->fetch_all('chromosome');
}

sub _load_registry {
  my ($self) = @_;
  if (!$self->config->{'skip_registry'}) {
    Bio::EnsEMBL::Registry->load_registry_from_db(
      -host => 'ensembldb.ensembl.org',
      -user => 'anonymous',
      -verbose => $self->{'debug'} );
  }
  Bio::EnsEMBL::Registry->set_disconnect_when_inactive();
  return;
}

sub _apply_override {
  my ($self) = @_;
  my $dbname = $self->config->{'dbname'};

  if ($dbname) {

    my ($species, $group) = $dbname =~ m/([a-z_]+)_([a-z]+)_\d+/mx;
    if ($species  eq 'ensembl') {
      $species = 'multi';
    }
    $self->{'_species'} ||= $species;
    $self->{'_group'}   ||= $group;

    if (!$self->{'_species'} || !$self->{'_group'}) {
      croak "Unable to parse database species and group: $dbname";
    }

    $self->{'debug'} && carp sprintf "Overriding database with %s (%s,%s)\n",
                                     $dbname, $self->{'_species'}, $self->{'_group'};

    # This is a map from group names to Ensembl DB adaptors.
    # Taken from Bio::EnsEMBL::Registry
    my %group2adaptor = (
      'blast'   => 'Bio::EnsEMBL::External::BlastAdaptor',
      'compara' => 'Bio::EnsEMBL::Compara::DBSQL::DBAdaptor',
      'core'    => 'Bio::EnsEMBL::DBSQL::DBAdaptor',
      'estgene' => 'Bio::EnsEMBL::DBSQL::DBAdaptor',
      'funcgen' => 'Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor',
      'haplotype' =>
        'Bio::EnsEMBL::ExternalData::Haplotype::DBAdaptor',
      'hive' => 'Bio::EnsEMBL::Hive::DBSQL::DBAdaptor',
      'lite' => 'Bio::EnsEMBL::Lite::DBAdaptor',
      'otherfeatures' => 'Bio::EnsEMBL::DBSQL::DBAdaptor',
      'pipeline' =>
        'Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor',
      'snp' =>
        'Bio::EnsEMBL::ExternalData::SNPSQL::DBAdaptor',
      'variation' =>
        'Bio::EnsEMBL::Variation::DBSQL::DBAdaptor',
      'vega' => 'Bio::EnsEMBL::DBSQL::DBAdaptor',
    );

    my $adaptorclass = $group2adaptor{ $self->{'_group'} } || croak 'Unknown database group: '.$self->{'_group'};
    # Creating a new connection will add it to the registry.
    eval "require $adaptorclass";
    $@ && die $@;
    $adaptorclass->new(
      -host    => $self->config->{'host'}     || 'localhost',
      -port    => $self->config->{'port'}     || '3306',
      -user    => $self->config->{'username'} || 'ensro',
      -pass    => $self->config->{'password'},
      -dbname  => $dbname,
      -species => $self->{'_species'},
      -group   => $self->{'_group'},
    );
  }
  return;
}

1;
