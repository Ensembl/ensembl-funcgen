#########
# Author:  Stefan Graf
# Maintainer: sg550@cam.ac.uk
# Created: 2009-05-11
# Last Modified: 2009-05-11
# Builds DAS features from eFG databases
#
package Bio::Das::ProServer::SourceAdaptor::Transport::efg;
use strict;
use warnings;
use Carp;
use Data::Dumper;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor;
use base qw(Bio::Das::ProServer::SourceAdaptor::Transport::generic);

our $VERSION  = do { my ($v) = (q$Revision: 1.1 $ =~ /\d+/mxg); $v; };

sub init {
	#warn "<<< transport::init >>>";
  my ($self) = @_;
  #warn($self->config->{'dbname'}) if $self->{'debug'};
  my ($basename, $species, $group, $data_version) 
	  = $self->config->{'dbname'} =~ m/([^_]+_)?([^_]+_[^_]+)_([^_]+)_([^_]+_[^_]+)$/;
  #print Dumper ($basename, $species, $group, $data_version);
  $self->{'_basename'}     = $basename;
  $self->{'_species'}      = $species;
  $self->{'_group'}        = $group;
  $self->{'_data_version'} = $data_version;
  $self->{'_db'}           = $self->adaptor();
#  $self->_apply_override;
#  $self->_load_registry;
  return;
}

sub adaptor {
  my ($self) = @_;
  my $cdb = Bio::EnsEMBL::DBSQL::DBAdaptor->new(
	  #-host => 'ensembldb.ensembl.org',
	  #-port => '5306',
	  #-user => 'anonymous',
	  -host => 'ens-staging',
	  -port => '3306',
	  -user => 'ensro',
	  -dbname => $self->{'_species'}.'_core_'.$self->{'_data_version'},
	  -species => $self->{'_species'}
	  );
  #warn("Got core database adaptor for ".$cdb->dbc->dbname."\n") if $self->{'debug'};

  my $db = Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor->new
	  (
	   -host    => $self->config->{'host'},
	   -port    => $self->config->{'port'},
	   -user    => $self->config->{'user'},
	   -pass    => $self->config->{'pass'} || undef,
	   -dbname  => $self->config->{'dbname'},
	   -species => $self->{'_species'},
	   -dnadb   => $cdb
     );
  #warn("Got funcgen database adaptor for ".$db->dbc->dbname."\n") if $self->{'debug'};

  return $db;
}

sub slice_adaptor {
  my ($self) = @_;
  return $self->{'_db'}->get_SliceAdaptor();
}

sub feature_set_adaptor {
  my ($self) = @_;
  return $self->{'_db'}->get_FeatureSetAdaptor();
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

sub gene_adaptor {
  my ($self, $species, $group) = @_;
  return $self->adaptor($species, $group)->get_GeneAdaptor();
}

sub chromosome_by_region { ## no critic
  my ($self, $chr, $start, $end, $species, $group) = @_;
  return $self->slice_adaptor($species, $group)->fetch_by_region('chromosome', $chr, $start, $end);
}

sub chromosomes {
  my ($self, $species, $group) = @_;
  return $self->slice_adaptor($species, $group)->fetch_all('chromosome');
}

sub gene_by_id {
  my ($self, $id, $species, $group) = @_;
  return $self->gene_adaptor($species, $group)->fetch_by_stable_id($id);
}

sub genes {
  my ($self, $species, $group) = @_;
  return $self->gene_adaptor($species, $group)->fetch_all();
}

sub version {
  my ($self) = @_;
  return Bio::EnsEMBL::Registry->software_version();
}

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

sub disconnect {
  my $self = shift;
  Bio::EnsEMBL::Registry->disconnect_all();
  $self->{'debug'} and carp "$self performed disconnect\n";
  return;
}

1;

=head1 NAME

Bio::Das::ProServer::SourceAdaptor::Transport::efg

=head1 AUTHOR

Stefan Graf <graef@ebi.ac.uk or sg550@cam.ac.uk>

=head1 LICENSE AND COPYRIGHT

Copyright (c) 2009 EMBL-EBI and CRUK-CRI

=cut
