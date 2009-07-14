=head1 NAME

efg_transport.pm

=head1 DESCRIPTION

Transport module to retrieve Annotated/ResultFeatures from an efg database

=head1 LICENCE

This code is distributed under an Apache style licence. Please see
http://www.ensembl.org/info/about/code_licence.html for details.

=head1 CONTACT

Please post comments/questions to the Ensembl development list
<ensembl-dev@ebi.ac.uk>

=cut

package Bio::Das::ProServer::SourceAdaptor::Transport::ensembl_funcgen;



#BEGIN {
#  my ($eroot) = $ENV{'ENS_ROOT'}     =~ m|([a-zA-Z0-9_/\.\-]+)|;
#  my ($broot) = $ENV{'BIOPERL_HOME'} =~ m|([a-zA-Z0-9_/\.\-]+)|;
#
#  unshift(@INC,"$eroot/ensembl/modules");
#  unshift(@INC,"$broot");
#}

use strict;
use Carp;
use Data::Dumper;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor;
#use Bio::EnsEMBL::Registry;
use Bio::Das::ProServer::SourceAdaptor::Transport::generic;

use vars qw(@ISA);
@ISA = qw(Bio::Das::ProServer::SourceAdaptor::Transport::generic);



#No init as data_version now either explicilty set as part of dnadb config
#Or autoselect in the DBAdaptor

sub adaptor {
    my $self = shift;

    unless($self->{'_adaptor'}) {
		  
        
	  #Can we add actaully ini file to this error?

	  $self->{'_adaptor'} ||= Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor->new
		(
		 -host    => $self->config->{'host'}     || 'localhost',
		 -user    => $self->config->{'username'} || 'ensro',
		 -dbname  => $self->config->{'dbname'}   || carp('No dbname defined for '.$self->{'dsn'}.' in your das_source.ini file'),
		 -species => $self->config->{'species'}  || carp('No species defined for '.$self->{'dsn'}.' in your das_source.ini file'),
		 -pass    => $self->config->{'password'} || '',
		 -port    => $self->config->{'port'}     || '3306',
		 -type    => 'funcgen',
		 -dnadb_assembly => $self->config('dnadb_assembly'),
		 -dnadb_host => $self->config('dnadb_host'),
		 -dnadb_port => $self->config('dnadb_post'),
		 -dnadb_user => $self->config('dnadb_user'),
		 -dnadb_pass => $self->config('dnadb_pass'),
		 


		);
	  print Dumper $self->{'_adaptor'} if ($self->{'debug'});

	  #Test the connection
	  $self->{'_adaptor'}->dbc->db_handle;
	  return $self->{'_adaptor'};
    }
}

sub chromosome_by_region {
  my ($self, $chr, $start, $end) = @_;
  my $slice = $self->adaptor->dnadb->get_SliceAdaptor->fetch_by_region
      (
       'chromosome', $chr, $start, $end
       );
  return $slice;
}

sub disconnect {
  my $self = shift;
  return unless (exists $self->{'_adaptor'});
  $self->{'_adaptor'}->disconnect();
  delete $self->{'_adaptor'};
  $self->{'debug'} and print STDERR "$self PERFORMED DBA DISCONNECT\n";
}

sub DESTROY {
  my $self = shift;
  $self->disconnect();
}


1;

__END__

#Are any of these used?
#Or is this hangover from ensembl adaptor?
#Using the registry will create problems
#as we may had multiple DBs with the same species
#which will cause the registry to redefine the species name

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
