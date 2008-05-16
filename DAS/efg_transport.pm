=head1 NAME

efg_transport.pm

=head1 DESCRIPTION

Transport module to retrieve REsult features from efg database

=head1 LICENCE

This code is distributed under an Apache style licence. Please see
http://www.ensembl.org/info/about/code_licence.html for details.

=head1 AUTHOR

Stefan Graf, Ensembl Functional Genomics

=head1 CONTACT

Please post comments/questions to the Ensembl development list
<ensembl-dev@ebi.ac.uk>

=cut

package Bio::Das::ProServer::SourceAdaptor::Transport::efg_transport;



#BEGIN {
#  my ($eroot) = $ENV{'ENS_ROOT'}     =~ m|([a-zA-Z0-9_/\.\-]+)|;
#  my ($broot) = $ENV{'BIOPERL_HOME'} =~ m|([a-zA-Z0-9_/\.\-]+)|;
#
#  unshift(@INC,"$eroot/ensembl/modules");
#  unshift(@INC,"$broot");
#}

use strict;
use Data::Dumper;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Registry;
use Bio::Das::ProServer::SourceAdaptor::Transport::generic;

use vars qw(@ISA);
@ISA = qw(Bio::Das::ProServer::SourceAdaptor::Transport::generic);

sub adaptor {
    my $self = shift;
    unless($self->{'_adaptor'}) {
        
        my $host     = $self->config->{'host'}     || 'localhost';
        my $port     = $self->config->{'port'}     || '3306';
        my $dbname   = $self->config->{'dbname'};
        my $username = $self->config->{'username'} || 'ensro';
        my $password = $self->config->{'password'} || '';
        my $species  = $self->config->{'species'}  || 'homo_sapiens';
        my $data_version  = $self->config->{'data_version'};
        
        $self->{'_adaptor'} ||= Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor->new
            (
             -host    => $host,
             -user    => $username,
             -dbname  => $dbname,
             -species => $species,
             -pass    => $password,
             -port    => $port,
             );
        print Dumper $self->{'_adaptor'} if ($self->{'debug'});
        return $self->{'_adaptor'};
    }
}

sub chromosome_by_region {
  my ($self, $chr, $start, $end) = @_;
  my $slice = $self->adaptor->get_SliceAdaptor->fetch_by_region
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
