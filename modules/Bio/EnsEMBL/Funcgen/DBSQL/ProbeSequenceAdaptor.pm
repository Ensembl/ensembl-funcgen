=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2018] EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

  Questions may also be sent to the Ensembl help desk at
  <http://www.ensembl.org/Help/Contact>.


=head1 NAME

  Bio::EnsEMBL::Funcgen::DBSQL::ProbeFeatureTranscriptMappingAdaptor

=head1 SYNOPSIS

=cut

package Bio::EnsEMBL::Funcgen::DBSQL::ProbeSequenceAdaptor;

use strict;
use warnings;
use Bio::EnsEMBL::Utils::Exception qw( throw warning );
use DBI qw(:sql_types);

use vars '@ISA';
@ISA    = qw(Bio::EnsEMBL::DBSQL::BaseAdaptor);

sub _tables {
  return ['probe_seq', 'ps'];
}

sub _columns {
  my $self = shift;
  
  return qw(
    ps.probe_seq_id
    ps.sequence
  );
}

sub _objs_from_sth {
  my ($self, $sth) = @_;

  my (
    $sth_fetched_dbID,
    $sth_fetched_sequence,
  );

  $sth->bind_columns (
    \$sth_fetched_dbID,
    \$sth_fetched_sequence,
  );
  
  use Bio::EnsEMBL::Funcgen::ProbeSequence;

  my @return_objects;
  ROW: while ( $sth->fetch() ) {
  
    my $probe_sequence = Bio::EnsEMBL::Funcgen::ProbeSequence->new(
      -dbID         => $sth_fetched_dbID,
      -sequence     => $sth_fetched_sequence,
      -adaptor      => $self->db,
    );
    push @return_objects, $probe_sequence
  }
  return \@return_objects;
}

sub fetch_by_sequence {
  my $self     = shift;
  my $sequence = shift;

  my $constraint = "ps.sequence_upper_sha1 = sha1(upper(?))";
  $self->bind_param_generic_fetch($sequence, SQL_VARCHAR);
  
  my $probe_sequence = $self->generic_fetch($constraint);
  
  if (!$probe_sequence || @$probe_sequence==0) {
    return;
  }
  if (@$probe_sequence!=1) {
    throw("Found ". @$probe_sequence ." entries in the probe sequence table with the same probe sequence!");
  }
  return $probe_sequence->[0];
}

sub store {
  my $self = shift;
  my $probe_sequence = shift;
  
  if (! defined $probe_sequence->sequence) {
    use Carp;
    confess(
      "Probe sequence not set:\n"
      . Dumper($probe_sequence)
    );
  }

  # Supress flurry of error messages that can appear when trying to store 
  # probes with the same sequence. The error is caught and handled below.
  #
  my $saved_printerror_value = $self->db->dbc->db_handle->{PrintError};
  $self->db->dbc->db_handle->{PrintError} = 0;
  
  my $probe_seq_sth = $self->prepare('insert into probe_seq (sequence, sequence_upper, sequence_upper_sha1) values (?, upper(?), cast(sha1(upper(?)) as char))');
  my $probe_seq_id;
  
  eval {
    $probe_seq_sth->bind_param(1, $probe_sequence->sequence, SQL_VARCHAR);
    $probe_seq_sth->bind_param(2, $probe_sequence->sequence, SQL_VARCHAR);
    $probe_seq_sth->bind_param(3, $probe_sequence->sequence, SQL_VARCHAR);

    $probe_seq_sth->execute;
    $probe_seq_id = $probe_seq_sth->{mysql_insertid};
  };
  
  if ($@) {
    my $error_msg = $@;

    # Check, if the exception was triggered by the index on the probe_sha1 
    # column.
    #
    if ($error_msg=~/DBD::mysql::st execute failed: Duplicate entry/) {
      my $probe_sequence_from_db = $self->fetch_by_sequence($probe_sequence->sequence);
      $probe_seq_id = $probe_sequence_from_db->dbID;
    } else {
      use Carp;
      confess($@);
    }
  }
  $probe_sequence->dbID($probe_seq_id);
  $probe_sequence->adaptor($self);
  
  $self->db->dbc->db_handle->{PrintError} = $saved_printerror_value;
  
  return $probe_sequence;
}

1;
