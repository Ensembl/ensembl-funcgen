package Bio::EnsEMBL::Funcgen::RunnableDB::ProbeMapping::DumpUnmappedSeqs;

use strict;
use base ('Bio::EnsEMBL::Hive::Process');
use Hash::Util qw( lock_keys );

sub run {
    my $self = shift;
    
    my $unmapped_sequences_file = $self->param('unmapped_sequences_file');
    my $species                 = $self->param('species');
    
    my $funcgen_adaptor         = Bio::EnsEMBL::Registry->get_DBAdaptor($species, 'funcgen');
    my $dbc_tracking = $funcgen_adaptor->dbc;

    my $helper = Bio::EnsEMBL::Utils::SqlHelper->new(
      -DB_CONNECTION => $dbc_tracking 
    );

    my $unmapped_sequences_count =
      $helper->execute_single_result(
	-SQL => 'select count(*) from probe_seq',
    );

    use Bio::EnsEMBL::Utils::Logger;
    my $logger = Bio::EnsEMBL::Utils::Logger->new();
    $logger->info("There are $unmapped_sequences_count unmapped probe sequences in the database.\n");

    if ($unmapped_sequences_count==0) {
      die("There are no unmapped probes in the database.");
    }

    my $sth = $dbc_tracking->prepare('select probe_seq_id, sequence from probe_seq');
    $sth->execute;

    use Bio::Seq;
    use Bio::SeqIO;

    my $out = Bio::SeqIO->new(
      -file   => '>' . $unmapped_sequences_file,
      -format => 'Fasta'
    );

    my $progressbar_id = $logger->init_progress($unmapped_sequences_count, 100);
    $logger->info("Writing unmapped sequences to $unmapped_sequences_file\n");

    my $skipped_sequences = 0;
    
    my $num_sequence_written=0;
    PROBE_SEQUENCE:
    while (my $data = $sth->fetchrow_hashref) {
    
      lock_keys(%$data);
      
      if (!$data->{sequence}) {
        use Data::Dumper;
        warn("This probe had no sequence!\n" . Dumper($data));
        $skipped_sequences++;
        next PROBE_SEQUENCE;
      }

      my $seq_obj = Bio::Seq->new(
	-id       => $data->{probe_seq_id},
	-seq      => $data->{sequence},
	-alphabet => 'dna',
      );
      $out->write_seq($seq_obj);
      
      $num_sequence_written++;
      $logger->log_progressbar($progressbar_id, $num_sequence_written);
    }
    $logger->info("Done writing $num_sequence_written sequences.");
    
    if ($num_sequence_written + $skipped_sequences != $unmapped_sequences_count) {
      my $msg = "$unmapped_sequences_count sequences had to be written, but only wrote $num_sequence_written sequences and skipped $skipped_sequences sequences!";
      $logger->error($msg);
      die($msg)
    }
}

1;
