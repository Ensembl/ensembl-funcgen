package Bio::EnsEMBL::Funcgen::RunnableDB::ChIPSeq::QcStoreProportionOfReadsInPeaks;

use warnings;
use strict;
use Data::Dumper;

use base 'Bio::EnsEMBL::Hive::Process';

sub run {
  my $self = shift;

  my $peak_file          = $self->param_required('peak_file');
  my $temp_dir           = $self->param_required('tempdir');
  my $peak_calling_name  = $self->param_required('peak_calling');
  my $bam_file           = $self->param_required('bam_file');
  my $species            = $self->param_required('species');
  
  my $num_reads_in_peaks           = $self->param_required('num_reads_in_peaks');
  my $num_reads_in_total           = $self->param_required('num_reads_in_total');
  my $proportion_of_reads_in_peaks = $self->param_required('proportion_of_reads_in_peaks');

  use Bio::EnsEMBL::Utils::Logger;
  my $logger = Bio::EnsEMBL::Utils::Logger->new();
  $logger->init_log;
  
  my $peak_calling_adaptor = Bio::EnsEMBL::Registry
    ->get_adaptor(
        $species, 
        'funcgen', 
        'PeakCalling'
    );
  my $peak_calling = $peak_calling_adaptor->fetch_by_name($peak_calling_name);

  my $logic_name = 'Proportion of reads in peaks';
  my @proportion_of_reads_in_peaks_analysis_details = (
    -logic_name      => $logic_name,
    -program         => 'proportion_of_reads_in_peaks.pl',
    -parameters      => undef,
    -description     => 'Computation of the proportion of reads in peaks',
    -display_label   => 'Proportion of reads in peaks',
    -displayable     => undef,
  );

  my $analysis_adaptor = Bio::EnsEMBL::Registry->get_adaptor( $species, 'Funcgen', 'Analysis' );
  my $analysis = $analysis_adaptor->fetch_by_logic_name($logic_name);

  if (! $analysis ) {
    $logger->info("No analysis with logic name $logic_name found. Creating one.\n");
    $analysis = Bio::EnsEMBL::Analysis->new(@proportion_of_reads_in_peaks_analysis_details);
    $analysis_adaptor->store($analysis);
  }

  my $analysis_id = $analysis->dbID;
  $logger->info("The analysis_id is: $analysis_id.\n");

  my $dbc = $analysis_adaptor->db->dbc;
  #print Dumper($dbc);  die;

  my $sql_processor = sub {
    my $sql = shift;
    $logger->info("Running: " . $sql . "\n");
    
    eval {
      $dbc->do($sql);
    };
    if ($@) {
      
      my $has_already_been_stored = $@ =~ /Duplicate entry/;
      
      if (! $has_already_been_stored) {
        $self->throw($@);
      }
      
    }
  };

  create_table({ 
    sql_processor => $sql_processor,
  });

  create_insert_sql({
    analysis_id         => $analysis_id,
    peak_calling_id     => $peak_calling->dbID,
    proportion_of_reads_in_peaks => $proportion_of_reads_in_peaks,
    num_reads_in_total         => $num_reads_in_total,
    sql_processor       => $sql_processor,
  });

  $logger->finish_log;
  return;
}

=head2 create_insert_sql
=cut
sub create_insert_sql {

  my $param = shift;
  
  my $analysis_id         = $param->{analysis_id};
  my $peak_calling_id     = $param->{peak_calling_id};
  my $proportion_of_reads_in_peaks = $param->{proportion_of_reads_in_peaks};
  my $num_reads_in_total         = $param->{num_reads_in_total};
  my $sql_processor       = $param->{sql_processor};  

  my $sql = "INSERT INTO peak_calling_qc_prop_reads_in_peaks ("
    . "analysis_id, "
    . "peak_calling_id, "
    . "prop_reads_in_peaks, "
    . "total_reads"
    . ")  VALUES ("
  . (
    join ', ', (
      $analysis_id,
      $peak_calling_id,
      $proportion_of_reads_in_peaks,
      $num_reads_in_total
      )
    )
  . ");";
  $sql_processor->($sql);
}

=head2 create_table
=cut
sub create_table {

  my $param = shift;
  my $sql_processor = $param->{sql_processor};

my $sql = <<SQL
CREATE TABLE if not exists peak_calling_qc_prop_reads_in_peaks (
  peak_calling_qc_prop_reads_in_peaks_id int(10) unsigned NOT NULL AUTO_INCREMENT,
  analysis_id                            int(10) unsigned,
  peak_calling_id                        int(10) unsigned unique,
  prop_reads_in_peaks                    double default NULL,
  total_reads                            int(10),
  PRIMARY KEY ( peak_calling_qc_prop_reads_in_peaks_id )
);
SQL
;
  $sql_processor->($sql);
}

1;
