package Bio::EnsEMBL::Funcgen::Hive::Ftp::JobFactoryAnnotatedFeatures;

use strict;
use Data::Dumper;
use base ('Bio::EnsEMBL::Hive::Process');

sub run {
  my $self = shift;
#   my $db_conn  = $self->param('regulation_database');
  my $temp_dir = $self->param('temp_dir');
  my $species  = $self->param('species');


  my $annotated_feature_adaptor = Bio::EnsEMBL::Registry->get_adaptor( $species, 'Funcgen', 'AnnotatedFeature' );

#   use Data::Dumper;
#   print Dumper($self->data_dbc);
  use File::Path qw(make_path remove_tree);
  make_path($temp_dir);

  my $batch_size = 10000;
  my $sql = 'select floor(annotated_feature_id / ?) as batch_number, min(annotated_feature_id) as min_id, max(annotated_feature_id) as max_id from annotated_feature group by batch_number order by batch_number asc';
  
#   my $dbc = $self->data_dbc;
  my $funcgen_adaptor = Bio::EnsEMBL::Registry->get_DBAdaptor( $species, 'Funcgen' );
  my $dbc = $funcgen_adaptor->dbc;
  
  my $sth = $dbc->prepare($sql);
  $sth->bind_param(1, $batch_size);
  $sth->execute();
  
  my $gff_files_from_batches = $temp_dir . '/gff_files_from_batches.txt';
  
  open OUT, ">$gff_files_from_batches";
  
  while (my $current_batch = $sth->fetchrow_hashref) {
  
    my $batch_number = $current_batch->{batch_number};
    
    my @directory_components = $temp_dir;
    push @directory_components, reverse split '', $batch_number;
    
    my $min_id = $current_batch->{min_id};
    my $max_id = $current_batch->{max_id};
          
    my $directory = join '/', @directory_components;
    my $file      = join '_', (
      'annotated_features',
      $min_id,
      $max_id,
    );
    $file .= '.gff';
    
    $current_batch->{directory} = $directory;
    $current_batch->{file}      = $file;
    
    print OUT "$directory/$file\n";
    
    $self->dataflow_output_id($current_batch, 2);
  }
  close(OUT);
  
  $self->dataflow_output_id({
    gff_files_from_batches => $gff_files_from_batches,
    merged_gff             => $temp_dir . '/AnnotatedFeatures.gff',
  }, 1);
}

1;
