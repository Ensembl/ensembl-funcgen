=pod 
=head1 NAME

Bio::EnsEMBL::Funcgen::RunnableDB::Segmentation::GenerateSegmentationParameterFile

=head1 DESCRIPTION
=cut

package Bio::EnsEMBL::Funcgen::RunnableDB::Segmentation::GenerateSegmentationParameterFile;

use warnings;
use strict;

use base 'Bio::EnsEMBL::Hive::Process';
use Hash::Util qw( lock_keys );

sub run {
  my $self = shift;
  
  my $file               = $self->param_required('file');
  my $segmentation_lists = $self->param_required('segmentation_lists');
  
  use Data::Dumper;
  print Dumper($segmentation_lists);

  print "Writing output to $file\n";
  
  open my $fh, '>', $file || die("Can't open $file for writing!");

  my @superclasses = keys %$segmentation_lists;
  
  foreach my $superclass (@superclasses) {
  
    my $class_to_segmentation_meta_data = $segmentation_lists->{$superclass};
    my @classes = keys %$class_to_segmentation_meta_data;
    
    foreach my $class (@classes) {
    
        my $segmentation_meta_data = $class_to_segmentation_meta_data->{$class};
        
        lock_keys(%$segmentation_meta_data);

        my $segmentation_name     = $segmentation_meta_data->{segmentation_name};
        my $learn_model_directory = $segmentation_meta_data->{learn_model_directory};
        
        my $line = 
            join 
                "\t", 
                'segmentation',
                $segmentation_name,
                'ChromHMM',
                $learn_model_directory
                ;
        $fh->print($line . "\n");
        print($line . "\n");
    }
  }
  $fh->close;
  return;
}

1;
