package Bio::EnsEMBL::Funcgen::RunnableDB::ProbeMapping::JobFactory;

use strict;
use base ('Bio::EnsEMBL::Hive::Process');

sub run {
    my $self = shift;
    my $probe_directories = $self->param('probe_directories');

    use Bio::EnsEMBL::Utils::Logger;
    my $logger = Bio::EnsEMBL::Utils::Logger->new();
    
    use Data::Dumper;
    use File::Spec;

    opendir(D, $probe_directories) || die "Can't open directory $probe_directories: $!\n";
    my @list = grep !/^\.\.?$/, readdir(D);
    closedir(D);

    my $input_id = [];
    DIR: foreach my $subdirectory (@list) {
    
      my $full_path = File::Spec->catfile($probe_directories, $subdirectory);
      
      # Skip any files that may be in there, only directories are processed.
      #
      next DIR if (! -d $full_path);
      next DIR if ($full_path =~ /GENOMICSEQS/);
      next DIR if ($full_path =~ /TRANSCRIPTSEQS/);

      push @$input_id, {    
	array_format => $subdirectory,
	QUERYSEQS    => File::Spec->catfile($full_path, 'arrays.'.$subdirectory.'.fasta'),
      };
    }    
    foreach my $current_input_id (@$input_id) {
      $self->dataflow_output_id( $current_input_id, 1 );
    }
}

1;
