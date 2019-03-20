package Bio::EnsEMBL::Funcgen::RunnableDB::Ftp::CreatePeakCallingDescription;

use strict;
use Data::Dumper;
use Bio::EnsEMBL::Registry;
use base ('Bio::EnsEMBL::Hive::Process');

use constant {
  BRANCH_JOB_OUTPUT => 2,
};

sub run {
  my $self = shift;
  my $species   = $self->param('species');
  my $file_name = $self->param('file_name');
  my $registry  = $self->param('reg_conf');
  
  my $epigenome_id    = $self->param('epigenome_id');
  my $feature_type_id = $self->param('feature_type_id');
  
  my $epigenome_adaptor    = Bio::EnsEMBL::Registry->get_adaptor($species, 'funcgen', 'epigenome');
  my $feature_type_adaptor = Bio::EnsEMBL::Registry->get_adaptor($species, 'funcgen', 'featuretype');
  my $peak_calling_adaptor = Bio::EnsEMBL::Registry->get_adaptor($species, 'funcgen', 'peakcalling');
  
  my $epigenome    = $epigenome_adaptor->fetch_by_dbID($epigenome_id);
  my $feature_type = $feature_type_adaptor->fetch_by_dbID($feature_type_id);
  
  my $peak_callings = $peak_calling_adaptor
    ->fetch_all_by_Epigenome_FeatureType(
        $epigenome, 
        $feature_type
    );

  $self->say_with_header($file_name, 1);
  
  use File::Path qw( make_path );
  use File::Basename;
  my $directory = dirname($file_name);
  make_path($directory);
  
  $self->say_with_header("Creating description $file_name", 1);
  
  use Bio::EnsEMBL::Funcgen::Report::PeakCalling;
  my $report = Bio::EnsEMBL::Funcgen::Report::PeakCalling->new(
    -species      => $species,
    -registry     => $registry,
    -output_file  => $file_name,
    -peak_calling => $peak_callings->[0],
  );

  $report->generate_report;
  
  return
}


1;
