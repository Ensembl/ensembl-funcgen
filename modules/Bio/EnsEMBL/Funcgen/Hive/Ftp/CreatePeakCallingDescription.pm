package Bio::EnsEMBL::Funcgen::Hive::Ftp::CreatePeakCallingDescription;

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
  
  my $output;
  
  use Template;
  my $tt = Template->new;

  use Bio::EnsEMBL::Funcgen::Template::PeakCallingDescription qw( 
    PEAK_CALLING_TXT_TEMPLATE
  );

  my $description_template = PEAK_CALLING_TXT_TEMPLATE;
  use Number::Format qw( :all );
  
  open my $fh, ">", $file_name;
  for my $peak_calling (@$peak_callings) {
$tt->process(
  \$description_template, 
  {
    peak_calling  => $peak_calling,
    
    canonpath => sub {
      my $path = shift;
      return File::Spec->canonpath($path)
    },
    
    bool_to_yes_no => sub {
      my $boolean = shift;
      if ($boolean) {
        return 'yes'
      }
      return 'no'
    },
    
    round_percent => sub {
      my $number = shift;
      return sprintf("%.2f", $number) . '%';
    },
    default_round => sub {
      my $number = shift;
      return sprintf("%.2f", $number);
    },
    scientific_notation => sub {
      my $number = shift;
      return sprintf("%.2e", $number);
    },
    format_number => sub {
      my $number = shift;
      return format_number($number);
    }
  },
  \$output
)
    || die $tt->error;
  }
  $fh->print($output);
  $fh->close;
  return
}

1;
