package Bio::EnsEMBL::Funcgen::RunnableDB::Ftp::JobFactoryFtpRegulatoryActivityFiles;

use strict;
use Data::Dumper;
use Bio::EnsEMBL::Registry;
use base ('Bio::EnsEMBL::Hive::Process');

use Bio::EnsEMBL::Funcgen::Utils::FtpUtils qw (
  get_all_regulatory_activity_file_infos_from_ftp
);

sub run {
  my $self = shift;

  my $species = $self->param('species');
  my $ftp_dir = $self->param('ftp_base_dir');

  my $funcgen_adaptor = Bio::EnsEMBL::Registry->get_DBAdaptor( $species, 'funcgen');
  my $dbc = $funcgen_adaptor->dbc;
  my $regulatory_activity_file_name_infos = get_all_regulatory_activity_file_infos_from_ftp({
      dbc     => $dbc,
      ftp_dir => $ftp_dir . '/' . $species,
      species => $species,
  });

  foreach my $regulatory_activity_file_name_info (@$regulatory_activity_file_name_infos) {
  
    if (! defined $regulatory_activity_file_name_info->{file_name}) {
      confess("Can't find file name in " . Dumper($regulatory_activity_file_name_info));
    }
  
    $self->dataflow_output_id({
      species => $species,
      %$regulatory_activity_file_name_info
    }, 2);
  }
}

1;
