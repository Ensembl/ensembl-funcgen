package Bio::EnsEMBL::Funcgen::Utils::FtpUtils;
=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2019] EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=cut


=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

  Questions may also be sent to the Ensembl help desk at
  <http://www.ensembl.org/Help/Contact>.

=cut
use strict;
use warnings;
use base qw( Exporter );
use vars qw( @EXPORT_OK );

@EXPORT_OK = qw(
  create_file_handle
);

sub create_file_handle {

  my $output_file = shift;

  my $output_fh;
  if ($output_file) {

    use File::Basename;
    my $ftp_dir = dirname($output_file);

    use File::Path qw(make_path);
    make_path($ftp_dir);

    use IO::File;
    $output_fh = IO::File->new(">$output_file");
  } else {
    $output_fh = *STDOUT;
  }
  return $output_fh;
}

1;