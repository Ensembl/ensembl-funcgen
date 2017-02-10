#
# Ensembl module for Bio::EnsEMBL::Funcgen::DBSQL::CellTypeAdaptor
#

=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2017] EMBL-European Bioinformatics Institute

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

Bio::EnsEMBL::Funcgen::DBSQL::CellTypeAdaptor - A database adaptor for fetching and
storing Funcgen CellType objects.

=head1 SYNOPSIS

my $ct_adaptor = $efgdba->get_CellTypeAdaptor();

my $cell_type = $ct_adaptor->fetch_by_name("HeLa-S3");


=head1 DESCRIPTION

The CellTypeAdaptor is a database adaptor for storing and retrieving
Funcgen CellType objects.

=cut

package Bio::EnsEMBL::Funcgen::DBSQL::CellTypeAdaptor;

use strict;
use warnings;
use Bio::EnsEMBL::Utils::Exception qw( deprecate );

use base qw(Bio::EnsEMBL::Funcgen::DBSQL::EpigenomeAdaptor);

deprecate(
    "Bio::EnsEMBL::Funcgen::DBSQL::CellTypeAdaptor has been deprecated and will be removed in Ensembl release 89."
        . " Please use Bio::EnsEMBL::Funcgen::DBSQL::EpigenomeAdaptor instead."
);

1;
