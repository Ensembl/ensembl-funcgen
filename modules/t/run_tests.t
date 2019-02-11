#!/usr/bin/perl

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

=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

  Questions may also be sent to the Ensembl help desk at
  <http://www.ensembl.org/Help/Contact>.

=cut

use strict;
use warnings;

use File::Basename;

# Load Tests
BEGIN {
    require Test::Class::Load;
    # Do not drop test databases, use the same for all Tests
    $ENV{'RUNTESTS_HARNESS'} = 1;

    # Load and run tests
    my ($filename, $dir, $suffix) = fileparse(__FILE__);
    my $test_dir = $dir . '/lib';
    Test::Class::Load->import($test_dir);
}

# Run Tests
# Test::Class->runtests();

# Clean up test databases and temp file(s)
delete $ENV{'RUNTESTS_HARNESS'};
my ($filename, $dir, $suffix) = fileparse(__FILE__);
my $db_conf = Bio::EnsEMBL::Test::MultiTestDB->get_db_conf($dir);

foreach my $species ( keys %{ $db_conf->{'databases'} } ) {
    # as soon as the $multi object goes out of scope, its DESTROY function
    # will take care of dropping the test databases
    my $multi = Bio::EnsEMBL::Test::MultiTestDB->new($species);
}

# CLEAN.t is created by MultiTestDB when environmental variable RUNTESTS_HARNESS
# is set to 1 and should be removed at the end of the run.
my $clean_file = $dir . 'CLEAN.t';
print "# Deleting $clean_file\n";
unlink $clean_file;