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
use Test::More;
use Test::Class::Load;

# Load all Test Classes
my ($filename, $dir, $suffix) = fileparse(__FILE__);
my $test_dir                  = $dir . '/lib';
Test::Class::Load->import($test_dir);

# Create test database
my $multi = Bio::EnsEMBL::Test::MultiTestDB->new('homo_sapiens');

# Run all available Test Classes or just the ones provided by the user
my @test_classes_to_run = @ARGV? @ARGV: Test::Class->_test_classes();

# Plan number of Tests
my $total_number_of_tests = 0;
my @tests;
for my $test_class (@test_classes_to_run){
    if ($test_class =~ m/Test::Class/){
        next;
    }
    my $test = $test_class->new($multi);
    push @tests, $test;
    $total_number_of_tests += $test->expected_tests();
}
plan tests => $total_number_of_tests;

# Run Tests
for my $t (@tests){
    $t->runtests();
}