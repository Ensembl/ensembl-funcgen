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

package Bio::EnsEMBL::Funcgen::Test;

use strict;
use warnings;

use File::Basename;

use Test::More;
use Test::Exception;

use Bio::EnsEMBL::Test::MultiTestDB;
use Bio::EnsEMBL::Test::TestUtils qw(test_getter_setter);

use parent qw(Test::Class Class::Data::Inheritable);

BEGIN {
    __PACKAGE__->mk_classdata('full_class');
    __PACKAGE__->mk_classdata('short_class');
}

INIT {Test::Class->runtests}

sub db_connect :Test(startup) {
    my $self = shift;

    my ($filename, $dir, $suffix) = fileparse(__FILE__);
    my $multitestdb_conf_dir      = $dir . '../../../t';

    my $multi =
        Bio::EnsEMBL::Test::MultiTestDB->new('homo_sapiens',
                                             $multitestdb_conf_dir);

    $self->{multi}      = $multi;
    $self->{funcgen_db} = $multi->get_DBAdaptor('funcgen');
    $self->{core_db}    = $multi->get_DBAdaptor('core');
}

sub create_class_methods :Test(startup) {
    my $self = shift;

    my $full_class = 'Bio::EnsEMBL::Funcgen::Test';
    my $short_class;

    my $package_name = ref($self);
    if ($package_name ne __PACKAGE__) {
        $package_name =~ s/Test//;
        $short_class = $package_name;
        my $prefix   = 'Bio::EnsEMBL::Funcgen::';
        if ($short_class =~ /adaptor/) {
            $prefix .= 'DBSQL::'
        }
        $full_class = $prefix . $short_class;
    }
    $self->full_class($full_class);
    $self->short_class($short_class);
}

sub _creation :Test(1) {
    my $self = shift;

    use_ok($self->full_class)
        or $self->FAIL_ALL('Can not use ' . $self->full_class);
}

sub parameters :Test(setup) {
    my $self = shift;

    $self->{mandatory_constructor_parameters} = {};
}

sub test_getters_setters :Test(no_plan) {
    my $self = shift;

    my @getters_setters = @{$self->getters_setters};
    my @getters         = @{$self->getters};
    $self->num_tests(scalar @getters_setters + scalar @getters);

    my $short_class = $self->short_class();

    for my $getter_setter (@getters_setters) {
        ok(test_getter_setter($self->{fetched}->{$short_class},
                              $getter_setter,
                              'placeholder_string'),
           'getter/setter for ' . $getter_setter);
    }

    for my $getter (@getters) {
        is_deeply($self->{fetched}->{$short_class}->$getter,
                  $self->{expected}->{$getter},
                  'getter for ' . $getter . ' works'
        );
    }
}

sub constructor :Test(no_plan) {
    my $self = shift;

    my $full_class = $self->full_class;
    my %mandatory_parameters = %{$self->{mandatory_constructor_parameters}};

    $self->num_tests(scalar(keys(%mandatory_parameters)) + 1);

    my $new_object = $full_class->new(%{$self->{constructor_parameters}});
    isa_ok($new_object, $self->full_class);

    my %incomplete_parameters;
    my $error_message;
    for my $parameter (keys %mandatory_parameters) {
        %incomplete_parameters = %mandatory_parameters;
        delete $incomplete_parameters{$parameter};
        $error_message = "Must supply .* parameter";
        throws_ok {
            $full_class->new(%incomplete_parameters)
        }
            qr/$error_message/,
            "... and exception is thrown when $parameter parameter is missing";
    }
}

sub summary_as_hash :Test(no_plan) {
    my $self = shift;

    if (exists $self->{expected}->{summary}) {
        $self->num_tests(1);

        is_deeply($self->{fetched}->{$self->short_class}->summary_as_hash,
                  $self->{expected}->{summary},
                  'summary_as_hash() works');
    }
}

sub db_disconnect :Test(shutdown) {
    my $self = shift;
}

sub getters_setters {return []};

sub getters {return []};


1;