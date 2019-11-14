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
    __PACKAGE__->mk_classdata('class');
}

sub new {
    my $class = shift;
    my $self = $class->SUPER::new(@_);

    if (! exists $self->{multi}){
        $self->db_connect();
    }
    $self->create_class_methods();
    $self->fetch_from_test_db();
    $self->parameters();
    $self->define_expected();

    my $i = scalar @{$self->getters_setters} + scalar @{$self->getters};
    $self->num_method_tests('test_getters_setters', $i);

    my $j = 1;
    if ($self->{mandatory_constructor_parameters}){
        my %mcps = %{$self->{mandatory_constructor_parameters}};
        $j += scalar(keys(%mcps));
    }
    $self->num_method_tests('constructor', $j);

    my $n = 0;
    if (exists $self->{expected}->{summary}){
        $n = 1;
    }
    $self->num_method_tests('summary_as_hash', $n);

    return $self;
}

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
    my $class;

    my $package_name = ref($self);
    if ($package_name ne __PACKAGE__) {
        $package_name =~ s/Test//;
        $class     = $package_name;
        my $prefix = 'Bio::EnsEMBL::Funcgen::';
        if ($class =~ /Adaptor/) {
            $prefix .= 'DBSQL::'
        }
        $full_class = $prefix . $class;
    }
    $self->full_class($full_class);
    $self->class($class);
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

sub define_expected :Test(setup){
    my $self = shift;
    $self->{expected} = {};
}

sub fetch_from_test_db :Test(setup) {
    my $self = shift;

    my $dbIDs = $self->dbIDs_to_fetch;
    my $class = $self->class();

    if (scalar @{$dbIDs} > 0) {
        $self->{fetched} =
            $self->{funcgen_db}->get_adaptor($class)
                 ->fetch_all_by_dbID_list($dbIDs);
    }

}

sub test_getters_setters :Test(no_plan) {
    my $self = shift;

    my @getters_setters = @{$self->getters_setters};
    my @getters         = @{$self->getters};
    $self->num_tests(scalar @getters_setters + scalar @getters);

    my $class = $self->class();

    for my $getter_setter (@getters_setters) {
        ok(test_getter_setter($self->{fetched}->[0],
                              $getter_setter,
                              'placeholder_string'),
           $getter_setter . '() getter/setter works'
        );
    }

    for my $getter (@getters) {
        is_deeply($self->{fetched}->[0]->$getter,
                  $self->{expected}->{$getter},
                  $getter . '() getter works'
        );
        1;
    }
}

sub constructor :Test(no_plan) {
    my $self = shift;

    my $full_class           = $self->full_class;

    if ($self->{mandatory_constructor_parameters}){
        my %mandatory_parameters = %{$self->{mandatory_constructor_parameters}};

        $self->num_tests(scalar(keys(%mandatory_parameters)) + 1);

        my $new_object = $full_class->new(%{$self->{constructor_parameters}});
        isa_ok($new_object, $self->full_class);

        my %missing_parameters;
        my $error_message;
        for my $parameter (keys %mandatory_parameters) {
            %missing_parameters = %mandatory_parameters;
            delete $missing_parameters{$parameter};
            $error_message = "Must supply a .* parameter";
            throws_ok {
                $full_class->new(%missing_parameters)
            }
                qr/$error_message/,
                "... and exception is thrown when $parameter parameter is missing";
        }
    }
}

sub summary_as_hash :Test(no_plan) {
    my $self = shift;

    if (exists $self->{expected}->{summary}) {
        $self->num_tests(1);

        is_deeply($self->{fetched}->[0]->summary_as_hash,
                  $self->{expected}->{summary},
                  'summary_as_hash() works');
    }
}

sub db_disconnect :Test(shutdown) {
    my $self = shift;
}

sub getters_setters {return []};

sub getters {return []};

sub dbIDs_to_fetch {return [];}

sub _quick_fetch {
    my ($self, $class, $dbID, $dbtype) = @_;

    my $default_dbtype = 'funcgen';
    if(! $dbtype){
        $dbtype = $default_dbtype;
    }
    $dbtype .= '_db';

    my $adaptor = $self->{$dbtype}->get_adaptor($class);
    return $adaptor->fetch_by_dbID($dbID);
}

sub _quick_fetch_all {
    my ($self, $class, $dbIDs, $dbtype) = @_;

    my $default_dbtype = 'funcgen';
    if(! $dbtype){
        $dbtype = $default_dbtype;
    }
    $dbtype .= '_db';

    my $adaptor = $self->{$dbtype}->get_adaptor($class);
    return $adaptor->fetch_all_by_dbID_list($dbIDs);
}

1;