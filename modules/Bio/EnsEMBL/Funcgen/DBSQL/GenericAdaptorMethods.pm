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


=head1 NAME

=head1 SYNOPSIS

=head1 DESCRIPTION

=cut

package Bio::EnsEMBL::Funcgen::DBSQL::GenericAdaptorMethods;

use strict;
use Role::Tiny;
use Bio::EnsEMBL::Utils::Exception qw( throw );
use Bio::EnsEMBL::Funcgen::GenericGetSetFunctionality qw(
  _generic_get_or_set
);

# Attributes
#
sub column_set  { return shift->_generic_get_or_set('column_set',  @_); }
sub primary_key { return shift->_generic_get_or_set('primary_key', @_); }
sub autoinc_id  { return shift->_generic_get_or_set('autoinc_id',  @_); }
sub sql_helper  { return shift->_generic_get_or_set('sql_helper',  @_); }

# Hooks
#
sub _columns {
  return []
}
sub _add_dependencies {
}
sub _load_dependencies  {
  return;
}
sub _store_dependencies {}
sub insertion_method {
    return 'insert'
}
sub input_column_mapping {
    return {};
}

sub fetch_by_dbID {
    my $self = shift;
    my $dbID = shift;

    my $primary_key = $self->primary_key;

    if (@$primary_key!=1) {
      die;
    }

    my $all = $self->fetch_all($primary_key->[0] . ' = ? ', [ $dbID ]);
    return $all->[0];
}

sub fetch_by_name {
    my $self = shift;
    my $name = shift;

    my $all = $self->fetch_all('name = ?', [ $name ]);
    return $all->[0];
}

requires 'object_class';

sub init_generic_adaptor {
  my $self = shift;

  $self->_table_info_loader;

  use Bio::EnsEMBL::Utils::SqlHelper;
  my $sql_helper = Bio::EnsEMBL::Utils::SqlHelper->new(
    -DB_CONNECTION => $self->dbc
  );
  $self->sql_helper($sql_helper);

  my $class = ref $self;

  # Otherwise the object can't be instantiated later.
  #
  eval "require " . $self->object_class;
};

sub _table_alias_pair {
    my $self   = shift;
    my @tables = $self->_tables;

    if (@tables > 1) {

        use Data::Dumper;
        throw("Only one table supported!" . Dumper(\@tables));
    }

    my $table_alias_pair = $tables[0];

    my $table = $table_alias_pair->[0];
    my $alias = $table_alias_pair->[1];

    return ($table, $alias);
}

sub table_name {
    my $self   = shift;
    (
      my $table_name,
      my $alias
    ) = $self->_table_alias_pair;
    return $table_name;
}

sub table_alias {
    my $self   = shift;
    (
      my $table_name,
      my $alias
    ) = $self->_table_alias_pair;
    return $alias;
}

sub prepare {
    my ( $self, $sql ) = @_;

    # Uncomment next line to cancel caching on the SQL side.
    # Needed for timing comparisons etc.
    #$sql =~ s/SELECT/SELECT SQL_NO_CACHE/i;

    return $self->db->dbc->prepare($sql);
}

sub objectify { # turn the hashref into an object
    my $self    = shift;
    my $hashref = shift;

    my $autoinc_id = $self->autoinc_id;
    
    my $object;
    eval {
    
      my %constructor_parameters = (
      
        # Providing both 'db' and 'adaptor' until all API objects have been
        # switched to the latter.
        #
        -db      => $self, 
        -adaptor => $self, 
        map {
          ( ($_ eq $autoinc_id) ? -dbID : '-' . $_ ) => $hashref->{$_} 
        } keys %$hashref 
      );
      
      $object = $self->object_class->new(
        %constructor_parameters
      );
    };
    if ($@) {
      my $error = $@;
      if ($error =~ /Can't locate object method "new" via package/) {
        my $error_message
          = "Can't instantiate object of type " . $self->object_class . ".\n"
          . "This can happen, if the object was misspelt or if the module"
          . " can't be compiled.\n"
          . "Try running\n\n"
          . "   perl -cw <Filename of the module ".$self->object_class.">\n\n"
          . "to find compilation problems."
        ;
        throw($error_message);
      }
      #throw($error);
      throw(
        "An unknown error occurred:\n\n"
        . "\t$error"
      );
    }
    return $object;
}

sub _table_info_loader {

    my $self = shift;

    my $table_name  = $self->table_name;
    my $db          = $self->db;

    if (! defined $db) {
      throw("db is undefined!");
    }
    if (! $db->isa('Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor')) {
      throw("Type error, got a " . (ref $db) );
    }

    my $dbc         = $db->dbc;
    my $dbh         = $dbc->db_handle;

    my %column_set  = ();
    my $autoinc_id  = '';
    my @primary_key = $dbh->primary_key(undef, undef, $table_name);
    my $sth         = $dbh->column_info(undef, undef, $table_name, '%');
    $sth->execute;

    while (my $row = $sth->fetchrow_hashref) {

        my ( $column_name, $column_type ) = @$row{'COLUMN_NAME', 'TYPE_NAME'};

        #warn "ColumnInfo [$table_name/$column_name] = $column_type\n";
        $column_set{$column_name}  = $column_type;

        if ($column_name eq $table_name.'_id') {
            $autoinc_id = $column_name;
        }
    }
    $sth->finish;

    $self->column_set(  \%column_set );
    $self->primary_key( \@primary_key );
    $self->autoinc_id(   $autoinc_id );
}

sub count_all {
    my ($self, $constraint, $key_list) = @_;

    my $table_name      = $self->table_name();
    my $count_col_name  = 'COUNT(*)';

    my $sql = "SELECT ".($key_list ? join(', ', @$key_list, '') : '')."COUNT(*) FROM $table_name";

    if($constraint) {
            # in case $constraint contains any kind of JOIN (regular, LEFT, RIGHT, etc) do not put WHERE in front:
        $sql .= (($constraint=~/\bJOIN\b/i) ? ' ' : ' WHERE ') . $constraint;
    }

    if($key_list) {
        $sql .= " GROUP BY ".join(', ', @$key_list);
    }
    #warn "SQL: $sql\n";

    my $sth = $self->prepare($sql);
    $sth->execute;

    my $result_struct;  # will be autovivified to the correct data structure

    while(my $hashref = $sth->fetchrow_hashref) {

        my $pptr = \$result_struct;
        if($key_list) {
            foreach my $syll (@$key_list) {
                $pptr = \$$pptr->{$hashref->{$syll}};   # using pointer-to-pointer to enforce same-level vivification
            }
        }
        $$pptr = $hashref->{$count_col_name};
    }

    unless(defined($result_struct)) {
        if($key_list and scalar(@$key_list)) {
            $result_struct = {};
        } else {
            $result_struct = 0;
        }
    }
    return $result_struct;
}

sub fetch_single_object {
    my $self = shift;
    my $constraint = shift;
    my $parameters = shift;

    my $result = $self->fetch_all($constraint, $parameters, @_);
    
    if (@$result == 0) {
      return;
    }
    if (@$result > 1) {
      throw(
        "Expected one result, but got " . scalar @$result . " "
        . "for the constraint:\n\n\t" . $constraint . "\n\n"
        . "and parameters:\n\n\t" . Dumper($parameters) . "\n"
      );
    }
    return $result->[0];
}

sub _generate_sql {
    my $self       = shift;
    my $constraint = shift;

    my $table_name  = $self->table_name;
    my $table_alias = $self->table_alias;
    
    my $input_column_mapping    = $self->input_column_mapping;

    my $sql = 'SELECT ' . join(', ', map { $input_column_mapping->{$_} // "$table_alias.$_" } keys %{$self->column_set}) . " FROM $table_name $table_alias";
    
    if($constraint) {
        # in case $constraint contains any kind of JOIN (regular, LEFT, RIGHT, etc) do not put WHERE in front:
        $sql .= (($constraint=~/\bJOIN\b/ or $constraint=~/^LIMIT|ORDER|GROUP/) ? ' ' : ' WHERE ') . $constraint;
    }
    return $sql;
}

sub fetch_all {
    my $self       = shift;
    my $constraint = shift;
    my $parameters = shift;

    my $table_name              = $self->table_name();
    my $input_column_mapping    = $self->input_column_mapping;
    
    my $sql = $self->_generate_sql($constraint);
    my $sth = $self->prepare($sql);

    eval {
      $sth->execute(@$parameters);
    };
    if ($@) {
      throw(
      "Error executing sql:\n\n"
      . "\t$sql\n\n"
      . "With these parameters:\n\n"
      . "\t(" . join(',', @$parameters) . ")\n\n"
      . "The error message is:\n\n"
      . "\t$@\n"
      );
    }

    my $result_list_ref = [];
    while (my $hashref = $sth->fetchrow_hashref) {
        my $object = $self->objectify($hashref);
        push @$result_list_ref, $object;
    }
    $sth->finish;

    if (! defined($result_list_ref)) {
        $result_list_ref = [];
    }
    foreach my $current_object (@$result_list_ref) {
      $self->_load_dependencies($current_object);
    }
    
    return $result_list_ref;
}

sub keys_to_columns {
    my ($self, $object) = @_;

    my $autoinc_id  = $self->autoinc_id();
    #my $sorted_keys = [ sort grep { ($_ ne $autoinc_id) and defined($object->$_()) } keys %{ $self->column_set } ];

    my $sorted_keys = [ sort grep { ($_ ne $autoinc_id) and defined($object->can($_)) } keys %{ $self->column_set } ];

    return ( $sorted_keys, join(', ', @$sorted_keys) );
}

sub slicer {    # take a slice of the object
    my ($self, $object, $fields) = @_;

    my $autoinc_id      = $self->autoinc_id();

    return [
      map {
        ($_ eq $autoinc_id)
            ? $object->dbID
            : eval {
                my $value = $object->$_;
                $value
            }
      } @$fields
    ];
}

sub mark_stored {
    my ($self, $object, $dbID) = @_;

    if($self->autoinc_id()) {
        $object->dbID($dbID);
    }
    $object->db($self);
}

sub update {
  my $self = shift;
  my $object = shift;

  my (
    $columns_being_stored, 
    $column_key
  ) 
    = $self->keys_to_columns($object);

  my @update_component;
  foreach my $current_column (@$columns_being_stored) {
    push 
      @update_component, 
      "$current_column = ?"
    ;
  }
  my $update_assignment_part = join ', ', @update_component;

  my $values_being_stored = [
    @{$self->slicer( $object, $columns_being_stored )},
    $object->dbID
  ];
  my $table_name = $self->table_name;
  my $primary_key = $self->primary_key->[0];
  
  if (! defined $primary_key) {
    throw("Can't update, primary key is undefined!");
  }
  
  my $sql = "update $table_name set $update_assignment_part where $primary_key = ?" ;
  
#   warn($sql);
#   warn(Dumper($values_being_stored));
  
  my $sth = $self->prepare( $sql );
  eval {
    $sth->execute(@$values_being_stored);
  };
  if ($@) {
    throw(
        "An error occurred runing this sql:\n"
        . $sql . "\n"
        . "\n"
        . "The error was:\n"
        . $@
    );
  }
  return;
}

sub _delete {
  my $self = shift;
  my $object = shift;

  my $table_name = $self->table_name;
  my $primary_key = $self->primary_key->[0];
  
  if (! defined $primary_key) {
    throw("Can't update, primary key is undefined!");
  }
  
  my $sql = "delete from $table_name where $primary_key = ?" ;
  
#   warn($sql);
#   warn($object->dbID);
  my $sth = $self->prepare( $sql );
  $sth->execute($object->dbID);
  return;
}

sub store {
    my ($self, $object_or_list) = @_;

    # Eensure we get an array of objects to store.
    #
    my $objects = (ref($object_or_list) eq 'ARRAY')
        ? $object_or_list
        : [ $object_or_list ];

    return ([], 0) unless(scalar(@$objects));

    my $table_name              = $self->table_name;
    my $autoinc_id              = $self->autoinc_id;
    my $all_storable_columns    = [ grep { $_ ne $autoinc_id } keys %{ $self->column_set } ];

    # INSERT, INSERT_IGNORE or REPLACE
    my $insertion_method        = $self->insertion_method;

    # Do not prepare statements until there is a real need.
    my %hashed_sth = ();

    my $stored_this_time        = 0;
    my $sql;

    foreach my $object (@$objects) {
        my ($columns_being_stored, $column_key) = $self->keys_to_columns($object);

        # If the dbID is defined, store this too. This is useful when
        # creating the test database.
        #
        if (defined $object->dbID) {
            push @$columns_being_stored, $self->primary_key->[0]
        }

        # warn "COLUMN_KEY='$column_key'\n";

        # only prepare (once!) if we get here:
        my $this_sth;
        $this_sth = $hashed_sth{$column_key};

        if (! defined $this_sth) {

            $sql = "$insertion_method INTO $table_name (".join(', ', @$columns_being_stored).') VALUES ('.join(',', (('?') x scalar(@$columns_being_stored))).')';
            
            $this_sth = $self->prepare( $sql );
            if (! defined $this_sth) {
                throw("Could not prepare statement: $sql");
            }
            $hashed_sth{$column_key} = $this_sth;
        }
        my $values_being_stored = $self->slicer( $object, $columns_being_stored );

#         warn "STORE: $sql\n";
#         warn "STORE: ". join(',', @$values_being_stored) ."\n";

        # using $return_code in boolean context allows to skip the
        # value '0E0' ('no rows affected') that Perl treats as zero
        # but regards as true:
        #
        my $return_code;
        eval {
          $return_code = $this_sth->execute( @$values_being_stored );
        };
        if ($@ || ! $return_code) {
           throw(
            "Error executing sql:\n\n"
            . "\t$sql\n\n"
            . "Could not store this data set:\n\n"
            . "\t(" . join(',', @$values_being_stored) . ")\n\n"
            . "into these columns:\n\n"
            . "\t{$column_key}\n\n"
            . "for " . (ref $self) . "\n\n"
            . "The error message is:\n\n"
            . "\t$@\n"
           );
        }

        # For the same reason we have to be explicitly numeric here
        #
        if($return_code > 0) {
            my $liid = $autoinc_id && $self->dbc->db_handle->last_insert_id(undef, undef, $table_name, $autoinc_id);
            $self->mark_stored($object, $liid );
            ++$stored_this_time;
        }
    }
    foreach my $sth (values %hashed_sth) {
        $sth->finish();
    }
    foreach my $current_object (@$objects) {
      $self->_store_dependencies($current_object);
    }
    return ($objects, $stored_this_time);
}

1;
