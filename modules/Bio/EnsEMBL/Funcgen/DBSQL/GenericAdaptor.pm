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

=head1 SYNOPSIS

=head1 DESCRIPTION

=cut

package Bio::EnsEMBL::Funcgen::DBSQL::GenericAdaptor;

use strict;
use warnings;
use Bio::EnsEMBL::Utils::Exception qw( throw warning );

use base (
  'Bio::EnsEMBL::DBSQL::BaseAdaptor',
  'Bio::EnsEMBL::Funcgen::GenericGetSetFunctionality',
);

sub new {
  my ($class, @args) = @_;
  
  # Calls constructor of Bio::EnsEMBL::DBSQL::BaseAdaptor
  my $self = $class->SUPER::new(@args);
  $self->init(@args);
  
  return $self;
}

sub init {
  my $self = shift;
  $self->SUPER::init($@);
  
  $self->_table_info_loader;
  
  use Bio::EnsEMBL::Utils::SqlHelper;
  my $sql_helper = Bio::EnsEMBL::Utils::SqlHelper->new(
    -DB_CONNECTION => $self->dbc
  );
  $self->sql_helper($sql_helper);
  
  # Otherwise the object can't be instantiated later.
  #
  eval "require " . $self->object_class;
  
  return $self;
}

sub _simple_accessors {
  return [
    { method_name => 'column_set',    hash_key => '_column_set',    },
    { method_name => 'primary_key',   hash_key => '_primary_key',   },
    { method_name => 'autoinc_id',    hash_key => '_autoinc_id',    },
    { method_name => 'sql_helper',    hash_key => '_sql_helper',    },
  ]
}

sub _columns {
  return []
}

sub _add_dependencies {}

sub default_input_column_mapping {
    return {
        # 'original_column1' => "original_column1*10 AS c1_times_ten",
        # 'original_column2' => "original_column2+1 AS c2_plus_one",
        # ...
    };
}

sub _load_dependencies  { return; }
sub _store_dependencies {}

sub insertion_method {
    return 'insert'
}

sub input_column_mapping {
    return {};
}

sub table_name {
    my $self   = shift;
    my @tables = $self->_tables;
    
    if (@tables > 1) {
        use Bio::EnsEMBL::Utils::Exception qw( throw );
        use Data::Dumper;
        throw("Only one table supported!" . Dumper(\@tables));
    }
    
    my $table_synonym_pair = $tables[0];
    
    my $table   = $table_synonym_pair->[0];
    my $synonym = $table_synonym_pair->[1];
    
    return $table;
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

    my $object = $self->object_class->new( 
      -db => $self, 
      map { 
        ( ($_ eq $autoinc_id) ? -dbID : '-' . $_ ) => $hashref->{$_} 
      } keys %$hashref 
    );
    return $object;
}

sub _table_info_loader {

    my $self = shift;

    my $table_name  = $self->table_name;
    my $dbc         = $self->db->dbc();
    my $dbh         = $dbc->db_handle();

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

sub fetch_all {
    my $self       = shift;
    my $constraint = shift;
    my $parameters = shift;

    my $table_name              = $self->table_name();
    my $input_column_mapping    = $self->input_column_mapping;
    
    my $sql = 'SELECT ' . join(', ', map { $input_column_mapping->{$_} // "$table_name.$_" } keys %{$self->column_set}) . " FROM $table_name";

    if($constraint) { 
        # in case $constraint contains any kind of JOIN (regular, LEFT, RIGHT, etc) do not put WHERE in front:
        $sql .= (($constraint=~/\bJOIN\b/i or $constraint=~/^LIMIT|ORDER|GROUP/) ? ' ' : ' WHERE ') . $constraint;
    }

    warn "SQL: $sql\n";

    my $sth = $self->prepare($sql);
    $sth->execute(@$parameters);

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

sub fetch_by_name {
    my $self = shift;
    my $name = shift;
    
    my $all = $self->fetch_all('name = ?', [ $name ]);
    return $all->[0];
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

sub keys_to_columns {
    my ($self, $object) = @_;

    my $autoinc_id  = $self->autoinc_id();
    my $sorted_keys = [ sort grep { ($_ ne $autoinc_id) and defined($object->$_()) } keys %{ $self->column_set } ];

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
            push $columns_being_stored, $self->primary_key->[0]
        }
        
        # warn "COLUMN_KEY='$column_key'\n";

        # only prepare (once!) if we get here:
        my $this_sth;
        $this_sth = $hashed_sth{$column_key};

        if (! defined $this_sth) {

            $sql = "$insertion_method INTO $table_name (".join(', ', @$columns_being_stored).') VALUES ('.join(',', (('?') x scalar(@$columns_being_stored))).')';
            #warn "STORE: $sql\n";
            $this_sth = $self->prepare( $sql );
            if (! defined $this_sth) {
                throw("Could not prepare statement: $sql");
            }
            $hashed_sth{$column_key} = $this_sth;
        }
        my $values_being_stored = $self->slicer( $object, $columns_being_stored );

        # using $return_code in boolean context allows to skip the 
        # value '0E0' ('no rows affected') that Perl treats as zero 
        # but regards as true:
        #
        my $return_code;
        
        eval {
          $return_code = $this_sth->execute( @$values_being_stored );
        };
        if ($@ || ! $return_code) {
           throw("Error executing sql:\n$sql\n\nCould not store fields\n\t{$column_key}\nwith data:\n\t(".join(',', @$values_being_stored).')');
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
