#
# EnsEMBL module for Bio::EnsEMBL::DBSQL::CoordSystemAdaptor
#
#

=head1 NAME

Bio::EnsEMBL::DBSQL::CoordSystemAdaptor

=head1 SYNOPSIS

  my $db = Bio::EnsEMBL::DBSQL::DBAdaptor->new(...);

  my $csa = $db->get_CoordSystemAdaptor();

  #
  # Get all coord systems in the database:
  #
  foreach my $cs (@{$csa->fetch_all()}) {
    print $cs->name, ' ',  $cs->version, "\n";
  }

  #
  # Fetching by name:
  #

  #use the default version of coord_system 'chromosome' (e.g. NCBI33):
  $cs = $csa->fetch_by_name('chromosome');

  #get an explicit version of coord_system 'chromosome':
  $cs = $csa->fetch_by_name('chromsome', 'NCBI34');

  #get all coord_systems of name 'chromosome':
  foreach $cs (@{$csa->fetch_all_by_name('chromosome')}) {
     print $cs->name, ' ', $cs->version, "\n";
  }

  #
  # Fetching by rank:
  #
  $cs = $csa->fetch_by_rank(2);

  #
  # Fetching the pseudo coord system 'toplevel'
  #

  #Get the default top_level coord system:
  $cs = $csa->fetch_top_level();

  #can also use an alias in fetch_by_name:
  $cs = $csa->fetch_by_name('toplevel');

  #can also request toplevel using rank=0
  $cs = $csa->fetch_by_rank(0);

  #
  # Fetching by sequence level:
  #

  #Get the coord system which is used to store sequence:
  $cs = $csa->fetch_sequence_level();

  #can also use an alias in fetch_by_name:
  $cs = $csa->fetch_by_name('seqlevel');

  #
  # Fetching by id
  #
  $cs = $csa->fetch_by_dbID(1);


=head1 DESCRIPTION

This adaptor allows the querying of information from the coordinate system
adaptor.

Note that many coordinate systems do not have a concept of a version
for the entire coordinate system (though they may have a per-sequence version).
The 'chromosome' coordinate system usually has a version (i.e. the
assembly version) but the clonal coordinate system does not (despite having
individual sequence versions).  In the case where a coordinate system does
not have a version an empty string ('') is used instead.

=head1 AUTHOR - Graham McVicker

=head1 CONTACT

Post questions to the EnsEMBL development list ensembl-dev@ebi.ac.uk

=head1 METHODS

=cut

use strict;
use warnings;

package Bio::EnsEMBL::DBSQL::CoordSystemAdaptor;

use Bio::EnsEMBL::DBSQL::BaseAdaptor;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::CoordSystem;

use vars qw(@ISA);

@ISA = qw(Bio::EnsEMBL::DBSQL::BaseAdaptor);


=head2 new

  Arg [1]    : See BaseAdaptor for arguments (none specific to this
               subclass)
  Example    : $cs = $db->get_CoordSystemAdaptor(); #better than new()
  Description: Creates a new CoordSystem adaptor and caches the contents
               of the coord_system table in memory.
  Returntype : Bio::EnsEMBL::DBSQL::CoordSystemAdaptor
  Exceptions : none
  Caller     :
  Status     : Stable

=cut

sub new {
  my $caller = shift;

  my $class = ref($caller) || $caller;

  my $self = $class->SUPER::new(@_);

  #
  # Cache the entire contents of the coord_system table cross-referenced
  # by dbID and name
  #

  #keyed on name, list of coord_system value
  $self->{'_name_cache'} = {};

  #keyed on id, coord_system value
  $self->{'_dbID_cache'} = {};

  #keyed on rank
  $self->{'_rank_cache'} = {};

  #keyed on id, 1/undef values
  $self->{'_is_sequence_level'} = {};
  $self->{'_is_default_version'} = {};

  my $sth = $self->prepare(
    'SELECT coord_system_id, name, rank, version, attrib ' .
    'FROM coord_system');
  $sth->execute();

  my ($dbID, $name, $rank, $version, $attrib);
  $sth->bind_columns(\$dbID, \$name, \$rank, \$version, \$attrib);

  while($sth->fetch()) {
    my $seq_lvl = 0;
    my $default = 0;
    if($attrib) {
      foreach my $attrib (split(',', $attrib)) {
        $self->{"_is_$attrib"}->{$dbID} = 1;
        if($attrib eq 'sequence_level') {
          $seq_lvl = 1;
        } elsif($attrib eq 'default_version') {
          $default = 1;
        }
      }
    }

    my $cs = Bio::EnsEMBL::CoordSystem->new
      (-DBID           => $dbID,
       -ADAPTOR        => $self,
       -NAME           => $name,
       -VERSION        => $version,
       -RANK           => $rank,
       -SEQUENCE_LEVEL => $seq_lvl,
       -DEFAULT        => $default);

    $self->{'_dbID_cache'}->{$dbID} = $cs;

    $self->{'_name_cache'}->{lc($name)} ||= [];
    $self->{'_rank_cache'}->{$rank} = $cs;
    push @{$self->{'_name_cache'}->{lc($name)}}, $cs;
  }
  $sth->finish();



  #
  # Retrieve a list of available mappings from the meta table.
  # this may eventually be moved a table of its own if this proves too
  # cumbersome
  #

  my %mapping_paths;
  my $mc = $self->db()->get_MetaContainer();

 MAP_PATH:
  foreach my $map_path (@{$mc->list_value_by_key('assembly.mapping')}) {
    my @cs_strings = split(/[|#]/, $map_path);

    if(@cs_strings < 2) {
      warning("Incorrectly formatted assembly.mapping value in meta " .
              "table: $map_path");
      next MAP_PATH;
    }

    my @coord_systems;
    foreach my $cs_string (@cs_strings) {
      my($name, $version) = split(/:/, $cs_string);
      my $cs = $self->fetch_by_name($name, $version);
      if(!$cs) {
        warning("Unknown coordinate system specified in meta table " .
                " assembly.mapping:\n  $name:$version");
        next MAP_PATH;
      }
      push @coord_systems, $cs;
    }

    # if the delimiter is a # we want a special case, multiple parts of the same
    # componente map to same assembly part. As this looks like the "long" mapping
    # we just make the path a bit longer :-)

    if( $map_path =~ /\#/ && scalar( @coord_systems ) == 2 ) {
      splice( @coord_systems, 1, 0, ( undef ));
    }

    my $cs1 = $coord_systems[0];
    my $cs2  = $coord_systems[$#coord_systems];

    my $key1 = $cs1->name().':'.$cs1->version();
    my $key2 = $cs2->name().':'.$cs2->version();

    if(exists($mapping_paths{"$key1|$key2"})) {
      warning("Meta table specifies multiple mapping paths between " .
              "coord systems $key1 and $key2.\n" .
              "Choosing shorter path arbitrarily.");

      next MAP_PATH if(@{$mapping_paths{"$key1|$key2"}} < @coord_systems);
    }

    $mapping_paths{"$key1|$key2"} = \@coord_systems;
  }

  #
  # Create the pseudo coord system 'toplevel' and cache it so that
  # only one of these is created for each db...
  #
  my $toplevel = Bio::EnsEMBL::CoordSystem->new(-TOP_LEVEL => 1,
                                                -NAME      => 'toplevel',
                                                -ADAPTOR   => $self);
  $self->{'_top_level'} = $toplevel;

  $self->{'_mapping_paths'} = \%mapping_paths;

  return $self;
}



=head2 fetch_all

  Arg [1]    : none
  Example    : foreach my $cs (@{$csa->fetch_all()}) {
                 print $cs->name(), ' ', $cs->version(), "\n";
               }
  Description: Retrieves every coordinate system defined in the DB.
               These will be returned in ascending order of rank. I.e.
               The highest coordinate system with rank=1 would be first in the
               array.
  Returntype : listref of Bio::EnsEMBL::CoordSystems
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub fetch_all {
  my $self = shift;

  my @coord_systems;

  #order the array by rank in ascending order
  foreach my $rank (sort {$a <=> $b} keys %{$self->{'_rank_cache'}}) {
    push @coord_systems, $self->{'_rank_cache'}->{$rank};
  }

  return \@coord_systems;
}



=head2 fetch_by_rank

  Arg [1]    : int $rank
  Example    : my $cs = $coord_sys_adaptor->fetch_by_rank(1);
  Description: Retrieves a CoordinateSystem via its rank. 0 is a special
               rank reserved for the pseudo coordinate system 'toplevel'.
               undef is returned if no coordinate system of the specified rank
               exists.
  Returntype : Bio::EnsEMBL::CoordSystem
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub fetch_by_rank {
  my $self = shift;
  my $rank = shift;

  throw("Rank argument must be defined.") if(!defined($rank));
  throw("Rank argument must be a non-negative integer.") if($rank !~ /^\d+$/);

  if($rank == 0) {
    return $self->fetch_top_level();
  }

  return $self->{'_rank_cache'}->{$rank};
}


=head2 fetch_by_name

  Arg [1]    : string $name
               The name of the coordinate system to retrieve.  Alternatively
               this may be an alias for a real coordinate system.  Valid
               aliases are 'toplevel' and 'seqlevel'.
  Arg [2]    : string $version (optional)
               The version of the coordinate system to retrieve.  If not
               specified the default version will be used.
  Example    : $coord_sys = $csa->fetch_by_name('clone');
               $coord_sys = $csa->fetch_by_name('chromosome', 'NCBI33');
               # toplevel is an pseudo coord system representing the highest
               # coord system in a given region
               # such as the chromosome coordinate system
               $coord_sys = $csa->fetch_by_name('toplevel');
               #seqlevel is an alias for the sequence level coordinate system
               #such as the clone or contig coordinate system
               $coord_sys = $csa->fetch_by_name('seqlevel');
  Description: Retrieves a coordinate system by its name
  Returntype : Bio::EnsEMBL::CoordSystem
  Exceptions : throw if no name argument provided
               warning if no version provided and default does not exist
  Caller     : general
  Status     : Stable

=cut

sub fetch_by_name {
  my $self = shift;
  my $name = lc(shift); #case insensitve matching
  my $version = shift;

  throw('Name argument is required.') if(!$name);

  $version = lc($version) if($version);


  if($name eq 'seqlevel') {
    return $self->fetch_sequence_level();
  } elsif($name eq 'toplevel') {
    return $self->fetch_top_level($version);
  }

  if(!exists($self->{'_name_cache'}->{$name})) {
    if($name =~ /top/) {
      warning("Did you mean 'toplevel' coord system instead of '$name'?");
    } elsif($name =~ /seq/) {
      warning("Did you mean 'seqlevel' coord system instead of '$name'?");
    }
    return undef;
  }

  my @coord_systems = @{$self->{'_name_cache'}->{$name}};

  foreach my $cs (@coord_systems) {
    if($version) {
      return $cs if(lc($cs->version()) eq $version);
    } elsif($self->{'_is_default_version'}->{$cs->dbID()}) {
      return $cs;
    }
  }

  if($version) {
    #the specific version we were looking for was not found
    return undef;
  }

  #didn't find a default, just take first one
  my $cs =  shift @coord_systems;
  my $v = $cs->version();
  warning("No default version for coord_system [$name] exists. " .
      "Using version [$v] arbitrarily");

  return $cs;
}


=head2 fetch_all_by_name

  Arg [1]    : string $name
               The name of the coordinate system to retrieve.  This can be
               the name of an actual coordinate system or an alias for a
               coordinate system.  Valid aliases are 'toplevel' and 'seqlevel'.
  Example    : foreach my $cs (@{$csa->fetch_all_by_name('chromosome')}){
                 print $cs->name(), ' ', $cs->version();
               }
  Description: Retrieves all coordinate systems of a particular name
  Returntype : listref of Bio::EnsEMBL::CoordSystem objects
  Exceptions : throw if no name argument provided
  Caller     : general
  Status     : Stable

=cut

sub fetch_all_by_name {
  my $self = shift;
  my $name = lc(shift); #case insensitive matching

  throw('Name argument is required') if(!$name);

  if($name eq 'seqlevel') {
    return [$self->fetch_sequence_level()];
  } elsif($name eq 'toplevel') {
    return [$self->fetch_top_level()];
  }

  return $self->{'_name_cache'}->{$name} || [];
}





=head2 fetch_by_dbID

  Arg [1]    : int dbID
  Example    : $cs = $csa->fetch_by_dbID(4);
  Description: Retrieves a coord_system via its internal
               identifier, or undef if no coordinate system with the provided
               id exists.
  Returntype : Bio::EnsEMBL::CoordSystem or undef
  Exceptions : thrown if no coord_system exists for specified dbID
  Caller     : general
  Status     : Stable

=cut

sub fetch_by_dbID {
  my $self = shift;
  my $dbID = shift;

  throw('dbID argument is required') if(!$dbID);

  my $cs = $self->{'_dbID_cache'}->{$dbID};

  return undef if(!$cs);

  return $cs;
}



=head2 fetch_top_level

  Arg [1]    : none
  Example    : $cs = $csa->fetch_top_level();
  Description: Retrieves the toplevel pseudo coordinate system.
  Returntype : a Bio::EnsEMBL::CoordSystem object
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub fetch_top_level {
  my $self = shift;

  return $self->{'_top_level'};
}


=head2 fetch_sequence_level

  Arg [1]    : none
  Example    : ($id, $name, $version) = $csa->fetch_sequence_level();
  Description: Retrieves the coordinate system at which sequence
               is stored at.
  Returntype : Bio::EnsEMBL::CoordSystem
  Exceptions : throw if no sequence_level coord system exists at all
               throw if multiple sequence_level coord systems exists
  Caller     : general
  Status     : Stable

=cut

sub fetch_sequence_level {
  my $self = shift;

  my @dbIDs = keys %{$self->{'_is_sequence_level'}};

  throw('No sequence_level coord_system is defined') if(!@dbIDs);

  if(@dbIDs > 1) {
    throw('Multiple sequence_level coord_systems are defined.' .
          'Only one is currently supported');
  }

  return $self->{'_dbID_cache'}->{$dbIDs[0]};
}




=head2 get_mapping_path

  Arg [1]    : Bio::EnsEMBL::CoordSystem $cs1
  Arg [2]    : Bio::EnsEMBL::CoordSystem $cs2
  Example    : foreach my $cs @{$csa->get_mapping_path($cs1,$cs2);
  Description: Given two coordinate systems this will return a mapping path
               between them if one has been defined.  Allowed Mapping paths are
               explicitly defined in the meta table.  The following is an
               example:

          mysql> select * from meta where meta_key = 'assembly.mapping';
          +---------+------------------+--------------------------------------+
          | meta_id | meta_key         | meta_value                           |
          +---------+------------------+--------------------------------------+
          |      20 | assembly.mapping | chromosome:NCBI34|contig             |
          |      21 | assembly.mapping | clone|contig                         |
          |      22 | assembly.mapping | supercontig|contig                   |
          |      23 | assembly.mapping | chromosome:NCBI34|contig|clone       |
          |      24 | assembly.mapping | chromosome:NCBI34|contig|supercontig |
          |      25 | assembly.mapping | supercontig|contig|clone             |
          +---------+------------------+--------------------------------------+

               For a one-step mapping path to be valid there needs to be
               a relationship between the two coordinate systems defined in
               the assembly table.  Two step mapping paths work by building
               on the one-step mapping paths which are already defined.

               The first coordinate system in a one step mapping path must
               be the assembled coordinate system and the second must be
               the component.

               Example of use:
               my $cs1 = $cs_adaptor->fetch_by_name('contig');
               my $cs2 = $cs_adaptor->fetch_by_name('chromosome');

               my @path = @{$cs_adaptor->get_mapping_path($cs1,$cs2)};

               if(!@path) {
                 print "No mapping path.";
               }
               elsif(@path == 2) {
                 print "2 step mapping path.";
                 print "Assembled = " . $path[0]->name() . "\n";
                 print "Component = " . $path[1]->name() . "\n";
               } else {
                 print "Multi step mapping path\n";
               }

  Returntype : reference to a list of Bio::EnsEMBL::CoordSystem objects

  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub get_mapping_path {
  my $self = shift;
  my $cs1 = shift;
  my $cs2 = shift;

  if(!ref($cs1) || !ref($cs2) ||
     !$cs1->isa('Bio::EnsEMBL::CoordSystem') ||
     !$cs2->isa('Bio::EnsEMBL::CoordSystem')) {
    throw('Two Bio::EnsEMBL::CoordSystem arguments expected.');
  }

  my $key1 = $cs1->name() . ":" . $cs1->version();
  my $key2 = $cs2->name() . ":" . $cs2->version();

  my $path = $self->{'_mapping_paths'}->{"$key1|$key2"};

  return $path if($path);

  $path = $self->{'_mapping_paths'}->{"$key2|$key1"};

  if(!$path) {
    # No path was explicitly defined, but we might be able to guess a
    # suitable path.  We only guess for missing 2 step paths.

    my %mid1;
    my %mid2;

    foreach my $path (values(%{$self->{'_mapping_paths'}})) {
      next if(@$path != 2);

      my $match = undef;

      if($path->[0]->equals($cs1)) {
        $match = 1;
      } elsif($path->[1]->equals($cs1)) {
        $match = 0;
      }

      if(defined($match)) {
        my $mid = $path->[$match];
        my $midkey = $mid->name() . ':' . $mid->version();

        # is the same cs mapped to by other cs?
        if($mid2{$midkey}) {
          my $path = [$cs1,$mid,$cs2];
          $self->{'_mapping_paths'}->{"$key1|$key2"} = $path;
          $key1 =~ s/\:$//;
          $key2 =~ s/\:$//;
          $midkey =~ s/\:$//;
          warning("Using implicit mapping path between '$key1' and '$key2' " .
                  "coord systems.\n" .
                  "An explicit 'assembly.mapping' entry should be added " .
                  "to the meta table.\nExample: " .
                  "'$key1|$midkey|$key2'\n");
          return $path;
        } else {
          $mid1{$midkey} = $mid;
        }
      }

      $match = undef;

      if($path->[0]->equals($cs2)) {
        $match = 1;
      } elsif($path->[1]->equals($cs2)) {
        $match = 0;
      }


      if(defined($match)) {
        my $mid = $path->[$match];
        my $midkey = $mid->name() . ':' . $mid->version();

        # is the same cs mapped to by other cs?
        if($mid1{$midkey}) {
          my $path = [$cs2,$mid,$cs1];
          $self->{'_mapping_paths'}->{"$key2|$key1"} = $path;

          $key1 =~ s/\:$//;
          $key2 =~ s/\:$//;
          $midkey =~ s/\:$//;
          warning("Using implicit mapping path between '$key1' and '$key2' " .
                  "coord systems.\n" .
                  "An explicit 'assembly.mapping' entry should be added " .
                  "to the meta table.\nExample: " .
                  "'$key1|$midkey|$key2'\n");

          return $path;
        } else {
          $mid2{$midkey} = $mid;
        }
      }
    }
  }

  return $path || [];
}


sub _fetch_by_attrib {
  my $self = shift;
  my $attrib = shift;
  my $version = shift;

  $version = lc($version) if($version);

  my @dbIDs = keys %{$self->{"_is_$attrib"}};

  throw("No $attrib coordinate system defined") if(!@dbIDs);

  foreach my $dbID (@dbIDs) {
    my $cs = $self->{'_dbID_cache'}->{$dbID};
    if($version) {
      return $cs if(lc($version) eq $cs->version());
    } elsif($self->{'_is_default_version'}->{$dbID}) {
      return $cs;
    }
  }

  #specifically requested attrib system was not found
  if($version) {
    throw("$attrib coord_system with version [$version] does not exist");
  }

  #coordsystem with attrib exists but no default is defined:
  my $dbID = shift @dbIDs;
  my $cs = $self->{'_dbID_cache'}->{$dbID};
  my $v = $cs->version();
  warning("No default version for $attrib coord_system exists. " .
          "Using version [$v] arbitrarily");

  return $cs;
}


sub _fetch_all_by_attrib {
  my $self = shift;
  my $attrib = shift;

  my @coord_systems = ();
  foreach my $dbID (keys %{$self->{"_is_$attrib"}}) {
    push @coord_systems, $self->{"_dbID_cache"}->{$dbID};
  }

  return \@coord_systems;
}


=head2 store

  Arg [1]    : Bio::EnsEMBL::CoordSystem
  Example    : $csa->store($coord_system);
  Description: Stores a CoordSystem object in the database.
  Returntype : none
  Exceptions : Warning if CoordSystem is already stored in this database.
  Caller     : none
  Status     : Stable

=cut

sub store {
  my $self = shift;
  my $cs = shift;

  if(!$cs || !ref($cs) || !$cs->isa('Bio::EnsEMBL::CoordSystem')) {
    throw('CoordSystem argument expected.');
  }

  my $db = $self->db();
  my $name = $cs->name();
  my $version = $cs->version();
  my $rank    = $cs->rank();

  my $seqlevel = $cs->is_sequence_level();
  my $default  = $cs->is_default();

  my $toplevel = $cs->is_top_level();

  if($toplevel) {
    throw("The toplevel CoordSystem cannot be stored");
  }

  #
  # Do lots of sanity checking to prevent bad data from being entered
  #

  if($cs->is_stored($db)) {
    warning("CoordSystem $name $version is already in db.\n");
    return;
  }

  if($name eq 'toplevel' || $name eq 'seqlevel' || !$name) {
    throw("[$name] is not a valid name for a CoordSystem.");
  }

  if($seqlevel && keys(%{$self->{'_is_sequence_level'}})) {
    throw("There can only be one sequence level CoordSystem.");
  }

  if(exists $self->{'_name_cache'}->{lc($name)}) {
    my @coord_systems = @{$self->{'_name_cache'}->{lc($name)}};
    foreach my $c (@coord_systems) {
      if(lc($c->version()) eq lc($version)) {
        warning("CoordSystem $name $version is already in db.\n");
        return;
      }
      if($default && $self->{'_is_default_version'}->{$c->dbID()}) {
        throw("There can only be one default version of CoordSystem $name");
      }
    }
  }

  if($rank !~ /^\d+$/) {
    throw("Rank attribute must be a positive integer not [$rank]");
  }
  if($rank == 0) {
    throw("Only toplevel CoordSystem may have rank of 0.");
  }

  if(defined($self->{'_rank_cache'}->{$rank})) {
    throw("CoordSystem with rank [$rank] already exists.");
  }

  my @attrib;

  push @attrib, 'default_version' if($default);
  push @attrib, 'sequence_level' if($seqlevel);

  my $attrib_str = (@attrib) ? join(',', @attrib) : undef;

  #
  # store the coordinate system in the database
  #

  my $sth = $db->dbc->prepare('INSERT INTO coord_system ' .
                         'SET name    = ?, ' .
                             'version = ?, ' .
                             'attrib  = ?,' .
                             'rank    = ?');

  $sth->bind_param(1,$name,SQL_VARCHAR);
  $sth->bind_param(2,$version,SQL_VARCHAR);
  $sth->bind_param(3,$attrib_str,SQL_VARCHAR);
  $sth->bind_param(4,$rank,SQL_INTEGER);
  $sth->execute();
  my $dbID = $sth->{'mysql_insertid'};
  $sth->finish();

  if(!$dbID) {
    throw("Did not get dbID from store of CoordSystem.");
  }

  $cs->dbID($dbID);
  $cs->adaptor($self);

  #
  # update the internal caches that are used for fetching
  #
  $self->{'_is_default_version'}->{$dbID} = 1 if($default);
  $self->{'_is_sequence_level'}->{$dbID} = 1 if($seqlevel);

  $self->{'_name_cache'}->{lc($name)} ||= [];
  push @{$self->{'_name_cache'}->{lc($name)}}, $cs;

  $self->{'_dbID_cache'}->{$dbID} = $cs;
  $self->{'_rank_cache'}->{$rank} = $cs;

  return $cs;
}


1;
