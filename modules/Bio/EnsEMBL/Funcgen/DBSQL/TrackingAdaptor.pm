#
# Ensembl module for Bio::EnsEMBL::DBSQL::Funcgen::TrackingAdaptor
#

=head1 LICENSE

Copyright [1999-2014] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute

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

Bio::EnsEMBL::DBSQL::Funcgen::TrackingAdaptor - Handles interactions with the funcgen tracking DB

=head1 SYNOPSIS

my $tdb = Bio::EnsEMBL::DBSQL::Funcgen::DBSQL::TrackingAdaptor->new
  (
   -user     => $tdb_user,
   -pass     => $tdb_pass,
   -dbname   => $tdb_name,
   -host     => $tdb_host,
   -species  => $species,
  );



=head1 DESCRIPTION

Handles interaction with the funcgen tracking DB is used to store meta data and tracking information
about data sets. This provides a staging area where analysis and QC can be performed prior to 
migration to the DB destined for a particular release. This adaptor methods to store and access
these information given the appropriately related objects e.g.
    InputSet, InputSubset, ResultSet, FeatureSet & DataSet

=cut

package Bio::EnsEMBL::Funcgen::DBSQL::TrackingAdaptor;

use strict;
use warnings;
use DateTime;
#use POSIX                                  qw( strftime );
use Bio::EnsEMBL::Utils::Exception         qw( throw warning );
use Bio::EnsEMBL::Utils::Scalar            qw( check_ref assert_ref );
use Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Funcgen::Utils::EFGUtils qw( get_month_number split_CamelCase);

use base qw( Bio::EnsEMBL::DBSQL::BaseAdaptor );
#Currently don't use and Funcgen BaseAdaptor methods in here
#and this would require _true_tables to be defined

#values are mandatory booleans
my %tracking_columns = 
 (input_subset => {availability_date => 1,
                   download_url      => 1,
                   download_date     => 0,
                   local_url         => 0,
                   md5sum            => 0,
                   notes             => 0},
  result_set   => {idr_max_peaks        => 0,
                   idr_peak_analysis_id => 0},
 );
                        

=head2 new

  Arg [-DNADB_HOST]     : String - Overrides defaults (ensembldb or registry)
  Arg [-DNADB_USER]     : String - Overrides defaults (ensembldb or registry)
  Arg [-DNADB_PASS]     : String - Overrides defaults (ensembldb or registry)
  Arg [-DNADB_PORT]     : String - Overrides defaults (ensembldb or registry)
  Arg [-DNADB_NAME]     : String - Overrides defaults (ensembldb or registry)
  Arg [-DNADB_ASSEMBLY] : String - Overrides defaults (ensembldb or registry)

  Arg [...]         : Other args are passed to superclass Bio::EnsEMBL::DBSQL::DBAdaptor
  Example1          : $db = new Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor
                              (
						       -user   => 'readonly', #No password
						       -dbname => 'pog',
						       -host   => 'caldy',
						    );

                      #If dnadb is not defined in registry, this will automatically
                      #set it from the default dnadb_host (e.g. ensembldb)

  Exmaple2          : $db = new Bio::EnsEMBL::DBSQL::Funcgen::DBAdaptor
                              (
						       -user           => 'write',
                               -pass           => 'password',
						       -dbname         => 'pog',
						       -host           => 'caldy',
                               -dnadb_assmebly => '36',
						      );
                       #This will specifically look for a dnadb with assembly version 36
                       #on the default dnadb_host

   Exmaple2          : $db = new Bio::EnsEMBL::DBSQL::Funcgen::DBAdaptor
                              (
						       -user           => 'write',
                               -pass           => 'password',
						       -dbname         => 'pog',
						       -host           => 'caldy',
                               #The following will over-ride the default dnadb setting
                               -dnadb_name     => 'my_homo_sapiens_core_67_37',
                               -dnadb_host     => 'my_host',
                               #Can add more dnadb params if required 
                               #-dnadb_host
                               #-dnadb_user
                               #-dnadb_port
						      );

  Description: Constructor for DBAdaptor. Will automatically set the dnadb based on dnadb params.
               This makes some assumptions about how the DBs name are defined i.e. the last part must
               conform to the _RELEASE_ASSEMBLY format e.g. 67_37 in homo_sapiens_funcgen_67_37
  Returntype : Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor
  Exceptions : Throws if conflicting dnadb params found
  Caller     : general
  Status     : Stable

=cut


# We don't want this to be available via $efgdb->get_TrackingAdaptor
# So have to handle setting db
# alternatively just pass DBAdaptor here? 

sub new {
  my ($class, @args) = @_; 
  return $class->SUPER::new(Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor->new(@args));
}

sub repository_path{
  my ($self, $rpath) = shift;

  #Use MetaContainer instead
  my $sql;

  if(! defined $self->{repository_path}){
	$sql = "SELECT meta_value from meta where meta_key='repository';";
	($self->{repository_path}) = @{$self->db->dbc->db_handle->selectrow_arrayref($sql)};
  }

  if(defined $rpath){

	if(defined $self->{repository_path} &&
	   ($self->{repository_path} ne $rpath ) ){
	  throw("Supplied repository path does not match meta_value: $rpath ne ".$self->{repository_path});
	}
	elsif(! defined $self->{repository_path}){
	  $sql = "INSERT into meta(meta_key, meta_value) values('repository_path', '${rpath}'";
	  $self->db->dbc->do($sql);
	}
  }

  return $self->{repository_path};
}


#Can take InputSet or InputSubset
#Change this to take sets rather than one?

#todo update this wrt local_url
#what are force_download, data, and skip_beds for
#surely that should be handled in a wrapper script?

#download_input_set_data will need revising


#Remove InputSet support from this to simplify it
#why are we still getting two records for a single InputSubset?
#we are still having to contend with merged input_subsets here
#we really need to fix this?

sub _columns{
  my $self       = shift;
  #unlikely table_name will be 0
  my $table_name = shift || throw('Must defined a table_name argument');
   
  if(! defined $table_name){
    throw('Must defined a table_name argument');  
  }  
  
  return keys %{$tracking_columns{$table_name}};  
}


sub _is_mandatory_column{
  my $self       = shift;
  #unlikely table_name or col will be 0
  my $table_name = shift || throw('Must defined a table_name argument');
  my $col        = shift || throw('Must defined a columns argument');
    
  if(! exists $tracking_columns{$table_name}{$col}){
    throw($col." is not a valid column, must be one of:\t".join(' ', $self->_columns($table_name)));  
  }
  
  return $tracking_columns{$table_name}{$col};
}


sub _is_column{
  my $self       = shift;
  #unlikely table_name or col will be 0
  my $table_name = shift || throw('Must defined a table_name argument');
  my $col        = shift || throw('Must defined a columns argument');
    
  return exists $tracking_columns{$table_name}{$col};  
}


#optionally add in class specific subs like is_embargoed and set_dowload_status?
#As these methods are not simple getter/setters it maybe cleaner not to inject them
#But this leaves some methods available via the object, and some by the adaptor only :/
#This is a case to have separate class specific _inject methods, which call the generic
#_inject methods(for the getter/setters), then inject the additionally inject the
#class specific ones
#One downside to this is error handling, where line numbers will reference anonymous subroutines
#and not the true line number in this module. This is likely fine for the getter/setters, but might
#be problemtic with more complicated methods.

sub _inject_tracking_methods{
  my $self       = shift;
  my $storable   = shift || throw('Must pass a Storable argument');
  my $table_name = $self->get_valid_stored_table_name($storable);   
    
  #We don't need to test for existing methods in the namespace here(CvGV_name_or_bust)
  #as we know exactly what methods are already present
  my @cols = $self->_columns($table_name);
  
  if(! $storable->can($cols[0])){ 
    my $ref = ref($storable);
    no strict 'refs';
  
    for my $col(@cols){
      *{$ref."::${col}"} = 
        sub {
          my $self = shift;
          my $val  = undef;
          
          #to prevent auto vivification
          #which may prevent fetch                   
          if(exists $self->{tracking_info}){
            $val = $self->{tracking_info}{$col};  
          } 
          return $val; 
        };
    }
   
    use strict;
    #avoids being able to work with symbolic reference
    #i.e. Your::Package::Name->$param_method
    #rather than $obj->$param_method  
  }
  
  return;
}


#Todo 
#Handle clashes between existing tracking_info attribute and fetched data

sub fetch_tracking_info{
  my $self       = shift;
  my $storable   = shift || throw('Must pass a Storable argument');
  my $table_name = $self->get_valid_stored_table_name($storable);   
  my $db         = $self->db;

  $self->_inject_tracking_methods($storable);
  my $hashref;

  if(exists $storable->{tracking_info}){
    throw("$table_name Storable(".$storable->dbID.
      ") already has tracking_info attribute. Please store/update before fetching\n");  
  }else{

    my $sql = 'SELECT '.join(', ', ($self->_columns($table_name))).
              " FROM ${table_name}_tracking WHERE ${table_name}_id =".$storable->dbID;
    my $sth = $self->prepare($sql);
    $sth->execute;
    $hashref = $sth->fetchrow_hashref; #Will this be undef if there is no data?
    $storable->{tracking_info} = $hashref;
    $sth->finish;#otherwise we get disconnect warnings
    
    #my %column;
    #pseudo array/hash? 
    #$sth->bind_columns( \( @column{ @{$sth->{NAME_lc} } } ));
    #while( $sth->fetch ){
    #  my $record = {%column}; #deref properly as %columns will be updated
    #  push @tracking_info, {%$record}; #Keep the dbID here
    #  #otherwise there will be now way to identify what the record refers to
    #  my $dbID = delete($record->{input_subset_id}); #Don't need this
    #  push @{$subset_cache{$dbID}->{tracking_info}}, $record;
    #}
  }
  
  return $hashref;
}


#TO DO
#This should also set the local_url?
#Or should we drop this in favour of a generic store/update method? 

sub set_download_status_by_input_subset_id{
  my ($self, $iss_id, $to_null) = @_;
  
  my $date = ($to_null) ? 'NULL' : 'NOW()';
  my $sql = "UPDATE input_subset_tracking SET download_date=${date} WHERE input_subset_id=${iss_id};";
  $self->db->dbc->do($sql);

  return;
}


#can probably remove this, in favour of calling fetch_input_subset_tracking_info
#and then using direct methods?
#can we also add is_embargoed to the injected methods dependant on the presence of an embargo field?

#sub is_InputSubset_downloaded {
#  my ($self, $isset) = @_;
#  $self->fetch_InputSubset_tracking_info($isset);  
#  return (defined $isset->download_date) ? 1 : 0;
#}


#Careful this is american format for some reason!!!

#Can we inject a version this for input_subsets?

sub is_InputSubset_embargoed {
  my ($self, $isset, $yyyy_mm_dd) = @_;
  
  my $now_date = DateTime->today();
  my $rel_date;
  
  if(! defined $yyyy_mm_dd){
    $rel_date = $now_date;
  }
  else{
  
    my ($y, $m, $d) = split(/-/, $yyyy_mm_dd);
  
    if(! (defined $m && defined $d)){ #we have just passed a month number or abbrv
      $m = $y;
      undef $y;
    
      if($m !~ /^[0-9]{1,2}$/){
        my $month_number = get_month_number($m);
   
        if(! $month_number){
          throw("Unable to identify month number from:\t$m\n".
            'Please specify a two digit number (e.g. 01) rf a three letter code (e.g. Jan)');
        }
    
        $m = $month_number;        
      }
    
    
      #we need to figure out whether the month is this year or next year
      $y = $now_date->year;
    
      if($m < $now_date->month){
        $y++;    
      }
    
      $rel_date = DateTime->new(month => $m, year  => $y);
    }
  }
  
  #InputSubsets can have more than one tracking record at the moment
  #as they are currently merged at the InputSubset level rather than the InputSet level
  $self->fetch_InputSubset_tracking_info($isset);
 
  #This was causing issues as DateTime cannot take the default null date of 0000-00-00
  my $avail_date = $isset->availability_date;
    
  if($avail_date eq '0000-00-00'){
    warn "Found $avail_date availability date for InputSubset:\t".$isset->name.
      "\nDefaulting to 0001-01-01\n";
    $avail_date = '0001-01-01';
  }
 
  my ($year, $month, $day) = split(/-/, $avail_date);
  
  
  my $isset_date = DateTime->new(day   => $day,
                                 month => $month,
                                 year  => $year  );
                                                                  
  #Nice operator overloading!
  return ($isset_date > $rel_date) ? 1 : 0;  
}



#todo handle update safely
#INSERT IGNORE based on update flag or always INSERT_IGNORE?
#Will we always know to update?
#Maybe not, but we should know when not to update?
#Will it always be safe to update? There may be old info lurking in the table
#could have an overwrite function which will set the NULLs?
#tricky as this may bork good data.
#Probably need to pull back existing data?
#test return of INSERT IGNORE, if 0E0, then we know it has failed and should do an update
#certainly shouldn't update by default, as this will always silently overwrite existing data
#so we probably need allow_update which will update only those tracking attributes which
#have data, and a purge/overwrite flag, which will also update the NULL attributes to NULL
#have these as flags or params will be safer in case of a flag ordering issue

sub store_tracking_info{
  my $self     = shift;
  my $storable = shift;
  my $info     = shift;
  my $params   = shift;
  #update flag? or use separate methods for dates, urls md5s and ting?   
  
  my $table_name = $self->get_valid_stored_table_name($storable);   
  assert_ref($info, 'HASH', 'Tracking info HASH');
  
  my @cols       = ($table_name.'_id');
  my @values     = ($storable->dbID);
  my @valid_cols = $self->_columns($table_name);

  #Test for unexpect info items  
  foreach my $col(keys %$info){
    
    if(! $self->_is_column($table_name, $col) ){
      #Could also grep this from @valid_cols?
      throw("Found unexpected parameter in tracking info hash:\t".$col.
        "Must be one of:\t@valid_cols");  
    }  
  }
  
  #Test for mandatory info items, and build cols/values
  foreach my $col(@valid_cols){
    
    if($self->_is_mandatory_column($table_name, $col) &&
       ((! exists $info->{$col}) || (! defined $info->{$col}))){
      throw("Mandatory tracking column must be defined:\t$col\n$table_name:\t".$storable->name);     
    }       
    
    if(defined $info->{$col}){
      push @cols,   $col; 
      push @values, $info->{$col};
    }
  }
  
  #Use SQLHelper::execute_update here?
  

  
  my $sql = "INSERT into ${table_name}_tracking(".join(', ', @cols).
    ') values("'.join('", "', @values).'")';
  
  #Although working with DBConnection here provides error handling  
  $self->db->dbc->do($sql);  
  
  #make this generic and pass table_name?
 
  $self->_inject_tracking_methods($storable);
  #could pass $table_name here too, but let the method validate
  $storable->{tracking_info} = $info;  
  
  return $storable;  
}


sub get_valid_stored_table_name{
  my $self     = shift;
  my $storable = shift || throw('Must pass a Storable object argument');
    
  (my $class = ref($storable)) =~ s/.*:://;
  my $table_name = get_table_name_from_class($class);
  
  if(! exists $tracking_columns{$table_name}){
    throw("$class table_name ($table_name) does not correspond to a valid tracking class:\n\t".
      join("\n\t", keys %tracking_columns));  
  }
  
  $self->db->is_stored_and_valid('Bio::EnsEMBL::Funcgen::'.$class, $storable);
  
  return $table_name;  
}


#Move this to DBAdaptor?
#Not BaseAdaptor as this is not a BaseAdaptor (is multi table adaptor)
#but probably shouldn't be DBAdaptor either?

sub get_table_name_from_class{
  my $class = shift || throw('Must provide a class argument');
  return join('_', split_CamelCase($class));
}



1;
