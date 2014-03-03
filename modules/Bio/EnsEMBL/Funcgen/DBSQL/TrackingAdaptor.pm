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
use Bio::EnsEMBL::Funcgen::Utils::EFGUtils qw( get_month_number );

use base qw( Bio::EnsEMBL::DBSQL::BaseAdaptor );
#Currently don't use and Funcgen BaseAdaptor methods in here
#and this would require _true_tables to be defined


my %mandatory_columns = (#'input_subset_id'   => 0,
                         availability_date => 1,
                         download_url      => 1,
                         downloaded        => 0,
                         local_url         => 0,
                         md5sum            => 0,
                         notes             => 0);
                        

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
  return keys %mandatory_columns;  
}

sub _is_mandatory_column{
  my $self = shift;
  my $col  = shift;
  
  if(! defined $col){
    throw('Must defined a columns argument');  
  }  
    
  if(! exists $mandatory_columns{$col}){
    throw($col." is not a valid column, must be one of:\t".join(' ', $self->_columns));  
  }
  
  return $mandatory_columns{$col};
}

sub _is_column{
  my $self = shift;
  my $col  = shift;
  
  if(! defined $col){
    throw('Must defined a columns argument');  
  }  
    
  return exists $mandatory_columns{$col};  
}

#Change this to set the tracking info in the subsets!
#and inject the relevant methods based on the columns

sub _inject_input_subset_tracking_methods{
  my $self = shift;
  my $iss  = shift;
  
  assert_ref($iss, 'Bio::EnsEMBL::Funcgen::InputSubset');  
  
  #We don't need CvGV_name_or_bust here as we know 
  #exactly what method are already present
  my @cols = $self->_columns;
  
  if(! $iss->can($cols[0])){ 
    no strict 'refs';
  
    for my $col(@cols){
      *{ref($iss)."::${col}"} = 
        sub {
          my $self = shift;
          my $val = undef;
          
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



#remove InputSet support from this!
#then we don't have to select input_subset_id
#although this also removes the ability to query >1 record at the same time
#this is fine for the tracking adaptor

sub fetch_InputSubset_tracking_info{
  my ($self, $set, $force_download, $date, $skip_beds) = @_;

  assert_ref($set, 'Bio::EnsEMBL::Funcgen::InputSubset');

  #can remove this when we update download_inputset_data.pl
  if($force_download || $date || $skip_beds){
    throw('fetch_InputSubset_tracking_info no longer supports thte force_dowload, date or skips_beds arguments');  
  }

  my $db = $self->db;

  #if(! $set){
	#  throw('You need to pass a valid stored InputSet or InputSubset');
  #}
  #elsif(check_ref($set, 'Bio::EnsEMBL::Funcgen::InputSet')){
  #  if($set->is_stored($db)){
  #	  push @sub_sets, @{$set->get_InputSubsets};
  #  }
  #  else{
  #    throw("InputSet is not stored in this DB:\t".$set->name);  
  #  }
  #}
  #elsif(check_ref($set, 'Bio::EnsEMBL::Funcgen::InputSubset')){
	#  if($set->is_stored($db)){
	#    push @sub_sets, $set;
  #  }
  #  else{
  #    throw("InputSubset is not stored in this DB:\t".$set->name);    
  #  }
  #}
  #else{
  #  throw("Set argument is not an InputSubset:\t$set");  
  #}


  $self->_inject_input_subset_tracking_methods($set);
  my $hashref;

  if(! exists $set->{tracking_info}){

    my $sql = 'SELECT '.join(', ', ($self->_columns)).
              ' FROM input_subset_tracking WHERE input_subset_id ='.$set->dbID;

    #warn $sql;
    #if ($date ne 'IGNORE'){
	  #  $date ||= "NOW()";
	  #  $sql .= " AND ( (isst.availability_date IS NULL) OR (isst.availability_date < '${date}')) ";
    #}

    #if(! $force_download){
	  #  $sql .= ' AND isst.downloaded IS NULL';
    #}

 
    my $sth = $self->prepare($sql);
    $sth->execute;
    $hashref = $sth->fetchrow_hashref; #Will this be undef if there is no data?
    $set->{tracking_info} = $hashref;
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

sub set_download_status_by_input_subset_id{
  my ($self, $iss_id, $to_null) = @_;
  
  my $date = ($to_null) ? 'NULL' : 'NOW()';
  my $sql = "UPDATE input_subset_tracking SET downloaded=${date} WHERE input_subset_id=${iss_id};";
  $self->db->dbc->do($sql);

  return;
}



#shouldn't this use the downloaded date?
sub is_InputSubset_downloaded {
  my ($self, $isset) = @_;
  $self->fetch_InputSubset_tracking_info($isset);  
  return (defined $isset->local_url) ? 1 : 0;
}


#Careful this is american format for some reason!!!

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
 
  my ($year, $month, $day) = split(/-/, $isset->availability_date);
  my $isset_date = DateTime->new(day   => $day,
                                 month => $month,
                                 year  => $year  );
                                                                  
  #Nice operator overloading!
  return ($isset_date > $rel_date) ? 1 : 0;  
}


sub store_input_subset_tracking_info{
  my ($self, $iss, $info) = @_;
  #update flag? or use separate methods for dates, urls md5s and ting?   
     
  $self->db->is_stored_and_valid('Bio::EnsEMBL::Funcgen::InputSubset', $iss);
  assert_ref($info, 'HASH', 'InputSubset tracking info HASH');
  
  my @cols       = ('input_subset_id');
  my @values     = ($iss->dbID);
  my @valid_cols = $self->_columns;

  #Test for unexpect info items  
  foreach my $col(keys %$info){
    
    if(! $self->_is_column($col) ){
      throw("Found unexpected parameter in tracking info hash:\t".$col.
        "Must be one of:\t@valid_cols");  
    }  
  }
  
  #Test for mandatory info items, and build cols/values
  foreach my $col(@valid_cols){
    
    if($self->_is_mandatory_column($col) &&
       ((! exists $info->{$col}) || (! defined $info->{$col}))){
      throw("Mandatory tracking column must be defined:\t$col\nInputSubset:\t".$iss->name);     
    }       
    
    if(defined $info->{$col}){
      push @cols,   $col; 
      push @values, $info->{$col};
    }
  }
  
  #Use SQLHelper::execute_update here?
  my $sql = 'INSERT into input_subset_tracking('.join(', ', @cols).
    ') values("'.join('", "', @values).'")';
  
  #Although working with DBConnection here provides error handling  
  $self->db->dbc->do($sql);  
  
  $self->_inject_input_subset_tracking_methods($iss);
  $iss->{tracking_info} = $info;  
  
  return $iss;  
}

1;
