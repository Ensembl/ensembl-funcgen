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
use Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor;
#use Bio::EnsEMBL::Funcgen::DBSQL::BaseAdaptor;#DBI sql_types import
use Bio::EnsEMBL::Utils::Exception         qw( throw warning );
use Bio::EnsEMBL::Funcgen::Utils::EFGUtils qw( get_month_number );

use base qw( Bio::EnsEMBL::DBSQL::BaseAdaptor );

#Currently don't use and Funcgen BaseAdaptor methods in here
#and this would require _true_tables to be defined


#Empty true tables method as required by funcgen BaseAdpator::new
#sub _true_tables{
#  return;  
#}

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

sub fetch_InputSubset_tracking_info{
  my ($self, $set, $force_download, $date, $skip_beds) = @_;

  #can remove this when we update download_inputset_data.pl
  if($force_download || $date || $skip_beds){
    throw('fetch_InputSubset_tracking_info no longer supports thte force_dowload, date or skips_beds arguments');  
  }


  my @sub_sets;

  if(! $set){
	  throw('You need to pass a valid stored InputSet or InputSubset');
  }
  elsif(ref($set) eq 'Bio::EnsEMBL::Funcgen::InputSet'){
	  $self->db->is_stored_and_valid('Bio::EnsEMBL::Funcgen::InputSet', $set);
	  push @sub_sets, @{$set->get_InputSubsets};
	
  }
  elsif(ref($set) eq 'Bio::EnsEMBL::Funcgen::InputSubset'){
	  $self->db->is_stored_and_valid('Bio::EnsEMBL::Funcgen::InputSet', $set);
	  push @sub_sets, $set;
  }


  my (%subset_cache, @tracking_info);

  foreach my $sset(@sub_sets){
  
    if(exists $sset->{tracking_info}){
      push @tracking_info, $sset->{tracking_info};
    }
    else{
      $sset->{tracking_info} = [];
      $subset_cache{$sset->dbID} = $sset;   
    }
  }
  
  if(keys %subset_cache){

    my $sql = 'SELECT iss.input_subset_id, isst.replicate, isst.download_url, isst.downloaded,
               isst.availability_date, isst.md5sum, isst.is_control, isst.not_pooled
               FROM input_subset iss, input_subset_tracking isst
               WHERE iss.input_subset_id=isst.input_subset_id 
               AND iss.input_subset_id IN('.join(',', (keys %subset_cache)) .')';


    #if ($date ne 'IGNORE'){
	  #  $date ||= "NOW()";
	  #  $sql .= " AND ( (isst.availability_date IS NULL) OR (isst.availability_date < '${date}')) ";
    #}

    #if(! $force_download){
	  #  $sql .= ' AND isst.downloaded IS NULL';
    #}


    #if($skip_beds){
	  #  $sql .= ' AND iss.name not like "%.bed%"';
    #}

    #warn $sql;
  
    #change this to return a hash with field key names
  
    my $sth = $self->prepare($sql);
    $sth->execute;
    my %column;
    #pseudo array/hash? 
    $sth->bind_columns( \( @column{ @{$sth->{NAME_lc} } } ));
  
  
  
    while( $sth->fetch ){
      my $record = {%column}; #deref properly as %columns will be updated
      my $dbID = $record->{input_subset_id};
      push @tracking_info, {%$record}; #Keep the dbID here
      #otherwise there will be now way to identify what the record refers to
      delete($record->{input_subset_id}); #Don't need this
      push @{$subset_cache{$dbID}->{tracking_info}}, $record;
    }
  }
  
  #return $self->db->dbc->db_handle->selectall_arrayref($sql);
  
  
  #This is all messed up as this handles InputSet as well as InputSubsets
  #and we want to cache the result and return the records?
  #should restrict this to InputSubsets only
  #where are the callers for this?
  
  return \@tracking_info;
}

sub set_download_status_by_input_subset_id{
  my ($self, $iss_id, $to_null) = @_;
  
  my $date = ($to_null) ? 'NULL' : 'NOW()';
  my $sql = "UPDATE input_subset_tracking SET downloaded=${date} WHERE input_subset_id=${iss_id};";
  $self->db->dbc->do($sql);

  return;
}



#These will both work for InputSets too as this is supported by 
#fetch_InputSubset_tracking_info

sub is_InputSubset_downloaded {
  my ($self, $isset) = @_;
 
  my $downloaded = 1;
 
  foreach my $tr_info(@{$self->fetch_InputSubset_tracking_info($isset)}){
    if(! defined $tr_info->{local_url}){
      $downloaded = 0;  
      last;
    }
  }
  
  return $downloaded;
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
  my $track_records = $self->fetch_InputSubset_tracking_info($isset);
  
  my @embargoed;
  
  #Need to change this as we don't have access to the dbID or the InputSubset name here
  
  foreach my $tr(@$track_records){
 
    #strip off time
    (my $avail_date = $tr->{availability_date}) =~ s/ .*//o;
    my ($year, $month, $day) = split(/-/, $avail_date);
   
   
    my $isset_date = DateTime->new(day   => $day,
                                   month => $month,
                                   year  => $year  );
                                  
                                  
    if($isset_date > $rel_date){ #Nice operator overloading!
      push @embargoed, $isset->name;
    }
  }
   
  return \@embargoed;  
}

1;
