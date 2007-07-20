#!/software/bin/perl -w


####!/opt/local/bin/perl -w


=head1 NAME

ensembl-efg stable_id_mapper.pl
  
=head1 SYNOPSIS

make_nr_probe_fasta.pl [options]

Options:

Mandatory


Optional


=head1 OPTIONS

=over 8

=item B<-name|n>

Mandatory:  Instance name for the data set, this is the directory where the native data files are located

=item B<-format|f>

Mandatory:  The format of the data files e.g. nimblegen

=over 8

=item B<-group|g>

Mandatory:  The name of the experimental group

=over 8

=item B<-data_root>

The root data dir containing native data and pipeline data, default = $ENV{'EFG_DATA'}

=over 8

=item B<-fasta>

Flag to turn on dumping of all probe_features in fasta format for the remapping pipeline

=item B<-norm>

Normalisation method, deafult is the Bioconductor vsn package which performs generalised log ratio transformations

=item B<-species|s>

Species name for the array.

=item B<-debug>

Turns on and defines the verbosity of debugging output, 1-3, default = 0 = off

=over 8

=item B<-help>

Print a brief help message and exits.

=item B<-man>

Prints the manual page and exits.

=back

=head1 DESCRIPTION

B<This program> takes a input redundant probe name fasta file and generates an NR probe dbID fasta file.

=cut


#add @INC stuff here, or leave to .bashrc/.efg?

BEGIN{
  if (! defined $ENV{'EFG_DATA'}) {
	if (-f "~/src/ensembl-functgenomics/scripts/.efg") {
	  system (". ~/src/ensembl-functgenomics/scripts/.efg");
	} else {
	  die ("This script requires the .efg file available from ensembl-functgenomics\n".
		   "Please source it before running this script\n");
	}
  }
}
	

#use Bio::EnsEMBL::Root; #Only used for rearrange see pdocs
#Roll own Root object to handle debug levels, logging, dumps etc.

### MODULES ###
use Getopt::Long;
#use Carp;#For dev only? cluck not exported by default Remove this and implement in Helper
use Pod::Usage;
#POSIX? File stuff
use File::Path;
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Utils::Exception qw( throw warning );
use Bio::EnsEMBL::Funcgen::Utils::EFGUtils qw (open_file run_system_cmd backup_file);
use Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Funcgen::DBSQL::AnnotatedFeatureAdaptor;
use strict;

$| = 1;							#autoflush
my ($dbname, $help, $man);
my ($ndbname, $from_file);
my $reg = "Bio::EnsEMBL::Registry";

my $user  = 'ensro';
my $nuser = $user;
my $pass = '';
my $npass = $pass;
my $host = 'ens-livemirror';
my $nhost = 'ens-genomics1';
my $port = '3306';
my $nport = $port;
my $dump = 1;
my $update = 0;
my $expand = 0;

#Definitely need some sort of Defs modules for each array?

$main::_debug_level = 0;
$main::_tee = 0;


#Use some sort of DBDefs for now, but need  to integrate with Register, and have put SQL into (E)FGAdaptor?
#Use ArrayDefs.pm module for some of these, class, vendor, format?
#ArrayDefs would also contain paths to data and vendor specific parse methods?

GetOptions (
			"odbpass=s"     => \$pass,
			"odbport=s"     => \$port,
			"ohost=s"     => \$host,
			"ouser=s"     => \$user,
			"odbname=s"   => \$dbname,
			"ndbpass=s"     => \$npass,
			"ndbport=s"     => \$nport,
			"nhost=s"     => \$nhost,
			"nuser=s"     => \$nuser,
			"ndbname=s"   => \$ndbname,
			"expand|e=s"  => \$expand,
			"dump"        => \$dump,
			"update"      => \$update,
			"from_file|f=s" =>\$from_file,
			"help|?"       => \$help,
			"man|m"        => \$man,
		   );


pod2usage(1) if $help;
pod2usage(-exitstatus => 0, -verbose => 2) if $man;


#check mandatory params here

my $ndb = Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor(
												  host   => $nhost,
												  user   => $nuser,
												  pass   => $npass,
												  port   => $nport,
												  dbname => $ndbname,
												 );

#Test connection, eval this?
$ndb->dbc();

my ($next_stable_id, %old_stable_ids, %new_stable_ids, %retired_stable_ids, %adaptors);
my $iterations = 0;
my $mappings = 0;

#we also need a hash with keys of new dbID
# or can we mapd the values of the old_stable_ids{'new'} and 
#jsut duplcaite the hashes for now for simplicity


if(! $from_file){

  my $odb = Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor(
												  host   => $ohost,
												  user   => $ouser,
												  pass   => $opass,
												  port   => $oport,
												  dbname => $odbname,
												 );

  #Test connection, eval this?
  $odb->dbc();

  $obj_cache{'OLD'}{'AF_ADAPTOR'}    = $odb->get_AnnotatedFeatureAdaptor();
  $obj_cache{'OLD'}{'SLICE_ADAPTOR'} = $odb->get_SliceAdaptor();
  $obj_cache{'OLD'}{'RF_TYPE'}       = $odb->get_FeatureTypeAdaptor->fetch_by_name('RegulatoryFeatures');
  $obj_cache{'NEW'}{'SLICE_ADAPTOR'} = $ndb->get_SliceAdaptor();
  $obj_cache{'NEW'}{'AF_ADAPTOR'}    = $ndb->get_AnnotatedFeatureAdaptor();
  $obj_cache{'NEW'}{'RF_TYPE'}       = $odb->get_FeatureTypeAdaptor->fetch_by_name('RegulatoryFeatures');

  $next_stable_id = $db->dbc->db_handle->do('SELECT display_label from annotated_feature af, feature_set fs where fs.feature_type_id='$oftype->dbID().' and fs.feature_set_id=af.feature_set_id order by display_label desc limit 1');
  $next_stable_id =~ s/\:.*//;
  $next_stable_id =~ s/ENSR[A-Z]*0*//;
  $next_stable_id ++;
  
  #we're going to have to have some level of redundancy here
  #as we can't depend on recording 
  
  #stable_id => {(
  #                old => [X, ..],
  #                new => [Y, ..],
  #                stable_ids   => [N, ..],
  #              )}
  #still not ideal, as we want all stable ids to be stored in once place

  #use redundant but probably faster three hashes
 



  #we need to do some sort of localised all vs all overlap
  #i.e. do the associated forward and reverse.

  #mapping logic shdould be very simple to start with as we're changing the data quite a lot
  #and we don't want to create a huge amount of useless messy data
  #what is the concept of a RegulatoryFeature?
  #It can cross cell/tissue types but have differing modifications, hence states
  #but maintain same stable_id
  #This makes mapping between assembly variable region extremely hard
  #we need to take account of multi cell/tissue types in some way.
  #we need to factor this into the schema



  foreach my $reg_feat(@{$opf_a->fetch_all_by_logic_name('RegulatoryRegion')}){
	my $from_db = 'OLD';
	#my $features_slices = []; #this has to be ref due to implementation in while
	my $feature = $reg_feat;
	my $mapping_info = {}; #this has to be ref due to implementation in while

	#o --  --  --------- -- -------------
	#n -----------   -----------   --------
	#how do we hande this?
	#middle old regf should be split
	#which stable id should be inherited for each new regf?
	#always take largest overlapping old regf stable_id
	#else if already taken just take next


	#what if next has even bigger overlap with another???
	#o--- ---------- --------------------- -----
	#n  ----- ---------------- --------------- ----
	#                                         ^

	#working from left to right, we're constantly shifting which feature inherits the stable_id
	
	#This is highly unlikely, but how are we going to handle it?
	#This is not really a problem which any other stable id mapping comes across
	#we should just log these chains and map manually?
	#cut off after 3 iteration i.e. forward, backward, forward
	#if we find a larger over lap on the last forward, then keep going until we find a smaller or no overlap
	#don't try and solve the id mappings, just log the details
	#or arguably, just create new stable IDs for all
	#or do full path through overlaps moving stable_ids between bins as we find new 
	#are more suited to new feature
	
	#basically we just need to store an array of old ids
	#until we reach the last one. Also store an array of the new features
	#Then we jus inherited from the count element, i.e. +1
	#then we need to handle end which may still have an overlap with another feature
	#we clean the arrays and pass back the last old feature ID so we don't include this in the next overlap anal
	#also have to consider reverse and middle feature
	
	#so we do then same thing backwards and if middle features does not donate it's stable ID, then it is recorded as a split.
	#if no overlap then it is recorded as a merge
	
	#NoNoNo!  We don't need to go backwards as we will have already seen any 5' old features and hence this stable feature 
	#will have been dealt with by the chaining algorithm

	#we only have to consider 5' new features



	#however we need to keep iterating until we com to the end of an overlapping sequence
	#or we find an overlap which is less than the previous
	#i.e. we assign the old stable_id with the largest overlap
	#so long as this doesn't overlap an adjacent feature more
	#we still then need to deal with the last old conjugated feature/stable_id
	#if we don't we could end up doing the same comparison but backwards
	#unles we ignore new dbIDs which have already been assigned?

	#do we want to expand on the new features too?



	#do we need a sub for the forward mapping to allow recursive mapping
	#during feature merging?


	#basically we want to keep doing the same comparison
	#but populating different hashes dependent on wether 
	#we're doing the forward or reverse mapping
	#until we don't return any features which over lap the end
	#This also needs to be bidrectional, ensuring we get the same
	#result is we start at either end of a contiguous chain of old and new features

	#so we need a tmp hash to store the current default stable_id
	#and the length it's overlap
	#all other stable IDs can be recorded in their appropriate hashes as we go along
	#we need to return pairs of current feature and overlappping feature, 
	#so we can exlude the current feature from the next overlap analysis

	#we need to supply list ref of current feature, last feature pairs
	#or just the start feature

	
	my $stable_id;
	($stable_id = $reg_feat->display_label()) =~ s/\:.*//;

	#$stable_db_ids{$stable_id} = {(
	#							   o => $reg_feat->dbID(),
	#							  )};

	#skip stable id if we've already seen it
	next if exists $old_stable_ids{$stable_id};
	

	#push @{$mapping_info{'stable_id'}} = $stable_id;
	#$mapping_info{'extended'}  = 0;
	#$mapping_info{'iteration'} = 0;
	#should also have a mappings element which is an array of all the full stable_id 
	#mappings  we've walked through, this will allows us to track back and change stable_ids if required
	#we also need to track feature slice pairs in each mapping?


	#when we're doing this by iteration the start point will be a new reg feature or could be slice
	
	#so maybe here is the point to call the first sub
	#$feature_slices = [$reg_feat->feature_Slice()];
	

	while(defined $feature){

	  #we don't need to sub this part do we?
	  #just the id assignment and id hash population

	  ($mapping_info, $from_db, $feature, $source_dbID) = &do_feature_map($mapping_info, $from_db, $feature, $source_dbID);

	}

	#do we need to deal with last here?  I don't think so
  

	#Now for each new feature which hasn't already been mapped, generate a new stable id



}
else{
  #Parse previously dumped file to avoid having to do mapping again
  
  
}


#THIS ASSUMES REGULATORY FEATURES NEVER OVERLAP WITHIN THE SAME BUILD
#This will currently give a direction bias to the inheritance
#if overlap are equal in % overlap and number of constituents, default to 5'
#then consider reverse overlap

#o  --- -- -    --------
#n   --------------  ----
#                   ^

#nfeat with more old features inherits last ofeat stable id

#o  --- ----
#n    ---

#nfeat inherits first ofeat stable ID as reverse overlap is greater

#The following will always give a direction bias, as we have no way to discriminate, so we plump for 5'
#o --- ---
#n   ---

sub do_feature_map{
  my ($mapping_info, $from_db, $feature, $source_dbID) = @_;

  my $to_db = ($from_db eq 'OLD') ? 'NEW' : 'OLD';
  my $slice = $feature->slice();

  #the source_rf_dbID is the dbID of the previous source feature
  #to ensure we don't redo the overlap
  #we use this to figure out whether this feat or the last inherits the stable ID
  #dependent on which direction we're going in
  
  #n-- ------
  #o -----  ------
  
  #In the above example the current feature is the 2nd nfeat
  #hence we must compare the overlap of the last source(which is the 1st ofeat i.e. the source_rf_dbID)
  #with the next ofeat.  Here the previous/source has the greater overlap so this nfeat inherits the previous 
  #ofeat stable ID
      
  #disabled default extend for now for simplicity
  # we should only do this if we have no mapping
  #if($expand && $mapping_info->{'iteration'} == 0){
  #	$mapping_info->{'extended'} = 1;
  #	$slice = $slice->expand($expand, $expand);
  #  }
  
  
  my $next_slice = $obj_cache{$to_db}{'SLICE_ADAPTOR'}->fetch_by_region('toplevel', $slice->seq_region_name(), $slice->start(), $slice->end());
  
  if(! $next_slice){
	#I don't think this should ever happen, can we remove?
	print $reg_feat->display_label()." is on a slice which is not present in the $to_db DB:\t".
	  $slice->seq_region_name().' '.$slice->start().' '.$slice->end."\n";
	next;
  }
  
  
  my @nreg_feats = @{$obj_cache{$to_db}{'AF_ADAPTOR'}->fetch_by_Slice_FeatureType($next_slice, 
																				  $obj_cache{$to_db}{'RF_TYPE'})};
  
  #disabled default extend for now for simplicity
  #only retry with default extend if we haven't extended already and we don't already have a mapping
  #if(! @nreg_feats && ! $extend && $mapping_info->{'iteration'} == 0){
  #  $next_slice = $next_slice->expand(200, 200);#Just enough to catch the next nucleosome
  #  $mapping_info->{'extended'} = 1;
  
  #  @nreg_feats =  @{$obj_cache{$to_db}{'AF_ADAPTOR'}->fetch_by_Slice_FeatureType($next_slice, 
  #$obj_cache{$to_db}{'RF_TYPE'})};
  
  #}
  
  
  if(! @nreg_feats){
	  
	#This is the first old to new mapping which has failed
	if(! defined $source_dbID){
	  print "No $to_db mapping possible for:\t".$mapping_info->{'stable_id'}."\n";

	  warn "we need to record retirement of stable_id here";

	  next;#there should only be one pair in this loop, so will exit loop
	}
  }
  else{
	#ignore source_dbID for comparison(include in overlap ids for logging)
	#store others associated with current dbID
	my ($overlap, $current_child, @nfeats, @overlap_ids);
	my $end_overlap = 0;
	
	foreach my $nfeat(@nreg_feats){
	  #if we're in the middle/end of a chain, the first will always be the source_dbID/dbID
	  push @overlap_ids, ($to_db eq 'NEW') ? $_->dbID() : $_->stable_id();
	  
	  if(defined $source_dbID && $nfeat == $source_dbID){
		#here we need to consider wether this over lap is greater than the last
		#o --     --   ----
		#n  ------------- ---  
		#                ^
		
		#comparing n to o
		#if it is not we can add the current nfeat to the mapping hash as a split
		#assign the stable id to the last nfeat
		#clean the mapping hash and carry on
		
				
		#do we next here?
		#yes, but we need to make sure we have recorded the split properly
		#aove example would spawn a new stable id for the last new feature
 
		#??????????????????????????????  IS THIS RIGHT
		#It's better for the assign_stable_ids alg
		#if this is the last in the chain i.e. not other nreg_feats
		#but do we not need this for overlap_ids?
		#can we no just populate this as we go along?
		#will this be enough to record all the relationships going on, merge, split etc...

		next;
		
	  }
	  
	  push @nfeats, $nfeat;
	  
	  if($nfest->start() < $slice->start()){
		#This can only be a 5' new feature
		#as we will have dealt with all 5' old features
		#o  -----
		#n --
		
		if ($to_db eq 'OLD'){
		  throw("Found 5' feature from $to_db, should have already dealt with this");
		  
		  #not quite true as we could be extending
		  #n--    -----
		  #o --
		  #we should sill never see this as the old feature will have been extended first and captured this
		  #new feature, hence we will have skipped it
		  
		}
		elsif(defined $source_dbID){
		  #found 5' new feature which wasn't the source feature we have skipped
		  #this should not be possible
		  throw("Found 5' feature from $to_db which does not match the slice source feature");
		}
		
		throw("Found more than one 5' feature") if(scalar(@nfeats) >1);
		$current_child = 0;
		
		#we need to consider extend here!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		
		if($nfeat->end() >= $slice->end()){
		  $overlap = $slice->length();
		}else{
		  $overlap = $slice->start() - $nfeat->end();
		}
		
	  } elsif ($nfeat->end() <= $slice->end()) {
		#all fully contained features;
		
		my $tmp_length = $nfeat->length();
		
		if ($tmp_length > $overlap) {
		  $overlap = $tmp_length;
		  $current_child = $#nfeats;
		}
	  } else {					#found end overlap
		  
		my $tmp_length = $nfeat->start() - $slice->end();
		$end_overlap = 1;
		  
		if ($tmp_length > $overlap) {
		  $overlap = $tmp_length;
		  $current_child = $#nfeats;
		}
	  }
	}
	  
	#so here we have to trap the case where we have a to nothing(only parent) mapping, i.e at the end of a chain
	#basically we don't want this transition as it doesn't really exist
	#o-1-  --3--
	#n  --2--

	#1>2
	#2>3(2)
	#3>nothing(2)
	
	#so basically if we don't have a current child, we have a nothing transition
	#so go straight to assigning stable_ids
	#else deal with building transition

	if(! $current_child){
	  &assign_stable_ids($mapping_info);		
	  undef $feature;#KEY!! while will exit and next reg feature retrieved!
	  $mapping_info = {};#just to be clean, altho this doesn't really matter
	  undef $source_dbID;#just to be clean, altho this doesn't really matter


	  #we can print here !!
	}
	else{
  
	  #now build transition hash and populate mapping info hash appropriately
	  #calling assign_stable_ids if we have found the end of a chain
	  
	  #3 situations
	  #i) Simple no end overlap
	  #b) End overlap, but not largest
	  #finally) End overlap and is largest of set
	  
	  
	  #mapping info
	  #overlap - required for reverse overlap comparison
	  #here we could add % last overlap also?
	  #array of hashes where each hash contitutes a transition?
	  #every other element would eventually boil down to one stable ID

	  my %transition = (
						source         => $from_db,
						overlaps       => \@overlap_ids;#used when overlap == to next
						current_child  => $current_child
					   );
	  
	  if ($to_db eq 'OLD') {
		push @{$mapping_info->{'stable_ids'}}, $tmp[$current_child]->stable_id();
		#we need to record dbID of feat we're mapping to
		push @{$mapping_info->{'dbIDs'}}, $feature->dbID();
	  } 
	  else {
		#do we need to record the dbIDs here too?
		push @{$mapping_info->{'dbIDs'}},  $tmp[$current_child]->dbID();
		push @{$mapping_info->{'stable_ids'}}, $feature->stable_id();
	  }
	
	  #reset source dbID
	  $source_dbID = $feature->dbID();
	  	  
	  if ($current_child == $#tmp) { 
		#end overlap is largest, need to update mapping info
		#o----  -----------  ------ ----------
		#n  -------  - -------------------
		
		$transition{'overlap'}   = $overlap;
		$transition{'overlap_factor'} = ($overlap/$feature->length());#used when overlaps ==
		$feature = $tmp[$#tmp];
		push @{$mapping_info->{'transitions'}}, \%transition;	
	  }
	  elsif($end_overlap) {	
		#end over lap, but not largest, do stable_id mapping and re-set mapping info?
		#we still need to consider last overlap even if last stable_id was not used by previosu feature
		#o ---- ----
		#n-------- ----		
		#here the 2nd stable id is not used for the first nfeat as there is a bigger overlap
		#but it overlaps the 2nd nfeat even less
		#this should result in a split, giving rise to an entirely new stable ID?
				
		$transition{'overlap'}   = $overlap;
		push @{$mapping_info->{'transitions'}}, \%tmp;	
		&assign_stable_ids($mapping_info);
		
		$mapping_info = {};
		$feature = $tmp[$#tmp];	#set the last feature
	  } 
	  else {					#no end overlap, do stable id mapping, clean mapping info
		push @{$mapping_info->{'transitions'}}, \%tmp;	
		&assign_stable_ids($mapping_info);
		
		undef $feature;#KEY!! while will exit and next reg feature retrieved!
		$mapping_info = {};#just to be clean, altho this doesn't really matter
		undef $source_dbID;#just to be clean, altho this doesn't really matter
		

		#we can print here!!
	  }
	}
  }


  #periodic printing would be dependent on having a tmp cache of currently 'active' stable_ids and dbids
  #as we can assign ids between split mappings, but split/merge info will not be complete until we finished the
  #next mapping. tricky

  
  #finally return the mapping info hash and switch the db
  
  return ($mapping_info, $to_db, $feature, $source_dbID);
	
}
 


sub assign_stable_ids{
  my $mapping_info = shift;
  

  #note 'nothing' transitions e.g. 5b>nothing(4b) below do not get passed
  #so we do not have to consider them

  #o-1--  ----3--------  --5a-- ------5b--
  #n  --2----  -4a ----4b-------------  xxxxx
  
  #xxx is possible further overlap, but we would split into two mappings including the xx in both,
  #and just pass 5b to the new mapings to ensure we don't migrate backwards

  #This would result in 5 transitions
  #1  > 2
  #2  > 3 (&1)
  #3  > 4a & 4b (&2)
  #4b > 5a & 5b (&3)
  #5b > xxxx which is a smaller overlap and therefore split(&4b)

   
  #we work backwards to define the mapping 
  #The results should be as follows
   
  #Transition 4b > 5a & 5b
  #gives
  #5b & 5a > 4b with 5b as the stable_id
  #due to larger overlap

  #we also need to record merge of 5a
  #but are we not missing the source_dbID of 3 in the merge?
  #make sure this is in the constituents
  #how are we going to record the split of 3 >2, 4a and 4b?!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  #so we just need to add this to the relevant hash when we get to the next transition
  #this problem could cross mapping_info boundaries when we reset due to a non-inheriting overlap!!!
  #this is okay?
  
  #Transition 3  > 4a & 4b
  #gives
  #so now we're working in reverse
  #we know the very last feature has inherited from it's repective position
  #so we inherit from the repective position too e.g. 2
  #and because we're going in the opposite direction we're projecting this st
  #we then has to check the dbID of the other constituents
  
  #what does 4a get?
  #we need to make sure we assign new stable_ids for all constituent dbIDs which are not primary overlap
  #4a new stable id
  

  #for all the other transition we know the biggest overlap is the conjugating overlap
  #at either end depending on inheritance of the last feature(current_child)

  #if we don't have an overlap in the very last transition then we know the biggest
  #is the current child i.e. the dbID if we're going to NEW
  #or the stable_id if we're going to OLD
  #but is this the previous conjugating overlap or something in the middle?
  #so basically we need to check whether current_chld for the last transition is >0
  #if so then we use the repesctive element of the stable_id array else
  #we inherit from the previous overlap, hence we inherit from respective position -1;
  #this is only true for stable_id and dbID arrays of same length!

  #we can have stable_id and dbID arrays off different lengths!
  #this is a function of whether we have a 5' new overlap and a 3' new overlap
  #so we can figure this out by just looking at the size of the stable_id and dbID arrays

  #if we have neither then the sizes from the above example are stable 3, dbID 2
  #then inheritance (OLD as we're working backwards) from stable_id array
  #is transition pos/elem + num seen OLD transitions - 1
  #for 4b>5b(5a & 3), 4b would inherit from    3 + 0 - 1 = 2 = last stable_id (5b)!!! CORRECT!!!
  #for 2>3, 2 would inherit from               1 + 1 - 1 = 1 = 2nd stable id(3) CORRECT!!! 
  #-2 if the last transition current_child == 0;

  #projection rules
  #transition pos/elem + num seen NEW transition - difference between stable_id and dbID lengths - 1
  #3>4b(2) would project to 2 + 0 - 1 - 1 = 0 = 1st dbID
  #1>2     would project to 0 + 1 - 1 - 1 = -1 = would be previous dbID if it existed
  #but doesn't project...so CORRECT !!!
  #remove -1 if last transition current_child == 0


  #so what if they are equal? (reverse will not give same number of transitions...see below)

  #o-1- --3---- ---5-------
  #n  -2--  --4----- --6-----

  #stable_ids == 3  dbIDs == 3
  #projection

  #5>6(4) 4%2 = 2 = last dbID...YAY!
  #3>4(2) 2%2 = 1 = 2nd dbID  YAY!
  #1>2    0%2 ? = 0 = 1st dbID YAY!!  Careful no to divide zero!!

  #inheritance == same rule?
  #we need to use modulus here to round down
  #4>3(5) 3%2 = 1 = 2nd stable_id  CORRECT BOY WONDER!!
  #2>1(3) 1%2 = 0 = 1st stable_id  KABLAMMMM!!


  #or..NO this is not equal!!

  #o  -1--  --3----- ---5----
  #n2a- --2b--- ---4-------

  #stable_ids = 3 and dbIDs = 2
  #remember at present 2a is not stored in the dbID array, this is probably wrong as it is a potential mapping

  #1>2a & 2b
  #2b>3 &1
  #3>4 & 2b
  #4> 5 & 3
 
  #this is same as first example

  #inheritance transition pos/elem + num seen OLD transitions - 1
  #4>5(3)   3 + 0 - 1 = 2 = 3rd stable_id...CORRECT!!
  #2b>3(1)  1 + 1 - 1 = 1 = 2nd stable_id...CORRECT!!
  #shift is same

  #projection transition pos/elem + num seen NEW transition - difference between stable_id and dbID lengths - 1
  #3>4(2b)  2 + 0 - 1 -1 = 0 = 1st dbID CORRECT!!
  #1>2b(2a) 0 + 1 - 1 -1 = -1 = ??? Is correct but element doesn't exist dbID

  #we need to add 2a to the dbIDs array, but then rules will change no?
  #using rules for equal sized arrays gives transition elem%2
  #inheritance 4>5(3)   3%2 = 1 = WRONG!!!
  
  #so we either need to change the rules to take account of which direction last transition is
  #or just store 2a dbID in the mapping hash and if we get -1 then we check it exists and use it
  #else we retire.
  

  #HERE!! WE NEED TO IMPLEMENT THIS BELOW!!!!!!!!!!!!!!!!!!!!!!!!!!
  #add 2a to mapping_info as ???? only if it is not the biggest overlap
  #as it would be added as a transiotion otherwise???
  #need to test this set up
  #o   -----    -----
  #n-----  --------

  #and test  if stable_ids < dbIDs


  #as we work backwards, we're dealing with every mapping twice due to nature of transitions
  #


}



1;
