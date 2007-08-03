#!/software/bin/perl -w


####!/opt/local/bin/perl -w


=head1 NAME

ensembl-efg stable_id_mapper.pl
  
=head1 SYNOPSIS

stable_id_mapper.pl [options]

Options:

Mandatory


Optional


=head1 OPTIONS

=over 8

=item B<-name|n>

Mandatory:  Instance name for the data set, this is the directory where the native data files are located


=item B<-debug>

Turns on and defines the verbosity of debugging output, 1-3, default = 0 = off

=over 8

=item B<-help>

Print a brief help message and exits.

=item B<-man>

Prints the manual page and exits.

=back

=head1 DESCRIPTION

B<This program> compares an existing set of Regulatory Features against a new set by performing
overlap analysis to project stable IDs onto the most appropriate new feature. As old and new features 
can overlap, there is always the possibilty of finding a better match than the first encountered.
This is overcome by building chains of overlapping features, taking the 5' most old feature and retrieving 
each new feature on this slice, then reciprocated using the new feature until we exhaust all overlapping 
old and new features in a given chain. Once we have built a chain and stored the overlap qualities of each 
transition from old to new or vice versa, then we can work backwards through the chain applying rules to 
each transition to denote which new feature inherits which stable ID, and which stable IDs are merged, 
split or retired.  Once this has been done any remaining non-overlapping new features are assigned a new 
stable ID.

Rule hierarchy for projecting or inheritance of stable IDs, given that the previous fails to resolve the mapping,
 are as follows:

1. Maximum overlap

old     --a--
new   -b-  --c-

c would inherit from a.
b would be given a new stable ID, but the split of a would be recorded as a > a(c) & b

2. Maximum coverage

old     --a--
new   -b-- --c----

b would inherit from a. Split behaviour is recorded as above.

3. Least displaced

old  ------a----
new    -b-  -c-

b would inherit from a as it is most central to a.

4. Maximum parents

old -x- ---a---
new   -b--   --c-

b would inherit from a as it has another parent regulatory feature.
x would get retired, with the split being recorded.
c would get a new stable id, again with the split of a being recorded

The above examples are just some very simplistic representations and there is obviously the possiblity to have many
combinations of overlap behaviour.  Due to the uncertainty of inheritance of linking overlaps in a chain, it is 
sometimes necessary to traverse 1 or two links before we can assign stable IDs and know the full split/merge behaviour.


????  In the case of a truely unresovled mapping do we take the 5', or do we split and retire ????

Some more examples:

old -x- ---a---
new   -b-     -c-

x>b due to coverage
a>c as b has already had an ID assigned

In reality this assignment is done in reverse, as when considering a, we don't know what might map to be, so we have 
to move the next transition before we discover x is a better match for b, then assign a to c.

???? Are we handling the assignment of a to c properly ???????????????????


example of normal first to previous bahviour here

This can get extremly complicated:

old --a--- c- --c1---- e- e1
new	-b- --b1-----  --d------

Here we have walked through a chain of ever extending or at least stable length linking overlaps.
The mapping for this chain would be as follows:

Walking forwards:

a > b b1
b is marked as the 'child' due to coverage.
b1 is an unknown so we have to wait until we have done the reverse walk to record the split of a properly

b1 > c c1
c1 is marked as 'child' (or more appropriately parent in the case of a new to old transition), however we don't
know whether c1 will have a better match, so we still need to move to the next transition before resovling this mapping.
c is marked as 2nd best or 'previous_child' due to coverage when compared with a, in case c1 is used for something else.

c1 > b1 d
waa waa oops, cannot resolve this yet as we don't know the number over coverage of parents of d, b1 get's marked as child,
but a flag is set to mark d as a twin to b1 in the next transition. This will be used in the reverse walk to make a 
comparison between the parents of d and b1.

?? should we not mark e as 2nd best so if b1 would inherit from c1, d would inherit from e??????????????????????????????????????????????


d > c1 e e1
c1 get marked as child(parent) due to overlap


Now we walk backwards through the chain to resolve the mappings

d > c1 e e1
Natural choice would be c1 due to the maximum overlap, but we know that d is a twin b1 with respect to c1.  So now we 
consider the b1 > a c c1 mapping against the current transition.  We see each has an equal number of parents but e1 has
better coverage than a, so we assign c1 to d and record the merge of e and e1 into d(c1). Potentially, if this were a truly 
symmetrical overlap neighbourhood, this would still be unresolved.  In this case we would split and retire c1 and take the
next best i.e. c > b1 and e > d.

??? We need to generate a sub to do this ??????

c1 > b1 d
We already know c1 has been assign to d, so we move on.  If this were not the case then we'd assign it to the recorded child
i.e. b1

b1 > a c c1
The child for this is c1, but we already know this has been taken so we use the 2nd best which is c

a > b b1
The child for this is b, we know this has not yet been assigned to, and we know there are no more transitions, so we assign
a to b.



things we need to implement on the back of this
recording of 2nd best mappings only immediately after a twin?
a sub loop to handle the bidrectional walking given a twin



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
my $out_dir = '.';


$main::_debug_level = 0;
$main::_tee = 0;


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
			"outdir=s"    => \$out_dir,
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

my ($next_stable_id, $mapping_file, $old_id_file, $new_id_file, %old_stable_ids, %new_stable_ids, %retired_stable_ids, %adaptors);
my $mappings = 0;

if($dump){
  #this contains the new dbID, it's old or new stable id, plus and merged/split stable IDs
  #new dbID/stable_id as focus + parents
  $mapping_file  = open_file($out_dir.'/db_id_stable_id.mappings');

  #this contains retired ids and the their child ids if any
  #old stable new/old stable id as focus + parents
  $old_id_file = open_file($out_dir.'/retired_stable_id.mappings');

  #old_id   new_id  parents
  #old_id1  old_id1 old_id2                 #merge
  #old_id2  old_id1                         #merge/death
  #old_id3  old_id3 old_id2 old_id4 old_id5 #merge
  #old_id2  old_id3 old_id4 old_id5         #merge/death
  #old_id4  old_id3 old_id2 old_id5         #merge/death
  #old_id5  old_id3 old_id2 old_id4         #merge/death
  #old_id7                                  #death
  #         new_id1                         #birth

  #do we need these merge/death instance recorded?
  #as information is encoded in the merge recods



  #Here!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  #consider the above
  #revise first block and tidy up annotation
  #then implement printing
  #NEARLY FINISHED!!! yeah right




  


  #do we need this last file can we not merge this with above as stable_ID to stable_id mapping?
  #contains only new stable ids and the parents
  #$new_id_file = open_file($out_dir.'/new_stable_id.mappings');#is this one needed?
}


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

  my ($o_start, $o_end, $overlap, $displacement, $num_parents, $current_child, $previous_child, @overlap_ids);
  my $overlap = 0;
  my $no_new_features = 0;
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
  
  

  ####Need to overhaul the above and below block!!!!!!

  if (! @nreg_feats) {
	
	#consider extend...as above
	

	#This is the first old to new mapping which has failed
	if (! defined $source_dbID) {
	  print "No $to_db mapping possible for:\t".$mapping_info->{'stable_id'}."\n";

	  warn "we need to record retirement of stable_id here";

	  next;						#there should only be one pair in this loop, so will exit loop
	}
  }


  #if(@nreg_feats){

  for my $i (0..$#nreg_feats) {
	#if we're in the middle/end of a chain, the first will always be the source_dbID/dbID
	
	if (defined $source_dbID && $#nreg_feats == 0) { #we have found the end of a chain;
	  $no_new_features = 1;
	  last;
	}
	
	push @overlap_ids, ($to_db eq 'NEW') ? $nreg_feats[$i]->dbID() : $nreg_feats[$i]->stable_id();
	
	#overlap length is always the difference between the two middle sorted values
	(undef, $o_start, $o_end, undef) = sort($slice->start(), $slice->end(), $nreg_feats[$i]->start(), $nreg_feats[$i]->end());
	my $tmp_length = $o_end - $o_start;
	my $tmp_displacement
	  my $new_child = 0;
	
	if ($tmp_length > $overlap) {
	  $new_child = 1;
	} elsif ($tmp_length == $overlap) {
	  #we need to account for %age overlap and number of consituents here?
	  #o-a- ---c--
	  #n -b---
	  #this would result in a>b(retire c)
	  
	  if ($nreg_feats[$i]->length() > $coverage) {
		$new_child = 1;
	  } elsif ($nreg_feats[$i]->length() == $coverage) {
		warn "Coverage is not resolving stable_id inheritance, need more rules, defaulting to 3' inheritance for \n";
		$new_child = 1;			#default to 3' inheritance!!!!!
		#We need to resolve this, or just opt for >= in above statement for now
		#this would give 3' presidence to inheritance
		
		#?????
		#would also need to acount for equal %age overlap and look at number of overlaps?
		#what would happen here?
		#o   -------  -------
		#n    --  ------  -
		#what if this does not resolve?
		#do we need another rule?
		#o------------
		#n    --    --
		#arguably the one in the middle should inherit


		#will we need to store coverage etc for resolving complex inheritance
		#in the reverse assign loop?
		#maybe not coverage, but maybe other info derive from the last transition 
		#which will not be available from the linking overlap, which will be available in the next transition
		#i.e. num of parents?
		#o -a c1 --c2----
		#n  --b-----  --d-----
		#Here b should in herit from c and not d
		#as their coverage is equal but b has other parental support
		#what about?
		#o -a c1  --c2--
		#n  --b-----  --d-----
		#b would inherit from c1 due to coverage
		#and this?
		#o  -a c1 --c2--
		#n  --b-----  --d-----
		#again b would inherit from c1, a, c1 and c2 have equal overlap, 
		#a and c1 have equal coverage, but c1 is more central to b.
		#is this the same in reverse?
		#o---a---- --c--
		#nb1 b2 -b3--
		#this is the same as we don't account for direction here!!??
		#this will take some funky logic in the assign block

		#do mid point comparison first as this will always find the best contained feature first
		#as these would fail the number of parents test as they will only have one
		#conversely the midpoint test for linking features would always fail as 
		#linking features will always score equally on the midpoint after scoring equally on the coverage
		
		if ($tmp_displacement < $displacement) {
		  $new_child = 1;
		} elsif ($tmp_displacement == $displacement) {
		  warn "Displacemnt is failing to resolve stable_id mapping, need more rules!!";
		  $new_child = 1;
		  #this is slightly different and may need to be worked out in the assign block!!!!!!!!!!!!!!!!
		  
		  #if($tmp_num_parents > $num_parents){

		  #}elsif($tmp_num_parents == $num_parents){
		  warn "Number of parents is failing to resolve stable_id mapping, need more rules";

		  #o -a c1 --c2---- x1 x2
		  #n --b------  --d------

		  #b should inherit from c2 due to increased number of children/parents
		  #so when considering the c2> b & d transition
		  #we currently don't know the number of parents of d or d
		  #so we can't pick a current child
		  #if they are truly equal then it's probably best to retire c2
		  #and assign c1 to b and x1 to d

		  #STOP HERE! THIS IS GETTNIG WAY MORE COMPLICATED THAN WE CAN HANDLE IN ONE GO!
		  #num_parent is essentially $#overlaps, so we will need to trace this backwards 
		  #through 3 transitions to be able to resolve this 
		  #d  > c2(x1 & x2) 
		  #c2 > b & d Don't know which is the child so can't assign an id yet
		  #b  > a > c1 & c2 can now compare num parents and assign if there is a difference
		  #else split c2 and use c1 and x1
		  
		  #what if parents have differing coverage?
		  #o--a c1 --c2---- x1 x2
		  #n --b------  --d------
		  #b now inherit from a
		  #and d from c2

		  #what is a has better match?
		  #o-----a c1 --c2---- x1 x2
		  #n-z- --b------  --d------
		  #we would have to work right back to z before assigning ids
		  #z should inherit from a
		  #b should inherit from previosu child i.e. c1
		  #d should inherit from c2

		  #That has to be it!!
		  #another extension here would be caught by preceding rules no?


		  #HERE!!!!!!!!!
		  #Right, no go back and update the first assign block in line with the bottom one, and
		  #add rules for non-maximal linking overlaps!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		  #THen maybe we can start logging these little fellas
		  #also make sure we're returning the right source_dbID and feature




		  #}

		


		}
	  }
	}
  
	#elsif($nreg_feats[$i]->end() > $slice->end()){
	#	#we have a non-maximal end overlap i.e. a potential true orphan
	#	#o----   -------
	#	#n  ------
	#	
	#	#we don't need to set true orphan
	#	#as we know if we are not using previous and the child isn't $# then we have a non-maximal link
	#
	#	  }


	#Found 3' overhanging feature
	if($nreg_feats[$i]->end() > $slice->end()){
	  #set source/parent dbID for next transtion
	  #and set new feature to overhang feature
	  $source_dbID = $feature->dbID();
	  $feature = $nreg_feats[$i];
	}else{
	  #make sure we don't pass back current values and get caught in loop
	  undef $source_dbID;
	  undef $feature;
	}

	#Found new child so update vars
	if($new_child){
	  $overlap = $tmp_length;
	  $coverage = $overlap / $nreg_feats[$i]->length();
	  $previous_child = $current_child if (defined $current_child);
	  $current_child = $i;
	}
		


	# $end_overlap = 1 if($nreg_feats[$i]->end() > $slice->end());
	#Don't need to do these check now as we are handling all nreg_feats/children/orphans
	#if we are going to consider a 5' and 3' extend we should do this as $i == 0 and $i == $# !!!!!!!!!!!!!!!!!!!!!!!

	#Some sensible checks
	#could we add more to check we aren't getting overlapping features within the same set
	#i.e. check current start > last end 
	#if ($nreg_feats[$i]->start() < $slice->start()) {
	#This can only be a 5' new feature
	#as we will have dealt with all 5' old features
	#o  -----
	#n --
	#	if ($to_db eq 'OLD') {
	#	  throw("Found 5' feature from $to_db, should have already dealt with this");
	#not quite true as we could be extending
	#n--    -----
	#o --
	#we should sill never see this as the old feature will have been extended first and captured this
	#new feature, hence we will have skipped it
	#	} elsif (defined $source_dbID) {
	#found 5' new feature which wasn't the source feature we have skipped
	#this should not be possible
	#	  throw("Found 5' feature from $to_db which does not match the slice source feature");
	#	}
	#Have bonafide 5' mapping
	#biggest overlap
	#these will both be stored in the dbID array
	#o -------
	#n--
	#o   ----
	#n ---- --
	#smallest overlap
	#would not be stored in array...can we just use the first overlap ID?
	#We don't know whether this is 5'or contained
	#we need to store dbID even if it is not the biggest overlap due to 
	#inheritance issue e.g.
	#o ---- ------
	#n-- -------
	#ug we also nee to consider other overlaps
	#o  -------------- ----------
	#n----  ---   ----------
	#this is only an issue fo rthe first transition, due to the nature of identifying chained overlaps?
	#o--a----- ---c----------  --e-------
	#n      -b--- -d1- -----d2---------
	###NO !! GASH!!! WE NEED TO CONSIDER THIS IN ALL MODELS!!!!!!!!!!!!!!!
	#Actually this is another case for splitting the mappings
	#We should always split this into
	#a>b we need to keep info about b > c but cannot have this as the longest overlap
	#as despite what might happen with inheritance further down the chain, d1 will 
	#always inherit over b
	#do we just add b>c, but do not add longest overlap?
	#mmm no I think we need, current child and previous child?
	#as we can't assume it will be previous chaining overlap
	#	throw("Found more than one 5' feature") if(scalar(@new_feats) >1);
	#  }
	#}
	

	#need to finish off revising block below
	#we're still not accounting for the following
	#o-a-  --c----      ---e-
	#n  --b--    ---d----  xxxx
	
	#we would split on c/d here
	#but this would result in not knowing what the stable_id for d is when building
	#the c>b & d mapping
	#in the case of e mapping to d we assign e, but if e maps to xx then we assign a
	#new stable id
	#likewise with the following?
	
	#o-a--  --c-----------   ---e-
	#n  --b--  -d1-  ---d2----
	
	#here we would split on b>c
	#but we know what the stable_id of c is, so this is not an issue
	
	
	#This is that same issue as when we work backwards when assigning..but obviously in reverse
	#so we should not split chain at all, but simply account for this in the assignment rules!!!!!!!!!
	#ARRRGGGG WILL THIS THING EVER END?
	#we need to mark these features in the mapping_info hash 
	#as we won't know it is a a true orphan(i.e.has two parents which want nothign to d with it)
	#rather than an single parent orphan
  }



  #Found end of chain so assign ids
  if($no_new_features || ! @nreg_features){
	#no_new_features ensures we don't add the last transition which has only parent as nreg_feat
	&assign_stable_ids($mapping_info);		
  }else{
	#Found some new overlapping features i.e. a transition

	push @{$mapping_info->{'transitions'}},	
	  {(source_id        => ($to_db eq 'NEW') ? $feature->stable_id() : $feature->dbID(),
		source           => $from_db, #current feature feature_set, OLD or NEW
		#overlap          => undef, #bp overlap with current child
		#coverage         => $coverage, #% of current child overlapping feature
		overlap_ids      => \@overlap_ids, #dbID or stableIDs of consitutent features of transition
		current_child    => $current_child, #index of current child in tmp new_feats array
		previous_child   => $previous_child, #used to identify inheritance if it get's shifted
		true_orphan      => ??????,
	   )};
  }

  
  #finally return the mapping info hash and switch the db
  return ($mapping_info, $to_db, $feature, $source_dbID);
}
 


sub assign_stable_ids{
  my $mapping_info = shift;
  
  my ($i);
  my $num_trans = $#{$mapping_info{'transitions'}};
  my $stable_or_dbid;
  my $use_previous = 0;
  my $child = 0; #default to 0 to avoid undef warning below
  #This is also making the assumption that the first or last transition, should use the current child 
  #as there are no further inheritance/projection issues


  #are we splitting above when previous child !=0 and shift == 1?
  #do we need to do this now, or can we work through a chain using the current/previous_child values
  #taking into account what the last(opposite) transition did?
  #i.e. if the child used was 0 we know that the current child should be used
  #else we know the last transition has a larger overlap somewhere and so will not inherit/project
  #to this feature, so we use the previous child

  #We're not splitting on multiple overlap so we need to account for the following
  #this is handled by previous child

  #o  --a-- --c1 ---c2-------------- xxxxxxxxxxxx
  #n     -----b-------        ---------d-------


  #basically we need to set use_previous if last child != 0

  #If the last child was 0 then the next child will be $# or current_child?
  #This only true is we maintain the splitting of non-maximal overlapping features

  #d>c(xxx) current_child = 0(1) 
  #therefore next time last_child will be 0(1) hence
  #c2>d & b = current_child = 1 or is xxx then previous child = 0


  #we could accomodate non-maximal overlaps by seeing if current_child == $#
  #if use_previous and current_child != $# then we would still use current
  my ($last_stable_id, %split_merge_cache);
  #split merge cache has old stable ids as keys
  #values are an array of new or inherited stable ids
  #first of which is always the id this has mapped to, the rest are other ids this feature has been split across
  #so first id maybe the same as the key, or a different inherited stable id
  #and the maybe >= 0 other ids.



  for($i = $num_trans; $i >= 0; $i--){
	my ($orphan, @new_stable_ids);
	my $last_orphan = $#{$mapping_info{'transitions'}[$i]{'overlap_ids'}};

	$use_previous = 1 if($child != 0);

	#reset child for current transition
	$child = ($use_previous) ? $mapping_info{'transitions'}->[$i]->{'previous_child'} :
	  $mapping_info{'transitions'}->[$i]->{'current_child'};

	$child_id = $mapping_info{'transitions'}[$i]{'overlap_ids'}[$child];
	
		
	if ($mapping_info{'transitions'}[$i]{'source'} eq 'NEW'){#inheritance
	  #assign stable_id
	  #record other merged stable_ids


	  #no no no, this might not be true, as the previous transition
	  #may have projected it's stable_id to this new feature
	  
	  if(exists $dbid_mappings{$mapping_info{'transitions'}[$i]{'source_id'}}){
		throw("Found duplicate inheritance:\tdbID ".$mapping_info{'transitions'}[$i]{'source_id'}." > ".
			  $dbid_mappings{$mapping_info{'transitions'}[$i]{'source_id'}}." & $child_id");
	  }
	
	  #key is dbID
	  $mappings++;
	  
	  #we don't need to cache this unless we want to check if we're getting duplicate inheritance
	  #therefore we only need to cache the dbID => stable_id
	  
	  $dbid_mappings{$mapping_info{'transitions'}[$i]{'source_id'}} = $child_id;
	  
	  #{(
	  #	stable_id  => $mapping_info{'transitions'}[$i]{'overlap_ids'}[$child],
	  #	merged_ids => $mapping_info{'transitions'}[$i]{'overlap_ids'},
	  #   )};
	  
	  
	  #print to file/log in DB here?
	  if($dump){
		
		#this is: dbid   mapped_stable_id  merged_stable_ids
		
		#print $mapping_file $dbid_mappings{$mapping_info{'transitions'}[$i]{'source_id'}}."\t".
		#  $child_id."\t".join("\t", @{$mapping_info{'transitions'}[$i]{'overlap_ids'}})."\n";#NR for mapped stable_id
		
		#oxxx  ---a-----
		#n  -b1- -b2- -b3-
		
		#we need to push the mapped_stable_id(x>b1) onto the previous split cache
		#as we won't know the b1 stable id until this point
		#previous split/merge may not always be valid

		#what about?
		#o ---a--   --c-
		#n     --b---
		#or infact omit c
		#then cache will be absent as this will be the first transition
		#so we can test whether the split/merge cache exists
		#it should have at least one entry otherwise we have the above
		#we need access to the last stable_id e.g. a.
		#which will always be the 2nd element in the overlap IDs
		#otherise this transition would not exist
	
		#my 
	
		#shoudl we loop through retiring stable id which are merged or split but not inherited?
		#is this not done in the other loop too?
		#no this is creating new stable ids not retiring old ones
		#so we do a similar thing
		#but we need to account for the last 5' transition here?

		#ox1xx x2x ---a---------
		#n  -b1------- -b2- -b3-


		#we know that if we are not the last transition then we shoud always have a split_merge_cache
		#we retire all non-child ids
		
		#now do very similar loop to below but for retiring IDs, no creating them
		#only populate previous cache with current stable_id if current transition is not $#
		#retire all contained stable ids, except the child/parent
		#
		#Here we see that if there is another transition, this would capture all the merge split info for 1x
		#

		#b1 may already have been assigned a stable_id, by the previous trainsition?????
		#nope! as we would be using previous!!!
		#however we do need to make sure that the last overlap stable_id is retrieved from the dbID cache
		#unless we are at the last transition!
		
   
			
		#deal with new stable_ids e.g. db
		
		#now populate the previous split_merge_cache
		


		########################IS THIS RIGHT?????
		#$start_orhpan = ($i == 0) ? 0 : 1;
		#should always be 0?????????????????????????????????????????????????
		#as we always know the stable ID
		#are we missing the reveser case with respet to dbIDs????


		#we need to push the child stable_id onto the split cache for the last orphan
		#only if this is not the last transition 
		#(and the current_child = $#overlaps
		#this will current always be the case as we're splitting chains if the link is not the maximal overlap)
		
		
		foreach my $orphan_cnt(0..$last_orphan){#these are stable_ids!
		  $orphan = $mapping_info{'transitions'}[$i]{'overlap_ids'}[$orphan_cnt];
		  #orphan is stable_id
		  
		  #next if $orphan == $child_id;#this is causing problems here!!!
		  #we need to be able to populate the previosu cache in this loop
		  #is this causing similar probs below??

		  
		  #merge all, with inherited id first
		  if($orphan == $child_id){
			unshift @stable_ids, $orphan;
		  }else{
			push @stable_ids, $orphan;
		  }



		  #do we really want to do this?
		  #o1xx  ---a--------- xxxx
		  #n  -b1--- -b2-  -b3------
		  
		  #use previously assigned stable_id for linking overlap if it exists
		  
		  #we need to create a new stable_id for b1 if it is not the child!!
		  #as we could have following
		  #o 1xxx  ---a--------- xxxx
		  #nxxx -b1-- -b2- -b3------
		  
		  #so here we are b1, but it is not child
		  #this would never happen as it would get split!!
		  

		  if($orphan_cnt == $last_orphan && exists $split_merge_cache{$orphan}){
			   #exists $dbid_mappings{$mapping_info{$orphan}}){
			  #this is analgous to checking that child != $#overlap

			  #record this split from
			  #o1xx  ---a--------- xxxx
			  #n  -b1-- -b2- -b3------

			  
			#we need to add b1 i.e. 1x or current_child to the split_merge_cache of a

			#only if b1 wasn't the child of a
			#we know this as we should be using_previous if this was the case

			push @{$split_merge_cache{$orphan}}, $child_id if ! ($use_previous);

		  }else{
			#in fact do we need the above block?
			#as we're essentially doing the same for all the orphans
			#adding the child_id
			#The difference is that they should be at the start of all the non-inherited orphans arrays
			#altho' there should only ever be one member, except for the potential overlap at the 5' end
			#which we will deal with in the next loop
			@{$split_merge_cache{$orphan}} = $child_id;
		  }

		  #now add individual contained stable_ids to split merge cache
		  #here we are, this is we're we need to cosider the next transition
		  #as the 5' end may be a split rather than a simple merge, hence we need to catch this in the next block



		  #####################HERE WE NEED TO NOW CHECK IN THE NEXT TRANSITION
		  #WHETHER TO PUSH OR UNSHIFT THE next stable ID depending on inheritance
		  #o xxxxxx  ---a----- xxxx
		  #nxx x -b1--- -b2- -b3-

		  #THis is already done by including the child_id as the first in the stable_id array
		  #Then skipping the child in the orhpan loop
		  #and using a pre-existing 5' stable or if present else generating a new one




		}


		print $mapping_file $dbid_mappings{$mapping_info{'transitions'}[$i]{'source_id'}}."\t".
		  join("\t", @stable_ids)."\n";#no longer NR for mapped stable_id


	  }
	  #else{
	#	#do the db stuff here
	#	my $sql = 'UPDATE annotated_feature';
	#	#would it not be better to dump to tab delimited text and load as one file?
	#  }
	  


	  #set last old stable_id here
	  #this is to enable updating of the overlapping new stable_id in the split merge cache
	  $last_stable_id = $orphan;

	}
	else{#projection
		
	  #stable_id
		my $outline = $mapping_info{'transitions'}[$i]{'source_id'};
		
		if(! defined $child){
		  #this can only be
		  #o --a-- --c----
		  #n     -b---
		  #we're a has no previous_child
		  #therefore retire!


		#or could be
		#o ---
		#n

		#Is this not caught above?
		#Maybe we should just pass it anyway and let this bblock take care of the logic


		$outline .= "\t".$child_id if (defined $child_id);



		#we need to record this as a retired stable ID and children or no children
		#The 'merge' behaviour will be captured by the b>c(a) inheritance transition


		#we need access to the stable_id of b

		#we only know this after we have done the previous transition assignment, so we need to store previous stable_id

	  }else{
		#o -a--- -d---
		#n-b- -c----


		#so this needs to record split of
		#a>b & c  where a will project it's stable id to b

		#this will not be caught by above block


		#will this d>e?

		#o -a--- -d-----------  xxxxxx
		#n-b- -c----  -e1- --e2--
		#n                  ---e3-
		
		#we still need to record the split here

		
		#we also need to record the dbid_mappings
		#as the next transition c>d & a will not capture this
		#as c will inherit from a due to coverage of e2 being > c
		#e1 should get a new stable_id
		# we need to account for this here too

		#however if we consider e3 over e2
		#then we will record the d to c inheritance here
		#and when we get to the next transition
		#we can't record a&d > c fully here as we don't know about a yet
		#so if current_child == 0 then we need to leave it until the next transition?
		#and just record d > c e1 e3

		#else it would be d > e2 c e1
		#This would need recording the dbID hash only if we haven't already done
		#e2> d xxx


		#so basically we need to record in the dbID hash here if the child is not 0 or $#
		#this would be a fully contained dbID, therefore we would not have any other merged stable_ids


		#also need to record if
		#1 we are at the last(3') transition and child is $# and $# !=0
		#2 we are at the first(5') transition and the child == 0 and $# != 0

	
		#we also need to record split
	
		#don't need this if now as the loop accounts for position of transition
		
		#if($child != 0 && $child != $#{$dbid_mappings{$mapping_info{'transitions'}[$i]{'overlap_ids'}}}){
		 
		$mappings++;
		#oxxx  ---a----- xxxxxxxxx
		#n  -b1- -b2- -b3--- xxx xx
		
		#a>b2 (b1, b3)
		
	 
		#this will not fail here
		#child will never be b3(which would already have a stable_id) as we will have specified use_previous
		#otherwise the chain would be split as the linking overlap would not be maximal
		
		if(exists $dbid_mappings{$child_id}){
		  throw("Found duplicate inheritance:\tdbID $child_id > ".$mapping_info{'transitions'}[$i]{'source_id'}.
				" & ".$dbid_mappings{$child_id});
		}
		
		#record dbID mapping first only done to check for duplicate inheritance
		$dbid_mappings{$child_id} = $mapping_info{'transitions'}[$i]{'source_id'};
		
		#then record splits?
		#can not record split until we have dealt with new stable id generation
		#and potential inherited stable_id from x>b1
		
		
		
		#can we sub this??
		
		#make sure inherited stable_id is first

		#@new_stable_ids = ($dbid_mappings{$child_id});
		
		#deal with new stable_ids e.g. db
		
		#my $start_orhpan = ($i == 0) ? 0 : 1;#last/5' transtion so deal with b1!
		
		#we normally deal with b1 in the next transition
		#but if b1 is the current child then will this get duplicated
		#we need to check that we don't add it again to the split merge cache of a
		#if we already have a stable_id


		#we need to deal with b1 if it is the current child as we are not adding this by deafult
		#so we always need to consider 0, but skip if it is not the current child or is i > 0




		#what if b1 inherits here due t use_previous
		#oxxx  ---a--------- xxxx
		#n  -b1--- -b2- -b3------
		#This will always be added as we include the child_id first!!
		
		foreach my $orphan_cnt(0..$last_orphan){
		  #orphan is dbID
		  my $orphan = $mapping_info{'transitions'}[$i]{'overlap_ids'}[$orphan_cnt];

		  #Rules for dealing with 5' position
		  #we need to account for 5' link not having stable_id assign yet
		  #b1 at 0 in the above case will inherit from the next transition or xxx
		  #so if we are the first transition(5') we assign stable_id
		  #else we populate the split merge cache for a with b2 and b3, and finish in the next transtion with xxx
		  
		  #Rules for accounting for 3' position
		  #if we know 3'($#) orphan has already been assigned a stable_id
		  #hence we know we should be using previous and will not be current_child
		  #so we use the previously assigned stable_id in our current stable_id array
		  #else if this current child, it could be mid chain or the last 3' transition
		  #so we project the stable_id and update the split_merge cache for the previous old stable_id
		  #else if it is not the current child it must be the 5' position
		  #otherwise the chain would have been split, so we create a new stable_id
		  
	
		  #Rules for all others: we assign new stable_ids
		
		  #3' rules
		  if($orphan_cnt == $last_orphan){

			if(exists $dbid_mappings{$mapping_info{$orphan}}){
			  #This will have already been added to the previous split cache as it will be the previous child id
			  $new_sid = $dbid_mappings{$mapping_info{$orphan}};
			}
			elsif($orphan == $child_id){
			  unshift @stable_ids, $dbid_mappings{$child_id};
			  
			  #we need add this to the previous split_merge_cache
			  if(defined $last_stable_id){

				if(exists $split_merge_cache{$last_stable_id}){
				  push @{$split_merge_cache{$last_stable_id}}, $dbid_mappings{$child_id};
				}
				else{
				  throw("Failed to find split_merge_cache for old $last_stable_id when updating with linking new feature ".
						$dbid_mappings{$child_id});
				}
			  }
			 

			  next; #to avoid adding the last new_sid again
			}
			elsif($i == $num_transitions){
			  $new_sid = &get_new_stable_id();
			  $dbid_mappings{$orphan} = $new_sid;
			}
			else{
			  throw('This should never happen as this means a link has been made with a none maximal overlap, as link has not been used as child for either transition');
			}		
		  }
		  elsif($orphan == $child_id){
			#this will also account for the child being the 5' orphan

			unshift @stable_ids, $dbid_mappings{$child_id};
			#don't need to add this to the split_merge_cache as it is not overlapping
			next;# to avoid adding the last new_sid

		  }
		  elsif($i > 0 || ($i == 0 && $orphan_cnt == 0)){#generate a new stable_id
			$new_sid = &get_new_stable_id();
			$dbid_mappings{$orphan} = $new_sid;
			print $mapping_file $orphan."\t".$new_sid;

			#we should have contributing stable_ids here too
			#is this not duplcaiting the split_merge_cache?
			
			#we need to log the dbID and stable_id associations, simply to be able to update the features from file?
			#what if we are just going to load a new set from file, need to dump the correct tab formatted data, but don't require
			#the dbID...DO REQUIRE THE NEW FEATURE/S SO WE CAN DUMP/LOAD THE INFO CORRECTLY

			
			#need to do correct printing and updating here
			#or store all in tmp caches and print in one go after w've finished all the transitions?
		  }
		  else{
			next;
			#5' orphan in the middle fo a chain
			#don't now stable id, update split mrege cache in next transition
		  }
				
		  #print to new stable_id file here
		  #also record in dbID hash
		  push @new_stable_ids, $new_sid;
			
		}

		#set the cache, may not contain the most 5' stable_id which may need to be assigned in the next transition
		@{$split_merge_cache{$mapping_info{'transitions'}[$i]{'source_id'}}} = @new_stable_ids;
	  }
	}
  }


  #so we should prit out split_merge_cache here
  #we should also print out dbID > stable id mappings too, in loop or here, and do we need other constributing/merged stable_ids in this file
  #lastly we need to print out retired stable ids? shouldn't this go in the split merge cache?
  


  #ignore all this old id array logic!!!!!!!!



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

  #store 2a dbID in the mapping hash(as previous_child index) and if we get -1 then we check 
  #it exists and use it else we retire.

  #we're now using the previous_child element to enable spliting over mappings based on a child which is not 
  #a joining overlap.
  


  #need to test this set up
  #o   --a--    --c--
  #n-b1---  --b2-----
  #this would get split on b2

  #and test  if stable_ids < dbIDs

  #o  -a--  -c----   ---e---
  #n     ---b- ---d------
  #inheritance  transition pos/elem + num seen OLD transitions - 1
  #3 + 0  -1 == 2  CORRECT!!!
  #d>e(&c)  >  e == 2
  #shifted -2 == CORRECT

  #b>c(a) == 1
  #1 + 1 - 1  CORRECT!!
  
  #projection
  #transition pos/elem + num seen NEW transition - difference between stable_id and dbID lengths - 1
  #a>b  nothing / -1
  #0 + 1 - 1 -1 = -1 CORRECT

  #shifted
  #0 + 1 - 1 = 0 CORRECT!!!
  
  
  


  #as we work backwards, we're dealing with every mapping twice due to nature of transitions
  #can we avoid this some how?
  #can we only store forward transition, and rely on this to record all mappings?
  #this would mess up the logic above?


  #Works, but overkill and messy, use simple previous child rules above


}



1;
