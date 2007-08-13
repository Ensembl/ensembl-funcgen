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

what is the concept of a RegulatoryFeature?
It can cross cell/tissue types but have differing modifications, hence states
but maintain same stable_id
This makes mapping between assembly variable region extremely hard
we need to take account of multi cell/tissue types in some way.
we need to factor this into the schema

=cut


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
my ($dbname, $help, $man, $slice_name);
my ($odbname, $ndbname, $npass, $nuser, $from_file, $ohost, $oport);
my $reg = "Bio::EnsEMBL::Registry";

my $ouser  = 'ensro';
my $opass = '';
#my $ohost = 'ens-livemirror';
my $nhost = 'ens-genomics1';
my $old_fset_name = 'RegulatoryFeatures';
my $new_fset_name = 'newRegulatoryFeatures';
my $stable_fset_name = 'stableRegulatoryFeatures';
my $nport = '3306';
my $no_load = 0;
my $from_file = 0;
#my $update = 0;
my $expand = 0;
my $out_dir = '.';


$main::_debug_level = 0;
$main::_tee = 0;


GetOptions (
			"odbpass=s"     => \$opass,
			"odbport=s"     => \$oport,
			"ohost=s"     => \$ohost,
			"ouser=s"     => \$ouser,
			"odbname=s"   => \$odbname,
			"ndbpass=s"     => \$npass,
			"ndbport=s"     => \$nport,
			"nhost=s"     => \$nhost,
			"nuser=s"     => \$nuser,
			"ndbname=s"   => \$ndbname,
			"expand|e=s"  => \$expand,
			"no_load"        => \$no_load,
			"slice_name=s" => \$slice_name,
			"outdir=s"    => \$out_dir,
			"update"      => \$update,
			"from_file|f=s" =>\$from_file,
			"help|?"       => \$help,
			"man|m"        => \$man,
			#add opt for old, new & stable fset name
		   );


pod2usage(1) if $help;
pod2usage(-exitstatus => 0, -verbose => 2) if $man;


#check mandatory params here
throw("You cannot specify load '-from_file' in conjunction with '-no_load'\n") if($no_load && $from_file);


#we need to resolve odb ndb relationship, defulat to same if only one specified?
#also need to add defaults for 

if(! ($nhost && $nuser && $npass && $ndbname)){
  throw("You must provide new DB details: -ndbhost -ndbuser -ndbpass -ndbname\n");
}

my $ndb = Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor(
												  host   => $nhost,
												  user   => $nuser,
												  pass   => $npass,
												  port   => $nport,
												  dbname => $ndbname,
												 );

#Test connection, eval this?
$ndb->dbc();

if(! ($odbname && $ohost && $ouser && $oport)){
  warn "Defaulting some old DB setting to new DB values\n";
  $oport ||= $nport;
  $odbname ||= $ndbname;
  $ouser ||= $nuser;
  $ohost ||= $nhost;
}

my $odb = Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor(
													host   => $ohost,
													user   => $ouser,
													pass   => $opass,
													port   => $oport,
													dbname => $odbname,
												   );
$odb->dbc();#Test connection, eval this?


my ($next_stable_id, $stable_id, $new_reg_feat_handle);
my ($new_id_handle, $stable_id_handle, %dbid_mappings, @slices);
my (%mapping_cache, %obj_cache, $mappings, $new_mappings);#%new_id_cache?
my $total_mappings = 0;

###Set up new stable ID mapped FeatureSet
my $fset_a = $odb->get_FeatureSetAdaptor();
my $stable_fset = $fset_a->fetch_by_name($stable_fset_name);

if(defined $stable_fset){
  
  if(! defined $clobber){
	throw('The FeatureSet $stable_fset_name already exists in the DB,'.
		  ' please redefine the -feature_set or specify -clobber to overwrite');
  }else{
	$odb->dbc->db_handle->do('DELETE from annotated_feature where feature_set_id='.$stable_fset->dbID);
  }
}
else{
  my $reg_ftype = $odb->get_FeatureTypeAdaptor->fetch_by_name('RegulatoryFeatures');
  my $reg_atype = $odb->get_AnalysisAdaptor->fetch_by_name('RegulatoryRegion');
  
  $stable_fset = Bio::EnsEMBL::Funcgen::FeatureSet->new(
														name         => $stable_fset_name,
														analysis     => $reg_atype,
														feature_type => $reg_ftype,
													   );
  
  ($stable_fset) =  @{$fset_a->store($stable_fset)};
}



  
##Set start stable_id
$next_stable_id = $odb->dbc->db_handle->do('SELECT display_label from annotated_feature af, feature_set fs where fs.feature_type_id='$oftype->dbID().' and fs.feature_set_id=af.feature_set_id order by display_label desc limit 1');
$next_stable_id =~ s/\:.*//;
$next_stable_id =~ s/ENSR[A-Z]*0*//;
$next_stable_id ++;

###Get Old and New FeatureSets and set adaptors
$obj_cache{'OLD'}{'AF_ADAPTOR'}    = $odb->get_AnnotatedFeatureAdaptor();
$obj_cache{'OLD'}{'SLICE_ADAPTOR'} = $odb->get_SliceAdaptor();
$obj_cache{'OLD'}{'FSET'}          = $odb->get_FeatureSetAdaptor->fetch_by_name($old_fset_name);
$obj_cache{'NEW'}{'SLICE_ADAPTOR'} = $ndb->get_SliceAdaptor();
$obj_cache{'NEW'}{'AF_ADAPTOR'}    = $ndb->get_AnnotatedFeatureAdaptor();
$obj_cache{'NEW'}{'FSET'}          = $ndb->get_FeatureSetAdaptor->fetch_by_name($new_fset_name);



if($slice_name){
  @slices = ( $obj_cache{'OLD'}{'SLICE_ADAPTOR'}->fetch_by_name($slice_name) );
}  
else{
  @slices = @{$obj_cache{'OLD'}{'SLICE_ADAPTOR'}->fetch_all('toplevel')};
}

#Process each top level seq_region
foreach my $slice (@slices){
  my $seq_name = $slice->seq_region_name();

  #Output files
  my $new_reg_feat_file = $out_dir.'/annotated_feature.'.$seq_name.'.txt';
  #tab txt dump of new stable features

  my $new_id_file = $out_dir.'/new_stable_ids.'.$seq_name;
  #This is slightly redundant? This will be present in the db_id_stable_id file 
  #and absent in the stable_id file
  #         new_id1                         #birth

  my $stable_id_file = $out_dir.'/stable_id_mappings.'.$seq_name;
  #this contains retired ids and the their child ids if any
  #old stable new/old stable id as focus + parents
  #old_id   new_id/split  parents/merges
  #old_id1  old_id1       old_id2                 #merge
  #old_id2  old_id1                         #merge/death
  #old_id3  old_id3       old_id2 old_id4 old_id5 #merge
  #old_id2  old_id3       old_id4 old_id5         #merge/death
  #old_id4  old_id3       old_id2 old_id5         #merge/death
  #old_id5  old_id3       old_id2 old_id4         #merge/death
  #old_id7                                        #death
 
  #Start the stable ID mapping
  if(! $from_file){
	#open files
	$new_reg_feat_handle  = open_file($new_reg_feat_file);
	$new_id_handle = open_file($new_id_file);
	
	#clean vars
	%mapping_cache = {};
	%dbid_mappings = {};
	$mappings = 0;
	$new_mappings = 0;

	foreach my $reg_feat(@{$obj_cache{'OLD'}{'AF_ADAPTOR'}->fetch_all_by_Slice_FeatureSet($slice, $obj_cache{'OLD'}{'FSET'})}){
	  my $from_db = 'OLD';
	  my $feature = $reg_feat;

	  #Can we change this to an array as it's only containing transitions array and the orignal stable_id
	  my $transitions = []; #this has to be ref due to implementation in while



	  my $source_dbID; #The source dbID this feature was fetched from, set to null for start feature.
	  ($stable_id = $reg_feat->display_label()) =~ s/\:.*//;
	  next if exists $mapping_cache{$stable_id};
	  
	  #keep building transition whilst we have a 3' overhanging feature returned
	  while(defined $feature){		
		#we don't need to sub this part do we?
		($transitions, $from_db, $feature, $source_dbID) = &build_transitions($transitions, $from_db, $feature, $source_dbID);
	  }
	}

	print 'Mapped $mappings old stable IDs for toplevel seq_region '.$seq_name."\n";
	$total_mappings += $mappings;

	#Now for this slice get each new feature which hasn't already been mapped, generate a new stable id
	my $new_slice = $obj_cache{'NEW'}{'SLICE_ADAPTOR'}->fetch_by_region('toplevel', $seq_name);

	foreach my $nreg_feat(@{$obj_cache{'NEW'}{'AF_ADAPTOR'}->fetch_all_by_Slice_FeatureSet($slice, $obj_cache{'NEW'}{'FSET'})}){
	  
	  if(! exists $dbid_mappings{$nreg_feat->dbID()}){
		#New regfeat has not been mapped to, so generate new stable_id
		&assign_and_log_new_stable_id($nreg_feat->dbID());		
	  }
	}

	print 'Generated $new_mappings new stable IDs for toplevel seq_region '.$seq_name."\n";
	print 'Total stable IDs for toplevel seq_region '.$seq_name.":\t".($mappings + $new_mappings)."\n";
	$total_mappings += $new_mappings;	

	#print mapping cache
	$stable_id_handle = open_file($stable_id_file);

	foreach my $old_stable_id(keys %mapping_cache){
	  print $stable_id_handle $old_stable_id."\t".(join "\t", @{$mapping_cache{$old_stable_id}})."\n";
	}

	close($stable_id_handle);

	#should we cache other data and print here too?
  }

  #Import from file
  if(! $dump){
	#recover mode should tidy stable id mapping tables
	#can this be done with a large IN statement for all stable IDs on one chromosome
	#we can't recover this per chromosome unless we then get all unassigned new RegFeats on a per chr basis
	#Need to implement this!!!

	my $cmd = "mysqlimport -d${ndbname} -P$nport} -h${nhost} -p{$npass} ".$new_reg_feat_file;
	system($cmd);#this will exit on fail?
	
	#load mappings
	if($from_file){
	  #parse mappings and populate mapping_cache
	  throw("-from_file parse not ye implemented\n");
	}

	foreach my $old_stable_id(keys %mapping_cache){
	  warn "Stable ID event import not yet implemented\n";
	  last;
	  #load each mapping event here
	}
  }
}

print "Finished mapping a total of $total_mappings old and new stable IDs\n" if( ! $from_file);




sub build_transitions{
  my ($transitions, $from_db, $feature, $source_dbID) = @_;
  #the source_dbID is the dbID of the previous source feature
  #to ensure we don't create a recursive transition, hence enabling detection of the end transition

  my ($o_start, $o_end, $overlap, $displacement, $num_parents, $current_child, $previous_child, @overlap_ids);
  my $overlap = 0;#length of overlap
  my $no_new_features = 0;
  my $to_db = ($from_db eq 'OLD') ? 'NEW' : 'OLD';

  my $slice = $feature->slice();

  #disabled default extend for now for simplicity
  # we should only do this if we have no mapping
  #if($expand && $mapping_info->{'iteration'} == 0){
  #	$mapping_info->{'extended'} = 1;
  #	$slice = $slice->expand($expand, $expand);
  #  }
    
  my $next_slice = $obj_cache{$to_db}{'SLICE_ADAPTOR'}->fetch_by_region('toplevel', $slice->seq_region_name(), $slice->start(), $slice->end());
  
  #if(! $next_slice){#I don't think this should ever happen, can we remove?
  #print $reg_feat->display_label()." is on a slice which is not present in the $to_db DB:\t".
  #	  $slice->seq_region_name().' '.$slice->start().' '.$slice->end."\n";
  #	next;
  #  }
    
  my @nreg_feats = @{$obj_cache{$to_db}{'AF_ADAPTOR'}->fetch_by_Slice_FeatureSet($next_slice, 
																				  $obj_cache{$to_db}{'FSET'})};
  
  #disabled default extend for now for simplicity
  #only retry with default extend if we haven't extended already and we don't already have a mapping
  #if(! @nreg_feats && ! $extend && $mapping_info->{'iteration'} == 0){
  #  $next_slice = $next_slice->expand(200, 200);#Just enough to catch the next nucleosome
  #  $mapping_info->{'extended'} = 1;
  #  @nreg_feats =  @{$obj_cache{$to_db}{'AF_ADAPTOR'}->fetch_by_Slice_FeatureSet($next_slice, 
  #$obj_cache{$to_db}{'FSET'})};
  #}
  
  ####Need to overhaul the above and below block!!!!!!

  if (! @nreg_feats) {
	#consider extend...as above
	
	#This is the first old to new mapping which has failed
	if (! defined $source_dbID) {
	  print "No $to_db mapping possible for:\t".$feature->stable_id()."\n";
	  $mapping_cache{$feature->stable_id()} = [];
	  next;						#there should only be one pair in this loop, so will exit loop
	}
	else{
	  throw("this should never happen!");#REMOVE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	}
  }

  #Loop through features assigning best overlap as child
  #recording next best as previous child(?or next best in the case of twins?)
  my $midpoint =  $slice->end() -  $slice->start();
  my $split_chain = 0;

  for my $i (0..$#nreg_feats) {
	#if we're in the middle/end of a chain, the first will always be the source_dbID/dbID
	

	if (defined $source_dbID && $#nreg_feats == 0) { #only found source i.e. end of the chain
	  #we could sanity check here that the source_dbID == the $nreg_feats[0] dbID
	  $no_new_features = 1;
	  last;
	}
		
	#overlap length is always the difference between the two middle sorted values
	(undef, $o_start, $o_end, undef) = sort($slice->start(), $slice->end(), $nreg_feats[$i]->start(), $nreg_feats[$i]->end());
	my $tmp_length = $o_end - $o_start;
	my $tmp_displacement = $midpoint - ($nreg_feats[$i]->end() -  $nreg_feats[$i]->start());
	my $new_child = 0;
	
	#Length test
	if ($tmp_length > $overlap) {
	  $new_child = 1;
	} 
	elsif ($tmp_length == $overlap) {

	  #Coverage test
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
	

		#Displacement test
		$tmp_displacement * -1 if ($tmp_displacement < 0);
	
		if ($tmp_displacement < $displacement) {
		  $new_child = 1;
		  $diplacement = $tmp_displacement;
		} 
		elsif ($tmp_displacement == $displacement) {
		  warn "Displacement is failing to resolve stable_id mapping, need num parent rules!! Defaulting to 5' mapping";
		  $new_child = 1;
		  $diplacement = $tmp_displacement;


		  #STOP HERE! THIS IS GETTNIG WAY MORE COMPLICATED THAN WE CAN HANDLE IN ONE GO!
		  #this is slightly different and may need to be worked out in the assign block!!!!!!!!!!!!!!!!
		  #We need to compare num of parent of parents, but we won't know this until we traverse the 
		  #next transition, hence this needs to be dealt with in the assign block when we walk backwards
		  #or can we do it here as we move along?
		  #as these 'twins' may be exact matches over many transitions
		  #so we need to set a 'twin' flag, and record which transition it occured at
		  #then walk forward comparing to the mirrored transition and seeing which wins, then track back to the
		  #transition where the twin occured assigning child and previous child accordingly
		  #can we leave this walk back until the assign step?
		  #do we also need to record 2nd child instead of previous child, as the next best child might be 3' rather than
		  #the default which is the previous seen i.e. 5'
	  
		  #if($tmp_num_parents > $num_parents){

		  #}elsif($tmp_num_parents == $num_parents){
		  #warn "Number of parents is failing to resolve stable_id mapping, need more rules";

		  #o -a c1 --c2---- x1 x2
		  #n --b------  --d------

		  #b should inherit from c2 due to increased number of children/parents
		  #so when considering the c2> b & d transition
		  #we currently don't know the number of parents of d or d
		  #so we can't pick a current child
		  #if they are truly equal then it's probably best to retire c2
		  #and assign c1 to b and x1 to d

		
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
		  #For walking forward then we find a twin, previous child has to be next best child i.e. 3' to current
		  #we would have to have a nested loop here to compare back to the mirrored 'twins' to define
		  #child/previous child for each before we pass to assign_ids?  Is this necessary, can we not resolves this in assign_ids?
		  #}
		}
	  }
	}

	#non-maximal 3' overlap will be handled by ! use_previous in assign block
	#will never get a unseen 5' new feature when going from new to old, 
	#as we are walking along the chromosome, hence we never have follow transitions 3' to 5'
	#as we will have already seen then


	#Need to split the chain to avoid reverse shift problem
	#If child  is 5'
	#we can't assign this yet until we do the previous transition!!!!
	#o   ---a----  --c--
	#nxxxxxxx ---b--  -d-
	#c > b & d
	#b is child but we don't know whether previous transition will be better match!
	#this could iterate ad nauseum, due to other non-maximal links i.e. ab in context of x being a better overlap

	

	#Found 3' overhanging feature
	if($nreg_feats[$i]->end() > $slice->end()){
	  #set source/parent dbID for next transtion
	  #and set new feature to overhang feature
	  $source_dbID = $feature->dbID();
	  $feature = $nreg_feats[$i];

	  #Found non-maximal link
	  $split_chain = 1 if(! $new_child);
	  
	}
	else{
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
  }


  #Found end of chain so assign ids
  if($no_new_features){
	#no_new_features ensures we don't add the last transition which has only parent as nreg_feat


	#This may be an empty hash as we can potentially break a chain on a non-maximal overlap
	#which gives rise to no more features, will find source

	if(@$transitions){
	  &assign_stable_ids($transitions);
	}
	else{
	  warn "No transitions to assign to, must have found broken non-maximal 3' overhang(dbID:".
		$feature->dbID().") with no more overlapping features\n";
	}
  }
  elsif(@nreg_features){
	#Found some new overlapping features i.e. a transition

	push @$transitions,
	  {(
		#feature_id        => ($to_db eq 'NEW') ? $feature->stable_id() : $feature->dbID(),
		feature          -> $feature,
		source           => $from_db, #current feature feature_set, OLD or NEW
		#overlap         => undef, #bp overlap with current child
		#coverage        => $coverage, #% of current child overlapping feature
		overlap_features => \@nreg_features,
		current_child    => $current_child, #index of current child in tmp new_feats array
		previous_child   => $previous_child, #used to identify inheritance if it get's shifted
		#true_orphan     => ??????,
		#twin?? this could be figured out if previous child is > current_child
	   )};


	#we want to split chain here if last overlap was non-maximal
	if($split_chain){
	  &assign_stable_ids($transitions);
	  #only clean the mapping info as we want to start from the right feature
	  #but we don't want to redo the assignment for the first part of the chain
	  undef $transitions;
	}
	#don't clean source_dbID and feature as we need them to start the next part of the chain
  }
  else{#no mapping/features
	undef $feature;
  }

  #finally return the mapping info hash and switch the db
  return ($transitions, $to_db, $feature, $source_dbID);
}

=head2 assign_stable_ids

  Arg [1]    : Arrayref of transition hashes
  Example    : &assign_stable_ids($transitions);
  Description: Traverses mapping transitions in reverse order assigning
               or creating new stable IDs and caching stable ID history 
  Returntype : None
  Exceptions : Throws if duplicate mappings are found
  Status     : At Risk

=cut

sub assign_stable_ids{
  my $transitions = shift;

  my ($i, $child_id, $last_stable_id);
  my $num_trans = $#{$transitions};
  my $child = 0; #default to 0 to avoid undef warning below

  #Walk backwards through transitions
  for($i = $num_trans; $i >= 0; $i--){
	my ($feature, $orphan, $orphan_id, @new_stable_ids);
	my $last_orphan = $#{$transitions->[$i]{'overlap_ids'}};

	#would need to add stuff here to account for the twin situation
	#would need next_child, rather than previous child
	#more

	
	#reset child for current transition
	if($child == 0 && 
	   ($transitions->[$i]->{'current_child'} == $#{$transitions->[$i]->{'overlap_features'}})){
	  #children clash
	  $child =  $transitions->[$i]->{'previous_child'}; 
	}
	else{#no clash
	  $transitions->[$i]->{'current_child'};
	}


	#INHERITANCE
	if ($transitions->[$i]->{'source'} eq 'NEW'){
	  $child_id = $transitions->[$i]->{'overlap_features'}->[$child]->stable_id();
	  $feature_id = $transitions->[$i]->{'feature'}->dbID();

	  #This should be true as will be using previous if last transition projected to this feature
	  if(exists $dbid_mappings{$feature_id}){
		throw("Found duplicate inheritance:\tdbID ".$feature_id." > ".
			  $dbid_mappings{$feature_id}." & $child_id");
	  }
	  
	  
	  #only one per transition here as we're only dealing with one new feature
	  $mappings++;
	  
	  #Assign mapping
	  $feature->stable_id($child_id);
	  $dbid_mappings{$feature_id} = $feature;
	  	  	  
	  
	  #record all stable_id splits
	  foreach my $orphan_cnt(0..$last_orphan){#these are stable_ids!
		$orphan = $transitions->[$i]->{'overlap_ids'}->[$orphan_cnt];
		$orphan_id = $orphan->stable_id();
		
		#Rules 
		#we build a new mapping cache entry for each old stable id
		#mapping it to the new child_id
		#if this is 3 prime then we need to check if old_stable_id entry exists
		#and unshift
		#else can only be one known new stable_id(this feature) 
		
		#consituent hash is not recorded, but can be infered from split cache
		#e.g. new_stable_id [old_stable_ids]
		#mapping cache is
		#old_stable_id [[old_stable_id] new_stable_ids]
				  
		#can't look at dbID cache here, as that is specifically for new dbIDs

		if($orphan_cnt == $last_orphan && exists $mapping_cache{$orphan_id}){
			#This must be a linking overlap
			#we either need to push or unshift depending on child
			unshift @{$mapping_cache{$orphan_id}}, $child_id;
		}
		else{

		  #logic sanity check  REMOVE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		  if(exists $mapping_cache{$orphan_id}){
			throw("Have found duplicate mappings for $orphan_id\n");
		  }

		  #is either entirely contained
		  #or is 3' with no linking overlap
		  #or is 5' with an unknown overlap
		  #i.e. has only one known new stable_id(this feature)
		   @{$mapping_cache{$orphan_id}}, $child_id;
		}
	  }

	  #set last old stable_id here enables updating of split cache in next transition
	  $last_stable_id = $orphan_id;
	}
	else{#PROJECTION i.e. source == OLD

	  #this should be dbID!!

	  $child_id = $transitions->[$i]->{'overlap_features'}->[$child]->dbID();
	  $feature_id = $transitions->[$i]->{'feature'}->stable_id();


	  #my $outline = $mapping_info{'transitions'}[$i]{'source_id'};	  #stable_id
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
	  
	  
	  $mappings++;
	  #oxxx  ---a----- xxxxxxxxx
	  #n  -b1- -b2- -b3--- xxx xx
	  
	  #a>b2 (b1, b3)
	  
	  
	  #this will not fail here
	  #child will never be b3(which would already have a stable_id) as we will have specified use_previous
	  #do we need to revise this in light of encountering a twin?	  
	  if(exists $dbid_mappings{$child_id}){
		throw("Found duplicate projection:\tdbID $child_id > ".$feature_id.
			  " & ".$dbid_mappings{$child_id});
	  }
	  
	  #record dbID mapping first only done to check for duplicate inheritance
	  #$dbid_mappings{$child_id} = $mapping_info{'transitions'}[$i]{'source_id'};
	  #can't do this here as we're testing for it in the loop
	  
	  #then record splits?
	  #can not record split until we have dealt with new stable id generation
	  #and potential inherited stable_id from x>b1
		
	  #can we sub this??
				
	  #what if b1 inherits here due t use_previous
	  #oxxx  ---a--------- xxxx
	  #n  -b1--- -b2- -b3------
	  #This will always be added as we include the child_id first!!
	  
	  foreach my $orphan_cnt(0..$last_orphan){
		$orphan = $transitions->[$i]->{'overlap_features'}->[$orphan_cnt];
		$orphan_id = $orphan->dbID();#new dbID


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
	

		#don't change the order of this block!
	
		#3' rules
		if($orphan_cnt == $last_orphan){
				
		  if(exists $dbid_mappings{$orphan_id}){
			#already assigned from last transition(cannot be child)
			#This will have already been added to the previous split cache as it will be the previous child id
			$new_sid = $dbid_mappings{$orphan_id}->stable_id();
			#added to split cache below
		  }
		  elsif($orphan_id == $child_id){
			#last orphan is child
			$orphan->stable_id($feature_id);
			$dbid_mappings{$orphan_id} = $orphan;
			$new_sid = $dbid_mappings{$orphan_id}->stable_id();
					
			#we need add this to the previous split_merge_cache
			if(defined $last_stable_id){
			  
			  if(exists $mapping_cache{$last_stable_id}){
				push @{$mapping_cache{$last_stable_id}}, $dbid_mappings{$orphan_id}->stable_id();
			  }
			  else{
				throw("Failed to find split_merge_cache for old $last_stable_id when updating with linking new feature ".
					  $dbid_mappings{$orphan_id}->stable_id());
			  }
			}
		  }
		  else{
			#can only be last in chain
			$new_sid = &assign_and_log_new_stable_id($orhpan_id);

			#But may have been a split link i.e. a non-maximal lining overlap
			#so check if we have to add it to the previous split mappings cache

			if(exists $mapping_cache{$last_stable_id}){
			  push @{$mapping_cache{$last_stable_id}}, $dbid_mappings{$orphan_id}->stable_id();
			}
			else{
			  throw("Failed to find mapping cache for old $last_stable_id when updating with non-maximal linking new feature ".
					$dbid_mappings{$orphan_id}->stable_id());
			}
		  }		
		}
		elsif($orphan_id == $child_id){
		  #this will also account for the child being the 5' orphan
		  $orphan->stable_id($feature_id);
		  $dbid_mappings{$orphan_id} = $orphan;
		  $new_sid = $dbid_mappings{$orphan_id}->stable_id();;
		  #don't need to add this to the split_merge_cache as it is not overlapping


		  #If this i s5'
		  #we can't assign this yet until we do the previous transition!!!!
		  
		  #o   ---a----  --c--
		  #nxxxxxxx ---b--  -d-

		  #working backwards
		  #c > b & d
		  #b is child but we don't know whether previous transition will be better match!
		  
		  #b > a & c
		  #a is child but again we can't assign as there might be a better match
		  #This is the same shift prolem but in reverse!

		  #but we already know the details for the next transition
		  #we would have to traverse the transitions until we found which way to shift
		  #can we just reverse the use_previous logic if we find a non-maximal link?
		  #or should we just ignore this, split the chain and make sure we assign b to c in the split merge cache
		  #Yes!!!!!!
		  #split merge cache has to persist between chains. i.e. instances of this sub



		}
		elsif($i > 0 || ($i == 0 && $orphan_cnt == 0)){
		  #not child and contained orphan or 5' orphan at start of chain
		  $new_sid = &assign_and_log_new_stable_id($orphan_id);
		}
		else{#5' orphan in the middle of a chain
		  #don't assign stable id yet, update mapping cache in next transition
		  next;	#to ensure we populate this transition's mapping cache
		}
		
   		if ($orphan_id == $child_id){
		  unshift @new_stable_ids, $new_sid;
		}else{
		  push @new_stable_ids, $new_sid;
		}

		#dump this file in the correct tab txt format to enable file import
		#rather than updating every single record individually
		#this would require passing the features rather than just their dbIDs
		#need to create new feature_set at start of process and use new feature_set id
		#we're going to migrate reg_feats to another table, so do we want to do this
		#we need a way of defining the old and new set within the reg_feature table
		#have build field?
		#adaptor would default to fetching all from the current schema version specified in meta.
		#optional argument for schema
		#so we would have 3 sets
		#old                       e.g. 46
		#new with no stable IDs    e.g. 47a
		#new with stable           e.g. 47final
		#would then delete 47a if flag specified and update table to shift dbIDs to start from last +1 dbID?
		#would need to add this to update_db_for release script, which is essentially complex healthcheck script
		#rather than simply mysql healthchecks


		#print $mapping_file $orphan."\t".$new_sid."\n";
		



		#print unresolved 5' orphans in other block

	  }
	  
	  #set the cache, may not contain the most 5' stable_id which may need to be assigned in the next transition
	  @{$mapping_cache{$feature_id}} = @new_stable_ids;
	}
  }
  

  #so we should print out split_merge_cache here
  #we should also print out dbID > stable id mappings too, in loop or here, and do we need other constributing/merged stable_ids in this file
  #lastly we need to print out retired stable ids? shouldn't this go in the split merge cache?

}

=head2 assign_and_log_new_stable_id

  Arg [1]    : Bio::EnsEMBL:Funcgen::AnnotatedFeature
  Example    : my $new_sid = &assign_and_log_new_stable_id($new_reg_feature);
  Description: Creates and assigns a new stable ID to a given RegulatoryFeature.
               New dbID to stable ID mappings are logged in the new_ids file.
  Returntype : string - Stable ID
  Exceptions : Throws if no valid AnnotatedFeature is passed
  Status     : At Risk - Use Bio::EnsEMBL::Funcgen::RegulatoryFeature

=cut


sub assign_and_log_new_stable_id{
  my $new_reg_feat = shift;

  if(! (defined $new_reg_feat && $new_reg_feat->isa('Bio::EnsEMBL::Funcgen::AnnotatedFeature'))){
	throw('You must provide a valid Bio::EnsEMBL::Funcgen::AnnotatedFeature');
  }

  #Need to add species intrafix here

  my $new_sid = sprintf("ENSR%011d", $next_stable_id);
  $new_reg_feat->stable_id($new_sid);
  $dbid_mappings{$new_reg_feat->dbID()} = $new_reg_feat;
  $new_mappings ++;
  $next_stable_id ++;

  print $new_ids_handle $new_reg_feat->dbID()."\t".$new_sid."\n";

  return $new_sid;
}

sub dump_new_stable_feature{
  my ($new_reg_feat, $stable_id) = @_;



}


#implement this?
#sub add_to_split_merge_cache{
#  my ($old_stable_id, $new_stable_id) = @_;

#  if(exists $split_merge_cache{$last_stable_id}){
#	push @{$split_merge_cache{$last_stable_id}}, $dbid_mappings{$child_id};
#  }
#  else{
#	throw("Failed to find split_merge_cache for old $last_stable_id when updating with linking new feature ".
#		  $dbid_mappings{$child_id});
#  }
  
#  return;
#}

1;
