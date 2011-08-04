#!/usr/bin/env perl

=head1 LICENSE


  Copyright (c) 1999-2011 The European Bioinformatics Institute and
  Genome Research Limited.  All rights reserved.

  This software is distributed under a modified Apache license.
  For license details, please see

    http://www.ensembl.org/info/about/code_licence.html

=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <ensembl-dev@ebi.ac.uk>.

  Questions may also be sent to the Ensembl help desk at
  <helpdesk@ensembl.org>.

=head1 NAME

ensembl-efg stable_id_mapper.pl
  
=head1 SYNOPSIS

stable_id_mapper.pl [options]

Options:

Mandatory


Optional


=head1 OPTIONS

			"odbpass=s"          => \$opass,
			"odbport=s"          => \$oport,
			"ohost=s"            => \$ohost,
			"ouser=s"            => \$ouser,
			"odbname=s"          => \$odbname,
			"npass=s"          => \$npass,
			"nport=s"          => \$nport,
			"nhost=s"            => \$nhost,
			"nuser=s"            => \$nuser,
			"ndbname=s"          => \$ndbname,
			"dnadb_host=s"       => \$dnadb_host,
			"dnadb_pass=s"       => \$dnadb_pass,
			"dnadb_port=s"       => \$dnadb_port,
 			"dnadb_user=s"       => \$dnadb_user,
			"dnadb_name=s"       => \$dnadb_name,

			"old_fset_name=s"    => \$old_fset_name,
			"new_fset_name=s"    => \$new_fset_name,
			'old_assembly=s'     => \$old_assm,
			'new_assembly=s'     => \$new_assm,
			'coord_system=s'     => \$cs_level,
			"species=s"          => \$species,
			"expand|e=s"         => \$expand,
			'no_load'            => \$no_load,
			'clobber'            => \$clobber,
            'recover'            => \$recover,

 -slice_name        The name of a slice to run. This must not be run in parallel as stable IDs are 
                    incremented within this script, not by MySQL. Running in parallel will result in
                    duplicate new stable IDs, re-run with clobber to correct this.

 -assign_all_nulls  If running in slice mode with a new assembly, this ensures that all remaining 
                    RegualtoryFeatures which have not been processed will be assigned a stable ID.  
                    This should only be used when processing the last slice. Turning this on before
                    processing the final slice will result in all remaning features being assigned 
                    new rather than mapped stable IDs(this can be fixed by running the problem 
                    slices with -clobber, but will obfuscate which seq_regions have been 
                    successfully mapped).

	"stable_id=s"        => \$stable_id,
			"out_dir=s"           => \$out_dir,
			#"update"             => \$update,
			#"from_file|f=s"      =>\$from_file,

			'log_file=s'        => \$main::_log_file,
			'tee'               => \$main::_tee,
			

			"help|?"             => \$help,
			"man|m"              => \$man,


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


#To do

#1 Run on toplevel non-ref, now added, but what impact will this have?
#2 remove empty files
#3 Implement HealthChecker/Helper for logging and check_stable_ids sub
#4 Implement assembly mapping
#5 Implement recover properly, cannot do for sub slices due to chaining
#  Also mapped feature count may not always match old feature cnt
#  Currently just redoes any slice with mismatched counts
#  Mapped count is based on whole slice not just feature slice
#  Should get same results due to chaining, but need to check
#6 Verbose level logging?
#7 Blast based mapping for none assembly mapped regions?
#8 Do we rollback new stable IDs when rerunning i.e. are we reusing new stable IDs?
#9 Add mapping reason to output
#10 Copy stable_ids to other feature_sets no ending in _v[0-9]+

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
use Bio::EnsEMBL::Funcgen::Utils::Helper;
use Bio::EnsEMBL::Funcgen::Utils::EFGUtils qw (open_file run_system_cmd backup_file);
use Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Funcgen::DBSQL::RegulatoryFeatureAdaptor;
use strict;

$| = 1;							#autoflush
my ($dbname, $help, $man, @slice_names, @skip_slices, $clobber, $no_load, $odb, $recover, $no_dump);
my ($odbname, $ndbname, $npass, $nuser, $ohost, $oport, $from_file, $stable_id, $next_stable_id, $assign_nulls);
my ($dnadb_name, $dnadb_pass, $dnadb_user, $dnadb_host, $dnadb_port, $old_assm, $new_assm, $species);
my $reg = "Bio::EnsEMBL::Registry";

my $ouser  = 'ensro';
my $opass = '';
my $nhost = 'ens-genomics1';
my $old_fset_name = 'RegulatoryFeatures';
my $new_fset_name = 'newRegulatoryFeatures';
my $nport = '3306';
my $expand = 0;
my $out_dir = '.';
my $cs_level = 'chromosome';


$main::_debug_level = 0;
$main::_tee = 0;



print "Args are:\t@ARGV\n";

GetOptions (
			"odbpass=s"          => \$opass,
			"odbport=s"          => \$oport,
			"ohost=s"            => \$ohost,
			"ouser=s"            => \$ouser,
			"odbname=s"          => \$odbname,
			"npass=s"          => \$npass,
			"nport=s"          => \$nport,
			"nhost=s"            => \$nhost,
			"nuser=s"            => \$nuser,
			"ndbname=s"          => \$ndbname,
			"dnadb_host=s"       => \$dnadb_host,
			"dnadb_pass=s"       => \$dnadb_pass,
			"dnadb_port=s"       => \$dnadb_port,
 			"dnadb_user=s"       => \$dnadb_user,
			"dnadb_name=s"       => \$dnadb_name,

			"old_fset_name=s"    => \$old_fset_name,
			"new_fset_name=s"    => \$new_fset_name,
			'old_assembly=s'     => \$old_assm,
			'new_assembly=s'     => \$new_assm,
			'coord_system=s'     => \$cs_level,
			"species=s"          => \$species,
			"expand|e=s"         => \$expand,
			'no_load'            => \$no_load,
			'no_dump'            => \$no_dump,
			'clobber'            => \$clobber,
            'recover'            => \$recover,

			#These modes will give different results
			#as mapping is conxtext dependent i.e. what is immediately before
			#this slice/stable_id may or may not have an effect on the mapping.
			
			'assign_all_nulls'   => \$assign_nulls,
			'slices=s{,}'        => \@slice_names,
			'skip_slices=s{,}'   => \@skip_slices,
			'next_stable_id=i'   => \$next_stable_id,
			"stable_id=s"        => \$stable_id,
			"out_dir=s"          => \$out_dir,
			#"update"            => \$update,
			#"from_file|f=s"      =>\$from_file,

			'log_file=s'         => \$main::_log_file,
			'tee'                => \$main::_tee,
            'debug_level=i'      => \$main::_debug_level,			

			"help|?"             => \$help,
			"man|m"              => \$man,
			#add opt for old, new & stable fset name
		   ) || die('Specified invalid option, pod2usage needed here');


pod2usage(1) if $help;
pod2usage(-exitstatus => 0, -verbose => 2) if $man;


#check mandatory params here
die("You cannot specify load '-from_file' in conjunction with '-no_load'\n") if($no_load && $from_file);


#we need to resolve odb ndb relationship, defulat to same if only one specified?
#also need to add defaults for 

if(! ($nhost && $nuser && $npass && $ndbname)){
  die("You must provide new DB details: -ndbhost -ndbuser -ndbpass -ndbname\n");
}
#else{
#  print "Using new DB:\n\t".join("\n\t", ($nhost, $nuser, $npass, $ndbname))."\n";
#}


if((@slice_names || @skip_slices) &&
   ! $next_stable_id){
  die("To run with slice subsets you must specify the next_stable_id, as there is currently no way to tell whether this should come form the old or new DB");
}


#This should always be the handover db on staging
#Do we need to add old dnadb too?
#Too account for unversioned slices which are not present in the new DB?

my $cdb;

if($dnadb_name){
  $cdb =  Bio::EnsEMBL::DBSQL::DBAdaptor->new(
											  -host    => $dnadb_host || 'ens-staging',
											  -user    => $dnadb_user || 'ensro',
											  -pass    => $dnadb_pass || undef,
											  -port    => $dnadb_port || 3306,
											  -dbname  => $dnadb_name,
											  -species => $species,
											 );
}


my $ndb = Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor->new(
													   -host   => $nhost,
													   -user   => $nuser,
													   -pass   => $npass,
													   -port   => $nport,
													   -dbname => $ndbname,
													   -species => $species,
													   -dnadb => $cdb,
													  );

#Test connections, eval this?
$ndb->dbc;
$ndb->dnadb->dbc;

$main::_log_file ||= $ENV{'HOME'}."/logs/stable_id_mapper.$ndbname.$$.log";
print "Writing log to:\t".$main::_log_file."\n";
my $helper = new Bio::EnsEMBL::Funcgen::Utils::Helper;


if(! ($odbname)){
  #Should test for any old vars set in case we have forgotten to specify ohost
  $helper->log("Setting old DB to new DB");
  $odb = $ndb;
}
else{
  $odb = Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor->new(
													  -host   => $ohost || $nhost,
													  -user   => $ouser,
													  -pass   => $opass,
													  -port   => $oport,
													  -dbname => $odbname,
													  -species => $species,
													 );

  $odb->dbc();#Test connection, eval this?
}

my ($new_reg_feat_handle, $split_id);
my ($new_id_handle, $stable_id_handle, %dbid_mappings);
my (%mapping_cache, %obj_cache, $mappings, $new_mappings);#%new_id_cache?
my $total_stable_ids = 0;
my $total_new_stable_ids = 0;

#Do not change the order of this array unless you know what you're doing!!!!!
my @comparators = ('overlap_length', 'coverage', 'displacement');
my %comparatees;


#Define referenced methods here
#Comparator calculation and compare methods
#$tmp_comparatees{$comparator} = &$method($nstart, $nend, $midpoint, $nfstart, $nfend, $nflength);
#these have to be generated in this order as the are dependent on each other

sub calculate_overlap_length{
  #overlap length is always the difference between the two middle sorted start/end value pairs
  my (undef, $overlap_start, $overlap_end, undef) = sort{$a<=>$b} ($_[0], $_[1], $_[3], $_[4]);
  return ($overlap_end - $overlap_start +1);
}

sub calculate_coverage{
  return ($comparatees{'overlap_length'}{'tmp'}/ $_[4]);
}

sub calculate_displacement{
  #source feature midpoint - new feature midpoint *-1 id < 0
  my $disp = ($_[2] - ($_[3] + ( ($_[4] - $_[3]) /2 )));
  return ($disp < 0 ) ? ($disp * -1) : $disp;
}

sub is_greater_than{
  return ($_[0] > $_[1]);
}

sub is_less_than{
  return ($_[0] < $_[1]);
}


#is this too obfuscated?
#redundancy of subs here
my %comparison_methods = 
  (
   'overlap_length' => 
   {
	compare   => \&is_greater_than,
	calculate => \&calculate_overlap_length,
   },

   'coverage' => 
   { 
	compare   => \&is_greater_than,
	calculate => \&calculate_coverage,
   },

   'displacement' => 
   {
	compare   => \&is_less_than,
	calculate => \&calculate_displacement,
   },
  );




###Get Old and New FeatureSets and set adaptors
#get FeatureAdaptors from Fset, this will allow use of any feature type
$obj_cache{'OLD'}{'FSET'}          = $odb->get_FeatureSetAdaptor->fetch_by_name($old_fset_name);
die("Could not find OLD RegulatoryFeature FeatureSet:\t$old_fset_name") if ! $obj_cache{'OLD'}{'FSET'};
$obj_cache{'OLD'}{'SLICE_ADAPTOR'} = $odb->get_SliceAdaptor();
$obj_cache{'NEW'}{'SLICE_ADAPTOR'} = $ndb->get_SliceAdaptor();
$obj_cache{'NEW'}{'FSET'}          = $ndb->get_FeatureSetAdaptor->fetch_by_name($new_fset_name);
die("Could not find NEW RegulatoryFeature FeatureSet:\t$new_fset_name") if ! $obj_cache{'NEW'}{'FSET'};


$helper->log("Old FeatureSet:\t".$odb->dbc->dbname.':'.$old_fset_name);
$helper->log("New FeatureSet:\t".$ndb->dbc->dbname.':'.$new_fset_name);



#Now we need to test the CoordSystem versions if they have been specified
#This will require a mapping path.
#To use the project method on the old features the mapping path must be in the olddb
#i.e. in the olddb and the newdb must be the same
#This could be done without the project method, but is much easier like this.



###Need this to store new features?
my $nrf_adaptor = $ndb->get_RegulatoryFeatureAdaptor();
my $orf_adaptor = $odb->get_RegulatoryFeatureAdaptor();
##Set start stable_id
#This looks accross all reg feat sets to avoid using a previously retired stable ID
#Or duplicating if we are running in recovery mode and some have already been assigned for the new set


#This is dependent on clobber or recover.

warn "Need to implement clobber/recover in next_stable_id lookup\n";

my ($cmd, $on_sid, $nn_sid);

if($next_stable_id){#check it is valid
  $helper->log("Using specified next stable:\t$next_stable_id");
  #Test both dbs else fail
  #should we rollback the stable_ids before we do this?
  #Yes! Otherwise we may end up skipping lot's of stable_ids

}else{
  warn "WARNING:\tFetching next_stable_id from old DB\n";

  $cmd = 'SELECT stable_id from regulatory_feature rf where feature_set_id='.$obj_cache{'OLD'}{'FSET'}->dbID().' order by stable_id desc limit 1';
  ($next_stable_id) = @{$odb->dbc->db_handle->selectrow_arrayref($cmd)};
  $next_stable_id ++;
  $helper->log("Next stable id:\t$next_stable_id");
}


die('Cannot specify both a slice and a stable_id to test run on') if($stable_id && @slice_names);


#Check the old and new coord_systems

if($old_assm || $new_assm){

  if(! $old_assm || ! $new_assm){
	die('Must provide both a -new_assm and an -old_assm if you wish to stable id map between assemblies');
  }

  my $cs = $odb->dnadb->get_coordSystemAdaptor->fetch_by_name($cs_level, $old_assm);
	
  if(! $cs){
	die("Cannot access $cs_level $old_assm assembly from $dnadb_name");
  }

  $odb->dnadb->get_coordSystemAdaptor->fetch_by_name($cs_level, $new_assm);
  
  if(! $cs){
	die("Cannot access $cs_level $old_assm assembly from $dnadb_name");
  }

  $helper->log("Projecting stable IDs between $old_assm and $new_assm");

}




#Grab the old slices
my (@top_level_slices, @slices, %new_slices);

#Need to implement EFGUtils::generate_slices_from_names here

if(@slice_names){
  $helper->log("Running only on slices:\t".join(', ', @slice_names));
  warn "WARNING:\tThis script cannot be run in parallel as stable ID increment is performed by the script, not by MySQL\nRunning in parallel will result in duplicate new stable IDs\n";

  foreach my $slice_name(@slice_names){
	
	if($old_assm && $slice_name !~ /$old_assm/){
	  die("You must map from an old slice($old_assm) if you wish to map between assemblies:\t$slice_name");
	}


	#Try region first
	my $slice;

	
	$slice = $obj_cache{'OLD'}{'SLICE_ADAPTOR'}->fetch_by_region(undef, $slice_name);
	
	if(! $slice){
	  eval{ $slice = $obj_cache{'OLD'}{'SLICE_ADAPTOR'}->fetch_by_name($slice_name) };
	  warn "Caught Exception:\n$@" if $@;
	}
	
	if(! $slice){
	  warn("You have specified a slice name which cannot be found in the database:\t".$slice_name);
	}
	else{
	  push @slices, $slice;
	}
  }
}
elsif($stable_id){
  #This will automatically use the old assembly if the old fset only has features on the old assembly
  $no_dump = 1;
  @slices = ($odb->get_RegulatoryFeatureAdaptor->fetch_by_stable_id($stable_id, $obj_cache{'OLD'}{'FSET'})->feature_Slice());
  $helper->log("Running only on $stable_id ".$slices[0]->name);
}  
else{
  #fetch top_level, including non-ref
  #This could actually be the new toplevel if we have set the olddb to the newdb
  #But we're only using this for the unversioned slices, which will will have to be present in the new DB
  #to be able to map them. This approach may miss logging of stable IDs present on old/expired seqs
  #Hence why we should force usage of old dna DB when projecting between assemblies
  @top_level_slices = @{$obj_cache{'OLD'}{'SLICE_ADAPTOR'}->fetch_all('toplevel', undef, 1)};
  
  #Fetch old assembled level
  #What about > 1 assembled/versioned level?
  #These would not be caught below
  @slices = @{$obj_cache{'OLD'}{'SLICE_ADAPTOR'}->fetch_all($cs_level, $old_assm)};

  #Add all toplevel nonref slices
  #So we have the old toplevel approximately(new non versioned slices may have been added/removed)
  #No mapping available for slices which are not present in both DBs
  #Would have to do this with blast

  #Non-ref seqs may have version(haps, mouse NT supercontigs?)
  #Only reference seqs are assembly mapped i.e. chromosomes(not inc haps)
  #some supercontigs have assembly versions, but can assume they are the same if they 
  #have the same name and have the same length. Should really have a name change if they are different.
  
  
  foreach my $slice(@top_level_slices){

	#supercontigs can be reference!
	#But how can they be reference if they are not versioned?
	my $cs_name = $slice->coord_system->name;

	if(($cs_name ne 'chromosome') ||
	   ! $slice->is_reference){

	  #Move these test to loop below
	  #So we can log the unmapped stable IDs in the cache

	  if($old_assm){
		#Validate old seq is same as new seq or has mapping

		my $sr_name = $slice->seq_region_name;

		if($cs_name ne 'lrg'){
		  #lrgs never change assembly

		  my $new_slice = $obj_cache{'NEW'}{'SLICE_ADAPTOR'}->fetch_by_region($cs_name, $sr_name, undef, undef, undef, $old_assm);

		  if(! defined $new_slice){
			#Do we add this to the array and let the rest of the do the logging 
			#of the unmappable stable IDs?
			warn "Could not find new Slice for old seq_region:\t".$slice->name."\n";
			next;
		  }
		  elsif(($cs_name eq 'supercontig') ||
				($cs_name eq 'chromosome'))	{
			#Do we have a new supercontig/hap of same name & length?
			
			if($slice->length ne $new_slice->length){
			  warn "Found different lengths for ${cs_name}:${old_assm}:${sr_name} and ${cs_name}:${new_assm}:${sr_name}.\n".
				"Skipping mapping for ${cs_name}:${old_assm}:${sr_name} as no assembly projection available\n";
			  next;
			  #Should we next here or will this be handled by the code?
			  #i.e. can we include this so we get the logs
			}
		  }
		  else{
			warn "Cannot handle assembly differences/projection for ${cs_name}:${sr_name}";
			next;
		  }
		}
	  }
	  
	  push @slices, $slice;
	}
  }


  map $new_slices{$_->name} = $_, @{$obj_cache{'NEW'}{'SLICE_ADAPTOR'}->fetch_all('toplevel', undef, 1)};
}


if(! %new_slices){ #Running with subset of slices

  if($assign_nulls){

	#should really have force_assign_nulls here

	if(@slice_names){
	  die("It is currently unsafe to -assign_all_nulls when running with a subset of slices");
	}

	%new_slices  = map {$new_slices{$_->name} = $_} @{$obj_cache{'NEW'}{'SLICE_ADAPTOR'}->fetch_all('toplevel', undef, 1)};
  }
  else{

	foreach my $slice(@slices){
	  $new_slices{$slice->name} = $obj_cache{'NEW'}{'SLICE_ADAPTOR'}->fetch_by_region($slice->coord_system->name, $slice->seq_region_name);
	  #May no have a new slice? Deal with this later as it may have a mapping
	}
  }
}



$helper->log("Number of slices to map:\t".scalar(@slices));


#Need to explicitly build seq_region cache as were using the 'private' seq_region methods out of context
$nrf_adaptor->build_seq_region_cache();
$orf_adaptor->build_seq_region_cache();
my $total_old_feats    = 0;
my $total_failed_proj  = 0;
my $total_new_feats    = 0;
my $total_mapped_feats = 0;


#Set disconnect_when_inactive on old DB if not the same as new DB
#This should stop server timeouts from causing uncaught lost connections

#if($ndb ne $odb){
  #Test obj ref to ensure this is the same connection
$helper->log("Setting disconnect_when_inactive on old dnadb:\t".$odb->dnadb->dbc->dbname);
$odb->dnadb->dbc->disconnect_when_inactive(1);
$helper->log("Setting disconnect_when_inactive on new dnadb:\t".$ndb->dnadb->dbc->dbname);
$ndb->dnadb->dbc->disconnect_when_inactive(1);


#}

#Process each top level old seq_region




foreach my $slice (@slices){

  #Need to change this to slice_name
  my $seq_name = $slice->name();

  #Output files
  my $new_reg_feat_file = $out_dir.'/regulatory_feature.'.$seq_name.'.txt';
  #tab txt dump of new stable features

  my $new_id_file = $out_dir.'/new_stable_ids.'.$seq_name.'.txt';
  #This is slightly redundant? This will be present in the db_id_stable_id file 
  #and absent in the stable_id file
  #         new_id1                         #birth

  my $stable_id_file = $out_dir.'/stable_id_mappings.'.$seq_name.'.txt';
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

  my $sr_id = $orf_adaptor->get_seq_region_id_by_Slice($slice);

  if(! defined $sr_id){
	#This should never happen as we are using new Slices
	#remove after testing
	warn 'Slice '.$slice->name." does not exist in the old database.\n";
	next;
  }

  #This reliably exited with success, when an undef sr_id was used in this query

  $cmd = 'select count(regulatory_feature_id) from regulatory_feature where feature_set_id='.
	$obj_cache{'OLD'}{'FSET'}->dbID().' and seq_region_id='.$sr_id;

  my $failed_proj_cnt = 0;
  my ($old_feat_cnt)   = $odb->dbc->db_handle->selectrow_array($cmd);
  $total_old_feats   += $old_feat_cnt;
  
  #Need to test for seq_region_id here!
  $sr_id = $nrf_adaptor->get_seq_region_id_by_Slice($new_slices{$slice->name});

  if(! defined $sr_id){
	#This should never happen as we are using new Slices
	#remove after testing
	warn 'Slice '.$slice->name." does not exist in the new database.\n".
	  "This should not affect the mappings. Run update_DB_for_relase.pl to fix this\n";
  }


  $cmd = 'select count(regulatory_feature_id) from regulatory_feature where feature_set_id='.$obj_cache{'NEW'}{'FSET'}->dbID().
	' and seq_region_id='.$sr_id.' and stable_id is not NULL';
 
  my ($mapped_feature_cnt) = @{$ndb->dbc->db_handle->selectrow_arrayref($cmd)}; 
  $helper->log_header("Processing slice $seq_name with $mapped_feature_cnt mapped features from a total of $old_feat_cnt");

  
  if($mapped_feature_cnt != 0){
	
	if($mapped_feature_cnt == $old_feat_cnt){
	  $helper->log("All RegulatoryFeatures have already been stable ID mapped for $seq_name") 
		if($mapped_feature_cnt == $old_feat_cnt);

	  if($recover){
		die("Skipping $seq_name in recover mode\n");#??? surely warn?
		next;
	  }
	}

	if(! $clobber){
	  die("$mapped_feature_cnt/$old_feat_cnt RegulatoryFeatures have already been stable ID mapped for $seq_name, specify -clobber to overwrite");
	}
	elsif($no_load){
	  $helper->log("There are features which have already been stable ID mapped for $seq_name in the DB, you have chosen leave these and dump a summary of the new ones");
	}
  }
  

  #Start the stable ID mapping
  if(! $from_file){
	#open files
	warn 'Need to implement recovery here for backing up files';
	
	$new_reg_feat_handle  = open_file($new_reg_feat_file, '>');
	$new_id_handle = open_file($new_id_file, '>');
	
	#clean vars
	%mapping_cache = ();
	%dbid_mappings = ();
	$mappings = 0;
	$new_mappings = 0;
	my $proj_slice;

	$helper->log("Fetching old ".$obj_cache{'OLD'}{'FSET'}->name." features for:\t".$slice->name);
	my @old_reg_feats = @{$obj_cache{'OLD'}{'FSET'}->get_Features_by_Slice($slice)};
	$odb->dnadb->dbc->disconnect_if_idle();
	#warn "Connected status is now ".	$odb->dnadb->dbc->connected;
	#This should called automatically bby the fetch above, but does not seem to work


	#If new_assm, test if we can project or exists in new DB
	#if not we need to log the unmappable stable_ids
	#Can we project whole chromosomes like this?
	
	$helper->log("Processing ".scalar(@old_reg_feats)." RegulatoryFeatures for:\t\t\t".$slice->name);

  REGFEAT: foreach my $reg_feat(@old_reg_feats){
	  my $from_db = 'OLD';
	  my $feature = $reg_feat;
	  my $use_next_child = 0;
	  #Can we change this to an array as it's only containing transitions array and the orignal stable_id
	  #my $transitions = []; #this has to be ref due to implementation in while

	  if(! defined $reg_feat->stable_id){
		die("Found OLD RegulatoryFeature without a stable_id, dbID is ".$reg_feat->dbID);
	  }


	  $helper->debug(1, "STARTING CHAIN:\tOLD ".$reg_feat->stable_id().' '.$reg_feat->feature_Slice->name, 1);

	  my $source_dbID; #The source dbID this feature was fetched from, set to null for start feature.
	  #($stable_id = $reg_feat->display_label()) =~ s/\:.*//;
	  $stable_id = $reg_feat->stable_id;
	  next if exists $mapping_cache{$stable_id};
	  
	  #keep building transition whilst we have a 3' overhanging feature returned
	  while(defined $feature){		
		#we don't need to sub this part do we?
		#can we not have this all in one hash and just test for the feature hash element?

		#Assembly map feature first if required
		if($new_assm && $feature->slice->is_reference){
		  my $fail;

		  my @segments = @{$feature->project($cs_level, $new_assm)};
		  # do some sanity checks on the projection results:
		  # discard the projected feature if
		  #   1. it doesn't project at all (no segments returned)
		  #   2. the projection is fragmented (more than one segment)
		  #   3. the projection doesn't have the same length as the original
		  #      feature

		  # this tests for (1) and (2)
		  if (scalar(@segments) == 0) {
			$failed_proj_cnt++;#Is this correct cnt?
			#$no_projection_cnt++;
			$fail = 1;
		  } 
		  elsif (scalar(@segments) > 1) {
			$failed_proj_cnt++;
			#$multi_segment_cnt++;
			$fail = 1;
		  }
		  else{
    
			# test (3)
			$proj_slice = $segments[0]->to_Slice;
		  
			if ($feature->length != $proj_slice->length) {
			  #$wrong_length_cnt++;
			  #if(! $ignore_length)){
			  $failed_proj_cnt++;
			  $fail = 1;
			  #}
			}
		  }

		  if($fail){
			#$helper->log("No assembly mapping possible for:\t".$feature->stable_id());
			$mapping_cache{$feature->stable_id()} = [];
			#Do we need anythign else in here?
			next REGFEAT;
		  }
    
		  # everything looks fine, so adjust the coords of your feature
		  #Have to generate new_slice here as we are not sure it is going to be 
		  #on the same slice as the old assembly
		  my $new_full_slice = $obj_cache{'NEW'}{'SLICE_ADAPTOR'}->fetch_by_region($cs_level, $proj_slice->seq_region_name);
		  $feature->start($proj_slice->start);
		  $feature->end($proj_slice->end);
		  $feature->slice($new_full_slice);
		}

		($from_db, $feature, $source_dbID, $use_next_child, $split_id) = &build_transitions($from_db, $feature, $source_dbID, $use_next_child, $split_id);
	  }
	}

	$helper->log_header("Mapped $mappings old stable IDs from toplevel seq_region $seq_name");
	

	if($new_assm){
	  $total_failed_proj+=$failed_proj_cnt;
	  $helper->log("Failed to project ${failed_proj_cnt}/${old_feat_cnt}");
	}



	#Now for this slice get each new feature which hasn't already been mapped, generate a new stable id
	#my $new_slice = $obj_cache{'NEW'}{'SLICE_ADAPTOR'}->fetch_by_region('toplevel', $seq_name);
	#what was I going to do with this???????????????????????????????????????????????????????????????????????????????????????

						 
	$helper->log("Assigning new stable IDs");

	#Here we need to use the new assembly to fetch remaining feats
	#These may not match directly and we may have new seq_regions in a new assembly
	#So we need to use the new toplevel?
	#Can't do this unless we deal with all seq_regions at once
	#unless we know we are mapping the last slice and force?

	
	#This needs to use the new slice
	#If we are assembly mapping
	#May not have new seq if we have lost a medium level assembled seq_region i.e. supercontig etc.
	#We may also have extra new assembled slices which aren't represented in the old slices
	#This may be a case for assembly/stable_id mapping using blast

	my $new_slice;



	if($slice->coord_system->version){
	  #Default version will be new assembly
	  #or just the same asssembly if we are not mapping
	  $new_slice = $obj_cache{'NEW'}{'SLICE_ADAPTOR'}->fetch_by_region($slice->coord_system->name, $slice->seq_region_name);
	}
	else{
	  $new_slice = $slice;
	}

	#Remove slice from new_slice hash
	delete $new_slices{$seq_name};
	my $new_seq_name = $new_slice->name;
	

	my @nregs = @{$obj_cache{'NEW'}{'FSET'}->get_Features_by_Slice($new_slice)};
	my $new_feat_cnt = scalar(@nregs);
  
	#&print_problem_stable_ID;

	foreach my $nreg_feat(@nregs){
	  
	  if(! exists $dbid_mappings{$nreg_feat->dbID()}){
		#New regfeat has not been mapped to, so generate new stable_id
		&assign_and_log_new_stable_id($nreg_feat);		
	  }
	}

	#&print_problem_stable_ID;
	
	$helper->log('Stable IDs mapped for seq_region '.$seq_name.":\t${mappings}/${old_feat_cnt}");
	$helper->log('New stable IDs for seq_region '.$new_seq_name.":\t${new_mappings}/${new_feat_cnt}");
	$total_new_feats      += $new_feat_cnt;
	$total_new_stable_ids += $new_mappings;
	$total_mapped_feats   += $mappings;
	$total_stable_ids     += ($mappings + $new_mappings);	



	#only do this for testing and comparison?
	#as it will be difficult to load regulatory feature from file once we have implemented the attributes table
	#or we could map dbID undef and set the new feature set for all the features and store directly?
	my $action = ($no_load) ? 'Dumping' : 'Updating';

	if( ($action eq 'Dumping') && 
		($no_dump) ){
	  $action = undef;
	};


	if(! defined $action){
	  $helper->log('Skipping dump/load RegulatoryFeatures for Slice '.$seq_name);
	}
	else{
	  $helper->log($action.' RegulatoryFeatures for Slice '.$seq_name);
	  
	  my $sql_file    =  $out_dir.'/stable_id_patch.'.$seq_name.'.sql';
	  my $update_file = open_file($sql_file, '>');
	  my @io_buffer;
	  
	  
	  #Dump update sql to file
	  foreach my $f(values %dbid_mappings){
		
		
		if($no_load){  #dump a summary.. Remove this?
		  print $new_reg_feat_handle join("\t", 
										  ($f->dbID, $f->slice->get_seq_region_id, $f->seq_region_start, 
										   $f->seq_region_end, $f->seq_region_strand, 
										   $f->display_label, $obj_cache{'NEW'}{'FSET'}->dbID, $f->stable_id))."\n";
		}
		
		#Strip prefix
		(my $sid = $f->{'stable_id'}) =~ s/ENS[A-Z]*0*//;

		#Appear to be much faster to dump the sql to file and run on cmdline
		#Also have patch log of mapping
		#DB connections(with discnnect_when_inactive_set) versus perl IO & cmdline patch
		#$cmd = 'UPDATE regulatory_feature set stable_id='.$sid.' where regulatory_feature_id='.$f->dbID();
		#$ndb->dbc->db_handle->do($cmd);


		push @io_buffer, 'UPDATE regulatory_feature set stable_id='.$sid.' where regulatory_feature_id='.$f->dbID().";\n";
		#Can we do this directly using a print join map
		#rather than using an io_buffer?
		#io_buffer will take less memory but more IO

		if(scalar(@io_buffer) >5000){
		  print $update_file join('', @io_buffer);
		  @io_buffer = ();
		}
	  }

	  #Flush IO buffer to update sql file
	  if(scalar(@io_buffer) >0){
		print $update_file join('', @io_buffer);
		@io_buffer = ();
	  }
	  close($update_file);


	  #Update stable IDs via cmdline
	  if(! $no_load){
		$helper->log($action.' RegulatoryFeatures for Slice '.$seq_name);
		my $cmd = "mysql -u${nuser} -P${nport} -h${nhost} -p${npass} ${ndbname}<".$sql_file;
		run_system_cmd($cmd);#will exit on fail
	  }


	  #Print stable ID mappings
	  $stable_id_handle = open_file($stable_id_file, '>');
	  
	  foreach my $old_stable_id(keys %mapping_cache){
		
		push @io_buffer, $old_stable_id."\t".(join "\t", @{$mapping_cache{$old_stable_id}})."\n";

		if(scalar(@io_buffer) >5000){
		  print $stable_id_handle join('', @io_buffer);
		  @io_buffer = ();
		}
	  }
	  
	  #Flush IO buffer to stable ID mapping file
	  if(scalar(@io_buffer) >0){
		print $stable_id_handle join('', @io_buffer);
		@io_buffer = ();
	  }
	  close($stable_id_handle);
	  
	  #should we cache other data and print here too?
	}
  }

  $helper->log("Finished mapping stable IDs for:\t".$slice->name."\n");
}



  #}
  
  #&print_problem_stable_ID;


  #This is supposed to be the block which populates the stable_id_mapping_event tables
  #Currently not supported
  #if(! $no_load){
  #
  #foreach my $old_stable_id(keys %mapping_cache){
  #	  warn "Stable ID event import not yet implemented\n";
  #	  last;
  #	  #load each mapping event here
  #	}
  #  }
#}



#sub print_problem_stable_ID{
#  my $tmp = $ndb->dbc->db_handle->selectrow_arrayref('select "DBID STABLE_ID", regulatory_feature_id, stable_id from regulatory_feature where regulatory_feature_id=4987772');
#  $helper->log(Data::Dumper::Dumper($tmp));
#}


#This no longer checks for NULL for all as we do that above
#

$new_mappings = 0;
my $new_slice_cnt = 0;

foreach my $slice_name(keys %new_slices){

  if(! exists $new_slices{$slice_name}){
	#This can only happen if we are mapping an old slice using -slice_name
	#Which does not appear in the new assembly
	warn "The slice you have specified does not appear in the new assembly:\t".$slice_name.
	  "\nMaybe we need to implement blast based mapping?";
  }
  else{


	#This should now only happen if we are assembly mapping and 
	#we have a new seq_region in the new assembly
	#This will happen when running without -slice_name
	#or with -slice_name and -assign_all_nulls
	$new_slice_cnt++;

	my (undef, undef, $sr_name, $sr_start, $sr_end) = split/:/, $slice_name;

	$helper->log("Found a $new_assm slice which was not present in $old_assm:\t".$slice_name);

	#This is only counting the core regions, will this be different from the mapping process?
	
	my $sql = "select count(stable_id) from regulatory_feature rf, seq_region sr where rf.seq_region_start <= $sr_start and rf.seq_region_end >= $sr_end and feature_set_id=".$obj_cache{'NEW'}{'FSET'}->dbID.' and stable_id is NULL';
	my ($null_count) = $ndb->dbc->db_handle->selectrow_array($sql);
	
	
	if($null_count){#>1
	  $helper->log("Slice $slice_name still has $null_count RegulatoryFeatures without a stable ID assignment");
	  
	  if($assign_nulls){
		warn "You are running with -assign_all_null and -slice_name.  This will assign new stable IDs to all other unmapped RegualtoryFeatures. If this has been done in error, you need to rerun all or each unmapped slice with the -clobber option\n";
	  }
	  
	  foreach my $nreg_feat(@{$obj_cache{'NEW'}{'FSET'}->get_Features_by_Slice($new_slices{$slice_name})}){
		
		
		#Still have to check as we may have had some project???
		#Not sure this is true but leave anyway for safety
		
		if(! exists $dbid_mappings{$nreg_feat->dbID()}){
		  #New regfeat has not been mapped to, so generate new stable_id
		  &assign_and_log_new_stable_id($nreg_feat);		
		}
	  }
	}
  }
}

$total_new_stable_ids += $new_mappings;

if( ! $from_file){
  $helper->log_header("Total slices mapped(inc new & old):\t".(scalar(@slices) + $new_slice_cnt));
  $helper->log("Total failed to project:\t\t\t${total_failed_proj}/${total_old_feats}") if $new_assm;
  $helper->log("Total mapped Stable IDs:\t\t\t${total_mapped_feats}/${total_old_feats}");
  $helper->log("Total new stable IDs:\t\t\t${total_new_stable_ids}/${total_new_feats}");
  $helper->log("Total stable IDs(inc new & mapped):\t$total_stable_ids");
}



=head2 build_transitions

  Arg [1]    : string - 'OLD' or 'NEW' denoting the DB source of the current feature 
  Arg [2]    : Bio::EnsEMBL::Funcgen::RegulatoryFeature
  Arg [3]    : int - The dbID of the source feature from which this feature was retrieved by overlap
  Example    : 	($from_db, $feature, $source_dbID) = &build_transitions($from_db, $feature, $source_dbID);
  Description: Traverses mapping transitions in reverse order assigning
               or creating new stable IDs and caching stable ID history 
  Returntype : None
  Exceptions : Throws if duplicate mappings are found
  Status     : At Risk - incorporate into loop above.

=cut

sub build_transitions{
  my ($from_db, $feature, $source_dbID, $use_next_child, $split_id) = @_;
  #the source_dbID is the dbID of the previous source feature
  #to ensure we don't create a recursive transition, hence enabling detection of the end transition
  
  #Will this join avoid undef warning?
  #print "Building transitions:\t".join(', ', ($from_db, $feature, $source_dbID, $use_next_child, $split_id))."\n";

  my ($o_start, $o_end, $previous_child, $current_child, $displacement, $num_parents, $coverage);
  my ($next_child, @transitions, @overlap_ids);
  #my $overlap = 0;#length of overlap
  $comparatees{'overlap_length'}{'child'} = 0;
  $comparatees{'overlap_length'}{'next_child'} = 0;
    

  my $no_new_features = 0;
  my $to_db = ($from_db eq 'OLD') ? 'NEW' : 'OLD';
  my $end_chain = 0;
  my $split_chain = 0;
  my $overhang = 0;

  #disabled default extend for now for simplicity
  # we should only do this if we have no mapping
  #if($expand && $mapping_info->{'iteration'} == 0){
  #	$mapping_info->{'extended'} = 1;
  #	$slice = $slice->expand($expand, $expand);
  #  }
   

  $helper->debug(1, "BUILDING TRANSITIONS:\tFetching $from_db Feature(".$feature->dbID.") as $to_db Slice\t".$feature->seq_region_name.':'.$feature->seq_region_start.':'.$feature->seq_region_end);

  #can we replace this with feature_Slice?
  my $next_slice = $obj_cache{$to_db}{'SLICE_ADAPTOR'}->fetch_by_region('toplevel', $feature->seq_region_name(), $feature->seq_region_start(), $feature->seq_region_end());
  

  #if(! $next_slice){#I don't think this should ever happen, can we remove?
  #print $reg_feat->display_label()." is on a slice which is not present in the $to_db DB:\t".
  #	  $slice->seq_region_name().' '.$slice->start().' '.$slice->end."\n";
  #	next;
  #  }
    

  #this needs to be replaced with some sort of intelligent sub which implements syntenic/seq comparisons
  #if we can't map dynamically using the assemlby mapper
  #this would also implement a prioritised code array of extension rules?

  my @nreg_feats = @{$obj_cache{$to_db}{'FSET'}->get_Features_by_Slice($next_slice)};
  my $id = ($from_db eq 'OLD') ? $feature->stable_id() : $feature->dbID();#Is this safe?
  
  #Surely this is always OLD to NEW?
  #print "\nMapping $from_db $id  to ".scalar(@nreg_feats)." ($to_db) RegFeats\n";
  #print "split_id is $split_id" if $split_id;
  #print "use next_child is $use_next_child" if $use_next_child;

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
	  $helper->debug(1, "No $to_db mapping possible for:\t".$feature->stable_id());
	  $mapping_cache{$feature->stable_id()} = [];
	  #next;					#there should only be one pair in this loop, so will exit loop
	}
	else{
	  die("this should never happen!");#REMOVE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	}
  }
  else{
	#set a few local vars to avoid accessors and redundancy
	my $fstart =  $feature->seq_region_start();
	my $fend =  $feature->seq_region_end();
	my $midpoint =  $fstart + (($fend -  $fstart)/2);#this is different to slice->centrepoint

	
	#Loop through features assigning best overlap as child
	#recording next best as previous child(?or next best in the case of twins?)
	
	for my $i (0..$#nreg_feats) {
	  #if we're in the middle/end of a chain, the first will always be the source_dbID/dbID
	  
	  
	  if (defined $source_dbID && $#nreg_feats == 0) { #only found source i.e. end of the chain
		#we could sanity check here that the source_dbID == the $nreg_feats[0] dbID
		$no_new_features = 1;
		last;
	  }
	  
	  $id = (defined $nreg_feats[$i]->stable_id()) ? $nreg_feats[$i]->stable_id() : '';
	  
	  $helper->debug(2, "FOUND $to_db:\t".$nreg_feats[$i]->dbID().":".$id, 1);
	  
	  #set a few more local vars
	  my $nfstart = $nreg_feats[$i]->seq_region_start();
	  my $nfend = $nreg_feats[$i]->seq_region_end();
	  my $nflength = $nreg_feats[$i]->length();
	  my %new_children = (
						  child => 0,
						  next_child => 0,
						 );
						 

	  #overlap length is always the difference between the two middle sorted values
	  #(undef, $o_start, $o_end, undef) = sort($feature->seq_region_start(), $feature->seq_region_end(), $nreg_feats[$i]->seq_region_start(), $nreg_feats[$i]->seq_region_end());
	  #my $tmp_overlap = $o_end - $o_start;
	  #my $tmp_displacement = $midpoint - ($nreg_feats[$i]->seq_region_end() -  $nreg_feats[$i]->seq_region_start());
	  #my $tmp_coverage = $tmp_overlap / $nreg_feats[$i]->length();
	  
	  
	  
	  #Calculate all comparatees
	  foreach my $comparator(@comparators){
		#warn "Calculating $comparator comparator from ($fstart, $fend, $midpoint, $nfstart, $nfend, $nflength)";

		$comparatees{$comparator}{'tmp'} = $comparison_methods{$comparator}->{'calculate'}->
		  ($fstart, $fend, $midpoint, $nfstart, $nfend, $nflength);
	  }
	  
	  #Loop through comparators child first then next child

	OFF_SPRING: foreach my $off_spring('child', 'next_child'){
		#Only compare new features for next child if they don't qualify for new_child
		last if $new_children{'child'} && $off_spring eq 'next_child';

		foreach my $comparator(@comparators){
		
		  $helper->debug(3, "Testing new $comparator(".$comparatees{$comparator}{'tmp'}.
						 ") vs $off_spring(".$comparatees{$comparator}{$off_spring}.")");


		  if($comparison_methods{$comparator}->{'compare'}->
			 ($comparatees{$comparator}{'tmp'}, $comparatees{$comparator}{$off_spring})){
			
			$helper->debug(2, "Found new $comparator(".$comparatees{$comparator}{'tmp'}.
						   ") vs $off_spring(".$comparatees{$comparator}{$off_spring}.")");
			
			$new_children{$off_spring} = 1;
		  }
		  elsif($comparatees{$comparator}{'tmp'} ==  $comparatees{$comparator}{$off_spring}){
			
			$helper->debug(3, "Found new equal $comparator(".$comparatees{$comparator}{'tmp'}.
						   ") vs $off_spring(".$comparatees{$comparator}{$off_spring}.")");

			if($comparator eq 'displacement'){#must be ==
			  #only required as comparators are not exhaustive
			  
			  warn "Displacement is failing to resolve stable_id mapping, need num parent rules!! Defaulting to 5' mapping";
			  $new_children{$off_spring} = 1;
			}
			
			#only call next if comparators are equal for both features
			next;
		  }
		
		  #default to exit comparisons as we know if new comparatee is
		  #a winner or not
		  last;#comparator
		}
	  }

	  #left for docs on symmetrical inheritance/twins etc
	  #	#Length test
	  #	if ($tmp_overlap > $overlap) {
	  #	  $new_child = 1;
	  #	} 
	  #	elsif ($tmp_overlap == $overlap) {
	  #
	  #	  #Coverage test
	  #	  if ($tmp_coverage > $coverage) {
	  #		$new_child = 1;
	  #	  } elsif ($tmp_coverage == $coverage) {
	  #		#warn "Coverage is not resolving stable_id inheritance, need more rules, defaulting to 3' inheritance for \n";
	  #		#$new_child = 1;			#default to 3' inheritance!!!!!
	  #		#We need to resolve this, or just opt for >= in above statement for now
	  #		#this would give 3' presidence to inheritance
	  #		
	  #		#?????
	  #		#would also need to acount for equal %age overlap and look at number of overlaps?
	  #		#what would happen here?
	  #		#o   -------  -------
	  #		#n    --  ------  -
	  #		#what if this does not resolve?
	  #		#do we need another rule?
	  #		#o------------
	  #		#n    --    --
	  #		#arguably the one in the middle should inherit
	  #
	  #
	  #		#will we need to store coverage etc for resolving complex inheritance
	  #		#in the reverse assign loop?
	  #		#maybe not coverage, but maybe other info derive from the last transition 
	  #		#which will not be available from the linking overlap, which will be available in the next transition
	  #		#i.e. num of parents?
	  #		#o -a c1 --c2----
	  #		#n  --b-----  --d-----
	  #		#Here b should in herit from c and not d
	  #		#as their coverage is equal but b has other parental support
	  #		#what about?
	  #		#o -a c1  --c2--
	  #		#n  --b-----  --d-----
	  #		#b would inherit from c1 due to coverage
	  #		#and this?
	  #		#o  -a c1 --c2--
	  #		#n  --b-----  --d-----
	  #		#again b would inherit from c1, a, c1 and c2 have equal overlap, 
	  #		#a and c1 have equal coverage, but c1 is more central to b.
	  #		#is this the same in reverse?
	  #		#o---a---- --c--
	  #		#nb1 b2 -b3--
	  #		#this is the same as we don't account for direction here!!??
	  #		#this will take some funky logic in the assign block
	  #
	  #		#do mid point comparison first as this will always find the best contained feature first
	  #		#as these would fail the number of parents test as they will only have one
	  #		#conversely the midpoint test for linking features would always fail as 
	  #		#linking features will always score equally on the midpoint after scoring equally on the coverage
	  #	
	  #
	  #		#Displacement test
	  #		$tmp_displacement *= -1 if ($tmp_displacement < 0);
	  #	
	  #		if ($tmp_displacement < $displacement) {
	  #		  $new_child = 1;
	  #		} 
	  #		elsif ($tmp_displacement == $displacement) {
	  #		  warn "Displacement is failing to resolve stable_id mapping, need num parent rules!! Defaulting to 5' mapping";
	  #		  $new_child = 1;
	  #		
	  #
	  #		  #STOP HERE! THIS IS GETTNIG WAY MORE COMPLICATED THAN WE CAN HANDLE IN ONE GO!
	  #		  #this is slightly different and may need to be worked out in the assign block!!!!!!!!!!!!!!!!
	  #		  #We need to compare num of parent of parents, but we won't know this until we traverse the 
	  #		  #next transition, hence this needs to be dealt with in the assign block when we walk backwards
	  #		  #or can we do it here as we move along?
	  #		  #as these 'twins' may be exact matches over many transitions
	  #		  #so we need to set a 'twin' flag, and record which transition it occured at
	  #		  #then walk forward comparing to the mirrored transition and seeing which wins, then track back to the
	  #		  #transition where the twin occured assigning child and previous child accordingly
	  #		  #can we leave this walk back until the assign step?
	  #		  #do we also need to record 2nd child instead of previous child, as the next best child might be 3' rather than
	  #		  #the default which is the previous seen i.e. 5'
	  #	  
	  #		  #if($tmp_num_parents > $num_parents){
	  #
	  #		  #}elsif($tmp_num_parents == $num_parents){
	  #		  #warn "Number of parents is failing to resolve stable_id mapping, need more rules";
	  #
	  #		  #o -a c1 --c2---- x1 x2
	  #		  #n --b------  --d------
	  #
	  #		  #b should inherit from c2 due to increased number of children/parents
	  #		  #so when considering the c2> b & d transition
	  #		  #we currently don't know the number of parents of d or d
	  #		  #so we can't pick a current child
	  #		  #if they are truly equal then it's probably best to retire c2
	  #		  #and assign c1 to b and x1 to d
	  #
	  #		
	  #		  #num_parent is essentially $#overlaps, so we will need to trace this backwards 
	  #		  #through 3 transitions to be able to resolve this 
	  #		  #d  > c2(x1 & x2) 
	  #		  #c2 > b & d Don't know which is the child so can't assign an id yet
#		  #b  > a > c1 & c2 can now compare num parents and assign if there is a difference
	  #		  #else split c2 and use c1 and x1
	  #		  
	  #		  #what if parents have differing coverage?
	  #		  #o--a c1 --c2---- x1 x2
	  #		  #n --b------  --d------
	  #		  #b now inherit from a
	  #		  #and d from c2
	  #
	  #		  #what is a has better match?
	  #		  #o-----a c1 --c2---- x1 x2
	  #		  #n-z- --b------  --d------
	  #		  #we would have to work right back to z before assigning ids
	  #		  #z should inherit from a
	  #		  #b should inherit from previosu child i.e. c1
	  #		  #d should inherit from c2
	  #
	  #		  #That has to be it!!
	  #		  #another extension here would be caught by preceding rules no?
	  #		  #For walking forward then we find a twin, previous child has to be next best child i.e. 3' to current
	  #		  #we would have to have a nested loop here to compare back to the mirrored 'twins' to define
	  #		  #child/previous child for each before we pass to assign_ids?  Is this necessary, can we not resolves this in assign_ids?
	  #		  #}
	  #		}
	  #	  }
	  #	}
	  
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
	  if($nreg_feats[$i]->seq_region_end() > $feature->seq_region_end()){
		$overhang = 1;
		

		#Found non-maximal link
		$split_chain = 1 if(! $new_children{'child'});
	  }
	  elsif($i == $#nreg_feats){
		#Last feature but no overhang
		$end_chain = 1;
	  }
	  
	  #Found new child so update vars
	  if($new_children{'child'}){
		map $comparatees{$_}{'child'} = $comparatees{$_}{'tmp'}, keys %comparatees; 
		$previous_child = $current_child if (defined $current_child);
		$current_child = $i;
	  }elsif($new_children{'next_child'}){
		map $comparatees{$_}{'next_child'} = $comparatees{$_}{'tmp'}, keys %comparatees; 
		$next_child = $i;
	  }
	}
  }


  #undef $feature if we're ending a chain
  #no need to undef source_dbID as we're not testing it in the caller

  #Found end of chain or no overlap features
  if($no_new_features || ! @nreg_feats){
	#no_new_features ensures we don't add the last transition which has only parent as nreg_feat
	#or no features????
	#can we ever get !@nreg_features? due to source
	#this is handled above in extend block

	#This may be an empty hash as we can potentially break a chain on a non-maximal overlap
	#which gives rise to no more features, will find source

	if(@transitions){
	  #print "2 passing split_id $split_id\n";

	  $split_id = &assign_stable_ids(\@transitions, $split_chain, $split_id);
	  
	  #this should never happen
	  die("split chain and id with id $split_id with no new feature") if $split_id || $split_chain;

	  #split chain will always be false here
	}
	#else{
	#  warn "No transitions to assign to, must have found broken non-maximal 3' overhang(dbID:".
	#	$feature->dbID().") with no more overlapping features\n";
	#}

	undef $feature;
	
  }
  elsif(@nreg_feats){
	#Found some new overlapping features i.e. a transition

	#warn "pushing transition";

	#right then
	#here we decide which is the current_child based on whether
	#we have found a split
	#if this is after a split then we only use the next child if the current_child is 0
	#i.e. child class and the previous link/transition has already been assigned
	#then for futher iterations we either revert to using current child if this
	#next_child is not a linking overlap i.e. $#nreg_feats
	#we should pass a next child flag to catch cases where we try and use previous in the assign step
	#as this should logically never happen

	#what happens when we have no next child
	#we need to assign a new stable_id
	#is this caight in the assign block?


	if($use_next_child){
	  
	  #warn "Using next child $next_child ".$nreg_feats[$next_child]->stable_id().
	#	" instead of current_child $current_child ".$nreg_feats[$current_child]->stable_id();
	  
	  #Only use next if current is a clash with the one used in the previous split link
	  $current_child = $next_child if $current_child == 0;
	  
	 # warn "current child is now $current_child";
	  #now we have lost current_child, and previous child is invalid?

	  #reset ot normal behaviour if current child is not overlap and has no chance of clash
	  #we also need to check whether this is an overlap??
	
	}


	push @transitions,
	  {(
		feature          => $feature,
		source           => $from_db, #current feature feature_set, OLD or NEW
		#overlap         => undef, #bp overlap with current child
		#coverage        => $coverage, #% of current child overlapping feature
		overlap_features => \@nreg_feats,
		current_child    => $current_child, #index of current child in tmp new_feats array
		previous_child   => $previous_child, #used to identify inheritance if it get's shifted


		#this is being immediately overwritten above
		no_use_previous  => $use_next_child,#don't really need this either



		#overhang         => $overhang,
		#next_child       => $next_child,
		#true_orphan     => ??????,
		#twin?? this could be figured out if previous child is > current_child
		#pass %comparatees?
	   )};

	if($use_next_child){

	  if($current_child != $#nreg_feats || ($current_child = $#nreg_feats  && (! $overhang))){
		$use_next_child = 0;
	  }
	}


	#warn "use next child for next transition will be $use_next_child";

	#we want to split chain here if last overlap was non-maximal
	if($end_chain || $split_chain){
	 # warn "end chain $end_chain split chain $split_chain";
	 # warn "passing split_id $split_id" if $split_id;

	  $split_id = &assign_stable_ids(\@transitions, $split_chain, $split_id);#, $rest_of_chain);

	  
	  if($end_chain){
		#make sure we don't pass back current values and get caught in loop
		undef $feature;
		#undef $source_dbID;
	  }
	  else{#split
		#set source/parent dbID for next transtion
		#and set new feature to overhang feature
		$source_dbID = $feature->dbID();
		$feature = $nreg_feats[$#nreg_feats];
		$use_next_child = 1; #this needs to be returned
	  }
	}
  }
  else{#no mapping/features???
	warn "not catching this!!!";
  }
  
  #is split_id performing the same task as use_next_child?
  #split_id is the last overhanging feature id in a split chain


  #finally return the mapping info hash and switch the db

  #warn "returning split id to build $split_id" if $split_id;

  return ($to_db, $feature, $source_dbID, $use_next_child, $split_id);
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
  my ($transitions, $split_chain, $split_id) = @_;#, $rest_of_chain) = @_;

  my ($i, $child_id, $feature_id, $last_stable_id, $orphan_id);
  my $num_trans = $#{$transitions};

  #Set child to 1, so avoid use previous at transition end e.g. 1 to 1 mapping
  my $child = 1;

  $helper->debug(1, "Assigning stable IDs to ".($num_trans+1)." transitions");
  $helper->debug(2, "split chain!!") if $split_chain;
  $helper->debug(2, "previous split_id of $split_id") if $split_id;
  

  #Walk backwards through transitions
  for($i = $num_trans; $i >= 0; $i--){
	my ($feature, $orphan, $orphan_id, @new_stable_ids);
	my $last_orphan = $#{$transitions->[$i]{'overlap_features'}};
	#warn "Processing transtion $i\n";

	#would need to add stuff here to account for the twin situation
	#would need next_child, rather than previous child
	#more

	
	#reset child for current transition
	##Only do this if we are not at the end of the chain??
	#why?


	#if($i != $num_trans){


	#This currently fails on a one to one mapping
	#we only set previous if we are not at the end
	#as if we are at the end then either it will be an overlap with no further transitions
	#or a non-maximal overlap with more transitions
	#so child here is the last child
	#this was the reason we set child to 1 as default;

	if($child == 0 && 
	   ($transitions->[$i]->{'current_child'} == $#{$transitions->[$i]->{'overlap_features'}})
	   && ! $transitions->[$i]->{'no_use_previous'}){
	  #child clash
	  #warn "current child and length of overlap features is ".$transitions->[$i]->{'current_child'};
	  $child =  $transitions->[$i]->{'previous_child'};
	  #warn "child clash setting to previous child $child\n";
	}
	#elsif(){
	#NOW HANDLED in the build_transitions block
	##This is for the use next child case if we are dealing with the >1st part of a split chain
	#and only dealing with the first link/transition in that chain
	#This may have a knock on effect to the rest of the chain
	#or will this have been resolved in the define transitions step?
	#should we handle all of this in the build_transitions step
	#as we don't know what impact there will be when working backwards through the transitions here
	#but we can deal with that when we actually work forwards when we build the transitions
	#we also generate this problem when working forwards by splitting the chain
	#so we only actually need next child when we split
	#and so long as that is not the last feature then we can carry on as noraml
	#other wise we have to check the there is not a child clash on the overlapping feature
	#else we need to keep checking the next child, as the normal child has already been taken by the previous
	#link/transition
	
	#we shouldn't ever have the problem of then using previous here as this should never occur given the prior context of knowing at least part of the projection/inheritance
	
				
	#}
	else{#no clash
	  $child = $transitions->[$i]->{'current_child'};
	  #warn "no clash setting to current_child $child\n";
	}

	# we need a case for next child also!!!!!!!!!!!!!!!!!!!!!!!!!!
	#for the case where we split a chain, then find the next bit of the chain has a first child which is shared with
	#the previous link.
	#we still need to add the stable_id to the previous mapping cache
	#can't avoid the next child!!
	
	#}

	#warn "child is $child\n";


	#INHERITANCE
	if ($transitions->[$i]->{'source'} eq 'NEW'){

	  $child_id = $transitions->[$i]->{'overlap_features'}->[$child]->{'stable_id'};
	  $feature_id = $transitions->[$i]->{'feature'}->dbID();
	  
	  $helper->debug(1, "Assigning IDs by Inheritance. Source is NEW dbID $feature_id stable ID $child_id from child $child");


	  #This should be true as will be using previous if last transition projected to this feature
	  if(exists $dbid_mappings{$feature_id}){
		die("Found duplicate inheritance:\tdbID ".$feature_id." > ".
			  $dbid_mappings{$feature_id}->{'stable_id'}." & $child_id");
	  }
	  
	  
	  #only one per transition here as we're only dealing with one new feature
	  $mappings++;
	  
	  #Assign mapping????
	  #Should this be done here 


	  #warn "Assign stable id $child_id to dbID $feature_id";
	  $transitions->[$i]->{'feature'}->stable_id($child_id);
	  $dbid_mappings{$feature_id} = $transitions->[$i]->{'feature'};
	  $helper->debug(1, "Assigned dbid_mapping $feature_id => ".$transitions->[$i]->{'feature'}->stable_id());
	  
	  #record all stable_id splits
	  foreach my $orphan_cnt(0..$last_orphan){#these are stable_ids!
		$orphan = $transitions->[$i]->{'overlap_features'}->[$orphan_cnt];
		$orphan_id = $orphan->{'stable_id'};
		
		#warn "Found $orphan_cnt orphan with stable_id $orphan_id\n";

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

			if(! ($split_id && ($i == $num_trans))){
			  die("Have found duplicate mappings on very 5' feature in a transition chain for $orphan_id: previous ".
				  join(',', @{$mapping_cache{$orphan_id}})." new $child_id");


		

			}

			#This is because chain has been split and this transition knows nothing about
			#this shared child between the two transitions
			#this is a case for not splitting
			#but can we overcome this by simply passing a flag to use the next child
			#if we have split the chain? bu tonly for the first transition!
			#can this not happen for the last transition also as with Projection? and vice versa?

			#Populate mapping cash from split chain with new inherited ID
			#warn "Assigning child $child_id to previous split mapping_cache $split_id";
			push @{$mapping_cache{$split_id}}, $child_id;
			undef $split_id;
		  }
		  else{
			#is either entirely contained
			#or is 3' with no linking overlap
			#or is 5' with an unknown overlap
			#i.e. has only one known new stable_id(this feature)
			@{$mapping_cache{$orphan_id}} = ($child_id);
		  }
		}
	  }

																										 
	  #now we need to add the new stabke_id to the old mapping cache
		#																								 print "mapping cache for $orphan_id is ".join(',', @{$mapping_cache{$orphan_id}});
	  

	  #set last old stable_id here enables updating of split cache in next transition
	  $last_stable_id = $orphan_id;
	}
	else{#PROJECTION i.e. source == OLD
	  $child_id = $transitions->[$i]->{'overlap_features'}->[$child]->dbID();
	  $feature_id = $transitions->[$i]->{'feature'}->{'stable_id'};

	  $helper->debug(1, "Assigning IDs by Projection. Source is OLD stable id $feature_id to child dbID $child_id");

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
		die("Found duplicate projection:\tdbID $child_id < stable_ids ".$feature_id.
			  " & ".$dbid_mappings{$child_id}->stable_id());
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
		my $new_sid;
		$orphan = $transitions->[$i]->{'overlap_features'}->[$orphan_cnt];
		$orphan_id = $orphan->dbID();#new dbID

		#warn "Found orphan $orphan_cnt with dbID $orphan_id\n";
		#warn "new_sids are @new_stable_ids";

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
		  #warn "\n\n\n\nis last orhpan $last_orphan";
		  
		  if(exists $dbid_mappings{$orphan_id}){
			#already assigned from last transition(cannot be child)
			#This will have already been added to the previous split cache as it will be the previous child id

			#warn "already have mapping";

			$new_sid = $dbid_mappings{$orphan_id}->{'stable_id'};
			#added to split cache below
		  }
		  #elsif($orphan_cnt == $child){
		  else{

			#warn "orphan count is $orphan_cnt and child is $child";

			if($orphan_cnt == $child){
			  #last orphan is child
			  $orphan->stable_id($feature_id);
			  #warn "1 setting dbid_mappings for $orphan_id to $orphan";
			  $dbid_mappings{$orphan_id} = $orphan;
			  $new_sid = $orphan->{'stable_id'};
			}
			else{
			  #can only be last in chain as we don't have a link

			  #so we only assign if this is not an overlapping child
			  #if it is then we have to account for the fact that the next part of the chain 
			  #this child may want to inherit from the current feature
			  #and that we need to add this to the mapping cash for this feature now
			  #or in the next assign processes
			  #this requires know what the last mapping cache was between assign calls
			  
			  #we also need to know whether there is any of the chain left otherwise we will skip the assignment
			  #do we need to know overlap? Just that there is no more links?
			  
			  #don't really need num_trans test here as it is implicit through the logic

			  if(($i == $num_trans) &&
				 (! $split_chain)){# || (! $transitions->[$i]->{'overlap'})){
				$new_sid = &assign_and_log_new_stable_id($orphan);
				#warn "Assigned $new_sid to last feature";
			  }else{
				#this is mid chain overlap which should either be the child
				#or this is the last transition and it is a split chain 

				die('Found mid chain overlapping orphan!') if($i != $num_trans);

				#split_chain
				#so we need to wait for inheritance in next assinment and then populate 
				#this mapping cache with the new inherited ID

				#set split_id to pass back for next assignment call
				$split_id = $feature_id;


				#how do we account for split_id when passed?!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
				#let's deal with it in the Inheritance block first

			  }

			  #warn "last in chain, assigned $new_sid" if $new_sid;
			  #warn "last in split chain, returning split_id $split_id" if $split_id;
			}

			#we need add this to the previous split_merge_cache
			if(defined $last_stable_id){
			  
			  if(exists $mapping_cache{$last_stable_id}){
				push @{$mapping_cache{$last_stable_id}}, $orphan->{'stable_id'};
			  }
			  else{
				die("Failed to find split_merge_cache for old $last_stable_id ".
					  "when updating with linking new feature ".$orphan->stable_id());
			  }
			}
		  }
		}
		elsif($orphan_cnt == $child){
		  #this will also account for the child being the 5' orphan
		  $orphan->stable_id($feature_id);
		  $dbid_mappings{$orphan_id} = $orphan;
		  $new_sid = $orphan->{'stable_id'};
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
		  #warn "not child and contained orphan or 5' orphan at start of chain";
		  $new_sid = &assign_and_log_new_stable_id($orphan);
		}
		else{#5' orphan in the middle of a chain
		  #this is nor really an orphan as this wil be the child of the next transition
		  #warn "waa next";
		  #don't assign stable id yet, update mapping cache in next transition
		  next;	#to ensure we populate this transition's mapping cache
		}
		
   		if ($orphan_cnt == $child){
		  unshift @new_stable_ids, $new_sid;
		}elsif(! $split_id){
		  #don't do this if we haven't assigned due to split chain
		  #warn "pushing non child stable_id $new_sid";
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

	  #warn "keys are now ".join(',', keys %dbid_mappings);
	  #warn "Mapping cache for $feature_id is @new_stable_ids\n";

	}
  }
  
  #warn "dbid_mappings:\n".join(',',keys %dbid_mappings);
  #warn "returning split is $split_id" if $split_id;
  #return the ID which needs it's mapping cash populating with
  #the ambiguous as yet unassinged feature
  #This will be picked up in the first transition(last in loop)
  return $split_id;
}

=head2 assign_and_log_new_stable_id

  Arg [1]    : Bio::EnsEMBL:Funcgen::RegulatoryFeature
  Example    : my $new_sid = &assign_and_log_new_stable_id($new_reg_feature);
  Description: Creates and assigns a new stable ID to a given RegulatoryFeature.
               New dbID to stable ID mappings are logged in the new_ids file.
  Returntype : string - Stable ID
  Exceptions : Throws if no valid RegulatoryFeature is passed
  Status     : At Risk - Use Bio::EnsEMBL::Funcgen::RegulatoryFeature

=cut


sub assign_and_log_new_stable_id{
  my $new_reg_feat = shift;

  if(! (ref($new_reg_feat) && $new_reg_feat->isa('Bio::EnsEMBL::Funcgen::RegulatoryFeature'))){
	die('You must provide a valid Bio::EnsEMBL::Funcgen::RegulatoryFeature');
  }

  #Need to add species intrafix here
  $new_reg_feat->stable_id($next_stable_id);
  #warn "Assigned new stable ID ".$new_reg_feat->stable_id()." to RegulatoryFeature ".$new_reg_feat->dbID()."\n"; 

  $dbid_mappings{$new_reg_feat->dbID()} = $new_reg_feat;
  $new_mappings ++;
  $next_stable_id ++;

  #Really need to add this to an IO buffer to improve output efficiency
  print $new_id_handle $new_reg_feat->dbID()."\t".$new_reg_feat->stable_id()."\n";

  return $new_reg_feat->{'stable_id'};
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
