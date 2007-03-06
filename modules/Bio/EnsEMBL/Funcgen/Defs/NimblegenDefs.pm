#
# EnsEMBL module for Bio::EnsEMBL::Funcgen::Defs::NimblegenDefs
#

=head1 NAME

Bio::EnsEMBL::Funcgen::Defs::NimblegenDefs

=head1 SYNOPSIS

  my $defs_type = "Bio::EnsEMBL::Funcgen::Defs::NimblegenDefs";
  push @INC, $defs_type;
  my $imp = $class->SUPER::new(@_);


=head1 DESCRIPTION

This is a definitions class which should not be instatiated directly, it 
normally inherited by the Importer.  NimblegenDefs contains meta data and 
methods specific to NimbleGen arrays to aid parsing and importing of 
experimental data.

=head1 AUTHOR

This module was written by Nathan Johnson.

=head1 CONTACT

Post questions to the EnsEMBL development list ensembl-dev@ebi.ac.uk

=head1 METHODS

=cut

package Bio::EnsEMBL::Funcgen::Defs::NimblegenDefs;

use Bio::EnsEMBL::Funcgen::Array;
use Bio::EnsEMBL::Funcgen::ProbeSet;
use Bio::EnsEMBL::Funcgen::Probe;
use Bio::EnsEMBL::Funcgen::ProbeFeature;
use Bio::EnsEMBL::Funcgen::FeatureType;
use Bio::EnsEMBL::Funcgen::ExperimentalChip;
use Bio::EnsEMBL::Funcgen::ArrayChip;
use Bio::EnsEMBL::Funcgen::Channel;
use Bio::EnsEMBL::Utils::Exception qw( throw warning deprecate );
use Bio::EnsEMBL::Funcgen::Utils::EFGUtils qw(species_chr_num open_file);
use Bio::EnsEMBL::Utils::Argument qw( rearrange );
use Bio::EnsEMBL::Funcgen::Helper;
use strict;

use vars qw(@ISA);
@ISA = qw(Bio::EnsEMBL::Funcgen::Helper);

=head2 new

  Example    : my $self = $class->SUPER::new(@_);
  Description: Constructor method for NimblegenDefs class
  Returntype : Bio::EnsEMBL::Funcgen::Defs::NimblegenDefs
  Exceptions : throws if Experiment name not defined or if caller is not Importer
  Caller     : Bio::EnsEMBL::Funcgen::Importer
  Status     : at risk

=cut


sub new{
  my $caller = shift;

  my $class = ref($caller) || $caller;
  my $self  = $class->SUPER::new();

  throw("This is a skeleton class for Bio::EnsEMBL::Importer, should not be used directly") if(! $self->isa("Bio::EnsEMBL::Funcgen::Importer"));
	
  my ($name) = rearrange(['NAME'], @_);

  throw('Must provide an Experiment name for a Nimblegen import') if ! defined $name;

  #should we provide args override for all of these?

  $self->{'defs'} =  
    {(
      #order of these data arrays is important!
      array_data   => ['experiment'],
      probe_data   => ["probe"],
      results_data => ["results"],
      sample_key_fields => ['DESIGN_ID', 'CHIP_ID', 'DYE', 'PROMOT_SAMPLE_TYPE', 'SAMPLE_LABEL'],
      # 'SAMPLE_DESCRIPTION removed due to naming disparities
      ndf_fields      => ['CONTAINER', 'PROBE_SEQUENCE', 'MISMATCH', 'FEATURE_ID', 'PROBE_ID'],
      pos_fields      => ['CHROMOSOME', 'PROBE_ID', 'POSITION', 'COUNT'],
      result_fields   => ['PROBE_ID', 'PM'],
      notes_fields   => ['DESIGN_ID', 'DESIGN_NAME', 'DESCRIPTION'],
      norm_method => 'VSN_GLOG',
      dye_freqs => {(
		     Cy5 => 635,
		     Cy3 => 532,
		    )},
     )};
	
  return $self;
}


=head2 set_defs

  Example    : my $self->set_defs;
  Description: Sets attribute dependent defs
  Returntype : None
  Exceptions : None
  Caller     : Bio::EnsEMBL::Funcgen::Importer
  Status     : at risk

=cut


sub set_defs{
  my $self = shift;

  #dir are not set in defs to enable generic get_dir method access

  $self->{'design_dir'} = $self->get_dir("data")."/input/".
    $self->vendor()."/".$self->name()."/DesignFiles";
    
  $self->{'defs'}{'chip_file'} = $self->get_dir("data")."/input/".
    $self->vendor()."/".$self->name()."/SampleKey.txt";
    
  $self->{'defs'}{'notes_file'} = $self->get_dir("data")."/input/".
    $self->vendor()."/".$self->name()."/DesignNotes.txt";
  
  $self->{'results_dir'} = $self->get_dir("data")."/input/".
    $self->vendor()."/".$self->name()."/PairData";

  return;
}




=head2 read_array_data

  Example    : $imp->read_array_data();
  Description: Parses NimbleGen DesignNotes.txt files to create and store new Arrays
  Returntype : none
  Exceptions : None
  Caller     : general
  Status     : At risk - Can this be generic? Can we force the creation of a DesignNotes file on other formats?

=cut


sub read_array_data{
  my ($self, $notes_file) = @_;

  my ($line, $array, $array_chip, @data, %hpos);
  my $oa_adaptor = $self->db->get_ArrayAdaptor();
  my $ac_adaptor = $self->db->get_ArrayChipAdaptor();

  #Slurp file to string, sets local delimtter to null and subs new lines
  my $fh = open_file("<", $notes_file);
  #($design_desc = do { local ($/); <$fh>;}) =~ s/\r*\n$//;
  #close($fh);
  
  while ($line = <$fh>){
    
    $line =~ s/\r*\n//;#chump
    @data =  split/\t/o, $line;

    #We need to have a DESIGN vendor type?
    #also need to be able to set file path independently of defs
    
    if($. == 1){
      %hpos = %{$self->set_header_hash(\@data, $self->get_def('notes_fields'))};
      next;
    }

    ### CREATE AND STORE Array and ArrayChips  
    if(! defined $array ){
      #This is treating each array chip as a separate array, unless arrayset is defined
      #AT present we have no way of differentiating between different array_chips on same array???!!!
      #Need to add functionality afterwards to collate array_chips into single array
           
      #This will use a stored array if present

      $array = Bio::EnsEMBL::Funcgen::Array->new
	(
	 -NAME        => $self->array_name() || $data[$hpos{'DESIGN_NAME'}],
	 -FORMAT      => uc($self->format()),
	 -VENDOR      => uc($self->vendor()),
	 -TYPE        => 'OLIGO',
	 -DESCRIPTION => $data[$hpos{'DESCRIPTION'}],#need to trim the array chip specific description here
	);

      ($array) = @{$oa_adaptor->store($array)};

      $array_chip = Bio::EnsEMBL::Funcgen::ArrayChip->new(
							  -ARRAY_ID  => $array->dbID(),
							  -NAME      => $data[$hpos{'DESIGN_NAME'}],
							  -DESIGN_ID => $data[$hpos{'DESIGN_ID'}],
							  #add description?
							 );

      #This will use a stored array_chip if present
      ($array_chip) = @{$ac_adaptor->store($array_chip)};
      $array->add_ArrayChip($array_chip);
        
    }
    elsif((! $array->get_ArrayChip_by_design_id($data[$hpos{'DESIGN_ID'}])) && ($self->array_set())){
      
      $self->log("Generating new ArrayChip(".$data[$hpos{'DESIGN_NAME'}].". for same Array ".$array->name()."\n");
      
      $array_chip = Bio::EnsEMBL::Funcgen::ArrayChip->new(
							  -ARRAY_ID  => $array->dbID(),
							  -NAME        => $data[$hpos{'DESIGN_NAME'}],
							  -DESIGN_ID => $data[$hpos{'DESIGN_ID'}],
							 );
      
      ($array_chip) = @{$ac_adaptor->store($array_chip)};
      $array->add_ArrayChip($array_chip);
      
    }
    elsif(! $array->get_ArrayChip_by_design_id($data[$hpos{'DESIGN_ID'}])){
      throw("Found experiment with more than one design without -array_set");
    }
  }
  

  $self->add_Array($array);

  close($fh);
  
  return;

}


=head2 read_experiment_data

  Example    : $imp->read_array_chip_data();
  Description: Parses and imports array & experimental chip meta data/objects
  Returntype : none
  Exceptions : throws if more than one array/design found and not an "array set"
  Caller     : Importer
  Status     : At risk - likely to be replaced/moved with/to VendorDefs

=cut

sub read_experiment_data{
  my $self = shift;

  
  $self->read_array_data($self->get_def('notes_file'));

  my ($design_desc, $line, $tmp_uid, $channel, $echip);
  my ($sample_desc, %hpos, @data);
  my $ec_adaptor = $self->db->get_ExperimentalChipAdaptor();
  my $chan_adaptor = $self->db->get_ChannelAdaptor();
  
  #Currently 1 design = 1 chip = 1 array /DVD
  #Different designs are not currently collated into a chip_set/array in any ordered manner
  #Register each design as an array and an array_chip
  #May want to group array_chips into array/chip sets by association though the API
  
 
  warn("Harcoded for one array(can have multiple chips from the same array) per experiment\n");
  my $fh = open_file("<", $self->get_def("chip_file"));
  $self->log("Reading chip data");

  warn("We need to change Array/ArrayChip retrieval to ensure that we have IMPORTED status, so avoid having an incomplete arraychip");
  warn "Do we need to validate each line here against the header array?";
  while ($line = <$fh>){
    next if $line =~ /^\s+\r*\n/;
    $line =~ s/\r*\n//;#chump
    @data =  split/\t/o, $line;
    #we could validate line against scalar of header array
     #ORD_ID  CHIP_ID DYE     DESIGN_NAME     DESIGN_ID       SAMPLE_LABEL    SAMPLE_SPECIES  SAMPLE_DESCRIPTION      TISSUE_TREATMENT        PROMOT_SAMPLE_TYPE
    if ($. == 1){
      %hpos = %{$self->set_header_hash(\@data, $self->get_def('sample_key_fields'))};

      #we need to set the sample description field name, as it can vary :(((
      @data = grep(/SAMPLE_DESCRIPTION/, keys %hpos);
      $sample_desc = $data[0];
      throw("More than one sample description(@data) in ".$self->get_def("chip_file")."\n") if(scalar @data >1);
      next;
    }
    
    #Need to handle array class here i.e. two channel arrays will have two lines
    #validate species here
    #look up alias from registry and match to self->species
    #registry may not be loaded for local installation


    ### CREATE AND STORE ExperimentalChips
    if ((! $tmp_uid) || ($data[$hpos{'CHIP_ID'}] ne $tmp_uid)){

      #Test both channels are available, i.e. the SampleKey has two TOTAL channels
      if($echip){

	for my $type('TOTAL', 'EXPERIMENTAL'){

	  my $test_chan =  $chan_adaptor->fetch_by_type_experimental_chip_id($type, $echip->dbID());
	  throw("ExperimentalChip(".$echip->unique_id().
		") does not have a $type channel, please check the SampleKey.txt file") if ! $test_chan;

	}
      }
      
      $tmp_uid = $data[$hpos{'CHIP_ID'}];
      $echip = $ec_adaptor->fetch_by_unique_id_vendor($data[$hpos{'CHIP_ID'}], 'NIMBLEGEN');

      if($echip){

	if(! $self->recovery()){
	  throw("ExperimentalChip(".$echip->unqiue_id().
		" already exists in the database\nMaybe you want to recover?");
	}
      }else{
      
   	$echip =  Bio::EnsEMBL::Funcgen::ExperimentalChip->new
	  (
	   -EXPERIMENT_ID  => $self->experiment->dbID(),
	   -DESCRIPTION    => $data[$hpos{$sample_desc}],
	   -ARRAY_CHIP_ID  => $self->arrays->[0]->get_ArrayChip_by_design_id($data[$hpos{'DESIGN_ID'}])->dbID(),
	   -UNIQUE_ID      => $data[$hpos{'CHIP_ID'}],
	  );
      
	($echip) = @{$ec_adaptor->store($echip)};
	$self->experiment->add_ExperimentalChip($echip); #need a contains method here, #need to do this as the probe cache init 'experimental_chips'
      }
    }


    ### CREATE AND STORE Channels
    $channel = $chan_adaptor->fetch_by_type_experimental_chip_id(uc($data[$hpos{'PROMOT_SAMPLE_TYPE'}]), $echip->dbID());

    if($channel){
      if(! $self->recovery()){
	throw("Channel(".$echip->unqiue_id().":".uc($data[$hpos{'PROMOT_SAMPLE_TYPE'}]).
	      " already exists in the database\nMaybe you want to recover?");
	}
    }else{
      #Handles single/mutli
      $channel =  Bio::EnsEMBL::Funcgen::Channel->new
	(
	 -EXPERIMENTAL_CHIP_ID => $echip->dbID(),
	 -DYE                  => $data[$hpos{'DYE'}],
	 -SAMPLE_LABEL         => $data[$hpos{'SAMPLE_LABEL'}],
	 -TYPE                 => uc($data[$hpos{'PROMOT_SAMPLE_TYPE'}]),
      );
    
      #-SPECIES              => $self->species(),#on channel/sample to enable multi-species chip/experiment
      #would never happen on one chip?  May happen between chips in one experiment

      ($channel) = @{$chan_adaptor->store($channel)};
      #$echip->add_Channel($channel);
    }    
  }
  
  #Set the array now we have stored all details

  close($fh);

  return;
}



=head2 read_probe_data

  Example    : $imp->read_probe_data();
  Description: Parses and imports probes, probe sets and features of a given array
  Returntype : none
  Exceptions : none
  Caller     : Importer
  Status     : Medium - Move parts to "Vendor"Defs.pm, should function the same

=cut


#Assumes one chip_design per experimental set.
sub read_probe_data{
  my ($self) = shift;
  
  my ($fh, $line, @data, %hpos, %probe_pos);#, %duplicate_probes);

  $self->log("Parsing ".$self->vendor()." probe data (".localtime().")");
  #can we do a reduced parse if we know the array chips are already imported");
  

  ### Read in
  #ndf file: probe_set, probe and probe_feature(.err contains multiple mappings)
  #pos file: probe chromosome locations
  #ndf file: build, chr, start stop,
  
  #Need to test whether array exists and is validated before importing probe data.
  #If already present, then can skip and just import results using probe_name to get probe_id and link to probe_feature.
  #just needs to create a mapping between probe_feature_ids and results
  #What if an array has been imported with only a few chips present?
  #Need check on size in import_array
  #e.g. someone has done an experiment with only a few chips from an array?
  #This would mean absent probes
  
  #also need to change how probe_names are generated for nimblegen?
  #native probe_ids may not be unique, but should be when combined with the seq_id which is currently being used as the xref_id
  #Handle with API!!
  
  #READ REGION POSITIONS
  #We need to handle different coord systems and possibly different assmemblies
  my $slice_a = $self->db->get_SliceAdaptor();
  my $cs = $self->db->get_FGCoordSystemAdaptor()->fetch_by_name_schema_build_version(
										     'chromosome', 
										     $self->db->_get_schema_build($self->db->dnadb())
										    );

  #should really test and store here if not valid, but this will be handled by adaptor store methods
  #This should really focus on the experiemental chips, 
  #altho' this amounts to the same thing if the chip has already been imported as




  

  
  foreach my $array(@{$self->arrays()}){
    
    foreach my $achip(@{$array->get_ArrayChips()}){

      my (@log);


      #warn "We need to account for different cs feature amppings here


      if($achip->has_status('IMPORTED')){
	$self->log("Skipping fully imported ArrayChip:\t".$achip->design_id());
	next;
      }elsif($self->recovery()){
	$self->log("Rolling back ArrayChip:\t".$achip->design_id());
	$self->db->rollback_ArrayChip($achip);
      }
      
      $self->log("Importing ArrayChip:".$achip->design_id());
      #mgd files are for a different product?
      #THIS BLOCK DOES NOT ACCOUNT FOR MULTIPLE ARRAYS PROPERLY, WOULD HAVE TO IMPLEMENT ARRAY SPECIFIC CACHES
      #all out files are generic, but are we converting to adaptor store?
      #      $fh = open_file("<", $self->get_def("design_dir")."/".$ac{'design_name'}.".ngd");
      #      my ($start, $stop, %regions, %probe_pos);
      #	  
      #      #May not have both ngd and pos file?
			#
      #      while ($line = <$fh>){
      #	$line =~ s/\r*\n//;#chump
      #	@data =  split/\||\t/o, $line;
      #	
      #	
      #	#SEQ_ID	SEQ_UNIQUE|BUILD|CHROMOSOME|LOCATION|DESCRIPTION|DATE_ENTERED|SOURCE_DB
      #	if ($. == 1){
      #	  %hpos = %{$self->set_header_hash(\@data)};
      #	  next;
      #	}
      #	
      #	#What about strand!!!!!!!!!!!
      #	$data[$hpos{'CHROMOSOME'}] =~ s/chr//;
      #	($start, $stop) = split/-/o, $data[$hpos{'LOCATION'}];
      #	
      #	#Do we need seq_id check here for validity?
      #	#overkill?
			#	if(exists $regions{$data[$hpos{'SEQ_ID'}]}){
      #	  croak("Duplicate regions\n");
      #	}else{
      #	  #$data[$hpos{'CHROMOSOME'}] = species_chr_num($self->species(), 	$data[$hpos{'CHROMOSOME'}]);
      #	  
      #	  #Set region hash for SEQ_ID
      #	  #Need to look up seq_region id here for given build
      #	  #Build should be manually specified as we can't guarantee it will be in the correct format
      #	  #or present at all
      #	  
      #	  $regions{$data[$hpos{'SEQ_ID'}]} = 
      #	    {
      #	     start => $start,
      #	     stop  => $stop,
      #	     seq_region_id => $self->get_chr_seq_region_id($data[$hpos{'CHROMOSOME'}], $start, $stop),
      #	     coord_system_id => $cs->dbID(),
      #	    };
      #	}
      #	
      #      }
      #      
      #      close($fh);
      
	
      
      #Always use pos file, ndf file cannot be guranteed to contain all location info
      #pos file also gives a key to which probes should be considered 'EXPERIMENTAL'
      
      #SLURP PROBE POSITIONS
      $fh = open_file("<", $self->get_dir("design")."/".$achip->name().".pos");
      
      #don't % = map ! Takes a lot longer than a while ;)
      while($line = <$fh>){
	$line =~ s/\r*\n//o;#Not using last element
	@data =  split/\t/o, $line;
	
	#SEQ_ID	CHROMOSOME	PROBE_ID	POSITION	COUNT
	if ($. == 1){
	  %hpos = %{$self->set_header_hash(\@data, $self->get_def('pos_fields'))};
	  next;
	}
	
	#Skip probe if there is a duplicate
	
	if(exists $probe_pos{$data[$hpos{'PROBE_ID'}]}){
	  throw("Found duplicate mapping for ".$data[$hpos{'PROBE_ID'}]." need implement duplicate logging/cleaning");
	  #need to build duplicate hash to clean elements from hash
	  # $duplicate_probes{$data[$hpos{'PROBE_ID'}]} = 1;
	  #next;
	}
	
	#
	if(! $self->cache_slice($data[$hpos{'CHROMOSOME'}])){
	  push @log, "Skipping probe ".$data[$hpos{'PROBE_ID'}]." with non-standard region ".$data[$hpos{'CHROMOSOME'}];
	}
	
	$probe_pos{$data[$hpos{'PROBE_ID'}]} = {(
						 chr => $data[$hpos{'CHROMOSOME'}],
						 start => $data[$hpos{'POSITION'}],
						)};
	
      }
      
      
      #Remove duplicate probes 

      
      close($fh);
      
      
      #OPEN PROBE IN/OUT FILES
      $fh = open_file("<", $self->get_dir("design")."/".$achip->name().".ndf");
      #Need to set these paths in each  achip hash, file names could be tablename.chip_id.txt
      #my $p_out = open_file(">", $self->get_dir("import")."/probe.".$ac{'design_name'}."txt");
      #my $ps_out = open_file(">", $self->get_dir("import")."/probe_set.".$ac{'design_name'}.".txt");
      #my $pf_out = open_file(">", $self->get_dir("import")."/probe_feature.".$ac{'design_name'}."txt");
      my $f_out = open_file(">", $self->get_dir("output")."/probe.".$achip->name()."fasta")	if($self->{'_dump_fasta'});
      my ($length, $ops, $op, $of, %pfs);

      #should define mapping_method arg to allows this to be set to LiftOver/EnsemblMap
      my $anal = $self->db->get_AnalysisAdaptor()->fetch_by_logic_name("VendorMap");
			
		
      my $strand = 0;	#default for nimblegen, should be defs hash?
      my $cig_line = "50M";	#default for nimblegen, should be defs hash?
      my $fasta = "";
      
      #$self->Timer()->mark("Starting probe loop");
      
      while($line = <$fh>){
	$line =~ s/\r*\n//;
	@data =  split/\t/o, $line;
	my $loc = "";
	my $class = "EXPERIMENTAL";
	
	#PROBE_DESIGN_ID	CONTAINER	DESIGN_NOTE	SELECTION_CRITERIA	SEQ_ID	PROBE_SEQUENCE	MISMATCH	MATCH_INDEX	FEATURE_ID	ROW_NUM	COL_NUM	PROBE_CLASS	PROBE_ID	POSITION	DESIGN_ID	X	Y
				#2067_0025_0001  BLOCK1          0       chrX    TTAGTTTAAAATAAACAAAAAGATACTCTCTGGTTATTAAATCAATTTCT      0       52822449        52822449        1       25      experimental    chrXP10404896   10404896        2067    25      1
				
	if ($. == 1){	
	  %hpos = %{$self->set_header_hash(\@data, $self->get_def('ndf_fields'))};
	  next;
	}
		  
	
	if (! exists $probe_pos{$data[$hpos{'PROBE_ID'}]}){
	  push @log, "Skipping non-experimental probe:\t".$data[$hpos{'PROBE_ID'}];
	  next;
	}
	
	#Which non-experimental probes might we want to store?
				#if($data[$hpos{'CONTAINER'}] =~ /control/io){
	#	$class = "CONTROL";
	#}
	#elsif($data[$hpos{'CONTAINER'}] =~ /random/io){
	#	$class = "RANDOM";
	#}
	#elsif($data[$hpos{'PROBE_CLASS'}] !~ /experimental/io){
	#	$class = "OTHER";
	#}
	#elsif(! exists $probe_pos{$data[$hpos{'PROBE_ID'}]}){	#HACKY HACKY HACK HACK!! Needed for valid region retrival
	#  $class = "OTHER";
	#}
	#SPIKE INS?
	
	

	#This assumes all probes in feature/probeset are next to each other!!!!!!!!!!

	
	if($data[$hpos{'FEATURE_ID'}] != $data[$hpos{'MATCH_INDEX'}]){#Probe set data
	  #print "Generating new probeset:\tFeature id:\t".$data[$hpos{'FEATURE_ID'}]."\tmatchindex:\t".$data[$hpos{'MATCH_INDEX'}]."\n";
	  
	  if($ops && ($data[$hpos{'FEATURE_ID'}] ne $ops->name())){
	    #THis is where we chose to update/validate
	    #Do we need to pass probes if they're already stored..may aswell to reduce mysql load?
	    #No point as we have to query anyway
	    $self->store_set_probes_features($achip->dbID(), \%pfs, $ops);
	    throw("ops still defined in caller") if defined $ops;
	  }
	  
	  $ops = Bio::EnsEMBL::Funcgen::ProbeSet->new(
						      -NAME => $data[$hpos{'FEATURE_ID'}],
						      -SIZE => undef,
						      -FAMILY  => $data[$hpos{'CONTAINER'}],
						      #xref_id => $data[$hpos{'SEQ_ID'}],#Need to populate xref table
															   );
					
	  #should we store straight away or build a probeset/probe/feature set, and then store and validate in turn?
	  #Store directly have separate method to validate and update?
	  #would need to check if one exists before storing anyway, else we could potentially duplicate the same probe/probeset from a different array
	  #remember for affy we need duplicate probe records with identical probe ids, probeset records unique across all arrays
	  
	  undef %pfs
	}
	elsif($. > 2){#may have previous ops set, but next has no ops, or maybe just no ops's at all
	  $self->store_set_probes_features($achip->dbID(), \%pfs, $ops);
	  throw("ops still defined in caller") if defined $ops;
	}
	
				
	###PROBES
	#should we cat $xref_id to $probe_id here to generate unique id?
	#would be messy to handle in the code, but would have to somewhere(in the retrieval code)
	
	$length = length($data[$hpos{'PROBE_SEQUENCE'}]);	
	#$probe_string .= "\t${psid}\t".$data[$hpos{'PROBE_ID'}]."\t${length}\t$ac_id\t${class}\n";
	#print "Generating new probe with $array ".$array->dbID()." and ac id ".$ac{'dbID'}."\n";
	
	$op = Bio::EnsEMBL::Funcgen::Probe->new(
						-NAME          => $data[$hpos{'PROBE_ID'}],
						-LENGTH        => $length,
						-ARRAY         => $array,
						-ARRAY_CHIP_ID => $achip->dbID(),
						-CLASS         => $class,
					       );
	
	
	%{$pfs{$data[$hpos{'PROBE_ID'}]}} = (
					     probe => $op,
					     features => [],
					    );
	
				
	###PROBE FEATURES
	#How can we be certain that we have the same mapping in the DB?
	#Put checks in here for build?
	#Need to handle controls/randoms here
	#won't have features but will have results!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	#The format of the pos file looks like it should have all the data required, but 
	# chromsome is missing, first undef :(
	#Need to use $pos here instead of .pos file
	#However have problems defining probe class, as not populated in test set
	#derive from container! :(
	#Ignore controls/random as they won't have a region
	#Also need to handle multiple mappings? 
	#As results reference probes, no features, can have multiple features on different builds
	
	#if($class eq "EXPERIMENTAL"){
	
	#if(exists $regions{$data[$hpos{'SEQ_ID'}]}){
	#  $fid++;
	#  $pf_string .= "\t".$regions{$data[$hpos{'SEQ_ID'}]}{'seq_region_id'}."\t".$data[$hpos{'POSITION'}]."\t".
	#    ($data[$hpos{'POSITION'}] + $length)."\t${strand}\t".$regions{$data[$hpos{'SEQ_ID'}]}{'coord_system_id'}.
	#	"\t${pid}\t${anal_id}\t".$data[$hpos{'MISMATCH'}]."\t${cig_line}\n";
	#$loc .= $regions{$data[$hpos{'SEQ_ID'}]}{'seq_region_id'}.":".$data[$hpos{'POSITION'}].
	#  "-".($data[$hpos{'POSITION'}] + $length).";" if ($self->{'_dump_fasta'});
	# }
	# else{ 
	#   die("No regions defined for ".$data[$hpos{'SEQ_ID'}]." ".$data[$hpos{'PROBE_ID'}].
	#" with family ".$data[$hpos{'CONTAINER'}]);
	#  }
	
	
				
	
	if ($self->{'_dump_fasta'}){
	  (my $chr = $probe_pos{$data[$hpos{'PROBE_ID'}]}->{'chr'}) =~ s/chr//;
	  $loc .= $chr.":".$probe_pos{$data[$hpos{'PROBE_ID'}]}->{'start'}."-".
	    ($probe_pos{$data[$hpos{'PROBE_ID'}]}->{'start'}+ $length).";";
	}
	
	#Hack!!!!!! Still importing probe (and result?)
	if(!  $self->cache_slice($probe_pos{$data[$hpos{'PROBE_ID'}]}->{'chr'})){
	  #warn("Skipping non standard probe (".$data[$hpos{'PROBE_ID'}].") with location:\t$loc\n");
	  next;
	}
	
				
	
	$of = Bio::EnsEMBL::Funcgen::ProbeFeature->new
	  (
	   -START         => $probe_pos{$data[$hpos{'PROBE_ID'}]}->{'start'},
	   -END           =>($probe_pos{$data[$hpos{'PROBE_ID'}]}->{'start'}  + $length),
	   -STRAND        => $strand,
	   -SLICE         => $self->cache_slice($probe_pos{$data[$hpos{'PROBE_ID'}]}->{'chr'}),
	   -ANALYSIS      => $anal,
	   -MISMATCHCOUNT => $data[$hpos{'MISMATCH'}],
	   -PROBE         => undef,#Need to update this in the store method
	  );
				
	
	push @{$pfs{$data[$hpos{'PROBE_ID'}]}{'features'}}, $of;
	
	if($self->{'_dump_fasta'}){			
	  #filter controls/randoms?  Or would it be sensible to see where they map
	  #wrap seq here?
	  $fasta .= ">".$data[$hpos{'PROBE_ID'}]."\t".$data[$hpos{'CHROMOSOME'}].
	    "\t$loc\n".$data[$hpos{'PROBE_SEQUENCE'}]."\n";
	}
      }
      
      #need to store last data here
      $self->store_set_probes_features($achip->dbID(), \%pfs, $ops);
      $self->log(join("\n", @log));
      $achip->adaptor->set_status("IMPORTED", $achip);
      $self->log("ArrayChip:\t".$achip->design_id()." has been IMPORTED");
      
      if ($self->{'_dump_fasta'}){
	print $f_out $fasta if($self->{'_dump_fasta'});
	close($f_out);
      }
    }
    
	
    #Should we build hash of probe_names:probe_feature_ids here for results import
    #Should we dump this as a lookup file for easier recoverability
    #This is the biggest step in the import
    #Building the hash would be fastest but least recoverable
    #Would have to write recover statments for each step i.e. build the most recent data structure required for the next import step
    #Individual queries for each result would take ages
    #This is all assuming there are no random records in the table i.e. ID series is linear with no gaps.
    
    
    #Could throw a random validation check in every X entries?
    #This would only work for non-parallel imports
    #periodic import when hit new probe_set, but no new data printed
		
  }
  
  #    #Need to print last probe_set here only if current and last probeset_id match
  #  if($probe_set{'name'} eq $data[$hpos{'FEATURE_ID'}]){
  #	$ps_string .= "\t".$probe_set{'name'}."\t".$probe_set{'size'}."\t".$probe_set{'family'}."\n";
  #      }
  
  
  #     print $p_out $probe_string;
  #     print $ps_out $ps_string;
  #     print $pf_out $pf_string;
  #     print $f_out $f_string if($self->{'_dump_fasta'});
  #$self->Timer()->mark("End of probe loop");
  
  #    close($fh);
  #    close($ps_out);
  #    close($p_out);
  #    close($pf_out);
  #close($f_out) if ($self->{'_dump_fasta'});
  
  
  
  
  $self->log("Finished parsing probe data");
  #Total probe_sets:\t$psid\n".
  #	       "Total probes:\t$pid\nTotal probe_features:\t$fid");
  
  return;
}


=head2 read_results_data

  Example    : $imp->read_results_data();
  Description: Parses and import raw results
  Returntype : none
  Exceptions : none
  Caller     : Importer
  Status     : Medium - move parts to VendorDefs, should still function as normal

=cut

 
sub read_results_data{
  my $self = shift;
  
  #TO DO
  #Recovery ? read feature_map.tmp into , query DB for last imported result(remove for safety?), restart from that point in file?
  #slurp here may require too much memory?
  #import files and do checks on ids to make sure they've imported properly
  #i.e. select the last entry based on the expected table id.
  
  $self->log("Parsing ".$self->vendor()." results");
  my (@header, @data, @design_ids, %hpos);#%channel_idx, %cc_id
  my ($fh, $pid, $line, $file);
  my $anal = $self->db->get_AnalysisAdaptor->fetch_by_logic_name("RawValue");
  my $result_set = $self->get_import_ResultSet('channel', $anal);

  
  if($result_set){#we have some new data to import
    
    #build channel_idx and chip_channel_id hash after we store
    #foreach my $echip(@{$self->experiment->get_ExperimentalChips()}){
    #  
    #  if( ! $echip->has_status('IMPORTED')){
    #
    #foreach my $chan (@{$echip->get_Channels()}){
    #	  my $field = $echip->unique_id()."_".$self->get_def('dye_freqs')->{$chan->dye()};
    ##	  $channel_idx{$field} = undef;
    #	  $cc_id{$field} = $result_set->get_chip_channel_id($chan->dbID());
    #	}
    #      }else{
    #	$self->log("ExperimentalChip (".$echip->unqiue_id().") already has status:\tIMPORTED");
    #      }
    #    }
    
    
    #foreach my $array(@{$self->arrays()}){
    #  #This should be accessing echips for this experiment, not design_ids.
    #  #we may have only a subset of designs/chips form the whole array
    #  #and we only want to import the ones we have in the context of this experiment
    #  
    #  @design_ids = @{$array->get_design_ids()};
      
    #  foreach my $design_id(@design_ids){
	
    foreach my $echip(@{$self->experiment->get_ExperimentalChips()}){

      if( ! $echip->has_status('IMPORTED')){

	foreach my $chan(@{$echip->get_Channels()}){
	  my ($probe_elem, $score_elem);
	  my $cnt = 0;
	  my $r_string = "";
	  my $chan_name = $echip->unique_id()."_".$self->get_def('dye_freqs')->{$chan->dye()};
	  my $cc_id = $result_set->get_chip_channel_id($chan->dbID());
	  
	  $self->log("Reading results for channel:\t${chan_name}");
	  #my (%tmp, @chan_fields);
	  #my $achip = $array->get_ArrayChip_by_design_id($design_id);
	  
	  #open/backup output
	  my $file = $self->get_dir("raw")."/result.".$chan_name.".txt";	
	  $self->backup_file($file);	#This may cause empty file backups
	  my $r_out = open_file(">", $file);
	  
	  
	  #open/slurp input
	  #$file_name = (scalar(@design_ids) > 1) ? $achip->name() : "All";  
	  (my $alt_chan_name = $chan_name) =~ s/\_/\_1h\_/;
	  my $found = 0;

	  FILE: foreach my $name($chan_name, $alt_chan_name){

	    foreach my $suffix("_pair.txt", ".pair", ".txt"){
	      
	      $file = $self->get_dir("results")."/".$name.$suffix;

	      if(-f $file){
		$found = 1;
		last FILE;
	      }
	    }
	  }

	  throw("Could not find result file for Channel(${chan_name}) in ".$self->get_dir('results')) if ! $found; 


	  $fh = open_file("<", $file);	
	  my @lines = <$fh>;
	  close($fh);
	
	  foreach $line(@lines){
	  
	    #can we preprocess effectively?
	    
	    next if $line =~ /^#/;
	    next if $line =~ /NGS_CONTROLS/;
	    next if $line =~ /V_CODE/;
	    next if $line =~ /H_CODE/;
	    next if $line =~ /RANDOM/;

	    $line =~ s/\r*\n//o;
	    @data = split/\t/o, $line;
	  
	    ###PROCESS HEADER
	    if ($line =~ /PROBE_ID/o){
	      
	      #GENE_EXPR_OPTION	SEQ_ID	PROBE_ID	POSITION	43827_532	43827_635	43837_532	43837_635	46411_532	46411_635	46420_532	46420_635	47505_532	47505_635	47525_532	47525_635
	  
	      
	      
	      #find probe column
	      #for my $i(0..$#header){
		
	      %hpos = %{$self->set_header_hash(\@data, $self->get_def('result_fields'))};

		#if($header[$i] eq "PROBE_ID"){
		#  $probe_elem = $i;
		#  #next;
		#}elsif($header[$i] eq "PM"){#found score field
		#  $score_elem = $i;
		#}	    
		
		#Populate channel idx for directly accessing unstored channel results
		#foreach my $field(keys %channel_idx){
		#$header[$i] =~  s/1h_//;
		#	
		#	if($header[$i] eq $field){
		#	  $channel_idx{$field} = $i;
		#	}
		#}
	      #}
	      
	      ##make tmp hash with remaining(other achip) channel fields
	      #foreach my $field(keys %channel_idx){
	      
	      #  if(! defined $channel_idx{$field}){
	      #	#warn here.
	      #	$tmp{$field} = undef;
	      #	delete $channel_idx{$field};
	      #      }
	      #    }
	      #    @chan_fields = keys %channel_idx;#do this once rather than having to call keys for everyline
	      #$self->log("Parsing results for channels:\t@chan_fields");
	      
	      next;#finished processing header
	    }
	    
	    ###PROCESS DATA
	    #foreach my $field(@chan_fields){
	 

	    #Is this string concat causing the slow down, would it befaster to use an array and print a join?
 	    
	    if($pid = $self->get_probe_id_by_name($data[$hpos{'PROBE_ID'}])){
	       $cnt ++;
	       $r_string .= "\t${pid}\t".$data[$hpos{'PM'}]."\t${cc_id}\n";
	     }else{
	      warn "Found unfiltered non-experimental probe in input $data[$hpos{'PROBE_ID'}]";
	    }
	  
	  
	    ###PRINT SOME RESULTS
	    if($cnt > 10000){
	      $cnt = 0;
	      print $r_out $r_string;
	      $r_string ="";
	      #could we fork here and import in the background?
	    }
	  	
	    #SET REMANING CHANNELS
	    #%channel_idx = %tmp;
	    #@chan_fields = keys %channel_idx;
	    
	    #warn "Finished parsing file";
	  } 
	  #PRINT/CLOSE Channel file
	  print $r_out $r_string;
	  close($r_out);
	  $self->log("Finished parsing $chan_name result");
	}
      }
      #if(%channel_idx){
      #	throw("Not all the Experiment Channel result fields were found\nAbsent channels:\t".(keys %channel_idx));
      #      }
    }
  }else{
    $self->log("Skipping results parse and import");
  }

  $self->log("Finished parsing results");


  return;
}



1;
