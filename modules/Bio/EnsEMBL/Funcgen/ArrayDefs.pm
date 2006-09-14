package Bio::EnsEMBL::Funcgen::ArrayDefs;

use Bio::EnsEMBL::Funcgen::OligoArray;
use Bio::EnsEMBL::Funcgen::OligoProbeSet;
use Bio::EnsEMBL::Funcgen::OligoProbe;
use Bio::EnsEMBL::Funcgen::OligoFeature;
use Bio::EnsEMBL::Funcgen::ExperimentalChip;
use Bio::EnsEMBL::Funcgen::Channel;
use Bio::EnsEMBL::Utils::Exception qw( throw warning );
use Bio::EnsEMBL::Funcgen::Utils::EFGUtils qw(species_chr_num open_file);
use strict;


my $reg = "Bio::EnsEMBL::Registry";

#define Array defs hashes here
#vendor specific paths etc
#method arrays for parse

#Need to use VendorDefs for vendor specific methods?
#We would need to reset them for each Import
#Create internal object and pass ref to Importer
#So importer calls generic methods on Defs object which then uses set ref to access


sub set_defs{
	my ($self) = @_;

	#Is this only valid unless somethign inherits from from Experiment ?
	throw("This is a skeleton class for Bio::EnsEMBL::Experiment, should not be used directly") if(! $self->isa("Bio::EnsEMBL::Funcgen::Importer"));



	#Can we do some funky check on $self to see if it is Experiment
	#Abstract these to ${vendor}Defs.pm modules
	#we can't use methods to build vars :(
	
	#could tidy up below but would disable over-ride for inpur_dir option
	#$self->{'_input_dir'} = $self->get_dir("data")."/native/".$self->vendor()."/".$self->instance();

	my %array_defs = (
					  nimblegen => {
									#order of these data arrays is important!
									array_data   => ["array_chip"],
									probe_data   => ["probe"],
									results_data => ["results"],
									#import_methods  => [],
									#data paths here?
									
									#This input dir is being recast as _input_dit
									input_dir      => $self->get_dir("data")."/raw/".$self->vendor()."/".$self->name(),
									#but this is vendor specific and remains an array_def
									design_dir     => $self->get_dir("data")."/raw/".$self->vendor()."/".$self->name()."/DesignFiles",

									chip_file        => $self->get_dir("data")."/raw/".$self->vendor()."/".$self->name()."/SampleKey.txt",
									array_file       => $self->get_dir("data")."/raw/".$self->vendor()."/".$self->name()."/DesignNotes.txt",
									#feature_map_file => $self->get_dir("import")."/feature_map.tmp",
									results_file     => $self->get_dir("data")."/raw/".$self->vendor()."/".$self->name()."/PairData/All_Pair.txt",

							results_dir     => $self->get_dir("data")."/raw/".$self->vendor()."/".$self->name()."/PairData",

									dye_freqs => {(
												   Cy5 => 635,
												   Cy3 => 532,
												  )},

									#PROBE_DESIGN_ID	CONTAINER	DESIGN_NOTE	SELECTION_CRITERIA	SEQ_ID	PROBE_SEQUENCE	MISMATCH	MATCH_INDEX	FEATURE_ID	ROW_NUM	COL_NUM	PROBE_CLASS	PROBE_ID	POSITION	DESIGN_ID	X	Y
									#header_fields => [


									#this won't work as design_name is undef at this point
									#need to make this a ref so it will update the values?
									#gene_description_file => $self->get_array_def("design_dir")."/".$self->array_data("design_name").".ngd",
								   },
					 );

	#print "Setting ".$self->vendor()." defs to ".$array_defs{$self->vendor()}." ".Data::Dumper::Dumper($array_defs{$self->vendor()})."\n";

	$self->{'array_defs'} = $array_defs{$self->vendor()};
	
	#Set mandatory defs here?



	return;

}



sub get_def{
	my ($self, $data_name) = @_;
	return $self->get_data('array_defs', $data_name);#will this cause undefs?
}





#How are we going to use different formats here?
#This is hardcoded for nimblegen at present
#Nimblgen data does not contain array data as such!

sub read_array_chip_data{
  my ($self) = shift;

  my ($design_desc);
  my $fh = open_file("<", $self->get_def("array_file"));

  #Slurp file to string, sets local delimtter to null and subs new lines
  ($design_desc = do { local ($/); <$fh>;}) =~ s/\r*\n$//;
  close($fh);
  
  #Currently 1 design = 1 chip = 1 array /DVD
  #Different designs are not currently collated into a chip_set/array in any ordered manner
  #Register each design as an array and an array_chip
  #May want to group array_chips into array/chip sets by association though the API
  
  $fh = open_file("<", $self->get_def("chip_file"));
  
  $self->log("Reading chip data");
  
  my ($line, $chip_uid, $dye, $design_name, $design_id, $sample_label,  $sample_desc);
  my ($array, $species, $sample_type, $channel, %tmp, %hpos, @data);
  my $tmp_uid = "FIRST";
  #Need some way of capturing other experimental variables?
  
  
  my ($echip, $ac_id);
  my $oa_adaptor = $self->db->get_OligoArrayAdaptor();
  my $ec_adaptor = $self->db->get_ExperimentalChipAdaptor();
  my $chan_adapt = $self->db->get_ChannelAdaptor();

  warn("Harcoded for one array per experiment\n");

  while ($line = <$fh>){
    $line =~ s/\r*\n//;#chump
    @data =  split/\t/o, $line;

    #ORD_ID  CHIP_ID DYE     DESIGN_NAME     DESIGN_ID       SAMPLE_LABEL    SAMPLE_SPECIES  SAMPLE_DESCRIPTION      TISSUE_TREATMENT        PROMOT_SAMPLE_TYPE
    #ORD_ID  CHIP_ID DYE     DESIGN_NAME     DESIGN_ID       SAMPLE_LABEL    SAMPLE_SPECIES  SAMPLE_DESCRIPTION-1    PROMOT_SAMPLE_TYPE
    if ($. == 1){
      %hpos = %{$self->set_header_hash(\@data)};

      #we need to set the sample description field name, as it can vary :(((
      @data = grep(/SAMPLE_DESCRIPTION/, keys %hpos);
      $sample_desc = $data[0];
      warn("More than one sample description(@data) in ".$self->get_def("chip_file")."\n") if(scalar @data >1);
      next;
    }
    
    
    #Need to handle array class here i.e. two channel arrays will have two lines
    #ignore tissue treatment, too wordy with typos.
    #Need to map sample_label to dye and description?
    #Some duplicate sample labels do not have the same description > different channels
    
    #(undef, $chip_uid, $dye, $design_name, $design_id, $sample_label,
    # $species, $sample_desc, undef, $sample_type) = split/\t/o, $line;

    #validate species here
    #look up alias from registry and match to self->species
    #registry may not be loaded for local installation
    
    if(! defined $array){
      warn("We need to accomodate arrays which have only been partially imported here");
      #Retrieve the array by the array_name, defined by the user or use the design_name
      #$array = $oa_adaptor->fetch_by_name_vendor($design_name, $self->vendor());

      #Then we can update the array chips as we go along

      %tmp = (
	      array_id    => undef,
	      dbID        => undef,
	      name => $data[$hpos{'DESIGN_NAME'}],
	     );
         
      #This is treating each array chip as a separate array, unless arrayset is defined
      #AT present we have no way of differentiating between different array_chips on same array???!!!
      #Need to add functionality afterwards to collate array_chips into single array
      
      
      #store now checks whether already stored and updates array chips accordingly
      $array = Bio::EnsEMBL::Funcgen::OligoArray->new
	(
	 -NAME        => $self->{'array_name'} || $data[$hpos{'DESIGN_NAME'}],
	 -FORMAT      => uc($self->format()),
	 -VENDOR      => uc($self->vendor()),
	 -DESCRIPTION => $design_desc,
	);
      $array->add_array_chip($data[$hpos{'DESIGN_ID'}], \%tmp);
      ($array) = @{$oa_adaptor->store($array)};

        
      #Need to reg array chip here (to differentiate from what is already in DB, handles subset of arrays)
      #This is a registry of array chips which have previously been stored for validation purpose
      #To avoid importing probes twice

    }elsif((! $array->get_array_chip_by_design_id($data[$hpos{'DESIGN_ID'}])) && ($self->{'array_set'})){

      print "generating new ac for same array ".$data[$hpos{'DESIGN_ID'}]."\n";

      %tmp = (
	      array_id    => undef,
	      dbID        => undef,
	      name => $data[$hpos{'DESIGN_NAME'}],
	     );
      $array->add_array_chip($data[$hpos{'DESIGN_ID'}], \%tmp);
      ($array) = @{$oa_adaptor->store($array)};
      
    }elsif(! $array->get_array_chip_by_design_id($data[$hpos{'DESIGN_ID'}])){
      throw("Found experiment with more than one design");
    }
    
    #Parse and Populate ExperimentalChip/Channels
    if ($tmp_uid eq "FIRST" || $data[$hpos{'CHIP_ID'}] != $tmp_uid){
      $tmp_uid = $data[$hpos{'CHIP_ID'}];
     
      $echip =  Bio::EnsEMBL::Funcgen::ExperimentalChip->new
	(
	 -DESIGN_ID      => $data[$hpos{'DESIGN_ID'}],
	 -EXPERIMENT_ID  => $self->experiment->dbID(),
	 -DESCRIPTION    => $data[$hpos{$sample_desc}],
	 -ARRAY_CHIP_ID  => $array->get_array_chip_by_design_id($data[$hpos{'DESIGN_ID'}])->{'dbID'},
	 -UNIQUE_ID      => $data[$hpos{'CHIP_ID'}],
	);
      
      ($echip) = @{$ec_adaptor->store($echip)};

      #should we nest these in the Experiment and 
      #don't need to add them to store, just have method which always retrieves all echips from db
      $self->{'echips'}{$data[$hpos{'CHIP_ID'}]} = $echip;#do we still need this?
      
    }
    
    
    #Handles single/mutli
    $channel =  Bio::EnsEMBL::Funcgen::Channel->new
      (
       -EXPERIMENTAL_CHIP_ID => $echip->dbID(),
       -DYE                  => $data[$hpos{'DYE'}],
       -SAMPLE_LABEL         => $data[$hpos{'SAMPLE_LABEL'}],
       -SPECIES              => $self->species(),#on channel/sample to enable multi-species chip/experiment
       -TYPE                 => uc($data[$hpos{'PROMOT_SAMPLE_TYPE'}]),
      );
    
        
    $self->set_channel($data[$hpos{'CHIP_ID'}], @{$chan_adapt->store($channel)});#do we need o set this anymore? or rename reg_channel?
    
  }
  
  
  #Set the array now we have stored all details
  $self->arrays($array);
  
  close($fh);
  return;
}


sub arrays{
	my ($self) = shift;

	my $module = 'Bio::EnsEMBL::Funcgen::OligoArray';

	if(@_ && ! $_[0]->isa($module)){
		throw("Must supply a $module");
	}
	
	if(@_){
		push @{$self->{'arrays'}}, @_;
	}


	throw("Does not yet support multiple array imports") if(scalar (@{$self->{'arrays'}}) > 1);
	#need to alter read_probe data at the very least

	return $self->{'arrays'};
		
}

sub echip_data{
	my ($self, $design_id, $data_type, $value) = @_;

	throw("Need to specify a design_id  and a data_type") if (! defined $data_type || ! defined $design_id);

	if(defined $value){
		${$self->get_data('echips', $design_id)}{$data_type} = $value;
	}
	else{
		return ${$self->get_data('echips', $design_id)}{$data_type};#will this cause undefs?
	}
}


sub get_echip{
	my ($self, $chip_uid) = @_;

	throw("Need to specify unique_id") if (! defined $chip_uid);

	return $self->{'echips'}{$chip_uid};
}




sub get_channels{
	my ($self, $chip_uid) = @_;
		
	return $self->get_echip($chip_uid)->{'channels'};
}


sub get_channel{
	my ($self, $chan_uid) = @_;
	
	#This may be different for other vendors!!		
	#my ($chip_uid);
	#($chip_uid = $chan_uid) =~ s/_[0-9]*//o;

	return $self->{'channels'}{$chan_uid};
}

#Nest the channels as a hash with key = dye?



#This should be in ExperimentalChip, get_all_Channels?
#defs issue, should this remain in the importer?


#do we need this now, can we not do this dynamically?
#or maybe just store the channel ids, rather than the channel
#then retrieve the relevant channels from the echip channel hash

sub set_channel{
	my ($self, $chip_uid, $channel) = @_;
	
	throw("Need to pass a chip_uid and a Channel object") if(! defined $chip_uid || ! $channel->isa('Bio::EnsEMBL::Funcgen::Channel'));
	$self->{'channels'}{"${chip_uid}_".${$self->get_def("dye_freqs")}{$channel->dye()}} = $channel;

	return;
}

sub channel_data{
	my ($self, $chan_uid, $data_type, $value) = @_;
	

	throw("Need to provide a channel uid and a data_type") if (! defined $chan_uid || ! defined $data_type);

	if(defined $value){
		throw("choudl now use Channel methods \n");
		${$self->get_channel($chan_uid)}{$data_type} = $value;
	}
	else{
		return ${$self->get_channel($chan_uid)}{$data_type};
	}
}



#this handle probe_set, probe, probe_feature
#Assumes no chip_design per experimental set.
sub read_probe_data{
  my ($self) = shift;

  my ($fh, $line, @data, %hpos);
  $self->log("Parsing probe data (".localtime().")...can we do a reduced parse if we know the array chips are already imported?");

  ### Read in
  #ndf file: probe_set, probe and probe_feature(.err contains multiple mappings)
  #pos file: probe locations(counts)
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
  #my $op_a = $self->db->get_OligoProbeAdaptor();
  #my $ops_a = $self->db->get_OligoProbeSetAdaptor();
  #my $opf_a = $self->db->get_OligoProbeFeatureAdaptor();


  my $slice_a = $self->db->get_SliceAdaptor();
  #we need to register a coord_system here, 
  #hardcode for chr for now
  #is this ever reset can this not be outside the loop?
  #do we need to handle multiple arrays which may be mapped to different builds?
  my $cs = $self->db->get_FGCoordSystemAdaptor()->fetch_by_name_schema_build_version('chromosome', $self->db->_get_schema_build($self->db->dnadb()));
  #should really test and store here if not valid, but this will be handled by adaptor store methods
  
  foreach my $array(@{$self->arrays()}){

    foreach my $design_id(@{$array->get_design_ids()}){
      my %ac = %{$array->get_array_chip_by_design_id($design_id)};

   
      #check status of array chip here
      if($self->db->fetch_status_by_name('array_chip', $ac{'dbID'}, 'IMPORTED')){
	warn("Skipping array chip ($design_id) already fully imported\n");
	next;
      }


      #THIS BLOCK DOES NOT ACCOUNT FOR MULTIPLE ARRAYS PROPERLY, WOULD HAVE TO IMPLEMENT ARRAY SPECIFIC CACHES
      #all out files are generic, but are we converting to adaptor store?

=pod

      $fh = open_file("<", $self->get_def("design_dir")."/".$ac{'design_name'}.".ngd");
      my ($start, $stop, %regions, %probe_pos);
	  
      #May not have both ngd and pos file?

      while ($line = <$fh>){
	$line =~ s/\r*\n//;#chump
	@data =  split/\||\t/o, $line;
	
	
	#SEQ_ID	SEQ_UNIQUE|BUILD|CHROMOSOME|LOCATION|DESCRIPTION|DATE_ENTERED|SOURCE_DB
	if ($. == 1){
	  %hpos = %{$self->set_header_hash(\@data)};
	  next;
	}
	
	#What about strand!!!!!!!!!!!
	$data[$hpos{'CHROMOSOME'}] =~ s/chr//;
	($start, $stop) = split/-/o, $data[$hpos{'LOCATION'}];
	
	#Do we need seq_id check here for validity?
	#overkill?
	if(exists $regions{$data[$hpos{'SEQ_ID'}]}){
	  croak("Duplicate regions\n");
	}else{
	  #$data[$hpos{'CHROMOSOME'}] = species_chr_num($self->species(), 	$data[$hpos{'CHROMOSOME'}]);
	  
	  #Set region hash for SEQ_ID
	  #Need to look up seq_region id here for given build
	  #Build should be manually specified as we can't guarantee it will be in the correct format
	  #or present at all
	  
	  $regions{$data[$hpos{'SEQ_ID'}]} = 
	    {
	     start => $start,
	     stop  => $stop,
	     seq_region_id => $self->get_chr_seq_region_id($data[$hpos{'CHROMOSOME'}], $start, $stop),
	     coord_system_id => $cs->dbID(),
	    };
	}
	
      }
      
      close($fh);
      
      
      #ONLY USE THIS FOR VALIDATION OF EXPERIMENTAL PROBES!!!! CAN REMOVE IF PROBE_CLASS POPULATED
      #SLURP PROBE POSITIONS
      $fh = open_file("<", $self->get_def("design_dir")."/".$ac{'design_name'}.".pos");
      
      #don't % = map ! Takes a lot longer than a while ;)
      while($line = <$fh>){
	#$line =~ s/\r*\n//;#Not using last element
	@data =  split/\t/o, $line;
	
	#SEQ_ID	CHROMOSOME	PROBE_ID	POSITION	COUNT
	if ($. == 1){
	  %hpos = %{$self->set_header_hash(\@data)};
	  next;
	}
	#($seq_id, undef, $probe_id, $lstart) = split/\t/o, $line;
	
	#can we remove this?
	throw("Found duplicate mapping for ".$data[$hpos{'PROBE_ID'}]) if(exists $probe_pos{$data[$hpos{'PROBE_ID'}]});
	
	$probe_pos{$data[$hpos{'PROBE_ID'}]} = {(
						 seq_id => $data[$hpos{'SEQ_ID'}],
						 lstart => $data[$hpos{'POSITION'}],
						)};
	
      }
      
      close($fh);
      
=cut
      
  #OPEN PROBE IN/OUT FILES
  

  
  $fh = open_file("<", $self->get_def("design_dir")."/".$ac{'name'}.".ndf");
      #Need to set these paths in each  achip hash, file names could be tablename.chip_id.txt
      #my $p_out = open_file(">", $self->get_dir("import")."/probe.".$ac{'design_name'}."txt");
      #my $ps_out = open_file(">", $self->get_dir("import")."/probe_set.".$ac{'design_name'}.".txt");
      #my $pf_out = open_file(">", $self->get_dir("import")."/probe_feature.".$ac{'design_name'}."txt");
      my $f_out = open_file(">", $self->get_dir("import")."/probe.".$ac{'name'}."fasta")	if($self->{'_dump_fasta'});
      my ($length);
      my ($ops, $op, $of, $chr, @probes, @features, %slices, %pfs);
      my $anal = $self->db->get_AnalysisAdaptor()->fetch_by_logic_name("VendorMap");
      ##HARDCODED!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      my $strand = 0;
      my $cig_line = "50M";
      my $f_string = "";
      #$self->Timer()->mark("Starting probe loop");
      
      while($line = <$fh>){
	$line =~ s/\r*\n//;
	@data =  split/\t/o, $line;
	my $loc = "";
	my $class = "EXPERIMENTAL";
	
	#PROBE_DESIGN_ID	CONTAINER	DESIGN_NOTE	SELECTION_CRITERIA	SEQ_ID	PROBE_SEQUENCE	MISMATCH	MATCH_INDEX	FEATURE_ID	ROW_NUM	COL_NUM	PROBE_CLASS	PROBE_ID	POSITION	DESIGN_ID	X	Y
	#2067_0025_0001  BLOCK1          0       chrX    TTAGTTTAAAATAAACAAAAAGATACTCTCTGGTTATTAAATCAATTTCT      0       52822449        52822449        1       25      experimental    chrXP10404896   10404896        2067    25      1

	if ($. == 1){	
	  %hpos = %{$self->set_header_hash(\@data)};
	  next;
	}

	
	if($data[$hpos{'CONTAINER'}] =~ /control/io){
	  $class = "CONTROL";
	}
	elsif($data[$hpos{'CONTAINER'}] =~ /random/io){
	  $class = "RANDOM";
	}
	elsif($data[$hpos{'PROBE_CLASS'}] !~ /experimental/io){
	  $class = "OTHER";
	}
	#elsif(! exists $probe_pos{$data[$hpos{'PROBE_ID'}]}){	#HACKY HACKY HACK HACK!! Needed for valid region retrival
	#  $class = "OTHER";
	#}
	
	
	#This assumes all probes in feature/probeset are next to each other!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	#Which non-EXPERIMENTAL probes do we want to store??!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	next if ($class ne "EXPERIMENTAL");#SPIKE?  Cannot have features tho' or would display

	if($data[$hpos{'FEATURE_ID'}] != $data[$hpos{'MATCH_INDEX'}]){#Probe set data
	  #print "Generating new probeset:\tFeature id:\t".$data[$hpos{'FEATURE_ID'}]."\tmatchindex:\t".$data[$hpos{'MATCH_INDEX'}]."\n";
	  
	  if($ops && ($data[$hpos{'FEATURE_ID'}] ne $ops->name())){
	    #THis is where we chose to update/validate
	    #Do we need to pass probes if they're already stored..may aswell to reduce mysql load?
	    #No point as we have to query anyway
	    #$self->store_set_probes_features($ac{'dbID'}, $ops, \@probes, \@features);
	    $self->store_set_probes_features($ac{'dbID'}, $ops, \%pfs);
	    throw("ops still defined in caller") if defined $ops;
	  }

       	  $ops = Bio::EnsEMBL::Funcgen::OligoProbeSet->new(
							   -NAME => $data[$hpos{'FEATURE_ID'}],
							   -SIZE => undef,
							   -FAMILY  => $data[$hpos{'CONTAINER'}],
							   #xref_id => $data[$hpos{'SEQ_ID'}],#Need to populate xref table
						       );

	  #should we store straight away or beuilf a probeset/probe/feature set, and then store and validate in turn?
	  #Store directly have separate method to validate and update?
	  #would need to check if one exists before storing anyway, else we could potentially duplicate the same probe/probeset from a different array
	  #remember for affy we need duplicate probe records with identical probe ids, probeset records unique across all arrays

	  #undef @probes;
	  undef %pfs
	}
	#elsif($ops){#Got ops, but new probe does not have ops, need to store old ops
	elsif($. > 2){#may have previous ops set, but next has no ops, or maybe just no ops's at all
	  #$self->store_set_probes_features($ac{'dbID'}, $ops, \@probes, \@features);
	  $self->store_set_probes_features($ac{'dbID'}, $ops, \%pfs);

	  throw("ops still defined in caller") if defined $ops;
					   
	  #undef $ops;
	  #undef @probes;
	}
	
	
	###PROBES
	#should we cat $xref_id to $probe_id here to generate unique id?
	#would be messy to handle in the code, but would have to somewhere(in the retrieval code)
      
	$length = length($data[$hpos{'PROBE_SEQUENCE'}]);	
	#$probe_string .= "\t${psid}\t".$data[$hpos{'PROBE_ID'}]."\t${length}\t$ac_id\t${class}\n";
	#print "Generating new probe with $array ".$array->dbID()." and ac id ".$ac{'dbID'}."\n";
      
	$op = Bio::EnsEMBL::Funcgen::OligoProbe->new(
						     -NAME          => $data[$hpos{'PROBE_ID'}],
						     -LENGTH        => $length,
						     -ARRAY         => $array,
						     -ARRAY_CHIP_ID => $ac{'dbID'},
						     -CLASS         => $class,
						    );
      
	
      	%{$pfs{$data[$hpos{'PROBE_ID'}]}} = (
					     probe => $op,
					     features => [],
					    );

		
	

	#push @probes, $op;
      
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
	

	#build slice hash here
	($chr = $data[$hpos{'SEQ_ID'}]) =~ s/chr//;
	if(! exists $slices{$chr}){
	  $slices{$chr} = $slice_a->fetch_by_region('chromosome', $chr);
	}
	
	
	$loc .= $chr.":".$data[$hpos{'POSITION'}].
	  "-".($data[$hpos{'POSITION'}] + $length).";" if ($self->{'_dump_fasta'});
	
	
	#Hack!!!!!!
	if(!  $slices{$chr}){
	  #warn("Skipping non standard probe (".$data[$hpos{'PROBE_ID'}].") with location:\t$loc\n");
	  #pop @probes;
	  next;
	}
	
	
	
	$of = Bio::EnsEMBL::Funcgen::OligoFeature->new(
						       -START         => $data[$hpos{'POSITION'}],
						       -END           =>($data[$hpos{'POSITION'}] + $length),
						       -STRAND        => $strand,
						       -SLICE         => $slices{$chr},
						       -ANALYSIS      => $anal,
						       -MISMATCHCOUNT => $data[$hpos{'MISMATCH'}],
						       -PROBE         => undef,#Need to update this in the store method
						      );
	

	#print "Pushing feature $of\n";
	
	#push @features, $of;
	push @{$pfs{$data[$hpos{'PROBE_ID'}]}{'features'}}, $of;
	
	
	  #}
	
	
	if($self->{'_dump_fasta'}){			
	  #filter controls/randoms?  Or would it be sensible to see where they map
	  #wrap seq here?
	  $f_string .= ">".$data[$hpos{'PROBE_ID'}]."\t".$data[$hpos{'SEQ_ID'}]."\t$loc\n".$data[$hpos{'PROBE_SEQUENCE'}]."\n";
	}
      #}
      }
    
    #need to store last data here
      #$self->store_set_probes_features($ac{'dbID'}, $ops, \@probes, \@features);
      $self->store_set_probes_features($ac{'dbID'}, $ops, \%pfs);
    $self->db->set_status('array_chip', $ac{'dbID'}, "IMPORTED");
    
    if ($self->{'_dump_fasta'}){
      print $f_out $f_string if($self->{'_dump_fasta'});
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

    
  
  
  #$self->log("Finished parsing probe data\nTotal probe_sets:\t$psid\n".
  #	       "Total probes:\t$pid\nTotal probe_features:\t$fid");
  
  return;
}


sub store_set_probes_features{
  #my ($self, $ac_id, $ops, $probes, $features) = @_;

  my ($self, $ac_id, $ops, $pf_hash) = @_;
  

  #just call these directly rather than setting each time?
  #my $op_a = $self->db->get_OligoProbeAdaptor();
  #my $opf_a = $self->db->get_OligoProbeFeatureAdaptor();


  #if(scalar(@$probes) != scalar(@$features)){
  #  my $t = scalar(@$probes);
  #  my $f = scalar(@$features);
  #  throw("Does not accomodate multiple features($f) per probe $t");
  #}

  #Maybe do some validation here against all probes for probeset and ac_id? 

  if($ops){
    $ops->size(scalar(keys %$pf_hash));
    ($ops) = $self->db->get_OligoProbeSetAdaptor->store($ops);
  }



  #If we're going to validate fully, we need to check for probes in this probeset on this array chip
  #Update size if we have any new probes
  #Overkill? Only do on recover? Do not read if array chip is IMPORTED 

  #This does not make any attempt to validate probes/set vs previously stored data


   
  for my $probe_id(keys %$pf_hash){
    
    #set probeset in probe and store
    #the process corresponding feature
    my $probe = $pf_hash->{$probe_id}->{'probe'};
    $probe->probeset($ops) if $ops;
    ($probe) = @{$self->db->get_OligoProbeAdaptor->store($probe)};

        
    foreach my $feature(@{$pf_hash->{$probe_id}->{'features'}}){
      $feature->probe($probe);
      ($feature) = @{$self->db->get_OligoFeatureAdaptor->store($feature)};
      
      
      #or can we keep probes?  Too much memmory?
      #do we need to key 
      #Can't use get_all_Arrays here as we can't guarantee this will only ever be the array we've generated
      #Might dynamically load array if non-present
      #This is allowing multiple dbIDs per probe???  Is this wrong?
      push @{$self->{'_probe_map'}{$probe->get_probename()}}, $probe->dbID();
    }
  }

  undef $ops;#Will this persist in the caller?
  undef %{$pf_hash};

# undef @$probes;
#  undef @$features,

  return;
}




#Remove array element to this?
sub get_probe_ids_by_name{
  my ($self, $name) = @_;

  my (@op_ids);

  #Should only ever be one pid per probe per array per n array_chips
  #i.e. duplicate records per array chip with same pid

  if(defined $self->{'_probe_map'}){
    @op_ids = @{$self->{'_probe_map'}->{$name}};
  }
  else{#get from db
    my $op= $self->db->get_OligoProbeAdaptor->fetch_by_array_probe_probeset_name($self->arrays->[0]->name(), $name);
    
    print "Got probe $op with dbid ".$op->dbID()."\n";
    push @op_ids, $op->dbID();

  }

  return \@op_ids;
}

sub read_results_data{
  my $self = shift;
  
  #TO DO
  #Recovery ? read feature_map.tmp into , query DB for last imported result(remove for safety?), restart from that point in file?
  #slurp here may require too much memory?
  #import files and do checks on ids to make sure they've imported properly
  #i.e. select the last entry based on the expected table id.
  
  $self->log("Parsing results(".localtime().")...");
  
   
  my ($i, $fh, $tmp, $line, $probe_elem, $first_result, $file_name, @header, @data, @design_ids);
  my $r_string = "";
  my $anal = $self->db->get_AnalysisAdaptor->fetch_by_logic_name("RawValue");
  #should check here if defined
  
  my $anal_id = $anal->dbID();
  
  foreach my $array(@{$self->arrays()}){

    @design_ids = @{$array->get_design_ids()};
   
    foreach my $design_id(@design_ids){
      my %ac = %{$array->get_array_chip_by_design_id($design_id)};
      $file_name = (scalar(@design_ids) > 1) ? $ac{'name'} : "All";  
      $fh = open_file("<", $self->get_def("results_dir")."/".$file_name."_pair.txt");

      while($line = <$fh>){
	$line =~ s/\r*\n//;
    
	if ($. == 1){
	  #GENE_EXPR_OPTION	SEQ_ID	PROBE_ID	POSITION	43827_532	43827_635	43837_532	43837_635	46411_532	46411_635	46420_532	46420_635	47505_532	47505_635	47525_532	47525_635
	  @header = split/\t/o, $line;
	  
	  for $i(0..$#header){
	    $probe_elem = $i if($header[$i] eq "PROBE_ID");
	    
	    if($header[$i] =~ /[0-9]+_[0-9]+/){
	      $first_result = $i;
	      last;#assume all fields after first result, are result feilds
	    }
	  }
	  next;
	}
	
	@data = split/\t/o, $line;
	
	#could validate here against reulsts vs. number of channels
	for $i($first_result..$#data){
	  
	  #multiple mappings????????????????????????????????????SHOULD ONLY HAVE ONE PID PER UNIQUE PROBE!!!!!
	  foreach my $pid(@{$self->get_probe_ids_by_name($data[$probe_elem])}){

	    #Need to change this get_channel call?
	    ($tmp = $header[$i]) =~ s/1h_//;

	    $r_string .= "\t${pid}\t".$data[$i]."\t$anal_id\t".$self->get_channel($tmp)->dbID()."\tchannel\n";
	  }
	}
      }
      
      my $r_out = open_file(">", $self->get_dir("import")."/result.".$ac{'name'}.".txt");
      print $r_out $r_string;
      close($r_out);
    }
  }

  $self->log("Finished parsing and results");
  return;
}



1;
