#
# EnsEMBL module for Bio::EnsEMBL::Funcgen::ArrayDefs
#

=head1 NAME

Bio::EnsEMBL::Funcgen::ArrayDefs

=head1 SYNOPSIS

  my $imp = Bio::EnsEMBL::Funcgen::Importer->new(%params);

  $imp->set_defs();


=head1 DESCRIPTION

This is a definitions class which should not be instatiated directly, it 
normally inherited from the Importer.  ArrayDefs contains meta data and methods specific 
to NimbleGen arrays to aid parsing and importing of experimental data.

This module is currently NimbleGen specific but is likely to be genericised even more, 
extracting vendor specific data and methods to a separate VendorDefs module e.g. NimbleGenDefs.


=head1 AUTHOR

This module was written by Nathan Johnson.

=head1 CONTACT

Post questions to the EnsEMBL development list ensembl-dev@ebi.ac.uk

=head1 METHODS

=cut

package Bio::EnsEMBL::Funcgen::ArrayDefs;

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
use strict;


my $reg = "Bio::EnsEMBL::Registry";


#Create internal object and pass ref to Importer
#So importer calls generic methods on Defs object which then uses set ref to access
#Do we need to "set" this or can we just have all of this in $self?

=head2 set_defs

  Example    : $imp->set_defs();
  Description: Sets a definitions hash for the vendor.
  Returntype : none
  Exceptions : Throws if not called by Importer
  Caller     : Importer
  Status     : At risk - likely to be replaced/moved with/to VendorDefs

=cut

sub set_defs{
	my ($self) = @_;

	throw("This is a skeleton class for Bio::EnsEMBL::Importer, should not be used directly") if(! $self->isa("Bio::EnsEMBL::Funcgen::Importer"));
	

	#Need to separate this into vendor defs and experiment defs?
	
	#these cause undefs if we're just setting sanger defs without exp context
	my %array_defs = (
			  NIMBLEGEN => {
					#order of these data arrays is important!
					array_data   => ["array_chip"],
					probe_data   => ["probe"],
					results_data => ["results"],
					sample_key_fields => ['DESIGN_ID', 'DESIGN_NAME', 'CHIP_ID', 'DYE',
							      'SAMPLE_DESCRIPTION', 'PROMOT_SAMPLE_TYPE'],
					ndf_fields      => ['CONTAINER', 'PROBE_SEQUENCE', 'MISMATCH', 'FEATURE_ID', 'PROBE_ID'],
					pos_fields      => ['CHROMOSOME', 'PROBE_ID', 'POSITION', 'COUNT'],
					result_fields   => ['PROBE_ID', 'PM'],
					#import_methods  => [],
					#data paths here?
													
					#but this is vendor specific and remains an array_def
					design_dir     => $self->get_dir("data")."/input/".$self->vendor()."/".$self->name()."/DesignFiles",
					chip_file        => $self->get_dir("data")."/input/".$self->vendor()."/".$self->name()."/SampleKey.txt",
					array_file       => $self->get_dir("data")."/input/".$self->vendor()."/".$self->name()."/DesignNotes.txt",
					results_file     => $self->get_dir("data")."/input/".$self->vendor()."/".$self->name()."/PairData/All_Pair.txt",
					
					results_dir     => $self->get_dir("data")."/input/".$self->vendor()."/".$self->name()."/PairData",
					norm_method => 'VSN_GLOG',
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


			  SANGER => {
	#order of these data arrays is important!
				     array_data   => [],#["array_chip"],
				     probe_data   => ["sanger_array_probe"],
				     results_data => ["sanger_result"],
				     #import_methods  => [],
				     #data paths here?
				     norm_method => undef,
				     #is this disabling -input_dir override option?
				    },


			 );
	
	
	$self->{'array_defs'} = $array_defs{$self->vendor()};

	#warn "Input dir is ".$self->input_dir()."\n";
	
	#Set mandatory defs here?
	return;
}

=head2 get_def

  Arg [1]    : mandatory - name of the data element to retrieve from the defs hash
  Example    : %dye_freqs = %{$imp->get_def('dye_freqs')};
  Description: returns data from the definitions hash
  Returntype : various
  Exceptions : none
  Caller     : Importer
  Status     : Medium - throw if no data_name?

=cut


sub get_def{
  my ($self, $data_name) = @_;
  return $self->get_data('array_defs', $data_name);#will this cause undefs?
}



=head2 read_array_chip_data

  Example    : $imp->read_array_chip_data();
  Description: Parses and imports array & experimental chip meta data/objects
  Returntype : none
  Exceptions : throws if more than one array/design found and not an "array set"
  Caller     : Importer
  Status     : At risk - likely to be replaced/moved with/to VendorDefs

=cut

sub read_array_chip_data{
  my $self = shift;

  my ($design_desc, $line, $tmp_uid, $array, $channel, $array_chip, $echip);
  my ($sample_desc, %hpos, @data);
  my $oa_adaptor = $self->db->get_ArrayAdaptor();
  my $ec_adaptor = $self->db->get_ExperimentalChipAdaptor();
  my $chan_adaptor = $self->db->get_ChannelAdaptor();
  my $ac_adaptor = $self->db->get_ArrayChipAdaptor();

  #Slurp file to string, sets local delimtter to null and subs new lines
  my $fh = open_file("<", $self->get_def("array_file"));
  ($design_desc = do { local ($/); <$fh>;}) =~ s/\r*\n$//;
  close($fh);
  
  #Currently 1 design = 1 chip = 1 array /DVD
  #Different designs are not currently collated into a chip_set/array in any ordered manner
  #Register each design as an array and an array_chip
  #May want to group array_chips into array/chip sets by association though the API
  
 
  warn("Harcoded for one array(can have multiple chips from the same array) per experiment\n");
  $fh = open_file("<", $self->get_def("chip_file"));
  $self->log("Reading chip data");

  warn("We need to change Array/ArrayChip retrieval to ensure that we have IMPORTED status, so avoid having an incomplete arraychip");

  while ($line = <$fh>){
    $line =~ s/\r*\n//;#chump
    @data =  split/\t/o, $line;

    #ORD_ID  CHIP_ID DYE     DESIGN_NAME     DESIGN_ID       SAMPLE_LABEL    SAMPLE_SPECIES  SAMPLE_DESCRIPTION      TISSUE_TREATMENT        PROMOT_SAMPLE_TYPE

    if ($. == 1){
      %hpos = %{$self->set_header_hash(\@data, $self->get_def('sample_key_fields'))};

      #we need to set the sample description field name, as it can vary :(((
      @data = grep(/SAMPLE_DESCRIPTION/, keys %hpos);
      $sample_desc = $data[0];
      warn("More than one sample description(@data) in ".$self->get_def("chip_file")."\n") if(scalar @data >1);
      next;
    }
    
    #Need to handle array class here i.e. two channel arrays will have two lines
    #validate species here
    #look up alias from registry and match to self->species
    #registry may not be loaded for local installation
    
	if(! defined $array){
      #This is treating each array chip as a separate array, unless arrayset is defined
      #AT present we have no way of differentiating between different array_chips on same array???!!!
      #Need to add functionality afterwards to collate array_chips into single array
           
      #This will use a stored array if present
      $array = Bio::EnsEMBL::Funcgen::Array->new
		(
		 -NAME        => $self->{'array_name'} || $data[$hpos{'DESIGN_NAME'}],
		 -FORMAT      => uc($self->format()),
		 -VENDOR      => uc($self->vendor()),
		 -TYPE        => 'OLIGO',
		 -DESCRIPTION => $design_desc,
		);

      ($array) = @{$oa_adaptor->store($array)};


      $array_chip = Bio::EnsEMBL::Funcgen::ArrayChip->new(
							  -ARRAY_ID  => $array->dbID(),
							  -NAME      => $data[$hpos{'DESIGN_NAME'}],
							  -DESIGN_ID => $data[$hpos{'DESIGN_ID'}],
							 );

      #This will use a stored array_chip if present
      ($array_chip) = @{$ac_adaptor->store($array_chip)};
      $array->add_ArrayChip($array_chip);
        
    }
    elsif((! $array->get_ArrayChip_by_design_id($data[$hpos{'DESIGN_ID'}])) && ($self->{'array_set'})){

      $self->log("Generating new ac for same array ".$data[$hpos{'DESIGN_ID'}]."\n");

      $array_chip = Bio::EnsEMBL::Funcgen::ArrayChip->new(
							  -ARRAY_ID  => $array->dbID(),
							  -NAME        => $data[$hpos{'DESIGN_NAME'}],
							  -DESIGN_ID => $data[$hpos{'DESIGN_ID'}],
							 );
   
      ($array_chip) = @{$ac_adaptor->store($array_chip)};
      $array->add_ArrayChip($array_chip);
      
    }
    elsif(! $array->get_ArrayChip_by_design_id($data[$hpos{'DESIGN_ID'}])){
      throw("Found experiment with more than one design");
    }
    
    #Parse and Populate ExperimentalChip/Channels
    if ((! $tmp_uid) || $data[$hpos{'CHIP_ID'}] != $tmp_uid){

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
	   -ARRAY_CHIP_ID  => $array->get_ArrayChip_by_design_id($data[$hpos{'DESIGN_ID'}])->dbID(),
	   -UNIQUE_ID      => $data[$hpos{'CHIP_ID'}],
	  );
      
	($echip) = @{$ec_adaptor->store($echip)};
	$self->experiment->add_ExperimentalChip($echip); #need a contains method here, #need to do this as the probe cache init 'experimental_chips'
      }
    }
   

 
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
  $self->add_array($array);

  close($fh);

  return;
}


=head2 add_array

  Arg [1]    : Bio::EnsEMBL::Funcgen::Array
  Example    : $self->add_array($array);
  Description: Setter for array elements
  Returntype : none
  Exceptions : throws if passed non Array or if more than one Array set
  Caller     : Importer
  Status     : Medium - Remove/Implement multiple arrays?

=cut

sub add_array{
  my $self = shift;

  warn("Move add_array to Experiment");


  #do we need to check if stored?
  if(! $_[0]->isa('Bio::EnsEMBL::Funcgen::Array')){
    throw("Must supply a Bio::EnsEMBL::Funcgen::Array");
  }elsif(@_){
    push @{$self->{'arrays'}}, @_;
  }
   
  throw("Does not yet support multiple array imports") if(scalar (@{$self->{'arrays'}}) > 1);
  #need to alter read_probe data at the very least
  
  return;
}


sub arrays{
  my $self = shift;
  return $self->{'arrays'};
}





sub read_sanger_array_probe_data{
  my ($self, $array_file) = @_;

  $array_file||= $self->array_file();
  my ($line, $fh, @list, $array_file_format, $cmd);
  my ($op, $of, $imported, $fimported);
  my $oa_adaptor = $self->db->get_ArrayAdaptor();
  my $op_adaptor = $self->db->get_ProbeAdaptor();
  my $of_adaptor = $self->db->get_ProbeFeatureAdaptor();
  my $ec_adaptor = $self->db->get_ExperimentalChipAdaptor();
  my $ac_adaptor = $self->db->get_ArrayChipAdaptor();
  my $slice_adaptor = $self->db->get_SliceAdaptor();
  my $analysis = $self->db->get_AnalysisAdaptor->fetch_by_logic_name("SangerPCR")->dbID();
  #have LiftOver? Could then use liftover in  pipeline to redo mappings

 #store now checks whether already stored and updates array chips accordingly
  my $array = Bio::EnsEMBL::Funcgen::Array->new
    (
     -NAME        => $self->{'array_name'},
     -FORMAT      => uc($self->format()),
     -VENDOR      => uc($self->vendor()),
     -TYPE        => 'PCR',
     -DESCRIPTION => "Sanger ENCODE PCR array 3.1.1",
    );

  ($array) = @{$oa_adaptor->store($array)};
         
  #This is treating each array chip as a separate array, unless arrayset is defined
  #AT present we have no way of differentiating between different array_chips on same array???!!!
  #Need to add functionality afterwards to collate array_chips into single array
  my $array_chip = Bio::EnsEMBL::Funcgen::ArrayChip->new(
							 -NAME      => $array->name(),
							 -DESIGN_ID => $array->name(),
							 -ARRAY_ID  =>$array->dbID(),
							);

  ($array_chip) = @{$ac_adaptor->store($array_chip)};
  $array->add_ArrayChip($array_chip);
  $self->add_array($array);


  #we also need to test wether the array as been imported as well as the mappings
  #THis needs to use coord_sys-id not schema_build!!  Duplcaite entries for different schema_builds 
  #with same assembly

  my $dnadb_cs = $self->db->dnadb->get_CoordSystemAdaptor->fetch_by_name('chromosome');
  my $fg_cs = $self->db->get_FGCoordSystemAdaptor->validate_coord_system($dnadb_cs);


  #This fails if we're pointing to an old DB during the release cycle.  Will be fine if we manage to cs mapping dynamically


  if($array_chip->has_status('IMPORTED')){
    $imported = 1;
    warn("Skipping ArrayChip probe import (".$array_chip->name().") already fully imported");
  }elsif($self->recovery()){
    $self->log("Rolling back partially imported ArrayChip:\t".$array_chip->name());
    $self->db->rollback_ArrayChip($array_chip);#This should really remove all CS imports too?
  }


  #should never really have CS imports if not IMPORTED
  #there is however the potential to trash a lot of data if we were to remove the CS importes by mistake
  #do we need to check whether any other sets are using the data?
  #we have to check for result using relevant cs_id and cc_id
  #no removal of probes is the key thing here as nothing is dependent on the feature_ids
  #get all result sets by array chip?  or get all ExperimentalChips by array chip
  #would have to be result set as we would find our own ecs.  May find our own rset
  
  

  if($array_chip->has_status('IMPORTED_CS_'.$fg_cs->dbID())){
    $fimported = 1;
    warn("Skipping ArrayChip feature import (".$array_chip->name().") already fully imported for ".$self->data_version());
  }elsif($self->recovery()){
    $self->log("Rolling back partially imported ArrayChip features:\t".$array_chip->name());
    $self->db->rollback_ArrayChip_features($array_chip, $fg_cs);
  }


  #need to check whether already imported on specified schema_build
  #check for appropriate file given format in input dir or take path

  if(! $fimported){

    if(! $array_file){

      if(! defined $self->get_dir('input')){
	throw("No input_dir defined, if you are running in a non Experiment context please use -array_file");
      }
      
      #hacky ..do better?
      for my $suffix("gff", "adf"){
	$cmd = $self->get_dir('input')."/".$self->{'array_name'}."*".$suffix;
	@list = `ls $cmd`;
	
	if((scalar(@list) == 1) && 
	   ($list[0] !~ /No such file or directory/o)){###this is only printed to STDERR?
	  
	  if(! defined $array_file){
	    $array_file = $list[0];
	  }else{
	    throw("Found more than one array file : $array_file\t$list[0]\nSpecify one with -array_file");
	  }
	}
      }
      
      throw("Cannot find array file. Specify one with -array_file") if (! defined $array_file);
    }
    
    
    if($array_file =~ /gff/io){
      $array_file_format = "GFF";
    }elsif($array_file =~ /adf/io){
      $array_file_format = "ADF";
      throw("Does not yet accomodate Sanger adf format");
    }else{
      throw("Could not determine array file format: $array_file");
    }

    my $fanal = $self->db->get_AnalysisAdaptor->fetch_by_logic_name(($array_file_format eq "ADF") ? "VendorMap" : "LiftOver");
   
    $self->log("Parsing ".$self->vendor()." array data (".localtime().")");
    warn("Parsing ".$self->vendor()." array data (".localtime().")");
    $fh = open_file("<", $array_file);
    

    my ($chr, $start, $end, $strand, $pid);

    while($line = <$fh>){
      $line =~ s/\r*\n//;
      
      #($chr, $start, $end, $ratio, $pid) = split/\t/o, $line;
      ($chr, undef, undef, $start, $end, undef, $strand, undef, $pid) = split/\t|\;/o, $line;
      $pid =~ s/reporter_id=//o;
      $chr  =~ s/chr//;
      $strand = ($strand eq "+") ? 0 : 1;

      #Hack!!!!!!  This is still maintaining the probe entry (and result?)
      if(!  $self->cache_slice($chr)){
	warn("Skipping non standard probe (".$pid.") with location:\t${chr}:${start}-${end}\n");
	next;
      }


      #need to parse dependant on file format 
      #also need to account for duplicate probes on grid
      
      if(! $self->get_probe_id_by_name($pid)){


	if(! $imported){

	  #when we utilise array coords, we need to look up probe cache and store again with new coords
	  #we're currently storing duplicates i.e. different ids with for same probe
	  #when we should be storing two records for the same probe/id
	  #the criteria for this will be different for each vendor, may have to check container etc for NimbleGen
	  
	  #$length = $start - $end;
	  #warn "length is $length";
	  
	  $op = Bio::EnsEMBL::Funcgen::Probe->new(
						  -NAME          => $pid,
						  -LENGTH        => ($end - $start),
						  -ARRAY         => $array,
						  -ARRAY_CHIP_ID => $array->get_ArrayChip_by_design_id($array->name())->dbID(),
						  -CLASS         => 'EXPERIMENTAL',
						 );
	  ($op) = @{$op_adaptor->store($op)};
	  $self->cache_name_id($op->get_probename(), $op->dbID);
	  
	}
      

	$of = Bio::EnsEMBL::Funcgen::ProbeFeature->new(
						       -START         => $start,
						       -END           => $end,
						       -STRAND        => $strand,
						       -SLICE         => $self->cache_slice($chr),
						       -ANALYSIS      => $fanal,
						       -MISMATCHCOUNT => 0,
						       -PROBE_ID     => $self->get_probe_id_by_name($pid),#work around to avoid cacheing probes
						      );
	
	$of_adaptor->store($of);
	
      }else{
	#warn "Need to accomdate duplicate probes here $pid";　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　
	#generate probe and add an x y position and store coordindates
      }
    }
    $array_chip->adaptor->set_status('IMPORTED_CS_'.$fg_cs->dbID(), $array_chip);
    $self->log("ArrayChip:\t".$array_chip->design_id()." has been IMPORTED_CS_".$fg_cs->dbID());
  }

  

  if(! $imported){
    $array_chip->adaptor->set_status('IMPORTED', $array_chip);
    $self->log("ArrayChip:\t".$array_chip->design_id()." has been IMPORTED");
  }
  
  $self->log("Finished parsing ".$self->vendor()." array/probe data (".localtime().")");
  #warn("Finished parsing ".$self->vendor()." array/probe data (".localtime().")");
  
  return;
}

=head2 read_sanger_result_data

  Example    : $imp->read_sanger_result_data();
  Description: Parses and imports result for the sanger PCR array platform
  Returntype : none
  Exceptions : none
  Caller     : Importer
  Status     : At risk - Move parts to "Vendor"Defs.pm, should function the same

=cut

sub read_sanger_result_data{
  my $self = shift;

  #change this to read_gff_chip_results
  #as opposed to gff channel results
  #This should also use the default logic names for the Vendor, or take a user defined list 
  $self->log("Parsing ".$self->vendor()." result data (".localtime().")");

  my ($file, $chip_uid, $line, $echip);
  my ($ratio, $pid, %chip_files);
  my $of_adaptor = $self->db->get_ProbeFeatureAdaptor();
  my $ec_adaptor = $self->db->get_ExperimentalChipAdaptor();
  my $chan_adaptor = $self->db->get_ChannelAdaptor();
  my $analysis = $self->db->get_AnalysisAdaptor->fetch_by_logic_name("SangerPCR");
  my $result_adaptor = $self->db->get_ResultSetAdaptor();
  #this is done to avoid having to self->array_name in loop, will make multiple array loop easier 
  my $array = ${$self->arrays()}[0];

  #This works a little differently as we're not parsing a meta file
  #so the echips haven't been added yet.
  #This is treating each array chip as a separate array, unless arrayset is defined
  #AT present we have no way of differentiating between different array_chips on same array???!!!
  #Need to add functionality afterwards to collate array_chips into single array

  #First add the echips to the Experiment
  
  if(! $self->result_files()){
    my $list = "ls ".$self->input_dir().'/[0-9]*-[0-9]*\.all\.*';
    my @rfiles = `$list`;
    $self->result_files(\@rfiles);
  }

  
  foreach $file(@{$self->result_files()}){
    chomp $file;
    $self->log("Found SANGER results file:\t$file");

    ($chip_uid = $file) =~ s/.*\///;
    $chip_uid =~ s/\..*//;
    $chip_files{$chip_uid} = $file;
    

    $echip = $ec_adaptor->fetch_by_unique_id_vendor($chip_uid, 'SANGER');

    #this should throw if not recovery
    #Nee to check Nimbelgen methods

    if($echip){
  
      if(! $self->recovery()){
	throw("ExperimentalChip(".$echip->unqiue_id().") already exists in the database\nMaybe you want to recover?");
      }
    }else{

      $echip =  Bio::EnsEMBL::Funcgen::ExperimentalChip->new
	(
	 -EXPERIMENT_ID  => $self->experiment->dbID(),
	 -ARRAY_CHIP_ID  => $self->arrays->[0]->get_ArrayChip_by_design_id($array->name())->dbID(),
	 -UNIQUE_ID      => $chip_uid,
	);
          
      ($echip) = @{$ec_adaptor->store($echip)};	
      $self->experiment->add_ExperimentalChip($echip); #if we need a contains method in  here , always add!!
    }

   
    #sub this passing the echip?
    foreach my $type('DUMMY_TOTAL', 'DUMMY_EXPERIMENTAL'){

      my $channel = $chan_adaptor->fetch_by_type_experimental_chip_id($type, $echip->dbID());
      
      if($channel){
	if(! $self->recovery()){
	  throw("Channel(".$echip->unique_id().":$type) already exists in the database\nMaybe you want to recover?");
	}
      }else{

	$channel =  Bio::EnsEMBL::Funcgen::Channel->new
	(
	 -EXPERIMENTAL_CHIP_ID => $echip->dbID(),
	 -TYPE                 => $type,
	);
	
	($channel) = @{$chan_adaptor->store($channel)};
      }
    }
  }


  
  #Now get rset using experiment echips
  my $rset = $self->get_import_ResultSet('experimental_chip', $analysis);

  if($rset){#we have some new data

    foreach my $echip(@{$self->experiment->get_ExperimentalChips()}){
      
      if($echip->has_status('IMPORTED_SangerPCR', $echip)){
	$self->log("ExperimentalChip(".$echip->unique_id().") has already been imported");
      }else{
	$self->log("Reading SANGER result file:\t".$chip_files{$echip->unique_id()});

	my $fh = open_file("<", $chip_files{$echip->unique_id()});
	my $rfile = open_file(">", $self->get_dir("norm")."/result.SangerPCR.".$echip->unique_id().".txt");
	my $r_string = "";
	my $cc_id = $rset->get_chip_channel_id($echip->dbID());
	
	while($line = <$fh>){
	  $line =~ s/\r*\n//o;
	  
	  ($ratio, $pid) = (split/\t/, $line)[3..4];
	  $pid =~ s/.*://o;

	  $ratio = '\N' if $ratio eq 'NA';#NULL is still useful info to store in result

	  
	  #this is throwing away the encode region which could be used for the probeset/family?	
	  $r_string .= "\t".$self->get_probe_id_by_name($pid)."\t${ratio}\t${cc_id}\n";
	}
	
	print $rfile $r_string;
	close($rfile);
      }
    }
  }else{
    $self->log("No new data, skipping result parse");
  }

  $self->log("Finished parsing ".$self->vendor()." probe data (".localtime().")");
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
	$self->log("Rolling back partially imported ArrayChip:\t".$achip->design_id());
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
      $fh = open_file("<", $self->get_def("design_dir")."/".$achip->name().".pos");
      
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
      $fh = open_file("<", $self->get_def("design_dir")."/".$achip->name().".ndf");
      #Need to set these paths in each  achip hash, file names could be tablename.chip_id.txt
      #my $p_out = open_file(">", $self->get_dir("import")."/probe.".$ac{'design_name'}."txt");
      #my $ps_out = open_file(">", $self->get_dir("import")."/probe_set.".$ac{'design_name'}.".txt");
      #my $pf_out = open_file(">", $self->get_dir("import")."/probe_feature.".$ac{'design_name'}."txt");
      my $f_out = open_file(">", $self->get_dir("output")."/probe.".$achip->name()."fasta")	if($self->{'_dump_fasta'});
      my ($length, $ops, $op, $of, @probes, @features, %pfs);

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
	    #$self->store_set_probes_features($achip->dbID(), $ops, \@probes, \@features);
	    $self->store_set_probes_features($achip->dbID(), $ops, \%pfs);
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
	  $self->store_set_probes_features($achip->dbID(), $ops, \%pfs);
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
      $self->store_set_probes_features($achip->dbID(), $ops, \%pfs);
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






=head2 store_set_probes_features

  Arg [1]    : mandatory - array chip id
  Arg [2]    : optional - Bio::EnsEMBL::Funcgen::ProbeSet
  Arg [3]    : mandatory - hashref of keys probe id, values are 
               hash of probe/features with values 
               Bio::EnsEMBL::Funcgen::Probe/Features for a given 
               probe set if defined.
  Example    : $self->store_set_probes_features($ac->dbID(), $ops, \%pfs);
  Description: Stores probe set, probes and probe features 
  Returntype : none
  Exceptions : none
  Caller     : self
  Status     : Medium

=cut


sub store_set_probes_features{
  #my ($self, $ac_id, $ops, $probes, $features) = @_;

  my ($self, $ac_id, $ops, $pf_hash) = @_;
  

  #just call these directly rather than setting each time?
  #my $op_a = $self->db->get_ProbeAdaptor();
  #my $opf_a = $self->db->get_ProbeFeatureAdaptor();


  #if(scalar(@$probes) != scalar(@$features)){
  #  my $t = scalar(@$probes);
  #  my $f = scalar(@$features);
  #  throw("Does not accomodate multiple features($f) per probe $t");
  #}

  #Maybe do some validation here against all probes for probeset and ac_id? 

  if($ops){
    $ops->size(scalar(keys %$pf_hash));
    ($ops) = $self->db->get_ProbeSetAdaptor->store($ops);
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
    ($probe) = @{$self->db->get_ProbeAdaptor->store($probe)};

      

    #Can't use get_all_Arrays here as we can't guarantee this will only ever be the array we've generated
    #Might dynamically load array if non-present
    #This is allowing multiple dbIDs per probe???  Is this wrong?
    $self->cache_name_id($probe->get_probename(), $probe->dbID());

        
    foreach my $feature(@{$pf_hash->{$probe_id}->{'features'}}){
      $feature->probe($probe);
      ($feature) = @{$self->db->get_ProbeFeatureAdaptor->store($feature)};
    }
  }

  undef $ops;#Will this persist in the caller?
  undef %{$pf_hash};

# undef @$probes;
#  undef @$features,

  return;
}

sub cache_slice{
  my ($self, $region_name, $cs_name) = @_;

  throw("Need to define a region_name to cache a slice from") if ! $region_name;

  $cs_name ||= 'chromosome';
  $self->{'slice_cache'} ||= {};
  $region_name =~ s/chr//;
  $region_name = "MT" if $region_name eq "M";
  
  #can we handle UN/random chromosomes here?
  
  
  if(! exists $self->{'slice_cache'}->{$region_name}){
    $self->{'slice_cache'}->{$region_name} = $self->slice_adaptor->fetch_by_region($cs_name, $region_name);
    warn("Could not generate a slice for ${cs_name}:$region_name") if ! defined $self->{'slice_cache'}->{$region_name};
  }
  
  return $self->{'slice_cache'}->{$region_name};
}


=head2 cache_name_id

  Arg [1]    : mandatory - probe name
  Arg [2]    : mandatory - probe dbID
  Example    : $self->cache_name_id("Probe1", $probe->dbID());
  Description: Setter for probe cache values
  Returntype : none
  Exceptions : throws is cache conflict encountered
  Caller     : self
  Status     : At risk - merge with following?

=cut


sub cache_name_id{
  my ($self, $pname, $pid) = @_;


  throw("Must provide a probe name and id") if (! defined $pname || ! defined $pid);


  if(defined $self->{'_probe_cache'}->{$pname} && ($self->{'_probe_cache'}->{$pname} != $pid)){
    throw("Found two differing dbIDs for $pname, need to sort out redundant probe entries");
  }

  $self->{'_probe_cache'}->{$pname} = $pid;
  return;
}

=head2 get_probe_id_by_name

  Arg [1]    : mandatory - probe name
  Example    : $pid = $self->get_probe_id_by_name($pname);
  Description: Getter for probe cache values
  Returntype : int
  Exceptions : none
  Caller     : self
  Status     : At risk - merge with previous, move to importer?

=cut


sub get_probe_id_by_name{
  my ($self, $name) = @_;
  
  #Should only ever be one pid per probe per array per n array_chips
  #i.e. duplicate records per array chip with same pid

  if(! defined $self->{'_probe_cache'}){ #|| (! defined $self->{'_probe_cache'}->{$name})){

    #this fails if we're testing for the probe_id
    #this is because we have no chips associated with the experiemnt yet.
  
    $self->{'_probe_cache'} = $self->db->get_ProbeAdaptor->fetch_probe_cache_by_Experiment($self->experiment());
    
  
    #my $op = $self->db->get_ProbeAdaptor->fetch_by_array_probe_probeset_name($self->arrays->[0]->name(), $name);
    #$self->{'_probe_map'}{$name} = $op->dbID() if $op;
  }
  

  return (exists $self->{'_probe_cache'}->{$name}) ? $self->{'_probe_cache'}->{$name} : undef;
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
	      
	      $file = $self->get_def("results_dir")."/".$name.$suffix;

	      if(-f $file){
		$found = 1;
		last FILE;
	      }
	    }
	  }

	  throw("Could not find result file for Channel(${chan_name}) in ".$self->get_def('results_dir')) if ! $found; 


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
