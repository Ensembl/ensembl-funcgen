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
use Bio::EnsEMBL::Funcgen::Feature;
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
					#import_methods  => [],
					#data paths here?
													
					#but this is vendor specific and remains an array_def
					design_dir     => $self->get_dir("data")."/input/".$self->vendor()."/".$self->name()."/DesignFiles",
					chip_file        => $self->get_dir("data")."/input/".$self->vendor()."/".$self->name()."/SampleKey.txt",
					array_file       => $self->get_dir("data")."/input/".$self->vendor()."/".$self->name()."/DesignNotes.txt",
					results_file     => $self->get_dir("data")."/input/".$self->vendor()."/".$self->name()."/PairData/All_Pair.txt",
					
					results_dir     => $self->get_dir("data")."/input/".$self->vendor()."/".$self->name()."/PairData",
					norm_method => 'vsn_norm',
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
	
	#warn "Setting ".$self->vendor()." defs to ".$array_defs{$self->vendor()}." ".Data::Dumper::Dumper($array_defs{$self->vendor()});
	
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
  warn ("Reading chip data");

  my ($line, $chip_uid, $dye, $design_name, $design_id, $sample_label,  $sample_desc);
  my ($array, $species, $sample_type, $channel, $array_chip, %hpos, @data);
  my $tmp_uid = "FIRST";
  #Need some way of capturing other experimental variables?
  
  
  my ($echip, $ac_id);
  my $oa_adaptor = $self->db->get_ArrayAdaptor();
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

      #%tmp = (
      #       array_id    => undef,
      #	      dbID        => undef,
      #	      name => $data[$hpos{'DESIGN_NAME'}],
      #	     );
         
      #This is treating each array chip as a separate array, unless arrayset is defined
      #AT present we have no way of differentiating between different array_chips on same array???!!!
      #Need to add functionality afterwards to collate array_chips into single array
      
      
      #store now checks whether already stored and updates array chips accordingly
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
							  -NAME      => $data[$hpos{'DESIGN_NAME'}],
							  -DESIGN_ID => $data[$hpos{'DESIGN_ID'}],
							 );


      $array->addArrayChip($array_chip);
      #($array) = @{$oa_adaptor->store($array)};

        
      #Need to reg array chip here (to differentiate from what is already in DB, handles subset of arrays)
      #This is a registry of array chips which have previously been stored for validation purpose
      #To avoid importing probes twice

    }elsif((! $array->get_ArrayChip_by_design_id($data[$hpos{'DESIGN_ID'}])) && ($self->{'array_set'})){

      warn "generating new ac for same array ".$data[$hpos{'DESIGN_ID'}]."\n";

      $array_chip = Bio::EnsEMBL::Funcgen::ArrayChip->new(
							  -NAME        => $data[$hpos{'DESIGN_NAME'}],
							  -DESIGN_ID => $data[$hpos{'DESIGN_ID'}],
							 );

   
      $array->add_ArrayChip($array_chip);
   
      
    }elsif(! $array->get_ArrayAhip_by_design_id($data[$hpos{'DESIGN_ID'}])){
      throw("Found experiment with more than one design");
    }
    
    #Parse and Populate ExperimentalChip/Channels
    if ($tmp_uid eq "FIRST" || $data[$hpos{'CHIP_ID'}] != $tmp_uid){
      $tmp_uid = $data[$hpos{'CHIP_ID'}];
     
      $echip =  Bio::EnsEMBL::Funcgen::ExperimentalChip->new
	(
	 #-DESIGN_ID      => $data[$hpos{'DESIGN_ID'}],
	 -EXPERIMENT_ID  => $self->experiment->dbID(),
	 -DESCRIPTION    => $data[$hpos{$sample_desc}],
	 -ARRAY_CHIP_ID  => $array->get_ArrayChip_by_design_id($data[$hpos{'DESIGN_ID'}])->dbID(),
	 -UNIQUE_ID      => $data[$hpos{'CHIP_ID'}],
	);
      
      $echip = $self->experiment->add_experimental_chip($echip);

      #should we nest these in the Experiment and 
      #don't need to add them to store, just have method which always retrieves all echips from db
      #$self->{'echips'}{$data[$hpos{'CHIP_ID'}]} = $echip;#do we still need this?
      
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


=head2 arrays

  Arg [1]    : optional - Bio::EnsEMBL::Funcgen::Array
  Example    : $self->arrays($array);
  Description: Getter/Setter for array element
  Returntype : list ref to Bio::EnsEMBL::Funcgen::Array objects
  Exceptions : throws if passed non Array or if more than one Array set
  Caller     : Importer
  Status     : Medium - Remove/Implement multiple arrays?

=cut

sub arrays{
  my ($self) = shift;

  if(@_ && ! $_[0]->isa('Bio::EnsEMBL::Funcgen::Array')){
    throw("Must supply a Bio::EnsEMBL::Funcgen::Array");
  }elsif(@_){
    warn "this is adding an array, not resetting";
    push @{$self->{'arrays'}}, @_;
  }
   
  throw("Does not yet support multiple array imports") if(scalar (@{$self->{'arrays'}}) > 1);
  #need to alter read_probe data at the very least
  
  return $self->{'arrays'};
}


=head2 echip_data

  Arg [1]    : mandatory - ExperimentalChip design_id
  Arg [2]    : mandatory - data type/name
  Arg [3]    : optional - value
  Example    : $self->echip_data($design_id, $data_type, $value));
  Description: Getter/Setter for ExperimentalChip data from cache
  Returntype : various
  Exceptions : throws if no design id or data type defined
  Caller     : Importer
  Status     : Deprecated

=cut


#sub echip_data{
 # my ($self, $design_id, $data_type, $value) = @_;

#  deprecate("Use ExperimentalChip methods directly");

#  throw("Need to specify a design_id  and a data_type") if (! defined $data_type || ! defined $design_id);
  
#  if(defined $value){
#    ${$self->get_data('echips', $design_id)}{$data_type} = $value;#can we deref with -> instead to set value?
#}
#else{
#  return ${$self->get_data('echips', $design_id)}{$data_type};#will this cause undefs?
#}
#}


=head2 get_channels

  Arg [1]    : mandatory - ExperimentalChip unique_id
  Example    : $chip = $self->get_echip($chip_uid);
  Description: Getter for Channel hash from cache
  Returntype : Hashref
  Exceptions : throws if no chip unique id defined
  Caller     : Importer
  Status     : At risk - use ExperimentalChip to cache

=cut


sub get_channels{
  my ($self, $chip_uid) = @_;
  throw("Need to specify a chip unique_id") if (! defined $chip_uid);

  #rename and move to Experiment?

  return $self->experiment->get_experimental_chip_by_unique_id($chip_uid)->{'channels'};
}

=head2 get_channel

  Arg [1]    : mandatory - channel "unique_id", ${chip_uid}_${dye_freq}
  Example    : $channel = $self->get_echip($chan_uid);
  Description: Getter for Channel hash from cache
  Returntype : Hashref
  Exceptions : throws if no channel unique id defined
  Caller     : Importer
  Status     : At risk - use ExperimentalChip to cache

=cut


sub get_channel{
  my ($self, $chan_uid) = @_;

  #rename and move to Experiment?

  throw('Need to specify a chan unique_id (${chip_uid}_${dye_freq})') if (! defined $chan_uid);	
  return $self->{'channels'}{$chan_uid};
}


=head2 set_channel

  Arg [1]    : mandatory - ExperimentalChip unique_id
  Arg [2]    : mandatory - Bio::EnsEMBL::Funcgen::Channel 
  Example    : $self->set_channel($chip_uid, $channel);
  Description: Setter for Channel cache
  Returntype : none
  Exceptions : throws if no chip unique id or channel defined
  Caller     : Importer
  Status     : At risk - Move to ExperimentalChip, hash on dye freq

=cut


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

=head2 channel_data

  Arg [1]    : mandatory - channel "unique_id", ${chip_uid}_${dye_freq}
  Arg [2]    : mandatory - data type/name
  Arg [3]    : optional - value
  Example    : $self->set_channel($chip_uid, $channel);
  Description: Setter for Channel cache
  Returntype : none
  Exceptions : throws if no chan unique id or data type defined
  Caller     : Importer
  Status     : Deprecated

=cut


sub channel_data{
  my ($self, $chan_uid, $data_type, $value) = @_;
	
  deprecate("Use Channel methods directly via echip cache");
  
  throw("Need to provide a channel uid and a data_type") if (! defined $chan_uid || ! defined $data_type);

  if(defined $value){
    throw("should now use Channel methods \n");
    ${$self->get_channel($chan_uid)}{$data_type} = $value;
  }
  else{
    return ${$self->get_channel($chan_uid)}{$data_type};
  }
}


sub read_sanger_array_probe_data{
  my ($self, $array_file) = @_;

  $array_file||= $self->array_file();
  my ($line, $fh, @list, $array_file_format, $cmd);
  my ($op, $of, $imported, $fimported, %slices);
  my $oa_adaptor = $self->db->get_ArrayAdaptor();
  my $op_adaptor = $self->db->get_ProbeAdaptor();
  my $of_adaptor = $self->db->get_ProbeFeatureAdaptor();
  my $ec_adaptor = $self->db->get_ExperimentalChipAdaptor();
  my $slice_adaptor = $self->db->get_SliceAdaptor();
  my $anal_id = $self->db->get_AnalysisAdaptor->fetch_by_logic_name("SangerPCR")->dbID();
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
							);
  $array->add_ArrayChip($array_chip);
 
  $self->arrays($array);
  my $ac_id = $array->get_ArrayChip_by_design_id($array->name())->dbID();


  #we also need to test wether the array as been imported as well as the mappings
  #THis needs to use coord_sys-id not schema_build!!  Duplcaite entries for different schema_builds 
  #with same assembly

  my $dnadb_cs = $self->db->dnadb->get_CoordSystemAdaptor->fetch_by_name('chromosome');
  my $fg_cs = $self->db->get_FGCoordSystemAdaptor->validate_coord_system($dnadb_cs);


  #This fails if we're pointing to an old DB during the release cycle.  Will be fine if we manage to cs mapping dynamically

  #if($self->db->fetch_status_by_name('array_chip', $ac_id, 'IMPORTED_CS_'.$fg_cs->dbID())){
  if($self->status_adaptor->has_status('IMPORTED_CS_'.$fg_cs->dbID(), $array_chip)){

    $fimported = 1;
    $imported = 1;
    warn("Skipping array chip feature import (".$array->name().") already fully imported for ".$self->data_version()."\n");
  }
  elsif($self->status_adaptor->has_status('IMPORTED', $array_chip)){
    $imported = 1;
    warn("Skipping array chip probe import (".$array->name().") already fully imported\n");
  }

  #need to check whether already imported on specified schema_build
  #check for appropriate file given format in input dir or take path

  if(! $fimported){

    if(! $array_file){

      throw("No input_dir defined, if you are running in a non Experiment context please use -array_file") if(! defined $self->get_dir('input'));
      
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
    }else{
      throw("Could not determine array file format: $array_file");
    }

    my $fanal = $self->db->get_AnalysisAdaptor->fetch_by_logic_name(($array_file_format eq "ADF") ? "VendorMap" : "LiftOver");
   
    $self->log("Parsing ".$self->vendor()." array data (".localtime().")");
    warn("Parsing ".$self->vendor()." array data (".localtime().")");
    
    $fh = open_file("<", $array_file);
    

    warn("Hardcoded for GFF file format");

    my ($chr, $start, $end, $strand, $pid);

    while($line = <$fh>){
      $line =~ s/\r*\n//;
      
      #($chr, $start, $end, $ratio, $pid) = split/\t/o, $line;
      ($chr, undef, undef, $start, $end, undef, $strand, undef, $pid) = split/\t|\;/o, $line;
      $pid =~ s/reporter_id=//o;

      #warn ("($chr, undef, undef, $start, $end, undef, $strand, undef, $pid)");
 

      #need to parse dependant on file format 
      #also need to account for duplicate probes on grid
      
      if(! $imported){
	#when we utilise array coords, we need to look up probe cache and store again with new coords
	#we're currently storing duplicates i.e. different ids with for same probe
	#when we should be storing two records for the same probe/id
	#the criteria for this will be different for each vendor, may have to check container etc for NimbleGen

	if(! $self->get_probe_id_by_name($pid)){
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

       	  #build slice hash cache
	  $chr  =~ s/chr//;
	  $strand = ($strand eq "+") ? 0 : 1;
	  if(! exists $slices{$chr}){
	    $slices{$chr} = $slice_adaptor->fetch_by_region('chromosome', $chr);
	  }
	  
	  #Hack!!!!!!
	  if(!  $slices{$chr}){
	    warn("Skipping non standard probe (".$pid.") with location:\t${chr}:${start}-${end}\n");
	    #pop @probes;
	    next;
	  }
	  
	
	  $of = Bio::EnsEMBL::Funcgen::ProbeFeature->new(
							 -START         => $start,
							 -END           => $end,
							 -STRAND        => $strand,
							 -SLICE         => $slices{$chr},
							 -ANALYSIS      => $fanal,
							 -MISMATCHCOUNT => 0,
							 -_PROBE_ID     => $self->get_probe_id_by_name($pid),#work around to avoid cacheing probes
							);

	  $of_adaptor->store($of);

	}else{
	  #warn "Need to accomdate duplicate probes here";¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡!!!!!!!
	}
      }
    }
    warn"setting ac cs status";

    #$self->db->set_status('array_chip', $ac_id, 'IMPORTED_CS_'.$fg_cs->dbID());
    $self->status_adaptor->set_status('IMPORTED_CS_'.$fg_cs->dbID(), $array_chip);
  }

  

  if(! $imported){
    warn"setting ac imp  status";
    #$self->db->set_status('array_chip', $ac_id, 'IMPORTED');
    $self->status_adaptor->set_status('IMPORTED', $array_chip);

    $imported = 1;
  }
  
  $self->log("Finished parsing ".$self->vendor()." array/probe data (".localtime().")");
  warn("Finished parsing ".$self->vendor()." array/probe data (".localtime().")");
  
  return;
}

=head2 read_sanger_probe_data

  Example    : $imp->read_sanger_probe_data();
  Description: Parses and imports probes, features and result for the sanger PCR array platform
  Returntype : none
  Exceptions : none
  Caller     : Importer
  Status     : At risk - Move parts to "Vendor"Defs.pm, should function the same

=cut

sub read_sanger_result_data{
  my $self = shift;


  warn("This should never import array/probes/feature.  Force adf/gff import and use this simply to import results");
  #shoudl also check wether experimental_chips have been previously imported
  

  my ($file, $result_set, $chip_uid, $line, $echip, $fh, $r_string, $rfile);
  my ($ratio, $pid, $imported, %tmp);
  my $of_adaptor = $self->db->get_ProbeFeatureAdaptor();
  my $ec_adaptor = $self->db->get_ExperimentalChipAdaptor();
  my $anal_id = $self->db->get_AnalysisAdaptor->fetch_by_logic_name("SangerPCR")->dbID();
  my $result_adaptor = $self->db->get_ResultSetAdaptor();
    			       

         
  #This is treating each array chip as a separate array, unless arrayset is defined
  #AT present we have no way of differentiating between different array_chips on same array???!!!
  #Need to add functionality afterwards to collate array_chips into single array
  
 
  #this is done to avoid having to self->array_name in loop, will make multiple array loop easier 
  my $array = ${$self->arrays()}[0];


  $self->log("Parsing ".$self->vendor()." probe data (".localtime().")");
  warn("Parsing ".$self->vendor()." probe data (".localtime().")");
  
  my $list = "ls ".$self->input_dir()."/*all*";
  
  foreach $file(`$list`){

    warn($file);

    ($chip_uid = $file) =~ s/.*\///;
    $chip_uid =~ s/\.all.*\n//;
    
    $echip =  Bio::EnsEMBL::Funcgen::ExperimentalChip->new
      (
       -EXPERIMENT_ID  => $self->experiment->dbID(),
       -DESCRIPTION    => "",
       -ARRAY_CHIP_ID  => $self->arrays->[0]->get_ArrayChip_by_design_id($array->name())->dbID(),
       -UNIQUE_ID      => $chip_uid,
      );
      
    $echip = $self->experiment->add_experimental_chip($echip);

    if($self->status_adaptor->has_status('IMPORTED', $echip)){
      $imported = 1;
      warn("Skipping experimental chip import (".$echip->unique_id().") already fully imported\n");
    }

    #Need to check fo rpartial import here


    
    



    #should we nest these in the Experiment and 
    #don't need to add them to store, just have method which always retrieves all echips from db
   
    
    if(! $imported){

      $fh = open_file("<", $file);
      $rfile = open_file(">", $self->get_dir("norm")."/result.".$echip->unique_id().".txt");
      $r_string = "";
      
      #as we don't know contig chips before hand set everything to the same result_set, then alter after wards
      if(! defined $result_set){

	$result_set = Bio::EnsEMBL::Funcgen::ResultSet->new
	  (
	   -analysis_id => $anal_id,
	   -table_name  => 'experimental_chip',
	   -table_id    => $echip->dbID(),
	  );
      }else{
	$result_set->add_table_id($echip->dbID());
      }
   
      while($line = <$fh>){
	$line =~ s/\r*\n//o;
	
	($ratio, $pid) = (split/\t/, $line)[3..4];


	$pid =~ s/.*://o;
	#warn"($ratio, $pid)";
	#this is throwing away the encode region which could be used for the probeset/family?
	
	#NA ratio imports as 0
	#$r_string .= "\t".$self->get_probe_id_by_name($pid)."\t${ratio}\t${anal_id}\t".$echip->dbID()."\texperimental_chip\n";
	$r_string .= "\t".$self->get_probe_id_by_name($pid)."\t${ratio}\n";
	
      
      }

      print $rfile $r_string;
      close($rfile);
    
      #should import drectly here?
    
      #need to do this for experimental_chip after results imported with logic name
      #if(! $imported){
      #  $self->db->set_status('_chip', $ac_id, "IMPORTED");
      #  $imported = 1;
      #}
    }
  }

  $result_adaptor->store($result_set);


  $self->log("Finished parsing ".$self->vendor()." probe data (".localtime().")");
  warn("Finished parsing ".$self->vendor()." probe data (".localtime().")");
  
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

  my ($fh, $line, @data, %hpos);
  $self->log("Parsing probe data (".localtime().")");
  warn("Parsing probe data (".localtime().")...can we do a reduced parse if we know the array chips are already imported");
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
      my $achip = $array->get_ArrayChip_by_design_id($design_id);

   
      #check status of array chip here
      #if($self->db->fetch_status_by_name('array_chip', $achip->dbID(), 'IMPORTED')){
      if($self->status_adaptor->has_status('IMPORTED', $achip)){
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
  

  
   $fh = open_file("<", $self->get_def("design_dir")."/".$achip->name().".ndf");
      #Need to set these paths in each  achip hash, file names could be tablename.chip_id.txt
      #my $p_out = open_file(">", $self->get_dir("import")."/probe.".$ac{'design_name'}."txt");
      #my $ps_out = open_file(">", $self->get_dir("import")."/probe_set.".$ac{'design_name'}.".txt");
      #my $pf_out = open_file(">", $self->get_dir("import")."/probe_feature.".$ac{'design_name'}."txt");
      my $f_out = open_file(">", $self->get_dir("output")."/probe.".$achip->name()."fasta")	if($self->{'_dump_fasta'});
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
	  $self->store_set_probes_features($achip->dbID(), $ops, \%pfs);

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
	
	
	
	$of = Bio::EnsEMBL::Funcgen::ProbeFeature->new(
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
      $self->store_set_probes_features($achip->dbID(), $ops, \%pfs);
    $self->db->set_status('array_chip', $achip->dbID(), "IMPORTED");
    
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

  if(defined $self->{'_probe_map'}->{$pname} && ($self->{'_probe_map'}->{$pname} != $pid)){
    throw("Found two differing dbIDs for $pname, need to sort out redundant oligo_probe entries");
  }

  $self->{'_probe_map'}->{$pname} = $pid;
  return;
}

=head2 get_probe_id_by_name

  Arg [1]    : mandatory - probe name
  Example    : $pid = $self->get_probe_id_by_name($pname);
  Description: Getter for probe cache values
  Returntype : int
  Exceptions : none
  Caller     : self
  Status     : At risk - merge with previous

=cut


#Remove array element to this?
sub get_probe_id_by_name{
  my ($self, $name) = @_;
  
  #Should only ever be one pid per probe per array per n array_chips
  #i.e. duplicate records per array chip with same pid

#  if((defined $self->{'_probe_map'}) && (defined $self->{'_probe_map'}->{$name})){
#    #@op_ids = @{$self->{'_probe_map'}->{$name}};
#  }
#  else{#get from db

  if((! defined $self->{'_probe_map'}) || (! defined $self->{'_probe_map'}->{$name})){ 
    my $op = $self->db->get_ProbeAdaptor->fetch_by_array_probe_probeset_name($self->arrays->[0]->name(), $name);
    #print "Got probe $op with dbid ".$op->dbID()."\n";
    #push @op_ids, $op->dbID();
    $self->{'_probe_map'}{$name} = $op->dbID() if $op;
  }
  
  #return \@op_ids;

 # warn "testing for pid $name";

  return (exists $self->{'_probe_map'}->{$name}) ? $self->{'_probe_map'}->{$name} : undef;
  
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
  
  $self->log("Parsing results(".localtime().")...");
  warn("Parsing results(".localtime().")...");
   
  my ($i, $fh, $tmp, $line, $r_string, $probe_elem, $first_result, $file_name, @header, @data, @design_ids);
  my $anal = $self->db->get_AnalysisAdaptor->fetch_by_logic_name("RawValue");
  #should check here if defined
  
  my $anal_id = $anal->dbID();
  

  #as we're importing all echip result in parallel, there is no point in checking for IMPORTED?
  #Or do we need to set flag for experiment or each chip and check somehow?
  my $cnt = 0;

  foreach my $array(@{$self->arrays()}){

    @design_ids = @{$array->get_design_ids()};
   
    foreach my $design_id(@design_ids){
     
      my $achip = $array->get_ArrayChip_by_design_id($design_id);
      $r_string = "";
      warn("Reading/Importing results for $design_id\n");
      my $r_out = open_file(">", $self->get_dir("raw")."/result.".$achip->name().".txt");


      $file_name = (scalar(@design_ids) > 1) ? $achip->name() : "All";  
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
	  #foreach my $pid(@{$self->get_probe_ids_by_name($data[$probe_elem])}){
	  #import directly here?
	  #print $r_out "\t${pid}\t".$data[$i]."\t$anal_id\t".$self->get_channel($tmp)->dbID()."\tchannel\n";	  
	  #$self->db->insert_table_row("result", $pid, $data[$i], $anal_id, $self->get_channel($tmp)->dbID(), "channel");
	  #Build multiple instert sql statement rather than executing once for each result?
	  

	  $cnt ++;
	  
	  #Need to change this get_channel call?
	  ($tmp = $header[$i]) =~ s/1h_//;

	  $r_string .= "\t".$self->get_probe_id_by_name($data[$probe_elem])."\t".$data[$i]."\t$anal_id\t".$self->get_channel($tmp)->dbID()."\tchannel\n";
	  #}
	}


	
	if($cnt > 10000){
	  $cnt = 0;
	  print $r_out $r_string;
	  $r_string ="";
	  #could we fork here and import in the background?

	}


      }
      
      #my $r_out = open_file(">", $self->get_dir("import")."/result.".$ac{'name'}.".txt");
      print $r_out $r_string;
      close($r_out);
    }
  }

  $self->log("Finished parsing and results");
  warn("Finished parsing and results");

  return;
}



1;
