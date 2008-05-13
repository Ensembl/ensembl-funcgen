#
# EnsEMBL module for Bio::EnsEMBL::Funcgen::Parsers::GFF
#

#Could this be based on a Generic Flat file parser?

=head1 NAME

Bio::EnsEMBL::Funcgen::Parsers::GFF

=head1 SYNOPSIS

  my $parser_type = "Bio::EnsEMBL::Funcgen::Parsers::GFF";
  push @INC, $parser_type;
  my $imp = $class->SUPER::new(@_);


=head1 DESCRIPTION

This is a definitions class which should not be instatiated directly, it 
normally set by the Importer as the parent class.  Bed contains meta 
data and methods specific to data in bed format, to aid 
parsing and importing of experimental data.

=head1 AUTHOR

This module was created by Stefan Graf.

=head1 CONTACT

Post questions to the EnsEMBL development list ensembl-dev@ebi.ac.uk

=head1 METHODS

=cut

package Bio::EnsEMBL::Funcgen::Parsers::Bed;

use Bio::EnsEMBL::Utils::Exception qw( throw warning deprecate );
use Bio::EnsEMBL::Utils::Argument qw( rearrange );
use strict;


use vars qw(@ISA);
@ISA = qw(Bio::EnsEMBL::Funcgen::Parsers::Simple);

=head2 new

  Example    : my $self = $class->SUPER::new(@_);
  Description: Constructor method for GFF class
  Returntype : Bio::EnsEMBL::Funcgen::Parsers::GFF
  Exceptions : None
  Caller     : Bio::EnsEMBL::Funcgen::Importer
  Status     : at risk

=cut


sub new{
  my $caller = shift;
  
  my $class = ref($caller) || $caller;

  #define default fields here and pass
  #We also need to be able to take custom attrs mappings

  #keys are array index of field, key are Feature paramter names
  #reverse this?
  #Unless we have a compound field which we name accordingly
  #And then call e.g. parse_attrs
  #Which will return a hash with the relevant Feature attributes

  #Is splitting this up simply going to make the parse slower due to acessor methods?

  #Pass or just set directly here?
  #<seqname> <source> <feature> <start> <end> <score> <strand> <frame> [attributes] [comments]


  #Some of these may be highly redundant due to the nature of the data.
  #We can hash things to lessen the impact but we're still going to be checking if exists for each one
  #No way around this?  Unless it is marginally faster to set a permanent type and then only check a boolean.
  #Yes there is, this is the exhaustive GFF definition, we can just redefine or delete some entries dynamically to
  #avoid ever considering a particular field index.

  my %fields = (
				0 => 'fetch_slice',
				1 => 'get_source',
				2 => 'get_feature_type',
				3 => '-start',
				4 => '-end',
				5 => '-strand',#Will most likely be , need to convert to -.+ > -1 0 1
				#6 => 'frame',#will most likely be .
				7 => 'get_attributes',
			   );

  my $self  = $class->SUPER::new(@_);
  

  #We need to define meta header method, starting with '##'
  #Also need to skip comments '#' at begining or end of line
  #Do we also need to skip field header? No methinks not.

  #Define result method
 # $self->{'file_ext'} => 'gff';#Could use vendor here?
  
  #define this if we want to override the generic method in Simple
  #$self->{'config'}{'results_data'} => ["and_import_gff"];  


  return $self;
}


=head2 set_config

  Example    : my $self->set_config;
  Description: Sets attribute dependent config
  Returntype : None
  Exceptions : None
  Caller     : Bio::EnsEMBL::Funcgen::Importer
  Status     : at risk

=cut


sub set_config{
  my $self = shift;

  $self->SUPER::set_config;

  #GFF specific stuff here.

  return;
}



 
sub read_and_import_bed_data{
    my $self = shift;
    
    $self->log("Reading and importing ".$self->vendor()." data");
    my (@header, @data, @design_ids, @lines);
    my ($anal, $fh, $file);
    
    my $eset_adaptor = $self->db->get_ExperimentalSetAdaptor();
    my $af_adaptor = $self->db->get_AnnotatedFeatureAdaptor();
    my $fset_adaptor = $self->db->get_FeatureSetAdaptor();
    my $dset_adaptor = $self->db->get_DataSetAdaptor();
	my $new_data = 0;
	my $eset = $eset_adaptor->fetch_by_name($self->experimental_set_name());
    
    if(! defined $eset){
        $eset = Bio::EnsEMBL::Funcgen::ExperimentalSet->new(
                                                            -name         => $self->experimental_set_name(),
                                                            -experiment   => $self->experiment(),
                                                            -feature_type => $self->feature_type(),
                                                            -cell_type    => $self->cell_type(),
                                                            -vendor       => $self->vendor(),
                                                            -format       => $self->format(),
															-analysis     => $self->feature_analysis,
                                                            );
        ($eset)  = @{$eset_adaptor->store($eset)};
    }

    #we need a way to define replicates on a file basis when we have no meta file!
    #can we make this generic for application to array imports?
    #currently we have to do a separate import for each replicate, specifying the result files each time
    #we need to add a experimental_set_name option
    #Actually ExperimentalSet is a little redundant as we can do the roll back which is exactly what this is designed to facilitate
    #It does however allow incremental addition of new subsets

    #Now define FeatureSet
    
    #shouldn't we be using the exp name?
    my $fset = $fset_adaptor->fetch_by_name($self->experiment->name());

    if(! defined $fset){

        #currently hardcoded, but we should probably add feature_analysis_name
	  $fset = Bio::EnsEMBL::Funcgen::FeatureSet->new(
													 -name         => $self->experiment->name(),
													 -feature_type => $self->feature_type(),
													 -cell_type    => $self->cell_type(),
													 -type         => 'annotated',
													 -analysis     => $self->feature_analysis,
													);
        ($fset)  = @{$fset_adaptor->store($fset)};
    }

    #Define DataSet
    my $dset = $dset_adaptor->fetch_by_name($self->experiment->name());
    
    if(! defined $dset){

        $dset = Bio::EnsEMBL::Funcgen::DataSet->new(
                                                    -name                => $self->experiment->name(),
                                                    -supporting_sets     => [$eset],
                                                    -feature_set         => $fset,
                                                    -displayable         => 1,
                                                    -supporting_set_type => 'experimental',
                                                    );
        ($dset)  = @{$dset_adaptor->store($dset)};
    }
    

    #Get file
    if (! @{$self->result_files()}) {
        my $list = "ls ".$self->input_dir().'/'.$self->name().'*.bed';
        my @rfiles = `$list`;
        $self->result_files(\@rfiles);
    }
    
    if (scalar(@{$self->result_files()}) >1) {
        warn("Found more than one bed file:\n".
             join("\n", @{$self->result_files()})."\nBed does not yet handle replicates.".
             "  We need to resolve how we are going handle replicates with random cluster IDs");
        #do we even need to?
    }

	#Here were are tracking the import of individual bed files by adding them as ExperimentalSubSets
	#Recovery would never know what to delete
	#So would need to delete all, Hence no point in setting status?



	### VALIDATE FILES ###
	#We need validate all the files first, so the import doesn't fall over half way through
	#Or if we come across a rollback halfway through
	my %new_data;
	my $roll_back = 0;

	foreach my $filepath(@{$self->result_files()}) {
	  chomp $filepath;
	  my $filename;
	  ($filename = $filepath) =~ s/.*\///;
	  my $sub_set;
	  $self->log("Found bed file\t$filename");
	  
	  if($sub_set = $eset->get_subset_by_name($filename)){
		
		if ($sub_set->has_status('IMPORTED')){
		  $self->log("ExperimentalSubset(${filename}) has already been imported");
		} 
		else {
		  $self->log("Found partially imported ExperimentalSubset(${filename})");
		  $roll_back = 1;

		  if ($self->recovery() && $roll_back) {
			$self->log("Rolling back results for ExperimentalSubset:\t".$filename);
			warn "Cannot yet rollback for just an ExperimentalSubset, rolling back entire set\n";
			
			$self->rollback_FeatureSet($fset);
			$self->rollback_ExperimentalSet($eset);
			last;
		  }
		  elsif($roll_back){
			throw("Found partially imported ExperimentalSubSet:\t$filepath\n".
				  "You must specify -recover  to perform a full roll back for this ExperimentalSet:\t".$eset->name);
		  }
		}
	  }
	  else{
		$self->log("Found new ExperimentalSubset(${filename})");
		$new_data{$filepath} = undef;
		$sub_set = $eset->add_new_subset($filename);
		$eset_adaptor->store_ExperimentalSubsets([$sub_set]);
	  }
	}


	### READ AND IMPORT FILES ###
	foreach my $filepath(@{$self->result_files()}) {
	  chomp $filepath;
	  my $filename;
	  my $count = 0;
	  ($filename = $filepath) =~ s/.*\///;

	  if($roll_back || exists $new_data{$filepath}){

		 $self->log("Reading bed file:\t".$filename);
		 my $fh = open_file($filepath);
		 my @lines = <$fh>;
		 close($fh);
		 
		 my ($line, $f_out);
		 my $fasta = '';
		 
		 #warn "we need to either dump the pid rather than the dbID or dump the fasta in the DB dir";
		 my $fasta_file = $ENV{'EFG_DATA'}."/fastas/".$self->experiment->name().'.'.$filename.'.fasta';

		 if($self->dump_fasta()){
		   $self->backup_file($fasta_file);
		   $f_out = open_file($fasta_file, '>');
		 }
		
		 foreach my $line (@lines) {
		   $line =~ s/\r*\n//o;
		   next if $line =~ /^\#/;	
		   next if $line =~ /^$/;
		   next unless $line =~ /^chr/i;
		   #next if $line =~ /^chr/i;#Mikkelson hack
		   
		   my ($chr, $start, $end, $pid, $score) = split/\t/o, $line;				  
		   #my ($chr, $start, $end, $score) = split/\t/o, $line;#Mikkelson hack				  
		   

		   if($self->ucsc_coords){
			 $start +=1;
			 $end +=1;
		   }
		   
		   if(!  $self->cache_slice($chr)){
			 warn "Skipping AnnotatedFeature import, cound non standard chromosome: $chr";
		   }
		   else{
			 
			 #this is throwing away the encode region which could be used for the probeset/family?	
			 my $feature = Bio::EnsEMBL::Funcgen::AnnotatedFeature->new
			   (
				-START         => $start,
				-END           => $end,
				-STRAND        => 1,
				-SLICE         => $self->cache_slice($chr),
				-ANALYSIS      => $fset->analysis,
				-SCORE         => $score,
				-DISPLAY_LABEL => $pid,
				-FEATURE_SET   => $fset,
			   );
			 
			 $af_adaptor->store($feature);
			 
			 $count++;

			 #dump fasta here
			 if ($self->dump_fasta()){
			   $fasta .= '>'.$pid."\n".$self->cache_slice($chr)->sub_Slice($start, $end, 1)->seq()."\n";
			 }
		   }
		 }
		 
		 
		 if ($self->dump_fasta()){
		   print $f_out $fasta;
		   close($f_out);
		 }


		 $self->log("Finished importing $count features:\t$filepath");
		 my $sub_set = $eset->get_subset_by_name($filename);
		 $sub_set->adaptor->store_status('IMPORTED', $sub_set);
	   }
    }

    $self->log("No new data, skipping result parse") if ! keys %new_data && ! $roll_back;   
    $self->log("Finished parsing and importing results");


	#Retrieve DataSet here as sanity check it is valid
	#We had one which had a mismatched feature_type between the feature_set and the experimental_set
	#How was this possible, surely this should have failed on generation/storage?

    
    return;
}
  


1;
