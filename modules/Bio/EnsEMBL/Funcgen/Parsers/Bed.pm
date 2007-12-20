#
# EnsEMBL module for Bio::EnsEMBL::Funcgen::Parsers::Bed
#

=head1 NAME

Bio::EnsEMBL::Funcgen::Parsers::Bed

=head1 SYNOPSIS

  my $parser_type = "Bio::EnsEMBL::Funcgen::Parsers::Bed";
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

use Bio::EnsEMBL::Funcgen::ExperimentalSet;
use Bio::EnsEMBL::Funcgen::FeatureSet;
use Bio::EnsEMBL::Funcgen::AnnotatedFeature;
use Bio::EnsEMBL::Utils::Exception qw( throw warning deprecate );
use Bio::EnsEMBL::Funcgen::Utils::EFGUtils qw(species_chr_num open_file);
use Bio::EnsEMBL::Utils::Argument qw( rearrange );
use Bio::EnsEMBL::Funcgen::Utils::Helper;
use strict;

use Data::Dumper;

use vars qw(@ISA);
@ISA = qw(Bio::EnsEMBL::Funcgen::Utils::Helper);

=head2 new

  Example    : my $self = $class->SUPER::new(@_);
  Description: Constructor method for Bed class
  Returntype : Bio::EnsEMBL::Funcgen::Parsers::Bed
  Exceptions : throws if Experiment name not defined or if caller is not Importer
  Caller     : Bio::EnsEMBL::Funcgen::Importer
  Status     : at risk

=cut


sub new{
    my $caller = shift;

  my $class = ref($caller) || $caller;
  my $self  = $class->SUPER::new();

  throw("This is a skeleton class for Bio::EnsEMBL::Importer, should not be used directly") 
      if(! $self->isa("Bio::EnsEMBL::Funcgen::Importer"));
	
  $self->{'config'} =  
    {(
      #order of these data arrays is important!
      array_data   => [],#['experiment'],
      probe_data   => [],#["probe"],
      results_data => ["and_import_bed"],
      norm_method => undef,

	  #Need to make these definable?
	  #have protocolfile arg and just parse tab2mage protocol section format
	  #SEE NIMBLEGEN FOR EXAMPLE
	  protocols => {()},
     )};
	
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

  throw('Must provide an ExperimentalSet name for a Bed import') if ! defined $self->experimental_set_name();
  #dir are not set in config to enable generic get_dir method access

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
   

	#Don't need this isa, tested in feature_analysis method already
    #if ($self->feature_analysis->isa("Bio::EnsEMBL::Analysis")) {
    #    $anal = $self->feature_analysis;
    #}

	#Commented this out as ExperimentalSets should not be normalised array data
	#rather processed data from any format i.e. a feature analysis
	#elsif (defined $self->norm_method){
	#        $anal = $self->db->get_AnalysisAdaptor->fetch_by_logic_name($self->norm_method);
	#    } 

	#else {
	  #We could set a default 'Experimental' placeholder analysis here and warn
	  #Or just remove this whole block as ExperimentalSet->new will barf if the feature_analysis is not valid
    #    throw("No analysis set.");
    #}

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

    warn "got fset $fset";
    
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

        warn "got fset $fset";
    }

    #Now define DataSet

    #shouldn't we be using the exp name?
    my $dset = $dset_adaptor->fetch_by_name($self->experiment->name());

    warn "got dset $dset";
    
    if(! defined $dset){

        $dset = Bio::EnsEMBL::Funcgen::DataSet->new(
                                                    -name                => $self->experiment->name(),
                                                    -supporting_sets     => [$eset],
                                                    -feature_set         => $fset,
                                                    -displayable         => 1,
                                                    -supporting_set_type => 'experimental',
                                                    );
        ($dset)  = @{$dset_adaptor->store($dset)};

        warn "got dset $dset";
    }
    

    #Get file
    if (! @{$self->result_files()}) {
        my $list = "ls ".$self->input_dir().'/'.$self->name().'*.bed';
        my @rfiles = `$list`;
        throw("Found more than one cluster file:\n@rfiles\nNeed to implement ExperimentalSubset rollback before removing this!")
            if (scalar(@rfiles) >1);
        
        $self->result_files(\@rfiles);
    }
    
    if (scalar(@{$self->result_files()}) >1) {
        warn("Found more than one bed file:\n".
             join("\n", @{$self->result_files()})."\nBed does not yet handle replicates\n".
             "we need to resolve how we are going handle replicates with random cluster IDs");
        #do we even need to?
    }


    #how are we going to track import of files if they are being directly imported into annotated_feature?
    #Current solution is to create dummy chips, but we want something neater
    #do we need to track import as closely as with chips?
    
    foreach my $filepath(@{$self->result_files()}) {
        chomp $filepath;
        my $filename;
        my $roll_back = 0;
        ($filename = $filepath) =~ s/.*\///;
        my $sub_set;

        $self->log("Found bed file\t$filename");

        if($sub_set = $eset->get_subset_by_name($filename)){
            $roll_back = 1;
        }else{
            $sub_set = $eset->add_new_subset($filename);
        }
        
        #store if not already, skips if stored
        $eset_adaptor->store_ExperimentalSubsets([$sub_set]);

        if ($sub_set->has_status('IMPORTED')){
            $self->log("ExperimentalSubset(${filename}) has already been imported");
        } 
        else {
            $new_data = 1;

            if ($self->recovery() && $roll_back) {
                $self->log("Rolling back results for ExperimentalSubset:\t".$filename);

                warn "Cannot yet rollback for just an ExperimentalSubset, rolling back entire set\n";
                warn ("Need to implement annotated_feature rollback!\n");
                #$self->db->rollback_results($cc_id);
            }
            
            $self->log("Reading bed file:\t".$filename);
            my $fh = open_file($filepath);
            my @lines = <$fh>;
            close($fh);
            
            #my $rfile_path = $self->get_dir("norm")."/result.Parzen.".$echip->unique_id().".txt";
            #my $rfile = open_file($rfile_path, '>');
            #my $r_string = "";
            my ($line, $f_out);
            my $fasta = '';
            
            #warn "we need to either dump the pid rather than the dbID or dump the fasta in the DB dir";
            my $fasta_file = $ENV{'EFG_DATA'}."/fastas/".$self->experiment->name().'.'.$filename.'.fasta';

            if($self->dump_fasta()){
                $self->backup_file($fasta_file);
                $f_out = open_file($fasta_file, '>');
            }
            
            $self->log("Parsing file:\t$filename");

            foreach my $line (@lines) {
                $line =~ s/\r*\n//o;
                next if $line =~ /^\#/;	
                next if $line =~ /^$/;
                next unless $line =~ /^chr/i;

                my ($chr, $start, $end, $pid, $score) = split/\t/o, $line;				  
                #change from UCSC to EnsEMBL coords
                $start +=1;
                $end +=1;
                
                if(!  $self->cache_slice($chr)){
                    warn "Skipping AnnotatedFeature import, cound non standard chromosome: $chr";
                }else{
                    
                    #this is throwing away the encode region which could be used for the probeset/family?	
                    my $feature = Bio::EnsEMBL::Funcgen::AnnotatedFeature->new
                        (
                         -START         => $start,
                         -END           => $end,
                         -STRAND        => 1,
                         -SLICE         => $self->cache_slice($chr),
                         -ANALYSIS      => $fset->anal,
                         -DISPLAY_LABEL => $pid,
                         -FEATURE_SET   => $fset,
                         );
                    
                    $af_adaptor->store($feature);
                    
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


            $self->log("Finished importing:\t$filepath");
            $sub_set->adaptor->set_status('IMPORTED', $sub_set);
        }
    }

    $self->log("No new data, skipping result parse") if ! $new_data;
    
    $self->log("Finished parsing and importing results");
    
    return;
}
  


1;
