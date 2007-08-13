#!/software/bin/perl -w


####!/opt/local/bin/perl -w


=head1 NAME

ensembl-efg convert_htilist_to_features.pl
  
=head1 SYNOPSIS

convert_hitlist_to_features.pl [options]

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

use strict;

### MODULES ###
use Getopt::Long;
#use Carp;#For dev only? cluck not exported by default Remove this and implement in Helper
use Pod::Usage;
#POSIX? File stuff
use File::Path;
use Data::Dumper;
#use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Utils::Exception qw( throw warning );
use Bio::EnsEMBL::Funcgen::Utils::EFGUtils qw (open_file run_system_cmd backup_file);
use Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Funcgen::FeatureType;
use Bio::EnsEMBL::Funcgen::FeatureSet;
use Bio::EnsEMBL::Funcgen::DataSet;
use Bio::EnsEMBL::Funcgen::AnnotatedFeature;
use Bio::EnsEMBL::Funcgen::Utils::Encode qw(get_encode_regions);
use Bio::EnsEMBL::Analysis;

use Bio::EnsEMBL::DBSQL::DBAdaptor;
use strict;

$| = 1;							#autoflush
my ($pass, $dbname, $help, $man, $ftname, $file, $species, @overlap_names, %overlap_sets, %overlap_names, %union_fsets);
my ($clobber, $write_features, $dump_features, $data_version, $host, $chr_name, $default_set);
my $max_length = 4000;
my $twindow = 2500;
my ($slice_name, @slices);
my $out_dir ='.';
my $user = "ensadmin";
my $port = '3306';

$main::_debug_level = 0;
$main::_tee = 0;


GetOptions (
			"pass|p=s"       => \$pass,
			"port=s"         => \$port,
			"host|h=s"       => \$host,
			"user|u=s"       => \$user,
			"dbname|d=s"     => \$dbname,
			"species=s"      => \$species,
			"help|?"         => \$help,
			"man|m"          => \$man,
		  	"overlap_sets=s"    => \@overlap_names,
			"default_set=s"   => \$default_set,
			"max_length=s"    => \$max_length,
			"transcript_window=s" => \$twindow,
			#"no_focus"       => \$no_focus,
			#"union_sets=s"    => \@union_set_names,
			"write_features" => \$write_features,
			"dump_features"  => \$dump_features,
			"clobber"        => \$clobber,
			"data_version=s" => \$data_version,
			"out_dir|o=s"    => \$out_dir,
			#"multiplex"      => \$multiplex,
			"slice_name=s"   => \$slice_name,
			"chr_name=s"     => \$chr_name,
		   );



#@union_set_names = split/,/, join(',', @union_set_names);
@overlap_names = split/,/, join(',', @overlap_names);

pod2usage(1) if $help;
pod2usage(-exitstatus => 0, -verbose => 2) if $man;



if (! ($pass && $host && $dbname && $data_version && $species)) {
    throw("Some mandatory parameters are not set, you must specify:\n".
          "-pass\t-port\t-host\t-dbname\t-data_version\t-species");
}

throw("You must provide a default, which is used if max length $max_length bps is exceeded") if(! defined $default_set);


run_system_cmd("mkdir -p $out_dir") if(! -d $out_dir);


print "You must supply a list of overlap sets" if (! @overlap_names);


if (! ($write_features || $dump_features)) {
    print "No output type specified turning on dump_features\n";
    $dump_features = 1;
}




my $cdb = Bio::EnsEMBL::DBSQL::DBAdaptor->new(
											  #-host => 'ensembldb.ensembl.org',
											  -host => 'ens-livemirror',
											  -user => 'ensro',
											  -dbname => $species."_core_".$data_version,
											  -species => $species,
											 );



my $db = Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor->new(
													  -dbname => $dbname,
													  -port   => $port,
													  -pass   => $pass,
													  -host   => $host,
													  -user   => $user,
													  -dnadb  => $cdb,
													 );


my $fset_a = $db->get_FeatureSetAdaptor();
my $dset_a = $db->get_DataSetAdaptor();
my $anal_a = $db->get_AnalysisAdaptor();
my $ft_adaptor = $db->get_FeatureTypeAdaptor();
my $slice_a = $db->get_SliceAdaptor();
my $afa = $db->get_AnnotatedFeatureAdaptor();


my $anal = Bio::EnsEMBL::Analysis->new(
									   -logic_name      => 'RegulatoryRegion',
									   -db              => 'NULL',
									   -db_version      => 'NULL',
									   -db_file         => 'NULL',
									   -program         => 'NULL',
									   -program_version => 'NULL',
									   -program_file    => 'NULL',
									   -gff_source      => 'NULL',
									   -gff_feature     => 'NULL',
									   -module          => 'NULL',
									   -module_version  => 'NULL',
									   -parameters      => 'NULL',
									   -created         => 'NULL',
									   -description     => 'Union of co-occurences',
									   -display_label   => 'RegulatoryRegion',
									   -displayable     => 1,
									  );

$anal = $anal_a->fetch_by_dbID($anal_a->store($anal));



foreach my $name (@overlap_names) {
        
  if (! ($overlap_sets{$name} = $fset_a->fetch_by_name($name))) {
	throw("One of your specified overlap FeatureSets does not exist:\t$name");
  } else {
	$overlap_names{$overlap_sets{$name}->dbID()} = $name;
  }
}

print "Overlap FeatureSets are:\t".join("\t", sort @overlap_names)."\n";


throw("Cannot specifiy a slice name and a chr name:\t$slice_name\t$chr_name") if($slice_name && $chr_name);


if ($slice_name) {

    if ($slice_name eq 'ENCODE') {
        my $encode_regions = &get_encode_regions($cdb);
        my @encode_region_names = sort keys %{$encode_regions};
        map {push @slices, $slice_a->fetch_by_name($encode_regions->{$_});} @encode_region_names;
        #push @slices, $slice_a->fetch_by_name($encode_regions->{ENr333});
        #print scalar(@slices), "\n";
        
    } else {

        @slices = ($slice_a->fetch_by_name($slice_name));
    }

    if (! @slices) {
        throw("-slice name did not retrieve a valid slice:\t$slice_name\n");
    }

} elsif ($chr_name) {
    
    @slices = ($slice_a->fetch_by_region('chromosome', $chr_name));
    
    if (! @slices) {
        throw("-chr_name did not retrieve a valid slice:\t$chr_name\n");
    }
    
} else {
    @slices = @{$slice_a->fetch_all('toplevel')};
}

print "Building co-occurrence features on slices:\n".join("\t", map $_->seq_region_name(), @slices)."\n";


my (%starts, %ends, %vectors, %union_ftypes, $fh, $transfer, $transfer_slice);
my (@start_ids, @names, @features, $feature, $slice, $reg_feature, $last_end, $first_start);
my $ensr_cnt = 1;



my $tx_adaptor = $cdb->get_TranscriptAdaptor();
my $g_adaptor = $cdb->get_GeneAdaptor();
my $union_fset_name = 'RegulatoryFeatures';
#'Union_'.join(':', (sort(@overlap_names), 'FivePrime','ThreePrime').;
my $union_fset = get_union_FeatureSet($union_fset_name);


if ($dump_features){
  $fh = open_file($out_dir."/".$union_fset_name.".out", '>');

  print $fh "# $union_fset_name\n";
  print $fh "# Default set is $default_set\n";
  print $fh "# Max length:\t$max_length\n";
  print $fh "# 5'/3' extension for Transcript start/end overlap:\t$twindow\n";
}
  
#144962236
foreach $slice (@slices) {
  %starts = ();
  %ends = ();
  %vectors = ();


  warn "Processing slice ".$slice->name."\n";

  foreach my $feature (@{$afa->fetch_all_by_Slice($slice)}) {
	
	#skip non union sets
	next if(! exists $overlap_sets{$feature->feature_set->name()});

	#get highest end value
	($last_end) = sort{$b <=> $a} map @$_, values %ends;
	
	print join("\t",$feature->start(), $feature->end(),$feature->feature_set->dbID(),
	          $feature->display_label()), "\n";
	
	
	#next feature start after last end
	
	if ((defined $last_end) && ($feature->start() > $last_end)) {
	  
	  ($first_start) =  sort{$a <=> $b} map @$_, values %starts;
	  
	  #Too long, just use default overlap set
	  #if(($last_end - $first_start) > $max_length){
	#	@features = &build_default_features($slice, \%starts, \%ends, \%vectors);
	#  }
	#  else{
	#	#build union feature with binary strings or'd
		@features = &build_union_features($slice);
	#  }     

	  &transfer_and_write_features($slice, @features);

	  #reset hashes
	  %starts = ();
	  %ends = ();
	  %vectors = ();
	}

	push @{$starts{$feature->feature_set->name()}}, $feature->start();
	push @{$ends{$feature->feature_set->name()}}, $feature->end();
	push @{$vectors{$feature->feature_set->name()}}, $feature->display_label();
	
  }


  #deal with last union
  ($last_end) = sort{$b <=> $a} map @$_, values %ends;
  ($first_start) =  sort{$a <=> $b} map @$_, values %starts;
  
  if($first_start){
	
  #Too long, just use default overlap set
	#if(($last_end - $first_start) > $max_length){
	#  
	#  @features = &build_default_features($slice, \%starts, \%ends, \%vectors);
	#}else{
	#  #build union feature with binary strings or'd
	  @features = &build_union_features($slice);
	#
	#}     
	
	warn "doing last with slice $slice";

	&transfer_and_write_features($slice, @features);
	
  }
}



sub transfer_and_write_features{

  my ($slice, @features) = @_;

  my $transfer=0;
  

  # get slice union features need to be transfered onto
  if( $slice->start != 1 || $slice->strand != 1) {
	$transfer=1;
	$transfer_slice = $slice_a->fetch_by_region
	  (
	   $slice->coord_system->name(),
	   $slice->seq_region_name(),
	   undef, #start
	   undef, #end
	   undef, #strand
	   $slice->coord_system->version()
	  );
	#print Dumper $transfer_slice;
  } 
  
  
  if($transfer){
	$slice = $transfer_slice;
		
	foreach my $reg_feature(@features){
	  
	  $reg_feature = $reg_feature->transfer($transfer_slice);

	}
  }
  

 
  #look at tss and 3' UTR
  #need to account for strandedness of transcripts
  
  foreach my $reg_feature(@features){
	my $five_prime = 0;
	my $three_prime = 0;

	

	my $window_start = $reg_feature->start() - $twindow;
	$window_start = 1 if $window_start <1;
	
	
	my $window_end = $reg_feature->end() + $twindow;
	my $window_length = $window_end - $window_start;
	my $fslice = $slice->sub_Slice($window_start, $window_end);



	#These bits should only be used for generating class patterns
	#They should not be used to assign the class as this would
	#exclude any un-annotated genes from this classification

	foreach my $gene(@{$g_adaptor->fetch_all_by_Slice($fslice)}){


	  #add stuff for genic bit here?
	
	  foreach my $tx(@{$tx_adaptor->fetch_all_by_Gene($gene)}){
	  
		if($tx->strand == 1){
		  
		  $five_prime = 1 if($tx->start() > 0);#TSS is in sub slice
		  $three_prime = 1 if($tx->end() <= $window_length);#3' region in sub slice
		}
		else{#opposite strand
		  $three_prime = 1 if($tx->start() >0);
		  $five_prime = 1 if($tx->end() <= $window_length);
		}
	  }
	}
 
	$reg_feature->display_label($reg_feature->display_label().$five_prime.$three_prime);
					 
	if($dump_features){
	  my ($ensr, $vector) = split/:/, $reg_feature->display_label();
	  print $ensr."\t".$reg_feature->slice->seq_region_name()."\t".$reg_feature->start()."\t".
		$reg_feature->end()."\t".$vector."\n";
	}
  }

  if($write_features){
	$afa->store(@features);
  }
}




#this does not account for features which may not contain default set features

sub build_default_features{
  my ($slice, $starts, $ends, $vectors) = @_;

  my $default = $default_set;

  if(! exists $starts{$default}){
	

	if(scalar(keys(%starts)) >1 ){
	  #$default = 'Overlap_Wiggle_H3K4me3::GM06990_DNASE_IMPORT:Nessie_NG_STD_2_ctcf_ren_BR1:Wiggle_H3K20me3:Wiggle_H3K27me3:Wiggle_H3K36me3:Wiggle_H3K4me3:Wiggle_H3K79me3:Wiggle_H3K9me3';
	  $default = 'Overlap_GM06990_DNASE_IMPORT::CD4_CTCF:CD4_H2AZ:CD4_H2BK5me1:CD4_H3K27me1:CD4_H3K27me2:CD4_H3K27me3:CD4_H3K36me1:CD4_H3K36me3:CD4_H3K4me1:CD4_H3K4me2:CD4_H3K4me3:CD4_H3K79me1:CD4_H3K79me2:CD4_H3K79me3:CD4_H3K9me1:CD4_H3K9me2:CD4_H3K9me3:CD4_H3R2me1:CD4_H3R2me2:CD4_H4K20me1:CD4_H4K20me3:CD4_H4R3me2:CD4_PolII:GM06990_DNASE_IMPORT:Nessie_NG_STD_2_ctcf_ren_BR1:Wiggle_H3K27me3:Wiggle_H3K36me3:Wiggle_H3K4me3:Wiggle_H3K79me3:Wiggle_H3K9me3:Wiggle_H4K20me3';
	}else{
	  ($default) = keys(%starts);
	}
  } 

  #my $starts = $starts{$default};

  #warn "starts are $starts";

  my @starts = @{$starts{$default}};
  my @ends = @{$ends{$default}};
  my @vectors = @{$vectors{$default}};
  my @features = ();  

  for my $i(0..$#starts){
	
	#ENSR00000184895
	my $ensr_id = sprintf("ENSR%011d", $ensr_cnt);
	
	warn "$ensr_id defaulted to $default\n";

	my $reg_feature = Bio::EnsEMBL::Funcgen::AnnotatedFeature->new
	  (
	   -slice => $slice,
	   -start => $starts[$i],
	   -end   => $ends[$i],
	   -feature_set => $union_fset,
	   -strand => 0,
	   -display_label => $ensr_id.":".$vectors[$i],
	   
	  );

	
	push @features, $reg_feature;
	
	$ensr_cnt ++;
  }

  return @features;
}


sub build_union_features{
  my ($slice) = @_;
  my (@features, %ord_vector);

	
	foreach my $vector(map @$_, values %vectors){
	  
	  
	  #warn "got donor vec $vector";
	  
	  my @tmp_vector = split//, $vector;
	  
	  foreach my $i(0..$#tmp_vector){
		$ord_vector{$i} ||= $tmp_vector[$i];
	  }
	}

	
	#sort values based on keys
	my $vector = join('', (map $ord_vector{$_}, sort (keys %ord_vector)));
	
	#warn "got vector $vector";
	
	
	
	my $ensr_id  = sprintf("ENSR%011d", $ensr_cnt);
	
	#warn "got ensr $ensr_id";
	
	my $reg_feature = Bio::EnsEMBL::Funcgen::AnnotatedFeature->new
	  (
	   -slice => $slice,
	   -start => $first_start,
	   -end   => $last_end,
	   -feature_set => $union_fset,
	   -strand => 0,
	   -display_label => $ensr_id.":".$vector,
	  );

	
	push @features, $reg_feature;
	
	$ensr_cnt ++;
		
	return @features;

}


sub get_union_FeatureSet{
    my $set_name = shift;
    
    if (! exists $union_fsets{$set_name}) {
        $union_fsets{$set_name} = $fset_a->fetch_by_name($set_name);
        
        if ($union_fsets{$set_name}) {
            
            if ($clobber && $write_features) {
			  my $cs_id = $db->get_FGCoordSystemAdaptor->fetch_by_name('chromosome')->dbID();

			  my $sql = 'DELETE from annotated_feature where feature_set_id='.
				$union_fsets{$set_name}->dbID().' and coord_system_id='.$cs_id;

			  $db->dbc->do($sql) 
				or throw('Failed to roll back annotated_features for feature_set_id'.
						 $union_fsets{$set_name}->dbID());
            } 
			elsif ($write_features) {
			  throw("Their is a pre-existing FeatureSet with the name '$set_name'\n".
                      'You must specify clobber is you want to delete and overwrite all'.
                      ' pre-existing PredictedFeatures');
            }
        } else {					#generate new fset
            
            #get ftype first
            if (! exists $union_ftypes{$set_name}) {
                
                $union_ftypes{$set_name} = $ft_adaptor->fetch_by_name($set_name);
                
                if (! $union_ftypes{$set_name}) {
                    
                    $union_ftypes{$set_name} = Bio::EnsEMBL::Funcgen::FeatureType->new
                        (
                         -name        => $set_name,
                         -description => $set_name,
                         );
                    
                    ($union_ftypes{$set_name}) = @{$ft_adaptor->store($union_ftypes{$set_name})} 
                    if $write_features;
                }
            }
            
            $union_fsets{$set_name} = Bio::EnsEMBL::Funcgen::FeatureSet->new
                (
                 -analysis     => $anal,
                 -feature_type => $union_ftypes{$set_name},
                 -name         => $set_name,
                 );
            
            
            if ($write_features) {
                ($union_fsets{$set_name}) = @{$fset_a->store($union_fsets{$set_name})};
                
                #generate data_set here too
                my $dset = Bio::EnsEMBL::Funcgen::DataSet->new
                    (
                     -feature_set => $union_fsets{$set_name},
                     -name        => $set_name,
                     );
                
                $dset_a->store($dset);
            }
        }
    }
    return $union_fsets{$set_name};
}

1;
