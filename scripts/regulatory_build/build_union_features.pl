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
use Bio::EnsEMBL::Funcgen::PredictedFeature;
use Bio::EnsEMBL::Funcgen::Utils::Encode qw(get_encode_regions);
use Bio::EnsEMBL::Analysis;

use Bio::EnsEMBL::DBSQL::DBAdaptor;
use strict;

$| = 1;							#autoflush
my ($pass, $dbname, $help, $man, $ftname, $file, $species, @focus_set_names, @union_set_names);
my ($clobber, $write_features, $dump_features, $data_version, $host, $multiplex, $chr_name);
my (%focus_sets, %focus_names, %union_sets, %union_names, $slice_name, $no_focus, @slices);
#my $reg = "Bio::EnsEMBL::Registry";
my $out_dir ='.';
#my $data_dir = $ENV{'EFG_DATA'};
my $user = "ensadmin";
#my $host = 'localhost';
my $port = '3306';

#Definitely need some sort of Defs modules for each array?

$main::_debug_level = 0;
$main::_tee = 0;

#Use some sort of DBDefs for now, but need  to integrate with Register, and have put SQL into (E)FGAdaptor?
#Use ArrayDefs.pm module for some of these, class, vendor, format?
#ArrayDefs would also contain paths to data and vendor specific parse methods?

GetOptions (
			"pass|p=s"       => \$pass,
			"port=s"         => \$port,
			"host|h=s"       => \$host,
			"user|u=s"       => \$user,
			"dbname|d=s"     => \$dbname,
			"species=s"      => \$species,
			"help|?"         => \$help,
			"man|m"          => \$man,
		  	"focus_sets=s"    => \@focus_set_names,
			"no_focus"       => \$no_focus,
			"union_sets=s"    => \@union_set_names,
			"write_features" => \$write_features,
			"dump_features"  => \$dump_features,
			"clobber"        => \$clobber,
			"data_version=s" => \$data_version,
			"out_dir|o=s"    => \$out_dir,
			"multiplex"      => \$multiplex,
			"slice_name=s"   => \$slice_name,
			"chr_name=s"     => \$chr_name,
		   );



@union_set_names = split/,/, join(',', @union_set_names);
@focus_set_names = split/,/, join(',', @focus_set_names);

pod2usage(1) if $help;
pod2usage(-exitstatus => 0, -verbose => 2) if $man;



if (! ($pass && $host && $dbname && $data_version && $species)) {
    throw("Some mandatory parameters are not set, you must specify:\n".
          "-pass\t-port\t-host\t-dbname\t-data_version\t-species");
}


run_system_cmd("mkdir -p $out_dir") if(! -d $out_dir);


if (! $no_focus) {
    throw('You must specificy some focus feature sets to build on') if(! @focus_set_names);
} elsif (@focus_set_names) {
    throw('You cannot specify -no_focus and -focus_sets');
}
print 'No union sets specified, using all feature set on $data_version' if (! @union_set_names);


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
my $pfa = $db->get_PredictedFeatureAdaptor();


my $anal = Bio::EnsEMBL::Analysis->new(
									   -logic_name      => 'Co-occurrence',
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
									   -description     => 'Co-occurrence of FeatureTypes',
									   -display_label   => 'Co-occurrence',
									   -displayable     => 1,
									  );

$anal_a->store($anal);


print "Running in multiplex mode\n" if $multiplex;

if (! $no_focus) {
    
    foreach my $name (@focus_set_names) {
        
        if (! ($focus_sets{$name} = $fset_a->fetch_by_name($name))) {
            throw("One of your sepcified focus FeatureSets does not exist:\t$name");
        } else {
            $focus_names{$focus_sets{$name}->dbID()} = $name;
        }
    }
} else {
    print "No focus, using all available feature sets\n";
    
    foreach my $fset (@{$fset_a->fetch_all()}) {
        push @focus_set_names, $fset->name();
        
        $focus_sets{$fset->name()} = $fset;
        $focus_names{$fset->dbID()} = $fset->name();
    }
}

print "Focus FeatureSets are:\t".join("\t", sort values %focus_names)."\n";

if ($no_focus) {
    print "Setting union sets to focus sets";
    %union_sets = %focus_sets;
    %union_names = %focus_names;
    @union_set_names = @focus_set_names;
} else {
    if (! @union_set_names) {
        print "No union FeatureSets specified, using all FeatureSets\n";
        
        foreach my $fset (@{$fset_a->fetch_all()}) {
            
            #next if(grep($fset->name(), @focus_set_names));
            push @union_set_names, $fset->name();
            
            $union_sets{$fset->name()} = $fset;
            $union_names{$fset->dbID()} = $fset->name();
            
        }
    } else {
        
        foreach my $name (@union_set_names) {
            
            if (! ($union_sets{$name} = $fset_a->fetch_by_name($name))){
                throw("One of your sepcified union FeatureSets does not exist:\t$name");
            } else {
                $union_names{$union_sets{$name}->dbID()} = $name;
            }
        }
        
        #print "now add focus sets to union sets\n";
        foreach my $name (@focus_set_names) {
            if (! exists $union_sets{$name}){
                $union_names{$focus_sets{$name}->dbID()} = $name;
                $union_sets{$name} = $focus_sets{$name};
            }
        }
    }
}

print "Union FeatureSets are:\t".join("\t", sort values %union_names)."\n";

throw("Cannot specifiy a slice name and a chr name:\t$slice_name\t$chr_name") if($slice_name && $chr_name);


if ($slice_name) {

    if ($slice_name eq 'ENCODE') {
        my $encode_regions = &get_encode_regions($cdb);
        my @encode_region_names = sort keys %{$encode_regions};
        map {push @slices, $slice_a->fetch_by_name($encode_regions->{$_});} @encode_region_names;
        #push @slices, $slice_a->fetch_by_name($encode_regions->{ENr333});
        #print scalar(@slices), "\n";
        
    } else {

        @slices = $slice_a->fetch_by_name($slice_name);
    }

    if (! @slices) {
        throw("-slice name did not retrieve a valid slice:\t$slice_name\n");
    }

} elsif ($chr_name) {
    
    @slices = $slice_a->fetch_by_region('chromosome', $chr_name);
    
    if (! @slices) {
        throw("-chr_name did not retrieve a valid slice:\t$chr_name\n");
    }
    
} else {
    @slices = @{$slice_a->fetch_all('toplevel')};
}

print "Building co-occurrence features on slices:\n".join("\t", map $_->seq_region_name(), @slices)."\n";


my (%starts, %ends, %current_unions, %union_fsets, %union_ftypes, %file_handles, %union_cnts);
my (@start_ids, @names);
#my $current_pos = 1;

foreach my $slice (@slices) {

    my $transfer=0;
    my $transfer_slice;
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

    foreach my $feature (@{$pfa->fetch_all_by_Slice($slice)}) {
        
        #skip non union sets
        next if(! exists $union_sets{$feature->feature_set->name()});

        my ($first_end) = sort{$b <=> $a} values %ends;

        print join("\t",$feature->start(), $feature->end(),$feature->feature_set->dbID(),
                   $feature->display_label()), "\n";

        

        if ((defined $first_end) && ($feature->start() > $first_end)) {
            #build co-occurence feature for feature with first_end and 
            #flush start end values appropriately
            print "***first_end ($first_end) < feature_start (".$feature->start().")\n";
            
            #warn "first end is $first_end and next start is ".$feature->start()." ".$feature->feature_set->name();
            
            foreach my $end_id (keys %ends) {
                print "\tend_id:\t$end_id\t$starts{$end_id}\t$ends{$end_id}\t$union_names{$end_id}\n";
                
                if ($ends{$end_id} == $first_end) { #find currently ending/ed features
                    #print "\tend{$end_id} $ends{$end_id} == first_end $first_end\n";

                    #now compare other start end vals of focus sets or all 
                    #if end_id is an focus set or no_focus sets specified
                    
                    if(scalar (keys %ends) > 1) {
                        #got over lap
                        
                        #warn "got overlap with ended feature ".$union_names{$end_id}." $first_end";
                        
                        @start_ids = sort { $starts{$a} <=> $starts{$b} } keys %starts;
                        #print "\tstart_ids by start: @start_ids \n";
                        
                        if (exists $focus_names{$end_id}) {
                            
                            print "\t*Focus feature $focus_names{$end_id}\n";
                            print "\tend_id:\t\t$end_id\t$starts{$end_id}\t$ends{$end_id}\t$union_names{$end_id}\n";

                            #compare vs all union sets next if end_id is same
                            
                            #we need to check whether we're duplicating last feature for a given co-occurence set
                            #will this not be solved by cleaning the start ends as we go move along the features?
                            #features should be named by sorting feature type names of co-occuring features
                            #need to build current co-occur feature hash to stop duplicates being formed
                            #we need to build these in order of starts to get multiplex co-occurences
                            
                            @names = ($focus_names{$end_id});
                            
                            foreach my $start_id (@start_ids) {
                                
                                next if $start_id == $end_id;
                                next if $starts{$end_id} > $ends{$start_id};
                                print "\t start_id:\t$start_id\t$starts{$start_id}\t$ends{$start_id}\t$union_names{$start_id}\n";
                                
                                if ($multiplex) {
                                    push @names, $union_names{$start_id};
                                } else {
                                    @names = ($focus_names{$end_id}, $union_names{$start_id});
                                }
                                
                                print "\tFound union of @names\n";
                                
                                @names = sort @names;
                                
                                #shift names to move thro' all possible co-occurences
                                
                                while (scalar(@names) > 1) {
                                    
                                    my $union_set_name = join(':', @names);
                                    
                                    my ($start,$end) = (sort {$a<=>$b} ($starts{$start_id},
                                                                        $ends{$start_id},
                                                                        $starts{$end_id},
                                                                        $ends{$end_id}))[1,2];
                                        
                                    my $union_feature = Bio::EnsEMBL::Funcgen::PredictedFeature->new
                                        (
                                         -slice  => $slice,
                                         -start  => $start,
                                         -end    => $end,
                                         -strand => 0,
                                         -feature_set => &get_union_FeatureSet($union_set_name),	
                                         );
                                    
                                    if ($transfer) {
                                        #warn("original uf:\t", join("\t", $union_feature->start, $union_feature->end), "\n");
                                        $union_feature = $union_feature->transfer($transfer_slice);
                                        #warn("transfered uf:\t", join("\t", $union_feature->start, $union_feature->end), "\n");
                                    }
                                    
                                    shift @names;
                                    
                                    if (! exists $current_unions{$union_set_name}) {
                                        $current_unions{$union_set_name} = $union_feature;
                                    } else {
                                        #this should happen as we're deleting the A start end values after AB, 
                                        #before we have a chance to build BA
                                        warn "Found duplicate union sets for $union_set_name at:\t".$slice->name;
                                    }
                                }
                            }
                            
                        } else {

                            print "\t*Union feature $union_names{$end_id}\n";
                            print "\tend_id:\t\t$end_id\t$starts{$end_id}\t$ends{$end_id}\t$union_names{$end_id}\n";

                            #just compare vs focus_sets
                            #no need to next as we know it's not present in focus_names
                            
                            @names = ($union_names{$end_id});
                            
                            
                            #warn "Non focus feature @names";
                            
                            
                            foreach my $start_id (@start_ids) {
                                #next if $start_id == $end_id;
                                
                                if ($focus_names{$start_id} && $ends{$start_id} > $starts{$end_id}) {
                                    print "\t start_id:\t$start_id\t$starts{$start_id}\t$ends{$start_id}\t$union_names{$start_id}\n";
                                    #print "focus_name $focus_names{$start_id} exists\n";
                                    
                                    if ($multiplex) {
                                        push @names, $focus_names{$start_id};
                                    } else {
                                        @names = ($union_names{$end_id}, $focus_names{$start_id});
                                    }
                                    
                                    @names = sort @names;
                                    
                                    print "\tFound union of @names\n";
                                    #warn "Found union of  @names";
                                    #exit;
                                    
                                    #shift names to move thro' all possible co-occurences
                                    
                                    while (scalar(@names) > 1) {
                                        
                                        my $union_set_name = join(':', @names);
                                        #print $union_set_name, "\n";

                                        my ($start,$end) = (sort {$a<=>$b} ($starts{$start_id},
                                                                            $ends{$start_id},
                                                                            $starts{$end_id},
                                                                            $ends{$end_id}))[1,2];
                                        
                                        my $union_feature = Bio::EnsEMBL::Funcgen::PredictedFeature->new
                                            (
                                             -slice  => $slice,
                                             -start  => $start,
                                             -end    => $end,
                                             -strand => 0,
                                             -feature_set => &get_union_FeatureSet(\@names),	
                                             );
                                        

                                        if ($transfer) {
                                            #warn("original uf:\t", join("\t", $union_feature->start, $union_feature->end), "\n");
                                            $union_feature = $union_feature->transfer($transfer_slice);
                                            #warn("transfered uf:\t", join("\t", $union_feature->start, $union_feature->end), "\n");
                                        }
                                        
                                        shift @names;
                                        
                                        if (! exists $current_unions{$union_set_name}) {
                                            $current_unions{$union_set_name} = $union_feature;
                                        } else {
                                            #this should happen as we're deleting the A start end values after AB,
                                            #before we have a chance to build BA
                                            warn "Found duplicate union sets for $union_set_name at:\t".$slice->name;
                                        }
                                    }
                                }
                            }
                        }
                        #end of if exists focus_names
                    }

                    #remove feature start ends for expired feature
                    delete $starts{$end_id};
                    delete $ends{$end_id};
                }
            }
            
            #print and load here as required
            #we could cache here by pusing the hash
            
            foreach my $uset_name (keys %current_unions) {

                $union_cnts{$uset_name} ||= 0;
                $union_cnts{$uset_name} ++;
                
                $pfa->store($current_unions{$uset_name}) if $write_features;
                
                if ($dump_features) {
                    
                    #warn "dumping $uset_name features";
                    
                    
                    #we should test for file and clobber here too
                    
                    $file_handles{$uset_name} ||= open_file($out_dir."/${uset_name}.hitlist", '>');
                    
                    my $fh = $file_handles{$uset_name};
                    
                    print $fh 'chr'.$current_unions{$uset_name}->slice->seq_region_name()."\t".
                        $current_unions{$uset_name}->start()."\t".$current_unions{$uset_name}->end()."\t\n";
                    print "chr".$current_unions{$uset_name}->slice->seq_region_name()."\t".
                        $current_unions{$uset_name}->start()."\t".$current_unions{$uset_name}->end()."\t$uset_name\n";
                    #score field empty
                    
                }
            }
            
            %current_unions = ();

        }
        
        $starts{$feature->feature_set->dbID()} = $feature->start();
        $ends{$feature->feature_set->dbID()} = $feature->end();

    }
}

foreach my $uset_name(keys %union_cnts){
    print 'Found '.$union_cnts{$uset_name}." union features of type $uset_name\n";
}


foreach my $fh (values %file_handles) {
    close($fh);
}


sub get_union_FeatureSet{
    my $set_name = shift;
    
    if (! exists $union_fsets{$set_name}) {
        $union_fsets{$set_name} = $fset_a->fetch_by_name($set_name);
        
        if ($union_fsets{$set_name}) {
            
            if ($clobber && $write_features) {
			  my $cs_id = $db->get_FGCoordSystemAdaptor->fetch_by_name('chromosome')->dbID();

			  my $sql = 'DELETE from predicted_feature where feature_set_id='.
				$union_fsets{$set_name}->dbID().' and coord_system_id='.$cs_id;

			  $db->dbc->do($sql) 
				or throw('Failed to roll back predicted_features for feature_set_id'.
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
                         -description => "Co-occurence of features $set_name",
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
				 -feature_class=> 'annotated',
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
