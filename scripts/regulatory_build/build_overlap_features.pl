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

build_overlap_features.pl

=head1 SYNOPSIS

build_overlap_features.pl -host host -P 3306 -u user -d password 
    -o /tmp -f testA -t testB,testC

=head1 DESCRIPTION

Generates regulatory features based on co-occurance between a 
given focus feature set and target sets. It calculates the
intersection between each focus-target-pair and generates for 
each focus features a binary vector representing the co-occurence.

=cut


### NJ added clobber based on slice seq_region_id and current coord_system_id
### NJ added various comments you may read or ignore as you see fit
### NJ Do we need to look at what and where the STDOUT is going
### There are log messages and data going to STDOUT, STDERR and OUT

use strict;
use warnings;
use Data::Dumper;
use Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor;
use Bio::EnsEMBL::DBSQL::DBAdaptor;


$| = 1;

### OLD ###
#my %opts;
#getopts('h:P:u:p:d:f:t:s:cwDo:iS:v:', \%opts);

#H was for ?
#f is focus set
#t is target sets comma separated list
#D is dump, but print to STDERR and STDOUT(actaully current default filehandle)

use Getopt::Long;
my ($pass,$port,$host,$user,$dbname,$species,$help,$man,
    $data_version,$outdir,$do_intersect,$write_features,$seq_name,$clobber,
    $focus,$target,$dump,$debug);
GetOptions (
            "pass|p=s"       => \$pass,
            "port=s"         => \$port,
            "host|h=s"       => \$host,
            "user|u=s"       => \$user,
            "dbname|d=s"     => \$dbname,
            "species=s"      => \$species,
            "help|?"         => \$help,
            "man|m"          => \$man,
            "data_version|v=s" => \$data_version,
            "outdir|o=s"     => \$outdir,
            "do_intersect|i=s" => \$do_intersect,
            "write_features|w" => \$write_features,
            "seq_name|s" => \$seq_name,
            "clobber" => \$clobber,
            "focus|f=s" => \$focus,
            "target|t=s" => \$target,
            "dump|d" => \$dump,
            "debug" => \$debug
            );

### defaults ###
$port = 3306 if !$port;
$species = 'homo_sapiens' if !$species;
$debug = 0;

### check options ###

throw("Must specify mandatory database hostname (-host).\n") if ! defined $host;
throw("Must specify mandatory database username. (-user)\n") if ! defined $user;
throw("Must specify mandatory database password (-pass).\n") if ! defined $pass;
throw("Must specify mandatory database name (-dbname).\n") if ! defined $dbname;
throw("Must specify mandatory database data version, like 47_36i (-data_version).\n") 
    if !$data_version;

throw("Must specify mandatory focus sets (-focus).\n") if ! defined $focus;
throw("Must specify mandatory target sets (-target).\n") if ! defined $target;

warn("write_features not set: do_intersect will be skipped\n") 
    if ($do_intersect && ! $write_features);

#throw("No output directory specified! Use -o option.") if (!$outdir);
if (defined $outdir && ! -d $outdir) {
    system("mkdir -p $outdir");
}
$outdir =~ s/\/$//;

### ChipSeq stuff
my %ChIPseq_cutoff = (
                      ### cutoff T/O <= 2
                      'CD4_CTCF'=>        5,
                      'CD4_H3K27me3'=>    8,
                      'CD4_H3K36me3'=>    4,
                      'CD4_H3K4me3'=>     6,
                      'CD4_H3K79me3'=>   22,
                      'CD4_H3K9me3'=>     7,
                      'CD4_H4K20me3'=>   17,
                      ### cutoff ~ <= 25000
                      ### see /lustre/work1/ensembl/graef/efg/input/SOLEXA/LMI/data/*.clstr.cutoff_25000.dat
                      'CD4_H2AZ'=>       15,
                      'CD4_H2BK5me1'=>   16,
                      'CD4_H3K27me1'=>    6,
                      'CD4_H3K27me2'=>    5,
                      'CD4_H3K36me1'=>    4,
                      'CD4_H3K4me1'=>    31,
                      'CD4_H3K4me2'=>    10,
                      'CD4_H3K79me1'=>    5,
                      'CD4_H3K79me2'=>    4,
                      'CD4_H3K9me1'=>    12,
                      'CD4_H3K9me2'=>     5,
                      'CD4_H3R2me1'=>     5,
                      'CD4_H3R2me2'=>     5,
                      'CD4_H4K20me1'=>   30,
                      'CD4_H4R3me2'=>     4,
                      'CD4_PolII'=>       8
                      );


use Bio::EnsEMBL::Utils::Exception qw(verbose throw warning info);
use Bio::EnsEMBL::Analysis::Tools::Logger qw(logger_verbosity logger_info);

my $utils_verbosity = 'WARNING';
my $logger_verbosity = 'OFF';
verbose($utils_verbosity);
logger_verbosity($logger_verbosity);

# databases and adaptors

# use ensembldb as we may want to use an old version

my $cdb = Bio::EnsEMBL::DBSQL::DBAdaptor->new(
                                              -host => 'ensembldb.ensembl.org',
                                              -port => 3306,
                                              -user => 'anonymous',
                                              #-host => '127.0.0.1',
                                              #-port => 33064,
                                              #-user => 'ensro',
                                              -dbname => $species.'_core_'.$data_version,
                                              -species => $species,
                                              );

my $db = Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor->new
    (
     -host   => $host,
     -user   => $user,
     -dbname => $dbname,
	   -species => $species,
     -pass   => $pass,
     -port   => $port,
     -dnadb  => $cdb
     );
#print Dumper $db;

my $fsa = $db->get_FeatureSetAdaptor();
my $fta = $db->get_FeatureTypeAdaptor();
my $afa = $db->get_AnnotatedFeatureAdaptor();
my $sa = $db->get_SliceAdaptor();
my $aa = $db->get_AnalysisAdaptor();

# parse focus and target sets and check that they exist
my @fset_ids;

my %focus_fsets;
map { $focus_fsets{$_} = $fsa->fetch_by_name($_); 
      throw("Focus set $_ does not exist in the DB") 
          if (! defined $focus_fsets{$focus}); 
      push(@fset_ids, $focus_fsets{$_}->dbID);
  } split(',', $focus);


my %target_fsets;
map { $target_fsets{$_} = $fsa->fetch_by_name($_); 
      throw("Target set $_ does not exist in the DB") 
          if (! defined $target_fsets{$_}); 
      push(@fset_ids, $target_fsets{$_}->dbID);
      } split(',', $target);

# make sure that target sets also contain focus sets (Do we really need this?)
map { $target_fsets{$_} = $fsa->fetch_by_name($_) } keys %focus_fsets;


# retrieve sequence to be analyzed 
my $slice;

if ($seq_name) {
  $slice = $sa->fetch_by_region('chromosome', $seq_name);
} elsif (defined $ENV{LSB_JOBINDEX}) {
  warn "Performing whole genome analysis on toplevel slices using the farm.\n";
    
  my $toplevel = $sa->fetch_all('toplevel');
  my @chr = sort (map $_->seq_region_name, @{$toplevel});
  #print Dumper @chr;
  
  my @slices;
  foreach my $chr (@chr) {
      
      next if ($chr =~ m/^NT_/);
      
      push @slices, $sa->fetch_by_region('chromosome', $chr);
  }
  #print Dumper @slices;

  $slice=$slices[$ENV{LSB_JOBINDEX}-1];      
  print Dumper ($ENV{LSB_JOBINDEX}, $slice->name);
  
} else {

    throw("Must specify mandatory chromosome name (-seq_name) or set\n ".
          "LSF environment variable LSB_JOBINDEX to perform whole\n ".
          "genome analysis on toplevel slices using the farm.\n");

}

#should this be printing to OUT or STDERR?
#or remove as we're printing this later?
print STDERR '# Focus set: ', join(" ", keys %focus_fsets), "\n",
    '# Target set(s): ', join(" ", sort keys %target_fsets), "\n",
    '# Species: ', $species, "\n",
    '# Chromosome: ', join(" ",$slice->display_id(),$slice->dbID), "\n";

if ($dump) {
    
    print STDERR "# Dumping annotated features from database to file.\n";
    print STDERR "# This will delete existing file dumps of that data version!\n";

    throw("Must specify directory to write the output (-outdir).\n") if ! defined $outdir;
    print STDERR "# Output goes to ", $outdir, "\n";


    my $sql = "select seq_region_id,seq_region_start,seq_region_end,".
        "seq_region_strand,feature_set_id from annotated_feature";
    my $command = "echo \"$sql\" ".
        " | mysql -h".$host." -P".$port." -u".$user." ".$dbname.
        " | gawk '{if (\$5==".join("||\$5==", @fset_ids).") print }'".
        " | sort -n -k 1,1 -k 2,2 -k3,3".
        " | gawk '{ print >> \"".$outdir."/".$dbname.".annotated_feature_\" \$1 \".txt\"}'";
        #." > ".$outdir."/".$dbname.".annotated_feature.txt";

    print STDERR "# Execute: ".$command."\n";
    
    # need to remove existing dump files, since we append to the file
    system("rm -f $outdir/$dbname.annotated_feature_*.txt") &&
        throw ("Can't remove files");

    system($command) &&
        throw ("Can't dump data to file $outdir/$dbname.annotated_feature.txt");

}

#my (@starts, @ends, @fset_ids, @features,
my (@starts, @ends, @features,
    $focus_start, $focus_end, $focus_fset_id, $focus_score, %cooc_fsets,
    @overlap_features, $overlap_fset, $intersect_fset, %ftypes);


#this next needs altering to accomodate multple focus sets
#i.e. we need to create multiple fsets, one for each focus

#set up new fsets here
if($write_features){
	my (@tfset_names, @ftype_names);
	
	#get output set string and feature types
	
	
	foreach my $target_fset(values %target_fsets){
		
		my $is_focus = 0;
		
		#foreach my $focus_fset(values %focus_fsets){
		#  $is_focus = 1 if $focus_fset->name() eq $target_fset->name();
		#}
		
		#next if $is_focus;
		
		my $ftype = $target_fset->feature_type();
		
		if(! defined $ftype){
			throw("You have a target FeatureSet with no associated FeatureType\n".
				  "You must associat a FeatureType to ensure valid FeatureSet generation");
		}
		
		push @ftype_names, $ftype->name();
		push @tfset_names, $target_fset->name();
	}
	
	
	warn "hardcoding for only one focus set!";
	
	if (scalar (keys %focus_fsets) > 1){
		throw('THis script does not yet accomodate multiple focus sets');
	}
	
	
	
	
	#these sorts may give different orders
	my $tfset_names = join('', keys(%focus_fsets))."::".join(':', sort(@tfset_names));
	my $tfset_ftype_names = join(':', sort(@ftype_names));
	
	
	#set up default overlap fset and analysis
	my $analysis = Bio::EnsEMBL::Analysis->new
		(
		 -logic_name      => 'Co-occurrence Overlap',
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
		 -description     => 'Co-occurrence of overlapping feature types',
		 -display_label   => 'Co-occurrence Overlap',
		 -displayable     => 1
		);
	
	$aa->store($analysis) if $write_features;
	
	$overlap_fset = get_FeatureSet('Overlap_'.$tfset_names, $tfset_ftype_names, $analysis);	
	
	
	#set up intersect analysis and fset
	if($do_intersect){
		
		$analysis = Bio::EnsEMBL::Analysis->new
			(
			 -logic_name      => 'Co-occurrence Intersect',
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
			 -description     => 'Intersect of co feature types',
			 -display_label   => 'Co-occurrence Intersect',
			 -displayable     => 1
			);
		
		$aa->store($analysis) if $write_features;;
		
		
		$intersect_fset = get_FeatureSet('Intersect_'.$tfset_names, $tfset_ftype_names, $analysis);	
		
	}
	
}

# get list of FeatureSet objects
my @target_fsets;
map {push @target_fsets, $target_fsets{$_} } keys %target_fsets;

# compare all target features against the focus features
foreach my $af (@{$afa->fetch_all_by_Slice_FeatureSets($slice, \@target_fsets)}) {
#foreach my $af (@{$afa->fetch_all_by_Slice($slice)}) {
    
    print STDERR join(" ", $af->start, $af->end, $af->score?$af->score:'',
					  $af->feature_set->dbID(),
					  $af->feature_set->name()), "\n" if ($debug);

    next if(! exists $target_fsets{$af->feature_set->name()});
    
    #shoud this print to STDERR?
    print STDERR join(" ", $af->start, $af->end, $af->score?$af->score:'',
					  $af->feature_set->dbID(),
					  $af->feature_set->name()), "\n" if ($debug);

    if (exists $focus_fsets{$af->feature_set->name()}) {
        
        ### fix revisting / counting feature twice here (focus testC vs. target testA) ###
        if (defined $focus_end){ # && !exists $regulatory_features{$focus_start}{$focus_end}) {
            
            push @overlap_features, @{&get_overlap_features(\@starts, \@ends, \@fset_ids, 
                                                            $focus_start, $focus_end, 
															$focus_fset_id, $focus_score)};
            #print Dumper $overlap_features;
            #&get_overlap_string(\%regulatory_features, \@overlap_features, $focus_start, $focus_end);
            
        }

        # focus feature set
        $focus_start = $af->start;
        $focus_end = $af->end;
        $focus_fset_id = $af->feature_set->dbID();
		$focus_score = $af->score?$af->score:'';

    } elsif (defined $focus_end && $af->start > $focus_end) {

        #warn ("focus: ".$focus_start." ".$focus_end." ".$focus_fset_id);

        # target feature sets
        push @overlap_features, @{&get_overlap_features(\@starts, \@ends, \@fset_ids, 
                                                        $focus_start, $focus_end,
														$focus_fset_id, $focus_score)};
        #print Dumper $overlap_features;
        #&get_overlap_string(\%regulatory_features, \@overlap_features, $focus_start, $focus_end);

    }
    
    # quick hack for 2nd reg. build version. need to disregard ChIPseq 
    # features below a certain threshold defined hardcoded in ChIPseq_cutoff hash.
    # This works here fir target set, focus set need to be filtered elsewhere.
    next if (exists $ChIPseq_cutoff{$af->feature_set->name()} &&
             $af->score() < $ChIPseq_cutoff{$af->feature_set->name()});
    
    print STDERR join(" ", $af->start, $af->end, $af->score?$af->score:'',
					  $af->feature_set->dbID(),
					  $af->feature_set->name()), "(filtered)\n" if ($debug);
    
    push(@starts, $af->start());
    push(@ends, $af->end());
    push(@fset_ids, $af->feature_set()->dbID());

}
if ($debug) {
    print STDERR join("\t", @starts), "\n";
    print STDERR join("\t", @ends), "\n";
    print STDERR join("\t", @fset_ids), "\n";
}


push @overlap_features, @{&get_overlap_features(\@starts, \@ends, \@fset_ids, 
                                                $focus_start, $focus_end, 
												$focus_fset_id, $focus_score)};
#print Dumper $overlap_features;
		  #&get_overlap_string(\%regulatory_features, \@overlap_features, $focus_start, $focus_end);
		  


# Dump/write output
my $regulatory_features = &get_overlap_strings(\@overlap_features);
#print Dumper $regulatory_features;

my $i = 0;

my $outfile = $focus.'_chr'.$slice->seq_region_name.'.overlap';

print "Output goes to $outdir/$outfile\n";
open(OUT, "> $outdir/$outfile")
    or throw("Can't open file $outdir/$outfile");
		  
print OUT '# Focus set: ', join(" ", keys %focus_fsets), "\n";
		  print OUT '# Target set(s): ', join(" ", sort keys %target_fsets), "\n";
		  print OUT '# ', $slice->name(), "\n";


		  


#sort starts
foreach my $s (sort {$a<=>$b} keys %{$regulatory_features}) {

  #for each end (e) of feature with start s
  foreach my $e (keys %{$regulatory_features->{$s}}) {
	my ($overlap_string, $count_string);
	
	
	foreach (sort keys %{$regulatory_features->{$s}->{$e}}) {
		next if (/^score$/);
		$overlap_string .= ($regulatory_features->{$s}->{$e}->{$_}>0) ? 1 : 0;
		$count_string .= $regulatory_features->{$s}->{$e}->{$_};
	}
	

	#should probably remove focus set name from stable_id
	#need to ensure we have unique ensr IDs
	#short term solution:
	#  generate based on (chr number - 1) * 100000?
	#we should get more than 100000 features per chr?
	#are these actaully the final ensr features???
	#maybe we should just load them as the final set, then do an 
	#update based on the feature_ids for the final ensr set

	print OUT join("\t",
				   sprintf("ENSR_%s_%06d",$focus,++$i),
				   $slice->seq_region_name, $s, $e,
				   $overlap_string, $regulatory_features->{$s}->{$e}->{'score'}
				   #$count_string
				  ), "\n";


	if($write_features){

	  my $feature = Bio::EnsEMBL::Funcgen::AnnotatedFeature->new
            (
             -slice  => $slice,
             -start  => $s,
             -end    => $e,
             -strand => 0,
             -feature_set => $overlap_fset,
             -display_label => $overlap_string,#omitted ensr id for now as we will generate these later
			);
	  
	  push @features, $feature;


	  #we should be able to generate the intersects here too!
	  #this would prevent having to loop through the features twice

	}
  }
}

$afa->store(@features) if $write_features;

close OUT;


#write intersect features
if($write_features && $do_intersect){
  @features = ();

  #build a hash to map the fset ids to the correct postion in the binary vector
  #based on a sort of the fset names
  

  #could do this with some funky map hash
  #but this is easier to understand

  my @target_fsets = sort($_->name(), values %target_fsets);
  my %fset_id_pos;

  #build lookup of dbID to position in vector
  foreach my $i(0..$#target_fsets){
	$fset_id_pos{$target_fsets[$i]->dbID()} = $i;
  }

  foreach my $ft (@overlap_features) {
        
	#build vector of 0 for appropriate length

	warn "check vector map is correct";

	my @vector = map(0, keys %fset_id_pos);

	my ($start, $end, $fset_id, 
		$focus_start, $focus_end, $focus_fset_id) = @$ft;
	
	next if ($focus_fset_id == $fset_id);
	
	#set the intersecting features to 1
	$vector[$fset_id_pos{$fset_id}] = 1;
	$vector[$fset_id_pos{$focus_fset_id}] = 1;


	#take middle values to get intersect
    my ($i_start,$i_end) = (sort {$a<=>$b} (
								$focus_start,
								$focus_end,
								$start,
								$end))[1,2];


	my $feature = Bio::EnsEMBL::Funcgen::PredictedFeature->new(
		-slice  => $slice,
		-start  => $i_start,
		-end    => $i_end,
		-strand => 0,
		-feature_set => $intersect_fset,
		-display_label => join('', @vector), #need to build binary string here
		);
	push @features, $feature;
  }

  $afa->store(@features) if $write_features;;

}






sub get_overlap_features{

  my ($starts, $ends, $fset_ids, 
	  $focus_start, $focus_end, 
	  $focus_fset_id, $focus_score) = @_;

  print STDERR $focus_start.':'.$focus_end."\t".$focus_score." (focus)\n" if ($debug);
  
  my @overlap_features = ();
  my ($start, $end, $fset_id);
  my (@s, @e, @f);
  while ($#starts>=0) {
	
	$start = shift @$starts;
        $end = shift @$ends;
        $fset_id = shift @$fset_ids;
        
        #don't compare the focus_set with itself
        #next if ($focus_fset_id == $fset_id);

        print STDERR $start.':'.$end."\n" if ($debug);
        
        if ($focus_start < $end && $focus_end > $start) {

            print STDERR "Found intersection: ".
                join(":", sort {$a<=>$b} ($focus_start, $focus_end, $start,$end)).
                " with ". $fsa->fetch_by_dbID($fset_id)->name()."\n" if ($debug);

            push @overlap_features, [ $start, $end, $fset_id, 
                                      $focus_start, $focus_end, $focus_fset_id, $focus_score];

        }
        if ($debug) {
            print STDERR join("\t", @$starts), "\n";
            print STDERR join("\t", @$ends), "\n";
            print STDERR join("\t", @$fset_ids), "\n";
        }

        if ($focus_end < $end) {
            push(@s, $start);
            push(@e, $end);
            push(@f, $fset_id);
        }

    }

    push @$starts, @s;
    push @$ends, @e;
    push @$fset_ids, @f;

    if ($debug) {
        print STDERR join("\t", @$starts), "\n";
        print STDERR join("\t", @$ends), "\n";
        print STDERR join("\t", @$fset_ids), "\n";
        print STDERR "=====\n";
    }

    return \@overlap_features;

}

##########################
sub get_overlap_strings()
##########################
{
    my ($overlap_features) = @_;
    #print Dumper $overlap_features;
    my (%regulatory_features, %overlap_string);

    foreach my $ft (@$overlap_features) {

		#print Dumper $ft;
        my ($start, $end, $fset_id, 
            $focus_start, $focus_end, $focus_fset_id, $focus_score) = @$ft;

        # initialize overlap string if not exists
        if (! exists $regulatory_features{$focus_start}{$focus_end}) {
            foreach my $name (keys %target_fsets) {
                #next if (exists $focus_fsets{$name});
                $regulatory_features{$focus_start}{$focus_end}{$name} = 0;
            }
        }

        my $fset = $fsa->fetch_by_dbID($fset_id);
        #print Dumper $fset->name;

		#count? will never be more than one for a given fset
        $regulatory_features{$focus_start}{$focus_end}{$fset->name}++;
		$regulatory_features{$focus_start}{$focus_end}{'score'} = $focus_score;
    }

    return \%regulatory_features;
}


###########################
#sub write_intersection_features()
###########################
#{
#    my ($overlap_features) = @_;
#    my (%i_fsets, %i_ftypes, @i_features);
#
#
#    my $analysis = Bio::EnsEMBL::Analysis->new
#        (
#         -logic_name      => 'Co-occurrence',
#         -db              => 'NULL',
#         -db_version      => 'NULL',
#         -db_file         => 'NULL',
#         -program         => 'NULL',
#         -program_version => 'NULL',
#         -program_file    => 'NULL',
#         -gff_source      => 'NULL',
#         -gff_feature     => 'NULL',
#         -module          => 'NULL',
#         -module_version  => 'NULL',
#         -parameters      => 'NULL',
#         -created         => 'NULL',
#         -description     => 'Co-occurrence of FeatureTypes',
#         -display_label   => 'Co-occurrence',
#         -displayable     => 1
#         );
#    $analysis = $aa->fetch_by_dbID($aa->store($analysis));
#
#    foreach my $ft (@$overlap_features) {
#        
#        my ($start, $end, $fset_id, 
#            $focus_start, $focus_end, $focus_fset_id) = @$ft;
#
#        next if ($focus_fset_id == $fset_id);
#
#        my ($i_start,$i_end) = (sort {$a<=>$b} ($focus_start,
#                                                $focus_end,
#                                                $start,
#                                                $end))[1,2];
#
#        my $i_ftype_name = 'Intersection';
#        my $i_fset_name = $i_ftype_name.'_'.
#            join(':',
#                 sort($fsa->fetch_by_dbID($focus_fset_id)->name(),
#                      $fsa->fetch_by_dbID($fset_id)->name()));
#
#        my $i_feature = Bio::EnsEMBL::Funcgen::PredictedFeature->new
#            (
#             -slice  => $slice,
#             -start  => $i_start,
#             -end    => $i_end,
#             -strand => 0,
#             -feature_set => &get_FeatureSet($i_fset_name, \%i_fsets, 
#                                             $i_ftype_name,\%i_ftypes, 
#                                             $analysis)	
#             );
#
#        push(@i_features, $i_feature);
#    }
#    
#    $afa->store(@i_features);
#}


sub get_FeatureSet{
    my ($fset_name, $ftype_name, $analysis) = @_;

	#whoa, this is setting the fset in the hash passed
	#careful we don't add the overlap/intersect fset to the focus_sets
	
	#The caching functionality is a bit redundant at the mo, as we're
	#only concerned with merged fsets i.e. we're just setting the fset in the main script
	#rather than accessing dynamically dependent on which multiplex fset we're dealing with
	#to implement this we really need to simplify access, remove some of the args?

	#this was originally written for handling a cache of multiplex fsets
	#hence the complexity, simplify, or keep for use with multiple focus sets

	#do we need to pass the $fsets hash? - removed from args
	#this is only necessary if we're populating more than one hash
	#removed $ftypes, as this can be global

	#this is really more suited to the old multiplex fsets
	#needs rewrite if we're going to maintain merged binary vector fsets


    if (! exists $cooc_fsets{$fset_name}) {
        
        if (defined ($cooc_fsets{$fset_name} = $fsa->fetch_by_name($fset_name))) {
            
            if ($clobber) {
			  my $cs_id = $db->get_FGCoordSystemAdaptor->fetch_by_name('chromosome')->dbID();

			  print "Deleting AnnotatedFeatures for $fset_name on chromosome ".
				  $slice->seq_region_name()."\n";

			  

                my $sql = 'DELETE from annotated_feature where feature_set_id='.
                    $cooc_fsets{$fset_name}->dbID().' and seq_region_id='.
					$slice->get_seq_region_id().
					' and coord_system_id='.$cs_id;
                
                $db->dbc->do($sql) 
                    or throw('Failed to roll back $fset_name PredictedFeatures on chromosome '.
							 $slice->seq_region_name());
            } 
			else {
			  throw("There is a pre-existing FeatureSet '$fset_name'.\n".
					'You must specify clobber (option -c) is you want to delete'.
					' and overwrite all pre-existing PredictedFeatures');
            }

        } else {
   
            warn("Feature set '$fset_name' not defined in database!\n".
                 " Generating new FeatureSet.");
            
            if (! exists $ftypes{$ftype_name}) {
                if (! defined ($ftypes{$ftype_name} = $fta->fetch_by_name($ftype_name))) {
                    
                    $ftypes{$ftype_name} = Bio::EnsEMBL::Funcgen::FeatureType->new
                        (
                         -name        => $ftype_name,
                         -description => "$ftype_name of $fset_name",
                         );
                    
                    ($ftypes{$ftype_name}) = @{$fta->store($ftypes{$ftype_name})};
                }
            }
            
            $cooc_fsets{$fset_name} = Bio::EnsEMBL::Funcgen::FeatureSet->new
			  (
			   -analysis     => $analysis,
			   -feature_type => $ftypes{$ftype_name},
			   -name         => $fset_name,
			   -feature_class=> 'annotated',
			  );
            
            ($cooc_fsets{$fset_name}) = @{$fsa->store($cooc_fsets{$fset_name})};
                
#            #generate data_set here too???
#            my $dset = Bio::EnsEMBL::Funcgen::DataSet->new
#                (
#                 -feature_set => $union_fsets{$set_name},
#                 -name        => $set_name,
#                 );
#            
#            $dsa->store($dset);

        }
            
    }

    return $cooc_fsets{$fset_name};

}




1;
