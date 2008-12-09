#!/software/bin/perl
###!/usr/bin/perl

=head1 NAME

build_regulatory_features.pl -- builds features for the "Ensembl 
Regulatory Build", the moral equivalent of the gene build

=head1 SYNOPSIS

run_build_regulatory_features.pl -host host -user user -pass password 
    -outdir output_directory -focus feature_setA,feature_setB  
    -target feature_setC,feature_setD -seq_region_name string


Options:

  Mandatory
    -host|h            Host for eFG DB
    -user|u            User for eFG DB
    -pass|p            Password for eFG DB
    -dbname|d          Name of eFG DB
    -data_version|v    Version of data in eFG DB (e.g. 51_36m)
    -outdir|o          Name of outputut directory
    
    -focus|f           Focus features
    -attrib|a          Attribute features
    

  Optional
    -port              Port for eFG DB, default is 3306
    -species           Latin species name e.g. homo_sapiens

    -seq_region_name|s 
    -gene_signature 

    -write_features|w
    -clobber
    -dump_annotated_features
    -dump_regulatory_features

    -stats

    -help              Does nothing yet
    -man               Does nothing yet


=head1 DESCRIPTION

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! STRATEGY CHANGED; DOCUMENTATION NEEDS TO BE UPDATED !!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

This script is the core to compute the regulatory features for the 
"Ensembl Regulatory Build". A regulatory feature consists of

 a) all features that directly overlap with a focus feature, and
 b) features that are contained in focus-feature-overlaping features

The following figure gives examples.

      |------|  |-------- F1 --------|   |----|

      |---------------------------------------|

  |--X--|  |------| |------|   |---------|  |--X---|

      |============= RegFeature ==============|

                                            |- F2 -|

      |=============== RegFeature =================|

[more documentation to be added]

=head1 LICENCE

This code is distributed under an Apache style licence. Please see
http://www.ensembl.org/info/about/code_licence.html for details.

=head1 AUTHOR

Stefan Graf, Ensembl Functional Genomics

=head1 CONTACT

Please post comments/questions to the Ensembl development list
<ensembl-dev@ebi.ac.uk>

=cut


#To do

# 1 There are warns and STDERR prints, but we don't specify an -e or a -o file when bsubing
#   This has the knock on effect of not knowing what has failed, as we're arbitrarily using the
#   job index to select a chromosome from an array.  As this does not necessarily map to the job name it's 
#   Could also separate the STDERR and STDOUT, not-essential.

# 2 Need to run using non-ref toplevel, currently just using toplevel

# 3 Fix LSB_JOBINDEX usage. We we're always getting one failure as we're subtracting 1 from the 
#   index which for LSB_JOBINDEX = 0 resulting in -1, which is not a valid array index. DONE.  
#   Still need to improve this so we're not hardcoding how many slices we're running with (now 
#   getting two failures from non-existant slices).  We need to convert the perl script to perl 
#   so we can grab the toplevel non-ref slices before bsubing.  Can we use LSB_JOBNAME instead 
#   and use the seq_region names directly? 
#   SG -- The run script is now written in perl and does facilitate running the build on reference
#   toplevel slices. For now we exclude haplotypes as they are no captured in the underlying data. 
#   However once we are using our own MAQ mappings strategy mightneed to change to reflect this (also 2).

# SG: Need to double-check that feature_set and data_set are stored and have status "displayable" as well as 
#     having supporting sets correctly associated with the data_set. In v52 there was still a problem.

use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;
$|=1;

my ($pass,$port,$host,$user,$dbname,$species,$help,$man,
    $data_version,$outdir,$do_intersect,$write_features,
    $dump_annotated_features,$dump_regulatory_features,
    $seq_region_name,$clobber,$focus_max_length,
    $focus_extend,$focus,$attrib,$dump,$gene_signature,$stats,
    $debug,$debug_start,$debug_end);

$host = $ENV{EFG_HOST};
$port = $ENV{EFG_PORT};
$user = $ENV{EFG_WRITE_USER};
$dbname = $ENV{EFG_DBNAME};
$species = $ENV{SPECIES};
$data_version = $ENV{DATA_VERSION};

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
            "dump_annotated_features"  => \$dump_annotated_features,
            "dump_regulatory_features"  => \$dump_regulatory_features,
            "seq_region_name|s=s" => \$seq_region_name,
            "clobber" => \$clobber,
            "focus|f=s" => \$focus,
            "focus_max_length=i" => \$focus_max_length,
            "focus_extend=i" => \$focus_extend,
            "attrib|a=s" => \$attrib,
            "dump" => \$dump,
            "gene_signature" => \$gene_signature,
            "stats" => \$stats,
            "debug" => \$debug,
            "debug_start=i" => \$debug_start,
            "debug_end=i" => \$debug_end
            );

#Can we catch unknown options here to avoid missing incorrect params
#Something like:
#my $opts_out = Getoptions();
#Then check $opts_out.

### defaults ###
$port = 3306 if !$port;
$species = 'homo_sapiens' if !$species;#NJ make this mandatory?

### check options ###

throw("Must specify mandatory database hostname (-host).\n") if ! defined $host;
throw("Must specify mandatory database username. (-user)\n") if ! defined $user;
throw("Must specify mandatory database password (-pass).\n") if ! defined $pass;
throw("Must specify mandatory database name (-dbname).\n") if ! defined $dbname;
throw("Must specify mandatory database data version, like 47_36i (-data_version).\n") 
     if !$data_version;

throw("Must specify mandatory output directory (-outdir).\n") 
     if !$outdir;

throw("Must specify mandatory focus sets (-focus).\n") if ! defined $focus;
throw("Must specify mandatory attribute sets (-attrib).\n") if ! defined $attrib;

$focus_max_length = 2000 if (! defined $focus_max_length);
$focus_extend = 2000 if (! defined $focus_extend);

#throw("No output directory specified! Use -o option.") if (!$outdir);
if (defined $outdir && ! -d $outdir) {
    system("mkdir -p $outdir");
}
$outdir =~ s/\/$//;


$| = 1;


use Bio::EnsEMBL::Utils::Exception qw(verbose throw warning info);
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Funcgen::Utils::EFGUtils qw(open_file);
use Bio::EnsEMBL::Funcgen::RegulatoryFeature;
#use Bio::EnsEMBL::Funcgen::Utils::RegulatoryBuild qw(is_overlap);
#NJ what was this being used for and can we use the RangeRegistry?


# use ensembldb as we may want to use an old version
#NJ Default should be staging, but add params for overriding

my $cdb = Bio::EnsEMBL::DBSQL::DBAdaptor->new
    (
     -host => 'ens-staging',
     -port => 3306,
     -user => 'ensro',
     #-host => 'ensembldb.ensembl.org',
     #-user => 'anonymous',
     #-port => 3306,
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
     -dnadb  => $cdb,
     );


my $fsa = $db->get_FeatureSetAdaptor();
my $dsa = $db->get_DataSetAdaptor();
my $fta = $db->get_FeatureTypeAdaptor();
my $afa = $db->get_AnnotatedFeatureAdaptor();
my $rfa = $db->get_RegulatoryFeatureAdaptor();
my $sa = $db->get_SliceAdaptor();
my $aa = $db->get_AnalysisAdaptor();
my $ga = $db->dnadb->get_GeneAdaptor();
#my $ta = $cdb->get_TranscriptAdaptor();


# parse focus and attribute sets and check that they exist
my (%focus_fsets, %attrib_fsets);
map { my $fset = $fsa->fetch_by_name($_);
      throw("Focus set $_ does not exist in the DB") 
          if (! defined $fset); 
      $focus_fsets{$fset->dbID} = $fset; 
  } split(',', $focus);
#print Dumper %focus_fsets;

map { 
    my $fset = $fsa->fetch_by_name($_);
    throw("Attribute set $_ does not exist in the DB") 
        if (! defined $fset); 
    $attrib_fsets{$fset->dbID()} = $fset; 
} split(',', $attrib);
#print Dumper %attrib_fsets;

# make sure that attribute sets also contain focus sets (Do we really need this?)
map { $attrib_fsets{$_} = $focus_fsets{$_} } keys %focus_fsets;

### build regulatory features

### ChipSeq stuff deactivated; filtering after peak calling with SWEmbl not ncessary any more
#my %ChIPseq_cutoff = (
#                      ### cutoff T/O <= 2
#                      'CD4_CTCF'=>        5,
#                      'CD4_H3K27me3'=>    8,
#                      'CD4_H3K36me3'=>    4,
#                      'CD4_H3K4me3'=>     6,
#                      'CD4_H3K79me3'=>   22,
#                      'CD4_H3K9me3'=>     7,
#                      'CD4_H4K20me3'=>   17,
#                      ### cutoff ~ <= 25000
#                      ### see /lustre/work1/ensembl/graef/efg/input/SOLEXA/LMI/data/*.clstr.cutoff_25000.dat
#                      ### original data now moved to /lustre/work1/ensembl/graef/efg/input/SOLEXA/LMI_methylation
#                      'CD4_H2AZ'=>       15,
#                      'CD4_H2BK5me1'=>   16,
#                      'CD4_H3K27me1'=>    6,
#                      'CD4_H3K27me2'=>    5,
#                      'CD4_H3K36me1'=>    4,
#                      'CD4_H3K4me1'=>    31,
#                      'CD4_H3K4me2'=>    10,
#                      'CD4_H3K79me1'=>    5,
#                      'CD4_H3K79me2'=>    4,
#                      'CD4_H3K9me1'=>    12,
#                      'CD4_H3K9me2'=>     5,
#                      'CD4_H3R2me1'=>     5,
#                      'CD4_H3R2me2'=>     5,
#                      'CD4_H4K20me1'=>   30,
#                      'CD4_H4R3me2'=>     4,
#                      'CD4_PolII'=>       8
#                      );

# retrieve sequence to be analyzed 
my $slice;

if ($seq_region_name) {
    eval {
        $slice = $sa->fetch_by_name($seq_region_name);
    };
    
    if ($@) {

        warn("Couldn't retrieve slice '$seq_region_name'");

        my $slices = $sa->fetch_all('toplevel');
        #print Dumper @slices;

        print STDERR ("Available toplevel slices are:\n");

        map { print STDERR ("\t", join (':', $_->coord_system_name, '', $_->seq_region_name ), "\n")} 
        sort {$a->seq_region_name cmp $b->seq_region_name} @$slices;

        throw("Select a toplevel slice from the above list");

    }

} elsif (defined $ENV{LSB_JOBINDEX}) {
  
    warn "Performing whole genome analysis on toplevel slices using the farm (LSB_JOBINDEX: ".$ENV{LSB_JOBINDEX}.").\n";
    
    my @slices = ();
    my $sa = $db->get_SliceAdaptor();
    #foreach my $s (sort {$a->name cmp $b->name} @{$sa->fetch_all('toplevel', undef, 1)}) {
    foreach my $s (sort {$a->name cmp $b->name} @{$sa->fetch_all('toplevel')}) {
        
        next if ($s->seq_region_name =~ m/^MT/);
        push @slices, $s;
        
    }
    #print Dumper @slices;

    throw ("LSB_JOBINDEX is too large. Use a a value between 1 and ".
           scalar @slices) if ($ENV{LSB_JOBINDEX} > scalar @slices);

    $slice=$slices[$ENV{LSB_JOBINDEX}-1];
    print Dumper ($ENV{LSB_JOBINDEX}, $slice->name);

	warn "Slice is:\t".$slice->name.".\n";
   
} else {

    throw("Must either specify mandatory chromosome name (-seq_name) or use the \n".
          "wrapper script 'run_build_regulatory_features.pl' to perform whole genome \n".
          "analysis on all toplevel slices using the farm (via LSF environment \n".
          "variable LSB_JOBINDEX).\n");

}

# get core and fg seq_region_id for slice
my ($core_sr_id, $fg_sr_id);
$core_sr_id = $slice->get_seq_region_id();

# need to build cache first
$afa->build_seq_region_cache();
$fg_sr_id = $afa->get_seq_region_id_by_Slice($slice);

throw("No eFG seq_region id defined for ".$slice->name." Almost certain there is no data".
      "available on this slice. You might want to run_update_DB_for_release.pl") 
    if ! defined $fg_sr_id;

#should this be printing to OUT or STDERR?
#or remove as we're printing this later?
print
    '# Focus set(s): ', join(", ", map {$_->name.' ('.$_->dbID.')'}
                             sort {$a->name cmp $b->name} values %focus_fsets), "\n",
    '# Attribute set(s): ', join(", ", map {$_->name.' ('.$_->dbID.')'} 
                              sort {$a->name cmp $b->name}  values %attrib_fsets), "\n",
    '# Species: ', $species, "\n",
    '# Seq region: ', join(" ",$slice->display_id(),$slice->get_seq_region_id), "\n",
    '#   core seq_region_id '.$core_sr_id." => efg seq_region_id ".$fg_sr_id."\n";


### Check whether analysis is already stored 

my $analysis = Bio::EnsEMBL::Analysis->new
    (
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
     -description     => 'Union of focus features, features overlapping focus features,'.
     ' and features that are contained within those',
     -display_label   => 'RegulatoryRegion',
     -displayable     => 1,
     );

my $logic_name = $analysis->logic_name();
my $ana = $aa->fetch_by_logic_name($logic_name);

if ( ! defined $ana ) { # NEW
    
    warn("Need to store new analysis with logic name $logic_name.");
    $aa->store($analysis);
    
} elsif ( $ana->compare($analysis) ) { # EXISTS, but with different options

    ### analysis compare
    # returns  1 if this analysis is special case of given analysis
    # returns  0 if they are equal
    # returns -1 if they are completely different
    
    throw('Analysis with logic name \''.$logic_name.'\' already exists, but '.
          "has different options! Use different logic_name for you analysis '$logic_name'!");

    #$self->efg_analysis->dbID($analysis->dbID);
    #$self->efg_analysis->adaptor($self->efgdb->get_AnalysisAdaptor);
    #$aa->update($self->efg_analysis);
 
} else { # EXISTS

    warn('Analysis with logic name \''.$logic_name.'\' already '.
         'exists.');
    
}

$analysis = $aa->fetch_by_logic_name($logic_name);

my $rfset = &get_regulatory_FeatureSet($analysis);
#print Dumper $rfset;


# Read from file and process sequentially in sorted by start, end order. Each 
# new feature is checked if it overlaps with the preceeding already seen 
# features. If yes we just carry on with the next one. Otherwise

#$dbname =~ s/sg_//;
my $af_file = $outdir.'/'.$dbname.'.annotated_features.'.$slice->seq_region_name.'.dat';

if (! -e $af_file  || $dump_annotated_features) {

    warn("File containing dumped annotated features doesn't exist. Need ".
         "to dump annotated features on slice '".$slice->name."' first!");
    
    #eval {
        &dump_annotated_features();
    #};

    #throw ("Couldn't dump annotated features on slice '".$slice->name."'.") if ($@);
    
} else {

    warn("File '$af_file' containing dumped annotated features exists.")

}

my $fh = open_file($af_file);

my (@rf,@af);
my ($af_id, $sr_id, $sr_name, $start, $end, $strand, $score, $fset_id);

### variables for statistics
my (%feature_count, %seen_af, %removed_af);

while (<$fh>) {

    next if (/^\#/);
    chomp;
    
    ($af_id, $sr_id, $sr_name, $start, $end, $strand, $score, $fset_id) = split (/\s+/, $_);
    #print Dumper ($af_id, $sr_id, $start, $end, $strand, $score, $fset_id);
    #print $af_id, "\n";

    if ($debug && $seq_region_name && $debug_start && $debug_end) {

        next if ($start < $debug_start);
        last if ($start > $debug_end);
        
    }
    
    # Quick hack for 2nd/3rd version of reg. build. Need to disregard ChIPseq 
    # features below a certain threshold defined hardcoded in ChIPseq_cutoff hash.
    #next if (exists $ChIPseq_cutoff{$attrib_fsets{$fset_id}->name()} 
    #         && $score < $ChIPseq_cutoff{$attrib_fsets{$fset_id}->name()});
    
    print $_, "\t", $attrib_fsets{$fset_id}->name, "\n" if ($debug);
    my  $length = $end-$start+1;
    
    
    # some stats
    $feature_count{$fset_id}++;
    
    if (exists $focus_fsets{$fset_id}) {
        
        # focus feature
        
        # focus features that are longer than a given threshold are only 
        # considered as attribute features
        if ($length > $focus_max_length) {
            print "focus feature ($af_id) longer than $focus_max_length; ",
            "added to attribute feature list\n" if ($debug);
            &add_feature();
            $removed_af{$fset_id}{$af_id} = 1;
            # need to check overlaps with current focus_features
            next;
        }
        
        # current feature is focus feature
        print "focus feature ($af_id)\n" if ($debug);
        
        
        # focus feature overlaps w/ rf core (0° of separation)
        if ( @rf && $start <= $rf[$#rf]{focus_end}) {
            
            &update_focus();
            
        # open new regulatory feature
        } else {
            
            &add_focus();
            &update_attributes();
            
        }
        
        &add_feature();
        
    } else {
        
        # attribute feature
        
        if ( @rf && $start <= $rf[$#rf]{focus_end} ) {
            
            print "overlap w/ focus feature; add (".$af_id.") to reg. feature\n" 
                if ($debug);
            
            # add annot. feature id to reg. feature
            $rf[$#rf]{annotated}{$af_id} = undef;
            $rf[$#rf]{fsets}{$fset_id}++;
            
            $seen_af{$fset_id}{$af_id} = 1;
            
            # update end of regulatory feature
            $rf[$#rf]{attribute_end} = $end
                if ($end > $rf[$#rf]{attribute_end});
            
        } elsif (@rf && $end <= $rf[$#rf]{attribute_end} &&
                 $start <= $rf[$#rf]{focus_end}+$focus_extend) {
            
            print "contained within reg. feature; add (".$af_id.") to reg. feature\n" 
                if ($debug);
            
            # add annot. feature id to reg. feature
            $rf[$#rf]{annotated}{$af_id} = undef;
            $rf[$#rf]{fsets}{$fset_id}++;
            
            $seen_af{$fset_id}{$af_id} = 1;
            
        }
        
        print "add to attribute feature list ($af_id)\n" if ($debug);
        &add_feature();
        
    }
    
}

# build binary string

map {
    
    $_->{binstring} = &build_binstring($_);
    
} @rf;

print "\n", Dumper @rf if ($debug);


if ($dump_regulatory_features) {
    
    my $outfile  = $outdir.'/'.$dbname.'.regulatory_features_'.$fg_sr_id.'.dat';
    my $out = open_file($outfile, ">");
    
    map {
        printf $out "%d\t%s\t%d\t%d\t%d\t%d\t%s\t%s\n", 
        $fg_sr_id, $slice->seq_region_name, 
        $_->{focus_start}, $_->{focus_end},
        $_->{attribute_start}, $_->{attribute_end}, 
        $_->{binstring}, join(",", sort {$a<=>$b} keys %{$_->{annotated}});
    } @rf;
    
    printf "# Regulatory features written to ".$outfile."\n";
    
}

if ($write_features) {
    
    my @f = ();
    
    map {
        push (@f, &get_regulatory_feature($_));
    } @rf;
    
    #print Dumper @f;
    #warn("STORING DEACTIVATED");
    $rfa->store(@f);
    
}


my $feature_count;
map {$feature_count+=$_} values %feature_count;
printf "# Number of read features: %10d\n", $feature_count;
printf "# Number of reg. features: %10d\n", scalar(@rf);

if ($stats) {
    my (%rf_count);
    map { 
        foreach my $k (keys %{$_->{fsets}}){
            #print join(" ", $k, $_->{fsets}->{$k}), "\n";
            $rf_count{$k} += $_->{fsets}->{$k};
        }
    } @rf;
    
    printf "# %-34s\t%8s\t%8s\t%8s\t%8s\n",
    'Number of feature sets', 'total', 'included', 'distinct', 'removed';
    map {
        printf "# %29s (%d)\t%8d\t%8d\t%8d\t%8d\n", $_->name, $_->dbID, 
        $feature_count{$_->dbID}||0, $rf_count{$_->dbID}||0,
        scalar(keys %{$seen_af{$_->dbID}})||0,
        scalar(keys %{$removed_af{$_->dbID}})||0
        } sort {$a->name cmp $b->name} values %attrib_fsets;
}

###############################################################################
# dump annotated features for given seq_region to file
sub dump_annotated_features () {

    my @fset_ids = keys %attrib_fsets;
    #print Dumper @fset_ids; 

    print STDERR "# Dumping annotated features from database to file.\n";
    print STDERR "# This will delete existing file dumps of that data version!\n";
    
    print STDERR "# Output goes to ", $outdir, "\n";
    
    my $sql = "select annotated_feature_id, af.seq_region_id,".
        "sr.name, seq_region_start, seq_region_end,seq_region_strand,".
        "score,feature_set_id from annotated_feature af, seq_region sr ".
        "where sr.seq_region_id=af.seq_region_id ".
        "and schema_build='$data_version' ".
        "and sr.name='".$slice->seq_region_name."' ".
        "and af.feature_set_id in (".join(',', @fset_ids).")";
    
    my $command = "echo \"$sql\" ".
        " | mysql -quick -N -h".$host." -P".$port." -u".$user." -p".$pass." ".$dbname.
        #" | gawk '{if (\$7==".join("||\$7==", @fset_ids).") print }'".
        " | sort -k3,3n -k4,4n -k5,5n".
        " | gawk '{ print >> \"".$outdir."/".$dbname.".annotated_features.".$slice->seq_region_name.".dat\" }'";
    
    warn("# Execute: ".$command."\n");# if ($debug);

    # need to remove existing dump files, since we append to the file
    #system("rm -f $outdir/$dbname.annotated_feature_*.dat") &&
    #        throw ("Can't remove files");
    
    system($command) &&
        throw ("Can't dump data to file in $outdir");
    
}

###############################################################################

sub add_focus ()
{

    print "add focus\n" if ($debug);

    push @rf, 
    {
        'focus_start' => $start,
        'focus_end' => $end,
        'attribute_start' => $start,
        'attribute_end' => $end,
        'annotated' => { 
            $af_id => undef
            },
                'fsets' => {
                    $fset_id => 1
                    }
    };

    $seen_af{$fset_id}{$af_id} = 1;
    
    #print Dumper @rf;
    
}

sub update_focus ()
{

    print "update focus_feature\n" if ($debug);

    $rf[$#rf]{focus_end} = $end
        if ($end > $rf[$#rf]{focus_end});
    $rf[$#rf]{attribute_end}=$rf[$#rf]{focus_end};

    $rf[$#rf]{fsets}{$fset_id}++;
    $rf[$#rf]{annotated}{$af_id}=undef;

    $seen_af{$fset_id}{$af_id} = 1;

    #print Dumper @rf if ($debug);

}

sub update_attributes ()
{

    print "update attributes\n" if ($debug);
    
    # look upstream for features overlaping focus feature
    for (my $i=0; $i<=$#af; $i++) {
        
        if ($af[$i]{end} >= $start) {

            # FIRST update attribute start and end of regulatory feature...
            $rf[$#rf]{attribute_start} = $af[$i]{start}
            if ($af[$i]{start} < $rf[$#rf]{attribute_start});
            $rf[$#rf]{attribute_end} = $af[$i]{end}
            if ($af[$i]{end} > $rf[$#rf]{attribute_end});
            
            # ... then remove features that are out of scope
            splice(@af, 0, $i);
            #print Dumper @af;

            # and add the remaining attribute features if applicable
            map {
                
                if ($_->{end} >= $start ||
                    $_->{start} >= $rf[$#rf]{focus_start}-$focus_extend) {

                    $rf[$#rf]{annotated}{$_->{af_id}} = undef;
                    $rf[$#rf]{fsets}{$_->{fset_id}}++;

                    $seen_af{$fset_id}{$af_id} = 1;

                }

            } @af;
            
            last;

        }

    }

}

sub add_feature ()
{

    push(@af, 
         {
             af_id => $af_id,
             start => $start,
             end => $end,
             strand => $strand,
             score => $score,
             fset_id => $fset_id
             }
         );

}

sub build_binstring()
{

    my ($rf) = @_;

    my $binstring = '';
    foreach (sort {$a->name cmp $b->name} values %attrib_fsets) {
        $binstring .= 
            (exists $rf->{fsets}->{$_->dbID})? 1 : 0;	
    }

    $binstring .= &get_gene_signature($rf) if ($gene_signature);
    
    return $binstring;

}

sub get_gene_signature()
{

    my ($rf) = @_;


    my $tss = 0;
    my $tts = 0;

    # get feature slice +/- 2.5kB
    my $expand = 2500;
    
    my ($feature_slice, $expanded_slice);
    $feature_slice = $sa->fetch_by_seq_region_id($core_sr_id,
                                                 $rf->{focus_start},
                                                 $rf->{focus_end});
    $expanded_slice = $feature_slice->expand($expand, $expand);

    foreach my $g (@{$ga->fetch_all_by_Slice($expanded_slice)}) {
        
        $tss = 1 if (($g->strand == 1  && $g->start > 0) ||
                     ($g->strand == -1 && $g->end <= $expanded_slice->length));
        $tts = 1 if (($g->strand == 1  && $g->end <= $expanded_slice->length) ||
                     ($g->strand == -1 && $g->start > 0));

        last if ($tss && $tts);

    }
    
    return $tss.$tts;

}

sub get_regulatory_feature{

    my ($rf) = @_;

    #print Dumper $rf;
    #print Dumper $rfset->feature_type;
    #print $rfset->feature_type;
    if ($slice->start != 1 || $slice->strand != 1) {

        warn("**** SLICE doesn't start at position 1; resetting start of slice to 1 ****");
        $slice->{'start'} = 1;
        #print "s_start:  ", $slice->start(), "\n";
        ###$regulatory_feature->transfer($s)
        
    }

    my $regulatory_feature = Bio::EnsEMBL::Funcgen::RegulatoryFeature->new
        (
         -slice            => $slice,
         -start            => $rf->{focus_start},
         -end              => $rf->{focus_end},
         -strand           => 0,
         -display_label    => $rf->{binstring} || '',
         -feature_set      => $rfset,
         -feature_type     => $rfset->feature_type,
         -_attribute_cache => {'annotated' => $rf->{annotated}},
         );


    return $regulatory_feature;



}




sub get_regulatory_FeatureSet{
    
    my $rfset = $fsa->fetch_by_name('RegulatoryFeatures');

    if ($rfset) {

        if ($write_features && $clobber) {
            
            my $sql = 
                'DELETE ra FROM regulatory_attribute ra, regulatory_feature rf'.
                ' WHERE ra.regulatory_feature_id=rf.regulatory_feature_id'.
                ' AND feature_set_id='.$rfset->dbID().' AND seq_region_id='.$fg_sr_id;
            $db->dbc->do($sql) 
                or throw('Failed to roll back regulatory_features for feature_set_id'.
                         $rfset->dbID());

            $sql = 
                ' DELETE FROM regulatory_feature WHERE feature_set_id='.$rfset->dbID().
                ' AND seq_region_id='.$fg_sr_id;
            $db->dbc->do($sql) 
                or throw('Failed to roll back regulatory_features for feature_set_id'.
                         $rfset->dbID());

        } elsif ($write_features) {
            throw("Their is a pre-existing FeatureSet with the name 'RegulatoryFeatures'\n".
                  'You must specify clobber if you want to delete and overwrite all'.
                  ' pre-existing RegulatoryFeatures');
        }
    } else {					#generate new fset

        my $ftype = $fta->fetch_by_name('RegulatoryFeature');
        
        if (! $ftype) {
            
            $ftype = Bio::EnsEMBL::Funcgen::FeatureType->new
                (
                 -name        => 'RegulatoryFeature',
                 -description => 'RegulatoryFeature',
                 );
            
            $ftype = @{$fta->store($ftype)} if ($write_features);

        }

        $rfset = Bio::EnsEMBL::Funcgen::FeatureSet->new
            (
             -analysis     => $analysis,
             -feature_type => $ftype,
             -name         => 'RegulatoryFeatures',
             -type         => 'regulatory'
             );

        $rfset->add_status('DISPLAYABLE');
        $rfset = @{$fsa->store($rfset)} if ($write_features);

        #generate data_set here too

        my $dset = $dsa->fetch_by_name('RegulatoryFeatures');

        if (! $dset) {

            my $dset = Bio::EnsEMBL::Funcgen::DataSet->new
                (
                 -feature_set => $rfset,
                 -name        => 'RegulatoryFeatures',
                 -supporting_set_type => 'feature'
                 );
            
            $dset->add_status('DISPLAYABLE');
            $dsa->store($dset) if ($write_features);
        }
    }

    return $rfset;


}

1;
