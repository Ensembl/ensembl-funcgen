#!/software/bin/perl
###!/usr/bin/perl

=head1 NAME

build_regulatory_features.pl -- builds features for the "Ensembl 
Regulatory Build", the moral equivalent of the gene build

=head1 SYNOPSIS

build_regulatory_features.pl -host host -user user -pass password 
    -outdir output_directory -focus feature_setA,feature_setB  
    -target feature_setC,feature_setD -seq_name chr


Options:

Mandatory
    -pass|p    Password for eFG DB
    -host|h    Host for eFG DB      => \$host,
            "user|u=s"       => \$user,
            "dbname|d=s"     => \$dbname,



            "data_version|v=s" => \$data_version,
            "outdir|o=s"     => \$outdir,
            "do_intersect|i=s" => \$do_intersect,
            "write_features|w" => \$write_features,
            "dump_features"  => \$dump_features,
            "seq_name|s=s" => \$seq_name,
            "clobber" => \$clobber,
            "focus|f=s" => \$focus,
            "focus_max_length=i" => \$focus_max_length,
            "focus_extend=i" => \$focus_extend,
            "target|t=s" => \$target,
            "dump" => \$dump,
            "gene_signature" => \$gene_signature,
            "stats" => \$stats,
            "debug" => \$debug,
            "debug_start=i" => \$debug_start,
            "debug_end=i" => \$debug_end


Optional
    -port      Port for eFG DB, default is 3306
    -species   Latin species name e.g. homo_sapiens

    -help    Does nothing yet.   
    -man     Does nothing yet.    


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

# 3 Fix LSB_JOBINDEX usage. We we're always getting one failure as we're subtracting 1 from the index which for LSB_JOBINDEX = 0 resulting in -1, which is not a valid array index. DONE.  Still need to improve this so we're not hardcoding how many slices we're running with(now getting two failures from non-existant slices).  We need to convert the perl script to perl so we can grab the toplevel non-ref slices before bsubing.  Can we use LSB_JOBNAME instead and use the seq_region names directly?

use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;
$|=1;

my ($pass,$port,$host,$user,$dbname,$species,$help,$man,
    $data_version,$outdir,$do_intersect,$write_features,
    $dump_features,$seq_name,$clobber,$focus_max_length,
    $focus_extend,$focus,$target,$dump,$gene_signature,$stats,
    $debug,$debug_start,$debug_end);

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
            "dump_features"  => \$dump_features,
            "seq_name|s=s" => \$seq_name,
            "clobber" => \$clobber,
            "focus|f=s" => \$focus,
            "focus_max_length=i" => \$focus_max_length,
            "focus_extend=i" => \$focus_extend,
            "target|t=s" => \$target,
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
throw("Must specify mandatory target sets (-target).\n") if ! defined $target;

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
     #-host => 'ensembldb.ensembl.org',
     #-user => 'anonymous',
     #-port => 3306,
     -host => 'ens-staging',
     -port => 3306,
     -user => 'ensro',
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


# parse focus and target sets and check that they exist
my (%focus_fsets, %target_fsets);
map { my $fset = $fsa->fetch_by_name($_);
      throw("Focus set $_ does not exist in the DB") 
          if (! defined $fset); 
      $focus_fsets{$fset->dbID} = $fset; 
  } split(',', $focus);
#print Dumper %focus_fsets;

map { 
    my $fset = $fsa->fetch_by_name($_);
    throw("Target set $_ does not exist in the DB") 
        if (! defined $fset); 
    $target_fsets{$fset->dbID()} = $fset; 
} split(',', $target);
#print Dumper %target_fsets;

# make sure that target sets also contain focus sets (Do we really need this?)
map { $target_fsets{$_} = $focus_fsets{$_} } keys %focus_fsets;

# dump data to files and exit
if ($dump) {
    
    my @fset_ids = keys %target_fsets;
    #print Dumper @fset_ids; 

    print STDERR "# Dumping annotated features from database to file.\n";
    print STDERR "# This will delete existing file dumps of that data version!\n";

    print STDERR "# Output goes to ", $outdir, "\n";

    my $sql = "select annotated_feature_id, af.seq_region_id,".
        "seq_region_start, seq_region_end,seq_region_strand,".
        "score,feature_set_id from annotated_feature af, seq_region sr ".
        "where sr.seq_region_id=af.seq_region_id and schema_build='$data_version'";
    my $command = "echo \"$sql\" ".
        " | mysql -quick -h".$host." -P".$port." -u".$user." -p".$pass." ".$dbname.
        " | gawk '{if (\$7==".join("||\$7==", @fset_ids).") print }'".
        " | sort -k2,2n -k3,3n -k4,4n".
        " | gawk '{ print >> \"".$outdir."/".$dbname.".annotated_feature_\" \$2 \".dat\" }'";

    print STDERR "# Execute: ".$command."\n" if ($debug);
    
    # need to remove existing dump files, since we append to the file
    system("rm -f $outdir/$dbname.annotated_feature_*.dat") &&
        throw ("Can't remove files");

    system($command) &&
        throw ("Can't dump data to file in $outdir");

    exit;

}


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

if ($seq_name) {
    $slice = $sa->fetch_by_region('chromosome', $seq_name);
} elsif (defined $ENV{LSB_JOBINDEX}) {
  
    warn "Performing whole genome analysis on toplevel slices using the farm (LSB_JOBINDEX: ".$ENV{LSB_JOBINDEX}.").\n";
    
    my $toplevel = $sa->fetch_all('toplevel');#This needs to be non-reference too!
    my @chr = sort (map $_->seq_region_name, @{$toplevel});
    #print Dumper @chr;
    
    my @slices;
    foreach my $chr (@chr) {
        
        next if ($chr =~ m/^NT_/);
        
        push @slices, $sa->fetch_by_region('chromosome', $chr);
    }
    #print Dumper @slices;

    $slice=$slices[$ENV{LSB_JOBINDEX}];#-1];      
    #print Dumper ($ENV{LSB_JOBINDEX}, $slice->name);


	warn "Slice is:\t".$slice->name.").\n";


    
} else {

    throw("Must specify mandatory chromosome name (-seq_name) or set\n ".
          "LSF environment variable LSB_JOBINDEX to perform whole\n ".
          "genome analysis on toplevel slices using the farm.\n");

}

# get core and fg seq_region_id for slice
my ($core_sr_id, $fg_sr_id);
$core_sr_id = $slice->get_seq_region_id();
#need to build cache first
$afa->build_seq_region_cache();
$fg_sr_id = $afa->get_seq_region_id_by_Slice($slice);

#should this be printing to OUT or STDERR?
#or remove as we're printing this later?
print
    '# Focus set(s): ', join(", ", map {$_->name.' ('.$_->dbID.')'}
                             sort {$a->name cmp $b->name} values %focus_fsets), "\n",
    '# Target set(s): ', join(", ", map {$_->name.' ('.$_->dbID.')'} 
                              sort {$a->name cmp $b->name}  values %target_fsets), "\n",
    '# Species: ', $species, "\n",
    '# Chromosome: ', join(" ",$slice->display_id(),$slice->get_seq_region_id), "\n",
    '#   core seq_region_id '.$core_sr_id." => fg seq_region_id ".$fg_sr_id."\n";



my $analysis = Bio::EnsEMBL::Analysis->new(
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
                                           ### display_label is going to be changed to "RegulatoryBuild" or so
                                           -display_label   => 'RegulatoryRegion',
                                           -displayable     => 1,
                                           );
$analysis = $aa->fetch_by_dbID($aa->store($analysis));# if ($write_features) ;


my $rfset = &get_regulatory_FeatureSet($analysis);
#print Dumper $rfset;

# Read from file and process sequentially in sorted by start, end order. Each 
# new feature is checked if it overlaps with the preceeding already seen 
# features. If yes we just carry on with the next one. Otherwise

#$dbname =~ s/sg_//;
my $fh = open_file($outdir.'/'.$dbname.'.annotated_feature_'.$fg_sr_id.'.dat');

my (@rf,@af);
my ($af_id, $sr_id, $start, $end, $strand, $score, $fset_id);

### variables for statistics
my (%feature_count, %seen_af, %removed_af);

while (<$fh>) {

    next if (/^\#/);
    chomp;
    
    ($af_id, $sr_id, $start, $end, $strand, $score, $fset_id) = split (/\s+/, $_);
    #print Dumper ($af_id, $sr_id, $start, $end, $strand, $score, $fset_id);
    #print $af_id, "\n";

    if ($debug && $seq_name && $debug_start && $debug_end) {

        next if ($start < $debug_start);
        last if ($start > $debug_end);
        
    }
    
    # Quick hack for 2nd/3rd version of reg. build. Need to disregard ChIPseq 
    # features below a certain threshold defined hardcoded in ChIPseq_cutoff hash.
    #next if (exists $ChIPseq_cutoff{$target_fsets{$fset_id}->name()} 
    #         && $score < $ChIPseq_cutoff{$target_fsets{$fset_id}->name()});
    
    print $_, "\t", $target_fsets{$fset_id}->name, "\n" if ($debug);
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


if ($dump_features) {
    
    my $outfile  = $outdir.'/'.$dbname.'.regulatory_feature_'.$fg_sr_id.'.dat';
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
        } sort {$a->name cmp $b->name} values %target_fsets;
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
    foreach (sort {$a->name cmp $b->name} values %target_fsets) {
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
    #print $rfset->feature_type;

    return Bio::EnsEMBL::Funcgen::RegulatoryFeature->new
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
            throw("Their is a pre-existing FeatureSet with the name 'Regulatory_Features'\n".
                  'You must specify clobber is you want to delete and overwrite all'.
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
            
            $dsa->store($dset) if ($write_features);
        }
    }

    return $rfset;


}

1;
