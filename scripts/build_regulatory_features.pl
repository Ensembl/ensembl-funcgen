#!/software/bin/perl

=head1 NAME

build_regulatory_features.pl

=head1 SYNOPSIS

build_regulatory_features.pl -H dbhost -P 3306 -u user -d password 
    -o /tmp -f testA -t testB,testC

=head1 DESCRIPTION

Generates regulatory features based on co-occurance between a 
given focus feature set and target sets. It calculates the
intersection between each focus-target-pair and generates for 
each focus features a binary vector representing the co-occurence.

=head1 LICENCE

This code is distributed under an Apache style licence. Please see
http://www.ensembl.org/info/about/code_licence.html for details.

=head1 AUTHOR

Stefan Graf <graef@ebi.ac.uk>, Ensembl Functional Genomics

=head1 CONTACT

Please post comments/questions to the Ensembl development list
<ensembl-dev@ebi.ac.uk>

=cut

use strict;
use warnings;
use Data::Dumper;
use Getopt::Std;

$| = 1;

my %opts;
getopts('hH:P:u:p:d:f:t:s:cwDo:', \%opts);

my $dbhost = $opts{H};
my $dbport = $opts{P};
my $dbuser = $opts{u};
my $dbpass = $opts{p};
my $dbname = $opts{d};

if (! ($dbhost && $dbport && $dbuser && $dbname)) {
    throw("Must specify mandatory database parameters, like:\n".
          " -H host -P 3306 -u XXXX -p XXXX -d dbname");
}

my $outdir = $opts{o};
throw("No output directory specified! Use -o option.") if (!$outdir);
if (! -d $outdir) {
    system("mkdir -p $outdir");
}

use Bio::EnsEMBL::Utils::Exception qw(verbose throw warning info);
use Bio::EnsEMBL::Analysis::Tools::Logger qw(logger_verbosity logger_info);

my $utils_verbosity = 'WARNING';
my $logger_verbosity = 'OFF';
verbose($utils_verbosity);
logger_verbosity($logger_verbosity);

# databases and adaptors
use Bio::EnsEMBL::Registry;
Bio::EnsEMBL::Registry->load_registry_from_db
    (
     -host => 'ens-livemirror',
     -user => 'ensro',
     #-verbose => "1" 
     );
my $cdb = Bio::EnsEMBL::Registry->get_DBAdaptor('human', 'core');

my $db = Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor->new
    (
     -host   => $dbhost,
     -user   => $dbuser,
     -dbname => $dbname,
     -pass   => $dbpass,
     -port   => $dbport,
     -dnadb  => $cdb
     );

my $fsa = $db->get_FeatureSetAdaptor();
my $fta = $db->get_FeatureTypeAdaptor();
my $pfa = $db->get_PredictedFeatureAdaptor();
my $sa = $db->get_SliceAdaptor();
my $aa = $db->get_AnalysisAdaptor();

# checking options
throw("No focus set specified! Use -f option.") if (!$opts{f});
my %focus_fsets;
$focus_fsets{$opts{f}} = $fsa->fetch_by_name($opts{f});
#map { $focus_fsets{$_} = $fsa->fetch_by_name($_) } split(',', $opts{f});
throw("Focus set $opts{f} is not defined!")if (! defined $focus_fsets{$opts{f}});

throw("No target set(s) specified! Use -t option.") if (!$opts{t});
my %target_fsets;
map { $target_fsets{$_} = $fsa->fetch_by_name($_) } split(',', $opts{t});
map { $target_fsets{$_} = $fsa->fetch_by_name($_) } keys %focus_fsets;

my $slice;
if ($opts{s}) {

    $slice = $sa->fetch_by_region('chromosome', $opts{s});

} else {

    warn("Perfoming whole genome analysis on toplevel slices using the farm.");

    throw("LSF environment variable LSB_JOBINDEX not defined.") 
        if (! defined $ENV{LSB_JOBINDEX});

    my $toplevel = $sa->fetch_all('toplevel');
    my @chr = sort (map $_->seq_region_name, @{$toplevel});
    #print Dumper @chr;

    my @slices;
    foreach my $chr (@chr) {
        
        next if ($chr =~ m/^NT_/);
        
        push @slices, $sa->fetch_by_region('chromosome', $chr);
    }

    $slice=$slices[$ENV{LSB_JOBINDEX}-1];      
    print Dumper ($ENV{LSB_JOBINDEX}, $slice->name);

}


print '# Focus set: ', join(" ", keys %focus_fsets), "\n";
print '# Target set(s): ', join(" ", sort keys %target_fsets), "\n";

my (@starts, @ends, @fset_ids, 
    $focus_start, $focus_end, $focus_fset_id,
    @overlap_features);


foreach my $pf (@{$pfa->fetch_all_by_Slice($slice)}) {
    
    #print join(" ", $pf->start, $pf->end, $pf->feature_set->name()), "\n";
    next if(! exists $target_fsets{$pf->feature_set->name()});
    
    print join(" ", $pf->start, $pf->end,
               $pf->feature_set->dbID(),
               $pf->feature_set->name()), "\n" if ($opts{D});

    if (exists $focus_fsets{$pf->feature_set->name()}) {
        
        ### fix revisting / counting feature twice here (focus testC vs. target testA) ###
        if (defined $focus_end){ # && !exists $regulatory_features{$focus_start}{$focus_end}) {
            
            push @overlap_features, @{&get_overlap_features(\@starts, \@ends, \@fset_ids, 
                                                            $focus_start, $focus_end, $focus_fset_id)};
            #print Dumper $overlap_features;
            #&get_overlap_string(\%regulatory_features, \@overlap_features, $focus_start, $focus_end);
            
        }

        # focus feature set
        $focus_start = $pf->start;
        $focus_end = $pf->end;
        $focus_fset_id = $pf->feature_set->dbID();

    } elsif (defined $focus_end && $pf->start > $focus_end) {

        #warn ("focus: ".$focus_start." ".$focus_end." ".$focus_fset_id);

        # target feature sets
        push @overlap_features, @{&get_overlap_features(\@starts, \@ends, \@fset_ids, 
                                                        $focus_start, $focus_end, $focus_fset_id)};
        #print Dumper $overlap_features;
        #&get_overlap_string(\%regulatory_features, \@overlap_features, $focus_start, $focus_end);

    }
    
    push(@starts, $pf->start());
    push(@ends, $pf->end());
    push(@fset_ids, $pf->feature_set()->dbID());

}
if ($opts{D}) {
    print STDERR join("\t", @starts), "\n";
    print STDERR join("\t", @ends), "\n";
    print STDERR join("\t", @fset_ids), "\n";
}
push @overlap_features, @{&get_overlap_features(\@starts, \@ends, \@fset_ids, 
                                                $focus_start, $focus_end, $focus_fset_id)};
#print Dumper $overlap_features;
#&get_overlap_string(\%regulatory_features, \@overlap_features, $focus_start, $focus_end);

#print Dumper %regulatory_features;


# Dump/write output
my $regulatory_features = &get_overlap_strings(\@overlap_features);

my $i = 0;

my $outfile = $opts{f}.'_chr'.$slice->seq_region_name.'.overlap';

print "Output goes to $outdir/$outfile\n";
open(OUT, "> $outdir/$outfile")
    or throw("Can't open file $outdir/$outfile");

print OUT '# Focus set: ', join(" ", keys %focus_fsets), "\n";
print OUT '# Target set(s): ', join(" ", sort keys %target_fsets), "\n";
print OUT '# ', $slice->name(), "\n";
foreach my $s (sort {$a<=>$b} keys %{$regulatory_features}) {
    foreach my $e (keys %{$regulatory_features->{$s}}) {
        
        my ($overlap_string, $count_string);
        foreach (sort keys %{$regulatory_features->{$s}->{$e}}) {
            $overlap_string .= ($regulatory_features->{$s}->{$e}->{$_}>0)?1:0;
            $count_string .= $regulatory_features->{$s}->{$e}->{$_};
        }
        
        print OUT join("\t",
                       sprintf("ENSR_%s_%06d",$opts{f},++$i),
                       $slice->seq_region_name, $s, $e,
                       $overlap_string,
                       #$count_string
                       ), "\n";
        
    }

}

close OUT;

#if ($opts{w}) {
#    
#    &write_intersection_features(\@overlap_features);
#
#}


##########################
sub get_overlap_features()
##########################
{
    my ($starts, $ends, $fset_ids, 
        $focus_start, $focus_end, $focus_fset_id) = @_;
    print STDERR $focus_start.':'.$focus_end." (focus)\n" if ($opts{D});

    my @overlap_features = ();
    my ($start, $end, $fset_id);
    my (@s, @e, @f);
    while ($#starts>=0) {
        
        $start = shift @$starts;
        $end = shift @$ends;
        $fset_id = shift @$fset_ids;
        
        #don't compare the focus_set with itself
        #next if ($focus_fset_id == $fset_id);

        print STDERR $start.':'.$end."\n" if ($opts{D});
        
        if ($focus_start < $end && $focus_end > $start) {

            print STDERR "Found intersection: ".
                join(":", sort {$a<=>$b} ($focus_start, $focus_end, $start,$end)).
                " with ". $fsa->fetch_by_dbID($fset_id)->name()."\n" if ($opts{D});

            push @overlap_features, [ $start, $end, $fset_id, 
                                      $focus_start, $focus_end, $focus_fset_id ];

        }
        if ($opts{D}) {
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

    if ($opts{D}) {
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

        my ($start, $end, $fset_id, 
            $focus_start, $focus_end, $focus_fset_id) = @$ft;

        # initialize overlap string if not exists
        if (! exists $regulatory_features{$focus_start}{$focus_end}) {
            foreach my $name (keys %target_fsets) {
                #next if (exists $focus_fsets{$name});
                $regulatory_features{$focus_start}{$focus_end}{$name} = 0;
            }
        }

        my $fset = $fsa->fetch_by_dbID($fset_id);
        #print Dumper $fset->name;
        $regulatory_features{$focus_start}{$focus_end}{$fset->name}++;
    }

    return \%regulatory_features;
}


##########################
sub write_intersection_features()
##########################
{
    my ($overlap_features) = @_;
    my (%i_fsets, %i_ftypes, @i_features);


    my $analysis = Bio::EnsEMBL::Analysis->new
        (
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
         -displayable     => 1
         );
    $analysis = $aa->fetch_by_dbID($aa->store($analysis));

    foreach my $ft (@$overlap_features) {
        
        my ($start, $end, $fset_id, 
            $focus_start, $focus_end, $focus_fset_id) = @$ft;

        next if ($focus_fset_id == $fset_id);

        my ($i_start,$i_end) = (sort {$a<=>$b} ($focus_start,
                                                $focus_end,
                                                $start,
                                                $end))[1,2];

        my $i_ftype_name = 'Intersection';
        my $i_fset_name = $i_ftype_name.'_'.
            join(':',
                 sort($fsa->fetch_by_dbID($focus_fset_id)->name(),
                      $fsa->fetch_by_dbID($fset_id)->name()));

        my $i_feature = Bio::EnsEMBL::Funcgen::PredictedFeature->new
            (
             -slice  => $slice,
             -start  => $i_start,
             -end    => $i_end,
             -strand => 0,
             -feature_set => &get_FeatureSet($i_fset_name, \%i_fsets, 
                                             $i_ftype_name,\%i_ftypes, 
                                             $analysis)	
             );

        push(@i_features, $i_feature);
    }
    
    $pfa->store(@i_features);
}

##########################
sub get_FeatureSet()
##########################
{
    my ($fset_name, $fsets, $ftype_name, $ftypes, $analysis) = @_;


    if (! exists $fsets->{$fset_name}) {
        
        if (defined ($fsets->{$fset_name} = $fsa->fetch_by_name($fset_name))) {
            
            if ($opts{c}) {
                my $sql = 'DELETE from predicted_feature where feature_set_id='.
                    $fsets->{$fset_name}->dbID();
                
                $db->dbc->do($sql) 
                    or throw('Failed to roll back predicted_features for feature_set_id'.
                             $fsets->{$fset_name}->dbID());
            } else {
                throw("There is a pre-existing FeatureSet '$fset_name'.\n".
                      'You must specify clobber (option -c) is you want to delete'.
                      ' and overwrite all pre-existing PredictedFeatures');
            }

        } else {
   
            warn("Feature set '$fset_name' not defined in database!\n".
                 " Generating new FeatureSet.");
            
            if (! exists $ftypes->{$ftype_name}) {
                if (! defined ($ftypes->{$ftype_name} = $fta->fetch_by_name($ftype_name))) {
                    
                    $ftypes->{$ftype_name} = Bio::EnsEMBL::Funcgen::FeatureType->new
                        (
                         -name        => $ftype_name,
                         -description => "$ftype_name of $fset_name",
                         );
                    
                    ($ftypes->{$ftype_name}) = @{$fta->store($ftypes->{$ftype_name})};
                }
            }
            
            $fsets->{$fset_name} = Bio::EnsEMBL::Funcgen::FeatureSet->new
                (
                 -analysis     => $analysis,
                 -feature_type => $ftypes->{$ftype_name},
                 -name         => $fset_name,
                 );
            
            ($fsets->{$fset_name}) = @{$fsa->store($fsets->{$fset_name})};
                
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

    return $fsets->{$fset_name};

}




1;
