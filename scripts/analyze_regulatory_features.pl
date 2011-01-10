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

analyze_regulatory_features.pl

=head1 SYNOPSIS

analyze_regulatory_features.pl \
    -host HOST \
    -port PORT \
    -user USER \
    -pass PASS \
    -dbname homo_sapiens_funcgen_47_36i \
    -data_version 46_36h \
    -focus FOCUS \
    -target TARGET \
    -indir INDIR \
    -outdir OUTDIR \
	-pattern PATTERN \
	-signature SIGNATURE

=head1 DESCRIPTION

reads in regulatory features from file and performs Chi-square tests based 
counted data by given pattern(s) (for binary string) and signature(s) 
(currently gene proximity only).

=cut

use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;

use Getopt::Long;

my ($pass,$port,$host,$user,$dbname,$species,$data_version,
	$help,$man,$focus,$target,$pattern,$signature,$indir,$outdir,
	$debug);

GetOptions (
            "pass|p=s"         => \$pass,
            "port=s"           => \$port,
            "host|h=s"         => \$host,
            "user|u=s"         => \$user,
            "dbname|d=s"       => \$dbname,
            "species=s"        => \$species,
            "data_version|v=s" => \$data_version,
            "help|?"           => \$help,
            "man|m"            => \$man,
            "focus|f=s"        => \$focus,
            "target|t=s"       => \$target,
            "pattern|r=s"      => \$pattern,
            "signature|s=s"    => \$signature,
            "indir|i=s"        => \$indir,
            "outdir|o=s"       => \$outdir,
            "debug"            => \$debug
            );

### defaults ###
$port = 3306 if !$port;
$species = 'homo_sapiens' if !$species;

### check options ###

throw("Must specify mandatory database hostname (-host).\n") if ! defined $host;
throw("Must specify mandatory database username. (-user)\n") if ! defined $user;
throw("Must specify mandatory database password (-pass).\n") if ! defined $pass;
throw("Must specify mandatory database name (-dbname).\n") if ! defined $dbname;
throw("Must specify mandatory database data version, like 47_36i (-data_version).\n") 
    if !$data_version;

throw("Must specify mandatory focus sets (-focus).\n") if ! defined $focus;
throw("Must specify mandatory target sets (-target).\n") if ! defined $target;

throw("Must specify mandatory directory containing file to read (-indir).\n") 
    if !$indir;
throw("Must specify mandatory directory containing file to read (-outdir).\n") 
    if !$outdir;

$| = 1;

use Bio::EnsEMBL::Utils::Exception qw(throw);
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Funcgen::Utils::EFGUtils qw(open_file);
use Bio::EnsEMBL::Funcgen::RegulatoryFeature;

# use ensembldb as we may want to use an old version

my $cdb = Bio::EnsEMBL::DBSQL::DBAdaptor->new
    (
     -host => 'ensembldb.ensembl.org',
     -port => 3306,
     -user => 'anonymous',
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

my $sa = $db->get_SliceAdaptor();
my $fsa = $db->get_FeatureSetAdaptor();
my $afa = $db->get_AnnotatedFeatureAdaptor();
my $rfa = $db->get_RegulatoryFeatureAdaptor();

# parse focus and target sets and check that they exist
my (%focus_fsets, %target_fsets);
map { my $fset = $fsa->fetch_by_name($_);
      $focus_fsets{$fset->dbID} = $fset; 
      throw("Focus set $_ does not exist in the DB") 
          if (! defined $focus_fsets{$fset->dbID}); 
  } split(',', $focus);
#print Dumper %focus_fsets;

map { 
    my $fset = $fsa->fetch_by_name($_);
    $target_fsets{$fset->dbID()} = $fset; 
    throw("Target set $_ does not exist in the DB") 
        if (! defined $target_fsets{$fset->dbID}); 
} split(',', $target);
#print Dumper %target_fsets;

# make sure that target sets also contain focus sets
map { $target_fsets{$_} = $focus_fsets{$_} } keys %focus_fsets;

my @fsets = sort {$a->name cmp $b->name} values %target_fsets;
#map { print Dumper $_->name } @fsets;

# get featrure set object
#my $fset = $fsa->fetch_by_name('RegulatoryFeatures');
#print Dumper $fset;

my (@rf);

# all toooo slow
#foreach my $slice (@{$sa->fetch_all('toplevel')})
#{
#   next if ($slice->seq_region_name =~ m/^NT/);
#   print Dumper $slice->name;
#   
#   push (@rf, $fset->get_Features_by_Slice($slice));
#
#
#}

my $R = "/usr/bin/R";

opendir(DIR, $indir) 
    or throw ("Can't open directory $indir");
my @infiles = grep /\.rf$/, readdir DIR;
closedir DIR;


my ($fg_sr_id, $sr_name, $focus_start, $focus_end, $attrib_start, 
	$attrib_end, $string, $attribute_ids, $sig);

my (@pattern, @signature, $total, %total_p,	%total_s, %total_not_s, %count);

# parse regular expressions to be analyzed
@pattern = split(",",$pattern) if ($pattern);
@signature = split(",",$signature) if ($signature);


foreach my $file (@infiles) {

    my $fh = open_file($indir.'/'.$file);

    while (<$fh>) {

        chomp;

        ($fg_sr_id, $sr_name, $focus_start, $focus_end, 
		 $attrib_start, $attrib_end, $string, $attribute_ids) = split(/\t/);

		#strip off last two genic bits
		$string =~ m/^(.+)(.{2})$/;
		$string = $1;
		$sig = $2;

		@pattern = ( $string ) if (! $pattern);
		@signature = ( $sig ) if (! $signature);

		$total++;

		foreach my $p (@pattern) {

			$total_p{$p}++ if ($string =~ /^$p$/);
			
			foreach my $s (@signature) {
				
				$total_s{$s}++ if ($sig =~ /^$s$/);
				
				
				if (! exists $count{$p}{$s}) {
					$count{$p}{$s} = [[0,0],[0,0]];
				}

				if ($string =~ /^$p$/ && $sig =~ /^$s$/) {
					$count{$p}{$s}[0][0]++;
				} elsif ($string !~ /^$p$/ && $sig =~ /^$s$/) {
					$count{$p}{$s}[0][1]++;
				} elsif ($string =~ /^$p$/ && $sig !~ /^$s$/) {
					$count{$p}{$s}[1][0]++;
				} elsif ($string !~ /^$p$/ && $sig !~ /^$s$/) {
					$count{$p}{$s}[1][1]++;
				} else {
					throw("Bugger!");
				}

			}

		}

	}
	
}

foreach my $p (keys %count) {
	
	foreach my $s (keys %{$count{$p}}) {
		if (!$pattern) {
			$count{$p}{$s}[0][1] = $total_s{$s} - $count{$p}{$s}[0][0];
			my $total_not_signature = $total - $total_s{$s};
			$count{$p}{$s}[1][1] = $total_not_signature - $count{$p}{$s}[1][0];
		}
		
		if (!$signature) {
			$count{$p}{$s}[1][0] = $total_p{$p} - $count{$p}{$s}[0][0];
			my $total_not_pattern = $total - $total_p{$p};
			$count{$p}{$s}[1][1] = $total_not_pattern - $count{$p}{$s}[0][1];
		}
		
	}
	
}

#print Dumper (%total_p,%total_s,%count);

#		my $feature_slice = $sa->fetch_by_region('chromosome',$sr_name,$start, $end);
#		my $expanded_slice = $feature_slice->expand(2500,2500);
#		my @genes = @{$expanded_slice->get_all_Genes()};
#
#		map { 
#			print join("\t", $_->display_id, $_->external_name, 
#					   $_->type, $_->biotype, $binstring), "\n"
#		} @genes;

print "total number of reg. features:\t", $total, "\n";


my $chisq;
foreach my $p (sort {$count{$b} <=> $count{$a}} keys %count){

	my @names = &get_fset_names($p);

	foreach my $s  (sort keys %{$count{$p}}) {
		
		my $R_cmd = sprintf 
			"chisq.test(matrix(c(%d,%d,%d,%d),2,2))",#"fisher.test(x), 
			$count{$p}{$s}[0][0], $count{$p}{$s}[0][1],
			$count{$p}{$s}[1][0], $count{$p}{$s}[1][1];
		
		###print $R_cmd, "\n";
		open(R, "echo \'$R_cmd\' | $R --no-save --slave 2>&1 |" )
				or die "Can't exec comand.";
		while (<R>) {
			chomp;
			$chisq = $_ if (/^X-squared/);
			print $_, "\n" if ($debug);
		}
		close R;

		if ($debug) {
			printf "%9d  %6d  %6d\n"x3,
			$count{$p}{$s}[0][0], $count{$p}{$s}[0][1], $total_s{$s},
			$count{$p}{$s}[1][0], $count{$p}{$s}[1][1], $total-$total_s{$s},
			$total_p{$p},         $total-$total_p{$p},  $total;
		}

		printf "%32s\t%s\t%s\t%s\n", $p, $s, $chisq, join(',', @names);

	}

	print "\n";

}


sub get_fset_names () 
{

	my ($p) = @_;
	my @bits = split('', $p);


	my @names;
	for (my $i=0; $i<=$#bits; $i++) {

		push(@names, $fsets[$i]->name) if $bits[$i] eq 1;

	}

	return @names;

}

exit;

__END__
map {

	my @bits = split('', $_);

	my @fset_names = ();
	for (my $i=0; $i<=$#bits; $i++) {

		push(@fset_names, $fsets[$i]->name) if $bits[$i] eq 1;

	}
	printf "%s\t%8d\t%s\n", $_, $total_r{$_},
	join(',', @fset_names), join("\t", keys %{$count{$_}});


} sort {$total_r{$b} <=> $total_r{$a}} keys %total_r;





1;
