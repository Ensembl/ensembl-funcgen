#!/usr/local/ensembl/bin/perl
use warnings;
use strict;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Funcgen::FeatureSet;

$| =1;

my $file = shift;
my $ftype_name = shift;
my $ctype_name = shift;
my @exp_id = @ARGV;



#we need to add functionality for feature set import
#do it by dbID? Should we change this to import DataSets too based on ResultSets?
#may need to manually hack this for this release

my $cdb = Bio::EnsEMBL::DBSQL::DBAdaptor->new
(
	-host => "ens-livemirror",
	-dbname => "homo_sapiens_core_42_36d",
	-species => "homo_sapiens",
	-user => "ensro",
	-pass => "",
#	-group => 'funcgen',
	-port => '3306',
);


my $db = Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor->new
(
	-host => "ens-genomics1",
	-dbname => "homo_sapiens_funcgen_43",
	-species => "homo_sapiens",
	-user => "ensadmin",
	-pass => "ensembl",
	-dnadb => $cdb,
	-port => '3306',
);

print "Getting adaptors\n";

my $anal_a = $db->get_AnalysisAdaptor();
my $anal = $anal_a->fetch_by_logic_name("Nessie");
my $ftype = $db->get_FeatureTypeAdaptor->fetch_by_name($ftype_name);
my $ctype = $db->get_CellTypeAdaptor->fetch_by_name($ctype_name);


die 'No valid cell or feature type name' if ! $ftype || ! $ctype; 

print "Got feature adn cell types $ftype $ctype\n";
print "Got analyssis $anal\n";

my $pfa = $db->get_PredictedFeatureAdaptor();

my $fset_adaptor = $db->get_FeatureSetAdaptor();

my $fset = Bio::EnsEMBL::Funcgen::FeatureSet->new
  (
   -CELL_TYPE => $ctype,
   -FEATURE_TYPE => $ftype,
   -ANALYSIS => $anal,
  );


print "Got feature set $fset\n";

($fset) = @{$fset_adaptor->store($fset)};

print "Got stored feature set $fset\n";

my @p_features;

open (FILE, $file) || die "No file\n";
while (my $line = <FILE>){

	chomp $line;
	my @tmp = split /\s+/, $line;
	my $start = $tmp[1];
	my $end = $tmp[2];
	my $score = $tmp[4];
	my $text = "enriched_site";
	my $chr = $tmp[0];

	$chr =~ s/chr//;

#	print STDERR "$cdb->get_SliceAdaptor()->fetch_by_region(\'chromosome\', $chr);\n";

	#cache this
	my $slice = $cdb->get_SliceAdaptor()->fetch_by_region('chromosome', $chr);

	my $pfeature = Bio::EnsEMBL::Funcgen::PredictedFeature->new
	  (
	   -SLICE         => $slice,
	   -START         => $start,
	   -END           => $end,
	   -STRAND        => 1,
	   -DISPLAY_LABEL => $text,
	   #-ANALYSIS      => $anal,
	   -SCORE         => $score,
	   -FEATURE_SET   => $fset,
	   #-EXPERIMENT_IDS => \@exp_id, 
	   #-FEATURE_TYPE => $ftype,
	);

	push @p_features, $pfeature;

}

$pfa->store(@p_features);
   							     	
