#!/usr/local/ensembl/bin/perl
use warnings;
use strict;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor;

my $file = shift;
my @exp_id = @ARGV;

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
	-host => "ecs3",
	-dbname => "homo_sapiens_funcgen_42_36d",
	-species => "homo_sapiens",
	-user => "ensadmin",
	-pass => "ensembl",
	-dnadb => $cdb,
	-port => '3307',
);

my $anal_a = $db->get_AnalysisAdaptor();
my $anal = $anal_a->fetch_by_logic_name("TilingHMM");
my $pfa = $db->get_PredictedFeatureAdaptor();
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


	my $slice = $cdb->get_SliceAdaptor()->fetch_by_region('chromosome', $chr);

	my $pfeature = Bio::EnsEMBL::Funcgen::PredictedFeature->new(
		-SLICE         => $slice,
		-START         => $start,
		-END           => $end,
		-STRAND        => 1,
		-DISPLAY_LABEL => $text,
		-ANALYSIS      => $anal,
		-SCORE         => $score,
		-EXPERIMENT_IDS => \@exp_id, 
		-FEATURE_TYPE_ID => 1
	);

	push @p_features, $pfeature;

}

$pfa->store(@p_features);
   							     	
