#!usr/bin/env perl
use strict;
use warnings;
use Test::More;
use Bio::EnsEMBL::Funcgen::AnnotatedFeature;
use Bio::EnsEMBL::Test::MultiTestDB;

BEGIN {
    $| = 1;

    #plan tests => 15;
}

#obtain Adaptors for dnabb and funcgen databases
my $multi = Bio::EnsEMBL::Test::MultiTestDB->new();
my $db    = $multi->get_DBAdaptor("funcgen");

# switch off the debug prints
our $verbose = 0;

#obtain Adaptors
my $fset_adaptor  = $db->get_adaptor('featureset');
my $af_adaptor    = $db->get_adaptor('annotatedfeature');
my $ftype_adaptor = $db->get_FeatureTypeAdaptor;
my $slice_adaptor = $db->get_adaptor("slice");

#create a test slice object
my $slice = $slice_adaptor->fetch_by_region( undef, 1 );

#Get all FeatureTypers from DB that have defined annotated features

my @histone_feature_types = @{ $ftype_adaptor->fetch_all_by_class('Histone') };
my @tf_feature_types =
  @{ $ftype_adaptor->fetch_all_by_class('Transcription Factor') };
my @open_feature_types =
  @{ $ftype_adaptor->fetch_all_by_class('Open Chromatin') };
my @pol_feature_types = @{ $ftype_adaptor->fetch_all_by_class('Polymerase') };


#Test annotated features from database for all defined FetureTypes
foreach my $ftype (
    @histone_feature_types, @tf_feature_types,
    @open_feature_types,    @pol_feature_types
  )
{
    my @anno_features =
      @{ $af_adaptor->fetch_all_by_Slice_FeatureType( $slice, $ftype ) };

    ok( @anno_features, $ftype->name . ' for annotated features' );

    next if !defined @anno_features;

    my $af1 = $anno_features[ rand @anno_features ];

    ok( defined $af1, $ftype->name . ' AnnotatedFeature object from DB' );

    is(
        $af1->display_label,
        $af1->{'set'}->display_label,
        'display_label for AnnotatedFeature object with dbID ' . $af1->dbID
    );

    ok(
        defined $af1->score,
        'score for '
          . $af1->display_label
          . ' with AnnotatedFeature dbID '
          . $af1->dbID
    );

    ok(
        defined $af1->summit,
        'summit for '
          . $af1->display_label
          . ' with AnnotatedFeature dbID '
          . $af1->dbID
    );

    can_ok( $af1, ('is_focus_feature') );

}

#Get FeatureSet object
my $fset1 =
  $fset_adaptor->fetch_by_name("K562_DNase1_ENCODE_Uw_SWEmbl_R0025_D150");

# create a test AnnotatedFeature object

my $anno_feature = Bio::EnsEMBL::Funcgen::AnnotatedFeature->new(
    -SLICE         => $slice,
    -START         => 1,
    -END           => 100,
    -STRAND        => -1,
    -DISPLAY_LABEL => "DNase1_K562",
    -SUMMIT        => 1900,
    -SCORE         => 1.02,
    -FEATURE_SET   => $fset1,
);

ok( defined $anno_feature, 'new object created' );

ok(
    $anno_feature->isa('Bio::EnsEMBL::Funcgen::AnnotatedFeature'),
    "the object belongs to the right class - AnnotatedFeature"
);

is( $anno_feature->score, 1.02, 'score () method works' );

is( $anno_feature->summit, 1900, 'summit () method works' );

is( $anno_feature->display_label,
    "DNase1_K562", 'display_label () method works' . "\n" );

#store the test AnnotatedFeature in database
#my @stored_afs=$af_adaptor->store($anno_feature);
#ok(@stored_afs,'store new AnnotatedFeature in DB');
