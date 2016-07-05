#!/bin/bash

export PERL5LIB=$PWD/bioperl-live-bioperl-release-1-2-3:$PWD/ensembl/modules:$PWD/ensembl-test/modules:$PWD/modules

echo "Running test suite"
echo "Using $PERL5LIB"

#skip these tests
SKIP_TESTS="--skip MultiTestDB.t,CoordSystem.t,Array_ArrayChip.t,Storable.t,RegulatoryFeature.t,InputSet_Set_BaseAdaptor.t,Channel.t,BaseFeatureAdaptor.t,BindingMatrix_MotifFeature.t,Annotated_SetFeatureAdaptor.t,DNAMethylationFeature.t"

# SetFeature.t,Set.t,SegmentationFeature.t,ResultFeature.t,ProbeSet.t,ProbeFeature.t,MirnaTargetFeature.t,InputSubset.t,ExternalFeature.t,Experiment.t,ExperimentalGroup.t,DNAMethylationFeature.t,DataSet.t,CellType.t,Array.t,ArrayChip.t,AnnotatedFeature.t,Probe.t,feature_class_Set.t,FeatureSet.t,FeatureType.t,ResultSet.t,

if [ "$COVERALLS" = 'true' ]; then
  PERL5OPT='-MDevel::Cover=+ignore,bioperl,+ignore,ensembl-test,+ignore,ensembl,+ignore,\.t$$,+ignore,Adaptor.pm$$,+ignore,EFGUtils.pm$$' perl $PWD/ensembl-test/scripts/runtests.pl -verbose $PWD/modules/t/ $SKIP_TESTS
else
  perl $PWD/ensembl-test/scripts/runtests.pl $PWD/modules/t/ $SKIP_TESTS
fi

rt=$?

if [ $rt -eq 0 ]; then
  if [ "$COVERALLS" = 'true' ]; then
    echo "Running Devel::Cover coveralls report"
    cover --nosummary -report coveralls
  fi
  exit $?
else
  exit $rt
fi