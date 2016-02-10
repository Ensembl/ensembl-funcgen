#!/bin/bash

export PERL5LIB=$PWD/bioperl-live-bioperl-release-1-2-3:$PWD/ensembl/modules:$PWD/ensembl-test/modules:$PWD/modules

echo "Running test suite"
echo "Using $PERL5LIB"

#skip these tests
SKIP_TESTS="--skip MultiTestDB.t,CoordSystem.t,BindingMatrix_MotifFeature.t,Array_ArrayChip.t,Storable.t,ResultSet.t,RegulatoryFeature.t,Probe.t,InputSet_Set_BaseAdaptor.t,FeatureType.t,FeatureSet.t,Channel.t,BaseFeatureAdaptor.t,Annotated_SetFeatureAdaptor.t"

# SetFeature.t Set.t SegmentationFeature.t ResultFeature.t ProbeSet.t ProbeFeature.t MirnaTargetFeature.t InputSubset.t ExternalFeature.t Experiment.t ExperimentalGroup.t DNAMethylationFeature.t DataSet.t CellType.t Array.t ArrayChip.t AnnotatedFeature.t

# skipping existing failing tests
# SKIP_TESTS="--skip Annotated_SetFeatureAdaptor.t,BaseFeatureAdaptor.t,BindingMatrix_MotifFeature.t,CoordSystem.t,FeatureSet.t,FeatureType.t,InputSet_Set_BaseAdaptor.t,MultiTestDB.t,RegulatoryFeature.t,ResultSet.t,Storable.t"

if [ "$COVERALLS" = 'true' ]; then
  PERL5OPT='-MDevel::Cover=+ignore,bioperl,+ignore,ensembl-test,+ignore,ensembl,+ignore,\.t$$' perl $PWD/ensembl-test/scripts/runtests.pl -verbose $PWD/modules/t/ $SKIP_TESTS
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