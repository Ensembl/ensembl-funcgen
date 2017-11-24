#!/bin/bash

export PERL5LIB=$PWD/bioperl-live-bioperl-release-1-2-3:$PWD/ensembl/modules:$PWD/ensembl-test/modules:$PWD/modules

echo "Running test suite"
echo "Using $PERL5LIB"

#skip these tests
SKIP_TESTS="--skip Set.t,SetFeature.t,Alignment.t,RegulatoryFeature.t,Probe.t,MirnaTargetFeature.t,ProbeFeature.t,FeatureSet.t,ReadFile.t,feature_class_Set.t,ExternalFeature.t,Peak.t,MultiTestDB.t,FeatureSet.t,Array_ArrayChip.t,Storable.t,InputSet_Set_BaseAdaptor.t,BaseFeatureAdaptor.t,BindingMatrix_MotifFeature.t,Annotated_SetFeatureAdaptor.t,DNAMethylationFeature.t,DataSet.t"


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