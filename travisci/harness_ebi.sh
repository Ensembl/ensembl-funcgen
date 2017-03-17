#!/bin/bash

rm -rf $PWD/modules/t/homo_sapiens.MultiTestDB.frozen.conf
SKIP_TESTS="--skip MultiTestDB.t,CoordSystem.t,Array_ArrayChip.t,Storable.t,RegulatoryFeature.t,InputSet_Set_BaseAdaptor.t,Channel.t,BaseFeatureAdaptor.t,BindingMatrix_MotifFeature.t,Annotated_SetFeatureAdaptor.t,DNAMethylationFeature.t"
perl $PWD/../ensembl-test/scripts/runtests.pl -clean $PWD/modules/t/ $SKIP_TESTS
rt=$?
exit $rt
