#!/bin/bash

export PERL5LIB=$PWD/bioperl-live-bioperl-release-1-2-3:$PWD/ensembl/modules:$PWD/ensembl-test/modules:$PWD/modules:$PWD/ensembl-hive/modules

echo "Running test suite"
echo "Using $PERL5LIB"

if [ "$COVERALLS" = 'true' ]; then
  PERL5OPT='-MDevel::Cover=+ignore,bioperl,+ignore,ensembl-test,+ignore,ensembl,+ignore,\.t$$,+ignore,\.t/lib$$,+ignore,EFGUtils.pm$$' perl $PWD/modules/t/run_tests.t
else
  perl $PWD/modules/t/run_tests.t
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