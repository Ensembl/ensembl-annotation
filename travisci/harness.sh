#!/bin/bash
export PERL5LIB=$PWD/bioperl-live-bioperl-release-1-2-3:$PWD/ensembl/modules:$PWD/modules:$PWD/scripts:$PWD/ensembl-compara/modules:$PWD/ensembl-killlist/modules:$PWD/ensembl-hive/modules:$PWD/ensembl-io/modules:$PWD/bioperl-live:$PWD/bioperl-run/lib:

export WORK_DIR=$PWD

echo "Running test suite"
echo "Using $PERL5LIB"
rt=0




if [ "$COVERALLS" = 'true' ]; then
  export PERL5LIB=$PERL5LIB:$PWD/ensembl-test/modules
  PERL5OPT='-MDevel::Cover=+ignore,bioperl,+ignore,ensembl-test' perl $PWD/ensembl-test/scripts/runtests.pl -verbose $PWD/modules/t $SKIP_TESTS
else

fi

if [ $rt -eq 0 ]; then
  if [ "$COVERALLS" = 'true' ]; then
    echo "Running Devel::Cover coveralls report"
    cover --nosummary -report coveralls
  fi
  exit $?
else
  exit $rt
fi
