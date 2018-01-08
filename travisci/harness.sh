#!/bin/bash

# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2018] EMBL-European Bioinformatics Institute
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

export PERL5LIB=$PWD/bioperl-live:$PWD/ensembl/modules:$PWD/modules:$PWD/scripts:$PWD/ensembl-compara/modules:$PWD/ensembl-killlist/modules:$PWD/ensembl-hive/modules:$PWD/ensembl-io/modules:$PWD/ensembl-test/modules

export TEST_AUTHOR=$USER

export WORK_DIR=$PWD

echo "Running test suite"
echo "Using $PERL5LIB"

echo "Test list"
pwd
ls -l t

echo "COVERALLS value=$COVERALLS"

if [ "$COVERALLS" = 'true' ]; then
  export PERL5LIB=$PERL5LIB:$PWD/ensembl-test/modules
  PERL5OPT='-MDevel::Cover=+ignore,bioperl,+ignore,ensembl-test' perl $PWD/ensembl-test/scripts/runtests.pl -verbose t $SKIP_TESTS
else
  perl $PWD/ensembl-test/scripts/runtests.pl t $SKIP_TESTS
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
