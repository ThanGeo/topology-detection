#! /bin/bash

set -e

# build the project
./build.sh
echo 

# go to build directory
cd build

echo "******** Scenario: O5EU-O6EU SCALABILITY OP2 ********"
./main -q find_relation -s SCALABILITY_OP2 -R O5EU -S O6EU
echo

echo "******** Scenario: O5EU-O6EU SCALABILITY OP3 ********"
./main -q find_relation -z -s SCALABILITY_OP3 -R O5EU -S O6EU
echo

# return to main dir
cd ..