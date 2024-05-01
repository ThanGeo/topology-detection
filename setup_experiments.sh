#! /bin/bash

set -e

# build the project
./build.sh

# dataets to run
declare -a dataset_pairs=("T1NA T2NA"
                         "T1NA T3NA"
                         "O5EU O6EU")

# go to build dir
cd build

# create the approximations (will perform some query but ignore that)
for i in "${dataset_pairs[@]}"
do
    pair=($i);
    echo Creating APRIL for datasets: ${pair[0]} ${pair[1]}
    ./main -q equal -c -s OP3 -R ${pair[0]} -S ${pair[1]}
    ./main -q equal -c -z -s OP3 -R ${pair[0]} -S ${pair[1]}
done




# return to main dir
cd ..