#! /bin/bash

set -e

# build the project
./build.sh
echo 

# dataets to run
declare -a dataset_pairs=("T1NA T2NA"
                         "T1NA T3NA"
                         "O5EU O6EU")

# go to build directory
cd build

# run each dataset with the 4 combinations of ST2, OP2, ST3, OP3
for i in "${dataset_pairs[@]}"
do
    pair=($i);
    echo "********" Scenario: ${pair[0]}"-"${pair[1]} "********"
    echo ST2:
    ./main -q find_relation -s ST2 -R ${pair[0]} -S ${pair[1]} 
    echo OP2:
    ./main -q find_relation -s OP2 -R ${pair[0]} -S ${pair[1]}
    echo ST3:
    ./main -q find_relation -z -s ST3 -R ${pair[0]} -S ${pair[1]}
    echo OP3:
    ./main -q find_relation -z -s OP3 -R ${pair[0]} -S ${pair[1]}
    echo
done




# return to main dir
cd ..