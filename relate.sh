#! /bin/bash

set -e

# build the project
./build.sh
echo 

# dataets to run
declare -a dataset_pairs=("T1NA T2NA"
                         "T1NA T3NA"
                         "O5EU O6EU")

# relate predicates
declare -a relate_predicates=("meet"
                              "equal"
                              "inside"
                              "contains"
                              "covers"
                              "covered_by"
                         )

# go to build directory
cd build

# run each dataset with the 4 combinations of ST2, OP2, ST3, OP3
for i in "${dataset_pairs[@]}"
do
    pair=($i);
    echo "********" Datasets: ${pair[0]} ${pair[1]} "********"

    echo "---- find relation ----"
    echo FR-OP3:
    ./main -q find_relation -s OP3 -R ${pair[0]} -S ${pair[1]}

    for j in "${relate_predicates[@]}"
    do
        relation=($j)
        echo "----" ${relation} "----"
        echo RE-OP3:
        ./main -q ${relation} -s OP3 -R ${pair[0]} -S ${pair[1]}
    done
done

# return to main dir
cd ..