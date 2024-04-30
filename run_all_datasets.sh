#! /bin/bash

# build the project
./build.sh

# set directories
mkdir -p experiment_results
csv_file="experiment_results/all_datasets.csv"
output_file="experiment_results/all_datasets.txt" 
build_dir_output_file="../${output_file}" 
build_dir_csv_file="../${csv_file}" 
echo -e "output \n" > $output_file

# dataset combinations to run
declare -a dataset_pairs=("T1NA T2NA"
                          "T1NA T3NA"
                          "O5AF O6AF"
                          "O5AS O6AS"
                          "O5EU O6EU"
                          "O5NA O6NA"
                          "O5OC O6OC"
                          "O5SA O6SA"
                          "T3NA T10NA"
                          "O3EU O6EU"
                          "O3NA O6NA"
                         )

# go to build directory
cd build

# run each dataset with the 4 combinations of filters (normal APRIL, no compression, 1 partition)
for i in "${dataset_pairs[@]}"
do
    pair=($i);
    echo ${pair[0]} ${pair[1]} >> $build_dir_output_file
    ./main -q find_relation -c -f APRIL_STANDARD -R ${pair[0]} -S ${pair[1]} -v $build_dir_csv_file -e >> $build_dir_output_file
    ./main -q find_relation -m MBR_FR -f APRIL_FR -R ${pair[0]} -S ${pair[1]} -v $build_dir_csv_file -e >> $build_dir_output_file
    ./main -q find_relation -R ${pair[0]} -S ${pair[1]} -v $build_dir_csv_file >> $build_dir_output_file -e
    ./main -q find_relation -m MBR_FR -R ${pair[0]} -S ${pair[1]} -v $build_dir_csv_file >> $build_dir_output_file -e
done

# run O5O6 with 1, 4 and 8 partitions
# echo O5 O6 >> $build_dir_output_file
# ./main -q find_relation -f APRIL_STANDARD -R ${pair[0]} -S ${pair[1]} -v $build_dir_csv_file -e >> $build_dir_output_file
# ./main -q find_relation -m MBR_FR -f APRIL_FR -R O5 -S O6 -p 1 -v $build_dir_csv_file -e >> $build_dir_output_file
# ./main -q find_relation -m MBR_FR -f APRIL_FR -R O5 -S O6 -p 4 -v $build_dir_csv_file -e >> $build_dir_output_file
# ./main -q find_relation -m MBR_FR -f APRIL_FR -R O5 -S O6 -p 8 -v $build_dir_csv_file -e >> $build_dir_output_file
# ./main -q find_relation -R O5 -S O6 -v $build_dir_csv_file >> $build_dir_output_file -e
# ./main -q find_relation -m MBR_FR -R O5 -S O6 -v $build_dir_csv_file >> $build_dir_output_file -e
    
# return to main dir
cd ..