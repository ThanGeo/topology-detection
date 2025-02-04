# Scalable Topology Detection paper - experiments code

## Datasets
Sample datasets to reproduce some of the experiments can be found in the following Google Drive link: 
https://drive.google.com/drive/folders/1AMOy2q9NGFnJ1eSXXJU82enFmTYXxryV?usp=drive_link

Download and place the following files in the **data/vector_data/** directory:
- T1NA_fixed_binary.dat
- T1NA_offset_map.dat
- T2NA_fixed_binary.dat
- T2NA_offset_map.dat
- O5_Europe_fixed_binary.dat
- O5_Europe_offset_map.dat
- O6_Europe_fixed_binary.dat
- O6_Europe_offset_map.dat

The files are in a custom binary format descibed in the 'GEOMETRY BINARY FORMAT.txt' file in the Google Drive link.
Each dataset is accompanied by an 'offset map' file that allows for faster retrieval of the geometry from the disk. 
Both the geometry and the offset map should be stored in the **data/vector_data/** directory.

To run your own dataset, transform your data into the described format and store them in the directory. Then, edit the scripts to include them.
All datasets are described by a keyword in the ```datasets.ini``` configuration file. The dataset's keyword is passes as argument in the program as follows:

```
./main -q find_relation -s OP3 -R dataset1 -S dataset2

```
Performs a find relation join between dataset1 and dataset2. Both keywords must exist in the datasets.ini and their respective fields must be set.
```-s OP3```: instructs the program to run the P+C filtering.
```-s ST3```: it runs the naive APRIL intersection filtering.
```-s ST2```: it runs the 2-step pipeline (no intermediate filters).
```-s OP3```: it runs the 2-step pipeline with an optimized MBR filter for find relation queries.


## Build
To manually build (not needed for tests, they build automatically), 
run ```source ./build.sh``` in main directory or 

```
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
cmake --build .
cd build
make
```

## Requirements
- Boost Geometry
- Enough RAM to load APRIL in-memory (see paper). If not enough memory, edit scripts and remove the dataset scenarios that it fails with and run the rest. 

## Experiments

To allow the bash scripts for the experiments to execute, run 
```
chmod +x *.sh
```

### create APRIL
First run the **setup_experiments.sh** script, to create and store on disk all APRIL approximations used in the experiments.
It will print some results because the creation of APRIL isn't decoupled from evaluating a query (yet), ignore. 

### reproduce experiments

#### To reproduce the numbers of Figure 7, run the following script:
```
./find_relation.sh
```

#### To run the experiment of Figure 8, run the following script:
```
./scalability.sh
```

#### To compare relate with find relation run:
```
./relate.sh
```
