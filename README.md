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
chmod 777 *.sh
```

### create APRIL
First run the **setup_experiments.sh** script, to create and store on disk all APRIL approximations used in the experiments.
It will print some results because the creation of APRIL isn't decoupled from evaluating a query (yet), ignore. 

### reproduce experiments

#### To reproduce the numbers of Figure 8, run the following script:
```
./find_relation.sh
```
#### To run the experiment of Table 6, run the following script:
```
./relate.sh
```

#### To run the experiment of Figure 9, run the following script:
```
./scalability.sh
```


