# Scalable Topology Detection paper - experiments code

## Datasets
Sample datasets to reproduce some of the experiments can be found in the following Google Drive link: 
https://drive.google.com/drive/folders/1McEj_TCbN8G2N8IVUJxG1FqyRduax53J?usp=sharing

Download and place all files in the **data/vector_data/** directory. 

## Build
In main directory, run ```source ./build.sh``` or 
```
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
cmake --build .
cd build
make
```

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
#### To run the experiment of Table 6 on the provided datasets (O5O6EU not included due to large size), run the following script:
```
./relate.sh
```

#### todo: scalability experiment


