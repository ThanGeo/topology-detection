# Topology Detection using APRIL
paper TBD...

## instructions
to run for specific predicate:
```
./sj -p 1000 -c -f -q -t <pred> <datasetA> <datasetB>
```
predicates are selected from the following list: 
{'crosses', 'meet', 'covers', 'intersect', 'covered_by', 'contains', 'within', 'equal', 'disjoint'}

to run for all predicates:
```
./sj -p 1000 -c -f -q <datasetA> <datasetB>
```

CROSSES predicate is not yet implemented
