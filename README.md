# CMSC818X-Project
CMSC818X Group Project

### Dependency
- Metis-5.1.0

### Installation
To install metis, simply
```bash
./INSTALL_METIS.sh
```

Run `make` in the `src` directory to build the project code into `./bin`.


### Running

Once built, the application can be run from command line with a graph input file and a chunk size (_s_ in the algorithm). For instance,

```bash
mpirun -np 4 ./bin/run ./metis-5.1.0/graphs/4elt.graph 10
```

for which the output might look like 

```
... random METIS jargon here ...
Rank 0 Duration:  0.140971
Min Duration:     0.14097
Max Duration:     0.140976
Average Duration: 0.140972
Colors Used:      7
Num Iterations    4
```

