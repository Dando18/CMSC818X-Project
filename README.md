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

You can run with a randomly generated graph using '-' as the filename. The the next parameters are s, the # of vertices per rank, max degree, and connectivity. 

```bash
mpirun -np 16 ./bin/run - 100 1000000 500 0.2
```

The above will generate a random graph with 1E6 vertices per rank. Each vertex will have a max possible of 500 edges. The connectivity, 0.2, is a number signifying roughly how many edges between ranks there will be. Closer to 0.0 is almost no inter-process edges and close to 1.0 is a lot of inter-process edges. Values between 0.1 and 0.3 are close to how metis would partition. Then the above algorithm is run with s=100 on 16 processes.
