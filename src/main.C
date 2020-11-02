/**
 *  Main code file for parallel graph partitioning
 *  @author Benjie Miao
 *  @author Yiheng Xu
 *  @author Daniel Nichols
 *  @author Mucong Ding
 * 
 *  @date November 2020
 */

// STL Headers
#include <iostream>
#include <string>
#include <ostream>

// 3rd party libs
#include "mpi.h"

/* forward declare types 
   TODO -- this can be removed once we include METIS headers */
struct graph_t;

/* forward declare functions */
graph_t *getGraph(std::string const& filename);
void colorGraph(graph_t *graph);
void outputColoring(graph_t *graph, std::ostream &out);

int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);

    /* primitive MPI info */
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    /* Read in graph structure and partition */
    graph_t *graph = getGraph("");

    /* color graph in parallel */
    double start = MPI_Wtime();
    colorGraph(graph);
    double end = MPI_Wtime();

    /* output coloring to file or stdout */ 
    outputColoring(graph, std::cout);

    /* aggregate timing info and output on rank 0 */
    double duration = end - start;
    double avgDuration, maxDuration, minDuration;
    MPI_Reduce(&duration, &avgDuration, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&duration, &minDuration, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
    MPI_Reduce(&duration, &maxDuration, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    avgDuration /= static_cast<double>(size);

    if (rank == 0) {
        std::cout << "Rank 0 Duration:  " << duration
                  << "Min Duration:     " << minDuration
                  << "Max Duration:     " << maxDuration
                  << "Average Duration: " << avgDuration
                  << "\n";
    }

    MPI_Finalize();
    return 0;
}


/** Get graph returns the graph partition for the current rank. 
 * It is intended to be called at the same time on all ranks.
 * @param[in] filename the file to read graph from
 * @return returns a pointer to a graph_t struct with the partitioned graph on this rank
 */
graph_t *getGraph(std::string const& filename) {

    /* [TODO] -- use ReadGraph and partition with MPI
                 or use ParMETIS */

    return static_cast<graph_t *>(nullptr);
}


/** Does parallel graph coloring.
 * @param[in,out] graph the graph to be colored
 */
void colorGraph(graph_t *graph) {
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size):

    /* [TODO] --implement graph algorithm */

}

/** Outputs the coloring info of a graph into a stream.
 * @param[in] graph the graph whose coloring will be outputted
 * @param[in] out the stream to output coloring info into
 */
void outputColoring(graph_t *graph, std::ostream &out) {

    /* [TODO] -- find a good way to output coloring so we can
                 check against other colorings */

}