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
#include <queue>
#include <ostream>
#include <algorithm>

// 3rd party libs
#include "mpi.h"
#include <metislib.h>


constexpr std::size_t MAX_COLOR = 1000;

/* forward declare functions */
graph_t *getGraph(std::string const& filename);
graph_t *tmpGetGraph();
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

    /* a temporary getGraph function for testing */
    graph = tmpGetGraph();

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


/** For testing until graph partitioning part is done.
 *  [TODO] -- remove this!
 *  @returns a graph_t object
 */
graph_t *tmpGetGraph() {
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (size != 3) {
        MPI_Abort(MPI_COMM_WORLD, 0);
    }

    graph_t *graph = new graph_t;
    constexpr idx_t NVERTS = 5;
    graph->nvtxs = 5;
    graph->vwgt = new idx_t[5];


    if (rank == 0) {

        constexpr idx_t NEDGES = 14;
        graph->nedges = NEDGES;
        graph->adjncy = new idx_t[NEDGES] {1, 3, 0, 4, 0, 3, 4, 0, 2, 4, 1, 2, 3, 0};
        graph->adjwgt = new idx_t[NEDGES] {0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 2};
        graph->xadj = new idx_t[NVERT+1] {0, };
        

    } else if (rank == 1) {

    } else if (rank == 2) {

    } else if (rank == 3) {

    }

    return static_cast<graph_t *>(nullptr);
}


/** Does parallel graph coloring.
 * @param[in,out] graph the graph to be colored
 */
void colorGraph(graph_t *graph) {
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    /* [TODO] --implement graph algorithm */
    graph_t graph = tmpGetGraph();
    idx_t U = //range(graph->nvtxs);
    int max_size_U = 1;
    while (max_size_U != 0) {
        // [U1, ..., Un] = partition(graph, s);
        
    }




}

/** Does the serial graph coloring on a subgraph.
 * Currently it's a naive breadth first search coloring method.
 * @param[in,out] graph the graph to be colored.
 */
// void colorGraphSerial_BFS(idx_t *xadj, idx_t *vcolors) {
//     std::queue<idx_t> q;
//     q.push(0);

//     bool neighbor_colors[MAX_COLOR];
//     while (!q.empty()) {
//         idx_t vertex = q.back();
//         q.pop();
//         std::fill(neighbor_colors + 0, neighbor_colors + MAX_COLOR, false);
//         for (int i = xadj[vertex]; i < xadj[vertex + 1]; i++) {
//             if (vcolors[i] != 0) { // not colored yet
//                 q.push(i);
//                 neighbor_colors[vcolors[i]] = true; // this color is used
//             }
//         }

//         for (std::size_t i = 1; i < MAX_COLOR; i++) {
//             if (!neighbor_colors[i]) { // first color not used
//                 vcolors[vertex] = i;
//                 break;
//             }
//         }
//     }
// }

void colorGraphSerial(idx_t *U, int size_U, graph_t graph) {
    bool neighbor_colors[MAX_COLOR];
    for (int i = 0; i < size_U; i++) {
        idx_t vertex = U[i];
        std::fill(neighbor_colors + 0, neighbor_colors + MAX_COLOR, false);
        for (int j = graph->xadj[vertex]; j < graph->xadj[vertex + 1]; j++) {
            idx_t neighbor = graph->adjncy[j];
            neighbor_colors[graph->vwgt[neighbor]] = true; // this color is used
        }
        for (std::size_t j = 1; j < MAX_COLOR; j++) {
            if (!neighbor_colors[j]) { // first color not used
                vcolors[vertex] = j;
                break;
            }
        }
    }
}

/** Outputs the coloring info of a graph into a stream.
 * @param[in] graph the graph whose coloring will be outputted
 * @param[in] out the stream to output coloring info into
 */
void outputColoring(graph_t *graph, std::ostream &out) {

    /* [TODO] -- find a good way to output coloring so we can
                 check against other colorings */

}