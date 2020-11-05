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
#include <vector>
#include <numeric>
#include <queue>
#include <ostream>
#include <map>
#include <algorithm>

// 3rd party libs
#include "mpi.h"
#include <metislib.h>


constexpr std::size_t MAX_COLOR = 1000;

/* forward declare functions */
graph_t *getGraph(std::string const& filename);
graph_t *tmpGetGraph();
void colorGraph(graph_t *graph, std::size_t s);
void colorGraphSerial(std::vector<idx_t> const& U, graph_t *graph);
std::vector<std::vector<idx_t>> partitionU(std::vector<idx_t> const& U, std::size_t s);
std::map<int, std::vector<std::tuple<int, int, int>>> getBoundaryVertices(graph_t *graph, std::vector<idx_t> const& U);
void outputColoring(graph_t *graph, std::ostream &out);

int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);

    /* primitive MPI info */
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    /* read in args */
    std::string fileName = std::string(argv[1]);
    std::size_t chunkSize = std::stoul(argv[2]);

    /* Read in graph structure and partition */
    graph_t *graph = getGraph(fileName);

    /* color graph in parallel */
    double start = MPI_Wtime();
    colorGraph(graph, chunkSize);
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

    /* clean up */
    delete graph;

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

    return tmpGetGraph();
    //return static_cast<graph_t *>(nullptr);
}


/** For testing until graph partitioning part is done.
 *  @returns a graph_t object
 */
graph_t *tmpGetGraph() {
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (size != 4) {
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
        graph->xadj = new idx_t[NVERTS+1] {0, 2, 5, 7, 10, 14};
        

    } else if (rank == 1) {

        constexpr idx_t NEDGES = 12;
        graph->nedges = NEDGES;
        graph->adjncy = new idx_t[NEDGES] {1, 1, 3, 0, 4, 3, 0, 2, 4, 1, 3, 0};
        graph->adjwgt = new idx_t[NEDGES] {0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 3};
        graph->xadj = new idx_t[NVERTS+1] {0, 3, 5, 6, 9, 12};

    } else if (rank == 2) {

        constexpr idx_t NEDGES = 13;
        graph->nedges = NEDGES;
        graph->adjncy = new idx_t[NEDGES] {4, 1, 2, 0, 2, 4, 0, 1, 3, 2, 1, 0, 4};
        graph->adjwgt = new idx_t[NEDGES] {0, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3};
        graph->xadj = new idx_t[NVERTS+1] {0, 3, 6, 9, 10, 13};

    } else if (rank == 3) {

        constexpr idx_t NEDGES = 11;
        graph->nedges = NEDGES;
        graph->adjncy = new idx_t[NEDGES] {4, 4, 1, 0, 2, 1, 3, 2, 4, 4, 3};
        graph->adjwgt = new idx_t[NEDGES] {1, 2, 3, 3, 3, 3, 3, 3, 3, 2, 3};
        graph->xadj = new idx_t[NVERTS+1] {0, 3, 5, 7, 9, 11};

    }

    return graph;
}


/** Does parallel graph coloring.
 * @param[in,out] graph the graph to be colored
 */
void colorGraph(graph_t *graph, std::size_t s) {
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // (4)
    std::vector<idx_t> U (graph->nvtxs);
    std::iota(U.begin(), U.end(), 0);

    int size_U = U.size();
    int max_size_U = 1;
    while (max_size_U != 0) {   // (5)
        std::map<int, std::vector<std::tuple<int, int>>> boundary_vertices_for_other_ranks;
        std::map<int, std::vector<std::tuple<int, int>>> boundary_vertices_from_other_ranks;
        std::map<int, std::vector<std::tuple<MPI_Request, MPI_Request>>> send_requests;
        std::map<int, std::vector<std::tuple<MPI_Request, MPI_Request>>> recv_requests;
        if (size_U != 0) {      // (6)
            // Partition U into subsets
            std::vector<std::vector<idx_t>> Us = partitionU(U, s);

            // Supersteps -- color each subset of U
            for (idx_t k = 0; k < graph->nvtxs / s; k++) {     // (8)
                colorGraphSerial(Us.at(k), graph); // (9-10)

                // Communicate with other processors
                // 1. find all boundary vertices
                std::vector<idx_t> boundaryColors (graph->nedges);
                for (idx_t i : Us.at(k)) {
                    for (int j = graph->xadj[i]; j < graph->xadj[i+1]; j++) {    // for every edge
                        idx_t neighborIdx = graph->adjncy[j];   // this is their local index
                        idx_t neighborRank = graph->adjwgt[j];

                        if (neighborRank != rank) { // is a boundary vertex
                            
                            

                            std::tuple<int, int> f = std::make_tuple(i, graph->vwgt[i]);
                            boundary_vertices_for_other_ranks.at(neighborRank).push_back(f);

                            MPI_Request send_request1, send_request2;
                            std::tuple<int, int> f = std::make_tuple(send_request1, send_request2)
                            send_requests.at(neighborRank).push_back(f);
                            MPI_Isend(boundary_vertices_for_other_ranks.at(neighborRank).at(0), 1, MPI_INT, 
                                neighborRank, tag, MPI_COMM_WORLD, send_requests.at(neighborRank).back().at(0));
                            MPI_Isend(boundary_vertices_for_other_ranks.at(neighborRank).at(1), 1, MPI_INT, 
                                neighborRank, tag, MPI_COMM_WORLD, send_requests.at(neighborRank).back().at(1));
                            
                            MPI_Request recv_request1, recv_request2;
                            std::tuple<int, int> f = std::make_tuple(recv_request1, recv_request2)
                            recv_requests.at(neighborRank).push_back(f);
                            MPI_Irecv(boundary_vertices_from_other_ranks.at(neighborRank).at(0), 1, MPI_INT, 
                                neighborRank, tag, MPI_COMM_WORLD, recv_requests.at(neighborRank).back().at(0));
                            MPI_Irecv(boundary_vertices_from_other_ranks.at(neighborRank).at(1), 1, MPI_INT, 
                                neighborRank, tag, MPI_COMM_WORLD, recv_requests.at(neighborRank).back().at(1));
                        }
                    }
                }
                
            }
            


            
        }

        size_U = U.size();
        MPI_Allreduce(&size_U, &max_size_U, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    }
}


/** Does the serial graph coloring on a subgraph U.
 * Currently it's a naive breadth first search coloring method.
 * @param[in,out] U the set of vertices to be colored in this superstep.
 * @param[in,out] graph the whole graph to be colored in this rank.
 */
void colorGraphSerial(std::vector<idx_t> const& U, graph_t *graph) {
    bool neighbor_colors[MAX_COLOR];
    for (int i = 0; i < U.size(); i++) {

        idx_t vertex = U.at(i);
        std::fill(neighbor_colors + 0, neighbor_colors + MAX_COLOR, false);

        for (int j = graph->xadj[vertex]; j < graph->xadj[vertex + 1]; j++) {
            idx_t neighbor = graph->adjncy[j];
            neighbor_colors[graph->vwgt[neighbor]] = true; // this color is used
        }
        for (std::size_t j = 1; j < MAX_COLOR; j++) {
            if (!neighbor_colors[j]) { // first color not used
                graph->vwgt[vertex] = j;
                break;
            }
        }
    }
}

std::vector<std::vector<idx_t>> partitionU(std::vector<idx_t> const& U, std::size_t s) {
    std::vector<std::vector<idx_t>> Us (U.size() / s);

    for (std::size_t i = 0; i < U.size(); i += s) {
        auto last = std::min(U.size(), i + s);
        Us.emplace_back(U.begin() + i, U.begin() + last);
    }

    return Us;
}

std::map<int, std::vector<std::tuple<int, int>>> getBoundaryVertices(graph_t *graph, std::vector<idx_t> const& U) {
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    std::map<int, std::vector<std::tuple<int, int, int>>> boundaryVertices;


    for (int i : U) {    // for every vertex
        for (int j = graph->xadj[i]; j < graph->xadj[i+1]; j++) {    // for every edge
            idx_t neighborIdx = graph->adjncy[j];   // this is their local index
            // idx_t neighborRank = graph->adjwgt[j];

            if (neighborRank != rank) { // is a boundary vertex
                std::tuple<int, int, int> f = std::make_tuple(i, graph->vwgt[i]);
                boundaryVertices.at(neighborRank).push_back(f);
            }
        }
    }
    return boundaryVertices;
}

/** Outputs the coloring info of a graph into a stream.
 * @param[in] graph the graph whose coloring will be outputted
 * @param[in] out the stream to output coloring info into
 */
void outputColoring(graph_t *graph, std::ostream &out) {

    /* [TODO] -- find a good way to output coloring so we can
                 check against other colorings */

}