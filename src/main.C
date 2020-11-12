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
#include <set>
#include <numeric>
#include <queue>
#include <ostream>
#include <map>
#include <algorithm>
#include <cstdio>

// 3rd party libs
#include "mpi.h"
#include <metislib.h>


constexpr std::size_t MAX_COLOR = 1000;

struct coloring_stats_t {
    std::size_t numIterations;
};

/* forward declare functions */
graph_t *getGraph(std::string const& filename);
graph_t *tmpGetGraph();
void deleteGraph(graph_t *graph);
void colorGraph(graph_t *graph, std::size_t s, coloring_stats_t &out);
void colorGraphSerial(std::vector<idx_t> const& U, graph_t *graph, std::vector<idx_t> const& boundaryColors);
void colorGraphSerialTest(std::vector<idx_t> const& U, graph_t *graph, int loop, int rank);
std::vector<std::vector<idx_t>> partitionU(std::vector<idx_t> const& U, std::size_t s);
idx_t numColorsUsed(graph_t *graph);
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
    coloring_stats_t stats;
    double start = MPI_Wtime();
    colorGraph(graph, chunkSize, stats);
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

    int totalColors = numColorsUsed(graph);

    if (rank == 0) {
        std::cout << "Rank 0 Duration:  " << duration            << "\n"
                  << "Min Duration:     " << minDuration         << "\n"
                  << "Max Duration:     " << maxDuration         << "\n"
                  << "Average Duration: " << avgDuration         << "\n"
                  << "Colors Used:      " << totalColors         << "\n"
                  << "Num Iterations    " << stats.numIterations << "\n";
    }

    /* clean up */
    deleteGraph(graph);
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
    graph->nvtxs = NVERTS;
    graph->vwgt = new idx_t[NVERTS];
    std::fill(graph->vwgt + 0, graph->vwgt + NVERTS, -1);


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

void deleteGraph(graph_t *graph) {
    delete graph->vwgt;
    delete graph->adjncy;
    delete graph->adjwgt;
    delete graph->xadj;
}


/** Does parallel graph coloring.
 * @param[in,out] graph the graph to be colored
 */
void colorGraph(graph_t *graph, std::size_t s, coloring_stats_t &out) {
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // (4)
    std::vector<idx_t> U (graph->nvtxs);
    std::vector<idx_t> boundaryColors (graph->nedges, -1);
    std::vector<idx_t> previousColors (graph->nvtxs, -1);
    std::iota(U.begin(), U.end(), 0);

    int size_U = U.size();
    int max_size_U = 1;
    int loop = 0;
    while (max_size_U != 0) {   // (5)
        
        if (size_U != 0) {      // (6)
            // Partition U into subsets
            std::vector<std::vector<idx_t>> Us = partitionU(U, s);

            // Communicate with other processors
            // 1. find all boundary vertices
            //std::vector<idx_t> boundaryColors (graph->nedges, -1);
            std::vector<MPI_Request> requests;
            int num_external_edges = 0;
            
            // count all of graph->adjwgt != rank
            for (idx_t i = 0; i < graph->nedges; i++) {
                num_external_edges += (graph->adjwgt[i] != rank) ? 1 : 0;
            }
            requests.reserve(num_external_edges);

            idx_t *recv_buf = new idx_t[3 * num_external_edges];
            std::fill(recv_buf + 0, recv_buf + (3*num_external_edges), -1);
            idx_t *send_buf = new idx_t[3 * num_external_edges];
            int counter = 0;
        
            std::cout << "loop " << loop << " [" << rank << "] boundary colors = ";
            for (auto const& b : boundaryColors) {
                std::cout << b << " ";
            }
            std::cout << std::endl;

            // Supersteps -- color each subset of U
            for (size_t k = 0; k < Us.size(); k++) {     // (8)
                std::copy(graph->vwgt + 0, graph->vwgt + graph->nvtxs, previousColors.begin());
                colorGraphSerial(Us.at(k), graph, boundaryColors); // (9-10)
                //colorGraphSerialTest(Us.at(k), graph, loop, rank); // (9-10)
                
                for (idx_t i : Us.at(k)) {  // (11-12)
                    idx_t myColor = graph->vwgt[i];

                    for (int j = graph->xadj[i]; j < graph->xadj[i+1]; j++) {    // for every edge
                        idx_t neighborIdx = graph->adjncy[j];   // this is their local index
                        idx_t neighborRank = graph->adjwgt[j];

                        bool isConflictPossible = (boundaryColors.at(j) == -1) || (boundaryColors.at(j) == previousColors.at(i));
                        if (neighborRank != rank && isConflictPossible) { // is a boundary vertex

                            send_buf[counter] = i; send_buf[counter + 1] = myColor; send_buf[counter + 2] = neighborIdx;
                        
                            MPI_Request recvRequest, sendRequest;
                            MPI_Isend(&(send_buf[counter]), 3, MPI_INT, neighborRank, 0, MPI_COMM_WORLD, &sendRequest);

                            printf("loop %d [%d] sent (%d, %d, %d) to %d\n", loop, rank, i, myColor, neighborIdx, neighborRank);
                            fflush(stdout);
                            
                            MPI_Irecv(&(recv_buf[counter]), 3, MPI_INT, neighborRank, 0, MPI_COMM_WORLD, &recvRequest);
                            requests.push_back(recvRequest);
                            counter += 3;
                        }
                        
                    }
                }
            }   // end of superstep
                
            MPI_Waitall(requests.size(), requests.data(), MPI_STATUS_IGNORE);
            for (size_t i = 0; i < 3 * requests.size(); i = i + 3) {
                idx_t neighborIdx = recv_buf[i];
                idx_t neighborColor = recv_buf[i + 1];
                idx_t myIdx = recv_buf[i + 2]; 

                printf("loop %d [%d] received (%d, %d, %d)\n", loop, rank, neighborIdx, neighborColor, myIdx);

                for (int j = graph->xadj[myIdx]; j < graph->xadj[myIdx+1]; j++) {
                    if (neighborIdx == graph->adjncy[j]) {
                        boundaryColors.at(j) = neighborColor;
                        break;
                    }
                }
            }
            delete[] recv_buf;
            delete[] send_buf;
            
            std::set<idx_t> R;
            for (idx_t i : U) {
                idx_t myColor = graph->vwgt[i];
                for (int j = graph->xadj[i]; j < graph->xadj[i+1]; j++) {    // for every edge
                    idx_t neighborIdx = graph->adjncy[j];   // this is their local index
                    idx_t neighborRank = graph->adjwgt[j];

                    if (neighborRank != rank) {
                        int neighborValue = neighborRank*graph->nvtxs + neighborIdx;
                        int myValue = rank*graph->nvtxs + i;

                        if (myColor == boundaryColors.at(j)) {
                            if (myValue < neighborValue) {
                                R.insert(i);
                            } else {
                                boundaryColors.at(j) = -1;
                                R.insert(i);
                            }
                        }
                    }
                }
            }
            //U = R;
            U.assign(R.begin(), R.end());
            std::cout << "loop " << loop << " Color of rank [" << rank << "] = [";
            for (idx_t i = 0; i < graph->nvtxs; i++) {
                std::cout << graph->vwgt[i] << " ";
            }
            std::cout << "]" << std::endl;

            std::cout << "loop " << loop << " R of rank [" << rank << "] = [";
            for (idx_t i : R) {
                std::cout << i << " ";
            }
            std::cout << "]" << std::endl;
        }   // end of "if U is not empty"

        size_U = U.size();
        MPI_Allreduce(&size_U, &max_size_U, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);

        fflush(stdout);
        if (rank == 0) {    
            std::cout << "\n\nEND OF WHILE\n\n";
            fflush(stdout);
        }
        loop++;
        //if (loop > 3) max_size_U = 0;
    }   // end of while loop

    out.numIterations = loop;
}


/** Does the serial graph coloring on a subgraph U.
 * Currently it's a naive breadth first search coloring method.
 * @param[in,out] U the set of vertices to be colored in this superstep.
 * @param[in,out] graph the whole graph to be colored in this rank.
 */
void colorGraphSerial(std::vector<idx_t> const& U, graph_t *graph, std::vector<idx_t> const& boundaryColors) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
    std::array<bool, MAX_COLOR> neighbor_colors;
    for (size_t i = 0; i < U.size(); i++) {

        idx_t vertex = U.at(i);
        std::fill(neighbor_colors.begin(), neighbor_colors.end(), false);

        for (int j = graph->xadj[vertex]; j < graph->xadj[vertex + 1]; j++) {
            idx_t neighbor = graph->adjncy[j];
            idx_t neighborRank = graph->adjwgt[j];
            if (neighborRank == rank && graph->vwgt[neighbor] != -1) {
                neighbor_colors.at(graph->vwgt[neighbor]) = true; // this color is used
            } else if (boundaryColors.at(j) != -1) {
                idx_t neighborColor = boundaryColors.at(j);
                neighbor_colors.at(neighborColor) = true; // this color is used
            }
        }
        for (std::size_t j = 0; j < MAX_COLOR; j++) {
            if (!neighbor_colors[j]) { // first color not used
                graph->vwgt[vertex] = j;
                break;
            }
        }
    }
}

void colorGraphSerialTest(std::vector<idx_t> const& U, graph_t *graph, int loop, int rank) {
    if (loop >= 1) {
        for (idx_t i : U) {
            graph->vwgt[i] = (i+1) + (rank % 2);
        }
    } else {
        for (idx_t i : U) {
            graph->vwgt[i] = i+1;
        }
    }
}

std::vector<std::vector<idx_t>> partitionU(std::vector<idx_t> const& U, std::size_t s) {
    std::vector<std::vector<idx_t>> Us;
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    for (std::size_t i = 0; i < U.size(); i += s) {
        auto last = std::min(U.size(), i + s);
        Us.emplace_back(U.begin() + i, U.begin() + last);
    }

    return Us;
}

idx_t numColorsUsed(graph_t *graph) {
    int localMaxColor = *std::max_element(graph->vwgt, graph->vwgt + graph->nvtxs);

    int maxColor = -1;
    MPI_Reduce(&localMaxColor, &maxColor, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);

    return maxColor + 1;
}

/** Outputs the coloring info of a graph into a stream.
 * @param[in] graph the graph whose coloring will be outputted
 * @param[in] out the stream to output coloring info into
 */
void outputColoring(graph_t *graph, std::ostream &out) {

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    /* [TODO] -- find a good way to output coloring so we can
                 check against other colorings */



    /* tmp output */
    for (int i = 0; i < graph->nvtxs; i++) {
        out << "rank[" << rank << "].vertex[" << i << "].color = " << graph->vwgt[i] << "\n";
    }

}