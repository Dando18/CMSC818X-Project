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
#include <fstream>
#include <sstream>
#include <map>
#include <algorithm>
#include <cstdio>

// 3rd party libs
#include "mpi.h"
// #include <metislib.h>
#include "metisbin.h"
#include <metis.h>

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

// Graph partitioning stuff
idx_t simpleReadGraph(std::string filename, graph_t *graph);
void setParams(params_t *params, graph_t *graph, idx_t nparts);
void prepareForPartition(idx_t options[], params_t *params, int starting_index = 0);
void printGraph(graph_t *graph);
void GPPrintInfo(params_t *params, graph_t *graph);

int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);

    /* check args to prevent undefined abortion */
    if (argc < 3){
        std::cout << "Usage:" << argv[0]  
                  << " [graph filename] [chunkSize]" << std::endl;
                //   << " [graph filename]" << std::endl;
        // MPI_Abort(MPI_COMM_WORLD, 1);
    }

    /* primitive MPI info */
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    /* read in args */
    std::string fileName = std::string(argv[1]);
    std::size_t chunkSize = std::stoul(argv[2]);

    /* Read in graph structure and partition */
    graph_t *graph = getGraph(fileName);
    graph->vwgt = new idx_t[graph->nvtxs];
    std::fill(graph->vwgt + 0, graph->vwgt + graph->nvtxs, -1);

    // For debug
    for (int proc=0; proc<size; proc++) {
        if (rank == proc) {
            std::cout << "=========================" << std::endl; 
            std::cout << "Print Graph in process " << proc << std::endl;
            //printGraph(graph);
            std::cout << "=========================" << std::endl; 
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }

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
    /* The rank 0 reads the graph file, calls the partition
        and then send to each rank */
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    
    graph_t *graph_local = new graph_t;
    graph_t *graph;
    if (rank == 0){
        graph = new graph_t;
        params_t *params = new params_t ;

        idx_t starting_index = simpleReadGraph(filename.c_str(), graph);
        std::cout << "Starting index" << starting_index << "\tSize" << size << std::endl;

        setParams(params, graph, size);
        GPPrintInfo(params, graph);
        idx_t options[METIS_NOPTIONS];
        prepareForPartition(options, params, starting_index);

        idx_t objval;
        idx_t *part = new idx_t[graph->nvtxs];
        int status = 0;
        getchar();
        status = METIS_PartGraphKway(&graph->nvtxs, &graph->ncon, graph->xadj, 
                graph->adjncy, graph->vwgt, graph->vsize, graph->adjwgt, 
                &params->nparts, params->tpwgts, params->ubvec, options, 
                &objval, part);
        for (idx_t i=0; i<graph->nvtxs; i++){
            part[i] -= starting_index;
            //std::cout << "Node:" << i << "\tPart:" << part[i] << std::endl; 
        }
        getchar();
        // broadcast
        // set up a mapping from global_idx to (rank, local_idx)
        std::map<idx_t, std::pair<idx_t, idx_t> > g_to_l;
        std::vector<idx_t> local_num(size, 0);
        for (idx_t i=0; i<graph->nvtxs; i++) {
            g_to_l[i] = std::make_pair(part[i], local_num[part[i]]++);
        }
        // for (auto x:g_to_l) {
        //     std::cout << x.first << " (" << x.second.first 
        //         << ", " << x.second.second << ")" << std::endl; 
        // }   

        graph_t **graph_local_all = new graph_t* [size];
        for (int i=0; i<size; i++) graph_local_all[i] = new graph_t;
        
        std::vector<std::vector<idx_t>> adjncy_all(size);
        std::vector<std::vector<idx_t>> adjwgt_all(size);
        std::vector<std::vector<idx_t>> xadj_all(size);
        for (int i=0; i<size; i++) xadj_all[i].push_back(0);
        // std::cout << "1 :" << adjncy_all.size() << std::endl;
        // std::cout << "2 :" << adjwgt_all.size() << std::endl;

        // for (int i=0; i<10000; i++) {
        //     adjncy_all[0].push_back(i); 
        //     adjwgt_all[0].push_back(i); 
        //     std::cout << i << "??" << std::endl;
        // }
        for (int proc=0; proc<size; proc++)
            graph_local_all[proc]->nvtxs = local_num[proc];
        for (int i=0; i<graph->nvtxs; i++)
        {
            idx_t lrank = g_to_l[i].first;
            for (int edg_idx=graph->xadj[i]; 
                     edg_idx<graph->xadj[i+1];
                     edg_idx++)
            {
                idx_t edg = graph->adjncy[edg_idx] - starting_index;
                // std::cout << "Edge: " << edg << "local:(" <<g_to_l[edg].first << "," 
                //     << g_to_l[edg].second << ")" <<  std::endl;
                // std::cout << "edg_idx " << g_to_l[edg].second << " "  << lrank << std::endl;
                // std::cout << adjncy_all[lrank].size() << "!" <<  adjwgt_all[lrank].size() << "!" << std::endl;
                // std::cout <<  g_to_l[edg].first << "!!" << g_to_l[edg].second << "!!" << std::endl;
                adjncy_all[lrank].push_back(g_to_l[edg].second);
                // std::cout << "Here" << std::endl;
                adjwgt_all[lrank].push_back(g_to_l[edg].first);
            
            }
            // std::cout << "Here Outside : " << adjncy_all[lrank].size() << std::endl;
            // std::cout << "xadj_all" << xadj_all[lrank][0] << std::endl;
            xadj_all[lrank].push_back(adjncy_all[lrank].size());
        }
        for (int proc=0; proc<size; proc++)
            graph_local_all[proc]->nedges = adjncy_all[proc].size();
        for (int proc=0; proc<size; proc++) { 
            graph_t *now_graph = graph_local_all[proc];
            now_graph->adjncy = new idx_t [adjncy_all[proc].size()];
            std::copy(adjncy_all[proc].begin(), adjncy_all[proc].end(), now_graph->adjncy);
            now_graph->adjwgt = new idx_t [adjwgt_all[proc].size()];
            std::copy(adjwgt_all[proc].begin(), adjwgt_all[proc].end(), now_graph->adjwgt);
            now_graph->xadj = new idx_t [xadj_all[proc].size()];
            std::copy(xadj_all[proc].begin(), xadj_all[proc].end(), now_graph->xadj);
        }


        // Now finally we are transporting the graph
        for (int proc=0; proc<size; proc++){
            graph_t *now_graph = graph_local_all[proc];
            if (proc == 0) {
                memcpy(graph_local, now_graph, sizeof(graph));
                graph_local->xadj = new idx_t[graph_local->nvtxs];
                memcpy(graph_local->xadj, now_graph->xadj, sizeof(idx_t) * (1+graph_local->nvtxs));
            
                graph_local->adjncy = new idx_t[graph_local->nedges];
                memcpy(graph_local->adjncy, now_graph->adjncy, sizeof(idx_t) * graph_local->nedges);
                graph_local->adjwgt = new idx_t[graph_local->nedges];
                memcpy(graph_local->adjwgt, now_graph->adjwgt, sizeof(idx_t) * graph_local->nedges);
            } else {

                MPI_Datatype datatype;
                if (sizeof(idx_t) == 4) datatype = MPI_INT;
                else if (sizeof(idx_t) == 8) datatype = MPI_LONG;
                else datatype = MPI_INT;
                MPI_Send(now_graph, sizeof(graph_t), MPI_CHAR, proc, 49, MPI_COMM_WORLD);
                MPI_Send(now_graph->xadj, now_graph->nvtxs+1, datatype, proc, 50, MPI_COMM_WORLD);
                MPI_Send(now_graph->adjncy, now_graph->nedges, datatype, proc, 51, MPI_COMM_WORLD);
                MPI_Send(now_graph->adjwgt, now_graph->nedges, datatype, proc, 52, MPI_COMM_WORLD);
            }
           
        }
        // (TODO) clean up
        // deleteGraph(graph);
        // for (int proc=0; proc<size; proc++) 
        //     deleteGraph(graph_local_all[proc]);
        // delete [] graph_local_all;
        // delete [] part; 
        // delete params;
    }
    else {
        MPI_Status status;
        MPI_Recv(graph_local, sizeof(graph_t), MPI_CHAR, 0, 49, MPI_COMM_WORLD, &status);
        MPI_Datatype datatype;
        if (sizeof(idx_t) == 4) datatype = MPI_INT;
        else if (sizeof(idx_t) == 8) datatype = MPI_LONG;
        else datatype = MPI_INT;
        graph_local->xadj = new idx_t[graph_local->nvtxs];
        MPI_Recv(graph_local->xadj, graph_local->nvtxs+1, datatype, 0, 50, MPI_COMM_WORLD, &status);
        graph_local->adjncy = new idx_t[graph_local->nedges];
        graph_local->adjwgt = new idx_t[graph_local->nedges];
        MPI_Recv(graph_local->adjncy, graph_local->nedges, datatype, 0, 51, MPI_COMM_WORLD, &status);
        MPI_Recv(graph_local->adjwgt, graph_local->nedges, datatype, 0, 52, MPI_COMM_WORLD, &status);
    }
    // return tmpGetGraph();
    // std::cout << "Final" << " " << rank << std::endl;
    return static_cast<graph_t *>(graph_local);
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
    delete[] graph->vwgt;
    delete[] graph->adjncy;
    delete[] graph->adjwgt;
    delete[] graph->xadj;
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
                    printf("Here %d\n", k);
                
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


// Graph partitioning stuff

idx_t simpleReadGraph(std::string filename, graph_t *graph)
{
    std::ifstream infile(filename);
    idx_t NVERTS;
    idx_t NEDGES;
    std::string line;
    // First line
    {
        std::getline(infile, line);
        std::stringstream ss(line);
        ss >> NVERTS >> NEDGES;
        graph->nvtxs = NVERTS;
        graph->nedges = NEDGES;
    }

    std::vector<idx_t> adjncy;
    std::vector<idx_t> xadj;
    xadj.push_back(0);
    idx_t starting_idx = NVERTS; //  
    for (idx_t i=0; i<NVERTS; i++)
    {
        std::getline(infile, line);
        std::stringstream ss(line);
        idx_t temp;
        while (ss >>temp) {
            adjncy.push_back(temp);
            if (temp < starting_idx) starting_idx = temp;
        }
        xadj.push_back(adjncy.size());
    }
    graph->adjncy = new idx_t [adjncy.size()];
    std::copy(adjncy.begin(), adjncy.end(), graph->adjncy);
    graph->xadj = new idx_t [xadj.size()];
    std::copy(xadj.begin(), xadj.end(), graph->xadj);

    // for (int i=0; i<xadj.size(); i++)
    //     std::cout << graph->xadj[i] << " ";
    // std::cout << std::endl;

    // for (int i=0; i<adjncy.size(); i++)
    //     std::cout << graph->adjncy[i] << " ";
    // std::cout << std::endl;
    
    // graph->adjncy = new idx_t[NEDGES] {1, 2, 0, 2, 0, 1, 3, 2, 4, 5, 3, 5, 3, 4};
    // graph->xadj = new idx_t[NVERTS+1] {0, 2, 4, 7, 10, 12, 14};

    // Setting those unrelated stuff
    graph->vwgt = new idx_t[NVERTS];
    std::fill(graph->vwgt + 0, graph->vwgt + NVERTS, 1);
    graph->adjwgt = new idx_t[NEDGES];
    std::fill(graph->adjwgt + 0, graph->adjwgt + NVERTS, 1);
    
    // graph->adjwgt = new idx_t[NEDGES] {0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 2};
    graph->ncon = 1;
    return starting_idx;
}

void setParams(params_t *params, graph_t *graph, idx_t nparts)
{
  memset((void *)params, 0, sizeof(params_t));

  /* initialize the params data structure */
  params->ptype         = METIS_PTYPE_KWAY;
  params->objtype       = METIS_OBJTYPE_CUT;
  params->ctype         = METIS_CTYPE_SHEM;
  params->iptype        = -1;
  params->rtype         = -1;

  params->no2hop        = 0;
  params->minconn       = 0;
  params->contig        = 0;

  params->nooutput      = 1;
  params->wgtflag       = 3;

  params->ncuts         = 1;
  params->niter         = 10;

  params->dbglvl        = 0;
  params->balance       = 0;
  params->seed          = -1;
  params->dbglvl        = 0;

  params->tpwgtsfile    = NULL;

  params->filename      = NULL;

  params->ufactor       = -1;

  params->ubvecstr      = NULL;
  params->ubvec         = NULL;

  params->nparts        = nparts;
    params->tpwgts = new real_t [graph->ncon * params->nparts];
    std::fill(params->tpwgts + 0, 
        params->tpwgts + graph->ncon * params->nparts,
         1.0/ params->nparts);
    // params->ubvec = new real_t [graph->ncon]{1.05};

//     /* initialize the params data structure */
//   params->gtype     = METIS_GTYPE_DUAL;
//   params->ncommon   = 1;
//   params->dbglvl    = 0;
//   params->filename  = NULL;
//   params->outfile   = NULL;
  /* Set the ptype-specific defaults */
  if (params->ptype == METIS_PTYPE_RB) {
    params->rtype   = METIS_RTYPE_FM;
  }
  if (params->ptype == METIS_PTYPE_KWAY) {
    params->iptype  = METIS_IPTYPE_METISRB;
    params->rtype   = METIS_RTYPE_GREEDY;
  }

         
  /* Setup iptype */
  if (params->iptype == -1) {
    if (params->ptype == METIS_PTYPE_RB) {
      if (graph->ncon == 1)
        params->iptype = METIS_IPTYPE_GROW;
      else
        params->iptype = METIS_IPTYPE_RANDOM;
    }
  } 
}

void prepareForPartition(idx_t options[], params_t *params, int starting_index)
{
    METIS_SetDefaultOptions(options);
    options[METIS_OPTION_OBJTYPE] = params->objtype;
    options[METIS_OPTION_CTYPE]   = params->ctype;
    options[METIS_OPTION_IPTYPE]  = params->iptype;
    options[METIS_OPTION_RTYPE]   = params->rtype;
    options[METIS_OPTION_NO2HOP]  = params->no2hop;
    options[METIS_OPTION_MINCONN] = params->minconn;
    options[METIS_OPTION_CONTIG]  = params->contig;
    options[METIS_OPTION_SEED]    = params->seed;
    options[METIS_OPTION_NITER]   = params->niter;
    options[METIS_OPTION_NCUTS]   = params->ncuts;
    options[METIS_OPTION_UFACTOR] = params->ufactor;
    options[METIS_OPTION_DBGLVL]  = params->dbglvl;
    options[METIS_OPTION_NUMBERING]  = starting_index;
}

void printGraph(graph_t *graph){
    std::cout << "N=" << graph->nvtxs << " m=" << graph->nedges << std::endl; 
    std::cout << "xadj:";
    for (int n=0; n<graph->nvtxs+1; n++) 
        std::cout << graph->xadj[n] << " ";
    std::cout << std::endl;

    std::cout << "adjncy:";
    for (int n=0; n<graph->nedges; n++) 
        std::cout << graph->adjncy[n] << " ";
    std::cout << std::endl;

    std::cout << "adjwgt:";
    for (int n=0; n<graph->nedges; n++) 
        std::cout << graph->adjwgt[n] << " ";
    std::cout << std::endl;
    
    std::cout << "vwgt:";
    for (int n=0; n<graph->nvtxs; n++) 
        std::cout << graph->vwgt[n] << " ";
    std::cout << std::endl;
    
}




/*************************************************************************/
/*! This function prints run parameters */
/*************************************************************************/
void GPPrintInfo(params_t *params, graph_t *graph)
{ 
  idx_t i;

  if (params->ufactor == -1) {
    if (params->ptype == METIS_PTYPE_KWAY)
      params->ufactor = KMETIS_DEFAULT_UFACTOR;
    else if (graph->ncon == 1)
      params->ufactor = PMETIS_DEFAULT_UFACTOR;
    else
      params->ufactor = MCPMETIS_DEFAULT_UFACTOR;
  }

  printf("******************************************************************************\n");
  printf("%s", METISTITLE);
  printf(" (HEAD: %s, Built on: %s, %s)\n", "", __DATE__, __TIME__);
  printf(" size of idx_t: %zubits, real_t: %zubits, idx_t *: %zubits\n", 
      8*sizeof(idx_t), 8*sizeof(real_t), 8*sizeof(idx_t *));
  printf("\n");
  printf("Graph Information -----------------------------------------------------------\n");
  printf(" Name: %s, #Vertices: %"PRIDX", #Edges: %"PRIDX", #Parts: %"PRIDX"\n", 
      params->filename, graph->nvtxs, graph->nedges/2, params->nparts);
  if (graph->ncon > 1)
    printf(" Balancing constraints: %"PRIDX"\n", graph->ncon);

  printf("\n");
  printf("Options ---------------------------------------------------------------------\n");
  printf(" ptype=%s, objtype=%s, ctype=%s, rtype=%s, iptype=%s\n",
      ptypenames[params->ptype], objtypenames[params->objtype], ctypenames[params->ctype], 
      rtypenames[params->rtype], iptypenames[params->iptype]);

  printf(" dbglvl=%"PRIDX", ufactor=%.3f, no2hop=%s, minconn=%s, contig=%s, nooutput=%s\n",
      params->dbglvl,
      I2RUBFACTOR(params->ufactor),
      (params->no2hop   ? "YES" : "NO"), 
      (params->minconn  ? "YES" : "NO"), 
      (params->contig   ? "YES" : "NO"),
      (params->nooutput ? "YES" : "NO")
      );

  printf(" seed=%"PRIDX", niter=%"PRIDX", ncuts=%"PRIDX"\n", 
      params->seed, params->niter, params->ncuts);

  if (params->ubvec) {
    printf(" ubvec=(");
    for (i=0; i<graph->ncon; i++)
      printf("%s%.2e", (i==0?"":" "), (double)params->ubvec[i]);
    printf(")\n");
  }

  printf("\n");
  switch (params->ptype) {
    case METIS_PTYPE_RB:
      printf("Recursive Partitioning ------------------------------------------------------\n");
      break;
    case METIS_PTYPE_KWAY:
      printf("Direct k-way Partitioning ---------------------------------------------------\n");
      break;
  }
}
