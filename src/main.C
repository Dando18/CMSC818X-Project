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
#include <cmath>
#include <cassert>

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
graph_t *getGraph(std::string const& filename, int seed);
graph_t *tmpGetGraph();
graph_t *tmpGetBipartieGraph();
graph_t *getRandomGraph(size_t vertsPerProcess, size_t maxDegree, size_t numBoundaries);
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
        std::cout << "Usage: \t" << argv[0]  
                  << " [graph filename] [chunkSize] [seed]\n" 
                  << "or if filename='-'\t" << argv[0] << " - [chunkSize] [vertsPerProcess] [maxDegree] [connectivity]\n" << std::endl;
    }

    /* primitive MPI info */
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    /* read in args */
    int seed = -1;
    std::string fileName = std::string(argv[1]);
    std::size_t chunkSize = std::stoul(argv[2]);
    if (argc == 4) {
        seed = std::stoi(argv[3]);
    }

    /* Read in graph structure and partition */
    graph_t *graph;
    if (fileName != "-") {
        graph = getGraph(fileName, seed);
        
    } else {
        std::size_t vertsPerProcess = (argc >= 4) ? std::stoul(argv[3]) : 1000;
        std::size_t maxDegree = (argc >= 5) ? std::stoul(argv[4]) : 50;
        double connectivity = (argc >= 6) ? std::stod(argv[5]) : 0.1;
        assert( connectivity > 0.0 && connectivity <= 1.0 );
        graph = getRandomGraph(vertsPerProcess, maxDegree, static_cast<size_t>(vertsPerProcess*size*connectivity*0.9));
    }
    graph->vwgt = new idx_t[graph->nvtxs];
    std::fill(graph->vwgt + 0, graph->vwgt + graph->nvtxs, -1);
    
    MPI_Barrier(MPI_COMM_WORLD);

    /* color graph in parallel */
    coloring_stats_t stats;
    double start = MPI_Wtime();
    if (size == 1) {
        std::vector<int> U (graph->nvtxs);
        std::iota(U.begin(), U.end(), 0);
        colorGraphSerial(U, graph, {});
        stats.numIterations = 1;
    } else {
        colorGraph(graph, chunkSize, stats);
    }
    double end = MPI_Wtime();

    /* output coloring to file or stdout */ 
    //outputColoring(graph, std::cout);

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
graph_t *getGraph(std::string const& filename, int seed) {

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

        if (size == 1) {
            delete graph_local;
            delete params;

            std::fill(graph->adjwgt + 0, graph->adjwgt + graph->nedges, 0);
            // re-index adjncy
            for (int i = 0; i < graph->nedges; i++) {
                graph->adjncy[i] -= starting_index;
            }
            return graph;
        }
        
        setParams(params, graph, size);
        params->seed = seed;
        GPPrintInfo(params, graph);
        idx_t options[METIS_NOPTIONS];
        prepareForPartition(options, params, starting_index);

        idx_t objval;
        idx_t *part = new idx_t[graph->nvtxs];
        int status = 0;

        status = METIS_PartGraphKway(&graph->nvtxs, &graph->ncon, graph->xadj, 
                graph->adjncy, graph->vwgt, graph->vsize, graph->adjwgt, 
                &params->nparts, params->tpwgts, params->ubvec, options, 
                &objval, part);
        
        switch (status) {
            case METIS_OK: std::cout << "METIS_OK" << std::endl; break;
            case METIS_ERROR_INPUT: std::cout << "METIS_ERROR_INPUT" << std::endl; break;
            case METIS_ERROR_MEMORY: std::cout << "METIS_ERROR_MEMORY" << std::endl; break;
            case METIS_ERROR: std::cout << "METIS_ERROR" << std::endl; break;
        }

        std::cout << "params->nparts = " << params->nparts << std::endl;

        std::set<idx_t> uniqueParts;

        for (idx_t i=0; i<graph->nvtxs; i++){
            uniqueParts.insert(part[i]);
            part[i] -= starting_index;
        }

        std::cout << "num unique parts = " << uniqueParts.size() << std::endl;
        if (uniqueParts.size() != size) {
            std::cerr << "For some odd reason 'METIS_PartGraphKway' did not generate the correct number of graph partitions." << std::endl;
            MPI_Abort(MPI_COMM_WORLD, 1);
        }

        // broadcast
        // set up a mapping from global_idx to (rank, local_idx)
        std::map<idx_t, std::pair<idx_t, idx_t> > g_to_l;
        std::vector<idx_t> local_num(size, 0);
        for (idx_t i=0; i<graph->nvtxs; i++) {
            g_to_l.insert(std::make_pair(i, std::make_pair(part[i], local_num.at(part[i])++)));
        }

        graph_t **graph_local_all = new graph_t* [size];
        for (int i=0; i<size; i++) graph_local_all[i] = new graph_t;
        
        std::vector<std::vector<idx_t>> adjncy_all(size);
        std::vector<std::vector<idx_t>> adjwgt_all(size);
        std::vector<std::vector<idx_t>> xadj_all(size);
        for (int i=0; i<size; i++) xadj_all.at(i).push_back(0);
        for (int proc=0; proc<size; proc++)
            graph_local_all[proc]->nvtxs = local_num[proc];
        for (int i=0; i<graph->nvtxs; i++)
        {
            idx_t lrank = g_to_l.at(i).first;
            for (int edg_idx=graph->xadj[i]; 
                     edg_idx<graph->xadj[i+1];
                     edg_idx++)
            {
                idx_t edg = graph->adjncy[edg_idx] - starting_index;
                adjncy_all[lrank].push_back(g_to_l.at(edg).second);
                adjwgt_all[lrank].push_back(g_to_l.at(edg).first);
            
            }
            xadj_all.at(lrank).push_back(adjncy_all.at(lrank).size());
        }
        for (int proc=0; proc<size; proc++)
            graph_local_all[proc]->nedges = adjncy_all.at(proc).size();
        
        idx_t runningTotal = 0;
        for (int proc=0; proc<size; proc++) { 
            graph_t *now_graph = graph_local_all[proc];
            now_graph->adjncy = new idx_t [adjncy_all.at(proc).size()];
            std::copy(adjncy_all.at(proc).begin(), adjncy_all.at(proc).end(), now_graph->adjncy);
            now_graph->adjwgt = new idx_t [adjwgt_all.at(proc).size()];
            std::copy(adjwgt_all.at(proc).begin(), adjwgt_all.at(proc).end(), now_graph->adjwgt);
            now_graph->xadj = new idx_t [xadj_all.at(proc).size()];
            std::copy(xadj_all.at(proc).begin(), xadj_all.at(proc).end(), now_graph->xadj);

            now_graph->ncon = runningTotal;
            runningTotal += now_graph->nvtxs;
        }


        // Now finally we are transporting the graph
        for (int proc=0; proc<size; proc++){
            graph_t *now_graph = graph_local_all[proc];
            if (proc == 0) {
                memcpy(graph_local, now_graph, sizeof(graph));
                graph_local->xadj = new idx_t[graph_local->nvtxs+1];
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
                //MPI_Send(now_graph, sizeof(graph_t), MPI_CHAR, proc, 49, MPI_COMM_WORLD);
                MPI_Send(&(now_graph->nvtxs), 1, datatype, proc, 49, MPI_COMM_WORLD);
                MPI_Send(&(now_graph->nedges), 1, datatype, proc, 50, MPI_COMM_WORLD);
                MPI_Send(&(now_graph->ncon), 1, datatype, proc, 51, MPI_COMM_WORLD);
                MPI_Send(now_graph->xadj, now_graph->nvtxs+1, datatype, proc, 52, MPI_COMM_WORLD);
                MPI_Send(now_graph->adjncy, now_graph->nedges, datatype, proc, 53, MPI_COMM_WORLD);
                MPI_Send(now_graph->adjwgt, now_graph->nedges, datatype, proc, 54, MPI_COMM_WORLD);
                
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

        MPI_Datatype datatype;
        if (sizeof(idx_t) == 4) datatype = MPI_INT;
        else if (sizeof(idx_t) == 8) datatype = MPI_LONG;
        else datatype = MPI_INT;
        MPI_Recv(&(graph_local->nvtxs), 1, datatype, 0, 49, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(&(graph_local->nedges), 1, datatype, 0, 50, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(&(graph_local->ncon), 1, datatype, 0, 51, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        graph_local->xadj = new idx_t[graph_local->nvtxs+1];
        MPI_Recv(graph_local->xadj, graph_local->nvtxs+1, datatype, 0, 52, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        graph_local->adjncy = new idx_t[graph_local->nedges];
        graph_local->adjwgt = new idx_t[graph_local->nedges];
        MPI_Recv(graph_local->adjncy, graph_local->nedges, datatype, 0, 53, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(graph_local->adjwgt, graph_local->nedges, datatype, 0, 54, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
    return static_cast<graph_t *>(graph_local);
}


/**
 * @param[in] vertsPerProcess the number of vertices to put on each process
 * @param[in] numBoundaries the number of external edges per rank ~ ish
 */
graph_t *getRandomGraph(size_t vertsPerProcess, size_t maxDegree, size_t numBoundaries) {
    int rank, size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int *buf = new int[numBoundaries*4];

    if (rank == 0) {
        /* Build table */
        std::set<std::tuple<int, int, int, int>> unique;
        while (unique.size() < numBoundaries) {
            int rank1 = rand() % size;
            int rank2;
            do { rank2 = rand() % size; } while (rank2 == rank1 && size != 1);
            int vertex1 = rand() % vertsPerProcess, vertex2 = rand() % vertsPerProcess;

            unique.insert(std::make_tuple(rank1, vertex1, rank2, vertex2));
        }

        int i = 0;
        for (std::tuple<int, int, int, int> const& vals : unique) {
            int r1, v1, r2, v2;
            std::tie(r1, v1, r2, v2) = vals;
            buf[i]=r1; buf[i+1]=v1; buf[i+2]=r2; buf[i+3]=v2;

            i += 4;
        }
    }

    /* distributed table */
    MPI_Bcast(buf, numBoundaries*4, MPI_INT, 0, MPI_COMM_WORLD);

    /* build graph locally on each rank */
    graph_t *graph = new graph_t;
    graph->nvtxs = vertsPerProcess;
    //graph->vwgt = new idx_t[vertsPerProcess];
    
    std::vector<std::vector<idx_t>> adjList(graph->nvtxs);
    std::vector<std::vector<idx_t>> adjListRanks(graph->nvtxs);
    for (size_t i = 0; i < numBoundaries*4; i += 4) {
        int rank1 = buf[i+0], v1 = buf[i+1], rank2 = buf[i+2], v2 = buf[i+3];

        if (rank1 == rank) {
            adjList.at(v1).push_back(v2);
            adjListRanks.at(v1).push_back(rank2);
        } else if (rank2 == rank) {
            adjList.at(v2).push_back(v1);
            adjListRanks.at(v2).push_back(rank1);
        }
    }

    graph->xadj = new idx_t[vertsPerProcess + 1];
    graph->xadj[0] = 0;
    for (int i = 0; i < vertsPerProcess; i++) {
        int n = adjList.at(i).size() + rand() % (maxDegree - adjList.at(i).size());
        graph->xadj[i + 1] = graph->xadj[i] + n;
    }

    for (size_t i = 0; i < adjList.size(); i++) {
        while (adjList.at(i).size() < graph->xadj[i + 1] - graph->xadj[i]) {
            int vertex;
            do { vertex = rand() % vertsPerProcess; } while (vertex == i);

            adjList.at(i).push_back(vertex);
            adjListRanks.at(i).push_back(rank);
        }
    }

    

    graph->nedges = graph->xadj[vertsPerProcess];
    graph->adjncy = new idx_t[graph->nedges];
    graph->adjwgt = new idx_t[graph->nedges];
    for (int i = 0; i < vertsPerProcess; i++) {
        for (int j = graph->xadj[i]; j < graph->xadj[i + 1]; j++) {
            graph->adjncy[j] = adjList.at(i).at(j - graph->xadj[i]);
            graph->adjwgt[j] = adjListRanks.at(i).at(j - graph->xadj[i]);
        }
    }

    delete[] buf;
    return graph;
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

graph_t *tmpGetBipartieGraph() {
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (size != 4) {
        MPI_Abort(MPI_COMM_WORLD, 0);
    }

    graph_t *graph = new graph_t;



    if (rank == 0) {
        constexpr idx_t NVERTS = 2;
        graph->nvtxs = NVERTS;
        graph->vwgt = new idx_t[NVERTS];
        graph->ncon = 0;
        std::fill(graph->vwgt + 0, graph->vwgt + NVERTS, -1);
        constexpr idx_t NEDGES = 3;
        graph->nedges = NEDGES;
        graph->adjncy = new idx_t[NEDGES] {0, 1, 0};
        graph->adjwgt = new idx_t[NEDGES] {1, 3, 1};
        graph->xadj = new idx_t[NVERTS+1] {0, 2, 3};
        

    } else if (rank == 1) {
        constexpr idx_t NVERTS = 3;
        graph->nvtxs = NVERTS;
        graph->vwgt = new idx_t[NVERTS];
        graph->ncon = 2;
        std::fill(graph->vwgt + 0, graph->vwgt + NVERTS, -1);
        constexpr idx_t NEDGES = 12;
        graph->nedges = NEDGES;
        graph->adjncy = new idx_t[NEDGES] {0, 1, 2, 0, 1, 0, 1, 0, 2, 0, 1, 0};
        graph->adjwgt = new idx_t[NEDGES] {3, 1, 1, 0, 0, 1, 3, 2, 3, 1, 3, 2};
        graph->xadj = new idx_t[NVERTS+1] {0, 5, 9, 12};

    } else if (rank == 2) {
        constexpr idx_t NVERTS = 2;
        graph->nvtxs = NVERTS;
        graph->vwgt = new idx_t[NVERTS];
        graph->ncon = 5;
        std::fill(graph->vwgt + 0, graph->vwgt + NVERTS, -1);
        constexpr idx_t NEDGES = 4;
        graph->nedges = NEDGES;
        graph->adjncy = new idx_t[NEDGES] {0, 1, 2, 0};
        graph->adjwgt = new idx_t[NEDGES] {3, 1, 1, 3};
        graph->xadj = new idx_t[NVERTS+1] {0, 3, 4};

    } else if (rank == 3) {
        constexpr idx_t NVERTS = 3;
        graph->nvtxs = NVERTS;
        graph->vwgt = new idx_t[NVERTS];
        graph->ncon = 7;
        std::fill(graph->vwgt + 0, graph->vwgt + NVERTS, -1);
        constexpr idx_t NEDGES = 11;
        graph->nedges = NEDGES;
        graph->adjncy = new idx_t[NEDGES] {0, 1, 0, 2, 1, 0, 1, 2, 0, 0, 1};
        graph->adjwgt = new idx_t[NEDGES] {1, 3, 2, 3, 2, 3, 1, 1, 0, 3, 1};
        graph->xadj = new idx_t[NVERTS+1] {0, 5, 9, 11};

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
    std::vector<idx_t> graphSizes(size+1, 0);
    std::iota(U.begin(), U.end(), 0);

    MPI_Allgather(&(graph->ncon), 1, MPI_INT, graphSizes.data(), 1, MPI_INT, MPI_COMM_WORLD);

    int size_U = U.size();
    int max_size_U = 1;
    int loop = 0;
    while (max_size_U != 0) {   // (5)

        if (size_U != 0) {      // (6)
            // Partition U into subsets
            std::vector<std::vector<idx_t>> Us = partitionU(U, s);

            // Communicate with other processors
            // 1. find all boundary vertices
            std::vector<MPI_Request> recvRequests, sendRequests;
            std::vector<int> recvRanks;
            int num_external_edges = 0;
            
            // count all of graph->adjwgt != rank
            for (idx_t i = 0; i < graph->nedges; i++) {
                num_external_edges += (graph->adjwgt[i] != rank) ? 1 : 0;
            }
            recvRequests.reserve(num_external_edges);
            sendRequests.reserve(num_external_edges);
            recvRanks.reserve(num_external_edges);

            idx_t *recv_buf = new idx_t[3 * num_external_edges];
            std::fill(recv_buf + 0, recv_buf + (3*num_external_edges), -1);
            idx_t *send_buf = new idx_t[3 * num_external_edges];
            int counter = 0;

            // Supersteps -- color each subset of U
            for (size_t k = 0; k < Us.size(); k++) {     // (8)
                std::copy(graph->vwgt + 0, graph->vwgt + graph->nvtxs, previousColors.begin());
                colorGraphSerial(Us.at(k), graph, boundaryColors); // (9-10)
                
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
                            sendRequests.push_back(sendRequest);
                            
                            MPI_Irecv(&(recv_buf[counter]), 3, MPI_INT, neighborRank, 0, MPI_COMM_WORLD, &recvRequest);
                            recvRequests.push_back(recvRequest);
                            recvRanks.push_back(neighborRank);
                            counter += 3;
                        }
                        
                    }
                }
            }   // end of superstep
            
            MPI_Waitall(sendRequests.size(), sendRequests.data(), MPI_STATUS_IGNORE);
            MPI_Waitall(recvRequests.size(), recvRequests.data(), MPI_STATUS_IGNORE);

            for (size_t i = 0; i < 3 * recvRequests.size(); i = i + 3) {
                idx_t neighborIdx = recv_buf[i];
                idx_t neighborColor = recv_buf[i + 1];
                idx_t myIdx = recv_buf[i + 2];
                idx_t neighborRank = recvRanks.at(i / 3);

                for (int j = graph->xadj[myIdx]; j < graph->xadj[myIdx+1]; j++) {
                    if (neighborIdx == graph->adjncy[j] && neighborRank == graph->adjwgt[j]) {
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
                        int neighborValue = graphSizes.at(neighborRank) + neighborIdx;
                        int myValue = graphSizes.at(rank) + i;

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
            // U = R
            U.assign(R.begin(), R.end());
        }   // end of "if U is not empty"

        size_U = U.size();
        MPI_Allreduce(&size_U, &max_size_U, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);

        loop++;
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
            if (neighborRank == rank) {
                if (graph->vwgt[neighbor] != -1) { 
                    neighbor_colors.at(graph->vwgt[neighbor]) = true; // this color is used
                }
            } else if (boundaryColors.at(j) != -1) {
                idx_t neighborColor = boundaryColors.at(j);
                neighbor_colors.at(neighborColor) = true; // this color is used
            }
        }
        for (std::size_t j = 0; j < MAX_COLOR; j++) {
            if (!neighbor_colors.at(j)) { // first color not used
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
    std::vector<std::vector<idx_t>> Us;// (ceil(U.size() / s));
    //int rank;
    //MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    for (std::size_t i = 0; i < U.size(); i += s) {
        auto last = std::min(U.size(), i + s);
        //Us.emplace_back(U.begin() + i, U.begin() + last);
        std::vector<idx_t> v;
        for (std::size_t j = i; j < last; j++) {
            v.push_back(U.at(j));
        }
        //std::vector<idx_t> v (U.begin() + i, U.begin() + last);
        Us.push_back(v);
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

    // Double the number of edges
    graph->nedges *= 2;

    // Setting those unrelated stuff
    graph->vwgt = new idx_t[NVERTS];
    std::fill(graph->vwgt + 0, graph->vwgt + NVERTS, 1);
    graph->adjwgt = new idx_t[graph->nedges];
    std::fill(graph->adjwgt + 0, graph->adjwgt + graph->nedges, 1);

    graph->vsize = new idx_t[NVERTS];
    
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
  params->seed          = 1;
  params->dbglvl        = 0;

  params->tpwgtsfile    = NULL;

  params->filename      = NULL;

  params->ufactor       = -1;

  params->ubvecstr      = NULL;
  params->ubvec         = NULL;

  params->nparts        = nparts;
  params->tpwgts = NULL;
    //params->tpwgts = new real_t [graph->ncon * params->nparts];
    //std::fill(params->tpwgts + 0, 
    //    params->tpwgts + graph->ncon * params->nparts,
    //     1.0/ params->nparts);
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
