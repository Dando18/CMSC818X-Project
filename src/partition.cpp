#include <iostream>
#include <fstream>
#include <sstream>
#include "metisbin.h"
#include <vector>
#include <map>
#include "mpi.h"
// #include <metis.h>
// #include <metislib.h>

void simpleReadGraph(std::string filename, graph_t *graph)
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
    for (idx_t i=0; i<NVERTS; i++)
    {
        std::getline(infile, line);
        std::stringstream ss(line);
        idx_t temp;
        while (ss >>temp) adjncy.push_back(temp);
        xadj.push_back(adjncy.size());
    }
    graph->adjncy = new idx_t [adjncy.size()];
    std::copy(adjncy.begin(), adjncy.end(), graph->adjncy);
    graph->xadj = new idx_t [xadj.size()];
    std::copy(xadj.begin(), xadj.end(), graph->xadj);

    for (int i=0; i<xadj.size(); i++)
        std::cout << graph->xadj[i] << " ";
    std::cout << std::endl;

    for (int i=0; i<adjncy.size(); i++)
        std::cout << graph->adjncy[i] << " ";
    std::cout << std::endl;
    
    // graph->adjncy = new idx_t[NEDGES] {1, 2, 0, 2, 0, 1, 3, 2, 4, 5, 3, 5, 3, 4};
    // graph->xadj = new idx_t[NVERTS+1] {0, 2, 4, 7, 10, 12, 14};

    // Setting those unrelated stuff
    graph->vwgt = new idx_t[NVERTS];
    std::fill(graph->vwgt + 0, graph->vwgt + NVERTS, 1);
    graph->adjwgt = new idx_t[NEDGES];
    std::fill(graph->adjwgt + 0, graph->adjwgt + NVERTS, 1);
    
    // graph->adjwgt = new idx_t[NEDGES] {0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 2};
    graph->ncon = 1;

}

void setParams(params_t *params, graph_t *graph)
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

  params->nparts = 2;
  params->ufactor       = -1;

  params->ubvecstr      = NULL;
  params->ubvec         = NULL;

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

void prepareForPartition(idx_t options[], params_t *params, int start_index=0)
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
    options[METIS_OPTION_NUMBERING]  = start_index;
}

void printGraph(graph_t *graph){
    std::cout << "N=" << graph->nvtxs << " m=" << graph->nedges << std::endl; 
    std::cout << "xadj:";
    for (int n=0; n<graph->nvtxs; n++) 
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
}

void cleanGraph(graph_t *graph){
    // (TODO)
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

int main(int argc, char *argv[])
{
    int ierr;
    ierr = MPI_Init(&argc, &argv);
    
    int nprocs, rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    graph_t *graph_local = new graph_t;
    graph_t *graph;
    if (rank == 0){
        graph = new graph_t;
        params_t *params = new params_t ;

        simpleReadGraph(argv[1], graph);

        setParams(params, graph);
        GPPrintInfo(params, graph);

        int start_index = 1;
        idx_t options[METIS_NOPTIONS];
        prepareForPartition(options, params, start_index);
        params->nparts = nprocs;

        idx_t objval;
        idx_t *part = new idx_t[graph->nvtxs];
        int status = 0;
        status = METIS_PartGraphKway(&graph->nvtxs, &graph->ncon, graph->xadj, 
                graph->adjncy, graph->vwgt, graph->vsize, graph->adjwgt, 
                &params->nparts, params->tpwgts, params->ubvec, options, 
                &objval, part);
        for (idx_t i=0; i<graph->nvtxs; i++){
            part[i] -= start_index;
            // std::cout << "Node:" << i << "\tPart:" << part[i] << std::endl; 
        }
    
        // broadcast
        // set up a mapping from global_idx to (rank, local_idx)
        std::map<idx_t, std::pair<idx_t, idx_t> > g_to_l;
        std::vector<idx_t> local_num(nprocs, 0);
        for (idx_t i=0; i<graph->nvtxs; i++) {
            g_to_l[i] = std::make_pair(part[i], local_num[part[i]]++);
        }
        // for (auto x:g_to_l) {
        //     std::cout << x.first << " (" << x.second.first 
        //         << ", " << x.second.second << ")" << std::endl; 
        // }   

        graph_t **graph_local_all = new graph_t* [nprocs];
        std::vector<std::vector<idx_t>> adjncy_all(nprocs);
        std::vector<std::vector<idx_t>> adjwgt_all(nprocs);
        std::vector<std::vector<idx_t>> xadj_all(nprocs, std::vector<idx_t>(1,0));

        for (int proc=0; proc<nprocs; proc++)
            graph_local_all[proc]->nvtxs = local_num[proc];
        for (int i=0; i<graph->nvtxs; i++)
        {
            idx_t lrank = g_to_l[i].first;
            for (int edg_idx=graph->xadj[i]; 
                     edg_idx<graph->xadj[i+1];
                     edg_idx++)
            {
                idx_t edg = graph->adjncy[edg_idx] - start_index;
                // std::cout << "Edge: " << edg << "local:(" <<g_to_l[edg].first << "," 
                //     << g_to_l[edg].second << ")" <<  std::endl;
                adjncy_all[lrank].push_back(g_to_l[edg].second);
                adjwgt_all[lrank].push_back(g_to_l[edg].first);
            }
            xadj_all[lrank].push_back(adjncy_all[lrank].size());
        }
        for (int proc=0; proc<nprocs; proc++)
            graph_local_all[proc]->nedges = adjncy_all[proc].size();

        for (int proc=0; proc<nprocs; proc++) { 
            graph_t *now_graph = graph_local_all[proc];
            now_graph->adjncy = new idx_t [adjncy_all[proc].size()];
            std::copy(adjncy_all[proc].begin(), adjncy_all[proc].end(), now_graph->adjncy);
            now_graph->adjwgt = new idx_t [adjwgt_all[proc].size()];
            std::copy(adjwgt_all[proc].begin(), adjwgt_all[proc].end(), now_graph->adjwgt);
            now_graph->xadj = new idx_t [xadj_all[proc].size()];
            std::copy(xadj_all[proc].begin(), xadj_all[proc].end(), now_graph->xadj);
        }

        for (int proc=0; proc<nprocs; proc++) {
            std::cout << "=========================" << std::endl; 
            std::cout << "Print Graph in process " << proc << std::endl;
            printGraph(graph_local_all[proc]);
            std::cout << "=========================" << std::endl; 
        }

        // Now finally we are transporting the graph
        for (int proc=0; proc<nprocs; proc++){
            graph_t *now_graph = graph_local_all[proc];
            if (proc == 0) {
                memcpy(graph_local, now_graph, sizeof(graph));
                graph_local->xadj = new idx_t[graph_local->nvtxs];
                memcpy(graph_local->xadj, now_graph->xadj, sizeof(idx_t) * graph_local->nvtxs);
            
                graph_local->adjncy = new idx_t[graph_local->nedges];
                memcpy(graph_local->adjncy, now_graph->adjncy, sizeof(idx_t) * graph_local->nedges);
                graph_local->adjwgt = new idx_t[graph_local->nedges];
                memcpy(graph_local->adjwgt, now_graph->adjwgt, sizeof(idx_t) * graph_local->nedges);
            } else {
                MPI_Datatype datatype;
                if (sizeof(idx_t) == 4) datatype = MPI_INT;
                else if (sizeof(idx_t) == 8) datatype = MPI_LONG;
                else datatype = MPI_INT;
                MPI_Send(now_graph, sizeof(graph_t), MPI_CHAR, 0, 49, MPI_COMM_WORLD);
                MPI_Send(now_graph->xadj, now_graph->nvtxs, datatype, 0, 50, MPI_COMM_WORLD);
                MPI_Send(now_graph->adjncy, now_graph->nedges, datatype, 0, 51, MPI_COMM_WORLD);
                MPI_Send(now_graph->adjwgt, now_graph->nedges, datatype, 0, 52, MPI_COMM_WORLD);
            }

            
        }

        // clean upp
        cleanGraph(graph);
        for (int proc=0; proc<nprocs; proc++) 
            cleanGraph(graph_local_all[proc]);
        delete [] part; 
        delete params;
    }
    else {
        MPI_Recv(graph_local, sizeof(graph_t), MPI_CHAR, 0, 49, MPI_COMM_WORLD, nullptr);
        MPI_Datatype datatype;
        if (sizeof(idx_t) == 4) datatype = MPI_INT;
        else if (sizeof(idx_t) == 8) datatype = MPI_LONG;
        else datatype = MPI_INT;
        graph_local->xadj = new idx_t[graph_local->nvtxs];
        MPI_Recv(graph_local->xadj, graph_local->nvtxs, datatype, 0, 50, MPI_COMM_WORLD, nullptr);
        graph_local->adjncy = new idx_t[graph_local->nedges];
        graph_local->adjwgt = new idx_t[graph_local->nedges];
        MPI_Recv(graph_local->adjncy, graph_local->nedges, datatype, 0, 51, MPI_COMM_WORLD, nullptr);
        MPI_Recv(graph_local->adjwgt, graph_local->nedges, datatype, 0, 52, MPI_COMM_WORLD, nullptr);
    }

    MPI_Finalize();
    return 0;
}

