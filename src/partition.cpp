
// #include <metis.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include "metisbin.h"
#include <vector>
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
    graph_t *graph = new graph_t;
    params_t *params = new params_t ;
    idx_t options[METIS_NOPTIONS];
    idx_t *part;
    idx_t objval;

    int status = 0;
    simpleReadGraph(argv[1], graph);


    part = new idx_t[graph->nvtxs];

    setParams(params, graph);

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
    options[METIS_OPTION_NUMBERING ]  = 1;
    GPPrintInfo(params, graph);

//   switch (params->ptype) {
    // case METIS_PTYPE_RB:
    // status = METIS_PartGraphRecursive(&graph->nvtxs, &graph->ncon, graph->xadj, 
    //             graph->adjncy, graph->vwgt, graph->vsize, graph->adjwgt, 
    //             &params->nparts, params->tpwgts, params->ubvec, options, 
    //             &objval, part);
    //   break;

    // case METIS_PTYPE_KWAY:
      status = METIS_PartGraphKway(&graph->nvtxs, &graph->ncon, graph->xadj, 
                   graph->adjncy, graph->vwgt, graph->vsize, graph->adjwgt, 
                   &params->nparts, params->tpwgts, params->ubvec, options, 
                   &objval, part);
    //   break;
    for (idx_t i=0; i<graph->nvtxs; i++){
        std::cout << "Node:" << i << "\tPart:" << part[i] << std::endl; 
    }
    delete [] part; 
    // clean_graph(graph);
    return 0;
  }

