#include "metis.h"
#include "metisbin.h"

int main(int argc, char* argv[])
{
    graph_t *graph;
    idx_t options[METIS_NOPTIONS];
    params_t *params;
    char *curptr, *newptr;
    int i;
    params = parse_cmdline(argc, argv);
    
    gk_startcputimer(params->iotimer);
    graph = ReadGraph(params);

    ReadTPwgts(params, graph->ncon);
    gk_stopcputimer(params->iotimer);
    
    /* Check if the graph is contiguous */
    if (params->contig && !IsConnected(graph, 0)) {
      printf("***The input graph is not contiguous.\n"
             "***The specified -contig option will be ignored.\n");
      params->contig = 0;
    }
    /* Get ubvec if supplied */
    if (params->ubvecstr) {
      params->ubvec = rmalloc(graph->ncon, "main");
      curptr = params->ubvecstr;
      for (i=0; i<graph->ncon; i++) {
        params->ubvec[i] = strtoreal(curptr, &newptr);
        if (curptr == newptr)
          errexit("Error parsing entry #%"PRIDX" of ubvec [%s] (possibly missing).\n",
              i, params->ubvecstr);
        curptr = newptr;
      }
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
    GPPrintInfo(params, graph);

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

