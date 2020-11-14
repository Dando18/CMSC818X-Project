#include "metis.h"
#include "metisbin.h"
#include <iostream>
#include <string>

graph_t *getGraph(std::string);

int main(int argc, char *argv[])
{
    /* read in args */
    std::string filename = std::string(argv[1]);
    std::size_t chunkSize = std::stoul(argv[2]);
    graph_t *graph;
    graph = getGraph(filename);
}

graph_t *getGraph(std::string filename)
{
    graph_t *graph;
    params_t *params;
    params = (params_t*)malloc(sizeof(params_t));
    memset((void*)params, 0, sizeof(params_t));
    params->gtype     = METIS_GTYPE_DUAL;
    params->ncommon   = 1;
    params->dbglvl    = 0;
    strcpy(params->filename, filename.c_str());
    params->outfile   = NULL;
    graph = ReadGraph(params);

}