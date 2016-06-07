#ifndef _GRAPH_COLOR_GREEDY_H
#define _GRAPH_COLOR_GREEDY_H

#include "graphColor.h"
#include "graphColorOrdering.h"

int greedy(sparseRowMajor<int,int>* _graph, unsigned int* _orderedVertices, int* _vertexColors);
int colorGraphGreedy(sparseRowMajor<int,int>* _graph, int* _vertexColors, std::string _ordering, unsigned int _seed = SEED_VALUE);

#endif _GRAPH_COLOR_GREEDY_H