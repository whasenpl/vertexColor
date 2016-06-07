// Copyright 2014 Will Hasenplaugh

#include <ctime>
#include <string>

#include "/afs/csail/proj/courses/6.172/cilkutil/include/cilktools/cilkview.h"

#include "graph.h"
#include "graphIO.h"
#include "graphUtils.h"
#include "graphColor.h"
#include "graphColorGreedy.h"
#include "graphColorJP.h"

double getTime() {
    struct timeval t;
    struct timezone tzp;
    gettimeofday(&t, &tzp);
    return t.tv_sec + t.tv_usec*1e-6;
}

int cilk_main(
  std::string _inputname,
  std::string _coloringAlgorithm,
  std::string _numWorkers,
  std::string _ordering) {
  __cilkrts_set_param("nworkers", _numWorkers.c_str());

  sparseRowMajor<int, int> inputGraph;
  // Parse the graph in .mtx or adjacency list representation.
  if (_inputname.find(".mtx") != string::npos) {
    graph<int> inter = graphFromEdges<int>(edgesFromMtxFile<int>(_inputname.c_str()), true);
    inputGraph = sparseFromGraph<int>(inter);
  } else if (_inputname.find(".edges") != string::npos) {
    graph<int> inter = graphFromEdges<int>(readEdgeArrayFromFile<int>((char*)_inputname.c_str()), true);
    inputGraph = sparseFromGraph<int>(inter);
  } else {
    graph<int> inter = graphFromEdges<int>(edgesFromGraph<int>(readGraphFromFile<int>((char*)_inputname.c_str())), true);
    inputGraph = sparseFromGraph<int>(inter);
  }


  // Initialize the vertex color array, and the permutation.
  // int* vertexColors = (int*) malloc(sizeof(int)*inputGraph.numRows);
  if (inputGraph.Values == NULL) {
    inputGraph.Values = (int*) malloc(sizeof(int)*inputGraph.numRows);
  }
  int* vertexColors = inputGraph.Values;

  unsigned int seed = 1234134;
  double beginTime = getTime();
  int numColors;
  if (_coloringAlgorithm.compare("Greedy") == 0) {
    numColors = colorGraphGreedy(&inputGraph, vertexColors, _ordering, seed);
  } else if (_coloringAlgorithm.compare("JP") == 0) {
    numColors = colorGraphJP(&inputGraph, vertexColors, _ordering, seed);
  } else {
    printf("Invalid coloring algorithm\n");
  }
  double endTime = getTime();

  for (int i = 0; i < inputGraph.numRows; i++) {
    for (int j = 0; j < getDegree(&inputGraph, i); j++) {
      int nbr = getNeighbor(&inputGraph, i, j);
      if (vertexColors[i] == vertexColors[nbr]) {
        printf("Error: vertex %d and vertex %d are neighbors and share color %d\n", i, nbr, vertexColors[i]);
        return -1;
      }
    }
  }
  printf("Total colors: %d\n", numColors);
  printf("Total time: %f\n", endTime - beginTime);

  return 0;
}

int main(int argc, char **argv) {
  if (argc < 5) {
    printf("ERROR: not enough arguments.\n");
    printf("Usage:\n>>>%s <file> <coloring algorithm {Greedy, JP}> ");
    printf("<numWorkers> <ordering {FF, R, LF, LLF, SL, SLL, ID, SD}>\n", argv[0]);
    printf("FF:\tFirst Fit / Input order\n");
    printf("R:\tRandom\n");
    printf("LF:\tLargest degree first\n");
    printf("LLF:\tLargest log-degree\n");
    printf("SL:\tAllwright's parallel Smallest last\n");
    printf("SLL:\tSmallest log-degree last\n");
    printf("ID:\tIncidence degree\n");
    printf("SD:\tSaturation degree\n");
  } else {
    cilk_main(argv[1], argv[2], argv[3], argv[4]);
    return 0;
  }

  return 0;
}

