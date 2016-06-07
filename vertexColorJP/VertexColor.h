#ifndef _GRAPH_COLOR_JP_H
#define _GRAPH_COLOR_JP_H

#include "graphColor.h"
#include "graphColorOrdering.h"

#define CILK_FOR_THRESHOLD 256
#define COUNTER_THRESHOLD 512
#define NUM_LEAF_MEMBERS 64

struct LeafClass
{
  volatile unsigned long bitColor;
  
  inline unsigned long updateLeafCounter(unsigned long _bitColor);
  volatile unsigned short mutex;
  volatile unsigned short counter;
};

class Vertex
{
 public:
  unsigned long finalColor;
  unsigned int numPredecessors;
  unsigned int numSuccessors; 
  int firstSuccessor;
  int secondSuccessor;
  unsigned int edgeIndex;
  unsigned short counter;
  unsigned char logSize;
  volatile unsigned char mutex;
  
  void vertexInitR(unsigned char *_tournamentArray, unsigned int _vertexID, unsigned int *_neighbors, unsigned int _neighborIndex, unsigned int _numNeighbors, unsigned int _mask, unsigned int _randVal);
  void vertexInitFF(unsigned char *_tournamentArray, unsigned int _vertexID, unsigned int *_neighbors, unsigned int _neighborIndex, unsigned int _numNeighbors);
  void vertexInitLF(sparseRowMajor<int,int>* _graph, unsigned char *_tournamentArray, unsigned int _vertexID, unsigned int *_neighbors, unsigned int _neighborIndex, unsigned int _numNeighbors, unsigned int _mask, unsigned int _randVal);
  void vertexInitLLF(sparseRowMajor<int,int>* _graph, unsigned char *_tournamentArray, unsigned int _vertexID, unsigned int *_neighbors, unsigned int _neighborIndex, unsigned int _numNeighbors, unsigned int _shiftAmount, unsigned int _randVal);
  void vertexInitSLL(unsigned char *_tournamentArray, unsigned int _vertexID, unsigned int *_neighbors, unsigned int _neighborIndex, unsigned int _numNeighbors, unsigned int *_priority, unsigned int _mask, unsigned int _randVal);
  void vertexInitSL(unsigned char *_tournamentArray, unsigned int _vertexID, unsigned int *_neighbors, unsigned int _neighborIndex, unsigned int _numNeighbors, unsigned int *_priority, unsigned int _mask, unsigned int _randVal);
  inline void scalarInit(unsigned int _neighborIndex);
  inline void tournamentInit(unsigned char *_tournamentArray, unsigned int _neighborIndex, unsigned int *_predecessor);
  inline unsigned long competeInTournament(unsigned int _hash, unsigned int _color, unsigned char *_tournamentArray);
  inline void colorVertex(Vertex *_vertices, unsigned int _vertexID, unsigned char *_tournamentArray, unsigned int *_neighbors, unsigned long _bitColor);
  inline unsigned long updateLocalCounter(unsigned long _bitColor);
};

inline void prefetchVertexInit(sparseRowMajor<int,int>* _graph, Vertex *_vertices, unsigned char *_tournamentArray, unsigned int _index, bool _fetchDegree);
void prefetchSuccessors(unsigned int *_successors,unsigned int _numSuccessors,Vertex *_vertices,int firstSuccessor,int secondSuccessor);
void processSuccessors(unsigned int *_successors,unsigned int _numSuccessors,Vertex *_vertices,unsigned int _vertexID,unsigned long _finalColor,unsigned char *_tournamentArray, unsigned int *_neighbors);

void colorGraphJPInit(sparseRowMajor<int,int>* _graph, Vertex *_vertices, unsigned char *_tournamentArray, std::string _ordering, unsigned int seed = SEED_VALUE);
void JP(sparseRowMajor<int,int>* _graph, Vertex *_vertices, unsigned char *_tournamentArray);
void colorGraphJP(sparseRowMajor<int,int>* _graph, Vertex *_vertices, unsigned char *_tournamentArray, std::string _ordering, unsigned int _seed = SEED_VALUE);
int colorGraphJP(sparseRowMajor<int,int>* _graph, int* _vertexColors, std::string _ordering, unsigned int _seed = SEED_VALUE);
 
#endif _GRAPH_COLOR_JP_H
