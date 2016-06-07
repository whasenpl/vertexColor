#ifndef _GRAPH_COLOR_ORDERING_H
#define _GRAPH_COLOR_ORDERING_H

#include "graphColor.h"
#include <assert.h>
#include <ext/hash_set>
#include <ext/hash_map>

struct IDStruct 
{
  unsigned int degree;
  unsigned int vertexID;
  IDStruct(unsigned int _degree, unsigned int _vertexID) : degree(_degree),vertexID(_vertexID) { }
};

inline bool operator < (const IDStruct &_a, const IDStruct &_b)
{
  if( _a.degree == _b.degree ) return _a.vertexID < _b.vertexID;
  else return _a.degree < _b.degree;
}

struct SDStruct {
  SDStruct *up;
  SDStruct *down;
  SDStruct *left;
  SDStruct *right;
  unsigned long bitColor;
  unsigned int effectiveDegree;
  unsigned int saturationDegree;
  unsigned int vertexID;
  int order;
  int color;
};

struct LFStruct {
  unsigned long bitColor;
  int color;
  int effectiveDegree;
  int saturationDegree;
  int maxVisibleColor;
  int order;
  int priority;
};

inline unsigned long packUnsignedInts(unsigned int a, unsigned int b) {
  return (((unsigned long) a) << 32L) + (unsigned long) b;
}

typedef __gnu_cxx::hash_map<unsigned int, SDStruct *> SDListType;
typedef __gnu_cxx::hash_map<unsigned long, SDStruct *> SDEDListType;

void checkAllSDVertex(SDStruct *_sdVertex, SDListType &_SDLists, SDEDListType &_SDEDLists, int _N);
void checkSDVertex(SDStruct *_sdVertex, unsigned int _vid, SDListType &_SDLists, SDEDListType &_SDEDLists);
int deleteSDVertex(SDStruct *_sdVertex, unsigned int _vid, unsigned int _maxSD, SDListType &_SDLists, SDEDListType &_SDEDLists);
int insertSDVertex(SDStruct *_sdVertex, unsigned int _vid, unsigned int _maxSD, SDListType &_SDLists, SDEDListType &_SDEDLists);

int calculateSLPriorities(sparseRowMajor<int,int>* _graph, unsigned int *_priorities);
int calculateSLLPriorities(sparseRowMajor<int,int>* _graph, unsigned int *_priorities);

void orderVertices(sparseRowMajor<int,int>* _graph, unsigned int* _orderedVertices, std::string _ordering, unsigned int _randVal = SEED_VALUE);


#endif _GRAPH_COLOR_ORDERING_H 
