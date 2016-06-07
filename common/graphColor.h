#ifndef GRAPH_COLOR_H_
#define GRAPH_COLOR_H_

#include <cilk/cilk.h>
#include "graph.h"

#define OVERFLOW_BIT_COLOR (1L << 63L)
#define MAX_BIT_COLOR 62
#define NO_BIT_COLOR_AVAILABLE 0x7FFFFFFFFFFFFFFF
#define SOFTWARE_PREFETCHING_FLAG1 true
#define PREFETCH_DISTANCE 6
#define HALF_PREFETCH_DISTANCE 2
#define MAX(a, b) ((a < b) ? b : a)
#define MIN(X, Y) ((X) < (Y) ? (X) : (Y))
#define SEED_VALUE 0xF1807D63

inline unsigned int hashVertexID(unsigned int _id, unsigned int _randVal = 0xF1807D63) {
  return _mm_crc32_u32(_randVal, _id);
}

inline bool isPowerOfTwo(int _num) {
  if (_num > 0) {
    return ( ( _num & (_num - 1)) == 0 );
  } else {
    return true;
  }
}

inline int logBaseTwoFloor(unsigned int _num) {
  return 31 - __builtin_clz(_num);
}

inline int logFloor(unsigned int _num) {
  return 31 - __builtin_clz(_num);
}

inline int getNeighbor(sparseRowMajor<int, int>* _graph, int vid, int i) {
  return _graph->ColIds[_graph->Starts[vid] + i];
}

inline int getDegree(sparseRowMajor<int, int>* _graph, int v) {
  return _graph->Starts[v+1] - _graph->Starts[v];
}

inline unsigned int getDegree(sparseRowMajor<int, int>* _graph, unsigned int _vertex_ID) {
  return (unsigned int) (_graph->Starts[_vertex_ID+1] - _graph->Starts[_vertex_ID]);
}

inline unsigned int hashR(unsigned int _vertexID, unsigned int _mask, unsigned int _randVal = SEED_VALUE) {
  return _mask & hashVertexID(_vertexID, _randVal);
}

inline unsigned int hashLF(sparseRowMajor<int, int>* _graph, unsigned int _vertexID) {
  return _graph->Starts[_vertexID+1] - _graph->Starts[_vertexID];
}

inline unsigned int hashLLF(sparseRowMajor<int, int>* _graph, unsigned int _vertexID, int _shiftAmount, unsigned int _randVal = SEED_VALUE) {
  unsigned int mask = (1 << _shiftAmount) - 1;
  return (logFloor(hashLF(_graph, _vertexID)+1) << _shiftAmount) + hashR(_vertexID, mask, _randVal);
}

inline unsigned int hashSLL(unsigned int _vertexID, unsigned int *_priority) {
  return _priority[_vertexID];
}

inline unsigned int hashSL(unsigned int _vertexID, unsigned int *_priority) {
  return _priority[_vertexID];
}

inline void setBit(unsigned long* _bitVector, unsigned int _index) {
  _bitVector[_index >> 6] |= (1L << (_index & 63L));
}

inline int getLowestUnsetBit(unsigned long* _bitVector) {
  const unsigned long ALL_ONES = 0xFFFFFFFFFFFFFFFF;
  int i = 0;
  while (_bitVector[i] == ALL_ONES) i++;
  return (i*64) + 63 - __builtin_clzl(_bitVector[i] ^ (_bitVector[i]+1));
}

inline bool isSuccessor(unsigned int _myHash, unsigned int _hisHash, unsigned int _myID, unsigned int _hisID) {
  if (_myHash != _hisHash) {
    return _myHash > _hisHash;
  } else {
    return _myID < _hisID;
  }
}

inline int findFirstBitSet(unsigned long _color) {
  return 63 - __builtin_clzl(_color);
}

#endif  // GRAPH_COLOR_H_
