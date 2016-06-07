#include "graphColorOrdering.h"

int calculateSLPriorities(sparseRowMajor<int,int>* _graph, unsigned int *_priorities) 
{
  unsigned int *degrees = (unsigned int *) malloc(sizeof(unsigned int)*_graph->numRows);
  unsigned int totalVertices = _graph->numRows;
  unsigned int *queuedUpSize = (unsigned int *) malloc(sizeof(unsigned int)*totalVertices);
  unsigned int *queuedUpIndex = (unsigned int *) malloc(sizeof(unsigned int)*totalVertices);
  unsigned int **queuedUp = (unsigned int **) malloc(sizeof(unsigned int *)*totalVertices);
  unsigned int *overflow = (unsigned int *) malloc(sizeof(unsigned int)*totalVertices);
  unsigned int overflowSize = totalVertices;
  unsigned int overflowIndex = 0;
  for( unsigned int i = 0; i < totalVertices; i++ ) {
    queuedUpSize[i] = 0;
    queuedUpIndex[i] = 0;
    _priorities[i] = 0;
  }
  unsigned int maxDegree = 1;
  for( int vertexID = 0; vertexID < totalVertices; vertexID++ ) {
    degrees[vertexID] = getDegree(_graph, vertexID);
    queuedUpSize[degrees[vertexID]]++;
    if( degrees[vertexID] == 0 ) {
      _priorities[vertexID] = 1;
    } 
    else {
      maxDegree = MAX(degrees[vertexID],maxDegree);
    }
  }
  unsigned int cumSize = 0;
  for( int degree = maxDegree; degree >= 0; degree-- ) {
    cumSize += queuedUpSize[degree];
    queuedUpSize[degree] = cumSize;
    //queuedUpSize[degree] = _graph->numRows;
    queuedUp[degree] = (unsigned int *) malloc(sizeof(unsigned int)*queuedUpSize[degree]);
  }
  for( int vertexID = 0; vertexID < totalVertices; vertexID++ ) {
    queuedUp[degrees[vertexID]][queuedUpIndex[degrees[vertexID]]] = vertexID;
    queuedUpIndex[degrees[vertexID]]++;
  }
  unsigned int priority = 0;
  unsigned int remainingVertices = _graph->numRows;
  unsigned int iterations = 0;
  for( unsigned int round = 0; round <= maxDegree; round++ ) {
    bool moreVerticesToProcess = true;
    //printf("%u: %u\n", round, iterations);
    iterations = 0;
    const unsigned int MAX_ROUND_ITERATIONS = 1 << 31;
    while( moreVerticesToProcess && ( iterations < MAX_ROUND_ITERATIONS ) ) {
      priority++;
      iterations++;
      unsigned int nextRound = (iterations == MAX_ROUND_ITERATIONS) ? round+1 : round;
      moreVerticesToProcess = false;
      unsigned int *tmpOverflow = overflow;
      unsigned int tmpOverflowSize = overflowSize;
      overflow = queuedUp[round];
      overflowSize = queuedUpSize[round];
      overflowIndex = queuedUpIndex[round];
      queuedUp[round] = tmpOverflow;
      queuedUpSize[round] = tmpOverflowSize;
      queuedUpIndex[round] = 0;
      const bool PRIORITY_PREFETCH = true;
      const unsigned int PRIORITY_PREFETCH_DISTANCE = 8;
      const unsigned int PRIORITY_HALF_PREFETCH_DISTANCE = 4;
      const unsigned int PRIORITY_QUARTER_PREFETCH_DISTANCE = 2;
      if( PRIORITY_PREFETCH ) {
        for( unsigned int j = 0; j < MIN(PRIORITY_PREFETCH_DISTANCE,overflowIndex); j++ ) {
          __builtin_prefetch(&_priorities[overflow[j]],0,0);
        }
      }
      for( unsigned int j = 0; j < overflowIndex; j++ ) {
        if( PRIORITY_PREFETCH ) {
          if( j + PRIORITY_PREFETCH_DISTANCE < overflowIndex ) {
            unsigned int prefetchVertexID = overflow[j+PRIORITY_PREFETCH_DISTANCE];
            __builtin_prefetch(&_priorities[prefetchVertexID],0,0);
            __builtin_prefetch(&_graph->Starts[prefetchVertexID],0,0);
          }
          if( j + PRIORITY_HALF_PREFETCH_DISTANCE < overflowIndex ) {
            unsigned int prefetchVertexID = overflow[j+PRIORITY_HALF_PREFETCH_DISTANCE];
            if( _priorities[prefetchVertexID] == 0 ) {
              __builtin_prefetch(&_graph->ColIds[_graph->Starts[prefetchVertexID]],0,0);
            }
          }
          if( j + PRIORITY_QUARTER_PREFETCH_DISTANCE < overflowIndex ) {
            unsigned int prefetchVertexID = overflow[j+PRIORITY_QUARTER_PREFETCH_DISTANCE];
            if( _priorities[prefetchVertexID] == 0 ) {
              int *nbrs = &_graph->ColIds[_graph->Starts[prefetchVertexID]];
              for( unsigned int k = 0; k < MIN(PRIORITY_PREFETCH_DISTANCE,getDegree(_graph, prefetchVertexID)); k++ ) {
                __builtin_prefetch(&degrees[nbrs[k]],1,1);
                __builtin_prefetch(&_priorities[nbrs[k]],0,0);
              }
            }                
          }
        }
        unsigned int vertexID = overflow[j];
        if( _priorities[vertexID] == 0 ) {
          _priorities[vertexID] = priority;
          remainingVertices--;
          if( remainingVertices == 0 ) return priority;
          int *nbrs = &_graph->ColIds[_graph->Starts[vertexID]];
          unsigned int vertexDegree = getDegree(_graph, vertexID);
          for( unsigned int k = 0; k < vertexDegree; k++ ) {
            if( PRIORITY_PREFETCH ) {
              if( k + PRIORITY_PREFETCH_DISTANCE < vertexDegree ) {
                __builtin_prefetch(&degrees[nbrs[k+PRIORITY_PREFETCH_DISTANCE]],1,1);
                __builtin_prefetch(&_priorities[nbrs[k+PRIORITY_PREFETCH_DISTANCE]],0,0);
              }
            }
            if( _priorities[nbrs[k]] == 0 ) {
              degrees[nbrs[k]]--;
              int newDegree = degrees[nbrs[k]];
              if( newDegree >= nextRound ) {
                if( newDegree == nextRound ) {
                  moreVerticesToProcess = true;
                }
                if( queuedUpIndex[newDegree] >= queuedUpSize[newDegree] ) {
                  if(true) {
                    printf("blammo: %u, %u, %d,%u,%u\n", k, nbrs[k], newDegree, priority, remainingVertices);
                    printf("boom: %u, %u, %u, %u\n", queuedUpSize[newDegree],queuedUpIndex[newDegree], maxDegree, nextRound);
                  }
                  queuedUpSize[newDegree] = 2*queuedUpIndex[newDegree]+16;
                  unsigned int *tmpBuffer = NULL;
                  tmpBuffer = (unsigned int *) realloc((void*)queuedUp[newDegree],sizeof(unsigned int)*queuedUpSize[newDegree]);
                  
                  if( tmpBuffer == NULL ) {
                    printf("blammo: %u, %u, %d, %u\n", k, nbrs[k], newDegree,queuedUpSize[newDegree]);
                  }
                  queuedUp[newDegree] = tmpBuffer;
                }
                queuedUp[newDegree][queuedUpIndex[newDegree]] = nbrs[k];
                queuedUpIndex[newDegree]++;
              }
            }
          }
        }
      }
    }
  }
  for( unsigned int j = 0; j <= maxDegree; j++ ) {
    free(queuedUp[j]);
  }
  free(queuedUp);
  free(queuedUpSize);
  free(queuedUpIndex);
  free(overflow);
  return priority;
}

int calculateSLLPriorities(sparseRowMajor<int,int>* _graph, unsigned int *_priorities) 
{
  unsigned int numEdges = (unsigned int) _graph->Starts[_graph->numRows];
  unsigned int N = _graph->numRows;
  unsigned int logV = logBaseTwoFloor(_graph->numRows);
  unsigned int mask = (1 << logV) - 1;
  unsigned int shiftAmount = logV - logBaseTwoFloor(logFloor(_graph->numRows + 1)+1) - 1;
  
  unsigned int maxLogDegree = logFloor(N + 1);
  const unsigned int logJobs = logBaseTwoFloor(1 << ((logV+1)>>1));
  unsigned int verticesPerJob = _graph->numRows >> logJobs;
  int *degrees = (int *) malloc(sizeof(int)*_graph->numRows);
  unsigned int **queuedUpSize = (unsigned int **) malloc(sizeof(unsigned int *)*(1 << logJobs));
  unsigned int **queuedUpIndex = (unsigned int **) malloc(sizeof(unsigned int *)*(1 << logJobs));
  unsigned int ***queuedUp = (unsigned int ***) malloc(sizeof(unsigned int **)*(1 << logJobs));
  unsigned int **overflow = (unsigned int **) malloc(sizeof(unsigned int *)*(1 << logJobs));
  unsigned int *overflowSize = (unsigned int *) malloc(sizeof(unsigned int)*(1 << logJobs));
  unsigned int *overflowIndex = (unsigned int *) malloc(sizeof(unsigned int)*(1 << logJobs));
  for( unsigned int i = 0; i < (1 << logJobs); i++ ) {
    unsigned int totalVertices = verticesPerJob;
    if( i == ((1 << logJobs) - 1) ) totalVertices = _graph->numRows - verticesPerJob*i;  
    queuedUp[i] = (unsigned int **) malloc(sizeof(unsigned int **)*(maxLogDegree+1));
    queuedUpSize[i] = (unsigned int *) malloc(sizeof(unsigned int *)*(maxLogDegree+1));
    queuedUpIndex[i] = (unsigned int *) malloc(sizeof(unsigned int *)*(maxLogDegree+1));
    overflow[i] = (unsigned int *) malloc(sizeof(unsigned int)*totalVertices);
    overflowSize[i] = totalVertices;
    overflowIndex[i] = 0;
    for( unsigned int j = 0; j < (maxLogDegree+1); j++ ) {
      queuedUp[i][j] = (unsigned int *) malloc(sizeof(unsigned int *)*totalVertices);
      queuedUpSize[i][j] = totalVertices;
      queuedUpIndex[i][j] = 0;
    }
  }
  for( int i = 0; i < (1 << logJobs); i++ ) {
    unsigned int totalVertices = verticesPerJob;
    if( i == (1 << logJobs) - 1 ) totalVertices = _graph->numRows - verticesPerJob*i;  
    for( int j = 0; j < totalVertices; j++ ) {
      int vertexID = i*verticesPerJob + j;
      degrees[vertexID] = getDegree(_graph, vertexID);
      _priorities[vertexID] = 0;
      if( degrees[vertexID] == 0 ) _priorities[vertexID] = 1;
      else {
        int logDegree = logFloor(degrees[vertexID] + 1);
        queuedUp[i][logDegree][queuedUpIndex[i][logDegree]] = vertexID;
        queuedUpIndex[i][logDegree]++;
      }
    }
  }
  unsigned int priority = 0;
  unsigned int iterations = 0;
  for( unsigned int round = 0; round <= maxLogDegree; round++ ) {
    bool moreVerticesToProcess = true;
    //printf("%u: %u\n", round, iterations);
    iterations = 0;
    const int MAX_ROUND_ITERATIONS = 64;
    while( moreVerticesToProcess && ( iterations < MAX_ROUND_ITERATIONS ) ) {
      priority++;
      iterations++;
      unsigned int nextRound = (iterations == MAX_ROUND_ITERATIONS) ? round+1 : round;
      moreVerticesToProcess = false;
      for( int i = 0; i < (1 << logJobs); i++ ) {
        unsigned int *tmpOverflow = overflow[i];
        unsigned int tmpOverflowSize = overflowSize[i];
        overflow[i] = queuedUp[i][round];
        overflowSize[i] = queuedUpSize[i][round];
        overflowIndex[i] = queuedUpIndex[i][round];
        queuedUp[i][round] = tmpOverflow;
        queuedUpSize[i][round] = tmpOverflowSize;
        queuedUpIndex[i][round] = 0;
        const bool PRIORITY_PREFETCH = true;
        const unsigned int PRIORITY_PREFETCH_DISTANCE = 8;
        const unsigned int PRIORITY_HALF_PREFETCH_DISTANCE = 4;
        const unsigned int PRIORITY_QUARTER_PREFETCH_DISTANCE = 2;
        if( PRIORITY_PREFETCH ) {
          for( unsigned int j = 0; j < MIN(PRIORITY_PREFETCH_DISTANCE,overflowIndex[i]); j++ ) {
            __builtin_prefetch(&_priorities[overflow[i][j]],0,0);
          }
        }
        for( unsigned int j = 0; j < overflowIndex[i]; j++ ) {
          if( PRIORITY_PREFETCH ) {
            if( j + PRIORITY_PREFETCH_DISTANCE < overflowIndex[i] ) {
              unsigned int prefetchVertexID = overflow[i][j+PRIORITY_PREFETCH_DISTANCE];
              __builtin_prefetch(&_priorities[prefetchVertexID],0,0);
              __builtin_prefetch(&_graph->Starts[prefetchVertexID],0,0);
            }
            if( j + PRIORITY_HALF_PREFETCH_DISTANCE < overflowIndex[i] ) {
              unsigned int prefetchVertexID = overflow[i][j+PRIORITY_HALF_PREFETCH_DISTANCE];
              if( _priorities[prefetchVertexID] == 0 ) {
                __builtin_prefetch(&_graph->ColIds[_graph->Starts[prefetchVertexID]],0,0);
              }
            }
            if( j + PRIORITY_QUARTER_PREFETCH_DISTANCE < overflowIndex[i] ) {
              unsigned int prefetchVertexID = overflow[i][j+PRIORITY_QUARTER_PREFETCH_DISTANCE];
              if( _priorities[prefetchVertexID] == 0 ) {
                int *nbrs = &_graph->ColIds[_graph->Starts[prefetchVertexID]];
                for( unsigned int k = 0; k < MIN(PRIORITY_PREFETCH_DISTANCE,getDegree(_graph, prefetchVertexID)); k++ ) {
                  __builtin_prefetch(&degrees[nbrs[k]],1,1);
                  __builtin_prefetch(&_priorities[nbrs[k]],0,0);
                }
              }                
            }
          }
          unsigned int vertexID = overflow[i][j];
          if( _priorities[vertexID] == 0 ) {
            _priorities[vertexID] = priority;
            int *nbrs = &_graph->ColIds[_graph->Starts[vertexID]];
            unsigned int vertexDegree = getDegree(_graph, vertexID);
            for( unsigned int k = 0; k < vertexDegree; k++ ) {
              if( PRIORITY_PREFETCH ) {
                if( k + PRIORITY_PREFETCH_DISTANCE < vertexDegree ) {
                  __builtin_prefetch(&degrees[nbrs[k+PRIORITY_PREFETCH_DISTANCE]],1,1);
                  __builtin_prefetch(&_priorities[nbrs[k+PRIORITY_PREFETCH_DISTANCE]],0,0);
                }
              }
              if( _priorities[nbrs[k]] == 0 ) {
                int newDegree = __sync_sub_and_fetch(&degrees[nbrs[k]],1);
                if( (logFloor(newDegree+1) != logFloor(newDegree)) 
                    && _priorities[nbrs[k]] == 0 ) {
                  int logDegree = logFloor(newDegree);
                  if( newDegree == 0 ) logDegree = 0; 
                  if( logDegree >= nextRound ) {
                    if( logDegree == nextRound ) {
                      //degrees[nbrs[k]] = 0;
                      moreVerticesToProcess = true;
                    }
                    if( queuedUpIndex[i][logDegree] >= queuedUpSize[i][logDegree] ) {
                      //printf("resize: %u/%u/%u/%u\n", i,logDegree,queuedUpIndex[i][logDegree],queuedUpSize[i][logDegree]);
                      unsigned int *tmpBuffer = (unsigned int *) malloc(sizeof(unsigned int)*queuedUpSize[i][logDegree]*2);
                      for( unsigned int tmp_i = 0; tmp_i < queuedUpIndex[i][logDegree]; tmp_i++ ) {
                        tmpBuffer[tmp_i] = queuedUp[i][logDegree][tmp_i];
                      }
                      free(queuedUp[i][logDegree]);
                      queuedUp[i][logDegree] = tmpBuffer;
                      queuedUpSize[i][logDegree] = 2*queuedUpSize[i][logDegree];
                    }
                    queuedUp[i][logDegree][queuedUpIndex[i][logDegree]] = nbrs[k];
                    queuedUpIndex[i][logDegree]++;
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  for( unsigned int i = 0; i < (1 << logJobs); i++ ) {
    for( unsigned int j = 0; j < (maxLogDegree+1); j++ ) {
      free(queuedUp[i][j]);
    }
    free(queuedUp[i]);
    free(queuedUpSize[i]);
    free(queuedUpIndex[i]);
    free(overflow[i]);
  }
  free(queuedUp);
  free(overflowSize);
  free(overflowIndex);
  return priority;
}

int insertSDVertex(SDStruct *_sdVertex, unsigned int _vid, unsigned int _maxSD, SDListType &_SDLists, SDEDListType &_SDEDLists) {

  assert(_sdVertex[_vid].left == NULL);
  assert(_sdVertex[_vid].right == NULL);
  assert(_sdVertex[_vid].up == NULL);
  assert(_sdVertex[_vid].down == NULL);

  unsigned int SD = _sdVertex[_vid].saturationDegree;
  unsigned int ED = _sdVertex[_vid].effectiveDegree;
  unsigned long SDED = packUnsignedInts(SD, ED);  
  if (_SDEDLists.count(SDED) > 0) {

    _sdVertex[_vid].up = _SDEDLists[SDED]->up;
    _sdVertex[_vid].down = _SDEDLists[SDED];
    _SDEDLists[SDED]->up->down = &_sdVertex[_vid];
    _SDEDLists[SDED]->up = &_sdVertex[_vid];
  } 
  else {
    _SDEDLists[SDED] = &_sdVertex[_vid];
    
    assert(_SDEDLists[SDED] == &_sdVertex[_vid]);
    
    _sdVertex[_vid].up = &_sdVertex[_vid];
    _sdVertex[_vid].down = &_sdVertex[_vid];
    if (_SDLists.count(SD) > 0) {
      _sdVertex[_vid].left = _SDLists[SD]->left;
      _sdVertex[_vid].right = _SDLists[SD];
      _SDLists[SD]->left->right = &_sdVertex[_vid];
      _SDLists[SD]->left = &_sdVertex[_vid];
      if (_SDLists[SD]->effectiveDegree < ED) {
        _SDLists[SD] = &_sdVertex[_vid];
        
        assert(_SDLists[SD] == &_sdVertex[_vid]);
        
      }
    } 
    else {
      _SDLists[SD] = &_sdVertex[_vid];

      assert(_SDLists[SD] == &_sdVertex[_vid]);

      _sdVertex[_vid].right = &_sdVertex[_vid];
      _sdVertex[_vid].left = &_sdVertex[_vid];
    }
  }
  return MAX(SD,_maxSD); 
}

int deleteSDVertex(SDStruct *_sdVertex, unsigned int _vid, unsigned int _maxSD, SDListType &_SDLists, SDEDListType &_SDEDLists) {

  assert(_sdVertex[_vid].down != NULL);
  assert(_sdVertex[_vid].up != NULL);

  unsigned int SD = _sdVertex[_vid].saturationDegree;
  unsigned int ED = _sdVertex[_vid].effectiveDegree;
  unsigned long SDED = packUnsignedInts(SD, ED);
  int newMaxSD = _maxSD;
  
  assert(_SDEDLists.count(SDED) > 0);
  assert(_SDLists.count(SD) > 0);
  
  if (&_sdVertex[_vid] == _sdVertex[_vid].down) { // singleton SDED list
    assert(_SDEDLists[SDED] == &_sdVertex[_vid]);

    _SDEDLists.erase(SDED);
    if (&_sdVertex[_vid] == _sdVertex[_vid].right) { // singleton SD list
      assert(_SDLists[SD] == &_sdVertex[_vid]);

      _SDLists.erase(SD);
      if (SD == _maxSD) {
        while ((_SDLists.count(newMaxSD) == 0) && (newMaxSD > 0)) newMaxSD--; 
      }
    } else {
      _sdVertex[_vid].left->right = _sdVertex[_vid].right;
      _sdVertex[_vid].right->left = _sdVertex[_vid].left;
      if (_SDLists[SD] == &_sdVertex[_vid]) { // current maxED
        _SDLists[SD] = _SDLists[SD]->right;
      }
    }
    if ((SD == _maxSD) && (_SDLists.count(newMaxSD) > 0)) { // this is _maxSD, so we need to find the next maxED
      SDStruct *sdPtr = _SDLists[newMaxSD];
      SDStruct *sdStart = _SDLists[newMaxSD];
      while (sdPtr->right != sdStart) {
        sdPtr = sdPtr->right;
        if (_SDLists[newMaxSD]->effectiveDegree < sdPtr->effectiveDegree) {
          _SDLists[newMaxSD] = sdPtr;
        }
      }
    }
  } else {
    if (_SDEDLists[SDED] == &_sdVertex[_vid]) { 
      // head of SDED list
      _SDEDLists[SDED] = _sdVertex[_vid].down;
      if (_SDLists[SD] == &_sdVertex[_vid]) {
        _SDLists[SD] = _SDEDLists[SDED];
      }
      _sdVertex[_vid].right->left = _SDEDLists[SDED];
      _SDEDLists[SDED]->right = _sdVertex[_vid].right;
      _sdVertex[_vid].left->right = _SDEDLists[SDED];
      _SDEDLists[SDED]->left = _sdVertex[_vid].left;
    } 
    _sdVertex[_vid].up->down = _sdVertex[_vid].down;
    _sdVertex[_vid].down->up = _sdVertex[_vid].up;    
  }
  _sdVertex[_vid].left = NULL;
  _sdVertex[_vid].right = NULL;
  _sdVertex[_vid].up = NULL;
  _sdVertex[_vid].down = NULL;
  
  return newMaxSD;
}

void checkSDVertex(SDStruct *_sdVertex, unsigned int _vid, SDListType &_SDLists, SDEDListType &_SDEDLists) {
  unsigned int SD = _sdVertex[_vid].saturationDegree;
  unsigned int ED = _sdVertex[_vid].effectiveDegree;
  unsigned long SDED = packUnsignedInts(SD, ED);
  if (_sdVertex[_vid].order == -1) {
    assert(_SDLists.count(SD) > 0);
    if( _SDLists[SD] == &_sdVertex[_vid] ) {
      assert(_sdVertex[_vid].right != NULL);
      assert(_sdVertex[_vid].left != NULL);
    }
    assert(_SDEDLists.count(SDED) > 0);
    if( _SDEDLists[SDED] == &_sdVertex[_vid] ) {
      assert(_sdVertex[_vid].right != NULL);
      assert(_sdVertex[_vid].left != NULL);
    }
    assert(_sdVertex[_vid].up->down == &_sdVertex[_vid]);
    assert(_sdVertex[_vid].down->up == &_sdVertex[_vid]);
    if ((_sdVertex[_vid].down == &_sdVertex[_vid])
        || (_sdVertex[_vid].right != NULL)
        || (_sdVertex[_vid].left != NULL)){ // singleton SDED
      assert(_SDEDLists[SDED] == &_sdVertex[_vid]);
      assert(_sdVertex[_vid].right->left == &_sdVertex[_vid]);
      assert(_sdVertex[_vid].left->right == &_sdVertex[_vid]);
    }
  } else {
    assert(_sdVertex[_vid].up == NULL);
    assert(_sdVertex[_vid].down == NULL);
    assert(_sdVertex[_vid].right == NULL);
    assert(_sdVertex[_vid].left == NULL);
  }
}

void checkAllSDVertex(SDStruct *_sdVertex, SDListType &_SDLists, SDEDListType &_SDEDLists, int _N) {
  for (int i = 0; i < _N; i++) {
    checkSDVertex(_sdVertex, i, _SDLists, _SDEDLists);
  }
}

void orderVertices(sparseRowMajor<int,int>* _graph, unsigned int* _orderedVertices, std::string _ordering, unsigned int _randVal)
{
  unsigned int* counts = (unsigned int*) calloc(sizeof(unsigned int), _graph->numRows);
  unsigned int logV = logBaseTwoFloor(_graph->numRows);
  if (_ordering.compare("R") == 0)
  {
    unsigned int mask = (1 << logV) - 1;
    for (int i = 0; i < _graph->numRows; i++) {
      counts[hashR(i,mask,_randVal)]++;
    }
    for( int i = 1; i < _graph->numRows; i++)
      counts[i] += counts[i-1];
    for (int i = 0; i < _graph->numRows; i++) {
      unsigned int hash = hashR(i,mask,_randVal);
      _orderedVertices[_graph->numRows - counts[hash]] = i;
      counts[hash]--;
    }    
  }
  else if (_ordering.compare("LF") == 0)
  {
    unsigned int* orderedVerticesTmp = (unsigned int*) calloc(sizeof(unsigned int), _graph->numRows);
    unsigned int mask = (1 << logV) - 1;
    for (int i = 0; i < _graph->numRows; i++) {
      counts[hashR(i,mask,_randVal)]++;
    }
    for( int i = 1; i < _graph->numRows; i++)
      counts[i] += counts[i-1];
    for (int i = 0; i < _graph->numRows; i++) {
      unsigned int hash = hashR(i,mask,_randVal);
      orderedVerticesTmp[_graph->numRows - counts[hash]] = i;
      counts[hash]--;
    }    
    for (int i = 0; i < _graph->numRows; i++) {
      counts[i] = 0;
    }
    for (int i = 0; i < _graph->numRows; i++) {
      counts[hashLF(_graph, i)]++;
    }
    for( int i = 1; i < _graph->numRows; i++) {
      counts[i] += counts[i-1];
    }
    for (int i = 0; i < _graph->numRows; i++) {
      unsigned int hash = hashLF(_graph, orderedVerticesTmp[i]);
      _orderedVertices[_graph->numRows - counts[hash]] = orderedVerticesTmp[i];
      counts[hash]--;
    }    
  }
  else if (_ordering.compare("LLF") == 0)
  {
    unsigned int shiftAmount = logV - logBaseTwoFloor(logFloor(_graph->numRows + 1)+1) - 1;
    
    for (int i = 0; i < _graph->numRows; i++) {
      counts[hashLLF(_graph, i,shiftAmount,_randVal)]++;
    }
    for( int i = 1; i < _graph->numRows; i++)
      counts[i] += counts[i-1];
    for (int i = 0; i < _graph->numRows; i++) {
      unsigned int hash = hashLLF(_graph, i,shiftAmount,_randVal);
      _orderedVertices[_graph->numRows - counts[hash]] = i;
      counts[hash]--;
    }    
  }
  else if (_ordering.compare("FF") == 0)
  {
    for(int i = 0; i < _graph->numRows; i++)
      _orderedVertices[i] = i;
  }
  else if (_ordering.compare("ID") == 0) 
  {
    int N = _graph->numRows;
    orderVertices(_graph, _orderedVertices, std::string("LF"));
    SDStruct *sdVertex = (SDStruct *) malloc(sizeof(SDStruct)*N);

    int numColored = 0;
    for (int i = 0; i < N; i++) {
      sdVertex[i].effectiveDegree = getDegree(_graph, i);
      sdVertex[i].saturationDegree = 0;
      sdVertex[i].vertexID = i;
      sdVertex[i].order = -1;
      sdVertex[i].up = NULL;
      sdVertex[i].down = NULL;
      sdVertex[i].left = NULL;
      sdVertex[i].right = NULL;
    }

    unsigned int maxColor = 0;
    unsigned int maxSD = 0;

    SDListType SDLists;
    SDEDListType SDEDLists;

    for (unsigned int i = 0; i < N; i++) {
      unsigned int vid = _orderedVertices[i];
      if (sdVertex[vid].order == -1) {
        maxSD = insertSDVertex(sdVertex, vid, maxSD, SDLists, SDEDLists);
      }
    }
        
    while (numColored < N) {
      int vid = SDLists[maxSD]->vertexID;
      
      sdVertex[vid].order = numColored;
      numColored++;
      maxSD = deleteSDVertex(sdVertex, vid, maxSD, SDLists, SDEDLists);
      
      for (int i = 0; i < getDegree(_graph, vid); i++) {
        unsigned int nbrID = getNeighbor(_graph, vid,i);

        if (sdVertex[nbrID].order == -1) {
          maxSD = deleteSDVertex(sdVertex, nbrID, maxSD, SDLists, SDEDLists);

          sdVertex[nbrID].saturationDegree++;
          // reinsert nbr w/ new saturation and effective degrees
          maxSD = insertSDVertex(sdVertex, nbrID, maxSD, SDLists, SDEDLists);
        }
      }
            
    }
    for (int i = 0; i < N; i++) {
      _orderedVertices[sdVertex[i].order] = i;
    }
  }
  else if (_ordering.compare("SD") == 0)
  {    
    int N = _graph->numRows;
    orderVertices(_graph, _orderedVertices, std::string("LF"));
    SDStruct *sdVertex = (SDStruct *) malloc(sizeof(SDStruct)*N);

    int numColored = 0;
    for (int i = 0; i < N; i++) {
      sdVertex[i].effectiveDegree = getDegree(_graph, i);
      sdVertex[i].saturationDegree = 0;
      sdVertex[i].bitColor = 0L;
      sdVertex[i].vertexID = i;
      sdVertex[i].color = -1;
      sdVertex[i].order = -1;
      sdVertex[i].up = NULL;
      sdVertex[i].down = NULL;
      sdVertex[i].left = NULL;
      sdVertex[i].right = NULL;
    }

    unsigned int maxColor = 0;
    unsigned int maxSD = 0;

    SDListType SDLists;
    SDEDListType SDEDLists;
    __gnu_cxx::hash_set<unsigned long> largeColors;

    for (unsigned int i = 0; i < N; i++) {
      unsigned int vid = _orderedVertices[i];
      if (sdVertex[vid].order == -1) {
        maxSD = insertSDVertex(sdVertex, vid, maxSD, SDLists, SDEDLists);
      }
    }
    
    while (numColored < N) {
      while (SDLists[maxSD]->order != -1) {
        maxSD = deleteSDVertex(sdVertex, SDLists[maxSD]->vertexID, maxSD, SDLists, SDEDLists);
      }
      int vid = SDLists[maxSD]->vertexID;
      
      sdVertex[vid].order = numColored;
      numColored++;
      unsigned int color = 0;
      unsigned long thisBitColor = 0L;
      const unsigned long ALL_ONES = 0xFFFFFFFFFFFFFFFF;
      if (sdVertex[vid].bitColor != ALL_ONES) {
        color = findFirstBitSet(sdVertex[vid].bitColor ^ (sdVertex[vid].bitColor + 1));
        thisBitColor = 1L << color;
      } else {
        unsigned char *colorArray = (unsigned char *) calloc(getDegree(_graph, vid) + 1, sizeof(unsigned char));
        for (int i = 0; i < getDegree(_graph, vid); i++) {
          int nbrID = getNeighbor(_graph, vid, i);
          if ((sdVertex[nbrID].color <= getDegree(_graph, vid)) 
              && (sdVertex[nbrID].color >= 0))
            colorArray[sdVertex[nbrID].color] = 1;
        }
        while (colorArray[color] != 0) color++;
      }
      sdVertex[vid].color = color;
      maxSD = deleteSDVertex(sdVertex, vid, maxSD, SDLists, SDEDLists);
      if (color > maxColor) {
        maxColor = color;
      }
      for (int i = 0; i < getDegree(_graph, vid); i++) {
        unsigned int nbrID = getNeighbor(_graph, vid,i);
        if (sdVertex[nbrID].order == -1) {
          maxSD = deleteSDVertex(sdVertex, nbrID, maxSD, SDLists, SDEDLists);

          sdVertex[nbrID].effectiveDegree--;
          // if nbrID wasn't deferred, update saturationDegree and re-insert

          if (color < 64) {
            if ((sdVertex[nbrID].bitColor & thisBitColor) == 0L) {
              sdVertex[nbrID].bitColor |= thisBitColor;
              sdVertex[nbrID].saturationDegree++;
            }
          } else {
            unsigned long vertexIDAndColor = packUnsignedInts(nbrID, color);
            if (largeColors.count(vertexIDAndColor) == 0) {
              largeColors.insert(vertexIDAndColor);
              sdVertex[nbrID].saturationDegree++;
            }
          }     
          // reinsert nbr w/ new saturation and effective degrees
          maxSD = insertSDVertex(sdVertex, nbrID, maxSD, SDLists, SDEDLists);
        }          
      }
    }
    for (int i = 0; i < N; i++) {
      _orderedVertices[sdVertex[i].order] = i;
    }
  }
  else if (_ordering.compare("SL") == 0) {
    const bool random_tie_breaking = true;
    unsigned int *priorities = (unsigned int *) malloc(sizeof(unsigned int)*_graph->numRows);
    int maxPriority = calculateSLPriorities(_graph, priorities);
    unsigned int* orderedVerticesTmp;
    if( random_tie_breaking ) {
      orderedVerticesTmp = (unsigned int*) calloc(sizeof(unsigned int), _graph->numRows);
      unsigned int mask = (1 << logV) - 1;
      for (int i = 0; i < _graph->numRows; i++) {
        counts[hashR(i,mask,_randVal)]++;
      }
      for( int i = 1; i < _graph->numRows; i++)
        counts[i] += counts[i-1];
      for (int i = 0; i < _graph->numRows; i++) {
        unsigned int hash = hashR(i,mask,_randVal);
        orderedVerticesTmp[_graph->numRows - counts[hash]] = i;
        counts[hash]--;
      }    
      for (int i = 0; i < _graph->numRows; i++) {
        counts[i] = 0;
      }
    }
    for (int i = 0; i < _graph->numRows; i++) {
      counts[priorities[i]]++;
    }
    for( int i = 1; i < _graph->numRows; i++) {
      counts[i] += counts[i-1];
    }
    for (int i = 0; i < _graph->numRows; i++) {
      unsigned int hash;
      if( random_tie_breaking ) {
        hash = priorities[orderedVerticesTmp[i]];
        _orderedVertices[_graph->numRows - counts[hash]] = orderedVerticesTmp[i];
      } else {
        hash = priorities[i];
        _orderedVertices[_graph->numRows - counts[hash]] = i;
      }
      counts[hash]--;
    }          
  }
  else if (_ordering.compare("SLL") == 0 ) {
    unsigned int *priorities = (unsigned int *) malloc(sizeof(unsigned int)*_graph->numRows);
    int maxPriority = calculateSLLPriorities(_graph, priorities);
    unsigned int* orderedVerticesTmp = (unsigned int*) calloc(sizeof(unsigned int), _graph->numRows);
    unsigned int mask = (1 << logV) - 1;
    for (int i = 0; i < _graph->numRows; i++) {
      counts[hashR(i,mask,_randVal)]++;
    }
    for( int i = 1; i < _graph->numRows; i++)
      counts[i] += counts[i-1];
    for (int i = 0; i < _graph->numRows; i++) {
      unsigned int hash = hashR(i,mask,_randVal);
      orderedVerticesTmp[_graph->numRows - counts[hash]] = i;
      counts[hash]--;
    }    
    for (int i = 0; i < _graph->numRows; i++) {
      counts[i] = 0;
    }
    for (int i = 0; i < _graph->numRows; i++) {
      counts[priorities[i]]++;
    }
    for( int i = 1; i < _graph->numRows; i++) {
      counts[i] += counts[i-1];
    }
    for (int i = 0; i < _graph->numRows; i++) {
      unsigned int hash = priorities[orderedVerticesTmp[i]];
      _orderedVertices[_graph->numRows - counts[hash]] = orderedVerticesTmp[i];
      counts[hash]--;
    }          
  }
  else
    printf("not a valid ordering\n");

  free(counts);
}

 
