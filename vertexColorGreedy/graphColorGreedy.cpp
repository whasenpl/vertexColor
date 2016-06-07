#include "graphColorGreedy.h"

int greedy(sparseRowMajor<int,int>* _graph, unsigned int* _orderedVertices, int* _vertexColors)
{
  for (int i = 0; i < _graph->numRows; i++) {
    _vertexColors[i] = -1;
  }
  int maxColor = 0;
  for(int k = 0; k < _graph->numRows; k++) {
    unsigned int vid = _orderedVertices[k];
    int* neighbors = &_graph->ColIds[_graph->Starts[vid]];
    unsigned int degree = (_graph->Starts[vid+1] - _graph->Starts[vid]);
    char* neighborColors;
    unsigned long bitColor = OVERFLOW_BIT_COLOR;
    if( degree > MAX_BIT_COLOR ) {
      neighborColors = (char*) calloc(degree + 1, sizeof(char));
    }
    int effectiveDegree = 0;
    for (int i = 0; i < degree; i++) {
      int neighbor = neighbors[i];
      int color = _vertexColors[neighbor];
      if (color >= 0 && color <= degree) {
        if( color <= MAX_BIT_COLOR ) {
          bitColor |= (1L << color);
        }
        else {
          neighborColors[color] = 1;
        }
      }
    }
    int chosenColor = -1;
    if( (bitColor & NO_BIT_COLOR_AVAILABLE) == NO_BIT_COLOR_AVAILABLE ) {
      for (int j = MAX_BIT_COLOR + 1; j < degree + 1; j++ ) {
        if (neighborColors[j] == 0){
          chosenColor = j;
          break;
        }
      }
    } else {
      chosenColor = findFirstBitSet(bitColor ^ (bitColor+1));
    }
    if( degree > MAX_BIT_COLOR ) {
      free(neighborColors);
    }
    _vertexColors[vid] = chosenColor;
    
    maxColor = MAX(maxColor, chosenColor);
  }
  return maxColor + 1;
}

int colorGraphGreedy(
  sparseRowMajor<int,int>* _graph, 
  int* _vertexColors, 
  std::string _ordering, 
  unsigned int _seed)
{
  unsigned int* orderedVertices = (unsigned int*) malloc(sizeof(unsigned int) * _graph->numRows);

  orderVertices(_graph, orderedVertices, _ordering, _seed);
  int returnVal = greedy(_graph, orderedVertices, _vertexColors);

  free(orderedVertices);
  return returnVal;
}


