/****************************************************************
 GridMap.H
 Copyright (C)2023 William H. Majoros (bmajoros@alumni.duke.edu).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#ifndef INCL_GridMap_H
#define INCL_GridMap_H
#include <iostream>
using namespace std;

class GridMap {
public:
  GridMap(int gridSize, float min, float max);

  int pToIndex(float p) const;
  float indexToP(int index) const;

  int getGridSize() const { return gridSize; }
  float getMin() const { return min; }
  float getMax() const { return max; }
  float getRange() const { return max - min; }
  float getStep() const { return getRange() / (gridSize - 1); }

private:
  const int gridSize;
  const float min;
  const float max;
};

#endif
