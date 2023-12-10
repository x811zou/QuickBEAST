/****************************************************************
 GridMap.C
 Copyright (C)2023 William H. Majoros (bmajoros@alumni.duke.edu).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#include "GridMap.H"
#include <iostream>

GridMap::GridMap(int gridSize, float min, float max)
    : gridSize(gridSize), min(min), max(max) {
  // ctor
}

int GridMap::pToIndex(float p) const {
  p = (p - min) / (max - min);
  return int(p * gridSize);
}

float GridMap::indexToP(int index) const {
  const float u = float(index) / float(gridSize - 1);
  return u * max + (1 - u) * min;
}
