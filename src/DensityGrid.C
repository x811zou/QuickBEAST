/****************************************************************
 DensityGrid.C
 Copyright (C)2023 William H. Majoros (bmajoros@alumni.duke.edu).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/

#include "DensityGrid.H"

#include "GridMap.H"

std::vector<float> fillGrid(const std::function<float(float)> &f,
                            const GridMap &gridMap) {
  std::vector<float> grid(gridMap.getGridSize());
  for (int i = 0; i < gridMap.getGridSize(); ++i) {
    const float p = gridMap.indexToP(i);
    grid[i] = f(p);
  }
  return grid;
}