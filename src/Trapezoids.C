/****************************************************************
 Trapezoids.C
 Copyright (C)2023 William H. Majoros (bmajoros@alumni.duke.edu).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#include "Trapezoids.H"
#include <iostream>

#include "DensityGrid.H"

std::vector<float> computeTrapezoidAreas(const std::vector<float> &grid) {
  std::vector<float> areas(grid.size() - 1);
  const int numTrapezoids = areas.size();

  const float dx = 1.0 / float(numTrapezoids);
  for (int i = 0; i < numTrapezoids; ++i) {
    areas[i] = dx * (grid[i] + grid[i + 1]) / 2.0;
  }

  // std::cout << "numTrapezoids=" << numTrapezoids << "\n";
  // std::cout << "dx=" << dx << "\n";
  // std::cout << "areas[0]=" << areas[0] << "\n";

  // // grid inputs 20 to a line
  // std::cout << "grid=\n";
  // for (int i = 0; i < grid.size(); ++i) {
  //   std::cout << grid[i] << " ";
  //   if (i % 20 == 19) {
  //     std::cout << "\n";
  //   }
  // }

  return areas;
}
