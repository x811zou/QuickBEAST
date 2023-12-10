/****************************************************************
 DensityGrid.H
 Copyright (C)2023 William H. Majoros (bmajoros@alumni.duke.edu).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#ifndef INCL_DensityGrid_H
#define INCL_DensityGrid_H

#include <functional>
#include <vector>

class GridMap;

std::vector<float> fillGrid(const std::function<float(float)> &f,
                            const GridMap &map);

#endif
