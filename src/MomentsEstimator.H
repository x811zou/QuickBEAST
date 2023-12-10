/****************************************************************
 MomentsEstimator.H
 Copyright (C)2023 William H. Majoros (bmajoros@alumni.duke.edu).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#ifndef INCL_MomentsEstimator_H
#define INCL_MomentsEstimator_H

#include <vector>

class GridMap;
class DensityFunction;

struct EstimatedMoments {
  float mean;
  float var;
  float mode;
};

EstimatedMoments estimateMoments(const std::vector<float> &trapezoidAreas,
                                 const GridMap &gridMap,
                                 const DensityFunction &densityFunction);

#endif
