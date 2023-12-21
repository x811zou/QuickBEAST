/****************************************************************
 MomentsEstimator.C
 Copyright (C)2023 William H. Majoros (bmajoros@alumni.duke.edu).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#include "MomentsEstimator.H"

#include "DensityFunction.H"
#include "DensityGrid.H"
#include "GridMap.H"
#include "Trapezoids.H"
#include <algorithm>
#include <iomanip>
#include <iostream>
#include <math.h>
#include <numeric>

float estimateMean(const MeanOptions &options,
                   const DensityFunction &densityFunction) {
  GridMap gridMap(options.numTrapezoids + 1, 0, 1);
  const auto &&grid = fillGrid(densityFunction, gridMap);
  const auto &&trapezoidAreas = computeTrapezoidAreas(grid);

  const int numTrap = trapezoidAreas.size();
  const float totalArea =
      std::accumulate(trapezoidAreas.begin(), trapezoidAreas.end(), 0.0f);
  const float partition = totalArea;

  if (totalArea == 0) {
    return 0;
  }

  float sumX = 0;
  for (int i = 0; i < numTrap; ++i) {
    const float PofX = trapezoidAreas[i] / partition;
    const float begin = gridMap.indexToP(i), end = gridMap.indexToP(i + 1);
    const float p = (begin + end) / 2;

    const float x = p;
    sumX += x * PofX; // expectation: sum of X*P(X)
  }

  return sumX;
}

float estimateMode(const ModeOptions &options,
                   const DensityFunction &densityFunction) {
  float startP = 0, endP = 1;
  float mode = 0;
  float maxDensity = 0;
  while (endP - startP > options.subgridThreshold) {
    GridMap subgrid(options.subgridSize, startP, endP);
    mode = subgrid.indexToP(options.subgridSize / 2);
    maxDensity = densityFunction.logLikelihood(mode);

    for (int i = 0; i < subgrid.getGridSize(); ++i) {
      const float p = subgrid.indexToP(i);
      const float density = densityFunction.logLikelihood(p);
      if (density > maxDensity) {
        maxDensity = density;
        mode = p;
      }
    }

    const float nextHalfWidth = subgrid.getStep() * options.subgridStepSlop;
    startP = mode - nextHalfWidth;
    endP = mode + nextHalfWidth;
  }

  return mode;
}

std::vector<std::pair<float, float>>
computeDistribution(const DistributionOptions &options,
                    const DensityFunction &densityFunction) {
  std::vector<std::pair<float, float>> distribution;
  distribution.reserve(options.numSamples);

  GridMap gridMap(options.numSamples, 0, 1);
  for (int i = 0; i < options.numSamples; ++i) {
    const float p = gridMap.indexToP(i);
    if (options.useLogLikelihood) {
      distribution.emplace_back(p, densityFunction.logLikelihood(p));
    } else {
      distribution.emplace_back(p, densityFunction(p));
    }
  }
  return distribution;
}