/****************************************************************
 MomentsEstimator.C
 Copyright (C)2023 William H. Majoros (bmajoros@alumni.duke.edu).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#include "MomentsEstimator.H"

#include "DensityFunction.H"
#include "GridMap.H"
#include <algorithm>
#include <iostream>
#include <math.h>
#include <numeric>

EstimatedMoments estimateMoments(const std::vector<float> &trapezoidAreas,
                                 const GridMap &gridMap,
                                 const DensityFunction &densityFunction) {
  const int numTrap = trapezoidAreas.size();
  const float totalArea =
      std::accumulate(trapezoidAreas.begin(), trapezoidAreas.end(), 0.0f);
  const float partition = totalArea;

  if (totalArea == 0) {
    return {
        .mean = 0,
        .var = 0,
        .mode = 0,
    };
  }

  // print trapezoids areas 20 to line
  // std::cout << "trapezoidAreas=\n";
  // for (int i = 0; i < numTrap; ++i) {
  //   std::cout << trapezoidAreas[i] << " ";
  //   if (i % 20 == 19) {
  //     std::cout << "\n";
  //   }
  // }

  // std::cout << "totalArea=" << totalArea << "\n";

  float sumX = 0, sumXX = 0, maxTrapI = 0, maxArea = -HUGE_VAL;
  for (int i = 0; i < numTrap; ++i) {
    const float PofX = trapezoidAreas[i] / partition;
    const float begin = gridMap.indexToP(i), end = gridMap.indexToP(i + 1);
    const float p = (begin + end) / 2;

    const float x = log(p / (1 - p));
    sumX += x * PofX; // expectation: sum of X*P(X)
    sumXX += x * x * PofX;

    if (trapezoidAreas[i] > maxArea) {
      maxArea = trapezoidAreas[i];
      maxTrapI = i;
    }
  }

  const float WIDTH_THRESHOLD = 0.000000001;
  float modeP =
      (gridMap.indexToP(maxTrapI) + gridMap.indexToP(maxTrapI + 1)) / 2;
  float width = gridMap.indexToP(maxTrapI + 1) - gridMap.indexToP(maxTrapI);
  while (width > WIDTH_THRESHOLD) {
    // std::cout << "modeP=" << modeP << " width=" << width << "\n";
    const float leftP = modeP - width / 4;
    const float rightP = modeP + width / 4;

    const float centerVal = densityFunction(modeP);
    const float leftVal = densityFunction(leftP);
    const float rightVal = densityFunction(rightP);

    if (centerVal > leftVal && centerVal > rightVal) {
      // modeP is unchanged
    } else if (leftVal > rightVal) {
      modeP = leftP;
    } else {
      modeP = rightP;
    }
    width /= 2;
  }
  const float modeX = modeP; // log(modeP / (1 - modeP));

  // std::cout << "sumX=" << sumX << "\n";
  // std::cout << "sumXX=" << sumXX << "\n";
  EstimatedMoments moments;
  // Mean = E[X]
  moments.mean = sumX;
  // Variance = E[X^2]=E[X]^2
  moments.var = sumXX - sumX * sumX;
  moments.mode = modeX;

  return moments;
}