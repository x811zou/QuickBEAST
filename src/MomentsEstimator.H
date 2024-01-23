/****************************************************************
 MomentsEstimator.H
 Copyright (C)2023 William H. Majoros (bmajoros@alumni.duke.edu).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#ifndef INCL_MomentsEstimator_H
#define INCL_MomentsEstimator_H

#include <utility>
#include <vector>

class DensityFunction;

struct MeanOptions {
  int numTrapezoids;
};
struct MeanAndVariance {
  float mean;
  float variance;
};
MeanAndVariance estimateMeanAndVariance(const MeanOptions &options,
                                        const DensityFunction &densityFunction);

struct ModeOptions {
  int subgridSize;
  float subgridThreshold;
  float subgridStepSlop;
};
float estimateMode(const ModeOptions &options,
                   const DensityFunction &densityFunction);

struct DistributionOptions {
  int numSamples;
  bool useLogLikelihood;
};
std::vector<std::pair<float, float>>
computeDistribution(const DistributionOptions &options,
                    const DensityFunction &densityFunction);

#endif
