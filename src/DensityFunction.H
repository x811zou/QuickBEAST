/****************************************************************
 DensityFunction.H
 Copyright (C)2023 William H. Majoros (bmajoros@alumni.duke.edu).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#ifndef INCL_DensityFunction_H
#define INCL_DensityFunction_H

#include "MatPat.H"
#include <iostream>
#include <vector>

class DensityFunction {
public:
  DensityFunction(const std::vector<MatPat> &counts,
                  const std::vector<float> &pis, float alpha, float beta,
                  bool fixMaxHetSite);

  float logLikelihood(float p) const;

  float operator()(float p) const; // This is *NOT* in log space!

private:
  const std::vector<MatPat> &counts;
  const std::vector<float> &pis; // switching error rate
  float alpha, beta;             // parameters of the beta prior on p

  int fixedSiteIndex;

  // cached data
  std::vector<std::pair<int, int>> i_piI_pairs;
  std::vector<std::pair<float, float>> logPi_logOneMinusPi_pairs;
};
#endif
