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
                  const std::vector<float> &pis, float alpha, float beta);

  float logLikelihood(float p) const;

  float operator()(float p) const; // This is *NOT* in log space!

private:
  const std::vector<MatPat> &counts;
  const std::vector<float> &pis; // switching error rate
  float alpha, beta;             // parameters of the beta prior on p
};
#endif
