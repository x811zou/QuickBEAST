/****************************************************************
 DensityFunction.C
 Copyright (C)2023 William H. Majoros (bmajoros@alumni.duke.edu).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#include "DensityFunction.H"

#include "LogBinomial.H"
#include <cmath>
#include <gsl/gsl_randist.h>
#include <utility>

const float NEGATIVE_INFINITY = -INFINITY;
float sumLogProbs(float logP, float logQ);

DensityFunction::DensityFunction(const std::vector<MatPat> &counts,
                                 const std::vector<float> &pis, float alpha,
                                 float beta)
    : counts(counts), pis(pis), alpha(alpha), beta(beta) {
  // ctor
}

float DensityFunction::operator()(float p) const {
  return exp(logLikelihood(p));
}

float DensityFunction::logLikelihood(float p) const {
  const int n = counts.size();

  // initial
  const int a0 = counts[0].mat, b0 = counts[0].pat;
  std::pair<float, float> mat_pat = std::make_pair(
      log(gsl_ran_beta_pdf(p, alpha, beta)) + logBinomialProb(p, a0, a0 + b0),
      NEGATIVE_INFINITY);

  // iterate
  for (int i = 1; i < n; ++i) {
    const int ai = counts[i].mat, bi = counts[i].pat;
    const float pi = pis[i - 1];
    const float logBinA = logBinomialProb(p, ai, ai + bi);
    const float logBinB = logBinomialProb(p, bi, ai + bi);
    float logPi = log(pi), logOneMinusPi = log(1 - pi);

    // std::cout << "i=" << i << " mat=" << mat_pat.first
    //           << " pat=" << mat_pat.second << " logBinA=" << logBinA
    //           << " logBinB=" << logBinB << std::endl;

    mat_pat =
        std::make_pair(sumLogProbs(mat_pat.first + logOneMinusPi + logBinA,
                                   mat_pat.second + logPi + logBinA),
                       sumLogProbs(mat_pat.first + logPi + logBinB,
                                   mat_pat.second + logOneMinusPi + logBinB));
  }

  // compute final result
  const float logDensity = sumLogProbs(mat_pat.first, mat_pat.second);
  return logDensity;
}

float sumLogProbs(float logP, float logQ) {
  if (logP == NEGATIVE_INFINITY) {
    return logQ;
  } else if (logQ == NEGATIVE_INFINITY) {
    return logP;
  } else if (logP > logQ) {
    return logP + log1p(exp(logQ - logP));
  } else {
    return logQ + log1p(exp(logP - logQ));
  }
}