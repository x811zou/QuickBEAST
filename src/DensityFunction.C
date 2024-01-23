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
                                 float beta, bool fixMaxHetSite)
    : counts(counts), pis(pis), alpha(alpha), beta(beta), fixedSiteIndex(0) {

  const int n = counts.size();

  if (fixMaxHetSite) {
    // find largest het count
    int maxHetCount = 0;
    int maxHetIndex = -1;
    for (int i = 0; i < n; ++i) {
      const int hetCount = counts[i].mat + counts[i].pat;
      if (hetCount > maxHetCount) {
        maxHetCount = hetCount;
        maxHetIndex = i;
      }
    }
    // fix largest read count
    fixedSiteIndex = maxHetIndex;
  }

  // precompute iteration indices
  // start with the fixed site, then go left to right to the right side
  // then finish right to left with the left side
  // a b c d e f
  // fix c => c | d e f | c b a
  for (int j = 1; j < n; ++j) {
    int i, piI;
    // right side first
    if (j < (n - fixedSiteIndex)) {
      i = fixedSiteIndex + j;
      piI = i - 1;
    } else {
      // left side reverse
      i = fixedSiteIndex - (j - (n - fixedSiteIndex)) - 1;
      piI = i;
    }
    i_piI_pairs.push_back({i, piI});
  }

  // std::cerr << "fixedSiteIndex=" << fixedSiteIndex << " n=" << n << '\n';
  // for (int i = 0; i < i_piI_pairs.size(); ++i) {
  //   std::cerr << "i=" << i_piI_pairs[i].first
  //             << " piI=" << i_piI_pairs[i].second << '\n';
  // }

  // precompute log(pi) and log(1-pi)
  for (int i = 0; i < pis.size(); ++i) {
    const float pi = pis[i];
    const float logPi = log(pi);
    const float logOneMinusPi = log(1 - pi);
    logPi_logOneMinusPi_pairs.push_back({logPi, logOneMinusPi});
  }
}

float DensityFunction::operator()(float p) const {
  return exp(logLikelihood(p));
}

float DensityFunction::logLikelihood(float p) const {
  const int n = counts.size();

  // initial
  const int a0 = counts[fixedSiteIndex].mat, b0 = counts[fixedSiteIndex].pat;
  float even =
      log(gsl_ran_beta_pdf(p, alpha, beta)) + logBinomialProb(p, a0, a0 + b0);
  float odd = NEGATIVE_INFINITY;

  // iterate
  for (int j = 0; j < i_piI_pairs.size(); ++j) {
    const int i = i_piI_pairs[j].first, piI = i_piI_pairs[j].second;
    const int ai = counts[i].mat, bi = counts[i].pat;
    const float logPi = logPi_logOneMinusPi_pairs[piI].first;
    const float logOneMinusPi = logPi_logOneMinusPi_pairs[piI].second;

    const float logBinA = logBinomialProb(p, ai, ai + bi);
    const float logBinB = logBinomialProb(p, bi, ai + bi);

    // std::cout << "i=" << i << " mat=" << mat_pat.first
    //           << " pat=" << mat_pat.second << " logBinA=" << logBinA
    //           << " logBinB=" << logBinB << std::endl;

    const float nextEven =
        sumLogProbs(even + logOneMinusPi + logBinA, odd + logPi + logBinA);
    const float nextOdd =
        sumLogProbs(even + logPi + logBinB, odd + logOneMinusPi + logBinB);
    even = nextEven;
    odd = nextOdd;
  }

  // compute final result
  const float logDensity = sumLogProbs(even, odd);
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