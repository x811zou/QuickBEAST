/****************************************************************
 LogBinomial.C
 Copyright (C)2023 William H. Majoros (bmajoros@alumni.duke.edu).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#include "LogBinomial.H"
#include <cmath>
#include <iostream>
#include <vector>
using namespace std;

const float LOG0 = log(0);
const float LOG1 = log(1);

float logBinomCoef(int n, int K);

float logBinomialProb(float p, int k, int n) {
  if (p == 0) {
    return k == 0 ? LOG1 : LOG0;
  } else if (p == 1) {
    return k == n ? LOG1 : LOG0;
  }
  const float coef = logBinomCoef(n, k);
  const float betaPart = k * log(p) + (n - k) * log(1 - p);
  const float prob = coef + betaPart;
  return prob;
}

std::vector<float> logFactorial(1);

void ensureLogFactorial(int n) {
  const size_t oldSize = logFactorial.size();
  if (oldSize >= n + 1) {
    return;
  }
  logFactorial.resize(n + 1);
  for (size_t i = oldSize; i < n + 1; i++) {
    logFactorial[i] = logFactorial[i - 1] + logf(i);
  }
}

float logBinomCoef(int n, int K) {
  const int k = K >= (n - K) ? K : n - K;

  ensureLogFactorial(n);
  return logFactorial[n] - logFactorial[k] - logFactorial[n - k];
}
