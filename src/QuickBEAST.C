/****************************************************************
 QuickBEAST.C
 Copyright (C)2023 William H. Majoros (bmajoros@alumni.duke.edu).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#include "DensityFunction.H"
#include "DensityGrid.H"
#include "GridMap.H"
#include "MatPat.H"
#include "MomentsEstimator.H"
#include "Trapezoids.H"
#include "cxxopts.hpp"
#include <fstream>
#include <iostream>
#include <math.h>
#include <sstream>

using namespace std;

#define XCC

struct Gene {
  std::string name;
  std::vector<MatPat> counts;
  std::vector<float> pis;
};

struct ComputedGeneData {
  float mean;
  float var;
  float mode;
};

Gene parseGeneFromTextLine(const std::string &line);
ComputedGeneData computeGeneData(const Gene &gene, float alpha, float beta,
                                 int numPoints);

class Application {
public:
  Application();
  int main(int argc, char *argv[]);
};

int main(int argc, char *argv[]) {
  try {
    Application app;
    return app.main(argc, argv);
  } catch (const char *p) {
    cerr << p << endl;
  } catch (const std::string &msg) {
    cerr << msg.c_str() << endl;
  } catch (const exception &e) {
    cerr << "STL exception caught in main:\n" << e.what() << endl;
  } catch (...) {
    cerr << "Unknown exception caught in main" << endl;
  }
  return -1;
}

Application::Application() {
  // ctor
}

int Application::main(int argc, char *argv[]) {
  cxxopts::Options options("QuickBEAST", "DESCRIPTION");
  options.add_options()("alpha", "Alpha for beta distribution",
                        cxxopts::value<float>())(
      "beta", "Beta for beta distribution", cxxopts::value<float>())(
      "numTrapezoids",
      "Number of trapezoids for numerical integration.  Must be even",
      cxxopts::value<int>())("v,verbose", "Verbose output")(
      "f,file", "Input file",
      cxxopts::value<std::string>()->default_value(""))("h,help", "Print help");
  options.parse_positional({"alpha", "beta", "numTrapezoids"});

  auto result = options.parse(argc, argv);

  if (result.count("help")) {
    std::cout << options.help() << std::endl;
    exit(0);
  }

  const float alpha = result["alpha"].as<float>();
  const float beta = result["beta"].as<float>();
  const int numTrapezoids = result["numTrapezoids"].as<int>();
  std::string inputFile = result["file"].as<std::string>();
  bool verbose = result.count("verbose") > 0;

  if (numTrapezoids % 2 != 0) {
    throw std::string("#trapezoids must be an even number");
  }
  const int numPoints = numTrapezoids + 1;

  if (verbose) {
    // print parsed args
    cerr << "alpha=" << alpha << endl;
    cerr << "beta=" << beta << endl;
    cerr << "#trapezoids=" << numTrapezoids << endl;
    cerr << "inputFile="
         << "'" << inputFile << "'" << endl;
  }

  std::istream *input = &std::cin;
  std::ifstream file;

  if (!inputFile.empty()) {
    if (verbose) {
      cerr << "Reading from file: " << inputFile << endl;
    }
    file.open(inputFile.c_str());
    if (!file.is_open()) {
      std::cerr << "Failed to open file: " << inputFile << std::endl;
      return 1;
    }
    input = &file;
  } else {
    if (verbose) {
      cerr << "Reading from stdin" << endl;
    }
  }

  // print header
  cout << "gene\tmean\tvar\tmode\n";
  std::string line;
  while (std::getline(*input, line)) {
    Gene gene = parseGeneFromTextLine(line);

    ComputedGeneData computedData =
        computeGeneData(gene, alpha, beta, numPoints);

    cout << gene.name << "\t" << computedData.mean << "\t" << computedData.var
         << "\t" << computedData.mode << endl;
  }

  return 0;
}

ComputedGeneData computeGeneData(const Gene &gene, float alpha, float beta,
                                 int numPoints) {
  DensityFunction f(gene.counts, gene.pis, alpha, beta);
  GridMap gridMap(numPoints + 1, 0, 1);
  const auto &&grid = fillGrid(f, gridMap);
  const auto &&trapezoidAreas = computeTrapezoidAreas(grid);

  EstimatedMoments moments = estimateMoments(trapezoidAreas, gridMap, f);
  return {.mean = moments.mean, .var = moments.var, .mode = moments.mode};
}

Gene parseGeneFromTextLine(const std::string &line) {
  // format is tab separated
  // geneID n_hets mat1 pat1 mat2 pat2 ... matn patn pi1 pi2 ... pin-1
  // where n_hets is the number of heterozygous sites
  // and pi1, pi2, ..., pin-1 are the switching error rates

  // split line into tokens
  std::vector<std::string> tokens;
  std::string token;
  std::istringstream tokenStream(line);
  while (std::getline(tokenStream, token, '\t')) {
    tokens.push_back(token);
  }

  // check that there are at least 4 tokens
  if (tokens.size() < 4) {
    std::stringstream ss;
    ss << "Unable to parse line (too few tokens) [" << tokens.size()
       << " < 4] : " << line;
    throw ss.str();
  }

  // validate length
  const int n_hets = std::atoi(tokens[1].c_str());
  const int n_tokens = tokens.size();
  const int expected_length =
      /* header */ 2 + /* read counts */ 2 * n_hets + /* pis */ (n_hets - 1);
  if (n_tokens != expected_length) {
    std::stringstream ss;
    ss << "Unable to parse line (wrong number of tokens) [found: " << n_tokens
       << " expected: " << expected_length << "] : " << line;
    throw ss.str();
  }

  Gene gene;
  gene.name = tokens[0];
  gene.counts.resize(n_hets);
  gene.pis.resize(n_hets - 1);

  // parse read counts
  for (int i = 0; i < n_hets; ++i) {
    gene.counts[i] = {std::atoi(tokens[2 + 2 * i].c_str()),
                      std::atoi(tokens[2 + 2 * i + 1].c_str())};
  }

  // parse switching error rates
  for (int i = 0; i < n_hets - 1; ++i) {
    gene.pis[i] = std::atof(tokens[2 + 2 * n_hets + i].c_str());
    if (gene.pis[i] < 0.0 || gene.pis[i] > 1.0) {
      std::stringstream ss;
      ss << "Unable to parse line (pi out of range) [pi=" << gene.pis[i]
         << "] : " << line;
      throw ss.str();
    }
  }

  // std::cout << "Parsed gene: " << gene.name << " n_hets=" << n_hets
  //           << " pis=[ ";
  // for (int i = 0; i < n_hets - 1; ++i) {
  //   std::cout << gene.pis[i] << " ";
  // }
  // std::cout << "]\n";

  return gene;
}