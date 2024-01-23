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
#include <iomanip>
#include <iostream>
#include <math.h>
#include <sstream>

using namespace std;

struct Gene {
  std::string name;
  std::vector<MatPat> counts;
  std::vector<float> pis;
};

Gene parseGeneFromTextLine(const std::string &line);

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
  // clang-format off
  options.add_options()
  ("alpha", "Alpha for beta distribution", cxxopts::value<float>()->default_value("8.789625"))
  ("beta", "Beta for beta distribution", cxxopts::value<float>()->default_value("8.789625"))

  ("f,file", "Input file, if omitted read from stdin", cxxopts::value<std::string>()->default_value(""))

  ("mode", "Compute mode")
  ("modeSubgridSize", "mode subgrid computation points", cxxopts::value<int>()->default_value("51"))
  ("modeSubgridThreshold", "final subgrid width for mode computation", cxxopts::value<float>()->default_value("1e-3"))
  ("modeSubgridStepSlop", "mode subgrid step slop", cxxopts::value<float>()->default_value("3"))

  ("mean", "Compute mean and variance")
  ("meanTrapezoids", "Number of trapezoids for mean computation. Must be even", cxxopts::value<int>()->default_value("100"))

  ("fixMaxHetSite", "Fix the site with the largest number of heterozygous reads")

  ("distributionFile", "File to optionally store posterior to", cxxopts::value<std::string>()->default_value(""))
  ("distributionLogLikelihood", "Compute distribution log likelihood instead of density")
  ("distributionPoints", "Number of points to use for distribution", cxxopts::value<int>()->default_value("10000"))

  ("v,verbose", "Verbose output")
  ("h,help", "Print help");
  // clang-format on

  auto result = options.parse(argc, argv);

  if (result.count("help")) {
    std::cout << options.help() << std::endl;
    exit(0);
  }

  const float alpha = result["alpha"].as<float>();
  const float beta = result["beta"].as<float>();
  const std::string inputFile = result["file"].as<std::string>();

  const bool computeMeanAndVariance = result.count("mean") > 0;
  const int meanTrapezoids = result["meanTrapezoids"].as<int>();

  const bool computeMode = result.count("mode") > 0;
  const int modeSubgridSize = result["modeSubgridSize"].as<int>();
  const float modeSubgridThreshold = result["modeSubgridThreshold"].as<float>();
  const float modeSubgridStepSlop = result["modeSubgridStepSlop"].as<float>();

  const bool fixMaxHetSite = result.count("fixMaxHetSite") > 0;

  const std::string distributionFile =
      result["distributionFile"].as<std::string>();
  const bool distributionLogLikelihood =
      result.count("distributionLogLikelihood") > 0;
  const int distributionPoints = result["distributionPoints"].as<int>();

  bool verbose = result.count("verbose") > 0;

  if (meanTrapezoids % 2 != 0) {
    std::cerr << "ERROR: meanTrapezoids must be an even number" << std::endl;
    return 1;
  }

  if (verbose) {
    // print parsed args
    cerr << "alpha=" << alpha << '\n';
    cerr << "beta=" << beta << '\n';

    cerr << "computeMode=" << computeMode << '\n';
    cerr << "modeSubgridSize=" << modeSubgridSize << '\n';
    cerr << "modeSubgridThreshold=" << modeSubgridThreshold << '\n';
    cerr << "modeSubgridStepSlop=" << modeSubgridStepSlop << '\n';
    cerr << "meanTrapezoids=" << meanTrapezoids << '\n';

    cerr << "computeMean=" << computeMeanAndVariance << '\n';
    cerr << "distributionPoints=" << distributionPoints << '\n';
    cerr << "distributionLogLikelihood=" << distributionLogLikelihood << '\n';
    cerr << "inputFile="
         << "'" << inputFile << "'" << '\n';
    cerr << "distributionFile="
         << "'" << distributionFile << "'" << '\n';
  }

  const bool writeOutput = computeMeanAndVariance || computeMode;
  const bool writeDistribution = !distributionFile.empty();

  if (!writeOutput && !writeDistribution) {
    std::cerr << "ERROR: No output requested" << std::endl;
    return 1;
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

  std::ofstream distributionFileStream;
  if (writeDistribution) {
    if (verbose) {
      cerr << "Writing distribution to file: " << distributionFile << endl;
    }
    distributionFileStream.open(distributionFile.c_str());
    if (!distributionFileStream.is_open()) {
      std::cerr << "Failed to open file: " << distributionFile << std::endl;
      return 1;
    }
  }

  // print header
  if (writeOutput) {
    cout << "gene";
    if (computeMeanAndVariance) {
      cout << "\tmean\tvariance";
    }
    if (computeMode) {
      cout << "\tmode";
    }
    cout << "\n";

    // configure float output
    cout << std::fixed << std::setprecision(12);
  }
  std::string line;
  while (std::getline(*input, line)) {
    Gene gene = parseGeneFromTextLine(line);

    DensityFunction f(gene.counts, gene.pis, alpha, beta, fixMaxHetSite);

    if (writeOutput) {
      std::cout << gene.name;

      if (computeMeanAndVariance) {
        const auto &meanAndVariance =
            estimateMeanAndVariance({meanTrapezoids}, f);
        std::cout << "\t" << meanAndVariance.mean << "\t"
                  << meanAndVariance.variance;
      }
      if (computeMode) {
        const float mode = estimateMode(
            {modeSubgridSize, modeSubgridThreshold, modeSubgridStepSlop}, f);
        std::cout << "\t" << mode;
      }

      std::cout << "\n";
    }

    if (writeDistribution) {
      distributionFileStream << gene.name << "\t";
      const auto &&distribution = computeDistribution(
          {distributionPoints, distributionLogLikelihood}, f);
      for (const auto &pair : distribution) {
        distributionFileStream << pair.first << "\t" << pair.second << "\t";
      }
      distributionFileStream << "\n";
    }
  }

  return 0;
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