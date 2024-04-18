## QuickBEAST: A comprehensive Tool for Estimating ASE Effect Size

#### QuickBEAST: A Comprehensive Tool for Estimating ASE Effect Size

##### Overview:
QuickBEAST is an open-source software tool designed for the precise estimation of allele-specific expression (ASE) effect size (θ) from gene data. The software employs advanced statistical methods to analyze allele counts, accounting for nuances such as switching error rates and heterogeneity in read counts.

##### Key Features:

Robust Statistical Framework: At its core, QuickBEAST utilizes a DensityFunction class to model the density of allele frequencies, offering functionalities to compute log-likelihoods and evaluate density functions.
Precise Estimations: The MomentsEstimator module extends the capabilities by providing accurate estimations of mean, variance, and mode, crucial for interpreting ASE effect sizes.
Advanced Numerical Methods: Implements trapezoidal integration (via the Trapezoids module) and efficient log-binomial probability calculations (through the LogBinomial module), ensuring high precision in numerical computations.
User-Friendly Command-Line Interface: QuickBEAST is accessible via a command-line interface, offering various parameters for detailed configuration and analysis. It supports input from files or standard input, and provides verbose output and file writing options for the results.
Methodology:
QuickBEAST performs an in-depth statistical analysis to estimate the binomial proportion p representing the ASE. The main output, qb_mode, is the mode of p estimated through a subgrid search algorithm, ensuring precision in identifying the most likely allele frequency. Subsequently, we can compute the ASE effect size (θ) downstream using the relation θ = p/(1-p), offering a clear quantification of the expression bias between alleles.

##### Usage:
The software is designed with flexibility in mind, allowing users to specify various options such as the alpha and beta parameters for the beta distribution, the choice of output (mean, variance, mode), and the ability to fix the site with the largest number of heterozygous reads, among others.

##### Conclusion:
QuickBEAST stands out for its robust statistical foundation, precision in computation, and user-centric design. It is a powerful tool for researchers and geneticists aiming to unravel the complexities of allele-specific expression and its implications in genomic studies.

##### Acknowledgements:
The initial version of QuickBEAST was developed by Bill (William H. Majoros). The current version has been enhanced and modified by myself, building upon Bill's original foundational work. We thank Bill Majoros for setting the high standards that continue to guide this project. For further information, inquiries, or contributions, please feel free to contact us: [Xue Zou](mailto:xz195@duke.edu), [William H. Majoros](mailto:bmajoros@alumni.duke.edu).
