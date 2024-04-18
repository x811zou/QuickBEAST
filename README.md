# QuickBEAST: A comprehensive Tool for Estimating ASE Effect Size

### Overview:
QuickBEAST is an open-source software designed to estimate allele-specific expression (ASE) effect size (θ) accurately. It leverages advanced statistical methods to analyze allele counts and addresses complications like switching error rates and read count heterogeneity.

### Methodology:
QuickBEAST performs an in-depth statistical analysis to estimate the binomial proportion p representing the ASE. The main output, "qb_mode", is the mode of p estimated through a subgrid search algorithm, ensuring precision in identifying the most likely allele frequency. Subsequently, we can compute the ASE effect size (θ) downstream using the relation θ = p/(1-p), offering a clear quantification of the expression bias between alleles.

### Methodology
QuickBEAST estimates the binomial proportion p representing ASE, using a subgrid search algorithm to find the most probable allele frequency. The effect size θ is then calculated as θ = p/(1-p).

### Usage:
```

```

### Citing:
Please cite this paper when using QuickBEAST for your publications.
```
[place holder]
```

### Acknowledgements:
The initial version of QuickBEAST was developed by Bill (William H. Majoros). The current version has been enhanced and modified by myself, building upon Bill's original foundational work. We thank Bill for setting the high standards that continue to guide this project. For further information, inquiries, or contributions, please feel free to contact us: [Xue Zou](mailto:xz195@duke.edu), [William H. Majoros](mailto:bmajoros@alumni.duke.edu).
