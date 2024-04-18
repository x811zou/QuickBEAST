# QuickBEAST: A comprehensive Tool for Estimating ASE Effect Size

### Overview:
QuickBEAST performs an in-depth statistical analysis to estimate the binomial proportion p representing the ASE. It leverages advanced statistical methods to analyze allele counts and addresses complications like switching error rates and read count heterogeneity. The main output, "qb_mode", is the mode of p estimated through a subgrid search algorithm, ensuring precision in identifying the most likely allele frequency. Subsequently, we can compute the ASE effect size (θ) downstream using the relation θ = p/(1-p), offering a clear quantification of the expression bias between alleles.

### Usage:
##### git clone QuickBEAST
```
git clone git@github.com:x811zou/QuickBEAST.git
cd QuickBEAST/
```
##### run python wrapper script with test data to obtain p values
```
in_file = test_data/bimodal_genes
out_qb_p_file = test_data/bimodal_genes_qb_p
python ./calculate_p_value_from_qb_mode.py --disable-cache $in_file $out_qb_p_file
```
##### run QuickBEAST directly to obtain qb estimates
```
in_file = test_data/bimodal_genes
out_qb_file = test_data/bimodal_genes_qb
./QuickBEAST --alpha 8.789625 --beta 8.789625 --mean --mode -f $in_file --fixMaxHetSite > $out_qb_file
```

### Citing:
Please cite this paper when using QuickBEAST for your publications.
```
[place holder]
```

### Acknowledgements:
The initial version of QuickBEAST was developed by Bill (William H. Majoros). The current version has been enhanced and modified by myself, building upon Bill's original foundational work. We thank Bill for setting the high standards that continue to guide this project. For further information, inquiries, or contributions, please feel free to contact us: [Xue Zou](mailto:xz195@duke.edu), [William H. Majoros](mailto:bmajoros@alumni.duke.edu).
