# QuickBEAST: A comprehensive Tool for Estimating ASE Effect Size

### Overview:
QuickBEAST employs the same statistical model as [BEASTIE](https://github.com/x811zou/BEASTIE), but it estimates the binomial proportion p representing the ASE. The main output, "qb_mode", is the MAP effect size (mode of p) estimated through a subgrid search algorithm. Subsequently, we can compute the p values from 1000 null simulation for each gene. In addition, we can also compute the ASE effect size (θ) downstream using the relation θ = p/(1-p), offering a clear quantification of the expression bias between alleles.

### Build QuickBEAST
```
place holder
```

### Usage:
##### run python wrapper script with test data to obtain p values
```
in_file = test_data/bimodal_genes
out_qb_p_file = test_data/bimodal_genes_qb_p
python ./calculate_p_value_from_qb_mode.py --disable-cache $in_file $out_qb_p_file
```
##### run QuickBEAST directly with test data to obtain qb estimates
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
The initial version of QuickBEAST was developed by Bill (William H. Majoros). The current version has been enhanced and modified by myself, building upon Bill's original foundational work. We thank Bill for setting the high standards that continue to guide this project. For further information, inquiries, or contributions, please feel free to contact [Xue Zou](mailto:xz195@duke.edu).
