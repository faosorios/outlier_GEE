# The gradient test statistic for outlier detection in generalized estimating equations

[![R](https://img.shields.io/badge/Made%20with-R%20under%20development-success)](https://cran.r-project.org/)
[![gee](https://img.shields.io/badge/gee-4.13--20-orange)](https://cran.r-project.org/package=gee)

Supplementary material to **The gradient test statistic for outlier detection in generalized estimating equations** by Felipe Osorio, Angelo Garate and Cibele Russo (Statistics & Probability Letters. DOI: [10.1016/j.spl.2024.110087](https://doi.org/10.1016/j.spl.2024.110087))

Code tested on:
- R under development (2022-01-27 r81578), running Linux Zorin 16 (64 bits)
- R version 3.3.0, running OS X 10.13.4 (64 bits)

Attached packages: gee 4.13-20

CONTENTS:
- case_study/case_GUIDE.R: R commands for the analysis of GUIDE dataset (described/analyzed at Section 4 from manuscript).
- case_study/case_LIMA.R: R commands for the analysis of Lima and Sanudo dataset (described/analyzed in Appendix E from supplementary material).
- code/influenceGEE_logit.R: R functions to compute diagnostic measures for GEE (specific for the logistic model).
- code/influenceGEE_norm.R: R functions to compute diagnostic measures for GEE (specific for the gaussian model, used in Appendix E).
- data/guide.rda: urinary incontinence in elderly patients (GUIDE) dataset, in RDA format.
- data/guide.txt: urinary incontinence in elderly patients (GUIDE) dataset, in text format. Extracted from http://www.bios.unc.edu/~preisser/personal/uidata/
- data/lima.rda: Lima and Sanudo dataset, in RDA format.
- data/lima.csv: Lima and Sanudo dataset, in CSV format.
- README.md: this file.
