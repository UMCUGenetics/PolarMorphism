# PolarMorphism
## Discovery of genetic variants with shared effect across multiple phenotypes from GWAS summary statistics

Preprint manuscript describing this package: https://www.biorxiv.org/content/10.1101/2022.01.14.476302v1

PolarMorphism uses genome-wide summary statistics from GWAS of multiple related phenotypes to output statistics per SNP that describe its degree of 'sharedness' across the phenotypes of interest and its overall (multivariate) effect size, as well as p-values indicating significance for each statistic. The method is based on a transform from Cartesian to polar coordinates; the 'sharedness' is indicated by the angle with the Cartesian x-axis *theta* and the overall effect size is indicated by the distance from the origin *r*.

See "template_for_github.Rmd" for an example R Notebook of how to use the package.

Necessary input:
- Genome-wide summary statistics from two or more GWAS. Even if you only want to know the sharedness of a few SNPs you need to supply genome-wide summary statistics because PolarMorhpism needs to calculate the covariance matrix from SNPs across the whole genome.

# Installation
```{r}
if(!require("devtools")){
  install.packages("devtools")
}
library(devtools)
install_github("UMCUgenetics/PolarMorphism")
library(PolarMorphism)
```

# Glossary
GWAS
: Genome-Wide Association Study/Studies

