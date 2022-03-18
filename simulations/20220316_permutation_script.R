library(readr)
library(tidyverse)

set.seed(293878)

cor <- 0.8
n_pleio <- 49317
n_trait_spec <- 49317
h2 <- 0.6

cor_matrix <- matrix(c(1, cor, cor, 1), 2)

for(i in 1:100){
      phenofile <- paste0(paste("100reps_1kgenomes_cor", cor, "npleio", n_pleio, "ntraitspec", n_trait_spec, "h2", h2, sep = "_"), "/phenos_", i, ".txt")
      if(file.exists(phenofile)){
	phenos <- read_delim(phenofile, col_names=TRUE, delim="\t")
     	ids <- phenos[,1:2]
      	permuted_phenos <- phenos[sample(x = 1:nrow(phenos), replace = FALSE),3:4]
        permuted_res <- ids %>% add_column(permuted_phenos)
        write_delim(permuted_res, file = paste0(paste("permuted_100reps_1kgenomes_cor", cor, "npleio", n_pleio, "ntraitspec", n_trait_spec, "h2", h2, sep = "_"), "/permuted_phenos_", i, ".txt"), delim = "\t")
	}
}
