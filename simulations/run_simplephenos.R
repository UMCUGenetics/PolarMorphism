library(simplePHENOTYPES)

set.seed(293878)
nsnps <- 493175 # from wc -l of the vcf - number of header lines
nreps <- 100
wanted_cor <- 0.8
prop_pleio <- 0.1
n_pleio <- floor(prop_pleio*nsnps)
prop_trait_spec <- 0.1
n_trait_spec <- floor(prop_trait_spec*nsnps)
h2 <- 0.6

cor_matrix <- matrix(c(1, wanted_cor, wanted_cor, 1), 2)
effects_pleio <- rnorm(n = n_pleio, mean = rnorm(n = n_pleio, mean = 0, sd = 1), sd = 1)
effects_traitspec <- rnorm(n = n_trait_spec, mean = rnorm(n = n_trait_spec, mean = 0, sd = 1), sd = 1)
effects <- c(effects_pleio, effects_traitspec)

create_phenotypes(
                                                geno_file = "../1Kgenomes/ALL.chip.omni_broad_sanger_combined.20140818.snps.genotypes_PASS_AFgt5perc_AFst95perc_snplist_nonFIN_EUR.vcf",
                                                ntraits = 2,
                                                pleio_a = n_pleio,
                                                trait_spec_a_QTN_num = c(n_trait_spec, n_trait_spec),
                                                h2 = c(h2, h2),
                                                add_effect = list(trait_1 = effects, trait_2 = effects),
                                                cor = cor_matrix,
                                                rep = nreps,
                                                output_dir = paste("100reps_1kgenomes_cor", wanted_cor, "npleio", n_pleio, "ntraitspec", n_trait_spec, "h2", h2, sep = "_"),
                                                output_format = "gemma",
                                                architecture = "partially",
                                                out_geno = "BED",
                                                to_r = FALSE,
                                                model = "A",
                                                home_dir = getwd(),
                                                quiet = T,
                                                sim_method = "custom")