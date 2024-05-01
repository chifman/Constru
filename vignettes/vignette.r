#!/usr/bin/env Rscript

library("devtools")
install("../", dependencies = TRUE)

library("CONSTRU")


gene_data=read.table("Ovarian_cancer_gene_expression.tsv")
survival_data=read.table("Ovarian_cancer_survival_data.tsv")
prognostic_variable=read.table("Prognostic_variable.tsv")

prognostic_variable_d=unlist(prognostic_variable[rownames(survival_data),])
ncores=parallel::detectCores()
cox_formula="Surv( OS_Time , OS_Event ) ~ prognostic_variable";

out=constru(survival_data, gene_data, prognostic_variable_d,cox_formula,ncores)

