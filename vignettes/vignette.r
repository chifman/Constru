#!/usr/bin/env Rscript

library("devtools")
install("../", dependencies = TRUE)

library("CONSTRU")

gene_data=read.table("Ovarian_cancer_gene_expression.tsv")
survival_data=read.table("Ovarian_cancer_survival_data.tsv")
prognostic_variable=read.table("Prognostic_variable.tsv")
cox_formula="Surv( OS_Time , OS_Event ) ~ prognostic_variable";
ncores=parallel::detectCores()

# format the input so that the order is the same
survival_data=survival_data[rownames(survival_data),]
prognostic_variable_d=unlist(prognostic_variable[rownames(survival_data),])


# Strategy 1 
output1=constru(survival_data, gene_data, prognostic_variable_d,cox_formula,ncores)

# Strategy 2
output2=constru_continuous(survival_data, gene_data, prognostic_variable_d,cox_formula,ncores)

# Strategy 3
n_kmeans_groups=3
output3=constru_kmeans(survival_data, gene_data, prognostic_variable_d,cox_formula,ncores,n_kmeans_groups)

# Strategy 4
output4=constru_Mclust(survival_data, gene_data, prognostic_variable_d,cox_formula,ncores)

