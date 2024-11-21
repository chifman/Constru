# CONSTRU - Computing Prognostic Marker Dependencies by Successive Testing of Gene-Stratified Subgroups

### Introduction

This algorithm is designed to enable the empirical discovery of groups of tumors that differ markedly with respect to the power and significance of a given prognostic variable in large tumor gene expression profiling datasets. 

The input is a tumor-gene expression data matrix, survival time and event of each sample, and the corresponding measures of the prognostic variable. To measure how the relative expression of each gene comprising a tumor gene expression data matrix enhances or antagonizes the association between a prognostic variable and patient survival, the normalized expression data for each gene is used to organize tumors into groups based the gene’s relative expression level.

In the current implementation, tumor grouping is based on gene expression tertiles, thereby representing a standardized measure of low, intermediate or high expression. For each gene tertile, a multivariable Cox proportional hazards regression model is fitted to the corresponding patient data, and the significance (p-value) and directionality (hazard ratio) of the association between the prognostic variable and survival outcome is computed. The Cox p-values and hazard ratios is used to generate a combined measure of significance and effect size for each tertile (T1 and T3), called the parity score. Higher parity scores equate with genes that have greater variable-survival associations in T1 (termed LowerT genes), while lower parity scores reflect genes with greater variable-survival associations in T3 (termed UpperT genes).

Parity scores are calculated as follows according to three different computational options:

```
PS1=
      (-log10("cox regression Pvalue of gene Tertile 1")/"cox regression Hazard ratio of gene Tertile 1") - 
      (-log10("cox regression Pvalue of gene Tertile 3")/"cox regression Hazard ratio of gene Tertile 3")
 
PS2=
      (log("cox regression Pvalue of gene Tertile 1")*("cox regression Hazard ratio of gene Tertile 1"-1)) - 
      (log("cox regression Pvalue of gene Tertile 3")*("cox regression Hazard ratio of gene Tertile 3"-1))
 
PS3= 
      (log("cox regression Pvalue of gene Tertile 1")*log("cox regression Hazard ratio of gene Tertile 1"))- 
      (log("cox regression Pvalue of gene Tertile 3")*log("cox regression Hazard ratio of gene Tertile 3"))
```

The parity score is used to rank the genes to distinguish the largest difference between lower (T1) or upper (T3) tertiles with respect to the statistical power of the prognostic variable. The ranked genes can then be assigned percentile ranks as a function of all genes used in the analysis, thus allowing a standardized approach for comparing genes across different datasets. A larger the difference between ranks, the greater the difference between T1 and T3 with respect to prognostic variable significance and effect size.
 
Parity score PS1 is the algorithm used to rank genes in the published article: _A patient stratification signature mirrors the immunogenic potential of high grade serous ovarian cancers. [Journal of Translational Medicine](https://doi.org/10.1186/s12967-024-05846-9)_. See the article for methods on selecting genes by reproducibility of parity scores, and the establishment of gene signatures derived from selected genes.

**Caveats**

Genes that have expression patterns correlated with the variable are poor candidates as the goal is to stratify patients based on relationships between gene expression levels and the prognostic variable of interest. In such cases, the tumor groups defined by the resulting gene tertiles will exhibit a compressed distribution of the prognostic variable values.

Genes with near zero correlation with the variable are preferred. Therefore, either before or after running CONSTRU, it is advisable to filter genes based on their inherent positive or negative correlations with the prognostic variable. In the published article, an empirical threshold of Pearson’s correlation coefficient of <-0.15 or >0.15 was used to exclude genes from analysis.

**Alternative method**

The central function of CONSTRU is akin to the method of modeling each gene’s effect on the prognostic variable as an interaction term in a Cox regression model. We have included this approach in the algorithm for comparison with CONSTRU gene parity score rankings.

In the original method, the calculation of the parity score is based on the comparison of Cox results for T1 and T3, and any effect relative to T2 is not considered. The alternative method does not exclude the T2 data and may result in a markedly different ranking for a given gene, particularly in cases where a significant association observed in CONSTRU T3 extends into T2.

Depending on the preferred gene classification parameters, such genes may be more desirable for signature construction.

### Installation

```{r}

# Install from local repository
# Step 1: download the package
# Step 2: install local package

library("devtools")
install(path_to_repository, dependencies = TRUE)

# We plan to make the package available on a public Github repository.
# When it is available, it can be installed with the commands:

install_github("chifman/Constru")

```

### Input

1. Patient survival data
    A data frame with survival data with the samples as rows and survival time and event as columns.

2. Gene expression data. 
    A data frame with normalized gene expression measurements or counts (eg, CPM/RPKM), where the columns are sample names and rows are genes.

3. Prognostic variable
    A vector with the prognostic variable values. It should be of the same length and order as the row of the survival data.

4. Formula used for cox regression. Example: 
    Surv( OS_TIME , OS_Event ) ~ prognostic_variable 

If additional prognostic factors need to be added to the model, add it as a column in the survival data and use its column name in the formula. Example:
    Surv( OS_TIME , OS_Event ) ~ prognostic_variable + AGE 

### Vignette

```{r}

# The sample data and vignette code is located in the vignette folder. 

library("CONSTRU")

gene_data=read.table("Ovarian_cancer_gene_expression.tsv")
survival_data=read.table("Ovarian_cancer_survival_data.tsv")
prognostic_variable=read.table("Prognostic_variable.tsv")
cox_formula="Surv( OS_Time , OS_Event ) ~ prognostic_variable";
ncores=parallel::detectCores()

# format the input so that the order is the same
survival_data=survival_data[rownames(survival_data),]
prognostic_variable_d=unlist(prognostic_variable[rownames(survival_data),])


# Original method
output1=constru(survival_data, gene_data, prognostic_variable_d,cox_formula,ncores)

# Alternative method
output2=constru_continuous(survival_data, gene_data, prognostic_variable_d,cox_formula,ncores)

```

### Citation

**_Laurel K. Berry, Ashok K. Pullikuth, Kristen L. Stearns, Yuezhu Wang, Calvin J. Wagner, Jeff W. Chou, Janelle P. Darby, Michael G. Kelly, Raghvendra Mall, Ming Leung, Julia Chifman and Lance D. Miller. 2024. A patient stratification signature mirrors the immunogenic potential of high grade serous ovarian cancers. [Journal of Translational Medicine](https://doi.org/10.1186/s12967-024-05846-9). PMID: [39568014](https://pubmed.ncbi.nlm.nih.gov/39568014/)_**

