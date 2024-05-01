# CONSTRU - Computing Prognostic Marker Dependencies by Successive Testing of Gene-Stratified Subgroups

### Introduction

This algorithm measures how the relative expression of each gene enhances or antagonizes the association between a prognostic marker and patient survival.

Three different strategies were used.

**Strategy 1: "constru_continuous"**

The simplest method is to include the association as a interaction term in the cox regression. 

**Strategy 2: "constru"**

Cox regression assumes that continuous predictors have a linear relationship with the outcome. However, multiple studies have shown that the prognostic and predictive power of immune signatures may be restricted to distinct tumor subgroups with favorable immunogenic properties. The discrete partitions of the subgroups often result in non-continuous relationships.

As such, the second method separates the samples by gene expression tertile before performing cox regression. The results of the cox regression of the first and third tertile is used to calculate a parity score to quantify the effect. The three parity scores calculated are as follows:

PS1=\
      (-log10("cox regression Pvalue of gene Tertile 1")/"cox regression Hazard ratio of gene Tertile 1") - \
      (-log10("cox regression Pvalue of gene Tertile 3")/"cox regression Hazard ratio of gene Tertile 3")
 
PS2=\
      (log("cox regression Pvalue of gene Tertile 1")*("cox regression Hazard ratio of gene Tertile 1"-1)) - \
      (log("cox regression Pvalue of gene Tertile 3")*("cox regression Hazard ratio of gene Tertile 3"-1))
 
PS3= \
      (log("cox regression Pvalue of gene Tertile 1")*log("cox regression Hazard ratio of gene Tertile 1"))- \
      (log("cox regression Pvalue of gene Tertile 3")*log("cox regression Hazard ratio of gene Tertile 3"))
 
Parity score PS1 is used to rank the genes by its effect on the prognostic marker in the paper _A patient stratification signature mirrors the immunogenic potential of high grade serous ovarian cancers. [Journal](link)_. Higher parity scores equate with genes that have greater Prognostic variable - survival associations in the first tertile, whereas lower parity scores equate with genes that have greater Prognostic variable - survival associations in the third tertile.

**Strategy 3: "constru_kmeans"**

Separating by tertiles is a simple way to categorize genes by high and low expression. However, the gene's behavior might now follow these specific divisions. 

To remove this assumption, kmeans clustering is used to determine the division between a set number of groups. The exact number of groups can be specified using the option ngroups.

**Strategy 4: "constru_Mclust"**

In strategies 2 and 3, the programs assumed that a fixed number of groups.

To remove this assumption, R package MClust (Model-based clustering based on parameterized finite Gaussian mixture models) is used to categorize the prognostic variable and gene expression into clusters. The parity scores above are calculated by selecting two of the groups identified by MClust; the selection prioritizes maximizing the parity score. 

Note: This algorithm is very computationally intensive compared to the other methods.

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
    A data frame with normalized gene counts (CPM/RPKM), where the columns are sample names and rows are genes.

3. Prognostic variable
    A vector with the prognostic variable data. It should be of the same length and order as the row of the survival data.

4. Formula used for cox regression. Example: 
    Surv( OS_TIME , OS_Event ) ~ prognostic_variable 

If additional prognostic factors need to be added to the model, add it as a column in the survival data and use its column name in the formula. Example:
    Surv( OS_TIME , OS_Event ) ~ prognostic_variable + AGE 

### Vignette

```{r}

# The sample data and vignette code is located in the vignette folder. 

library(CONSTRU)

```

### Citation

**_Laurel K. Berry, Ashok K. Pullikuth, Kristen L. Stearns, Yuezhu Wang, Calvin J. Wagner, Jeff W. Chou, Janelle P. Darby, Michael G. Kelly, Raghvendra Mall, Ming Leung, Julia Chifman and Lance D. Miller. 2024. A patient stratification signature mirrors the immunogenic potential of high grade serous ovarian cancers. [Journal](link). PMID: [number](link)_**

