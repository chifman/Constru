% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CONSTRU_function.r
\name{constru}
\alias{constru}
\title{constru}
\usage{
constru(survival_data, gene_data, prognostic_variable_d, cox_formula, ncores)
}
\arguments{
\item{survival_data}{A data frame with the samples as rows and survival time and event as columns.}

\item{gene_data}{A data frame sample names as columns and gene expression as rows.}

\item{prognostic_variable_d}{A vector with the prognostic variable data. It should be of the same length and order as the row of the survival data.}

\item{cox_formula}{Formula used for cox regression. Example: 
\itemize{
  \item Surv( OS_TIME , OS_Event ) ~ prognostic_variable \cr
}
  If additional prognostic factors need to be added to the model, add it as a column in the survival data and use its column name in the formula. Example:
\itemize{
  \item Surv( OS_TIME , OS_Event ) ~ prognostic_variable + AGE \cr
}}

\item{ncores}{The number of cores used during multithreading.}
}
\value{
a table with columns:
\itemize{
  \item mean - gene expression mean
  \item P5\% - 5\% confidence interval of the gene expression
  \item P95\% - 95\% confidence interval of the gene expression
  \item Diff P95\%_5\% - difference between the 95\% and 5\% confidence interval
  \item R - pearson correlation between the gene expression and prognostic variable
  \item PS1 - parity score 1 \cr
(-log10("Pvalue gene Tertile 1")/"HR gene Tertile 1") - \cr
(-log10("Pvalue gene Tertile 3")/"HR gene Tertile 3")   \cr
  \item PS2 - parity score 2 \cr
(log("Pvalue gene Tertile 1")*("HR gene Tertile 1"-1)) - \cr
(log("Pvalue gene Tertile 3")*("HR gene Tertile 3"-1))   \cr
  \item PS3 - parity score 3 \cr
(log("Pvalue gene Tertile 1")*log("HR gene Tertile 1"))- \cr
(log("Pvalue gene Tertile 3")*log("HR gene Tertile 3"))  \cr
  \item gT1 Pval - Pvalue of gene Tertile 1
  \item gT1 HR - Hazard ratio of gene Tertile 1
  \item gT1sT1 - number of samples in \cr
  \item gT1sT1 3 year - 3 year survival \cr
      Kaplan-Meier survival curve of the group above
  \item gT1sT1 6 year - 6 year survival \cr
      Kaplan-Meier survival curve of the group above
  \item gT1sT2 - number of samples in \cr
  \item gT1sT2 3 year - 3 year survival \cr
      Kaplan-Meier survival curve of the group above
  \item gT1sT2 6 year - 6 year survival \cr
      Kaplan-Meier survival curve of the group above
  \item gT1sT3 - number of samples in \cr
  \item gT1sT3 3 year - 3 year survival \cr
      Kaplan-Meier survival curve of the group above
  \item gT1sT3 6 year - 6 year survival \cr
      Kaplan-Meier survival curve of the group above
  \item gT2 Pval - Pvalue of gene Tertile 2
  \item gT2 HR - Hazard ratio of gene Tertile 2
  \item gT2sT1 - number of samples in \cr
  \item gT2sT1 3 year - 3 year survival \cr
      Kaplan-Meier survival curve of the group above
  \item gT2sT1 6 year - 6 year survival \cr
      Kaplan-Meier survival curve of the group above
  \item gT2sT2 - number of samples in \cr
  \item gT2sT2 3 year - 3 year survival \cr
      Kaplan-Meier survival curve of the group above
  \item gT2sT2 6 year - 6 year survival \cr
      Kaplan-Meier survival curve of the group above
  \item gT2sT3 - number of samples in \cr
  \item gT2sT3 3 year - 3 year survival \cr
      Kaplan-Meier survival curve of the group above
  \item gT2sT3 6 year - 6 year survival \cr
      Kaplan-Meier survival curve of the group above
  \item gT3 Pval - Pvalue of gene Tertile 3
  \item gT3 HR - Hazard ratio of gene Tertile 3
  \item gT3sT1 - number of samples in \cr
  \item gT3sT1 3 year - 3 year survival \cr
      Kaplan-Meier survival curve of the group above
  \item gT3sT1 6 year - 6 year survival \cr
      Kaplan-Meier survival curve of the group above
  \item gT3sT2 - number of samples in \cr
  \item gT3sT2 3 year - 3 year survival \cr
      Kaplan-Meier survival curve of the group above
  \item gT3sT2 6 year - 6 year survival \cr
      Kaplan-Meier survival curve of the group above
  \item gT3sT3 - number of samples in \cr
  \item gT3sT3 3 year - 3 year survival \cr
      Kaplan-Meier survival curve of the group above
  \item gT3sT3 6 year - 6 year survival \cr
      Kaplan-Meier survival curve of the group above
  \item warnings - if the cox regress fails, the warning should be recorded here
}
}
\description{
Cox regression after separating prognistic variable by gene expression tertiles
}
