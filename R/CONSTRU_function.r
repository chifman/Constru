#!/usr/bin/env Rscript


can.be.numeric <- function(x) {
    stopifnot(is.atomic(x) || is.list(x)) # check if x is a vector
    numNAs = sum(is.na(x))
    numNAs_new = suppressWarnings(sum(is.na(as.numeric(x))))
    return(numNAs_new == numNAs)
}

check_numeric=function(myDF){
	myDF2 = as.data.frame(lapply(as.data.frame(myDF), function(col) {
		if (can.be.numeric(col)) { as.numeric(col) } else { col }
	}))
	colnames(myDF2)=colnames(myDF)
	rownames(myDF2)=rownames(myDF)
	return(as.data.frame(myDF2))
}

categorize_tertiles=function(x){
	ox=order(x)
	y=rep(3,length(x))
	y[ox[1:ceiling(length(x)*1/3)]]=1
	y[ox[(ceiling(length(x)*1/3)+1):(ceiling(length(x)*2/3))]]=2
	return(y)
}

constru_single_continuous=function(gene, survival_data, gene_data, prognostic_variable_d,cox_formula){
	suppressMessages(require(survival))
	suppressMessages(require(survminer))
	# output format
	column_names<-c("mean","P5%","P95%","Diff P95%_5%","R",
	"prognostic_variable_coef",
	"prognostic_variable_HR",
	"prognostic_variable_Pval",
	"gene_data_coef",
	"gene_data_HR",
	"gene_data_Pval",
	"prognostic_variable:gene_data_coef",
	"prognostic_variable:gene_data_HR",
	"prognostic_variable:gene_data_Pval",
	"warnings"
	)
	output=rep(NA,length(column_names))
	names(output)=column_names
	tryCatch({
		gene_single=unlist(gene_data[gene,])
		data_subset=cbind(survival_data,prognostic_variable=as.vector(prognostic_variable_d),gene_data=gene_single)
		# 
		output['mean']=mean(gene_single)
		output['P5%']=quantile(gene_single, c(0.05))
		output['P95%']=quantile(gene_single, c(0.95))
		output["Diff P95%_5%"]=quantile(gene_single, c(0.95)) - quantile(gene_single, c(0.05))
		output["R"]=cor(gene_single, prognostic_variable_d, method="pearson")
		gene_tertiles=categorize_tertiles(gene_single)
		prognostic_variable_d_tertiles=categorize_tertiles(prognostic_variable_d)
		cox_formula2=gsub("prognostic_variable","prognostic_variable+gene_data+prognostic_variable:gene_data",cox_formula)
		res.cox = survival::coxph(as.formula(cox_formula2), data = data_subset)
		sres.cox= summary(res.cox)
		output["prognostic_variable_coef"]=sres.cox$coef["prognostic_variable","coef"]
		output["prognostic_variable_HR"]=sres.cox$coef["prognostic_variable","exp(coef)"]
		output["prognostic_variable_Pval"]=sres.cox$coef["prognostic_variable","Pr(>|z|)"]
		output["gene_data_coef"]=sres.cox$coef["gene_data","coef"]
		output["gene_data_HR"]=sres.cox$coef["gene_data","exp(coef)"]
		output["gene_data_Pval"]=sres.cox$coef["gene_data","Pr(>|z|)"]
		output["prognostic_variable:gene_data_coef"]=sres.cox$coef["prognostic_variable:gene_data","coef"]
		output["prognostic_variable:gene_data_HR"]=sres.cox$coef["prognostic_variable:gene_data","exp(coef)"]
		output["prognostic_variable:gene_data_Pval"]=sres.cox$coef["prognostic_variable:gene_data","Pr(>|z|)"]
#		res.cox = survival::survfit(as.formula(cox_formula2), data = data_subset)
	},error = function(cond){output["warnings"]=paste(output["warnings"],conditionMessage(cond))}) 
	return(output)
}

constru_single_Mclust=function(gene, survival_data, gene_data, prognostic_variable_d,cox_formula){
	suppressMessages(require(survival))
	suppressMessages(require(survminer))
	suppressMessages(require(mclust))
	# output format
	column_names<-c("mean","P5%","P95%","Diff P95%_5%","R",
	"nGroups",
	"PS1","PS1_max_group","PS1_min_group",
	"PS2","PS2_max_group","PS2_min_group",
	"PS3","PS3_max_group","PS3_min_group",
	"All_Pval","All_HR",
	gsub("^","Group_",colnames(gene_data)),
	"warnings")
	output=rep(NA,length(column_names))
	names(output)=column_names
	gene_single=unlist(gene_data[gene,])
	data_subset=cbind(survival_data,prognostic_variable=as.numeric(prognostic_variable_d),gene_data=gene_single)
	#
	gene_tertiles=categorize_tertiles(gene_single)
	prognostic_variable_d_tertiles=categorize_tertiles(prognostic_variable_d)
	mydata=data_subset[,c("prognostic_variable","gene_data")]
	# 
	output['mean']=mean(gene_single)
	output['P5%']=quantile(gene_single, c(0.05))
	output['P95%']=quantile(gene_single, c(0.95))
	output["Diff P95%_5%"]=quantile(gene_single, c(0.95)) - quantile(gene_single, c(0.05))
	output["R"]=cor(gene_single, prognostic_variable_d, method="pearson")
	#Model-based clustering based on parameterized finite Gaussian mixture models
	fit <- mclust::Mclust(mydata,verbose=F)
	output[gsub("^","Group_",colnames(gene_data))]=fit$classification
	nGroups=max(fit$classification)
	output["nGroups"]=nGroups
	#
	All_Pval=c()
	All_HR=c()
	for(group in 1:nGroups){
	tryCatch({
		keep= (output[gsub("^","Group_",rownames(data_subset))]==group)
		res.cox = survival::coxph(as.formula(cox_formula), data = data_subset[keep,])
		sres.cox= summary(res.cox)
		All_Pval=c(All_Pval,sres.cox$coef["prognostic_variable","Pr(>|z|)"])
		All_HR=c(All_HR,sres.cox$coef["prognostic_variable","exp(coef)"])
	},error = function(cond){output["warnings"]=paste(output["warnings"],conditionMessage(cond))}) 
	}
	output["All_Pval"]=paste(All_Pval,sep=",",collapse=",")
	output["All_HR"]=paste(All_HR,sep=",",collapse=",")
	tryCatch({
	PS1_part=sapply(1:length(All_Pval),function(g){
		(-log10(All_Pval[g])/All_HR[g])
	})
	PS1_max_group=which(PS1_part==max(PS1_part,na.rm=T))[1]
	PS1_min_group=which(PS1_part==min(PS1_part,na.rm=T))[1]
	output['PS1_max_group']=PS1_max_group
	output['PS1_min_group']=PS1_min_group
	output['PS1']=
		(-log10(All_Pval[PS1_max_group])/All_HR[PS1_max_group])-
		(-log10(All_Pval[PS1_min_group])/All_HR[PS1_min_group])
	PS2_part=sapply(1:length(All_Pval),function(g){
		(log(All_Pval[g])*(All_HR[g]-1))
	})
	PS2_max_group=which(PS2_part==max(PS2_part,na.rm=T))[1]
	PS2_min_group=which(PS2_part==min(PS2_part,na.rm=T))[1]
	output['PS2_max_group']=PS2_max_group
	output['PS2_min_group']=PS2_min_group
	output['PS2']=
		(log(All_Pval[PS2_max_group])*(All_HR[PS2_max_group]-1)) -
		(log(All_Pval[PS2_min_group])*(All_HR[PS2_min_group]-1))
	PS3_part=sapply(1:length(All_Pval),function(g){
		(log(All_Pval[g])*log(All_HR[g]))
	})
	PS3_max_group=which(PS3_part==max(PS3_part,na.rm=T))[1]
	PS3_min_group=which(PS3_part==min(PS3_part,na.rm=T))[1]
	output['PS3_max_group']=PS3_max_group
	output['PS3_min_group']=PS3_min_group
	output['PS3']=
		(log(All_Pval[PS3_max_group])*log(All_HR[PS3_max_group])) -
		(log(All_Pval[PS3_min_group])*log(All_HR[PS3_min_group]))
	},error = function(cond){output["warnings"]=paste(output["warnings"],conditionMessage(cond))}) 
	return(output)
}

constru_single=function(gene, survival_data, gene_data, prognostic_variable_d,cox_formula){
	suppressMessages(require(survival))
	suppressMessages(require(survminer))
	# output format
	column_names<-c("mean","P5%","P95%","Diff P95%_5%","R","PS1","PS2","PS3",
	"gT1 Pval","gT1 HR",
	"gT1sT1","gT1sT1 3 year","gT1sT1 6 year",
	"gT1sT2","gT1sT2 3 year","gT1sT2 6 year",
	"gT1sT3","gT1sT3 3 year","gT1sT3 6 year",
	"gT2 Pval","gT2 HR",
	"gT2sT1","gT2sT1 3 year","gT2sT1 6 year",
	"gT2sT2","gT2sT2 3 year","gT2sT2 6 year",
	"gT2sT3","gT2sT3 3 year","gT2sT3 6 year",
	"gT3 Pval","gT3 HR",
	"gT3sT1","gT3sT1 3 year","gT3sT1 6 year",
	"gT3sT2","gT3sT2 3 year","gT3sT2 6 year",
	"gT3sT3","gT3sT3 3 year","gT3sT3 6 year","warnings")
	output=rep(NA,length(column_names))
	names(output)=column_names
	# subset
	tryCatch({
		hold=unlist(gene_data[gene,])
		data_subset=cbind(survival_data,prognostic_variable=as.vector(prognostic_variable_d))
		# 
		output['mean']=mean(hold)
		output['P5%']=quantile(hold, c(0.05))
		output['P95%']=quantile(hold, c(0.95))
		output["Diff P95%_5%"]=quantile(hold, c(0.95)) - quantile(hold, c(0.05))
		output["R"]=cor(hold, prognostic_variable_d, method="pearson")
		#
		gene_tertiles=categorize_tertiles(hold)
		prognostic_variable_d_tertiles=categorize_tertiles(prognostic_variable_d)
		# 
		gt_cutoff=1
		keep=(gene_tertiles == gt_cutoff)
		res.cox = survival::coxph(as.formula(cox_formula), data = data_subset[keep,])
		sres.cox= summary(res.cox)
		output["gT1 Pval"]=sres.cox$coef["prognostic_variable","Pr(>|z|)"]
		output["gT1 HR"]=sres.cox$coef["prognostic_variable","exp(coef)"]
		# 
		gt_cutoff=2
		keep=(gene_tertiles == gt_cutoff)
		res.cox = survival::coxph(as.formula(cox_formula), data = data_subset[keep,])
		sres.cox= summary(res.cox)
		output["gT2 Pval"]=sres.cox$coef["prognostic_variable","Pr(>|z|)"]
		output["gT2 HR"]=sres.cox$coef["prognostic_variable","exp(coef)"]
		# 
		gt_cutoff=3
		keep=(gene_tertiles == gt_cutoff)
		res.cox = survival::coxph(as.formula(cox_formula), data = data_subset[keep,])
		sres.cox= summary(res.cox)
		output["gT3 Pval"]=sres.cox$coef["prognostic_variable","Pr(>|z|)"]
		output["gT3 HR"]=sres.cox$coef["prognostic_variable","exp(coef)"]
		#
		output["PS1"]=
			(-log10(output["gT1 Pval"])/output["gT1 HR"]) - 
			(-log10(output["gT3 Pval"])/output["gT3 HR"])
		output["PS2"]=
			(log(output["gT1 Pval"])*(output["gT1 HR"]-1)) - 
			(log(output["gT3 Pval"])*(output["gT3 HR"]-1))
		output["PS3"]= 
			(log(output["gT1 Pval"])*log(output["gT1 HR"]))- 
			(log(output["gT3 Pval"])*log(output["gT3 HR"]))
		#
	},error = function(cond){output["warnings"]=paste(output["warnings"],conditionMessage(cond))}) 
	tryCatch({
		data_subset_tertile=cbind(survival_data,prognostic_variable=as.vector(prognostic_variable_d_tertiles))
		#
		gt_cutoff=1
		keep=(gene_tertiles == gt_cutoff)
		res.cox = survival::survfit(as.formula(cox_formula), data = data_subset_tertile[keep,])
		output["gT1sT1"]=res.cox['prognostic_variable=1']$n
		output["gT1sT1 3 year"]=summary(res.cox['prognostic_variable=1'],time=c(3),extend=T)$surv
		output["gT1sT1 6 year"]=summary(res.cox['prognostic_variable=1'],time=c(6),extend=T)$surv
	},error = function(cond){output["warnings"]=paste(output["warnings"],conditionMessage(cond))}) 
	tryCatch({
		output["gT1sT2"]=res.cox['prognostic_variable=2']$n
		output["gT1sT2 3 year"]=summary(res.cox['prognostic_variable=2'],time=c(3),extend=T)$surv
		output["gT1sT2 6 year"]=summary(res.cox['prognostic_variable=2'],time=c(6),extend=T)$surv
	},error = function(cond){output["warnings"]=paste(output["warnings"],conditionMessage(cond))}) 
	tryCatch({
		output["gT1sT3"]=res.cox['prognostic_variable=3']$n
		output["gT1sT3 3 year"]=summary(res.cox['prognostic_variable=3'],time=c(3),extend=T)$surv
		output["gT1sT3 6 year"]=summary(res.cox['prognostic_variable=3'],time=c(6),extend=T)$surv
		#
	},error = function(cond){output["warnings"]=paste(output["warnings"],conditionMessage(cond))}) 
	tryCatch({
		gt_cutoff=2
		keep=(gene_tertiles == gt_cutoff)
		res.cox = survival::survfit(as.formula(cox_formula), data = data_subset_tertile[keep,])
	},error = function(cond){output["warnings"]=paste(output["warnings"],conditionMessage(cond))}) 
	tryCatch({
		output["gT2sT1"]=res.cox['prognostic_variable=1']$n
		output["gT2sT1 3 year"]=summary(res.cox['prognostic_variable=1'],time=c(3),extend=T)$surv
		output["gT2sT1 6 year"]=summary(res.cox['prognostic_variable=1'],time=c(6),extend=T)$surv
	},error = function(cond){output["warnings"]=paste(output["warnings"],conditionMessage(cond))}) 
	tryCatch({
		output["gT2sT2"]=res.cox['prognostic_variable=2']$n
		output["gT2sT2 3 year"]=summary(res.cox['prognostic_variable=2'],time=c(3),extend=T)$surv
		output["gT2sT2 6 year"]=summary(res.cox['prognostic_variable=2'],time=c(6),extend=T)$surv
	},error = function(cond){output["warnings"]=paste(output["warnings"],conditionMessage(cond))}) 
	tryCatch({
		output["gT2sT3"]=res.cox['prognostic_variable=3']$n
		output["gT2sT3 3 year"]=summary(res.cox['prognostic_variable=3'],time=c(3),extend=T)$surv
		output["gT2sT3 6 year"]=summary(res.cox['prognostic_variable=3'],time=c(6),extend=T)$surv
		#
	},error = function(cond){output["warnings"]=paste(output["warnings"],conditionMessage(cond))}) 
	tryCatch({
		gt_cutoff=3
		keep=(gene_tertiles == gt_cutoff)
		res.cox = survival::survfit(as.formula(cox_formula), data = data_subset_tertile[keep,])
	},error = function(cond){output["warnings"]=paste(output["warnings"],conditionMessage(cond))}) 
	tryCatch({
		output["gT3sT1"]=res.cox['prognostic_variable=1']$n
		output["gT3sT1 3 year"]=summary(res.cox['prognostic_variable=1'],time=c(3),extend=T)$surv
		output["gT3sT1 6 year"]=summary(res.cox['prognostic_variable=1'],time=c(6),extend=T)$surv
	},error = function(cond){output["warnings"]=paste(output["warnings"],conditionMessage(cond))}) 
	tryCatch({
		output["gT3sT2"]=res.cox['prognostic_variable=2']$n
		output["gT3sT2 3 year"]=summary(res.cox['prognostic_variable=2'],time=c(3),extend=T)$surv
		output["gT3sT2 6 year"]=summary(res.cox['prognostic_variable=2'],time=c(6),extend=T)$surv
	},error = function(cond){output["warnings"]=paste(output["warnings"],conditionMessage(cond))}) 
	tryCatch({
		output["gT3sT3"]=res.cox['prognostic_variable=3']$n
		output["gT3sT3 3 year"]=summary(res.cox['prognostic_variable=3'],time=c(3),extend=T)$surv
		output["gT3sT3 6 year"]=summary(res.cox['prognostic_variable=3'],time=c(6),extend=T)$surv
	},error = function(cond){output["warnings"]=paste(output["warnings"],conditionMessage(cond))}) 
	return(output)
}

#constru_single_mode=function(gene, survival_data, gene_data, prognostic_variable_d,cox_formula,mode){
#	if(mode=="tertile"){
#		return(constru_single(gene, survival_data, gene_data, prognostic_variable_d,cox_formula))
#	}
#	if(mode=="continuous"){
#		return(constru_single_continuous(gene, survival_data, gene_data, prognostic_variable_d,cox_formula))
#	}
#	if(mode=="Mclust"){
#		return(constru_single_Mclust(gene, survival_data, gene_data, prognostic_variable_d,cox_formula))
#	}
#}

survival_data_u=NULL
gene_data_u=NULL
prognostic_variable_d_u=NULL
cox_formula_u=NULL

#' constru
#'
#' Cox regression after separating prognistic variable by gene expression tertiles
#' @param survival_data A data frame with the samples as rows and survival time and event as columns.
#' @param gene_data A data frame sample names as columns and gene expression as rows.
#' @param prognostic_variable_d A vector with the prognostic variable data. It should be of the same length and order as the row of the survival data.
#' @param cox_formula Formula used for cox regression. Example: 
#' \itemize{
#'   \item Surv( OS_TIME , OS_Event ) ~ prognostic_variable \cr
#' }
#'   If additional prognostic factors need to be added to the model, add it as a column in the survival data and use its column name in the formula. Example:
#' \itemize{
#'   \item Surv( OS_TIME , OS_Event ) ~ prognostic_variable + AGE \cr
#' }
#' @param ncores The number of cores used during multithreading.
#' @return a table with columns:
#' \itemize{
#'   \item mean - gene expression mean
#'   \item P5\% - 5\% confidence interval of the gene expression
#'   \item P95\% - 95\% confidence interval of the gene expression
#'   \item Diff P95\%_5\% - difference between the 95\% and 5\% confidence interval
#'   \item R - pearson correlation between the gene expression and prognostic variable
#'   \item PS1 - parity score 1 \cr
#'      (-log10("Pvalue of gene Tertile 1")/"HR of gene Tertile 1") - \cr
#'      (-log10("Pvalue of gene Tertile 3")/"HR of gene Tertile 3")   \cr
#'   \item PS2 - parity score 2
#'      (log("Pvalue of gene Tertile 1")*("HR of gene Tertile 1"-1)) - \cr
#'      (log("Pvalue of gene Tertile 3")*("HR of gene Tertile 3"-1))   \cr
#'   \item PS3 - parity score 3
#'      (log("Pvalue of gene Tertile 1")*log("HR of gene Tertile 1"))- \cr
#'      (log("Pvalue of gene Tertile 3")*log("HR of gene Tertile 3"))  \cr
#'   \item gT1 Pval - Pvalue of gene Tertile 1
#'   \item gT1 HR - HR of gene Tertile 1
#'   \item gT1sT1 - 
#'   \item gT1sT1 3 year - 
#'   \item gT1sT1 6 year - 
#'   \item gT1sT2 - 
#'   \item gT1sT2 3 year - 
#'   \item gT1sT2 6 year - 
#'   \item gT1sT3 - 
#'   \item gT1sT3 3 year - 
#'   \item gT1sT3 6 year - 
#'   \item gT2 Pval - 
#'   \item gT2 HR - 
#'   \item gT2sT1 - 
#'   \item gT2sT1 3 year - 
#'   \item gT2sT1 6 year - 
#'   \item gT2sT2 - 
#'   \item gT2sT2  - 
#'   \item 3 year - 
#'   \item gT2sT2 6 year - 
#'   \item gT2sT3 - 
#'   \item gT2sT3 3 year - 
#'   \item gT2sT3 6 year - 
#'   \item gT3 Pval - 
#'   \item gT3 HR - 
#'   \item gT3sT1 - 
#'   \item gT3sT1 3 year - 
#'   \item gT3sT1 6 year - 
#'   \item gT3sT2 - 
#'   \item gT3sT2 3 year - 
#'   \item gT3sT2 6 year - 
#'   \item gT3sT3 - 
#'   \item gT3sT3 3 year - 
#'   \item gT3sT3 6 year - 
#'   \item warnings - 
#' }
#' @examples 
#' todo
#' @export
constru<-function(survival_data, gene_data, prognostic_variable_d,cox_formula,ncores){
	suppressMessages(require(parallel))
	gi=rownames(gene_data)
	oo=NULL
	if( Sys.info()[['sysname']] == 'Windows' ){
		survival_data_u=survival_data
		gene_data_u=gene_data
		prognostic_variable_d_u=prognostic_variable_d
		cox_formula_u=cox_formula
		cl <- parallel::makeCluster(ncores)
		parallel::clusterExport(cl, varlist=c("categorize_tertiles","constru_single","survival_data_u","gene_data_u","prognostic_variable_d_u","cox_formula_u"))
		parallel::clusterEvalQ(cl, { library(survival); library(survminer); library(parallel); library(mclust); })
		oo=parallel::parLapply(cl,gi,function(x){ constru_single(x,survival_data_u,gene_data_u,prognostic_variable_d_u,cox_formula_u) })
		parallel::stopCluster(cl)
	} else {
		oo=parallel::mclapply(gi,function(x){constru_single(x,survival_data,gene_data,prognostic_variable_d,cox_formula)} ,mc.cores=ncores)
	}
	oo=t(as.data.frame(oo))
	rownames(oo)=gi
	return(oo)
}

