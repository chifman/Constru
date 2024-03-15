#!/usr/bin/env Rscript
 
suppressMessages(library(parallel))
suppressMessages(library(survival))
suppressMessages(library(survminer))

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

#categorize_tertiles=function(x){
#	y=rep(3,length(x))
#	y[x <= quantile(x,2/3)]=2
#	y[x <= quantile(x,1/3)]=1
#	return(y)
#}

categorize_tertiles=function(x){
	ox=order(x)
	y=rep(3,length(x))
	y[ox[1:ceiling(length(x)*1/3)]]=1
	y[ox[(ceiling(length(x)*1/3)+1):(ceiling(length(x)*2/3))]]=2
	return(y)
}

constru_single=function(gene, clinical, gene_data, metagene_mean,cox_formula){
	require(survival)
	require(survminer)
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
		data_subset=cbind(clinical,metagene=as.vector(metagene_mean))
		# 
		output['mean']=mean(hold)
		output['P5%']=quantile(hold, c(0.05))
		output['P95%']=quantile(hold, c(0.95))
		output["Diff P95%_5%"]=quantile(hold, c(0.95)) - quantile(hold, c(0.05))
		output["R"]=cor(hold, metagene_mean, method="pearson")
		#
		gene_tertiles=categorize_tertiles(hold)
		metagene_mean_tertiles=categorize_tertiles(metagene_mean)
		# 
		gt_cutoff=1
		keep=(gene_tertiles == gt_cutoff)
		res.cox = coxph(as.formula(cox_formula), data = data_subset[keep,])
		sres.cox= summary(res.cox)
		output["gT1 Pval"]=sres.cox$coef["metagene","Pr(>|z|)"]
		output["gT1 HR"]=sres.cox$coef["metagene","exp(coef)"]
		# 
		gt_cutoff=2
		keep=(gene_tertiles == gt_cutoff)
		res.cox = coxph(as.formula(cox_formula), data = data_subset[keep,])
		sres.cox= summary(res.cox)
		output["gT2 Pval"]=sres.cox$coef["metagene","Pr(>|z|)"]
		output["gT2 HR"]=sres.cox$coef["metagene","exp(coef)"]
		# 
		gt_cutoff=3
		keep=(gene_tertiles == gt_cutoff)
		res.cox = coxph(as.formula(cox_formula), data = data_subset[keep,])
		sres.cox= summary(res.cox)
		output["gT3 Pval"]=sres.cox$coef["metagene","Pr(>|z|)"]
		output["gT3 HR"]=sres.cox$coef["metagene","exp(coef)"]
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
		data_subset_tertile=cbind(clinical,metagene=as.vector(metagene_mean_tertiles))
		#
		gt_cutoff=1
		keep=(gene_tertiles == gt_cutoff)
		res.cox = survfit(as.formula(cox_formula), data = data_subset_tertile[keep,])
		output["gT1sT1"]=res.cox['metagene=1']$n
		output["gT1sT1 3 year"]=summary(res.cox['metagene=1'],time=c(3),extend=T)$surv
		output["gT1sT1 6 year"]=summary(res.cox['metagene=1'],time=c(6),extend=T)$surv
	},error = function(cond){output["warnings"]=paste(output["warnings"],conditionMessage(cond))}) 
	tryCatch({
		output["gT1sT2"]=res.cox['metagene=2']$n
		output["gT1sT2 3 year"]=summary(res.cox['metagene=2'],time=c(3),extend=T)$surv
		output["gT1sT2 6 year"]=summary(res.cox['metagene=2'],time=c(6),extend=T)$surv
	},error = function(cond){output["warnings"]=paste(output["warnings"],conditionMessage(cond))}) 
	tryCatch({
		output["gT1sT3"]=res.cox['metagene=3']$n
		output["gT1sT3 3 year"]=summary(res.cox['metagene=3'],time=c(3),extend=T)$surv
		output["gT1sT3 6 year"]=summary(res.cox['metagene=3'],time=c(6),extend=T)$surv
		#
	},error = function(cond){output["warnings"]=paste(output["warnings"],conditionMessage(cond))}) 
	tryCatch({
		gt_cutoff=2
		keep=(gene_tertiles == gt_cutoff)
		res.cox = survfit(as.formula(cox_formula), data = data_subset_tertile[keep,])
	},error = function(cond){output["warnings"]=paste(output["warnings"],conditionMessage(cond))}) 
	tryCatch({
		output["gT2sT1"]=res.cox['metagene=1']$n
		output["gT2sT1 3 year"]=summary(res.cox['metagene=1'],time=c(3),extend=T)$surv
		output["gT2sT1 6 year"]=summary(res.cox['metagene=1'],time=c(6),extend=T)$surv
	},error = function(cond){output["warnings"]=paste(output["warnings"],conditionMessage(cond))}) 
	tryCatch({
		output["gT2sT2"]=res.cox['metagene=2']$n
		output["gT2sT2 3 year"]=summary(res.cox['metagene=2'],time=c(3),extend=T)$surv
		output["gT2sT2 6 year"]=summary(res.cox['metagene=2'],time=c(6),extend=T)$surv
	},error = function(cond){output["warnings"]=paste(output["warnings"],conditionMessage(cond))}) 
	tryCatch({
		output["gT2sT3"]=res.cox['metagene=3']$n
		output["gT2sT3 3 year"]=summary(res.cox['metagene=3'],time=c(3),extend=T)$surv
		output["gT2sT3 6 year"]=summary(res.cox['metagene=3'],time=c(6),extend=T)$surv
		#
	},error = function(cond){output["warnings"]=paste(output["warnings"],conditionMessage(cond))}) 
	tryCatch({
		gt_cutoff=3
		keep=(gene_tertiles == gt_cutoff)
		res.cox = survfit(as.formula(cox_formula), data = data_subset_tertile[keep,])
	},error = function(cond){output["warnings"]=paste(output["warnings"],conditionMessage(cond))}) 
	tryCatch({
		output["gT3sT1"]=res.cox['metagene=1']$n
		output["gT3sT1 3 year"]=summary(res.cox['metagene=1'],time=c(3),extend=T)$surv
		output["gT3sT1 6 year"]=summary(res.cox['metagene=1'],time=c(6),extend=T)$surv
	},error = function(cond){output["warnings"]=paste(output["warnings"],conditionMessage(cond))}) 
	tryCatch({
		output["gT3sT2"]=res.cox['metagene=2']$n
		output["gT3sT2 3 year"]=summary(res.cox['metagene=2'],time=c(3),extend=T)$surv
		output["gT3sT2 6 year"]=summary(res.cox['metagene=2'],time=c(6),extend=T)$surv
	},error = function(cond){output["warnings"]=paste(output["warnings"],conditionMessage(cond))}) 
	tryCatch({
		output["gT3sT3"]=res.cox['metagene=3']$n
		output["gT3sT3 3 year"]=summary(res.cox['metagene=3'],time=c(3),extend=T)$surv
		output["gT3sT3 6 year"]=summary(res.cox['metagene=3'],time=c(6),extend=T)$surv
	},error = function(cond){output["warnings"]=paste(output["warnings"],conditionMessage(cond))}) 
	return(output)
}

clinical_u=NULL
gene_data_u=NULL
metagene_mean_u=NULL
cox_formula_u=NULL

#' constru
#'
#' Cox regression after separating metagene_mean's tertile category
#' @param clinical, gene_data, metagene_mean,cox_formula,ncores
#' @return a table detailing ...
#' @examples 
#' todo
#' @export
constru<-function(clinical, gene_data, metagene_mean,cox_formula,ncores){
	# remove rows where more than 1/3 the sample have zero expression
	keep=apply(gene_data,1,function(x){sum(x==0)/length(x) < 1/3})
	gene_data=gene_data[keep,]
	# 
	gi=rownames(gene_data)
	oo=NULL
	if( Sys.info()[['sysname']] == 'Windows' ){
		clinical_u=clinical
		gene_data_u=gene_data
		metagene_mean_u=metagene_mean
		cox_formula_u=cox_formula
		cl <- makeCluster(ncores)
		clusterExport(cl, varlist=c("categorize_tertiles","constru_single","clinical_u","gene_data_u","metagene_mean_u","cox_formula_u"))
		clusterEvalQ(cl, { library(survival); library(survminer); library(parallel); })
		oo=parLapply(cl,gi,function(x){ constru_single(x,clinical_u,gene_data_u,metagene_mean_u,cox_formula_u) })
		stopCluster(cl)
	} else {
		oo=mclapply(gi,function(x){constru_single(x,clinical,gene_data,metagene_mean,cox_formula)} ,mc.cores=ncores)
	}
	oo=t(as.data.frame(oo))
	rownames(oo)=gi
	return(oo)
}
