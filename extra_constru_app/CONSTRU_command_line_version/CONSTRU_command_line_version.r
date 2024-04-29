#!/usr/bin/env Rscript
 
suppressMessages(library(optparse))
suppressMessages(library(parallel))
suppressMessages(library(survival))
suppressMessages(library(survminer))

# parse arguments
option_list = list(
make_option(c("-m", "--metadata"), type="character", default=NULL, 
	help="Metadata file name. This file should have a csv table with sample names as column and metadata as rows.", metavar="character"),
make_option(c("-d", "--data"), type="character", default=NULL, 
	help="Gene expression data file name. This file should have a csv table with sample names as column and gene expression as rows.", metavar="character"),
make_option(c("-r", "--metadata_row"), type="character", default=NULL, 
	help="Specify which row of the metadata should be used as the metagene values", metavar="character"),
make_option(c("-f", "--formula"), type="character", default=NULL, 
	help=" Formula used for cox regression. Example: 
		If the metadata file stores 
			the survival time in row 'OS_TIME ',
			the survival event in row 'OS_Event', and
			the covariate age in row 'AGE',
		the formula for cox regression formula is:
		'Surv( OS_TIME , OS_Event ) ~ metagene + AGE'", metavar="character"),
make_option(c("-t", "--threads"), type="numeric", default=detectCores(), 
	help="Number of threads [default= %default]", metavar="numeric"),
make_option(c("-o", "--out"), type="character", default="out.txt", 
	help="Output file name [default= %default]", metavar="character")
) 

description="
(Note: all non-alphanumeric characters except for underscore in the row and column names will be converted to a dot)
"

opt_parser = OptionParser(option_list=option_list,description=description)
opt = parse_args(opt_parser)

if (is.null(opt$metadata)){
	print_help(opt_parser)
	stop("At least one argument must be supplied for input metadata\n")
}
metadata_file=opt$metadata

if (is.null(opt$data)){
	print_help(opt_parser)
	stop("At least one argument must be supplied for input data\n")
}
data_file=opt$data

if (is.null(opt$metadata_row)){
	print_help(opt_parser)
	stop("At least one argument must be supplied for input metadata_row\n")
}
metagenes_value_name=opt$metadata_row

if (is.null(opt$formula)){
	print_help(opt_parser)
	stop("At least one argument must be supplied for input formula\n")
}
custom_formula=opt$formula

number_of_threads=opt$threads
output_filename=opt$out

# example inputs
#metadata_file="Test_COMBAT_Combined_forMing3_part1.csv"
#data_file="Test_COMBAT_Combined_forMing3_part2.csv"
#metagene_mean_name='CYT.Tertiles'
#custom_formula='Surv( OS.TIME..yrs. , OS.Event..0.censored..1.death.) ~ metagene'
#number_of_threads=1
#output_filename="Test_COMBAT_Combined_forMing3_CONSTRU.csv"

#####

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

######

# read in data
metadata=read.csv(metadata_file)
rownames(metadata)=metadata[,1]
metadata=check_numeric(t(metadata[,-1]))
colnames(metadata)=make.names(colnames(metadata))

data=read.csv(data_file)
rownames(data)=data[,1]
data=data[,-1]

# make sure that the dimensions match

# make sure that the order of data and metadata is the same.
metadata=metadata[colnames(data),]

# make sure the metagene column exists
# check column names
metagenes_value=metadata[,metagenes_value_name]

# run
output=constru(metadata, data, metagenes_value,custom_formula,number_of_threads)
write.csv(output,file=output_filename)

