#!/usr/bin/env Rscript

library(parallel)
library(survival)
library(survminer)

load("test.Rdata")

categorize_tertiles=function(x){
	y=rep(3,length(x))
	y[x <= quantile(x,2/3)]=2
	y[x <= quantile(x,1/3)]=1
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
	"gT3sT3","gT3sT3 3 year","gT3sT3 6 year")
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
	},error = function(cond){}) 
	tryCatch({
		data_subset_tertile=cbind(clinical,metagene=as.vector(metagene_mean_tertiles))
		#
		gt_cutoff=1
		keep=(gene_tertiles == gt_cutoff)
		res.cox = survfit(as.formula(cox_formula), data = data_subset_tertile[keep,])
		output["gT1sT1"]=res.cox['metagene=1']$n
		output["gT1sT1 3 year"]=summary(res.cox['metagene=1'],time=c(3),extend=T)$surv
		output["gT1sT1 6 year"]=summary(res.cox['metagene=1'],time=c(6),extend=T)$surv
	},error = function(cond){}) 
	tryCatch({
		output["gT1sT2"]=res.cox['metagene=2']$n
		output["gT1sT2 3 year"]=summary(res.cox['metagene=2'],time=c(3),extend=T)$surv
		output["gT1sT2 6 year"]=summary(res.cox['metagene=2'],time=c(6),extend=T)$surv
	},error = function(cond){}) 
	tryCatch({
		output["gT1sT3"]=res.cox['metagene=3']$n
		output["gT1sT3 3 year"]=summary(res.cox['metagene=3'],time=c(3),extend=T)$surv
		output["gT1sT3 6 year"]=summary(res.cox['metagene=3'],time=c(6),extend=T)$surv
		#
	},error = function(cond){}) 
	tryCatch({
		gt_cutoff=2
		keep=(gene_tertiles == gt_cutoff)
		res.cox = survfit(as.formula(cox_formula), data = data_subset_tertile[keep,])
	},error = function(cond){}) 
	tryCatch({
		output["gT2sT1"]=res.cox['metagene=1']$n
		output["gT2sT1 3 year"]=summary(res.cox['metagene=1'],time=c(3),extend=T)$surv
		output["gT2sT1 6 year"]=summary(res.cox['metagene=1'],time=c(6),extend=T)$surv
	},error = function(cond){}) 
	tryCatch({
		output["gT2sT2"]=res.cox['metagene=2']$n
		output["gT2sT2 3 year"]=summary(res.cox['metagene=2'],time=c(3),extend=T)$surv
		output["gT2sT2 6 year"]=summary(res.cox['metagene=2'],time=c(6),extend=T)$surv
	},error = function(cond){}) 
	tryCatch({
		output["gT2sT3"]=res.cox['metagene=3']$n
		output["gT2sT3 3 year"]=summary(res.cox['metagene=3'],time=c(3),extend=T)$surv
		output["gT2sT3 6 year"]=summary(res.cox['metagene=3'],time=c(6),extend=T)$surv
		#
	},error = function(cond){}) 
	tryCatch({
		gt_cutoff=3
		keep=(gene_tertiles == gt_cutoff)
		res.cox = survfit(as.formula(cox_formula), data = data_subset_tertile[keep,])
	},error = function(cond){}) 
	tryCatch({
		output["gT3sT1"]=res.cox['metagene=1']$n
		output["gT3sT1 3 year"]=summary(res.cox['metagene=1'],time=c(3),extend=T)$surv
		output["gT3sT1 6 year"]=summary(res.cox['metagene=1'],time=c(6),extend=T)$surv
	},error = function(cond){}) 
	tryCatch({
		output["gT3sT2"]=res.cox['metagene=2']$n
		output["gT3sT2 3 year"]=summary(res.cox['metagene=2'],time=c(3),extend=T)$surv
		output["gT3sT2 6 year"]=summary(res.cox['metagene=2'],time=c(6),extend=T)$surv
	},error = function(cond){}) 
	tryCatch({
		output["gT3sT3"]=res.cox['metagene=3']$n
		output["gT3sT3 3 year"]=summary(res.cox['metagene=3'],time=c(3),extend=T)$surv
		output["gT3sT3 6 year"]=summary(res.cox['metagene=3'],time=c(6),extend=T)$surv
	},error = function(cond){}) 
	return(output)
}

constru<-function(clinical, gene_data, metagene_mean,cox_formula,ncores){
	gi=rownames(gene_data)
	oo=mclapply(gi,function(x){constru_single(x,clinical,gene_data,metagene_mean,cox_formula)} ,mc.cores=ncores)
	oo=t(as.data.frame(oo))
	rownames(oo)=gi
	return(oo)
}

constru_old<-function(clinical, gene_data, metagene_mean)
{
  # make surv object
  clinical_surv<-Surv(clinical$OS_Time, clinical$OS_Event)
  
  # create a categorical vector for metagene holding 3 values High, Med and Low
  metagene_tertiles<-as.numeric(metagene_mean)
  metagene_tertiles[metagene_mean <= quantile(metagene_mean,1:3/3)[[3]]]<-"High"
  metagene_tertiles[metagene_mean <= quantile(metagene_mean,1:3/3)[[2]]]<-"Med"
  metagene_tertiles[metagene_mean <= quantile(metagene_mean,1:3/3)[[1]]]<-"Low"
  
  
  # Define empty vector to hold tertiles for each gene in gene_data
  gene_factors<-c()
  
  # Define column names for the table that will hold results of constru
  column_names<-c("mean", "P5%", "P95%", "Diff P95%_5%", "R", "PS1", "PS2","PS3", 
       "gT1 Pval","gT1 HR", "gT1sT1", "gT1sT1 3 year", "gT1sT1 6 year",
       "gT1sT2", "gT1sT2 3 year", "gT1sT2 6 year", "gT1sT3", "gT1sT3 3 year", 
       "gT1sT3 6 year", "gT2 Pval","gT2 HR", "gT2sT1", "gT2sT1 3 year", 
       "gT2sT1 6 year", "gT2sT2", "gT2sT2 3 year", "gT2sT2 6 year", "gT2sT3", 
       "gT2sT3 3 year", "gT2sT3 6 year", "gT3 Pval","gT3 HR", "gT3sT1", 
       "gT3sT1 3 year", "gT3sT1 6 year", "gT3sT2", "gT3sT2 3 year", 
       "gT3sT2 6 year", "gT3sT3", "gT3sT3 3 year", "gT3sT3 6 year")
  
  # Define a data frame that will hold the results of the constru and 
  # assign column names
  df <- data.frame(row.names = NULL, 
                   matrix(ncol = length(column_names), nrow = 0))
  colnames(df) <- column_names
  
  # main for loop to compute the results
  for (k in 1:nrow(gene_data)) 
  { 
    # Tells users which gene is being processed 
#    incProgress(1/nrow(gene_data), detail = paste("Getting results for gene", k, "out of", nrow(gene_data), "genes"))
    
    # Create a vector to hold values for gene k
    hold <- as.numeric(gene_data[k,])
    
    # do constru only on genes that have less than 33% of zero entries
    if (100*length(which(hold==0))/length(hold) < 33)
    {
      # create tertiles for gene k
      gene_factors[gene_data[k,] <= quantile(hold, 1:3/3)[[3]]] <- "High"
      gene_factors[gene_data[k,] <= quantile(hold, 1:3/3)[[2]]] <- "Med"
      gene_factors[gene_data[k,] <= quantile(hold, 1:3/3)[[1]]] <- "Low"
      
      # subset metagene names/values based on gene k tertiles
      # subset clinical data based on gene k tertiles 
      # and run coxph() on each subset
      metagene_low <- metagene_mean[gene_factors %in% c("Low")]
      metagene_med <- metagene_mean[gene_factors %in% c("Med")]
      metagene_high <- metagene_mean[gene_factors %in% c("High")]
      
      clinical_low <- clinical_surv[gene_factors%in%c("Low")]
      clinical_med <- clinical_surv[gene_factors%in%c("Med")]
      clinical_high <- clinical_surv[gene_factors%in%c("High")]
      
      cox_low <- coxph(clinical_low ~ as.numeric(metagene_low))
      cox_med <- coxph(clinical_med ~ as.numeric(metagene_med))
      cox_high <- coxph(clinical_high ~ as.numeric(metagene_high))
      
      # now take metagene categorical vector consisting of -- Low, Med and High
      # this vector was created earlier
      # and subset based on gene k tertiles
      # run survfit on each subset
      metagene_tertiles_low<-metagene_tertiles[gene_factors %in% c("Low")]
      metagene_tertiles_med<-metagene_tertiles[gene_factors %in% c("Med")]
      metagene_tertiles_high<-metagene_tertiles[gene_factors %in% c("High")]
      
      surv_fit_low<-survfit(clinical_low~metagene_tertiles_low, 
                            type="kaplan-meier", conf.type="log")
      surv_fit_med<-survfit(clinical_med~metagene_tertiles_med, 
                            type="kaplan-meier", conf.type="log")
      surv_fit_high<-survfit(clinical_high~metagene_tertiles_high, 
                             type="kaplan-meier", conf.type="log")
      
      # write results for gene k 
      df[k,] <- c(mean(hold), quantile(hold, c(0.05)), quantile(hold, c(0.95)), 
                  quantile(hold, c(0.95)) - quantile(hold, c(0.05)),
                  cor(hold, as.numeric(metagene_mean), method="pearson"),
                  
                  # next three blocks of lines are parity scores 
                  # and they need some additional work
                  # they include only P-val and HR of the low/high tertiles
                  # I think gene-k variability and gene-k correlation to the
                  # metagene are also should be part of the parity score,
                  # if we want to extract new metagenes for the user without 
                  # the user going through the table themselves. Some users
                  # might prefer to have new metagene being pre-selected. 
                  (-log10(summary(cox_low)$coef[5])/summary(cox_low)$coef[2])-
                    (-log10(summary(cox_high)$coef[5])/summary(cox_high)$coef[2]),
                  
                  (log(summary(cox_low)$coef[5])*(summary(cox_low)$coef[2]-1))-
                    (log(summary(cox_high)$coef[5])*(summary(cox_high)$coef[2]-1)),
                  
                  (log(summary(cox_low)$coef[5])*log(summary(cox_low)$coef[2]))-
                    (log(summary(cox_high)$coef[5])*log(summary(cox_high)$coef[2])),
                  
                  summary(cox_low)$coef[5],summary(cox_low)$coef[2], 
                  
                  if(tryCatch(surv_fit_low["metagene_tertiles_low=Low"]$n, 
                              error=function(err) 0) !=0)
                  {
                    c(surv_fit_low["metagene_tertiles_low=Low"]$n,
                      summary(surv_fit_low["metagene_tertiles_low=Low"],
                              time=c(3,6), 
                              extend = TRUE)$surv)
                  }
                  else
                  {
                    c(0,0,0)
                  }
                  ,
                  
                  if(tryCatch(surv_fit_low["metagene_tertiles_low=Med"]$n, 
                              error=function(err) 0) !=0)
                  {
                    c(surv_fit_low["metagene_tertiles_low=Med"]$n,
                      summary(surv_fit_low["metagene_tertiles_low=Med"],
                              time=c(3,6), 
                              extend = TRUE)$surv)
                  }
                  else
                  {
                    c(0,0,0)
                  }
                  ,
                  
                  if(tryCatch(surv_fit_low["metagene_tertiles_low=High"]$n, 
                              error=function(err) 0) !=0)
                  {
                    c(surv_fit_low["metagene_tertiles_low=High"]$n,
                      summary(surv_fit_low["metagene_tertiles_low=High"],
                              time=c(3,6), 
                              extend = TRUE)$surv)
                  }
                  else
                  {
                    c(0,0,0)
                  }
                  ,
                  
                  summary(cox_med)$coef[5],summary(cox_med)$coef[2], 
                  
                  if(tryCatch(surv_fit_med["metagene_tertiles_med=Low"]$n, 
                              error=function(err) 0) !=0)
                  {
                    c(surv_fit_med["metagene_tertiles_med=Low"]$n,
                      summary(surv_fit_med["metagene_tertiles_med=Low"],
                              time=c(3,6),
                              extend = TRUE)$surv)
                  }
                  else
                  {
                    c(0,0,0)
                  }
                  ,
                  
                  if(tryCatch(surv_fit_med["metagene_tertiles_med=Med"]$n, 
                              error=function(err) 0) !=0)
                  {
                    c(surv_fit_med["metagene_tertiles_med=Med"]$n,
                      summary(surv_fit_med["metagene_tertiles_med=Med"],
                              time=c(3,6),
                              extend = TRUE)$surv)
                  }
                  else
                  {
                    c(0,0,0)
                  }
                  ,
                  
                  if(tryCatch(surv_fit_med["metagene_tertiles_med=High"]$n, 
                              error=function(err) 0) !=0)
                  {
                    c(surv_fit_med["metagene_tertiles_med=High"]$n,
                      summary(surv_fit_med["metagene_tertiles_med=High"],
                              time=c(3,6),
                              extend = TRUE)$surv)
                  }
                  else
                  {
                    c(0,0,0)
                  }
                  ,
                  
                  summary(cox_high)$coef[5],summary(cox_high)$coef[2], 
                  
                  if(tryCatch(surv_fit_high["metagene_tertiles_high=Low"]$n, 
                              error=function(err) 0) !=0)
                  {
                    c(surv_fit_high["metagene_tertiles_high=Low"]$n,
                      summary(surv_fit_high["metagene_tertiles_high=Low"],
                              time=c(3,6),
                              extend = TRUE)$surv)
                  }
                  else
                  {
                    c(0,0,0)
                  }
                  ,
                  
                  if(tryCatch(surv_fit_high["metagene_tertiles_high=Med"]$n, 
                              error=function(err) 0) !=0)
                  {
                    c(surv_fit_high["metagene_tertiles_high=Med"]$n,
                      summary(surv_fit_high["metagene_tertiles_high=Med"],
                              time=c(3,6),
                              extend = TRUE)$surv)
                  }
                  else
                  {
                    c(0,0,0)
                  }
                  ,
                  
                  if(tryCatch(surv_fit_high["metagene_tertiles_high=High"]$n, 
                              error=function(err) 0) !=0)
                  {
                    c(surv_fit_high["metagene_tertiles_high=High"]$n,
                      summary(surv_fit_high["metagene_tertiles_high=High"],
                              time=c(3,6),
                              extend = TRUE)$surv)
                  }
                  else
                  {
                    c(0,0,0)
                  }
                  ) 
     }
     else # genes with more than 33% of zero's are not part of the computations
     {
       df[k,]<-c(mean(hold), quantile(hold, c(0.05)), quantile(hold, c(0.95)), 
                 quantile(hold, c(0.95)) - quantile(hold, c(0.05)),
                 cor(hold, metagene_mean, method="pearson"),NA, NA, NA,
                 NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, 
                 NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, 
                 NA, NA, NA, NA, NA)
     }
    
  }  
  rownames(df)<-rownames(gene_data)
  return(df)
}

print("Multi-threaded 6:")
system.time(constru(clinical, gene_data, metagene_mean,"Surv( OS_Time, OS_Event) ~ metagene",6))
print("Multi-threaded 4:")
system.time(constru(clinical, gene_data, metagene_mean,"Surv( OS_Time, OS_Event) ~ metagene",4))
print("Multi-threaded 2:")
system.time(constru(clinical, gene_data, metagene_mean,"Surv( OS_Time, OS_Event) ~ metagene",2))
print("Single-threaded - Ming:")
system.time(constru(clinical, gene_data, metagene_mean,"Surv( OS_Time, OS_Event) ~ metagene",1))
print("Single-threaded - Julia:")
system.time(constru_old(clinical, gene_data, metagene_mean))


