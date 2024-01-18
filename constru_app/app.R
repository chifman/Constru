###################################
# Author: Julia Chifman
# Contact: chifman@american.edu
# Version: 1.0
# Year: 2024
###################################

library(shinydashboard)
library(shiny)
library(shinythemes)
library(doParallel)
library(ggplot2)
library(gplots)
library(DT)
library(RColorBrewer)
library(ComplexHeatmap)
library(survival)
library(survminer)
library(dplyr)
library(ggfortify)
library(gridExtra)
library(textshape)
library(patchwork)

# max size of file upload
options(shiny.maxRequestSize=1000*1024^2)
options(htmlwidgets.TOJSON_ARGS = list(na = 'string'))

# define main function constru()
# this function will take on three files that were uploaded by the user
# the files are: 
# (1) clinical data -- rows are patient ids and two of the columns must be 
# labeled for now as OS_Time and OS_Event. I will adjust this and make it smarter.
# (2) gene expression data -- columns should be patient ids and rows 
# gene identifiers (ids or gene names). Patient ids must match between clinical 
# and gene expression data.
# (3) metagene mean values. Users upload either gene ids/names and then app
# makes a metagene or upload metagene mean values. 

constru<-function(clinical, gene_data, metagene_mean)
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
    incProgress(1/nrow(gene_data), detail = paste("Getting results for gene", k, "out of", nrow(gene_data), "genes"))
    
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


# apps input
# ui <- fluidPage(
#     theme = bslib::bs_theme(bootswatch = "flatly"),
#     titlePanel("Constru"),
#     tags$style(HTML('table.dataTable tr.selected td, 
#                     table.dataTable td.selected 
#                     {background-color: #C6BAB8 !important;}')),
ui <- navbarPage("Constru",
                 theme =  shinytheme("flatly"),
                 fluid = T, 
                 tags$head(
                   tags$style(HTML("table {table-layout: auto;
                                    witdth: 80%;
                                    font-size: 80%;
                                    tab-size: 15px;
                                    }"
                                  )
                              )
                 ),
    #tabsetPanel(
        # Upload clinical data tab
        tabPanel("Import data", 
                 h4("Upload Clinical Data"),
                 p("Choose a file to upload clinical data. 
                   File should be tab delimited and have a header and row names. 
                   Row names should be patients ID. 
                   Column names should only use standard letters,
                   numbers and no spaces."),
                 fileInput("file1", "Upload file",
                           multiple = FALSE,
                           accept = c("text/csv",
                                      "text/comma-separated-values,text/plain",
                                      ".csv")),
                 h6("Data preview: clinical"),
            tableOutput("contents")
        ,
        
        # Upload second file: gene expression data
        # Horizontal line ----
        tags$hr(),
        h4("Upload Gene expression Data"),
        p("Choose a file to upload gene expression data in txt/csv format. 
                  Data must be tab separated and have a header. 
                  Make sure that patient IDs in
                  clinical and gene expression data are the same."),
        fileInput("file2", "Upload file",
                  multiple = FALSE,
                  accept = c("text/csv",
                             "text/comma-separated-values,text/plain",
                             ".csv")),
        
        h6("Data preview: gene expression"),
        tableOutput("contents2"),
        
        # Upload metagene
        # Horizontal line ----
        tags$hr(),
        h4("Upload Gene Signature"),
        p("Choose a file to upload metagene names or metagene mean values in 
            txt/csv format. Check the appropriate box below before uploading 
            a file."),
        p(strong("Gene names file:"), "one gene name/id per line. Gene names must match
          gene names/ids in your gene expression file. If using this option,
          then metagene can be visualized by clicking on the tab above labeled", 
          strong("Visualise metagene.")),
          
        p(strong("Mean values file:"), "a file containing one row of mean values 
          with columns labeled by patient IDs. Make sure that patient IDs in 
          clinical and metagene file are the same. Label the row contating 
          mean values as follows:", em("mean_values")),
        
        radioButtons("radio", h3("Choose one button"),
                     choices = list("Gene names file" = 1, "Mean values file" = 2),
                     selected = 1),
        
        fileInput("file3", "Upload file",
                  multiple = FALSE,
                  accept = c("text/csv",
                             "text/comma-separated-values,text/plain",
                             ".csv")),
        h6("Data preview: metagene"),
        
        
        # Output: Data file ----
        tableOutput("contents3"),
        )
    ,
    tabPanel("Visualise metagene",
             h4("Metagene heatmap"),
             p("This page displays a  heatmap of the metagene and also 
               its average value"),
             plotOutput("plot_metagene_mean")),
    
    tabPanel("Run Constru",
             h4("Table"),
             p("results of the analysis are displayed in the table below."),
             p("Column names will be defined here ..."),
             # Horizontal line ----
             tags$hr(),
             # actionButton("LowT", "Save lowerT genes", style="color: #135ADA;
             #                                                  background-color: #D7E4FD;
             #                                                  border-color: #135ADA;
             #                                                  font-size:120%;
             #                                                  border-width: 2px"),
             # actionButton("HighT", "Save upperT genes", style="color: #D35400;
             #                                                    background-color: #FDF2E9;
             #                                                    border-color: #D35400;
             #                                                    font-size:120%;
             #                                                    border-width: 2px"),
             actionButton("LowT", "Save lowerT genes"),
             actionButton("HighT", "Save upperT genes"),
             # Horizontal line ----
             tags$hr(),
             br(),
             DTOutput("constru_out")
    ),
    
  
  navbarMenu("Constru plots",
             tabPanel("LowerT Genes",
                     fluidRow(class="lowerT_plots", "Heatmap and KM-plot LowerT Genes",
                              plotOutput("plot_heatmap_lowT_2", height = "300px"),
                              plotOutput("plot_KM_lowT_2", height = "350px"))
                      ),
             tabPanel("UpperT Genes",
                      fluidRow(class="upperT_plots","Heatmap and KM-plot UpperT Genes",
                               plotOutput("plot_heatmap_highT_2", height = "300px"),
                               plotOutput("plot_KM_highT_2", height = "350px"))
                      ),
             tabPanel("UpperT - LowerT",
                      fluidRow(class="diff_plots","Heatmap and KM-plot (UpperT-LowerT) Genes",
                               plotOutput("plot_heatmap_diffT_2", height = "350px"),
                               plotOutput("plot_KM_diffT_2", height = "350px"))
                      )
             
  ),
  
    tabPanel("Documentation",
             h4("Constru Documentation"),
             p("Some text goes here ... "))
    #)
)

server <- function(input, output, session) {
    thematic::thematic_shiny()
  
    clinical<-reactive({
        req(input$file1)
        read.table(input$file1$datapath, header = TRUE, sep="\t")
    })
    
    gene_data<-reactive({
        req(input$file2)
        read.table(input$file2$datapath, header = TRUE)
    })
    
  
    metagene_mean<-reactive({
      req(input$file3)
        if (input$radio == 1)
        {
         metagene<-read.table(input$file3$datapath, header = FALSE)
         meta_mean<-colMeans(gene_data()[rownames(gene_data()) %in% metagene[,1],])
         return(meta_mean)
        }
       if (input$radio == 2)
       {
         metagene<-read.table(input$file3$datapath, header = TRUE)
         return(t(metagene))
       }
      
    })
    
    constru_table<-reactive({
      withProgress(message = "Calculation in progress ...",
                   constru(clinical(), gene_data(), metagene_mean()))
    })
    
    output$contents <- renderTable({
        
            return(head(clinical()))
        
    },rownames = TRUE)
        
         output$contents2 <- renderTable({
        
             return(gene_data()[1:4,1:7])
           
    },rownames = TRUE)
        
        output$contents3 <- renderTable({
            
            return(head(metagene_mean()))
            
        },rownames = FALSE, colnames = FALSE)
        

        # a custom table container that displays custom column names 
        # for constru results
        custom_colnames = htmltools::withTags(table(
          #class = 'display',
          thead(
            tr(
              th(colspan = 9, ''),
              th(colspan = 11, 'Gene Low Tertile', class = "dt-center"),
              th(colspan = 11, 'Gene Medium Tertile', class = "dt-center"),
              th(colspan = 11, 'Gene High Tertile', class = "dt-center")
            ),
            tr(
              th(colspan = 1, ''),
              th(colspan = 5, 'Gene attributes', class = "dt-center"),
              th(colspan = 3, 'Parity scores', class = "dt-center"),
              th(colspan = 2, ''),
              th(colspan = 3, 'Signature Low', class = "dt-center"),
              th(colspan = 3, 'Signature Medium', class = "dt-center"),
              th(colspan = 3, 'Signature High', class = "dt-center"),
              th(colspan = 2, ''),
              th(colspan = 3, 'Signature Low', class = "dt-center"),
              th(colspan = 3, 'Signature Medium', class = "dt-center"),
              th(colspan = 3, 'Signature High', class = "dt-center"),
              th(colspan = 2, ''),
              th(colspan = 3, 'Signature Low', class = "dt-center"),
              th(colspan = 3, 'Signature Medium', class = "dt-center"),
              th(colspan = 3, 'Signature High', class = "dt-center")
            ),
             tr(
               lapply(c('ID', 'Mean', '5%', "95%", 'Diff', "R", "PS1", "PS2", "PS3", rep(c('Pval', 'HR', 'size', '3Y', '6Y','size', '3Y', '6Y','size', '3Y', '6Y'),3)), th)
             )
          )
        ))
        print(custom_colnames)
    
    # visualize metagene via heatmap if names were uploaded
    output$plot_metagene_mean<-renderPlot({
      
      expression_data_subset<-gene_data()[rownames(gene_data()) %in% read.table(input$file3$datapath, header = FALSE)[,1],]
      
      ht = Heatmap(t(scale(t(as.matrix(rbind(expression_data_subset, metagene_mean()))))), 
                   cluster_rows = FALSE,
                   cluster_columns = TRUE,
                   show_column_names = FALSE,
                   row_labels = c(rownames(expression_data_subset), "Mean value"),
                   heatmap_legend_param = list( 
                                title = "Color key", 
                                title_gp = gpar(fontsize = 16),
                                direction = "horizontal",
                                legend_width = unit(75, "mm")),
                   row_split = c(rep("1", nrow(expression_data_subset)), rep("2",1)),
                   row_gap = unit(5,"mm"),
                   show_row_dend = FALSE,
                   col=rev(brewer.pal(11,"BrBG")))
      draw(ht, heatmap_legend_side = "bottom")

    })
    
    # display constru results
    output$constru_out <- renderDT({
    datatable(constru_table(),
    extensions = c("Buttons", "Select"),
    selection = 'none',
    filter = 'top',
    container = custom_colnames,
    options = list(
      search = list(regex = TRUE, caseInsensitive = TRUE),
      pageLength = 10,
      select = list(style = 'multiple'),
      #fixedHeader = TRUE,
      #width = "85%",
      #columnDefs = list(list(width = '10px', targets = "_all")),
      dom = 'Bfrtip',
      buttons = list('pageLength', 'colvis', 'copy','selectNone', #'selectRows',
                  list(extend = 'collection',
                  buttons = c('csv', 'excel'),
                  text = 'Download',
                  header = T)),
      #editable = TRUE,
      lengthMenu = list(c(10, 15, 25, 100, -1), c('10', '15', '25', '100','All')))) %>%
        formatSignif(1:ncol(constru_table()), digits = 3)%>%
        formatRound(c("gT1sT1","gT1sT2","gT1sT3","gT2sT1","gT2sT2","gT2sT3","gT3sT1","gT3sT2","gT3sT3"), 
                    digits=0)%>%
        formatStyle(1:6, color = "#000000")%>%
        formatStyle(c("PS1"), fontWeight="bold", color="#214d08")%>%
        formatStyle(c("PS2"), fontWeight="bold", color="#214d08")%>%
        formatStyle(c("PS3"), fontWeight="bold", color="#214d08")%>%
        formatStyle(c("gT1 Pval"), fontWeight="bold", color="#135ADA")%>%
        formatStyle(c("gT1 HR"), fontWeight="bold", color="#135ADA")%>%
        formatStyle(c("gT2 Pval"), fontWeight="bold", color="#717D7E")%>%
        formatStyle(c("gT2 HR"), fontWeight="bold", color="#717D7E")%>%
        formatStyle(c("gT3 Pval"), fontWeight="bold", color="#D35400")%>%
        formatStyle(c("gT3 HR"), fontWeight="bold", color="#D35400")%>%
        formatStyle(c("gT1sT1"), backgroundColor = "#D7E4FD", fontWeight="bold", color="#1419BE")%>%
        formatStyle(c("gT1sT1 3 year", "gT1sT1 6 year"), backgroundColor = "#D7E4FD", fontWeight="bold", color="#000000")%>%
        formatStyle(c("gT1sT2"), backgroundColor = "#BDD4FE", fontWeight="bold", color="#1419BE")%>%
        formatStyle(c("gT1sT2 3 year", "gT1sT2 6 year"), backgroundColor = "#BDD4FE", fontWeight="bold", color="#000000")%>%
        formatStyle(c("gT1sT3"), backgroundColor = "#97BBFC", fontWeight="bold", color="#1419BE")%>%
        formatStyle(c("gT1sT3 3 year", "gT1sT3 6 year"), backgroundColor = "#97BBFC", fontWeight="bold", color="#000000")%>%
        formatStyle(c("gT2sT1"), backgroundColor = "#F4F6F6", fontWeight="bold", color="#1419BE")%>%
        formatStyle(c("gT2sT1 3 year", "gT2sT1 6 year"), backgroundColor = "#F4F6F6", fontWeight="bold", color="#000000")%>%
        formatStyle(c("gT2sT2"), backgroundColor = "#EAEDED", fontWeight="bold", color="#1419BE")%>%
        formatStyle(c("gT2sT2 3 year", "gT2sT2 6 year"), backgroundColor = "#EAEDED", fontWeight="bold", color="#000000")%>%
        formatStyle(c("gT2sT3"), backgroundColor = "#D5DBDB", fontWeight="bold", color="#1419BE")%>%
        formatStyle(c("gT2sT3 3 year", "gT2sT3 6 year"), backgroundColor = "#D5DBDB", fontWeight="bold", color="#000000")%>%
        formatStyle(c("gT3sT1"), backgroundColor = "#FDF2E9", fontWeight="bold", color="#1419BE")%>%
        formatStyle(c("gT3sT1 3 year", "gT3sT1 6 year"), backgroundColor = "#FDF2E9", fontWeight="bold", color="#000000")%>%
        formatStyle(c("gT3sT2"), backgroundColor = "#FAE5D3", fontWeight="bold", color="#1419BE")%>%
        formatStyle(c("gT3sT2 3 year", "gT3sT2 6 year"), backgroundColor = "#FAE5D3", fontWeight="bold", color="#000000")%>%
        formatStyle(c("gT3sT3"), backgroundColor = "#F5CBA7", fontWeight="bold", color="#1419BE")%>%
        formatStyle(c("gT3sT3 3 year", "gT3sT3 6 year"), backgroundColor = "#F5CBA7", fontWeight="bold", color="#000000")
    }, 
    server = FALSE
    
    )
    
      # get row names (gene ids/names) when user clicks on rows and then button
      # labeled "Save lowerT genes"
      lowT_data <- eventReactive(input$LowT, {
        gene_data()[rownames(gene_data()) %in% rownames(constru_table()[input$constru_out_rows_selected,]),]
      })
      
      # get row names (gene ids/names) when user clicks on rows and then button
      # labeled "Save upperT genes"
      highT_data <- eventReactive(input$HighT, {
        gene_data()[rownames(gene_data()) %in% rownames(constru_table()[input$constru_out_rows_selected,]),]
      })
    
      # The rest of the code will be simplified using function() definitions.
      # For now it is repetitive, but this way I can catch errors better.

    
    #------------------------
    # Display LowT genes via heatmap 
    output$plot_heatmap_lowT_2 <- renderPlot({
      
      heatmap_data_means<-data.frame(meta_mean = colMeans(lowT_data()))
      heatmap_data_means %>%
        mutate(newcol=NA)
      heatmap_data_means_tert<-heatmap_data_means %>%
        mutate(newcol = ntile(meta_mean, 3)) %>%
        mutate(newcol = if_else(newcol == 1, 'Low', if_else(newcol == 2, 'Medium', 'High'))) 
      
      heatmap_data_low<-lowT_data()[,colnames(lowT_data()) %in% 
                                      rownames(subset(heatmap_data_means_tert, newcol=="Low"))]
      heatmap_data_med<-lowT_data()[, colnames(lowT_data()) %in% 
                                      rownames(subset(heatmap_data_means_tert, newcol=="Medium"))]
      heatmap_data_high<-lowT_data()[, colnames(lowT_data()) %in% 
                                       rownames(subset(heatmap_data_means_tert, newcol=="High"))]
      
      heatmap_data_low_clustered<-cluster_matrix(heatmap_data_low, 
                                                 dim = 'col', method = "average")
      heatmap_data_med_clustered<-cluster_matrix(heatmap_data_med, 
                                                 dim = 'col', method = "average")
      heatmap_data_high_clustered<-cluster_matrix(heatmap_data_high, 
                                                  dim = 'col', method = "average")
      
      heatmap_data_tert<-cbind(heatmap_data_low_clustered,
                               heatmap_data_med_clustered,
                               heatmap_data_high_clustered)
      
      
      ht_constru = Heatmap(t(scale(t(heatmap_data_tert))),
                           cluster_rows = TRUE,
                           cluster_columns = FALSE,
                           show_column_names = FALSE,
                           column_split = c(rep("1: Low Tertile", ncol(heatmap_data_low)),
                                            rep("2: Medium Tertile", ncol(heatmap_data_med)),
                                            rep("3: High Tertile", ncol(heatmap_data_high))),
                           row_labels = c(rownames(lowT_data())),
                           column_gap = unit(10,"mm"), 
                           column_title_gp = gpar(fill=c("#BDD4FE","#F4F6F6","#FDF2E9"),
                                                  col=c("#0B3B90", "#717D7E","#D35400"),
                                                  fontsize=16),
                           show_heatmap_legend = FALSE,
                           show_row_dend = FALSE,
                           col=colorpanel(32, "#0650d1", "#dcdce6", "#ed4507") 
      )
      
      draw(ht_constru)
    }) 
    
    # Display HighT genes via heatmap 
    output$plot_heatmap_highT_2 <- renderPlot({
      
      heatmap_data_means_highT<-data.frame(meta_mean = colMeans(highT_data()))
      heatmap_data_means_highT %>%
        mutate(newcol=NA)
      heatmap_data_means_tert_highT<-heatmap_data_means_highT %>%
        mutate(newcol = ntile(meta_mean, 3)) %>%
        mutate(newcol = if_else(newcol == 1, 'Low', if_else(newcol == 2, 'Medium', 'High'))) 
      
      heatmap_data_low_highT<-highT_data()[,colnames(highT_data()) %in% 
                                             rownames(subset(heatmap_data_means_tert_highT, newcol=="Low"))]
      heatmap_data_med_highT<-highT_data()[, colnames(highT_data()) %in% 
                                             rownames(subset(heatmap_data_means_tert_highT, newcol=="Medium"))]
      heatmap_data_high_highT<-highT_data()[, colnames(highT_data()) %in% 
                                              rownames(subset(heatmap_data_means_tert_highT, newcol=="High"))]
      
      heatmap_data_low_clustered_highT<-cluster_matrix(heatmap_data_low_highT, 
                                                       dim = 'col', method = "average")
      heatmap_data_med_clustered_highT<-cluster_matrix(heatmap_data_med_highT, 
                                                       dim = 'col', method = "average")
      heatmap_data_high_clustered_highT<-cluster_matrix(heatmap_data_high_highT, 
                                                        dim = 'col', method = "average")
      
      heatmap_data_tert_highT<-cbind(heatmap_data_low_clustered_highT,
                                     heatmap_data_med_clustered_highT,
                                     heatmap_data_high_clustered_highT)
      
      
      ht_constru_highT = Heatmap(t(scale(t(heatmap_data_tert_highT))),
                                 cluster_rows = TRUE,
                                 cluster_columns = FALSE,
                                 show_column_names = FALSE,
                                 column_split = c(rep("1: Low Tertile", ncol(heatmap_data_low_highT)),
                                                  rep("2: Medium Tertile", ncol(heatmap_data_med_highT)),
                                                  rep("3: High Tertile", ncol(heatmap_data_high_highT))),
                                 row_labels = c(rownames(highT_data())),
                                 column_gap = unit(10,"mm"), 
                                 column_title_gp = gpar(fill=c("#BDD4FE","#F4F6F6","#FDF2E9"),
                                                        col=c("#0B3B90", "#717D7E","#D35400"),
                                                        fontsize=16),
                                 show_heatmap_legend = FALSE,
                                 show_row_dend = FALSE,
                                 col=colorpanel(32, "#0650d1", "#dcdce6", "#ed4507") 
      )
      
      draw(ht_constru_highT)
    }) 
    
    # Display LowT and HighT genes via heatmap 
    output$plot_heatmap_diffT_2 <- renderPlot({
      
      heatmap_data_diffT<-rbind(lowT_data(), highT_data())
      heatmap_data_means<-data.frame(meta_mean = colMeans(highT_data())-colMeans(lowT_data()))
      heatmap_data_means %>%
        mutate(newcol=NA)
      heatmap_data_means_tert<-heatmap_data_means %>%
        mutate(newcol = ntile(meta_mean, 3)) %>%
        mutate(newcol = if_else(newcol == 1, 'Low', if_else(newcol == 2, 'Medium', 'High'))) 
      
      heatmap_data_low<- heatmap_data_diffT[,colnames(heatmap_data_diffT) %in% 
                                              rownames(subset(heatmap_data_means_tert, newcol=="Low"))]
      heatmap_data_med<- heatmap_data_diffT[, colnames( heatmap_data_diffT) %in% 
                                              rownames(subset(heatmap_data_means_tert, newcol=="Medium"))]
      heatmap_data_high<- heatmap_data_diffT[, colnames( heatmap_data_diffT) %in% 
                                               rownames(subset(heatmap_data_means_tert, newcol=="High"))]
      
      heatmap_data_low_clustered<-cluster_matrix(heatmap_data_low, 
                                                 dim = 'col', method = "average")
      heatmap_data_med_clustered<-cluster_matrix(heatmap_data_med, 
                                                 dim = 'col', method = "average")
      heatmap_data_high_clustered<-cluster_matrix(heatmap_data_high, 
                                                  dim = 'col', method = "average")
      
      heatmap_data_tert<-cbind(heatmap_data_low_clustered,
                               heatmap_data_med_clustered,
                               heatmap_data_high_clustered)
      
      
      ht_constru = Heatmap(t(scale(t(heatmap_data_tert))),
                           cluster_rows = FALSE,
                           cluster_columns = FALSE,
                           show_column_names = FALSE,
                           column_split = c(rep("1: Low Tertile", ncol(heatmap_data_low)),
                                            rep("2: Medium Tertile", ncol(heatmap_data_med)),
                                            rep("3: High Tertile", ncol(heatmap_data_high))),
                           row_split = c(rep("lowerT",nrow(lowT_data())),
                                         rep("upperT",nrow(highT_data()))),
                           row_labels = c(rownames( heatmap_data_diffT)),
                           row_gap=unit(5,"mm"),
                           column_gap = unit(10,"mm"), 
                           row_title_gp = gpar(fill=c("#BDD4FE", "#FDF2E9"),
                                               col=c("#0B3B90","#D35400"),
                                               fontsize=14),
                           column_title_gp = gpar(fill=c("#BDD4FE","#F4F6F6","#FDF2E9"),
                                                  col=c("#0B3B90", "#717D7E","#D35400"),
                                                  fontsize=16),
                           heatmap_legend_param = list(
                             title = "Heatmap Color key",
                             title_gp = gpar(fontsize = 14),
                             direction = "horizontal",
                             legend_width = unit(75, "mm")),
                           show_row_dend = FALSE,
                           col=colorpanel(32, "#0650d1", "#dcdce6", "#ed4507") 
      )
      
      draw(ht_constru, heatmap_legend_side = "top")
    }) 
    
    # KM plots for LowT
    output$plot_KM_lowT_2 <-renderPlot({
      
      heatmap_data_means<-data.frame(meta_mean = colMeans(lowT_data())) 
      heatmap_data_means %>%
        mutate(newcol=NA)
      heatmap_data_means_tert<-heatmap_data_means %>%
        mutate(newcol = ntile(meta_mean, 3)) %>%
        mutate(newcol = if_else(newcol == 1, 'Low', if_else(newcol == 2, 'Medium', 'High')))
      
      #KM plots
      clinical_constru_low<-clinical()[rownames(clinical()) %in% 
                                         rownames(subset(heatmap_data_means_tert, newcol=="Low")),]
      clinical_constru_med<-clinical()[rownames(clinical()) %in% 
                                         rownames(subset(heatmap_data_means_tert, newcol=="Medium")),]
      clinical_constru_high<-clinical()[rownames(clinical()) %in% 
                                          rownames(subset(heatmap_data_means_tert, newcol=="High")),]
      
      
      # Split metagene into tertiles based on radio button selected -- issue with
      # reading single row -- it must be transposed if button 2 is selected
      # will fix that. Now it is just a quick way to display and see if 
      # everything makes sense 
      if (input$radio == 1)
      {
        metagene_mean_constru<-data.frame(meta_mean = metagene_mean())
        metagene_mean_constru %>%
          mutate(newcol=NA)
        metagene_mean_constru_tert<-metagene_mean_constru %>%
          mutate(newcol = ntile(meta_mean, 3)) %>%
          mutate(newcol = if_else(newcol == 1, 'Low', if_else(newcol == 2, 'Medium', 'High')))
      }
      
      if (input$radio == 2)
      {
        metagene_mean_constru<-data.frame(meta_mean = metagene_mean())
        metagene_mean_constru %>%
          mutate(newcol=NA)
        metagene_mean_constru_tert<-metagene_mean_constru %>%
          mutate(newcol = ntile(mean_values, 3)) %>%
          mutate(newcol = if_else(newcol == 1, 'Low', if_else(newcol == 2, 'Medium', 'High')))
      }
      
      
      metagene_mean_constru_low<-metagene_mean_constru_tert[rownames(metagene_mean_constru_tert) %in%
                                                              rownames(subset(heatmap_data_means_tert, newcol=="Low")),]
      metagene_mean_constru_med<-metagene_mean_constru_tert[rownames(metagene_mean_constru_tert) %in%
                                                              rownames(subset(heatmap_data_means_tert, newcol=="Medium")),]
      metagene_mean_constru_high<-metagene_mean_constru_tert[rownames(metagene_mean_constru_tert) %in%
                                                               rownames(subset(heatmap_data_means_tert, newcol=="High")),]
      
      
      data_for_fit_low<-merge(clinical_constru_low,metagene_mean_constru_low, by="row.names", all=TRUE)
      data_for_fit_med<-merge(clinical_constru_med,metagene_mean_constru_med, by="row.names", all=TRUE)
      data_for_fit_high<-merge(clinical_constru_high,metagene_mean_constru_high, by="row.names", all=TRUE)
      
      surv_fit_constru_low<-do.call(survfit, list(formula=Surv(OS_Time, OS_Event)~newcol, data=data_for_fit_low))
      surv_fit_constru_med<-do.call(survfit, list(formula=Surv(OS_Time, OS_Event)~newcol, data=data_for_fit_med))
      surv_fit_constru_high<-do.call(survfit, list(formula=Surv(OS_Time, OS_Event)~newcol, data=data_for_fit_high))
      
      print(surv_pvalue(surv_fit_constru_low))
      print(surv_pvalue(surv_fit_constru_med))
      print(surv_pvalue(surv_fit_constru_high))
      
      plot_km_1<-ggsurvplot(surv_fit_constru_low,conf.int = FALSE, pval=TRUE, legend = c(0.8, 0.8)) 
      plot_km_2<-ggsurvplot(surv_fit_constru_med,conf.int = FALSE, pval=TRUE, legend = c(0.8, 0.8))
      plot_km_3<-ggsurvplot(surv_fit_constru_high,conf.int = FALSE, pval=TRUE, legend = c(0.8, 0.8))
      
      arrange_ggsurvplots(list(plot_km_1, plot_km_2, plot_km_3),ncol=3)
      
    })
    
    # KM plots for HighT
    output$plot_KM_highT_2 <-renderPlot({
      
      heatmap_data_means_highT<-data.frame(meta_mean = colMeans(highT_data())) 
      heatmap_data_means_highT %>%
        mutate(newcol=NA)
      heatmap_data_means_tert_highT<-heatmap_data_means_highT %>%
        mutate(newcol = ntile(meta_mean, 3)) %>%
        mutate(newcol = if_else(newcol == 1, 'Low', if_else(newcol == 2, 'Medium', 'High')))
      
      #KM plots
      clinical_constru_low_highT<-clinical()[rownames(clinical()) %in% 
                                               rownames(subset(heatmap_data_means_tert_highT, newcol=="Low")),]
      clinical_constru_med_highT<-clinical()[rownames(clinical()) %in% 
                                               rownames(subset(heatmap_data_means_tert_highT, newcol=="Medium")),]
      clinical_constru_high_highT<-clinical()[rownames(clinical()) %in% 
                                                rownames(subset(heatmap_data_means_tert_highT, newcol=="High")),]
      
      
      # Split metagene into tertiles based on radio button selected -- issue with
      # reading single row -- it must be transposed if button 2 is selected
      # will fix that. Now it is just a quick way to display and see if 
      # everything makes sense 
      if (input$radio == 1)
      {
        metagene_mean_constru_highT<-data.frame(meta_mean = metagene_mean())
        metagene_mean_constru_highT %>%
          mutate(newcol=NA)
        metagene_mean_constru_tert_highT<-metagene_mean_constru_highT %>%
          mutate(newcol = ntile(meta_mean, 3)) %>%
          mutate(newcol = if_else(newcol == 1, 'Low', if_else(newcol == 2, 'Medium', 'High')))
      }
      
      if (input$radio == 2)
      {
        metagene_mean_constru_highT<-data.frame(meta_mean = metagene_mean())
        metagene_mean_constru_highT %>%
          mutate(newcol=NA)
        metagene_mean_constru_tert_highT<-metagene_mean_constru_highT %>%
          mutate(newcol = ntile(mean_values, 3)) %>%
          mutate(newcol = if_else(newcol == 1, 'Low', if_else(newcol == 2, 'Medium', 'High')))
      }
      
      metagene_mean_constru_low_highT<-metagene_mean_constru_tert_highT[rownames(metagene_mean_constru_tert_highT) %in%
                                                                          rownames(subset(heatmap_data_means_tert_highT, newcol=="Low")),]
      metagene_mean_constru_med_highT<-metagene_mean_constru_tert_highT[rownames(metagene_mean_constru_tert_highT) %in%
                                                                          rownames(subset(heatmap_data_means_tert_highT, newcol=="Medium")),]
      metagene_mean_constru_high_highT<-metagene_mean_constru_tert_highT[rownames(metagene_mean_constru_tert_highT) %in%
                                                                           rownames(subset(heatmap_data_means_tert_highT, newcol=="High")),]
      
      
      data_for_fit_low_highT<-merge(clinical_constru_low_highT,metagene_mean_constru_low_highT, by="row.names", all=TRUE)
      data_for_fit_med_highT<-merge(clinical_constru_med_highT,metagene_mean_constru_med_highT, by="row.names", all=TRUE)
      data_for_fit_high_highT<-merge(clinical_constru_high_highT,metagene_mean_constru_high_highT, by="row.names", all=TRUE)
      
      surv_fit_constru_low_highT<-do.call(survfit, list(formula=Surv(OS_Time, OS_Event)~newcol, data=data_for_fit_low_highT))
      surv_fit_constru_med_highT<-do.call(survfit, list(formula=Surv(OS_Time, OS_Event)~newcol, data=data_for_fit_med_highT))
      surv_fit_constru_high_highT<-do.call(survfit, list(formula=Surv(OS_Time, OS_Event)~newcol, data=data_for_fit_high_highT))
      
      print(surv_pvalue(surv_fit_constru_low_highT))
      print(surv_pvalue(surv_fit_constru_med_highT))
      print(surv_pvalue(surv_fit_constru_high_highT))
      
      plot_km_1_highT<-ggsurvplot(surv_fit_constru_low_highT,conf.int = FALSE, pval=TRUE, legend = c(0.8, 0.8)) 
      plot_km_2_highT<-ggsurvplot(surv_fit_constru_med_highT,conf.int = FALSE, pval=TRUE, legend = c(0.8, 0.8))
      plot_km_3_highT<-ggsurvplot(surv_fit_constru_high_highT,conf.int = FALSE, pval=TRUE, legend = c(0.8, 0.8))
      
      arrange_ggsurvplots(list(plot_km_1_highT, plot_km_2_highT, plot_km_3_highT),ncol=3)
      
    })
    
    # KM plots for diffT
    output$plot_KM_diffT_2 <-renderPlot({
      
      heatmap_data_means_diffT<-data.frame(meta_mean = colMeans(highT_data())-colMeans(lowT_data())) 
      heatmap_data_means_diffT %>%
        mutate(newcol=NA)
      heatmap_data_means_tert_diffT<-heatmap_data_means_diffT %>%
        mutate(newcol = ntile(meta_mean, 3)) %>%
        mutate(newcol = if_else(newcol == 1, 'Low', if_else(newcol == 2, 'Medium', 'High')))
      
      #KM plots
      clinical_constru_low_diffT<-clinical()[rownames(clinical()) %in% 
                                               rownames(subset(heatmap_data_means_tert_diffT, newcol=="Low")),]
      clinical_constru_med_diffT<-clinical()[rownames(clinical()) %in% 
                                               rownames(subset(heatmap_data_means_tert_diffT, newcol=="Medium")),]
      clinical_constru_high_diffT<-clinical()[rownames(clinical()) %in% 
                                                rownames(subset(heatmap_data_means_tert_diffT, newcol=="High")),]
      
      
      # Split metagene into tertiles based on radio button selected -- issue with
      # reading single row -- it must be transposed if button 2 is selected
      # will fix that. Now it is just a quick way to display and see if 
      # everything makes sense 
      if (input$radio == 1)
      {
        metagene_mean_constru_diffT<-data.frame(meta_mean = metagene_mean())
        metagene_mean_constru_diffT %>%
          mutate(newcol=NA)
        metagene_mean_constru_tert_diffT<-metagene_mean_constru_diffT %>%
          mutate(newcol = ntile(meta_mean, 3)) %>%
          mutate(newcol = if_else(newcol == 1, 'Low', if_else(newcol == 2, 'Medium', 'High')))
      }
      
      if (input$radio == 2)
      {
        metagene_mean_constru_diffT<-data.frame(meta_mean = metagene_mean())
        metagene_mean_constru_diffT %>%
          mutate(newcol=NA)
        metagene_mean_constru_tert_diffT<-metagene_mean_constru_diffT %>%
          mutate(newcol = ntile(mean_values, 3)) %>%
          mutate(newcol = if_else(newcol == 1, 'Low', if_else(newcol == 2, 'Medium', 'High')))
      }
      
      metagene_mean_constru_low_diffT<-metagene_mean_constru_tert_diffT[rownames(metagene_mean_constru_tert_diffT) %in%
                                                                          rownames(subset(heatmap_data_means_tert_diffT, newcol=="Low")),]
      metagene_mean_constru_med_diffT<-metagene_mean_constru_tert_diffT[rownames(metagene_mean_constru_tert_diffT) %in%
                                                                          rownames(subset(heatmap_data_means_tert_diffT, newcol=="Medium")),]
      metagene_mean_constru_high_diffT<-metagene_mean_constru_tert_diffT[rownames(metagene_mean_constru_tert_diffT) %in%
                                                                           rownames(subset(heatmap_data_means_tert_diffT, newcol=="High")),]
      
      
      data_for_fit_low_diffT<-merge(clinical_constru_low_diffT,metagene_mean_constru_low_diffT, by="row.names", all=TRUE)
      data_for_fit_med_diffT<-merge(clinical_constru_med_diffT,metagene_mean_constru_med_diffT, by="row.names", all=TRUE)
      data_for_fit_high_diffT<-merge(clinical_constru_high_diffT,metagene_mean_constru_high_diffT, by="row.names", all=TRUE)
      
      surv_fit_constru_low_diffT<-do.call(survfit, list(formula=Surv(OS_Time, OS_Event)~newcol, data=data_for_fit_low_diffT))
      surv_fit_constru_med_diffT<-do.call(survfit, list(formula=Surv(OS_Time, OS_Event)~newcol, data=data_for_fit_med_diffT))
      surv_fit_constru_high_diffT<-do.call(survfit, list(formula=Surv(OS_Time, OS_Event)~newcol, data=data_for_fit_high_diffT))
      
      print(surv_pvalue(surv_fit_constru_low_diffT))
      print(surv_pvalue(surv_fit_constru_med_diffT))
      print(surv_pvalue(surv_fit_constru_high_diffT))
      
      plot_km_1_diffT<-ggsurvplot(surv_fit_constru_low_diffT,conf.int = FALSE, pval=TRUE, legend = c(0.8, 0.8))
      plot_km_2_diffT<-ggsurvplot(surv_fit_constru_med_diffT,conf.int = FALSE, pval=TRUE, legend = c(0.8, 0.8))
      plot_km_3_diffT<-ggsurvplot(surv_fit_constru_high_diffT,conf.int = FALSE, pval=TRUE, legend = c(0.8, 0.8))
      
      arrange_ggsurvplots(list(plot_km_1_diffT, plot_km_2_diffT, plot_km_3_diffT),ncol=3)
      
    })
    #------------------
    
    session$onSessionEnded(function() {
       stopApp()
   })
    
}

# Run the application 
shinyApp(ui = ui, server = server)
