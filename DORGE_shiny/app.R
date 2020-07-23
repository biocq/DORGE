if(!require("data.table")) install.packages("data.table", repos = "http://cran.us.r-project.org")
if(!require("plotly")) install.packages("plotly", repos = "http://cran.us.r-project.org")
if(!require("shiny")) install.packages("shiny", repos = "http://cran.us.r-project.org")
if(!require("DT")) install.packages("DT", repos = "http://cran.us.r-project.org")
if(!require("plyr")) install.packages("plyr", repos = "http://cran.us.r-project.org")
if(!require("dplyr")) install.packages("dplyr", repos = "http://cran.us.r-project.org")
if(!require("reshape2")) install.packages("reshape2", repos = "http://cran.us.r-project.org")
if(!require("RColorBrewer")) install.packages("RColorBrewer", repos = "http://cran.us.r-project.org")
if(!require("stringr")) install.packages("stringr", repos = "http://cran.us.r-project.org")
if(!require("ggplot2")) install.packages("ggplot2", repos = "http://cran.us.r-project.org")
if(!require("heatmaply")) install.packages("heatmaply", repos = "http://cran.us.r-project.org")

#install.packages(c("data.table","plotly","DT","plyr","dplyr","reshape2","RColorBrewer","stringr","ggplot2","heatmaply"), repos = "http://cran.us.r-project.org")

suppressMessages(library(data.table))
suppressMessages(library(plotly))
suppressMessages(library(shiny))
suppressMessages(library(DT))
suppressMessages(library(plyr))
suppressMessages(library(dplyr))
suppressMessages(library(reshape2))
suppressMessages(library(RColorBrewer))
suppressMessages(library(stringr))
suppressMessages(library(ggplot2))
suppressMessages(library(heatmaply))

options(warn=-1)

rm(list = ls(), envir = environment())

#setwd("/Users/jlyu/Box Sync/TSGOG_Project/SA_sub/github/DORGE_paper/DORGE_shiny")

# Make an empty Google analytics file (for dev version - not for production)
if (!file.exists('google-analytics.js'))
{
  file.create('google-analytics.js')
}

H3K4me3Filename <- 'data/ENCODE_H3K4me3_peak_info.txt'
H3K4me3Filename <- normalizePath(H3K4me3Filename)

H3K4me3_sample_info <- 'data/H3K4me3_sample_information_for_website.txt'
H3K4me3_sample_info <- normalizePath(H3K4me3_sample_info)

canyonFilename <- 'data/Canyon_info.txt'
canyonFilename <- normalizePath(canyonFilename)

TCGAgeneexpressionFilename <- 'data/geneexpression_info.txt'
TCGAgeneexpressionFilename <- normalizePath(TCGAgeneexpressionFilename)

CCLEgeneexpressionFilename <- 'data/CCLE_RNASeq_matrix_melt.txt'
CCLEgeneexpressionFilename <- normalizePath(CCLEgeneexpressionFilename)

CCLEsampleannotationFilename <- 'data/CCLE_RNASeq_matrix_melt_anno.txt'
CCLEsampleannotationFilename <- normalizePath(CCLEsampleannotationFilename)

All_features <- 'data/All_features.csv'
All_features <- normalizePath(All_features)

gene_info <- 'data/gene_info.txt'
gene_info <- normalizePath(gene_info)

DORGE_prediction <- 'data/DORGE_prediction.txt'
DORGE_prediction <- normalizePath(DORGE_prediction)

fileInfo <- file.info(All_features)
modifiedDate <- strsplit(as.character(fileInfo$mtime), ' ')[[1]][1]

merged_features_table <- fread(All_features,sep=',',header=T,stringsAsFactors=TRUE)
DORGE_prediction_table <- fread(DORGE_prediction,sep='\t',header=T,stringsAsFactors=TRUE)
gene_info_table <- fread(gene_info,sep='\t',header=T,stringsAsFactors=TRUE)
merged_features_table<-cbind(gene_info_table,DORGE_prediction_table[,-1],merged_features_table[,-1])
TSG_threshold<-0.6233374 #FPR=0.01
OG_threshold<-0.6761319 #FPR=0.01
pvalues_TSG <- merged_features_table[,c("Gene","TSG_probability")]
pvalues_OG <- merged_features_table[,c("Gene","OG_probability")]
prioritized_genes<-unique(c(as.character(pvalues_TSG$Gene[pvalues_TSG$TSG_probability>TSG_threshold])),c(as.character(pvalues_OG$Gene[pvalues_OG$OG_probability>TSG_threshold])))
H3K4me3_peak_table <- fread(H3K4me3Filename,sep='\t',header=T,stringsAsFactors=TRUE)
H3K4me3_sample_info_table <- fread(H3K4me3_sample_info,sep='\t',header=T,stringsAsFactors=TRUE)
canyon_table <- fread(canyonFilename,sep='\t',header=T,stringsAsFactors=TRUE)
TCGAgeneexpression_table <- fread(TCGAgeneexpressionFilename,sep='\t',header=T,stringsAsFactors=TRUE)
CCLEgeneexpression_table <- fread(CCLEgeneexpressionFilename,sep='\t',header=T,stringsAsFactors=TRUE)
CCLE_sample_info_table <- fread(CCLEsampleannotationFilename,sep='\t',header=T,stringsAsFactors=TRUE)
CCLEgeneexpression_table_with_sample_info<-join(CCLEgeneexpression_table, CCLE_sample_info_table,type = "left",by="Sample_ID")

H3K4me3_peak_table_with_sample_info<-join(H3K4me3_peak_table, H3K4me3_sample_info_table,type = "left",by="Sample_name")
H3K4me3_peak_table_with_sample_info <- H3K4me3_peak_table_with_sample_info[order(-H3K4me3_peak_table_with_sample_info$Peak_length),]

H3K4me3_Explanation <- "<br /><br /><br /><b>H3K4me3 length (broad peaks) for top cell lines:</b><br />"
Canyon_Explanation <- "<br /><br /><br /><b>Methylation of Canyon genes for top cell lines:</b><br />"
TCGAgeneexpression_Explanation <- "<br /><b>Top expressed genes in TCGA Tumor samples:</b><br /><br />"
Geneexpression_CCLE_Explanation <- "<br /><b><a name='cellline_exp_barchart'>Top expressed genes in CCLE cell lines:</b><br />"
iPOT_table_Explanation <- "<br /><b>DORGE prediction:</b><br /><br />"
iPOT_customizing_Explanation <- "<br /><b>Prediction by predefined weights:</b><br /><br />"
TSG_feature_Explanation <- "<br /><b>TSG feature weights:</b><br /><br />"
OG_feature_Explanation <- "<br /><b>OG feature weights:</b><br /><br />"
feature_contribution_Explanation <- "<br /><b>Features:</b><br /><br /><br /><br />"

geneNames <- as.character(merged_features_table$Gene)
tissueNames <- unique(as.character(TCGAgeneexpression_table$Tissue))
celllineNames <- unique(as.character(CCLEgeneexpression_table_with_sample_info$Sample_illustration))

ui <- function(req) {
  fluidPage(
    tags$head(
      includeHTML("google-analytics.js"),
      includeHTML("metadata.html"),
      tags$style(".rightAlign{float:right; margin-left:5px; margin-bottom: 20px;}"),
      tags$style(".span1{ display: inline-block;vertical-align:top;width: 200px; }"),
      tags$style(".span2{ width: 1400px; }"),
      tags$style(".row-fluid{ display: inline-block;vertical-align:top;width: 1400px;}")
    ),
    titlePanel(
    windowTitle = "DORGE", 
	    fluidRow(
		    column(9, "Discovery of Oncogenes and Tumor SuppressoR Genes"), 
		    column(3, img(src = "DORGE_logo.svg",height = 50, width = 80,class = "pull-right"))
	  	)
    ),
    
    tabsetPanel(id="maintabs",
                type = "tabs",
                tabPanel("About", 
                         includeHTML("about.html"),
                         helpText(includeHTML("cc0.html"))
                ),
                tabPanel("DORGE prediction table", 
                         mainPanel(
                           HTML(iPOT_table_Explanation),
                           div(DT::dataTableOutput("iPOT_Table"), style = "font-size: 90%; width: 90%"),
                           verbatimTextOutput('clicked'),
                           helpText(paste("Last updated on:",modifiedDate)),
                           helpText(includeHTML("cc0.html"))
                         )
                ),tabPanel("by Tissue or Cell line", 
                         sidebarPanel(
                           selectizeInput("tissue_input", "Tissue", tissueNames, selected = 'Breast', multiple = FALSE, options = list(maxOptions = 2*length(tissueNames))),
                           selectizeInput("cellline_input", "Cell line", celllineNames, selected = 'AU565_BREAST', multiple = FALSE, options = list(maxOptions = 2*length(celllineNames))),
                           checkboxInput("topgenesonly", "show only DORGE-predicted genes", TRUE),
                           width=3
                         ),
                         mainPanel(
                           HTML(TCGAgeneexpression_Explanation),
                           plotlyOutput("tissue_exp_barchart"),
                           downloadButton("tissue_exp_download_collated_all", label = "Download All", class='rightAlign'),
                           #downloadButton("tissue_exp_download_collated_shown", label = "Download shown", class='rightAlign'),
                           DT::dataTableOutput("tissue_exp_table"),
                           
                           HTML(Geneexpression_CCLE_Explanation),
                           plotlyOutput("cellline_exp_barchart"),
                           downloadButton("cellline_exp_download_collated_all", label = "Download All", class='rightAlign'),
                           #downloadButton("cellline_exp_download_collated_shown", label = "Download shown", class='rightAlign'),
                           DT::dataTableOutput("cellline_exp_table"),
                           
                           helpText(paste("Last updated on:",modifiedDate)),
                           helpText(includeHTML("cc0.html"))
                         )
                ),
                tabPanel("by Gene", 
                         sidebarPanel(
                           selectizeInput("gene_input", "Gene", geneNames, selected = '', multiple = FALSE, options = list(maxOptions = 2*length(geneNames))),
                           HTML("Quick examples: HOXC6, PTEN, BCL2"),
                           HTML("<b>DORGE scores:</b>"),
                           plotlyOutput("gene_overview"),
                           width=3
                         ),
                         mainPanel(
                         	 htmlOutput("gene_information_by_gene_page"),
													 
                           HTML(H3K4me3_Explanation),
                           plotlyOutput("gene_barchart"),
                           downloadButton("gene_download_collated_all", label = "Download unaggregated H3K4me3 peak length information", class='rightAlign'),
                           downloadButton("gene_download_collated_shown", label = "Download mean H3K4me3 peak length for cell-lines", class='rightAlign'),
                           DT::dataTableOutput("gene_table"),
                           
                           HTML(Canyon_Explanation),
                           plotlyOutput("gene_barchart2"),
                           downloadButton("gene_download_collated_all2", label = "Download unaggregated gene-body methylation information", class='rightAlign'),
                           downloadButton("gene_download_collated_shown2", label = "Download mean gene-body methylation for cell-lines", class='rightAlign'),
                           DT::dataTableOutput("gene_table2"),
                           helpText(paste("Last updated on:",modifiedDate)),
                           helpText(includeHTML("cc0.html"))
                         )
                ),
                tabPanel("Help", 
                         includeHTML("help.html"),
                         helpText(includeHTML("cc0.html"))
                )
    ),
    helpText(includeHTML("subheader.html"))
  )
}

server <- function(input, output, session) {
    
  tissue_gene_exp_Table <- reactive({
    
    if (input$topgenesonly) {
    		table <- TCGAgeneexpression_table[TCGAgeneexpression_table$Tissue==input$tissue_input & TCGAgeneexpression_table$Gene %in% prioritized_genes,c(2,3)]
    		
    }else{
    		table <- TCGAgeneexpression_table[TCGAgeneexpression_table$Tissue==input$tissue_input,c(2,3)]
    }
    table
  })
  
  cellline_gene_exp_Table <- reactive({
    
    if (input$topgenesonly) {
    		table <- CCLEgeneexpression_table_with_sample_info[CCLEgeneexpression_table_with_sample_info$Sample_illustration==input$cellline_input & CCLEgeneexpression_table_with_sample_info$Gene %in% prioritized_genes,]
      	
    }else{
    		table <- CCLEgeneexpression_table_with_sample_info[CCLEgeneexpression_table_with_sample_info$Sample_illustration==input$cellline_input,]
    }
    table
  })
  
  geneTableProxy = dataTableProxy('gene_table')

  geneData_H3K4me3 <- reactive({
    table <- H3K4me3_peak_table_with_sample_info[H3K4me3_peak_table_with_sample_info$Gene==input$gene_input,]
    #table <- table[order(-table$Peak_length),]
    table
  })
  
  geneData_canyon <- reactive({
    table <- canyon_table[canyon_table$Gene==input$gene_input,]
    table <- table[order(-table$WGBS_methylation),]
    table
  })
  
  output$gene_information_by_gene_page <- renderUI({
	  specific_gene_info <- geneinfo()
	  	  	
	  if(nrow(specific_gene_info)>0){
		  table <- specific_gene_info[,1:20]

			genomelocation= ifelse(is.na(table$start),"none",paste0(table$chr,":", table$start,"-",table$end," (",table$chain,")"))
			
			prediction <- ""
			if(table$TSG_probability>TSG_threshold & table$OG_probability<OG_threshold){prediction<-"Predicted TSG"}
			if(table$OG_probability>OG_threshold & table$TSG_probability<TSG_threshold){prediction<-"Predicted OG"}
			if(table$TSG_probability>TSG_threshold & table$OG_probability>OG_threshold){prediction<-"Predicted TSG/OG"}
			CGC_status <- as.character(table$CGC_status)
			if(is.na(CGC_status)){prediction<-CGC_status}
			
			Tier_TSG <- "Not predicted"
			if(table$TSG_probability>TSG_threshold){Tier_TSG<-"5"}
			if(table$TSG_probability>0.7){Tier_TSG<-"4"}
			if(table$TSG_probability>0.8){Tier_TSG<-"3"}
			if(table$TSG_probability>0.9){Tier_TSG<-"2"}
			if(table$TSG_probability>0.95){Tier_TSG<-"1"}
			
			Tier_OG <- "Not predicted"
			if(table$OG_probability>OG_threshold){Tier_OG<-"5"}
			if(table$OG_probability>0.8){Tier_OG<-"4"}
			if(table$OG_probability>0.85){Tier_OG<-"3"}
			if(table$OG_probability>0.9){Tier_OG<-"2"}
			if(table$OG_probability>0.95){Tier_OG<-"1"}

			tagList(
		  	strong("Gene information:"),
		  	br(),
		  	strong(table$Gene),
		  	",",
			  table$Gene_illustration,
		  	strong(" Gene type: "),
		  	table$Gene_type,
		  	br(),
		  	strong("   Alias:"),
		  	table$Alias,
		  	strong("Position:"),
		  	genomelocation,
		  	strong(" Loci:"),
		  	table$Gene_loci,
		  	br(),
		  	strong("Ensembl Gene ID: "),
		  	table$Ensembl_gene_loci,
		  	strong(" RefSeq Gene ID: "),
		  	table$RefSeq,
		  	br(),
		  	strong("PubMed ID: "),
		  	table$Pubmed,
		  	br(),
		  	strong("DORGE Classification: "),
		  	prediction,
		  	strong("Tier TSG: "),
		  	Tier_TSG,
		  	strong(" Tier OG: "),
		  	Tier_OG,
		  	strong(" TSG-score: "),
		  	table$TSG_probability,
		  	strong(" OG-score: "),
		  	table$OG_probability,
		  	br(),
		  	strong("Other resources:"),
		  	br(),
			  "COSMIC (a resource exploring the impact of somatic mutations in human cancer):",
		  	a(input$gene_input,href=paste("https://cancer.sanger.ac.uk/cosmic/gene/analysis?ln=",input$gene_input, sep=""), target="_blank"),
		  	br(),
		  	"CancerMine (a literature database of drivers, TSGs and OGs in cancer):",
		  	a(input$gene_input,href=paste("http://bionlp.bcgsc.ca/cancermine/","", sep=""), target="_blank")
	  	)
	  }
	  
  })
  
  
  geneinfo <- reactive({
    specific_gene_info<- merged_features_table[merged_features_table$Gene == input$gene_input,]
    specific_gene_info
  })
  
	iPOT_Table_data <- reactive({
  	table <- merged_features_table[,c(1:20,50)]
  	table$log_gene_length<- round(2**as.numeric(table$log_gene_length),0)
		index_without_pos <- which(is.na(table$start))
		table <- mutate(table, genomelocation = paste0(chr,":", start,"-",end," (",chain,")"))
		table$genomelocation[index_without_pos]<-NA
		prediction <- rep("",length(table))
		table$TSG_probability<-round(as.numeric(table$TSG_probability),4)
		table$OG_probability<-round(as.numeric(table$OG_probability),4)
		pval_driver<-apply(cbind(table$TSG_probability,table$OG_probability),1,max)
		prediction[table$TSG_probability>TSG_threshold & table$OG_probability<OG_threshold]<-"Predicted TSG"
		prediction[table$OG_probability>OG_threshold & table$TSG_probability<TSG_threshold]<-"Predicted OG"
		prediction[table$TSG_probability>TSG_threshold & table$OG_probability>OG_threshold]<-"Predicted TSG/OG"

		merged_features_table$CGC_status <- as.character(merged_features_table$CGC_status)
		CGC_index <- which(!is.na(merged_features_table$CGC_status))

		#prediction[CGC_index]<-merged_features_table$CGC_status[CGC_index]
		table$prediction <- prediction
		table$prediction <- as.factor(table$prediction)
		quantile_table<-quantile(table$TSG_probability[table$TSG_probability>TSG_threshold], prob = seq(0, 1, length = 6), type = 7,na.rm=T)
		Tier_TSG <- ""
		Tier_TSG[table$TSG_probability>=TSG_threshold]<-"5"
		Tier_TSG[table$TSG_probability>=quantile_table[2]]<-"4"
		Tier_TSG[table$TSG_probability>=quantile_table[3]]<-"3"
		Tier_TSG[table$TSG_probability>=quantile_table[4]]<-"2"
		Tier_TSG[table$TSG_probability>=quantile_table[5]]<-"1"
		quantile_table<-quantile(table$OG_probability[table$OG_probability>OG_threshold], prob = seq(0, 1, length = 6), type = 7,na.rm=T)
		Tier_OG <- ""
		Tier_OG[table$OG_probability>OG_threshold]<-"5"
		Tier_OG[table$OG_probability>quantile_table[2]]<-"4"
		Tier_OG[table$OG_probability>quantile_table[3]]<-"3"
		Tier_OG[table$OG_probability>quantile_table[4]]<-"2"
		Tier_OG[table$OG_probability>quantile_table[5]]<-"1"
		
		table <- cbind(table,pval_driver,Tier_TSG,Tier_OG)
		table <- table[,c(1,8,22,6,23,24,14,15,25,26,7,9,10,11,12,13,21)]
		table<-table[order(-round(as.numeric(table$pval_driver),2)),]
  })
  

  observeEvent(input$iPOT_Table_rows_selected, {
  		table <- iPOT_Table_data()
      geneclicked <- as.character(table[input$iPOT_Table_rows_selected,1])
      updateTabsetPanel(session, 'maintabs', selected = "by Gene")
      updateSelectizeInput(session, "gene_input", selected = geneclicked)
  })
  
  output$iPOT_Table <- DT::renderDataTable({
      
			table <- iPOT_Table_data()
      DT::datatable(table,
      callback = JS("
				var tips = ['Gene','Gene full name','Genome location (hg19)','CGC annotation','Prediction status','Driver score','DORGE-TSG score','DORGE-OG','Tier TSG (1: score within 80-100 percentile; 2: score within 60-80 percentile; 3: score within 40-60 percentile; 4: score within 20-40 percentile; 5:score within 0-20 percentile of predicted TSGs)','Tier OG (1: score within 80-100 percentile; 2: score within 60-80 percentile; 3: score within 40-60 percentile; 4: score within 20-40 percentile; 5:score within 0-20 percentile of predicted OGs)','Alias','Gene type','Gene loci','Ensembl gene ID','RefSeq ID','Pubmed ID','CDS length'];
				header = table.columns().header();
				for (var i = 0; i < tips.length; i++) {
				  $(header[i]).attr('title', tips[i]);
				}
			"),
      selection = 'single',rownames = FALSE,colnames=c("Gene","Gene full name","Genome location",'CGC annotation',"Prediction status","Driver score","DORGE-TSG score","DORGE-OG score","Tier TSG","Tier OG","Alias","Gene type","Gene loci","Ensembl gene","RefSeq","Pubmed","CDS length"),filter = 'top',extensions = 'Buttons',
      options = list(
      								search = list(regex = TRUE, caseInsensitive = TRUE),
							      	dom = 'Bfrtip',buttons = list(c(I('colvis'),'copy', 'csv', 'excel')),
							      	pageLength = 20, lengthMenu = c(10, 20, 30), stateSave = TRUE
							      )
      ) %>% formatStyle('prediction',backgroundColor = styleEqual(c("Predicted OG","Predicted TSG","Predicted TSG/OG"), c('#4DAF4A','#984EA3','transparent'))) %>% formatStyle('CGC_status',backgroundColor = styleEqual(c("CGC OG","CGC TSG","CGC TSG/OG"), c('#E41A1C', '#377EB8', 'grey')))
      
  })
  
  output$tissue_exp_table <- DT::renderDataTable({
      table <- tissue_gene_exp_Table()
			table<-join(table, pvalues_TSG,type = "left",by="Gene")
			table<-join(table, pvalues_OG,type = "left",by="Gene")
			table<-table[,c(1,2,3,4)]
			table$Z_SCORE<-round(as.numeric(table$Z_SCORE),4)
			table$TSG_probability<-round(as.numeric(table$TSG_probability),4)
			table$OG_probability<-round(as.numeric(table$OG_probability),4)
      DT::datatable(table,selection = 'single',rownames = FALSE,colnames=c('Gene','Z SCORE',"TSG-score","OG-score"),filter = 'top',
      extensions = 'Buttons',options = list(dom = 'Bfrtip', buttons = c(I('colvis'),'copy', 'csv', 'excel'),pageLength = 20, lengthMenu = c(10, 20, 30), stateSave = TRUE))
  })
  
  output$cellline_exp_table <- DT::renderDataTable({
      table <- cellline_gene_exp_Table()
      table<-join(table, pvalues_TSG,type = "left",by="Gene")
			table<-join(table, pvalues_OG,type = "left",by="Gene")
			table<-table[,c(1,3,2,4,5,6)]
			table$TSG_probability<-round(as.numeric(table$TSG_probability),4)
			table$OG_probability<-round(as.numeric(table$OG_probability),4)
      DT::datatable(table,selection = 'single',rownames = FALSE,colnames=c('Gene','Expression TPM','Cell line',"Cell line illustration","TSG-score","OG-score"),filter = 'top',
      extensions = 'Buttons',options = list(dom = 'Bfrtip', buttons = c(I('colvis'),'copy', 'csv', 'excel'),pageLength = 20, lengthMenu = c(10, 20, 30), stateSave = TRUE))
  })
  
  output$tissue_exp_barchart <- renderPlotly({
  
    table <- tissue_gene_exp_Table()

    if(nrow(table)>0){
	    table$Z_SCORE<-as.numeric(as.character(table$Z_SCORE))
	    table$Gene<-as.character(table$Gene)
	   	table1 = aggregate(Z_SCORE ~ Gene,table,FUN = mean, simplify = TRUE, exclude = "", na.rm=TRUE)
	    table1 <- table1[order(-table1$Z_SCORE),]
	    genekept<-as.character(table1$Gene)[1:30]
	    
	    table1 <- table1[which(table1$Gene %in% genekept),]
	    table1 <- table1[order(-table1$Z_SCORE),]
	    levels <- as.character(table1$Gene)
	    table1$Gene <- factor(table1$Gene, levels = levels)
	    table <- table1
	    table$Z_SCORE <- as.numeric(table$Z_SCORE)
	    	    
	    p <- plot_ly(table, x=~Gene, y=~Z_SCORE, source='tissue_exp_barchart', type = 'bar', name = 'noname', marker=list(color="darkblue")) %>% layout(yaxis = list(title = 'Expression in TCGA Tumor tissues'), barmode = 'stack', margin = list(b = 100), xaxis=list(title = "", tickangle = 45))%>% config(displayModeBar = F)

    } else {
      p <- plotly_empty(type='pie') %>% config(displayModeBar = F)
    }
    
    p$elementId <- NULL
	  p
  })
  
  output$cellline_exp_barchart <- renderPlotly({
  
    table <- cellline_gene_exp_Table()
    if(nrow(table)>0){
    	
  		table$Exp<-as.numeric(as.character(table$Exp))
	    table$Gene<-as.character(table$Gene)
	   	table1 = aggregate(Exp ~ Gene,table,FUN = mean, simplify = TRUE, exclude = "", na.rm=TRUE)
	    table1 <- table1[order(-table1$Exp),]
	    genekept<-as.character(table1$Gene)[1:30]
	    
	    table1 <- table1[which(table1$Gene %in% genekept),]
	    table1 <- table1[order(-table1$Exp),]
	    levels <- as.character(table1$Gene)
	    table1$Gene <- factor(table1$Gene, levels = levels)
	    table <- table1
	    table$Exp <- as.numeric(table$Exp)
	    
	    p <- plot_ly(table, x=~Gene, y=~Exp, source='cellline_exp_barchart', type = 'bar', name = 'noname', marker=list(color="darkblue")) %>% layout(yaxis = list(title = 'Expression in CCLE cell lines'), barmode = 'stack', margin = list(b = 100), xaxis=list(title = "", tickangle = 45))%>% config(displayModeBar = F)
	   
    } else {
      p <- plotly_empty(type='pie') %>% config(displayModeBar = F)
    }
    
    p$elementId <- NULL
	  p
  })
  
  output$gene_table <- DT::renderDataTable({
      table <- geneData_H3K4me3()
      if(geneNames %in% table$Gene){
      	DT::datatable(table[,c('Gene','Peak_length','Sample_name','Cell_line')],selection = 'single',rownames = FALSE,colnames=c('Gene','Peak length','Sample name','Cell line'),options = list(pageLength = 20, lengthMenu = c(10, 20, 30), stateSave = TRUE))
      }
      
  })
  
  output$gene_table2 <- DT::renderDataTable({
      table <- geneData_canyon()
      if(geneNames %in% table$Gene){
      	DT::datatable(table[,c('Gene','Sample_name','WGBS_methylation','Tissue','Sample_ID')],selection = 'single',rownames = FALSE,colnames=c('Gene','Sample name','WGBS methylation (0-1)','Tissue','Sample ID'),options = list(pageLength = 20, lengthMenu = c(10, 20, 30), stateSave = TRUE))
      }
  })

  output$gene_overview <- renderPlotly({

    specific_gene_info <- geneinfo()
    
    if (nrow(specific_gene_info) > 0) {

	    role<-c("DORGE-TSG","DORGE-OG")
	    color<-c("darkblue","red")
	    role_pval <- data.frame(role=role, role_pval=c(round(as.numeric(specific_gene_info$TSG_probability),3),round(as.numeric(specific_gene_info$OG_probability),3)), color)

	    p <- plot_ly(role_pval, x=~role, y=~role_pval, source='gene_overview', type = 'bar',text = ~paste("TSG score: ", round(specific_gene_info$TSG_probability,3), '<br>OG score:', round(specific_gene_info$OG_probability,3)), marker=list(color=role_pval$color)) %>%
	      layout(yaxis = list(title = 'DORGE-TSG/OG score'), margin = list(b = 100), xaxis=list(title = "", tickangle = 80))%>% 
	      config(displayModeBar = F)
    } else {
      p <- plotly_empty(type='pie')%>% 
        config(displayModeBar = F)
    }
    p$elementId <- NULL
    p
  })
  
  gene_bar_event_val <- reactiveValues(count = 0)
  output$gene_barchart <- renderPlotly({
  
    table <- geneData_H3K4me3()
    
    if(nrow(table)>0){
    	if(nrow(table)>30){
	    	table <- table[order(-table$Peak_length),]
	    	table <- table[1:30,]
	    }
  		table$Peak_length<-as.numeric(as.character(table$Peak_length))
	    table$Gene<-as.character(table$Gene)
	   	table1 = aggregate(Peak_length ~ Cell_line,table,FUN = mean, simplify = TRUE, exclude = "", na.rm=TRUE)
	    table1 <- table1[order(-table1$Peak_length),]
	    levels <- as.character(table1$Cell_line)
	    
	    celllinekept<-as.character(table1$Cell_line)
	    if(nrow(table)>30){
	    	celllinekept<-as.character(table1$Cell_line)[1:30]
	    }
	    table1 <- table1[which(table1$Cell_line %in% celllinekept),]
	    table1 <- table1[order(-table1$Peak_length),]
	    table1$Cell_line <- factor(table1$Cell_line, levels = levels)
	    table1$Peak_length <- round(table1$Peak_length,0)
	    table <- table1
    
	    sourceName <- paste('gene_barchart_',gene_bar_event_val$count,sep='')
	    p <- plot_ly(table, x=~Cell_line, y=~Peak_length, source='gene_barchart', type = 'bar', name = 'noname', marker=list(color="darkblue")) %>% layout(yaxis = list(title = 'H3K4me3 peak length'), barmode = 'stack', margin = list(b = 100), xaxis=list(title = "", tickangle = 45))%>% config(displayModeBar = F)
	    
    } else {
      p <- plotly_empty(type='pie') %>% config(displayModeBar = F)
    }
    
    p$elementId <- NULL
	  p
  })
  
  output$gene_barchart2 <- renderPlotly({
  
  	table <- geneData_canyon()
     
    if (nrow(table) > 0) {
    	table$WGBS_methylation<-as.numeric(as.character(table$WGBS_methylation))
	    table$Gene<-as.character(table$Gene)
	   	table1 = aggregate(WGBS_methylation ~ Tissue,table,FUN = mean, simplify = TRUE, exclude = "", na.rm=TRUE)
	    table1 <- table1[order(-table1$WGBS_methylation),]
	    levels <- as.character(table1$Tissue)
	    
	    Tissuekept<-as.character(table1$Tissue)
	    if(nrow(table1)>30){
	    	Tissuekept<-as.character(table1$Tissue)[1:30]
	    }

	    table1 <- table1[which(table1$Tissue %in% Tissuekept),]
	    table1 <- table1[order(-table1$WGBS_methylation),]
	    table1$Tissue <- factor(table1$Tissue, levels = levels)
	    table <- table1
	    sourceName <- paste('gene_barchart_',gene_bar_event_val$count,sep='')
	    p <- plot_ly(table, x=~Tissue, y=~WGBS_methylation, source='gene_barchart2', type = 'bar', name = 'noname', marker=list(color="red")) %>% layout(yaxis = list(title = 'Methylation at Canyon'), barmode = 'stack', margin = list(b = 150), xaxis=list(title = "", tickangle = 45))%>% config(displayModeBar = F)
	    
    } else {
      p <- plotly_empty(type='pie') %>% config(displayModeBar = F)
    }
    
    p$elementId <- NULL
	  p
  })
  
  output$tissue_exp_download_collated_all <- downloadHandler(
    filename = function() {
      return("TCGA_gene_expression_info_for_selected_tissue.txt")
    },
    content = function(file) {
      table <- tissue_gene_exp_Table()
      write.table(table, file, row.names = FALSE, sep='\t', quote=F)
    }
  )
  
  output$cellline_exp_download_collated_all <- downloadHandler(
    filename = function() {
      return("CCLE_gene_expression_info_selected_cellline.txt")
    },
    content = function(file) {
      table <- cellline_gene_exp_Table()
      write.table(table, file, row.names = FALSE, sep='\t', quote=F)
    }
  )
  
  output$gene_download_collated_all <- downloadHandler(
    filename = function() {
      return("H3K4me3_peak_length_info_unaggregated.txt")
    },
    content = function(file) {
      table <- geneData_H3K4me3()
      write.table(table, file, row.names = FALSE, sep='\t', quote=F)
    }
  )
  
  output$gene_download_collated_shown <- downloadHandler(
    filename = function() {
      return("H3K4me3_mean_peak_length_info_for_cell_lines.txt")
    },
    content = function(file) {
	    table <- geneData_H3K4me3()
	    if(nrow(table)>0){
	    	
	  		table$Peak_length<-as.numeric(as.character(table$Peak_length))
		    table$Gene<-as.character(table$Gene)
		   	table1 = aggregate(Peak_length ~ Cell_line,table,FUN = mean, simplify = TRUE, exclude = "", na.rm=TRUE)
		    table1 <- table1[order(-table1$Peak_length),]
		    levels <- as.character(table1$Cell_line)
		    
		    celllinekept<-as.character(table1$Cell_line)

		    table1 <- table1[which(table1$Cell_line %in% celllinekept),]
		    table1 <- table1[order(-table1$Peak_length),]
		    table1$Cell_line <- factor(table1$Cell_line, levels = levels)
		    table1$Peak_length <- round(table1$Peak_length,0)
		    table <- table1
		  }
      write.table(table, file, row.names = FALSE, sep='\t', quote=F)
    }
  )
  
  output$gene_download_collated_all2 <- downloadHandler(
    filename = function() {
      return("Genebody_canyon_methylation_info_unaggregated.txt")
    },
    content = function(file) {
      table <- geneData_canyon()
      write.table(table, file, row.names = FALSE, sep='\t', quote=F)
    }
  )
  
  output$gene_download_collated_shown2 <- downloadHandler(
    filename = function() {
      return("Genebody_canyon_methylation_info_for_cell_lines.txt")
    },
    content = function(file) {
  		table <- geneData_canyon()
	    if (nrow(table) > 0) {
	    	table$WGBS_methylation<-as.numeric(as.character(table$WGBS_methylation))
		    table$Gene<-as.character(table$Gene)
		   	table1 = aggregate(WGBS_methylation ~ Tissue,table,FUN = mean, simplify = TRUE, exclude = "", na.rm=TRUE)
		    table1 <- table1[order(-table1$WGBS_methylation),]
		    levels <- as.character(table1$Tissue)
		    
		    Tissuekept<-as.character(table1$Tissue)

		    table1 <- table1[which(table1$Tissue %in% Tissuekept),]
		    table1 <- table1[order(-table1$WGBS_methylation),]
		    table1$Tissue <- factor(table1$Tissue, levels = levels)
		    table <- table1
		    
	    }
      write.table(table, file, row.names = FALSE, sep='\t', quote=F)
    }
  )
  
}

#shinyApp(ui, server)
shinyApp(ui=ui,server=server, options = list(port=7990))
#runApp("/Users/jlyu/Box Sync/TSGOG_Project/SA_sub/github/DORGE_paper/DORGE_shiny",launch.browser=T)