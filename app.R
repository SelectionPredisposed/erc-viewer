#!/usr/bin/Rscript

# Sets deploy flag if on the server node, small hack to have same repo for
# local testing and deployment
if(Sys.info()[["nodename"]]=='eserver'){ 
  deploy <- T
  setwd('/srv/shiny-server/ercviewer/')
}else{
  deploy <- F
  setwd('.')
}

#### ----- everything here is run only once ----- ####

library(shiny)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(data.table)
library(readxl)
library(httr)
library(jsonlite)
library(xml2)
library(biomaRt)
#library(Cairo) # Cairo graphics improves rendering on some Linux systems
options(stringsAsFactors=FALSE) #, shiny.usecairo=TRUE)


# file with gene ids and RefSeq positions:
allGenes = read.table("./resources/genelist.txt", h=F)
colnames(allGenes) = c("Gene", "Chr", "Start", "Stop")
allGenes = filter(allGenes, Chr %in% 1:23)
allGenes$Chr = as.numeric(allGenes$Chr)
allGenes$ypos = rep(1:5, length.out=nrow(allGenes))

# file with full SNP positions
if(!file.exists("chrTable")){
	d = fread("bmi3_additive", h=T, stringsAsFactors = F, data.table = F)
	d <- d[,c('chromosome','position','all_maf','frequentist_add_pvalue')]
	colnames(d) = c("Chr", "Pos", "MAF", "P")
	chrTable = group_by(d, Chr) %>%
	summarize(ChrL = as.numeric(max(Pos))) %>%
	 	arrange(Chr) %>%
		mutate(ChrStart = cumsum(ChrL)-ChrL)
	write.table(chrTable, "chrTable", quote=F, row.names=F, col.names=T)
} else {
	chrTable = read.table("chrTable", h=T)
}
allGenes = inner_join(allGenes, chrTable, by="Chr") %>%
	mutate(FullStart = Start + ChrStart, FullStop = Stop + ChrStart)

allRecomb = read.table(gzfile('./resources/recombination_rates.txt.gz'), header = T)
allRecomb.tmp = inner_join(allRecomb, chrTable, by="Chr") %>%
  mutate(AbsPos = Position.bp. + ChrStart)

# Load previous hits
# load markers extracted from zbmi3 rot1
b37pos <- fread('zcat rsid_chr_pos_b37.gz', h=T, data.table = F)

felix.hits = read.table("resources/felix_bmi43_snps_full_proxy_table.csv", h=T, sep='\t') %>% subset(USE==1)
felix <- subset(b37pos, rsid %in% felix.hits$SNP_PROXY) %>% dplyr::select(chromosome, position) %>% rename(Chr=chromosome, Pos=position)
 
locke.hits = read.table("resources/locke_bmi97_snps_full_proxy_table.csv", h=T, sep='\t') %>% subset(USE==1)
locke <- subset(b37pos, rsid %in% locke.hits$SNP_PROXY) %>% dplyr::select(chromosome, position) %>% rename(Chr=chromosome, Pos=position)
 
horikoshi.hits = read.table("resources/horikoshi_bw58_snps_full_proxy_table.csv", h=T, sep='\t') %>% subset(USE==1)
horikoshi <- subset(b37pos, rsid %in% horikoshi.hits$SNP_PROXY) %>% dplyr::select(chromosome, position) %>% rename(Chr=chromosome, Pos=position)
 
transt2d.hits = read.table("resources/t2d-transethnic-diagram-scott2017-t2d76_snps_full_proxy_table.csv", h=T, sep='\t') %>% subset(USE==1)
transt2d <- subset(b37pos, rsid %in% transt2d.hits$SNP_PROXY) %>% dplyr::select(chromosome, position) %>% rename(Chr=chromosome, Pos=position)

### -- input data files are loaded here -- ###
### if there's a lot of files and need to save RAM, can move this to 
### server.ui so that the user would wait for each file to load

# files must be named <geno>_<model> and be present in the ./data subfolder
# must have the following columns, with any names in the header:
# "chromosome position all_maf frequentist_add_pvalue"
if(deploy){
  batch = batch = c('HARVESTv9', 'ROTTERDAM1','META','HARVESTv10')
  genos = c("bmi0"="bmi0",
            "bmi1"="bmi1",
            "bmi2"="bmi2",
            "bmi3"="bmi3",
            "bmi4"="bmi4",
            "bmi5"="bmi5",
            "bmi6"="bmi6",
            "bmi8"="bmi8",
            "bmi9"="bmi9",
            "bmi10"="bmi10",
            "bmi11"="bmi11",
            "bmi12"="bmi12")
  models = c("additive"="additive")
} else {
  genos = c("bmi0"="bmi0", "bmi1"="bmi1")
  models = c("additive"="additive")
  batch = c('META','HARVEST', 'ROTTERDAM1')
}

#alldata = list()
# for(g in genos){
# 	for(m in models){
# 		model = paste(g, m, sep="_")
# 		infile = file.path("data", model)
# 		print(sprintf("looking for file %s", infile))
# 		if(file.exists(infile)) {
# 			f = fread(infile, h=T, data.table=F)
# 			      f <- f %>% select(chromosome,position,all_maf,frequentist_add_pvalue,rsid)
#         		colnames(f) = c("Chr", "Pos", "MAF", "P", "rsid")
# 			f$Chr = as.integer(f$Chr)
# 			alldata[[model]] = f
# 			print(sprintf("file %s is read", infile))
# 		}
# 	}
# }
gc()
### -- end data read -- ##

# JS functions for making log and scientific sliders
JS.logify = "function logifySlider (sliderId, sci = false) {
  if (sci) {
    // scientific style
    $('#'+sliderId).data('ionRangeSlider').update({
        'prettify': function (num) { return ('10<sup>'+num+'</sup>'); }
    })
  } else {
    // regular number style
    $('#'+sliderId).data('ionRangeSlider').update({
        'prettify': function (num) {
            return (parseFloat(Math.pow(5, num)).toFixed(4));
        }
    })
}}"

# call logifySlider for each sliderInput
JS.onload = "
// execute upon document loading
$(document).ready(function() {
  // wait a few ms to allow other scripts to execute
  setTimeout(function() {
    // include calls for each slider here:
    logifySlider('maf_cutoff', sci = false)
    logifySlider('p_cutoff', sci = true)
  }, 100)
})"
# functions for manhattan breaks and labels
breakmaker = function(x){
	seq(min(x), max(x), length.out = 10)
}
labelmaker = function(b){
       	breakchrs = findInterval(b, chrTable$ChrStart)
       	paste(chrTable$Chr[breakchrs], round(b - chrTable$ChrStart[breakchrs]), sep=":")
}

# function for creating link in biomart results
createLink <- function(val) {
  sprintf('<a href="https://www.ncbi.nlm.nih.gov/pubmed/%s" target="_blank" class="btn btn-primary">Article</a>',val)
}

#### ----- everything here is individually launched for each user ----- ####
ui <- fluidPage(
    tags$head(tags$script(HTML(JS.logify))),
    tags$head(tags$script(HTML(JS.onload))),

    titlePanel("GWAS Viewer"),

    # settings panel for the left side
    sidebarLayout(
      sidebarPanel(
          selectInput("select_batch",
                    label="Select analysis batch:",
                    choices=batch),
          selectInput("select_genomes",
                      label="Select the source of genomes to analyze:",
                      choices=genos),
          selectInput("select_model",
                      label="Select the type of model to analyze:",
                      choices=models),
	  # see scientific slider js for interpreation of values
          sliderInput("maf_cutoff",
                      "Select lower MAF cutoff:",
                      min = -5.7, max = -1.4307, value=-2.86, step=0.05, width="100%"),
          sliderInput("p_cutoff",
                      "Select the maximum -log10(P) cutoff to include SNPs in the genome viewer:",
                      min = -5, max = 0, value=-2.5, step=0.1, width="100%"),
          sliderInput("ld_window",
                      "Select the distance (in Kb) around selected marker to retrieve LD for:",
                      min=100, max=500, value=500, step=10, width="100%"),
          actionButton("load_file", "Load data"),
          checkboxInput("add_labels", "Add SNP IDs to the plot (when SNPs<50)",value = TRUE),
          checkboxInput("mark_locke", "Highlight SNPs from Locke (adult BMI) - orange"),
          checkboxInput("mark_felix", "Highlight SNPs from Felix (childhood BMI ~8-10y) - red"),
          checkboxInput("mark_horikoshi", "Highlight SNPs from Horikoshi (birth weight) - brown"),
          checkboxInput("mark_transt2d", "Highlight SNPs from Scott (transethnic T2D) - black"),
          checkboxInput("color_maf", "Add MAF color scale (not compatible with the highlighting above)")
    ),

    mainPanel(
	# status outputs
        fluidRow(
            column(8,
              textOutput("readout_maf"),
              textOutput("readout_p"),
              textOutput("readout_lzval")),
            column(4, wellPanel(textOutput("status")))
            
        ),

	# main plot
        plotOutput("manhattan", dblclick = "manh_click", click = "manh_single_click",
                  brush = brushOpts(direction="x", id="manh_brush", resetOnNew = T)),
          
	# plot adjustments
	fluidRow(
	      actionButton("make_snap", "Save current view as snapshot")
	),
	
	# additional plots/tables
	tabsetPanel(
	    tabPanel("Locus Zoom", 
	             plotOutput("locuszoom", click = "lz_single_click"), 
	             plotOutput("recombplot",height = 200),
	             plotOutput("geneplot2")),
	    tabPanel("Gene Context", plotOutput("geneplot")),
	    tabPanel("Snapshot", plotOutput("snapshot")),
	    tabPanel("Table", tableOutput("rawtable")),
	    tabPanel("biomaRt", tableOutput("biomarttable"))
	)
      )
   )
)

server <- function(input, output) {
    plot_maf = plot_p = plot_geno = plot_model = 0
    manh <<- NULL
    data = reactiveValues(d = NULL)
    ranges = reactiveValues(x = NULL)
    lzval <- reactiveValues(selmarker = NULL, lduse = NULL)
    
    output$readout_lzval <- reactive({
      sprintf("LD markers:", lzval$selmarker)
    })
    

    output$readout_maf <- reactive({
        sprintf("Showing only SNPs with MAF above %5.4f", 5^input$maf_cutoff)
    })
    output$readout_p <- reactive({
        sprintf("Showing only SNPs with P below %.3e", 10^input$p_cutoff)
    })
    
    # reactive check to see if any of the settings were changed since last plotting
    output$status <- renderText({
        if(input$maf_cutoff!=plot_maf | input$p_cutoff!=plot_p | input$select_genomes!=plot_geno | input$select_model!=plot_model) "plot is NOT current."
          else "plot is current."
    })
    
    observeEvent(input$load_file, {
	  # to check if plot is current:
	  # (should make reactives instead of globals here)
      plot_maf <<- input$maf_cutoff
      plot_p <<- input$p_cutoff
      plot_batch <<- input$select_batch
      plot_geno <<- input$select_genomes
      plot_model <<- input$select_model
      maf_cut = 5^input$maf_cutoff
      p_cut = 10^input$p_cutoff

    	# select the data source, filter MAF and P.
    	# could be moved out from observeEvent,
    	# so that the plot would update immediately
    	withProgress({
    	       
    	        # data loading happens on demand instead of loading all files upfront
    	        alldata = list()
              model = paste(plot_geno, plot_model, sep="_")
          	  infile = file.path("data",plot_batch, model)
          	  print(sprintf("looking for file %s", infile))
      	      if(file.exists(infile)) {
      	        f = fread(infile, h=T, data.table=F)
      	        f <- f %>% dplyr::select(chromosome,position,all_maf,frequentist_add_pvalue,rsid)
      	        colnames(f) = c("Chr", "Pos", "MAF", "P", "rsid")
      	        f$Chr = as.integer(f$Chr)
      	        alldata[[model]] = f
      	        print(sprintf("file %s is read", infile))
      	      }
    	  
           		data$d = alldata[[paste(plot_geno, plot_model, sep="_")]] 
           		print(data)
            	data$d = inner_join(data$d, chrTable, by="Chr") %>%
            	    mutate(FullPos = Pos + ChrStart, logP=-log(P,10))
    		data$d = filter(data$d, MAF>maf_cut, P<p_cut)
    	}, message="filtering data")
            dfull = data$d
    	# might be good to clean up:
    	gc()
        
	    # generate main plot
      output$manhattan = renderPlot({ withProgress({
		  # pre-filtering is faster than calling ggplot with all points
        	if(!is.null(ranges$x)){
        		data$d = filter(data$d, FullPos>ranges$x[1], FullPos<ranges$x[2])
        	} else {
        		data$d = dfull
        	}

		      # main plot call
        	manh <<- ggplot(data$d, aes(x=FullPos, y=logP)) +
			      geom_point() +
        		geom_hline(yintercept=5, col="orange") + geom_hline(yintercept=7.3, col="darkred") +
        		coord_cartesian(xlim = ranges$x) +
        		scale_x_continuous(breaks = chrTable$ChrStart+chrTable$ChrL/2,
        						   labels = chrTable$Chr) +
        		theme_bw()
        	
      		# change color scale to show MAF
      		if(input$color_maf){
      			manh <<- manh + aes(color = MAF) +
      				scale_color_continuous(trans="log10", low="lightyellow2", high="blue4", breaks=c(0.001, 0.005, 0.01, 0.02, 0.05, 0.1, 0.5), limits=c(0.001,0.5))
      			  #scale_color_gradientn(trans="log10", colours = rainbow(4), breaks=c(0.001, 0.005, 0.01, 0.02, 0.05, 0.1, 0.5), limits=c(0.001,0.5))
      		} else {
      			manh <<- manh + aes(color = Chr%%2==0) +
      				scale_color_manual(values=c("skyblue3", "lightgreen", "red"), guide=FALSE)
      		}
        	
		      # mark some predetermined candidate SNPs in plot
        	if(input$mark_locke){
        		locke_toplot = semi_join(data$d, locke, by=c("Chr", "Pos"))
        		if(nrow(locke_toplot)>0) {
        			manh <<- manh + geom_point(aes(col="X"), colour='orange', locke_toplot)
        		}
        	}
        	
        	if(input$mark_felix){
        	  felix_toplot = semi_join(data$d, felix, by=c("Chr", "Pos"))
        	  if(nrow(felix_toplot)>0) {
        	    manh <<- manh + geom_point(aes(col="X"),colour='red', felix_toplot)
        	  }
        	}
        	
        	if(input$mark_horikoshi){
        	  horikoshi_toplot = semi_join(data$d, horikoshi, by=c("Chr", "Pos"))
        	  if(nrow(horikoshi_toplot)>0) {
        	    manh <<- manh + geom_point(aes(col="X"),colour='brown', horikoshi_toplot)
        	  }
        	}
        	
        	if(input$mark_transt2d){
        	  transt2d_toplot = semi_join(data$d, transt2d, by=c("Chr", "Pos"))
        	  if(nrow(transt2d_toplot)>0) {
        	    manh <<- manh + geom_point(aes(col="X"),colour='black', transt2d_toplot)
        	  }
        	}
        	
		      # switch between chromosome/position labels depending on zoom
        	if(!is.null(ranges$x)){ 
        		manh <<- manh + scale_x_continuous(breaks = breakmaker,
        						   labels = labelmaker)
        		if(input$add_labels & nrow(data$d)<50){
        			manh <<- manh + geom_text_repel(aes(label=rsid))
        		}
        	}
        	manh
        }, message="plotting") })

      output$locuszoom = renderPlot({
        lz <<- ggplot(data$d, aes(x=FullPos, y=logP)) +
          geom_point() +
          geom_hline(yintercept=5, col="orange") + geom_hline(yintercept=7.3, col="darkred") +
          coord_cartesian(xlim = ranges$x) +
          scale_x_continuous(breaks = chrTable$ChrStart+chrTable$ChrL/2,
                             labels = chrTable$Chr) +
          theme_bw()
          
          if(!is.null(lzval$selmarker) & !is.null(lzval$lduse)){
            if(nrow(lzval$lduse) > 0){
              m <- merge(data$d, lzval$lduse, by.x='rsid', by.y='variation2', all.x=T)
              lz <<- lz + geom_point(data=m, aes(FullPos, -log10(P), colour=r2)) + 
                scale_colour_gradient(low = "red", high = "green", limits=c(0,1), breaks=seq(0,1,0.2))
              lz <<- lz + geom_point(data=subset(m, rsid==lzval$selmarker), aes(FullPos, -log10(P)),col='blue', size=4) +
                theme(legend.position = c(0.95,0.8))
            }
          }
        
          if(input$add_labels & nrow(data$d)<50){
            lz <<- lz + geom_text_repel(aes(label=rsid))
          }
          lz
        })
        
	      # make the gene-context plot
        output$geneplot = renderPlot({
            ggplot(allGenes) +
                geom_segment(aes(x=FullStart, xend=FullStop, y=ypos, yend=ypos), lwd=3, col="darkblue") +
                geom_text(aes(label=Gene, x=(FullStart+FullStop)/2, y=ypos+0.2), size=4, fontface="italic") +
                coord_cartesian(xlim = ranges$x) +
        		scale_x_continuous(breaks=breakmaker, labels=labelmaker) +
                theme_bw()
        })
        
        # Makge gene-context plot a second time due to it be using two different places and need unique ID
        output$geneplot2 = renderPlot({
          ggplot(allGenes) +
            geom_segment(aes(x=FullStart, xend=FullStop, y=ypos, yend=ypos), lwd=3, col="darkblue") +
            geom_text(aes(label=Gene, x=(FullStart+FullStop)/2, y=ypos+0.2), size=4, fontface="italic") +
            coord_cartesian(xlim = ranges$x) +
            scale_x_continuous(breaks=breakmaker, labels=labelmaker) +
            theme_bw()  
        })
        # 
        output$recombplot = renderPlot({
          if(!is.null(ranges$x)){
            data$recomb = filter(allRecomb.tmp, AbsPos>ranges$x[1], AbsPos<ranges$x[2])
            ggplot(data$recomb) +
              geom_line(aes(x=AbsPos, y=Rate.cM.Mb.)) +
              coord_cartesian(xlim = ranges$x) +
              theme_bw()
          } else {
            ggplot() + 
              annotate("text", x = 1, y = 1, label = "Recombination rates will be available on zoom") + 
              theme_bw()
          }
        })
        

	# check if any settings changed since last plot call
        output$status <- renderText({
            if(input$maf_cutoff!=plot_maf | input$p_cutoff!=plot_p) "plot is NOT current."
                else "plot is current."
        })
    }, ignoreNULL = TRUE, ignoreInit = TRUE)

    # When a double-click happens, check if there's a brush on the plot.
    # If so, zoom to the brush bounds; if not, reset the zoom.
    observeEvent(input$manh_click, {
        brush <- input$manh_brush
        if (!is.null(brush)) {
            ranges$x <- c(brush$xmin, brush$xmax)
        } else {
            ranges$x <- NULL
        }
    })

    # on single click, find nearest point
    observeEvent(input$manh_single_click, {
      withProgress({
        df <- nearPoints(data$d, input$manh_single_click, addDist = TRUE)
        if(nrow(df) > 0){
          selmarker <- df[1,]$rsid
          
          ensembl = useMart(biomart = 'ENSEMBL_MART_SNP',
                            host="www.ensembl.org",
                            dataset = 'hsapiens_snp',
                            path="/biomart/martservice")
          
          res <- getBM(attributes = c('refsnp_id','associated_gene','allele','minor_allele_freq','phenotype_name','pmid','title','year'),
                       filters = "snp_filter",
                       values = list(selmarker),
                       mart = ensembl)
          output$biomarttable = renderTable({
            # create link out of pmid
            res$pmid <- createLink(res$pmid)
            
            # to avoid duplicates where associated genes differs
            # the following adds all unique gene names to all entries
            g <- c(unique(res$associated_gene))
            g <- paste(g,collapse = ",") 
            g <- strsplit(g,split = ',')             
            g <- unlist(g)
            g <- unique(g)
            g <- paste(g, collapse = ',')
            res$associated_gene <- g
            
            # remove duplicates
            res <- res[!duplicated(res),]
            
            # order by year, newest first
            res <- arrange(res,-year)
            return(res)
          }, sanitize.text.function = function(x) x)
        }
      }, message="fetching GWAS-catalog info...")
    })
    
    # save current view for comparison
    observeEvent(input$make_snap, {
    	snap = manh
    	output$snapshot = renderPlot({ snap })
    })

    # print raw data file
    output$rawtable = renderTable({
	if(nrow(data$d)<100){
		out = data$d[,-c(6:7)]
		out$Chr = as.character(out$Chr)
		out$Pos = as.character(out$Pos)
		out
	}
	else { data.frame(warning="Table will appear when <100 SNPs are seen in the plot") }
    }, digits = -3)
    
  observeEvent(input$lz_single_click, {
    withProgress({
      df <- nearPoints(data$d, input$lz_single_click, addDist = TRUE)
      
      # To make sure a request is not sent if no marker has been selected
      if(nrow(df) > 0){
        selmarker <- df[1,]$rsid
        server <- "https://rest.ensembl.org"
        ext <- sprintf("/ld/human/%s/1000GENOMES:phase_3:CEU?window_size=%d",selmarker, input$ld_window)
        
        r <- GET(paste(server, ext, sep = ""), content_type("application/json"))
        ld <- data.frame(t(sapply(content(r),c)))
        
        # Make sure API responds with data
        if(!ncol(ld)==0){
          ld$r2 <- as.numeric(unlist(ld$r2))
          ld$variation2 <- unlist(ld$variation2)
          
          lzval$selmarker <- selmarker
          lzval$lduse <- ld %>% dplyr::select(variation2, r2)
        }
      }
    }, message="fetching LD-info...")
    
  })
}

### ----- Run the application ----- ###
shinyApp(ui = ui, server = server)
#runApp(list(ui = ui, server = server),port=52435, host="0.0.0.0", launch.browser = FALSE)
