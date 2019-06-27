library(shiny)
library(DT)
source("helpers.R")
#############Functions defined here
makePlotContainers <- function(pop_info,mh_info,all_sample_info, height=300, width="100%", ...) {
  ## Validate inputs
  #validateCssUnit(width)
  #validateCssUnit(height)
  #Render multiple plots given input$<id> and store them into different output$<id>
  sub_sample=all_sample_info$KenID[which(all_sample_info$population==pop_info)]  
  #
  ## Construct plotOutputs by converting renderPlot objects into html component through plotOutput()
  lst <- lapply(sub_sample, function(i)
    plotOutput(sprintf('%s_%s', mh_info, i), height=height, width=width))
  
  ## Make columns
  #lst <- lapply(split(lst, (seq.int(n)-1)%/%ncol), function(x) column(12/ncol, x))
  do.call(tagList, lst) #convert to tagList as html component for renderUI()
}

renderPlots <- function(pop_info, mh_info,read_info,amp_info,all_sample_info,all_mh_data,output) {
  #Render multiple plots given input$<id> and store them into different output$<id>
  sub_sample=all_sample_info$KenID[which(all_sample_info$population==pop_info)]
  mh_region=all_mh_data[which(all_mh_data$Amplicon==mh_info),]
  for (each_sample in sub_sample) {
    #use local to make sure loop is happening
    local({
      each_sample_t <- each_sample  # need to be evaluated here
      temp_sn=paste0(mh_info,"_",each_sample_t)
      #print(temp_sn)
      temp_mh_figure_data<-read.table(paste0("./www/data/pooled_vcf/mh_figure_data/",
                                             mh_info,"_",each_sample_t,
                                             "_plotdata.txt"),header=TRUE)
      temp_sample_total_read=read_info$Total[which(as.character(read_info$SampleID)==each_sample_t)]
      #print(temp_sample_total_read)
      #get amplicon count data for this sample
      amplicon_sample_idx=which(colnames(amp_info)==each_sample_t)
      amplicon_idx=which(amp_info$Amplicon==mh_info)
      temp_amplicon_total_read=amp_info[amplicon_idx,amplicon_sample_idx]
      ## These would be your plots instead
      output[[temp_sn]] <- renderPlot({
        #ggplot2 
        plot_microhaplotype_figure(converted_snp_for_plot = temp_mh_figure_data,
                                   mh_region = mh_region,
                                   pop_info = pop_info,
                                   each_sample_total_read = temp_sample_total_read,
                                   amplicon_total_read = temp_amplicon_total_read,
                                   each_sample_t)        
      })
    })
  }
}
convert_to_haplotype<-function(x){
  return(gsub(".$","",gsub("(.{1})", "\\1-", x)))
}
#############Pre-load the data 
#All microhaplotype region coordinate
microhaplotype_data=read.table("Bed_mMH_Miseq_08092018_P_remove_six_bad_amp.bed",header=FALSE)
colnames(microhaplotype_data)=c("CHROM","START","END","Amplicon","left_primer_start","left_primer_end","right_primer_start","right_primer_end")
all_sample_name=read.table("all_sample_name.txt",header=TRUE) #all sample name and population information 
#total_read=read.table(paste0("./www/data/pooled_vcf/qc/all_sample_total_read.txt"),header=TRUE)
#############Server function used in UI defined below
##Notes:
#All operations in the server function must be in reactive context
#every render* function is an implicit reactive variable
#reactive() function directly defines one reactive variable
#every input${id} and output$<id> is a reactive variable
#every *Output will produce an html component in the webpage and should be paired with render* variable 
function(input,output){
  #Get total read number for each sample
  total_read<-reactive({tr<-read.table(paste0("./www/data/pooled_vcf/qc/all_sample_total_read.txt"),header=TRUE)
                       return(tr)})
  #Get total amplicon count for each sample
  total_amplicon<-reactive({ta<-read.table(paste0("./www/data/pooled_vcf/qc/all_sample_uniformity_count.txt"),header=TRUE)
                            return(ta)})
  #########Draw mh plot of all samples with ggplot2 and save them into a list for further plotting
  
  #########Dynamic UI to show KenID belonging to all population
  output$selectSample<-renderUI({selectInput(inputId="sampleID",
                                             label="Sample ID",
                                             choices=all_sample_name$KenID[which(all_sample_name$population==input$population)],
                                             size=10,selectize = FALSE)})
  ##########Plot a single MH_sample vcf data and enable clicking to show variant public and QC information
  #Read the corresponding plotting data, the plot will change depending on input value
  figure_plot_data_single<-reactive({fd=read.table(paste0("./www/data/pooled_vcf/mh_figure_data/",
                                                   input$mh_region,"_",input$sampleID,
                                                   "_plotdata.txt"),header=TRUE)
  return(fd)
  })
  #ggplot2 and renderPlot() for the figure plot data of the single sample
  output$mh_figure_single<-renderPlot({
    #access sample coverage data
    each_sample_total_read=total_read()$Total[which(as.character(total_read()$SampleID)==input$sampleID)]
    #get amplicon count data for this sample
    amplicon_sample_idx=which(colnames(total_amplicon())==input$sampleID)
    amplicon_idx=which(total_amplicon()$Amplicon==input$mh_region)
    amplicon_total_read=total_amplicon()[amplicon_idx,amplicon_sample_idx]
    #population of this sample
    pop_info=input$population
    #get mh_region data of this specific mh_region
    mh_region=microhaplotype_data[which(microhaplotype_data$Amplicon==input$mh_region),]
    #ggplot2 
    plot_microhaplotype_figure(converted_snp_for_plot = figure_plot_data_single(),
                               mh_region = mh_region,
                               pop_info = pop_info,
                               each_sample_total_read = each_sample_total_read,
                               amplicon_total_read = amplicon_total_read,
                               input$sampleID)
  })
  #Click information 1: public database annotation
  output$info <- renderText({
    #based on input$plot_click$x and input$plot_click$y value, search and return details of the variant
    mh_region=microhaplotype_data[which(microhaplotype_data$Amplicon==input$mh_region),]#get mh_region data
    if (is.null(input$plot_click)){
      return()
    }
    if (input$plot_click$y<=2.2 && input$plot_click$y>=1.8 && input$plot_click$x>=mh_region$START && input$plot_click$x<=mh_region$END){
      query_haplotype="hap2"
      query=ceiling(input$plot_click$x)#get the coordinate
      temp_data=figure_plot_data_single()[figure_plot_data_single()$haplotype==query_haplotype,]#extract the data
      query_idx=which(temp_data$end==query)#get correct index
      if (length(query_idx)>0){
        query_result=temp_data[query_idx,]
        if (query_result$status!="Novel"){
          paste0("SNP ID:",query_result$snpid,"\n",
                 "REF:",query_result$REF," ALT:",query_result$ALT,"\n",
                 "Haplotype 2:",query_result$BaseCall,"\n",
                 "Position:",query_result$end,"\n",
                 "Variant Category:",query_result$status,"\n",
                 "Allele Frequency in 1000Genome\n",
                 "Eastern Asian:",query_result$EAS,"\n",
                 "American:",query_result$AMR,"\n",
                 "African:",query_result$AFR,"\n",
                 "European:",query_result$EUR,"\n",
                 "Southern Asian:",query_result$SAS)
        }else{
          paste0("SNP ID:",query_result$snpid,"\n",
                 "REF:",query_result$REF,"\n",
                 "Haplotype 2:",query_result$BaseCall,"\n",
                 "Position:",query_result$end,"\n",
                 "Variant Category:",query_result$status)
        }
      }else{
        paste0("Please Click on SNPs for details")
      }
    }else if (input$plot_click$y<=1.2 && input$plot_click$y>=0.8 && input$plot_click$x>=mh_region$START && input$plot_click$x<=mh_region$END){
      query_haplotype="hap1"
      query=ceiling(input$plot_click$x)
      temp_data=figure_plot_data_single()[figure_plot_data_single()$haplotype==query_haplotype,]
      query_idx=which(temp_data$end==query)
      if (length(query_idx)>0){
        query_result=temp_data[query_idx,]
        if (query_result$status!="Novel"){
          paste0("SNP ID:",query_result$snpid,"\n",
                 "REF:",query_result$REF," ALT:",query_result$ALT,"\n",
                 "Haplotype 1:",query_result$BaseCall,"\n",
                 "Position:",query_result$end,"\n",
                 "Variant Category:",query_result$status,"\n",
                 "Allele Frequency in 1000Genome\n",
                 "Eastern Asian:",query_result$EAS,"\n",
                 "American:",query_result$AMR,"\n",
                 "African:",query_result$AFR,"\n",
                 "European:",query_result$EUR,"\n",
                 "Southern Asian:",query_result$SAS)
        }else{
          paste0("SNP ID:",query_result$snpid,"\n",
                 "REF:",query_result$REF,"\n",
                 "Haplotype 1:",query_result$BaseCall,"\n",
                 "Position:",query_result$end,"\n",
                 "Variant Category:",query_result$status,"\n"
                 )
        }
      }else{
        paste0("Please Click on SNPs for details")
      }
    }else{
      paste0("Not valid search")
    }
  })
  #Click information 2: base calling quality control statistics
  output$info2<- renderText({
    #based on input$plot_click$x and input$plot_click$y value, search and return details of the variant
    mh_region=microhaplotype_data[which(microhaplotype_data$Amplicon==input$mh_region),]#get mh_region data
    if (is.null(input$plot_click)){
      return()
    }
    if (input$plot_click$y<=2.2 && input$plot_click$y>=1.8 && input$plot_click$x>=mh_region$START && input$plot_click$x<=mh_region$END){
      query_haplotype="hap2"
      query=ceiling(input$plot_click$x)#get the coordinate
      temp_data=figure_plot_data_single()[figure_plot_data_single()$haplotype==query_haplotype,]#extract the data
      query_idx=which(temp_data$end==query)#get correct index
      if (length(query_idx)>0){
        query_result=temp_data[query_idx,]
        paste0("Variant Calling QC Metric:\n",
               "Base Coverage:",query_result$ReadDepth,"\n",
               "Heterozygous Ratio:",format(query_result$HetRatio,digits = 2),"\n",
               "RefForward:",query_result$RefForward,"\n",
               "RefReverse:",query_result$RefReverse,"\n",
               "AltForward:",query_result$AltForward,"\n",
               "AltReverse:",query_result$AltReverse,"\n")
      }else{
        paste0("Please Click on SNPs for details")
      }
    }else if (input$plot_click$y<=1.2 && input$plot_click$y>=0.8 && input$plot_click$x>=mh_region$START && input$plot_click$x<=mh_region$END){
      query_haplotype="hap1"
      query=ceiling(input$plot_click$x)
      temp_data=figure_plot_data_single()[figure_plot_data_single()$haplotype==query_haplotype,]
      query_idx=which(temp_data$end==query)
      if (length(query_idx)>0){
        query_result=temp_data[query_idx,]
        paste0("Variant Calling QC Metric:\n",
               "Base Coverage:",query_result$ReadDepth,"\n",
               "Heterozygous Ratio:",format(query_result$HetRatio,digits = 2),"\n",
               "RefForward:",query_result$RefForward,"\n",
               "RefReverse:",query_result$RefReverse,"\n",
               "AltForward:",query_result$AltForward,"\n",
               "AltReverse:",query_result$AltReverse,"\n")
      }else{
        paste0("Please Click on SNPs for details")
      }
    }else{
      paste0("Not valid search")
    }
  })
  ##########Plot multiple MH_sample data simultaneously
  #render multiple plots if we observe change of mh region and population widget
  observeEvent(input$mh_region3,renderPlots(input$population2,input$mh_region3,total_read(),total_amplicon(),all_sample_name,microhaplotype_data,output))
  observeEvent(input$population2,renderPlots(input$population2,input$mh_region3,total_read(),total_amplicon(),all_sample_name,microhaplotype_data,output))
  #combine these plots together and display them on the webpage
  output$mh_figure_all <- renderUI({
    makePlotContainers(input$population2,input$mh_region3,all_sample_name)
  })
  ###########Show table of microhaplotype SNP data
  table_data=reactive({tbl=read.table(paste0("./www/data/pooled_vcf/mh_table/",
                                             input$mh_region2,
                                             "_snps_detail.txt"),
                                      header=TRUE)
  return(tbl)})
  output$mh_table<-DT::renderDataTable({datatable(table_data()) %>% formatStyle(columns='ReadDepthQCPass',
                                                                                color = styleEqual(c("Y","N"),
                                                                                                   c("black","red"))) %>% 
      formatStyle(columns='HetRatioQCPass',
                  color = styleEqual(c("Y","N"),
                                     c("black","red")))})
  ###########Show table of microhaplotype frequency table
  haplotype_data=reactive({freq_tbl=read.table(paste0("./www/data/pooled_vcf/mh_frequency/",
                                                      input$mh_region4,
                                                      "_haplotype.txt"),header=TRUE)
                          return(freq_tbl)})#read raw haplotype data
  output$mh_frequency_table<-DT::renderDataTable({
    temp_table=haplotype_data()[which(haplotype_data()$population==input$population3),]
    temp_summary=summary(as.factor(c(as.character(temp_table$Hap1),as.character(temp_table$Hap2))))
    
    out_table=data.frame("haplotype"=sapply(names(temp_summary),FUN = convert_to_haplotype),"count"=temp_summary,"frequency"=round(temp_summary/sum(temp_summary), 3))
    # no_sample_vec=c()
    # for (i in 1:dim(out_table)[1]){
    #   temp_query=as.character(out_table$haplotype[i])
    #   temp_sample_name=as.character(temp_table$sampleID[which(temp_table$Hap1==temp_query)])
    #   temp_sample_name=c(temp_sample_name,as.character(temp_table$sampleID[which(temp_table$Hap2==temp_query)]))
    #   no_sample_vec=c(no_sample_vec,length(unique(temp_sample_name)))
    # }
    #out_table=cbind(out_table,"sample_count"=no_sample_vec)
    out_table=rbind(out_table,data.frame("haplotype"="Total","count"=sum(temp_summary),"frequency"=sum(temp_summary/sum(temp_summary))))#"sample_count"=sum(no_sample_vec)))
    DT::datatable(out_table,rownames = FALSE,options=list(pageLength=20),width = 400)
  })
  #########QC figures for all three runs: use renderUI() to generate corresponding HTML tags or taglist  
  #output$qc_measure<-renderUI(
  #  #generate the tag to load pdf file
  #  tags$iframe(style="height:1000px; width:100%", 
  #              src=paste0("data/",input$dataset,"/qc/",input$qc_pdf,".pdf")))
  
}
#Load data
# data_set=reactive({ds=all_sample_name$run[which(as.character(all_sample_name$KenID)==as.character(input$sampleID))]
#                   return(ds)}
#                   )