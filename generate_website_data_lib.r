library(data.table)#contain fread() function
library(gdata)#read excel file
library(ggplot2)#plot function
library(reshape2)#contains melt() function
library(gridExtra)#for multiple graphs on one plot
library(RColorBrewer)
return_snp_summary_info<-function(micro_haplotype_coord,
                                  snp_data,
                                  known_snp_id,
                                  pubic_database_id){
  #annotate the SNPs:
  #micro_haplotype_coord: MH genomic interval
  #snp_data: all sites identified in our sample
  #known_snp_id: Ken's known snps on that MH
  #public_data_base_id:check if the SNP is recorded in 1000G
  output_list=list()
  matched_snp_id=data.frame()
  ken_provide=0
  new_snp=0
  in_1000genome=0#if the snp is in 1000genome
  in_dbsnp=0#if snp is in dbSNP
  not_recorded=0
  #Loop through all SNPs identified from our data
  for (i in 1:dim(snp_data)[1]){#determine if one SNP in the mH identified in our samples is in Ken's known list and 1000Genomes or not
    #check if the snp falls on this MH region
    if (as.character(snp_data$`#CHROM`[i])==as.character(micro_haplotype_coord$CHROM)){#on the same chromosome
      #The SNP does not fall on primer region
      if (as.numeric(snp_data$POS[i]) %in% seq(micro_haplotype_coord$START,micro_haplotype_coord$END)){
        #determine if this SNP is in Ken's list or not
        if (snp_data$ID[i] %in% known_snp_id$snpid){
          temp=data.frame(snpid=as.character(snp_data$ID[i]),
                          region=micro_haplotype_coord$Amplicon,
                          "known"="Ken",
                          "POS"=snp_data$POS[i])
          matched_snp_id=rbind(matched_snp_id,temp)
          ken_provide=ken_provide+1
        }else{
          #determine if the snp is in 1000genome or not
          if (snp_data$ID[i] %in% pubic_database_id$ID){
            in_1000genome=in_1000genome+1
            temp=data.frame(snpid=as.character(snp_data$ID[i]),region=micro_haplotype_coord$Amplicon,"known"="dbSNP","POS"=snp_data$POS[i])
          }else{
            if (as.character(snp_data$ID[i])=="."){
              temp=data.frame(snpid=as.character(snp_data$ID[i]),region=micro_haplotype_coord$Amplicon,"known"="Novel","POS"=snp_data$POS[i])
              not_recorded=not_recorded+1
            }else{
              temp=data.frame(snpid=as.character(snp_data$ID[i]),region=micro_haplotype_coord$Amplicon,"known"="dbSNP","POS"=snp_data$POS[i])
              in_dbsnp=in_dbsnp+1
            }
          }
          matched_snp_id=rbind(matched_snp_id,temp)
          new_snp=new_snp+1
        }
      }
    }
  }
  output_list[["full_list"]]=matched_snp_id
  output_list[["summary"]]=data.frame("Amplicon"=micro_haplotype_coord$Amplicon,
                                      "Ken"=ken_provide,
                                      "New"=new_snp,
                                      "in1000Genome"=in_1000genome,
                                      "indbSNP"=in_dbsnp,
                                      "NotRecorded"=not_recorded)
  return(output_list)
}
extract_race_frequency_1000genome<-function(vcf_info,col_no,max_allele_freq=TRUE){
  temp=unlist(strsplit(unlist(strsplit(vcf_info$INFO,split=";"))[col_no],split="="))
  freq_vec=unlist(strsplit(temp[2],split=","))
  if (length(freq_vec)==1){
    out=data.frame("population"=sub("_AF","",temp[1]),"allele_freq"=as.numeric(temp[2]))    
  }else{#multi-allelic SNP
    if (max_allele_freq==TRUE){
      out=data.frame("population"=sub("_AF","",temp[1]),"allele_freq"=max(as.numeric(freq_vec[1]),as.numeric(freq_vec[2])))      
    }else{
      out=data.frame("population"=sub("_AF","",temp[1]),"allele_freq"=paste(freq_vec[1],freq_vec[2],sep=","))
    }
  }
  return(out)
}
annotate_1000genome<-function(x,dir){
  #annotate each SNP variant with 1000Genome frequency
  vcf_file=paste0(dir,"mH_region_1000g/",as.character(x$region),".vcf")
  vcf_1000genome=fread(vcf_file,skip = "#CHROM")#obtain 1000Genome vcf data
  snp_1000genome_idx=which(as.character(vcf_1000genome$ID)==as.character(x$snpid))
  if (length(snp_1000genome_idx)>0){
    eas_af=extract_race_frequency_1000genome(vcf_1000genome[snp_1000genome_idx,],6,FALSE)
    amr_af=extract_race_frequency_1000genome(vcf_1000genome[snp_1000genome_idx,],7,FALSE)#American
    afr_af=extract_race_frequency_1000genome(vcf_1000genome[snp_1000genome_idx,],8,FALSE)#African
    eur_af=extract_race_frequency_1000genome(vcf_1000genome[snp_1000genome_idx,],9,FALSE)#Europe
    sas_af=extract_race_frequency_1000genome(vcf_1000genome[snp_1000genome_idx,],10,FALSE)#South asian
    af_combined=data.frame("EAS"=eas_af$allele_freq,"AMR"=amr_af$allele_freq,"AFR"=afr_af$allele_freq,
                           "EUR"=eur_af$allele_freq,"SAS"=sas_af$allele_freq)
    x=cbind(x,af_combined)
  }else{
    af_combined=data.frame("EAS"=NA,"AMR"=NA,"AFR"=NA,"EUR"=NA,"SAS"=NA)
    x=cbind(x,af_combined)
  }
  return(x)
}
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}
determine_variant_type<-function(vcf_line){
  if (nchar(vcf_line$REF)==1 && nchar(vcf_line$ALT)==1){
    return("SNP")
  }else if (nchar(vcf_line$REF)>1 || nchar(vcf_line$ALT)>1 ){#could be an indel or MNP
    if (grepl("^[ATCG]{1},[ATCG]{1}$",vcf_line$ALT)){
      return("MNP")
    }else{
      if (grepl(";STR$",vcf_line$INFO)){
        return("STRindel")
      }else{
        return("nonSTRindel")
      }
    }
  } 
  
}
extract_haplotype<-function(input_genotype,sample_ID){
  #Currently only SNP could be phased, so only extract phased SNP information, do not extract 
  output_haplotype=c()
  ref_out=c()
  alt_out=c()
  for (each_snp in 1:dim(input_genotype)[1]){#loop through each snp genotype in this MH region
    #print(each_snp)
    if (determine_variant_type(input_genotype[each_snp,])!="STRindel" && determine_variant_type(input_genotype[each_snp,])!="nonSTRindel"){
      ref_base=as.character(input_genotype$REF[each_snp])
      alt_base=unlist(strsplit(as.character(input_genotype$ALT[each_snp]),split=","))
      base_pool=c(ref_base,alt_base)
      ref_out=c(ref_out,as.character(input_genotype$REF[each_snp]))
      alt_out=c(alt_out,as.character(input_genotype$ALT[each_snp]))
      if (grepl("HP",as.character(input_genotype$FORMAT[each_snp]))){#which means this genotype is phased
        #Called genotype for this variant
        genotype=unlist(strsplit(as.character(input_genotype[,sample_ID][each_snp]),split=":"))[1]
        genotype=unlist(strsplit(genotype,split="/"))
        #convert the number to base character
        genotype_text=c()
        for (i in 1:length(genotype)){
          genotype_text=c(genotype_text,base_pool[as.numeric(genotype[i])+1])
        }
        #Haplotype for this variant
        temp_geno=unlist(strsplit(as.character(input_genotype[,sample_ID][each_snp]),split=":"))[5]
        temp_geno=unlist(strsplit(temp_geno,split=","))
        if (grepl("-1",temp_geno[1])){
          output_haplotype=c(output_haplotype,paste0(genotype_text[1],"|",genotype_text[2]))
        }else{
          output_haplotype=c(output_haplotype,paste0(genotype_text[2],"|",genotype_text[1]))
        }
      }else{
        genotype=unlist(strsplit(as.character(input_genotype[,sample_ID][each_snp]),split=":"))[1]
        genotype=unlist(strsplit(genotype,split="/"))
        if (genotype[1]!="."){#there is a genotype call at this site
          genotype_text=c()
          for (i in 1:length(genotype)){
            genotype_text=c(genotype_text,base_pool[as.numeric(genotype[i])+1])
          }           
        }else{#keep "./."
          genotype_text=genotype
        }
        output_haplotype=c(output_haplotype,as.character(paste0(genotype_text[1],"|",genotype_text[2])))
      }      
    }else{next()}

  }
  output_haplotype=data.frame(haplotype=output_haplotype,"REF"=ref_out,"ALT"=alt_out)
  return(output_haplotype)
}
convert_to_geom_rect_format<-function(snp_data_frame,mh_data){
  #for each snp data frame, determine its haplotype and convert it to geom_rect format for plotting later
  out=data.frame()
  for (each_line in 1:dim(snp_data_frame)[1]){#for each SNP record
    temp=snp_data_frame[each_line,]
    temp_hap=unlist(strsplit(as.character(temp$haplotype),split="\\|"))#get haplotype for this SNP site
    ref_base=as.character(temp$REF)
    status_hap1=""
    status_hap2=""
    #haplotype 1
    if (temp_hap[1]==ref_base && temp$known=="Ken"){#non variant from Ken
      status_hap1="KK(WT)"
    }else if (temp_hap[1]!=ref_base && temp$known=="Ken"){#variant from Ken
      status_hap1="KK(Variant)"
    }else if (temp_hap[1]!=ref_base && temp$known=="dbSNP"){
      status_hap1="dbSNP"
    }else if (temp_hap[1]!=ref_base && temp$known=="Novel"){
      status_hap1="Novel"
    }
    #haplotype 2
    if (temp_hap[2]==ref_base && temp$known=="Ken"){#non variant
      status_hap2="KK(WT)"
    }else if (temp_hap[2]!=ref_base && temp$known=="Ken"){#variant
      status_hap2="KK(Variant)"
    }else if (temp_hap[2]!=ref_base && temp$known=="dbSNP"){
      status_hap2="dbSNP"
    }else if (temp_hap[2]!=ref_base && temp$known=="Novel"){
      status_hap2="Novel"
    }
    if (status_hap1!=""){
      converted_hap1=data.frame("snpid"=temp$snpid,"haplotype"="hap1","start"=temp$POS-1,"end"=temp$POS,
                                "status"=status_hap1,"EAS"=temp$EAS,"AMR"=temp$AMR,"AFR"=temp$AFR,"EUR"=temp$EUR,
                                "SAS"=temp$SAS,"REF"=temp$REF,"ALT"=temp$ALT,
                                "BaseCall"=temp_hap[1],"haplotype"=temp$haplotype)
      out=rbind(out,converted_hap1)
    }
    if (status_hap2!=""){
      converted_hap2=data.frame("snpid"=temp$snpid,"haplotype"="hap2","start"=temp$POS-1,"end"=temp$POS,
                                "status"=status_hap2,"EAS"=temp$EAS,"AMR"=temp$AMR,"AFR"=temp$AFR,"EUR"=temp$EUR,
                                "SAS"=temp$SAS,"REF"=temp$REF,"ALT"=temp$ALT,
                                "BaseCall"=temp_hap[2],"haplotype"=temp$haplotype)      
      out=rbind(out,converted_hap2)
    }
  }
  return(out)
}
convert_to_text_for_plot<-function(snp_data_frame,
                                   mh_data,
                                   each_sample_mh_genotype,
                                   pop_info){
  #This script will convert the SNP record in SNP_data_frame into text format for visualization
  all_out=list()#output two data.frame()
  #out=data.frame()#output for text plotting
  out2=data.frame()#output for text table
  samplename=colnames(each_sample_mh_genotype)[dim(each_sample_mh_genotype)[2]]#Get sample name
  for (each_line in 1:dim(snp_data_frame)[1]){#read each snp record from data
    temp=snp_data_frame[each_line,]
    temp_geno=unlist(strsplit(as.character(temp$haplotype),split="\\|"))
    ref_base=as.character(temp$REF)
    #get het ratio information
    temp_full_info=as.character(each_sample_mh_genotype[each_line,dim(each_sample_mh_genotype)[2]])
    het_ratio=unlist(strsplit(temp_full_info,":"))[2]
    if (snp_data_frame$known[each_line]=="NewInDatabase"){# && snp_data_frame$known[each_line]!="LeftPrimer" && snp_data_frame$known[each_line] !="RightPrimer"){
      if (temp_geno[1]!=ref_base || temp_geno[2]!=ref_base){#The genotype contains non-reference base
        out2=rbind(out2,cbind(snp_data_frame[each_line,],
                              "SampleID"=samplename,
                              "Population"=pop_info,
                              "VariantCategory"="dbSNP"))
      }
    }else if (snp_data_frame$known[each_line]=="NewNotRecorded"){
      if (temp_geno[1]!=ref_base || temp_geno[2]!=ref_base){
        out2=rbind(out2,cbind(snp_data_frame[each_line,],
                              "SampleID"=samplename,
                              "Population"=pop_info,
                              "VariantCategory"="Novel"))          
      }
    }else if (snp_data_frame$known[each_line]=="Ken"){
      out2=rbind(out2,cbind(snp_data_frame[each_line,],
                            "SampleID"=samplename,
                            "Population"=pop_info,
                            "VariantCategory"="KK"))
    }
  }
  #all_out[["text_out"]]=out
  all_out[["table_out"]]=out2
  return(all_out)
}
get_low_het_dp<-function(each_sample_mh_genotype,dp_cutoff=20,het_ratio_cutoff=0.2,sample_name){
  #record and filter the variant site based on read depth and het ratio cutoff (default: DP<20||het_ratio<20%)
  output_df=data.frame()
  for (i in 1:dim(each_sample_mh_genotype)[1]){
    temp_genotype_data=unlist(strsplit(as.character(each_sample_mh_genotype[[sample_name]][i]),split=":"))
    read_depth=as.numeric(temp_genotype_data[3])
    temp_genotype=temp_genotype_data[1]
    temp_df=data.frame("ID"=each_sample_mh_genotype$ID[i])
    #Examine read depth for this site
    if (read_depth<dp_cutoff){#first check if the read depth is larger than cutoff
      temp_df=cbind(temp_df,data.frame("ReadDepthCheck"="N","ReadDepth"=read_depth))
    }else{
      temp_df=cbind(temp_df,data.frame("ReadDepthCheck"="Y","ReadDepth"=read_depth))
    }
    #Examine het ratio for this site
    if (temp_genotype!="0/0" && temp_genotype!="1/1" && temp_genotype!="2/2" & temp_genotype!="./."){#this is a het site, check het ratio
      het_count=min(as.numeric(unlist(strsplit(temp_genotype_data[2],split=","))))#which value is smaller
      het_ratio=het_count/read_depth#get non-ref/total rd ratio
      if (het_ratio<het_ratio_cutoff){#if the non-ref/ref ratio < cutoff, record this site as filtered
        temp_df=cbind(temp_df,data.frame("HetRatioCheck"="N","HetRatio"=het_ratio))
      }else{
        temp_df=cbind(temp_df,data.frame("HetRatioCheck"="Y","HetRatio"=het_ratio))        
      }
    }else{#homo site
      temp_df=cbind(temp_df,data.frame("HetRatioCheck"="Y","HetRatio"=0))      
    }
    output_df=rbind(output_df,temp_df)
  }
  return(output_df)
}
get_strand_allele<-function(each_sample_mh_genotype,sample_name){
  #Not applicable in our case
  output_df=data.frame()
  for (i in 1:dim(each_sample_mh_genotype)[1]){
    temp_genotype_data=unlist(strsplit(as.character(each_sample_mh_genotype[[sample_name]][i]),split=":"))
    # if (grepl(pattern = "SAC",as.character(each_sample_mh_genotype$FORMAT[i]))){#if this variant contains SAC record
    #   temp_SAC=unlist(strsplit(temp_genotype_data[length(temp_genotype_data)],split=","))#should be the last item
    #   temp_df=data.frame("RefForward"=as.numeric(temp_SAC[1]),
    #                      "RefReverse"=as.numeric(temp_SAC[2]),
    #                      "AltForward"=as.numeric(temp_SAC[3]),
    #                      "AltReverse"=as.numeric(temp_SAC[4]))
    # }else{
      temp_df=data.frame("RefForward"=NA,
                         "RefReverse"=NA,
                         "AltForward"=NA,
                         "AltReverse"=NA)      
    #}
    output_df=rbind(output_df,temp_df)
  }
  return(output_df)  
}
get_haplotype<-function(sample_haplotype,hap_no){
  output_string=c()
  for (i in 1:length(sample_haplotype$haplotype)){
    hap=unlist(strsplit(as.character(sample_haplotype$haplotype[i]),split="\\|"))[hap_no]
    output_string=paste0(output_string,hap)
  }
  return(output_string)
}
