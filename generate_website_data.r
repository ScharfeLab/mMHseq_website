library(xlsx)
library(sqldf)
#This script will scan through every filtered and phased vcf file for each sample in the directory
#Output necessary file to visualize on the mMHseq website
#############Set up necessary files
rm(list=ls())

# set input directory
work_directory = "" # local
setwd(work_directory)
source("generate_website_data_lib.r")

resource_dir= paste(work_directory, "MH_information_folder/", sep = "") #resource dir for MH and sample information
all_vcf_file= paste(work_directory,"Test_data/all_filtered.vcf", sep = "") #where phased sample is placed
individual_vcf_dir= paste(work_directory, "Test_data/filtered_and_phased/", sep = "") 
ken_mh_loc=read.xls(paste0(resource_dir,"Ken_MH_dbSNP.xlsx"))#contains the chromosome and gene information
all_sample_info=read.table(paste0(resource_dir,"all_sample_name.txt"),header=TRUE,stringsAsFactors = FALSE)
#Load 90 MH hg19 genomic coordinate (No primer)
microhaplotype_data=read.table(paste0(resource_dir,"mMHseq_90amp_no_primer.bed"),header=FALSE)
colnames(microhaplotype_data)=c("CHROM","START","END","Amplicon")#,"left_primer_start","left_primer_end","right_primer_start","right_primer_end")
#Load Ken's SNP information on each MH
ken_dbsnp=read.xlsx(paste0(resource_dir,"Ken_90MH_dbSNP_804.xlsx"), 1, header = T)#Ken's known snp ID on 96 MH region

# set output directory
# output_website_dir="/Users/iris/Documents/mMHseq_generated_data_test/"
output_website_dir= paste0(work_directory, "results/")

if (!dir.exists(output_website_dir)){#do a directory existence check
  dir.create(output_website_dir)#create this dir  
}

############ 1. Preprocessing the data, annotate each variant with 1000Genome frequency
#Read the vcf file which contains all filtered variants
all_snp_id=fread(all_vcf_file,skip = "#CHROM")
#############Annotate all identified SNP site with 1000Genome frequency in each MH region
snp_info_all_sample=all_snp_id[,1:5]#extract all snps on 90 amplicons found in our data
snp_summary_table=data.frame()#summarize the property of identified snps in all samples
snp_detail_table=data.frame()#annotate each snp with 1000genome data and ken's known snps information
for (each_mh in 1:dim(microhaplotype_data)[1]){#For each microhaplotype region, annotate each SNP within this region
  #Load 1000Genome frequnency data for this MH region
  mh_vcf_1000g=fread(paste0(resource_dir,"mH_region_1000g/",as.character(microhaplotype_data$Amplicon[each_mh]),".vcf"),skip = "#CHROM")#read 1000Genome data on this region
  #extract genomic interval of one microhaplotype region
  mh_region=microhaplotype_data[each_mh,]
  #Within each MH region, categorize each identifed SNP into "KK","dbSNP","Novel"
  temp_list=return_snp_summary_info(mh_region,snp_info_all_sample,ken_dbsnp,mh_vcf_1000g)
  snp_detail_table=rbind(snp_detail_table,temp_list[["full_list"]])#combine the mh summary into the table
  snp_summary_table=rbind(snp_summary_table,temp_list[["summary"]])#combine the mh detail into the table
}
snp_detail_table_1000genome=data.frame()#annotate the 1000G population frequency for all identified SNPs from all samples
for (each_snp in 1:dim(snp_detail_table)[1]){
  print(each_snp)
  snp_detail_table_1000genome=rbind(snp_detail_table_1000genome,annotate_1000genome(snp_detail_table[each_snp,],resource_dir))
}

write.table(snp_detail_table_1000genome,paste0(output_website_dir,"snp_detail_table_1000g.txt"),row.names=FALSE,quote=FALSE)
############END of pre-processing variant data

############ 2. Output data for website visualization ###################

#Output the annotated variant data to a user specified directory
snp_detail_table_1000genome=read.table(paste0(output_website_dir,"snp_detail_table_1000g.txt"),header=TRUE)

#read each filtered and phased vcf file of each sample into a list
all_sample_genotype=list()
for (each_sample in all_sample_info$KenID){
  mh_data=fread(paste0(individual_vcf_dir,each_sample,"_phased.vcf"),skip="#CHROM")#replace this line with fread(file,skip="#CHROM") later
  colnames(mh_data)[1]="CHROM"#,"POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT",each_sample)
  all_sample_genotype[[each_sample]]=mh_data
}
##########Output table for plotting purpose
#########CAUTION ###############
#mMHPoolV2-2: each_mh:15 and 31 has some GATK calling problems in S14, S30 due to low coverage of the amplicon
#some SNPs not called or called as ./.
#each_mh=1
#each_sample_ID=as.character(all_sample_info$KenID[1])
#counter <- 0
output_plot_data=data.frame()
for (each_mh in 1:dim(microhaplotype_data)[1]){#for each MH region
#for (each_mh in 9:9){#for each MH region
  print(paste0(each_mh," ",as.character(microhaplotype_data$Amplicon[each_mh])))
  mh_region=microhaplotype_data[each_mh,]#get a specific MH region
  amplicon_snp_output=data.frame()#output the SNP detail for each MH region
  mh_region_gene=ken_mh_loc[which(as.character(ken_mh_loc$MHapID)==as.character(mh_region$Amplicon)),]
  mh_haplotype_output=data.frame()#Output haplotype of each sample for this MH 
  #for each sample, get each SNP within the MH region and annotate them with 1000G allele frequency
  for (each_sample_ID in as.character(all_sample_info$KenID)){
    print(each_sample_ID)
    #get the genotype of all mH regions for this sample
    each_sample_genotype=all_sample_genotype[[each_sample_ID]]
    #determine which population this sample belongs to
    pop_info=as.character(all_sample_info$population[which(all_sample_info$KenID==each_sample_ID)])
    #get the identified snps annotation in 1000G in this MH region for this sample
    snp_in_mhregion=snp_detail_table_1000genome[which(snp_detail_table_1000genome$region==mh_region$Amplicon),]
    #extract the genotype data in this MH region from the sample
    snpid_in_mhregion_sampleidx=c()
    for (i in 1:dim(snp_in_mhregion)[1]){#for snp in each mH region within this specific sample
      if (snp_in_mhregion$snpid[i]!="."){#which means this is a SNP recorded in 1000G or dbSNP
        snpid_in_mhregion_sampleidx=c(snpid_in_mhregion_sampleidx,which(as.character(each_sample_genotype$ID)==snp_in_mhregion$snpid[i]))
      }else{#means this is a SNP not recorded in 1000Genome or dbSNP, use POS to find the index
        snpid_in_mhregion_sampleidx=c(snpid_in_mhregion_sampleidx,which(as.character(each_sample_genotype$POS)==snp_in_mhregion$POS[i]))
      }
    }
    ###########1. Extract genotype in this MH region for this sample from vcf file, but some variants in some sample will not be called
    each_sample_mh_genotype=each_sample_genotype[snpid_in_mhregion_sampleidx,]
    each_sample_mh_genotype=data.frame(each_sample_mh_genotype)
    
    #First check if the MH region of this sample contains any "./." call and remove it from the call set
    if (length(grep("^\\./\\.",each_sample_mh_genotype[,each_sample_ID])>0)){
      #remove no call loci in this sample
      each_sample_mh_genotype=each_sample_mh_genotype[-grep("^\\./\\.",each_sample_mh_genotype[,each_sample_ID]),]      
    }
    
    #contains no call in a specific sample, need to modify snp_in_mhregion a bit
    if (dim(each_sample_mh_genotype)[1]!=dim(snp_in_mhregion)[1]){
      temp_idx=c()
      for (i in 1:dim(each_sample_mh_genotype)[1]){
        temp_idx=c(temp_idx,which(as.character(snp_in_mhregion$POS)==as.character(each_sample_mh_genotype$POS[i])))
        # print("temp_idx")
        # print(temp_idx)
      }
      snp_in_mhregion=snp_in_mhregion[temp_idx,]#subset this 
    }
    
    ###########2. Filter/Tag sites with variants based on low DP|het_ratio, SNP only
    if (dim(snp_in_mhregion)[1]>0){#if there are snps identified in this region in this sample
      #First examine if there are indels in the call sets, currently the phasing only works in SNPs
      indel_idx=c()
      for (i in 1:dim(each_sample_mh_genotype)[1]){
        if (determine_variant_type(each_sample_mh_genotype[i,])=="STRindel" ||determine_variant_type(each_sample_mh_genotype[i,])=="nonSTRindel"){
          indel_idx=c(indel_idx,i)
        }
      }
      if (length(indel_idx)>0){
        #print("indel")
        each_sample_mh_genotype=each_sample_mh_genotype[-indel_idx,]
        snp_in_mhregion=snp_in_mhregion[-indel_idx,]
      }

      low_het_low_dp_check=get_low_het_dp(each_sample_mh_genotype,dp_cutoff=20,het_ratio=0.2,each_sample_ID)

      strand_allele_count=get_strand_allele(each_sample_mh_genotype,each_sample_ID)#Legacy, no such information from GATK UG

      snp_in_mhregion=cbind(snp_in_mhregion,low_het_low_dp_check[,-1])

      snp_in_mhregion=cbind(snp_in_mhregion,strand_allele_count)
      
      #convert genotype to haplotype in this sample
      sample_haplotype=extract_haplotype(each_sample_mh_genotype,sample_ID = each_sample_ID)
      sample_hap1=get_haplotype(sample_haplotype,1)
      sample_hap2=get_haplotype(sample_haplotype,2)

      mh_haplotype_output=rbind(mh_haplotype_output,data.frame("MHap"=mh_region$Amplicon,
                                                               "population"=pop_info,
                                                               "sampleID"=each_sample_ID,
                                                               "Hap1"=sample_hap1,
                                                               "Hap2"=sample_hap2))

      snp_in_mhregion=cbind(snp_in_mhregion,sample_haplotype)

      #first convert the data for plotting
      #convert each genotype/haplotype in this mH region to plotting format
      #Include base coverage and het ratio information
      converted_snp_for_plot=convert_to_geom_rect_format(snp_in_mhregion,mh_region)
      # print("converted_snp_for_plot")
      # print(converted_snp_for_plot)
      
      if(nrow(converted_snp_for_plot) > 0){
      check_annotation=data.frame()
      for (i in 1:dim(converted_snp_for_plot)[1]){
        if (converted_snp_for_plot$snpid[i]!="."){
          temp_idx=which(as.character(snp_in_mhregion$snpid)==converted_snp_for_plot$snpid[i])
        }else{
          temp_idx=which(snp_in_mhregion$POS==converted_snp_for_plot$end[i])
        }
        check_annotation=rbind(check_annotation,data.frame("ReadDepth"=snp_in_mhregion$ReadDepth[temp_idx],
                                                           "HetRatio"=snp_in_mhregion$HetRatio[temp_idx],
                                                           "RefForward"=snp_in_mhregion$RefForward[temp_idx],
                                                           "RefReverse"=snp_in_mhregion$RefReverse[temp_idx],
                                                           "AltForward"=snp_in_mhregion$AltForward[temp_idx],
                                                           "AltReverse"=snp_in_mhregion$AltReverse[temp_idx]))
      }
      converted_snp_for_plot=cbind(converted_snp_for_plot,check_annotation)
      converted_snp_for_plot=cbind(converted_snp_for_plot,
                                   "gene_or_locus"=rep(mh_region_gene$gene_or_locus,dim(converted_snp_for_plot)[1]),
                                   "Chr"=rep(mh_region_gene$Chr,dim(converted_snp_for_plot)[1]))
      #########Combine Output data for website plotting for mh region in this sample
      output_plot_data=rbind(output_plot_data,
                             data.frame("MHapID"=mh_region$Amplicon,
                                        "SampleID"=each_sample_ID,
                                        converted_snp_for_plot))
      #Output data for website table
      combined_output=convert_to_text_for_plot(snp_in_mhregion,
                                               mh_region,
                                               each_sample_mh_genotype,
                                               pop_info)
      amplicon_snp_output=rbind(amplicon_snp_output,combined_output[["table_out"]])}
    }else{#No SNPs identified in this region
      next()
    }
  }
  ##############Output table for the website
  amplicon_snp_output=cbind(amplicon_snp_output,
                            "gene_or_locus"=rep(mh_region_gene$gene_or_locus,dim(amplicon_snp_output)[1]),
                            "Chr"=rep(mh_region_gene$Chr,dim(amplicon_snp_output)[1]))
  amplicon_snp_output_reorder=amplicon_snp_output[,c("region","Chr","gene_or_locus","SampleID","Population",
                                                     "snpid","POS","REF","ALT","haplotype","VariantCategory",
                                                     "EAS","AMR","AFR","EUR","SAS",
                                                     "ReadDepthCheck","HetRatioCheck")]
  colnames(amplicon_snp_output_reorder)[1]="MHapID"#rename microhaplotype region
  colnames(amplicon_snp_output_reorder)[17]="ReadDepthQCPass"
  colnames(amplicon_snp_output_reorder)[18]="HetRatioQCPass"
  #amplicon_snp_output_reorder[,"HetRatio"]=format(amplicon_snp_output_reorder[,"HetRatio"],digits = 2)
  if (dir.exists(paste0(output_website_dir,"mh_table/"))){
    write.table(amplicon_snp_output_reorder,
                paste0(paste0(output_website_dir,"mh_table/"),as.character(mh_region$Amplicon),"_snps_detail.txt"),
                row.names = FALSE,
                col.names = c("MHapID", "Chr", "gene_or_locus", "SampleID", "Population", "snpid", "POS", "REF", "ALT", 
                              "genotype", "VariantCategory", "EAS", "AMR", "AFR", "EUR", "SAS", 
                              "ReadDepthQCPass", "HetRatioQCPass"),
                quote=FALSE,
                sep = "\t")    
  }else{
    dir.create(paste0(output_website_dir,"mh_table/"))
    write.table(amplicon_snp_output_reorder,
                paste0(output_website_dir,"mh_table/",as.character(mh_region$Amplicon),"_snps_detail.txt"),
                row.names = FALSE,
                col.names = c("MHapID", "Chr", "gene_or_locus", "SampleID", "Population", "snpid", "POS", "REF", "ALT", 
                              "genotype", "VariantCategory", "EAS", "AMR", "AFR", "EUR", "SAS", 
                              "ReadDepthQCPass", "HetRatioQCPass"),
                
                quote=FALSE,
                sep = "\t")       
  }
  
 # if (dir.exists(paste0(output_website_dir,"mh_frequency"))){
 #   write.table(mh_haplotype_output,paste0(output_website_dir,"mh_frequency/",
 #                                          as.character(mh_region$Amplicon),"_haplotype.txt"),
 #               row.names=FALSE,quote=FALSE,
 #               sep = "\t")
 # }else{
 #   dir.create(paste0(output_website_dir,"mh_frequency"))
 #   write.table(mh_haplotype_output,paste0(output_website_dir,"mh_frequency/",
 #                                          as.character(mh_region$Amplicon),"_haplotype.txt"),
 #               row.names=FALSE,
 #               quote=FALSE,
 #               sep = "\t")
 # }
}
if (dir.exists(paste0(output_website_dir,"mh_figure_data/"))){
  write.table(output_plot_data,
              paste0(output_website_dir,"mh_figure_data/website_plotdata.txt"),
              row.names=FALSE,
              quote=FALSE,
              sep = "\t")        
}else{
  dir.create(paste0(output_website_dir,"mh_figure_data/"))
  write.table(output_plot_data,
              paste0(output_website_dir,"mh_figure_data/website_plotdata.txt"),
              row.names=FALSE,
              quote=FALSE,
              sep = "\t")           
}
# missed_ken_mh_id=setdiff(unlist(strsplit(as.character(ken_mh_snpid$SNPs),split="/")),
# as.character(temp_list[["full_list"]]$snpid))#examine which Ken's SNP is not called in our data set
# missed_ken_snp=rbind(missed_ken_snp,data.frame("MH_name"=rep(mh_region$Amplicon,length(missed_ken_mh_id)),
# "SNP_ID"=missed_ken_mh_id))
# ken_mh_snpid=ken_dbsnp[which(as.character(ken_dbsnp$MHapID)==as.character(mh_region$Amplicon)),]

# bad_amplicon=read.table(paste0(resource_dir,"six_removed_amplicon.txt"))#6 bad amplicons on MH removed
# for (each in bad_amplicon$V1){
#   ken_dbsnp=ken_dbsnp[-which(ken_dbsnp$MHapID==as.character(each)),]
# }#remove six bad amplicons

#### novelCP ####
#### check for novels and add counterpart
#### if both hap1 and hap2 are novels at one position for a sample, skip to next ####

website_plotdata_original <- read.table(paste0(output_website_dir,"mh_figure_data/website_plotdata.txt"), header = TRUE, stringsAsFactors = FALSE)
out <- NULL

require(sqldf)
SampleID <- sqldf('SELECT DISTINCT SampleID FROM website_plotdata_original')
length(SampleID)
counter <- 0
counter_all_novelcp <- 0

for (i in 1:length(SampleID$SampleID)) {
  query <- paste0('SELECT * FROM website_plotdata_original WHERE snpid = \'.\' AND status = \'Novel\' AND SampleID = \'', SampleID$SampleID[i], '\'')
  novel_by_sample <- sqldf(query)
  print(SampleID$SampleID[i])
  
  if(nrow(novel_by_sample) > 0){
    for (j in 1:length(novel_by_sample$haplotype)) {
      query <- paste0('SELECT COUNT (*) FROM novel_by_sample WHERE end = \'', novel_by_sample$end[j], '\'')
      counter <- sqldf(query)
      print(paste0("counter = ", counter))
      if (counter > 1) { # if counterpart of novel hap is also novel, skip
        print(paste0(novel_by_sample$SampleID[j], " ", novel_by_sample$haplotype[j], " ", novel_by_sample$end[j]))
        counter <- 0
        temp <- NULL
      }else {
        temp <- NULL
        if(identical(novel_by_sample$haplotype[j],"hap1")) {
          print("1")
          
          temp <- novel_by_sample[j, ]
          temp$status <- "NovelCP"
          temp$haplotype <- "hap2"
          temp$BaseCall <- novel_by_sample[j, ]$REF
          counter <- 0
        }
        else if (identical(novel_by_sample$haplotype[j],"hap2")){
          print("2")
          
          temp <- novel_by_sample[j, ]
          temp$status <- "NovelCP"
          temp$haplotype <- "hap1"
          temp$BaseCall <- novel_by_sample[j, ]$REF
        }
        counter_all_novelcp <- counter_all_novelcp + 1
        print(counter_all_novelcp)
      }
      out <- rbind(out, temp)
    }
  }
  
}

print(paste0("counter_all_novelcp = ", counter_all_novelcp))
# print(out)
# write.table(as.data.frame(out), file = paste0(work_directory, "out.txt"), quote = F, sep = "\t", row.names = F, col.names = F)
new <- rbind(website_plotdata_original, out)
new <- sqldf("SELECT * FROM new ORDER BY MHapID, SampleID, end")
write.table(as.data.frame(new), file = paste0(work_directory, "out.txt"), quote = F, sep = "\t", row.names = F, col.names = T)
# out.txt file contains NovelCP, add to original website data (from last version) to proceed to the next step 