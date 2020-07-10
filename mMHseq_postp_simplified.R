rm(list=ls())

# set directories
work_directory = "" # local
post_processing_dir = paste0(work_directory, "post_processing/")
setwd(work_directory)
res_update_dir = paste0(work_directory, "res_update_info/")

ken_resource_dir = paste0(res_update_dir, "Ken_90MH_dbSNP_annotated.txt")
plotdata_complete_dir = paste0(post_processing_dir, "website_plotdata.txt")

final_output_dir = paste0(post_processing_dir, "additional_grey_bars_annotated.txt")
website_plotdata_complete_dir = paste0(post_processing_dir, "website_plotdata_complete.txt")

ken_resource <- read.table(ken_resource_dir, header = TRUE)
plotdata_complete <- read.table(plotdata_complete_dir, header = TRUE)

require(sqldf)
# get all mh IDs and sample IDs
MHapID <- sqldf('SELECT DISTINCT mhap FROM ken_resource')
length(MHapID$mhap)

SampleID <- sqldf('SELECT DISTINCT SampleID FROM plotdata_complete')
length(SampleID$SampleID)

# length(ken_resource$pos)
# length(ken_resource$snpid)

new_grey <- NULL
for (i in 1:length(SampleID$SampleID)) {
  query <- paste0('SELECT * FROM plotdata_complete WHERE SampleID = \'', SampleID$SampleID[i], '\'')
  # query <- paste0('SELECT * FROM plotdata_complete WHERE SampleID = \'JK5054\'')
  by_sample <- sqldf(query)
  for (j in 1:length(ken_resource$snpid)) {
    if (ken_resource$snpid[j] == '.') {
      if (!(ken_resource$pos[j] %in% by_sample$end)) {
        hap1 <- NULL
        hap2 <- NULL
        
        hap1 <- data.frame("MHapID" = ken_resource$mhap[j], "SampleID" = SampleID$SampleID[i], "snpid" = ken_resource$snpid[j], 
                           "haplotype" = "hap1", "start" = ken_resource$pos[j]-1, "end" = ken_resource$pos[j], 
                           "status" = "Purple_placeholder", "EAS" = ken_resource$EAS[j], "AMR" = ken_resource$AMR[j], 
                           "AFR" = ken_resource$AFR[j], "EUR" = ken_resource$EUR[j], "SAS" = ken_resource$SAS[j],
                           "REF" = ken_resource$REF[j], "ALT" = ken_resource$ALT[j], "BaseCall" = ken_resource$REF[j], 
                           "haplotype.1" = "NA", "ReadDepth" = "NA", "HetRatio" = "NA",
                           "RefForward" = "NA", "RefReverse" = "NA", "AltForward" = "NA", "AltReverse" = "NA", 
                           "gene_or_locus" = "NA", "Chr" = ken_resource$chr[j],
                           stringsAsFactors = FALSE)
        hap2 <- data.frame("MHapID" = ken_resource$mhap[j], "SampleID" = SampleID$SampleID[i], "snpid" = ken_resource$snpid[j], 
                           "haplotype" = "hap2", "start" = ken_resource$pos[j]-1, "end" = ken_resource$pos[j], 
                           "status" = "Purple_placeholder", "EAS" = ken_resource$EAS[j], "AMR" = ken_resource$AMR[j], 
                           "AFR" = ken_resource$AFR[j], "EUR" = ken_resource$EUR[j], "SAS" = ken_resource$SAS[j],
                           "REF" = ken_resource$REF[j], "ALT" = ken_resource$ALT[j], "BaseCall" = ken_resource$REF[j], 
                           "haplotype.1" = "NA", "ReadDepth" = "NA", "HetRatio" = "NA",
                           "RefForward" = "NA", "RefReverse" = "NA", "AltForward" = "NA", "AltReverse" = "NA", 
                           "gene_or_locus" = "NA", "Chr" = ken_resource$chr[j],
                           stringsAsFactors = FALSE)
        
        new_info <- NULL
        new_info <- rbind(hap1, hap2)
        new_grey <- rbind(new_grey, new_info)
      }
    } else if (!(ken_resource$snpid[j] %in% by_sample$snpid) & !(ken_resource$pos[j] %in% by_sample$end)) {
      hap1 <- NULL
      hap2 <- NULL
      
      hap1 <- data.frame("MHapID" = ken_resource$mhap[j], "SampleID" = SampleID$SampleID[i], "snpid" = ken_resource$snpid[j], 
                         "haplotype" = "hap1", "start" = ken_resource$pos[j]-1, "end" = ken_resource$pos[j], 
                         "status" = "Ancestral_placeholder", "EAS" = ken_resource$EAS[j], "AMR" = ken_resource$AMR[j], 
                         "AFR" = ken_resource$AFR[j], "EUR" = ken_resource$EUR[j], "SAS" = ken_resource$SAS[j],
                         "REF" = ken_resource$REF[j], "ALT" = ken_resource$ALT[j], "BaseCall" = ken_resource$REF[j], 
                         "haplotype.1" = "NA", "ReadDepth" = "NA", "HetRatio" = "NA",
                         "RefForward" = "NA", "RefReverse" = "NA", "AltForward" = "NA", "AltReverse" = "NA", 
                         "gene_or_locus" = "NA", "Chr" = ken_resource$chr[j],
                         stringsAsFactors = FALSE)
      hap2 <- data.frame("MHapID" = ken_resource$mhap[j], "SampleID" = SampleID$SampleID[i], "snpid" = ken_resource$snpid[j], 
                         "haplotype" = "hap2", "start" = ken_resource$pos[j]-1, "end" = ken_resource$pos[j], 
                         "status" = "Ancestral_placeholder", "EAS" = ken_resource$EAS[j], "AMR" = ken_resource$AMR[j], 
                         "AFR" = ken_resource$AFR[j], "EUR" = ken_resource$EUR[j], "SAS" = ken_resource$SAS[j],
                         "REF" = ken_resource$REF[j], "ALT" = ken_resource$ALT[j], "BaseCall" = ken_resource$REF[j], 
                         "haplotype.1" = "NA", "ReadDepth" = "NA", "HetRatio" = "NA",
                         "RefForward" = "NA", "RefReverse" = "NA", "AltForward" = "NA", "AltReverse" = "NA", 
                         "gene_or_locus" = "NA", "Chr" = ken_resource$chr[j],
                         stringsAsFactors = FALSE)
      
      new_info <- NULL
      new_info <- rbind(hap1, hap2)
      new_grey <- rbind(new_grey, new_info)
    }
  }
}

#### novelCP ####

# new_grey <- read.table(paste0(work_directory, "website_plotdata_JK5054_full.txt"), header = TRUE, stringsAsFactors = FALSE)
# out <- NULL
# 
# for (i in 1:length(SampleID$SampleID)) {
#   query <- paste0('SELECT * FROM new_grey WHERE snpid = \'.\' AND status = \'Novel\' AND SampleID = \'', SampleID$SampleID[i], '\'')
#   novel_by_sample <- sqldf(query)
#   
#   for (j in 1:length(novel_by_sample$haplotype)) {
#     temp <- NULL
#     if(novel_by_sample$haplotype[j] == "hap1") {
#       print("1")
#       
#       temp <- novel_by_sample[j, ]
#       temp$status <- "NovelCP"
#       temp$haplotype <- "hap2"
#     }
#     else if (novel_by_sample$haplotype[j] == "hap2"){
#       print("2")
# 
#       temp <- novel_by_sample[j, ]
#       temp$status <- "NovelCP"
#       temp$haplotype <- "hap1"
#     }
#     out <- rbind(out, temp)
#   }
#   
# }
# 
# print(out)
# write.table(as.data.frame(out), file = paste0(work_directory, "out.txt"), quote = F, sep = "\t", row.names = F, col.names = F)
write.table(as.data.frame(new_grey), file = final_output_dir, quote= F, sep= "\t", row.names= F)

# write.table(as.data.frame(new_grey), file = paste0(work_directory, "website_plotdata.txt"), quote= F, sep= "\t", row.names= F, append = TRUE)
plotdata_complete_new <- NULL
plotdata_complete_new <- rbind(plotdata_complete, new_grey)
plotdata_complete_new <- sqldf("SELECT * FROM plotdata_complete_new ORDER BY MHapID, SampleID, end")
#plotdata_complete_new <- rbind(plotdata_complete_new, out)
write.table(as.data.frame(plotdata_complete_new), file = website_plotdata_complete_dir, quote= F, sep= "\t", row.names= F)


######## generate freq table ########

freq_output_dir = paste0(post_processing_dir, "mh_frequency/")
sample_info_dir = paste0(work_directory, "MH_information_folder/all_sample_name.txt")

plotdata_complete <- data.frame()
plotdata_complete <- read.table(website_plotdata_complete_dir, header = TRUE, stringsAsFactors = FALSE)
sample_info <- read.table(sample_info_dir, header = TRUE, stringsAsFactors = FALSE)
# pop <- sqldf('SELECT DISTINCT population FROM sample_info')

require(sqldf)
sample_pop <- sqldf('SELECT KenID, population FROM sample_info')
basecall_hap1 <- sqldf('SELECT MHapID, SampleID, haplotype, end, BaseCall FROM plotdata_complete WHERE haplotype = \'hap1\'')
basecall_hap2 <- sqldf('SELECT MHapID, SampleID, haplotype, end, BaseCall FROM plotdata_complete WHERE haplotype = \'hap2\'')
MHapID <- sqldf('SELECT DISTINCT MHapID FROM basecall_hap1')

position_table <- data.frame()

if (!dir.exists(freq_output_dir)) {
  dir.create(freq_output_dir)
}

for (i in 1:length(MHapID$MHapID)) {
  freq_table <- data.frame()
  freq_table_position <- data.frame()
  
  print(paste0(i, " ", MHapID$MHapID[i]))
  
  query_hap1 <- paste0("SELECT MHapID, SampleID, haplotype, end, BaseCall FROM basecall_hap1 WHERE MHapID = \'",
                       MHapID$MHapID[i], "\'")
  query_hap2 <- paste0("SELECT MHapID, SampleID, haplotype, end, BaseCall FROM basecall_hap2 WHERE MHapID = \'",
                       MHapID$MHapID[i], "\'")
  basecall_by_MHapID_hap1 <- sqldf(query_hap1)
  basecall_by_MHapID_hap2 <- sqldf(query_hap2)
  
  for (j in 1:length(sample_pop$KenID)) {
    position <- NULL
    hap1 <- NULL
    hap2 <- NULL
    starts_at <- NULL
    ends_at <- NULL
    pop <- sample_pop$population[j]
    freq_entry <- data.frame()
    
    print(paste0(sample_pop$KenID[j], " ", pop))
    
    query_hap1 <- paste0("SELECT MHapID, SampleID, haplotype, end, BaseCall FROM basecall_by_MHapID_hap1 WHERE SampleID = \'", sample_pop$KenID[j], "\' ORDER BY end ASC")
    query_hap2 <- paste0("SELECT MHapID, SampleID, haplotype, end, BaseCall FROM basecall_by_MHapID_hap2 WHERE SampleID = \'", sample_pop$KenID[j], "\' ORDER BY end ASC")
    basecall_single_hap1 <- sqldf(query_hap1)
    basecall_single_hap2 <- sqldf(query_hap2)
    # pop <- sample_pop$population[j]
    
    for (k in 1:length(basecall_single_hap1$end)) {
      #hap1 <- c(hap1, basecall_hap1$BaseCall[k])
      #hap2 <- c(hap2, basecall_hap2$BaseCall[k])
      if (k == 1) {
        starts_at <- basecall_single_hap1$end[1]
      }
      if (k == length(basecall_single_hap1$end)) {
        ends_at <- basecall_single_hap1$end[k]
      }
      if (k == 1) {
        position <- basecall_single_hap1$end[k]
      } else {
        position <- paste0(position, "-", basecall_single_hap1$end[k])
      }
      hap1 <- paste0(hap1, basecall_single_hap1$BaseCall[k])
      hap2 <- paste0(hap2, basecall_single_hap2$BaseCall[k])
    }
    freq_entry <- data.frame("MHap" = MHapID$MHapID[i], "population"  = pop, "sampleID" = sample_pop$KenID[j], 
                             "Hap1" = hap1, "Hap2" = hap2,
                             # "starts_at" = starts_at, "ends_at" = ends_at,
                             stringsAsFactors = FALSE)
    freq_table <- rbind(freq_table, freq_entry)
  }
  
  # freq_table_position <- data.frame("MHap" = position, "population"  = "", "sampleID" = "", 
  #                          "Hap1" = "", "Hap2" = "")
  
  freq_table_position <- data.frame("MHap" = MHapID$MHapID[i], "position" = position, stringsAsFactors = FALSE)
  position_table <- rbind(position_table, freq_table_position)
  # freq_table <- rbind(freq_table_position, freq_table)
  write.table(as.data.frame(freq_table), file = paste0(freq_output_dir, MHapID$MHapID[i], "_haplotype.txt"), quote= F, sep= "\t", row.names= F)
}

write.table(as.data.frame(position_table), file = paste0(freq_output_dir, "positions.txt"), quote= F, sep= "\t", row.names= F)






