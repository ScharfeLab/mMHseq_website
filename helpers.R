#Helpers functions contains plotting functions for MH haplotype figures.
#can be combined into one function if time allowed
library(ggplot2)
library(grid)
library(stringr)
plot_microhaplotype_figure<-function(converted_snp_for_plot,
                                     mh_region,
                                     pop_info,
                                     each_sample_total_read,
                                     amplicon_total_read,
                                     each_sample_ID,
                                     spacing=45,
                                     spacing_2=45,
                                     fsize=5.5){
  haplotype_plot=data.frame("Haplotype"=c("Hap1","Hap2"),"start"=mh_region$START,"end"=mh_region$END)
  read_plot=data.frame(x1 = c(1.4,1.6), 
                       x2 = c(1.4,1.6), 
                       y1 = c(mh_region$END-250,mh_region$START), 
                       y2 = c(mh_region$END,mh_region$START+250))

  converted_snp_for_plot$status=factor(converted_snp_for_plot$status,
                                       levels=c("Ancestral",
                                                "Derived",
                                                "dbSNP",
                                                "Novel"))#set factor levels
  group.color=c("Ancestral"="#80796B99",
                "Derived"="red",
                "dbSNP"="#00A1D5FF",
                "Novel"="black")#color definition for each variant
  ######## Plot this MH region in this specific sample with geom_rect() in ggplot2
  #draw two empty rectangles as representations of haplotype
  g=ggplot(data=haplotype_plot)+geom_rect(aes(xmin=as.numeric(Haplotype)-0.2,
                                              xmax=as.numeric(Haplotype)+0.2,
                                              ymax=end,ymin=start),
                                              color="black",
                                              fill="white")+coord_flip()
  # #draw forward and reverse read
  #draw the variants if it is there
  g=g+geom_rect(data=converted_snp_for_plot,
                aes(xmin=as.numeric(haplotype)-0.2,xmax=as.numeric(haplotype)+0.2,
                    ymin=start-0.2,ymax=end+0.2,fill=status), color = NA)
  g=g+scale_fill_manual(name="Variant Category",values=group.color,
                                                    breaks=c("Derived",
                                                             "Ancestral",
                                                    "dbSNP",
                                                    "Novel"),drop=FALSE)
  g=g+geom_segment(aes(x=x1,y=y1,xend=x2,yend=y2),data=read_plot)
  g=g+theme(legend.position = "none",
            axis.text.y = element_blank(),
            axis.text.x = element_text(size=20),
            axis.ticks.y= element_blank(),
            axis.ticks.x = element_line(size=1, color="black"),
            axis.ticks.length = unit(0.5, "cm"),
            axis.title.y = element_blank(),
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(), 
            #panel.background = element_blank(),
            panel.background = element_blank(),
            # move x-axis label lower
            axis.title.x = element_text(margin = margin(t = 15)),
            text = element_text(size=20))
  # g = g+theme_bw()
  
  #Add x-axis and y-axis annoation
  g=g+ylab("Region(bp)")
  
  # g=g+scale_y_continuous(name = "Region(bp)",breaks=c(haplotype_plot$START,haplotype_plot$END),
  #                        labels=c(as.character(haplotype_plot$START,
  #                                              as.character(haplotype_plot$END))),limits = c(haplotype_plot$START,haplotype_plot$END))
  #Add sample info and QC information
  #g=g+annotate("text",x=2.5,y=haplotype_plot$start+18,label=mh_region$Amplicon,color="black",size=fsize)
  #g=g+annotate("text",x=2.5,y=haplotype_plot$start+3*spacing,label=paste0("Chr",converted_snp_for_plot$Chr[1]),color="black",size=fsize)
  #########g=g+annotate("text",x=2.5,y=haplotype_plot$start+2*spacing,label=converted_snp_for_plot$gene_or_locus[1],color="black",size=fsize)
  #g=g+annotate("text",x=2.5,y=haplotype_plot$start+4*spacing,label=each_sample_ID,color="black",size=fsize)
  #g=g+annotate("text",x=2.5,y=haplotype_plot$start+5*spacing,label=pop_info,color="black",size=fsize)
  
  # generate header
  label_header <- paste(mh_region$Amplicon, paste("Chr",converted_snp_for_plot$Chr[1], sep = ""), each_sample_ID, pop_info, sep = "       ")
  print(label_header)
  print(mh_region$Amplicon)
  g=g+annotate("text",x=1,y=haplotype_plot$start-9,label="H1",color="black",size=8)
  g=g+annotate("text",x=2,y=haplotype_plot$start-9,label="H2",color="black",size=8)
  g=g+annotate("text",x=2.5,y=haplotype_plot$start, hjust =0, label= label_header, color="black",size=fsize)
  print(converted_snp_for_plot$Chr[1])
  if (is.na(converted_snp_for_plot$Chr[1])) {
    g = g + annotate("text",x=2.5,y=haplotype_plot$end, hjust = 1, label= "Information not available", color="red",size=fsize)
  }
  print(each_sample_total_read)
  #if (each_sample_total_read<200000){
  #  g=g+annotate("text",x=2.5,y=haplotype_plot$start+6*spacing,label="QC-S",color="red",size=fsize)
  #}else{
  #  g=g+annotate("text",x=2.5,y=haplotype_plot$start+6*spacing,label="QC-S",color="blue",size=fsize)    
  #}
  #if (amplicon_total_read<100){
  #  g=g+annotate("text",x=2.5,y=haplotype_plot$start+6.7*spacing,label="QC-A",color="red",size=fsize)
  #}else{
  #  g=g+annotate("text",x=2.5,y=haplotype_plot$start+6.7*spacing,label="QC-A",color="blue",size=fsize)
  #}
  return(g)
}


# for download - simple version
# without forward, reversion read lines
# font and line thickness adjusted
plot_microhaplotype_figure_simple<-function(converted_snp_for_plot,
                                     mh_region,
                                     pop_info,
                                     each_sample_total_read,
                                     amplicon_total_read,
                                     each_sample_ID,
                                     spacing=45,
                                     spacing_2=45,
                                     fsize=5.5){
  
  extension <- 0.3
  haplotype_plot=data.frame("Haplotype"=c("Hap1","Hap2"),"start"=mh_region$START,"end"=mh_region$END)
  read_plot=data.frame(x1 = c(1.4,1.6), 
                       x2 = c(1.4,1.6), 
                       y1 = c(mh_region$END-250,mh_region$START), 
                       y2 = c(mh_region$END,mh_region$START+250))
  
  converted_snp_for_plot$status=factor(converted_snp_for_plot$status,
                                       levels=c("Ancestral",
                                                "Derived",
                                                "dbSNP",
                                                "Novel"))#set factor levels
  group.color=c("Ancestral"="#80796B99",
                "Derived"="red",
                "dbSNP"="#00A1D5FF",
                "Novel"="black")#color definition for each variant
  ######## Plot this MH region in this specific sample with geom_rect() in ggplot2
  #draw two empty rectangles as representations of haplotype
  g=ggplot(data=haplotype_plot)+geom_rect(aes(xmin=as.numeric(Haplotype)-extension,
                                              xmax=as.numeric(Haplotype)+extension,
                                              ymax=end,ymin=start),
                                          color="black",
                                          fill="white")+coord_flip()
  # #draw forward and reverse read
  #draw the variants if it is there
  g=g+geom_rect(data=converted_snp_for_plot,
                aes(xmin=as.numeric(haplotype)-extension,xmax=as.numeric(haplotype)+extension,
                    ymin=start-0.2,ymax=end+0.2,fill=status), color = NA)
  g=g+scale_fill_manual(name="Variant Category",values=group.color,
                        breaks=c("Derived",
                                 "Ancestral",
                                 "dbSNP",
                                 "Novel"),drop=FALSE)
  
  #g=g+geom_segment(aes(x=x1,y=y1,xend=x2,yend=y2),data=read_plot)
  
  g=g+theme(legend.position = "none",
            axis.text.y = element_blank(),
            axis.text.x = element_text(size=17),
            axis.ticks.y= element_blank(),
            axis.ticks.x = element_line(size=0.8, color="black"),
            axis.ticks.length = unit(0.5, "cm"),
            axis.title.y = element_blank(),
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(), 
            #panel.background = element_blank(),
            panel.background = element_blank(),
            # move x-axis label lower
            axis.title.x = element_text(margin = margin(t = 15)),
            text = element_text(size=20))
  # g = g+theme_bw()
  
  #Add x-axis and y-axis annoation
  g=g+ylab("Region(bp)")
  
  # generate header
  label_header <- paste0(str_pad(mh_region$Amplicon, 20, "right"))
  print(label_header)
  print(mh_region$Amplicon)
  g=g+annotate("text",x=1,y=haplotype_plot$start-8, label="H1",color="black",size=7)
  g=g+annotate("text",x=2,y=haplotype_plot$start-8, label="H2",color="black",size=7)
  g=g+annotate("text",x=2.6,y=haplotype_plot$start, hjust = 0, label= label_header, color="black", size=6)
  if (is.na(converted_snp_for_plot$Chr[1])) {
    g = g + annotate("text",x=2.6,y=haplotype_plot$end, hjust = 1, label= "Information not available", color="red",size=6)
  }
  
  
  return(g)
}

# for download - complete version
plot_microhaplotype_figure_complete<-function(converted_snp_for_plot,
                                     mh_region,
                                     pop_info,
                                     each_sample_total_read,
                                     amplicon_total_read,
                                     each_sample_ID,
                                     spacing=45,
                                     spacing_2=45,
                                     fsize=5.5){
  haplotype_plot=data.frame("Haplotype"=c("Hap1","Hap2"),"start"=mh_region$START,"end"=mh_region$END)
  read_plot=data.frame(x1 = c(1.4,1.6), 
                       x2 = c(1.4,1.6), 
                       y1 = c(mh_region$END-250,mh_region$START), 
                       y2 = c(mh_region$END,mh_region$START+250))
  
  converted_snp_for_plot$status=factor(converted_snp_for_plot$status,
                                       levels=c("Ancestral",
                                                "Derived",
                                                "dbSNP",
                                                "Novel"))#set factor levels
  group.color=c("Ancestral"="#80796B99",
                "Derived"="red",
                "dbSNP"="#00A1D5FF",
                "Novel"="black")#color definition for each variant
  ######## Plot this MH region in this specific sample with geom_rect() in ggplot2
  #draw two empty rectangles as representations of haplotype
  g=ggplot(data=haplotype_plot)+geom_rect(aes(xmin=as.numeric(Haplotype)-0.2,
                                              xmax=as.numeric(Haplotype)+0.2,
                                              ymax=end,ymin=start),
                                          color="black",
                                          fill="white")+coord_flip()
  # #draw forward and reverse read
  #draw the variants if it is there
  g=g+geom_rect(data=converted_snp_for_plot,
                aes(xmin=as.numeric(haplotype)-0.2,xmax=as.numeric(haplotype)+0.2,
                    ymin=start-0.2,ymax=end+0.2,fill=status), color = NA)
  g=g+scale_fill_manual(name="Variant Category",values=group.color,
                        breaks=c("Ancestral",
                                 "Derived",
                                 "dbSNP",
                                 "Novel"),drop=FALSE)
  g=g+geom_segment(aes(x=x1,y=y1,xend=x2,yend=y2),data=read_plot)
  g=g+theme(legend.position = "none",
            axis.text.y = element_blank(),
            axis.text.x = element_text(size=17),
            axis.ticks.y= element_blank(),
            axis.ticks.x = element_line(size=0.8, color="black"),
            axis.ticks.length = unit(0.5, "cm"),
            axis.title.y = element_blank(),
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(), 
            #panel.background = element_blank(),
            panel.background = element_blank(),
            # move x-axis label lower
            axis.title.x = element_text(margin = margin(t = 15)),
            text = element_text(size=20))
  # g = g+theme_bw()
  
  #Add x-axis and y-axis annoation
  g=g+ylab("Region(bp)")
  
  # generate header
  label_header <- paste(mh_region$Amplicon, paste("Chr",converted_snp_for_plot$Chr[1], sep = ""), each_sample_ID, pop_info, sep = "       ")
  print(label_header)
  print(mh_region$Amplicon)
  g=g+annotate("text",x=1,y=haplotype_plot$start-8,label="H1",color="black",size=7)
  g=g+annotate("text",x=2,y=haplotype_plot$start-8,label="H2",color="black",size=7)
  g=g+annotate("text",x=2.4,y=haplotype_plot$start, hjust = 0, label= label_header, color="black",size=fsize)
  if (is.na(converted_snp_for_plot$Chr[1])) {
    g = g + annotate("text",x=2.4,y=haplotype_plot$end, hjust = 1, label= "Information not available", color="red",size=fsize)
  }
  
  #print(each_sample_total_read)
  return(g)
}