#Helpers functions contains plotting functions for MH haplotype figures.
library(ggplot2)
plot_microhaplotype_figure<-function(converted_snp_for_plot,
                                     mh_region,
                                     pop_info,
                                     each_sample_total_read,
                                     amplicon_total_read,
                                     each_sample_ID,
                                     spacing=45,
                                     spacing_2=45,
                                     fsize=5){
  haplotype_plot=data.frame("Haplotype"=c("Hap1","Hap2"),"start"=mh_region$START,"end"=mh_region$END)
  read_plot=data.frame(x1 = c(1.4,1.6), 
                       x2 = c(1.4,1.6), 
                       y1 = c(mh_region$END-250,mh_region$START), 
                       y2 = c(mh_region$END,mh_region$START+250))

  converted_snp_for_plot$status=factor(converted_snp_for_plot$status,
                                       levels=c("KK(WT)",
                                                "KK(Variant)",
                                                "dbSNP",
                                                "Novel"))#set factor levels
  group.color=c("KK(WT)"="grey",
                "KK(Variant)"="red",
                "dbSNP"="blue",
                "Novel"="yellow")#color definition for each variant
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
                    ymin=start,ymax=end,fill=status))
  g=g+scale_fill_manual(name="Variant Category",values=group.color,
                                                    breaks=c("KK(Variant)",
                                                    "KK(WT)",
                                                    "dbSNP",
                                                    "Novel"),drop=FALSE)
  g=g+geom_segment(aes(x=x1,y=y1,xend=x2,yend=y2),data=read_plot)
  g=g+theme(legend.position = "none",
            axis.text.y = element_blank(),
            axis.ticks.y= element_blank(),
            axis.title.y = element_blank(),
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(), 
            panel.background = element_blank(),
            text = element_text(size=20))
  #Add x-axis and y-axis annoation
  g=g+ylab("Region(bp)")
  # g=g+scale_y_continuous(name = "Region(bp)",breaks=c(haplotype_plot$START,haplotype_plot$END),
  #                        labels=c(as.character(haplotype_plot$START,
  #                                              as.character(haplotype_plot$END))),limits = c(haplotype_plot$START,haplotype_plot$END))
  #Add sample info and QC information
  g=g+annotate("text",x=2.5,y=haplotype_plot$start+15,label=mh_region$Amplicon,color="black",size=fsize)
  g=g+annotate("text",x=2.5,y=haplotype_plot$start+2*spacing,label=paste0("Chr",converted_snp_for_plot$Chr[1]),color="black",size=fsize)
  #g=g+annotate("text",x=2.5,y=haplotype_plot$start+2*spacing,label=converted_snp_for_plot$gene_or_locus[1],color="black",size=fsize)
  g=g+annotate("text",x=2.5,y=haplotype_plot$start+3*spacing,label=each_sample_ID,color="black",size=fsize)
  g=g+annotate("text",x=2.5,y=haplotype_plot$start+4*spacing,label=pop_info,color="black",size=fsize)
  g=g+annotate("text",x=1,y=haplotype_plot$start-5,label="H1",color="black",size=fsize)
  g=g+annotate("text",x=2,y=haplotype_plot$start-5,label="H2",color="black",size=fsize)
  if (each_sample_total_read<200000){
    g=g+annotate("text",x=2.5,y=haplotype_plot$start+5*spacing,label="QC-S",color="red",size=fsize)
  }else{
    g=g+annotate("text",x=2.5,y=haplotype_plot$start+5*spacing,label="QC-S",color="blue",size=fsize)    
  }
  if (amplicon_total_read<100){
    g=g+annotate("text",x=2.5,y=haplotype_plot$start+6*spacing,label="QC-A",color="red",size=fsize)
  }else{
    g=g+annotate("text",x=2.5,y=haplotype_plot$start+6*spacing,label="QC-A",color="blue",size=fsize)
  }
  return(g)
}