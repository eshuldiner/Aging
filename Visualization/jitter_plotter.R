library(ggplot2)
library(dplyr)
library(scales)
library(RColorBrewer)
library(tidyverse)
library(ggrastr)


args <- commandArgs(TRUE)
input <- as.character(args[1])
inerts<-unlist(strsplit(as.character(args[2]),","))
root<- as.character(args[3])
project <-as.character(args[4])
parameter <-as.character(args[5])

'%!in%' <- function(x,y)!('%in%'(x,y))
blank_with_x_y_ticks =   theme(panel.background = element_rect(fill = "white"),
                               legend.position = "None",
                               panel.grid = element_blank(),
                               axis.line.y = element_line(linewidth = 0.5, colour = "black", linetype="solid"),
                               axis.line.x = element_line(linewidth = 0.5, colour = "black", linetype="solid"))

data<-read.csv(input, sep=",", header=TRUE)
data$sgID<-as.factor(data$sgID)
print(head(data))
print("inerts before")
print(inerts)
for (i in inerts){
  inerts<-c(inerts,strsplit(i,"_")[[1]][1])
}

inerts=unique(inerts)
print("inerts after")
print(inerts)

sgids_not_inert = levels(data$sgID)[levels(data$sgID) %!in% inerts]
print("sgids not inert")
print(sgids_not_inert)
sgids_ordered <-c(inerts, sgids_not_inert)
print("sgids ordered")
print(sgids_ordered)


data$sgID<-factor(data$sgID, levels=sgids_ordered)

data <- data %>% separate(sgID, "_", remove=FALSE,
   into=c("Gene","Guide_Num"))

data$Gene<-as.factor(data$Gene)
genes_not_inert = levels(data$Gene)[levels(data$Gene) %!in% inerts]
print("genes not inert")
print(genes_not_inert)

genes_inert <-levels(data$Gene)[levels(data$Gene) %in% inerts]
print("genes_inert is")
print(genes_inert)
genes_ordered<-c(genes_inert, genes_not_inert)
data$Gene<-factor(data$Gene, levels=genes_ordered)

#sorting out the color palette
inert_colors<- c(rep("gray", length(genes_inert)))
not_inert_colors <-c(brewer.pal(8,"Dark2"), brewer.pal(12,"Paired"), brewer.pal(8,"Set2"),brewer.pal(8,"Dark2"), brewer.pal(12,"Paired"), brewer.pal(8,"Set2"))
all_colors<-c(inert_colors, not_inert_colors)


##########THIS PLOT HAS ALL SAMPLES, SPLIT BY GENOTYPE (VIA FACETING)
#Note to self: using geom_jitter because geom_beeswarm takes too long for this much data
##############################################################################################################################
w =1+length(unique(data$sgID))*0.4

##############################################################################################################################
##############################################################################################################################
#####################These plots will have each treatment group split by sample
##############################################################################################################################
print(head(data))
treatments = unique(data$Treatment)
print("treatments are\n")
print(treatments)
for (t in treatments){
  data_t = data[data$Treatment==t,]
  print("sgids are ")
  print(unique(data_t$sgID))
  plot_t<-ggplot(data_t, aes(x=sgID, y=CellNum, col=Gene, size=CellNum)) +
  blank_with_x_y_ticks +
  theme(axis.line = element_line(colour = "black"),
    panel.border = element_blank(),
    panel.background = element_blank())+ 
  rasterise(geom_jitter()) + scale_size_area("",
                                  breaks = trans_breaks("log10", function(x) 10^x),
                                  labels = trans_format("log10", math_format(10^.x))) +
  scale_y_log10("Tumor Size (# Cells)",
                breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)))+
  guides(col=FALSE) + 
  theme(axis.text.x=element_text(angle=90, size=12)) +
  scale_color_manual(values =all_colors) +
  facet_grid(Sample~.)
  plot_name=paste0(root,"/",project, "/", parameter,"/figures/jitter_plots/",project, "_",parameter, "_", t, "_jitterplot_by_mouse.pdf")
  
  h = 2 + length(unique(data_t$Sample))*1
  
  cat(sprintf("Attempting to print a plot with dimensions %f by %f\n",w,h))
ggsave(plot_name,
       plot = plot_t,
       useDingbats=FALSE,
       width=w, height=h, limitsize=FALSE,
       units = "in",
       dpi = 500)

}

