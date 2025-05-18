library(ggplot2)
library(scales)
library(RColorBrewer)
library(tidyverse)
#library(reshape2)

blank_with_x_y_ticks =   theme(panel.background = element_rect(fill = "white"),
                               legend.position = "None",
                               panel.grid = element_blank(),
                               axis.line.y = element_line(linewidth = 0.5, colour = "black", linetype="solid"),
                               axis.line.x = element_line(linewidth = 0.5, colour = "black", linetype="solid"))

args <- commandArgs(TRUE)
input <- as.character(args[1])
project <- as.character(args[2])
parameter <- as.character(args[3])
root <- as.character(args[4])

data <- read.csv(input, header = T, sep = ",")
plot_n_str = paste0(root,project,"/",parameter,"/figures/QC_plots/",project,"_",parameter)#reminder: data is sgid_bc,sgid,bc,Nsamples,Nreads

print("header of data is\n")
print(head(data))

#plot 1: # of BC occuring in different numbers of samples
rejection_threshold = min(data[data$Rejected == "Rejected",]$N_samples)
sample_counts <- data %>% group_by(N_samples)%>%tally() %>%
mutate(Status=ifelse(N_samples>=rejection_threshold,"Rejected","Kept"))

p1<-ggplot(sample_counts,aes(x=N_samples,n,fill=Status))+
blank_with_x_y_ticks +
geom_bar(stat="identity")+
scale_y_log10("Number of barcodes",
              breaks = trans_breaks("log10", function(x) 10^x),
              labels = trans_format("log10", math_format(10^.x)))+
xlab("Number of samples barcode is found in")
theme(axis.text.x=element_text(angle=90, size=12))


plot_name=paste0(plot_n_str,"_NBC_by_NSample_full.png")

ggsave(plot_name,
       plot = p1,
       width = 7, height = 5, limitsize=FALSE,
       units = "in",
       dpi = 400)

####################################################
###### # of samples barcodes occur in my sgid ######
####################################################
get_threshold_per_sgid<-function(Rejected,N_samples){
  dat <- tibble(Rejected,N_samples)
  return(min(dat[dat$Rejected=="Rejected",]$N_samples))
}

get_thresholds <- data %>% group_by(sgID) %>% summarize(Threshold=get_threshold_per_sgid(Rejected, N_samples))

#Determining the # of rows and columns for the plot. This is somewhat complicated because different projects
# have very different #s of sgIDs
ncol=6
nsgid=length(unique(data$sgID))
if (nsgid%%ncol==0){
  nrow=nsgid/ncol
}else{
  nrow=1+round(nsgid/ncol)
}

width = ncol * 1.5
height = nrow * 1.5


sample_counts_by_sgid <- data %>% 
  group_by(N_samples,sgID) %>% 
  tally() %>% 
  filter(sgID!="Spi") %>% 
  left_join(get_thresholds) %>%
  mutate(Status = ifelse(N_samples>=Threshold,"Rejected","Kept")) 

p_sgid <- ggplot(sample_counts_by_sgid, aes(x=N_samples, y=n, fill=Status) )+
  blank_with_x_y_ticks +
  geom_bar(stat="identity")+
  scale_y_continuous(trans="log10", breaks = trans_breaks('log10', function(x) 10^x),
                     labels = trans_format('log10', math_format(10^.x)))+
  ylab("Number of barcodes") +
  xlab("Number of samples barcode is found in") +
  theme(axis.text.x=element_text(size=8))+
  facet_wrap(.~sgID, ncol=ncol)

plot_name=paste0(plot_n_str,"_Nsamples_by_sgID.png")

ggsave(plot_name,
       plot = p_sgid,
       width = width, height = height, limitsize=FALSE,
       units = "in",
       dpi = 400)
####################################
######plotting rejected data######
####################################
rejected_data <- data%>%filter(Rejected=="Rejected")
n_bc_removed<-rejected_data %>% group_by(sgID)%>%tally()


#plot 2: # of removed barcodes by sgid

p2 <- ggplot(n_bc_removed, aes(x=sgID,n))+
blank_with_x_y_ticks +
geom_bar(stat="identity")+
scale_y_continuous(label=comma) +
ylab("Number of barcodes removed") +
theme(axis.text.x=element_text(angle=90,size=12,hjust=0.95,vjust=0.2))

plot_name=paste0(plot_n_str,"_contamination_removal_nBC.png")

ggsave(plot_name,
       plot = p2,
       width = 7, height = 5, limitsize=FALSE,
       units = "in",
       dpi = 400)


n_reads_removed<-rejected_data %>% group_by(sgID)%>%summarize(totalReads=sum(N_reads))

#plot 2: # of removed reads removed by sgid

p3<-ggplot(n_reads_removed,aes(x=sgID,totalReads))+
blank_with_x_y_ticks +
scale_y_continuous(label=comma) +
geom_bar(stat="identity")+
ylab("Number of reads removed") +
theme(axis.text.x=element_text(angle=90,size=12,hjust=0.95,vjust=0.2))


plot_name=paste0(plot_n_str,"_contam_removal_nreads.png")

ggsave(plot_name,
       plot = p3,
       width = 7, height = 5, limitsize=FALSE,
       units = "in",
       dpi = 400)
