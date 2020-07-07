library(ggplot2)
library(dplyr)
library(reshape2)
library(stringr)
library(gridExtra)
library(cowplot)
setwd("F:/data/ucsc/query/")
pafdir="syntenies/"
paf_target_query_filename=paste0(pafdir,"target.vs.query.paf")
query_sizes_filename=paste0(pafdir,"query.scaffold.sizes")
target_sizes_filename=paste0(pafdir,"target.scaffold.sizes")
paf_names=c(
  "qname","qlen","qstart","qend","strand","tname",
  "tlen","tstart","tend","res_matches","aln_block_len",
  "mapq","tag1","tag2","tag3","tag4","tag5","tag6")
###################################################################################################
# target sizes
target_sizes=read.table(target_sizes_filename, stringsAsFactors = F)
colnames(target_sizes)=c("chr","size")
# Define order of the target chrom/scaffold orders
target_scaffold_order=target_sizes$chr
target_sizes$scaffold=factor(target_sizes$chr, levels = target_scaffold_order)
target_sizes$cum_size=cumsum(as.numeric(target_sizes$size))
# query sizes
query_sizes=read.table(query_sizes_filename, stringsAsFactors = F)
colnames(query_sizes)=c("scaffold","size")
query_sizes=query_sizes[with(query_sizes, order(size, decreasing = T)),]
# organize commulative sum of sizes (for plotting)
query_sizes$cum_size=cumsum(as.numeric(query_sizes$size))
# Define order of the query chrom/scaffold orders
query_scaffold_order=query_sizes$scaffold
query_sizes$scaffold=factor(query_sizes$scaffold, levels = query_scaffold_order)
#------------------------------------------------------------------------------------------------------
aln=read.table(paf_target_query_filename, stringsAsFactors = F, fill=T)
colnames(aln)=paf_names
aln$x1=aln$qstart
aln$x2=aln$qend
aln$y1=aln$tstart
aln$y2=aln$tend
# Filtering alignments per quality
aln=subset(aln, mapq>40)
# need to factorize qname and tname to be able to get an order to tranlate x,y coordinates per chromosome
aln$tname=factor(aln$tname, levels = target_scaffold_order)
aln$qname=factor(aln$qname, levels = query_scaffold_order)
# update target coordinates to add previous chrom sizes
for(chrIndex in 2:nrow(target_sizes)){
  # first chr does not need to update coordinates
  chr=target_sizes[chrIndex,]$chr
  cumsize=target_sizes[chrIndex-1,]$cum_size
  indices=which(aln$tname==chr)
  if(length(indices)>0){
    aln[indices,]$y1=aln[indices,]$tstart + cumsize
    aln[indices,]$y2=aln[indices,]$tend + cumsize
  }
}
# updating query coordinates to add previous chrom sizes and modify segments according to strand
for(chrIndex in 2:nrow(query_sizes)){
  # first chr does not need to update coordinates
  chr=query_sizes[chrIndex,]$scaffold
  cumsize=query_sizes[chrIndex-1,]$cum_size
  indices=which(aln$qname==chr)
  if (length(indices)>0){
    aln[indices,]$x1=aln[indices,]$qstart + cumsize
    aln[indices,]$x2=aln[indices,]$qend + cumsize
  }
  indices2=which(aln$qname==chr & aln$strand == "-")
  if(length(indices2)>0){
    aln[indices2,]$x1=aln[indices2,]$qend + cumsize
    aln[indices2,]$x2=aln[indices2,]$qstart + cumsize
  }
}

query_labels=query_sizes$cum_size-(query_sizes$size/2)
target_labels=target_sizes$cum_size-(target_sizes$size/2)
aln$percent=aln$res_matches/aln$aln_block_len

aln_plot=ggplot(aln) +
  geom_segment(aes(x = x1, y = y1, xend = x2, yend = y2, colour = percent),size=2, data = aln) +
  geom_hline(data=target_sizes,   aes(yintercept=cum_size)) +
  geom_vline(data=query_sizes, aes(xintercept=cum_size)) +
  scale_colour_gradientn(colours = rainbow(5)) +
  scale_x_continuous(labels = query_sizes$scaffold, breaks=query_labels) + #,sec.axis = dup_axis())+
  scale_y_continuous(labels = target_sizes$chr, breaks=target_labels) + #,sec.axis = dup_axis())+
  # coord_cartesian(xlim=c(0,223616942), ylim=c(0, 248956422)) +
  theme(legend.position="None",   
        axis.ticks = element_blank(),
        axis.text.x=element_text(vjust=0.5, hjust=1,angle = 90, size = 10),
        axis.text.y=element_text(vjust=0.5, hjust=1, size = 15),
        axis.title = element_text(size=20))+
  labs(x="query",y="target", title="Synteny map", subtitle="target (target) vs. query (query)")

pdf("aln.pdf", height =12, width= 18)
print(aln_plot)
dev.off()