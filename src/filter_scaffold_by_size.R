#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
# Args
# 1 = input
# 2 = number of bases to be considered as gaps (condition is >=)
# 3 = output

usageMessage=paste0(
    "\t 1 - Input bed file\n",
    "\t 2 - Number of bases to be considered as gaps\n",
    "\t     (Condition GapSize>=N)\n",
    "\t 3 - Output file\n",
    "\t     (TSV, format: <scaffold_id>     <gap_id>)"
)

if (length(args)<3) {
  stop(paste0("ERROR: 3 arguments must be supplied.\n", usageMessage), call.=FALSE)
} 

inputfile=args[1]
outputfile=args[3]
Nsize=as.numeric(args[2])

data=read.table(inputfile, header=F)
data$size=data$V3-data$V2
data$id=0
slist=unique(data$V1)
for (sc in slist){
    inds=which(data$V1==sc)
    data[inds,]$id=0:(length(inds)-1)
}
colnames(data)=c("scaffold","start","end", "size","id")
data2=data[data$size >= Nsize,]
write.table(file=outputfile,
    x=data2[,c("scaffold","id")],
    sep="\t",
    quote=F,
    row.names=F,
    col.names=F)

