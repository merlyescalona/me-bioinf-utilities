#!/usr/bin/env Rscript
# args = commandArgs(trailingOnly=TRUE)
##########################################################################################
suppressMessages(library(rlang));
suppressMessages(library(stringi));   
suppressMessages(library(stringr));   
suppressMessages(library(gdata));
suppressMessages(library(egg));       
suppressMessages(library(optparse));
suppressMessages(library(logger));
##########################################################################################
option_list = list(
    make_option(c("-d", "--discordant"), type="character", default=NULL, 
                help="Path of the BED file describing the discordant regions.",
                metavar="character"),
    make_option(c("-s", "--sizes"), type="character", default=NULL, 
                help="Path of the TSV file with the scaffold sizes.\n\t\tFormat: scaffold    size", 
                metavar="character"),
    make_option(c("-g", "--gaps"), type="character", default=NULL, 
                help="Path of the BED file describing the gaps in the scaffolds.\n\t\tFormat: scaffold   start   end", 
                metavar="character"),
    make_option(c("-t", "--threshold"), type="integer", default=50000, 
                help="Distance threshold to define closeness of a discordant region to a gap\n\t\t[default= %default]", 
                metavar="number"),
    make_option(c("-p", "--prefix"), type="character", default="out", 
                help="Output prefix  [default= %default]", 
                metavar="character"),
    make_option(c("-o", "--output"), type="character", default="~/", 
                help="Path to the output folder [default= %default]", 
                metavar="character")
); 
##########################################################################################
opt_parser = OptionParser(
    prog="discordant_regions_filter",
    usage="usage: %prog -d regions.bed -g gaps.bed -s scaffold.sizes\n\t\t[-t THRESHOLD -p PREFIX -o OUTPUT_FOLDER]",
    option_list=option_list, 
    description="",
    add_help_option=T,
    epilogue="(c) 2019 Merly Escalona <mmescalo@ucsc.edu>\n\t Paleogenomics Lab.\n\t Department of Biomolecular Engieering.\n\t University of California Santa Cruz");
opt = parse_args(opt_parser);
##########################################################################################
# Checking required options
if ( is.null(opt$discordant) | is.null(opt$sizes) | is.null(opt$gaps) ){
  print_help(opt_parser)
  stop("Required options -d/--discordant, -s/--sizes, -g/-gaps.", call.=F)
}
# Verifying existence of given arguments
if ( !file.exists(opt$discordant) ){
  stop("The BED file describing the discordant regions does not exist.", call.=F)
}
if ( !file.exists(opt$sizes) ){
  stop("The TSV file describing the sizes of the scaffolds does not exist.", call.=F)
}
if ( !file.exists(opt$gaps) ){
  stop("The BED file describing the gaps in the scaffolds does not exist.", call.=F)
}
##########################################################################################
# Genome information
DISTANCE_THRESHOLD=opt$threshold
sizesFile=opt$sizes
conflictsFile=opt$discordant
bedFile=opt$gaps
output=opt$output
prefix=opt$prefix
##########################################################################################
# Checking content of the files
log_info(paste0("Reading files..."))
gaps=try(read.table(bedFile, header=F, stringsAsFactors=F))
if ( !inherits(gaps, 'try-error')){
    colnames(gaps)=c("scaffold","start","end")
    gaps$start=as.numeric(gaps$start)
    gaps$end=as.numeric(gaps$end)
    gaps$size=gaps$end-gaps$start
    log_info(paste0("Number of gaps: ", nrow(gaps)))
}else{
    stop("There's a problem with the BED file describing the gaps in the scaffolds. Please verify and try again.", call.=F)
}
scafsizes=try(read.table(sizesFile, header=F,row.names=1))
if (!inherits(scafsizes, 'try-error')){
    log_info(paste0("Number of scaffolds: ", nrow(scafsizes)))
}else{
    stop("There's a problem with the TSV file describing the sizes of the scaffolds. Please verify and try again.", call.=F)
}
conflicts1=try(read.table(conflictsFile,stringsAsFactors=F,header=F, fill=T))
if ( !inherits(conflicts1, 'try-error')){
    colnames(conflicts1)=c("scaffold","start","end","type")
    conflicts1$start=as.numeric(as.character(conflicts1$start))
    conflicts1$end=as.numeric(as.character(conflicts1$end))
    log_info(paste0("Number of regions: ", nrow(conflicts1)))
    conflicts1=unique(conflicts1)
    log_info(paste0("Number of (unique) regions: ", nrow(conflicts1)))
    conflicts1[conflicts1$start == -1,]$start=conflicts1[conflicts1$start == -1,]$end
    conflicts1[conflicts1$end == -1,]$end=conflicts1[conflicts1$end == -1,]$start
    conflicts1=unique(conflicts1)
    log_info(paste0("Number of (unique - 2) regions: ", nrow(conflicts1)))
    conflicts1$size=conflicts1$end-conflicts1$start
}else{
    stop("There's a problem with the BED file describing the discordant regions. Please verify and try again.", call.=F)
}
##########################################################################################
# Conflitcs expansion
log_info(paste0("Widening regions..."))
conflicts1$lineStart=apply(conflicts1[,c(1,2)], 1, function(x){ stri_join_list(list(sort(unique(x))), sep=";")})
conflicts1$lineEnd=apply(conflicts1[,c(1,3)], 1, function(x){ stri_join_list(list(sort(unique(x))), sep=";")})
conflicts1$mod=F
tableLinesStart1=table(conflicts1$lineStart)
tableLinesEnd1=table(conflicts1$lineEnd)
conflicts2=duplicate(conflicts1)
# print(head(conflicts1))
# stop("conflicts1",call.=F)
for(subsetEnd in names(tableLinesEnd1[tableLinesEnd1>1])){
  allstart=conflicts1[conflicts1$lineEnd==subsetEnd,]$start
  allends=conflicts1[conflicts1$lineEnd==subsetEnd,]$end
  conflicts2[conflicts1$lineEnd==subsetEnd & conflicts1$start %in% allstart,]$start=min(allstart,na.rm=T)
  conflicts2[conflicts1$lineEnd==subsetEnd & conflicts1$start %in% allstart,]$mod=T
}
conflicts2$lineEnd=apply(conflicts2[,c(1,3)], 1, function(x){ stri_join_list(list(sort(unique(x))), sep=";")})
conflicts2$lineStart=apply(conflicts2[,c(1,2)], 1, function(x){ stri_join_list(list(sort(unique(x))), sep=";")})

log_info(paste0("Number of regions (after widening - side 1): ", nrow(conflicts2)))
tableLinesStart2=table(conflicts2$lineStart)
tableLinesEnd2=table(conflicts2$lineEnd)
conflicts3=duplicate(conflicts2)
for(subsetStart in names(tableLinesStart2[tableLinesStart2>1])){
  allstart=conflicts2[conflicts2$lineStart==subsetStart,]$start
  allends=conflicts2[conflicts2$lineStart==subsetStart,]$end
  indices=conflicts2$lineStart==subsetStart & conflicts2$end %in% allends
  conflicts3[indices,]$end=max(allends,na.rm=T)
  conflicts3[indices,]$mod=T
  if (length(unique(list(conflicts3[indices,]$type)[[1]]))>1){
    conflicts3[indices,]$type=stri_join_list(list(sort(unique(conflicts3[indices,]$type))), sep=";")
  }else{
    conflicts3[indices,]$type=unique(conflicts3[indices,]$type)
  }
}
conflicts3$lineStart=apply(conflicts3[,c(1,2)], 1, function(x){ stri_join_list(list(sort(unique(x))), sep=";")})
conflicts3$lineEnd=apply(conflicts3[,c(1,3)], 1, function(x){ stri_join_list(list(sort(unique(x))), sep=";")})
log_info(paste0("Number of regions (after widening - side 2): ", nrow(conflicts3)))
##########################################################################################
# Expanding conflicts
# Generating double conflict records for the onew that have size > 1
log_info(paste0("Expanding conflicts..."))
conflicts4=subset(conflicts3, conflicts3$type != "conflict")
singleconflicts=subset(conflicts3, conflicts3$type == "conflict" & conflicts3$size ==1 )
onlyconfs=subset(conflicts3, conflicts3$type == "conflict" &  conflicts3$size > 1 )
newrecords=c()
for (index in 1:nrow(onlyconfs)){
  newrecord1=onlyconfs[index,]
  newrecord2=onlyconfs[index,]
  newrecord1$end = newrecord1$start + 1 
  newrecord2$start = newrecord2$end - 1

  newrecord1$lineStart=stri_join_list(list(sort(unique(newrecord1[1,c("scaffold", "start")]))), sep=";")
  newrecord2$lineStart=stri_join_list(list(sort(unique(newrecord2[1,c("scaffold", "start")]))), sep=";")
  newrecord1$lineEnd=stri_join_list(list(sort(unique(newrecord1[1,c("scaffold", "end")]))), sep=";")
  newrecord2$lineEnd=stri_join_list(list(sort(unique(newrecord2[1,c("scaffold", "end")]))), sep=";")
  newrecords=rbind(newrecords, newrecord1,newrecord2)
  
}
newrecords$size=newrecords$end-newrecords$start
newrecords$mod=T
conflicts=rbind(unique(conflicts4), singleconflicts, newrecords)
conflicts=unique(conflicts)
log_info(paste0("Number of regions (after expansion): ", nrow(conflicts)))
##########################################################################################
log_info(paste0("Getting full data frame..."))
out=c()
for (indexD in 1:nrow(conflicts)){
  sc=conflicts[indexD,]$scaffold
  gapsSC=subset(gaps, gaps$scaffold==sc)
  start=which(gapsSC$start<=conflicts[indexD,]$start)
  end=which(gapsSC$end>=conflicts[indexD,]$end)
  startBlock=data.frame(scaffold="",start=0,end=0,size=0)
  endBlock=data.frame(scaffold="",start=0,end=0,size=0)
  if (length(start)>0){startBlock=gapsSC[max(start,na.rm=T),]}
  if (length(end)>0){endBlock=gapsSC[min(end,na.rm=T),]}
  out=rbind(out, cbind(
      conflicts[indexD,], 
      startBlock, 
      endBlock)
      )
}
out=data.frame(out)
out=out[,c(1:5,8,10:12,14:16)]
colnames(out)=c("scaffold","disc_start","disc_end","type","disc_size","mod",
                "right_gap_start", "right_gap_end", "right_gap_size",
                "left_gap_start","left_gap_end","left_gap_size")
##########################################################################################
# Checking for tags multiple tims
out$type=str_replace(out$type, "inversion;inversion", "inversion")
out$type=str_replace(out$type, "inversion_paired;inversion_paired", "inversion_paired")
out$type=str_replace(out$type, "inversion_partial;inversion_partial", "inversion_partial")
out$type=str_replace(out$type, "conflict;conflict", "conflict")
##########################################################################################

# log_info(str(out))
# stop("bu",call.=F)
# log_info(paste0("nrow(gapsSC): ", nrow(gapsSC)," length(start): ", length(start)," length(end):", length(end)))
out$dist_right_start_disc_start=abs(out$right_gap_start-out$disc_start)
out$dist_right_start_disc_end=abs(out$right_gap_start-out$disc_end)
out$dist_right_end_disc_start=abs(out$right_gap_end-out$disc_start)
out$dist_right_end_disc_end=abs(out$right_gap_end-out$disc_end)
out$dist_right=apply(
    cbind(
        out$dist_right_start_disc_start,out$dist_right_start_disc_end,
        out$dist_right_end_disc_start,out$dist_right_end_disc_end),
        1, min, na.rm = T)
out$dist_left_start_disc_start=abs(out$left_gap_start-out$disc_start)
out$dist_left_start_disc_end=abs(out$left_gap_start-out$disc_end)
out$dist_left_end_disc_start=abs(out$left_gap_end-out$disc_start)
out$dist_left_end_disc_end=abs(out$left_gap_end-out$disc_end)
out$dist_left=apply(
    cbind(
        out$dist_left_start_disc_start,out$dist_left_start_disc_end,
        out$dist_left_end_disc_start,out$dist_left_end_disc_end), 
        1,min, na.rm = T)
out$assembly_right=""
out$assembly_left=""
out$assembly=""

out[!is.na(out$right_gap_size) & out$right_gap_size==100 ,]$assembly_right="hirise"
out[!is.na(out$right_gap_size) & out$right_gap_size!=100 ,]$assembly_right="bionano"
out[!is.na(out$left_gap_size) & out$left_gap_size==100 ,]$assembly_left="hirise"
out[!is.na(out$left_gap_size) & out$left_gap_size!=100 ,]$assembly_left="bionano"
out[out$assembly_right=="hirise" & out$assembly_left=="hirise",]$assembly="hirise"
out[out$assembly_right=="hirise" & out$assembly_left!="hirise",]$assembly="bionano_hirise"
out[out$assembly_right!="hirise" & out$assembly_left=="hirise",]$assembly="hirise_bionano"
out[out$assembly_right!="hirise" & out$assembly_left!="hirise",]$assembly="bionano"
out$within_distance=F
##########################################################################################
log_info(paste0("Filtering..."))
out[out$dist_right < DISTANCE_THRESHOLD | out$dist_left < DISTANCE_THRESHOLD  ,]$within_distance=T
out_curated=out[, c(
  c("scaffold","disc_start","disc_end","disc_size","type","mod",
    "assembly","assembly_right","assembly_left","dist_right","dist_left",
    "right_gap_start", "right_gap_end","left_gap_start","left_gap_end", "within_distance")
)]
breaks_right=out[!is.infinite(out$dist_right) & out$dist_right< DISTANCE_THRESHOLD ,
    c("scaffold","right_gap_start", "right_gap_end", "type", "assembly_right")]
breaks_left=out[!is.infinite(out$dist_left) &  out$dist_left< DISTANCE_THRESHOLD,
    c("scaffold","left_gap_start", "left_gap_end", "type", "assembly_left")]
colnames(breaks_right)=c("scaffold","start","end","type", "assembly")
colnames(breaks_left)=c("scaffold","start","end","type", "assembly")
breaks=unique(rbind(breaks_right, breaks_left))
rownames(breaks)=1:nrow(breaks)
breaks$action="-"
breaks$gapid=paste0(breaks$scaffold,breaks$start, breaks$end, sep="|")
splitbreak=split(breaks, breaks$gapid)
##########################################################################################
tostore=c()
for (gapid in names(splitbreak)){
  newrec=splitbreak[gapid][[1]][1,]
  if( length(unique(splitbreak[gapid][[1]]))>1){
    uniquetype=sort(unique(splitbreak[gapid][[1]]$type))
    newrec$type=stri_join_list(list(uniquetype), sep=";")
  }
  if( length(unique(splitbreak[gapid][[1]]))>1){
    uniqueassembly=unique(splitbreak[gapid][[1]]$assembly)
    uniqueassembly=sort(uniqueassembly)
    newrec$assembly=stri_join_list(list(uniqueassembly), sep=";")
  }
  tostore=rbind(tostore,newrec)
}

tostore=unique(tostore)
dftostore=data.frame(
  type=tostore$type,
  assembly=tostore$assembly
)
tostore$type=as.character(tostore$type)
tostore$assembly=as.character(tostore$assembly)
##########################################################################################
log_info(paste0("Preparing data for output..."))
tostore=tostore[,c("scaffold","start","end","action","type","assembly")]
tostore=tostore[tostore$assembly == "hirise",]
brconflicts=subset(tostore, tostore$type=="conflict")
br9=  subset(tostore, tostore$type=="conflict;inversion")
br3=  subset(tostore, tostore$type=="conflict;inversion_partial")
br2=  subset(tostore, tostore$type=="conflict;inversion_paired")
br11= subset(tostore, tostore$type=="conflict;inversion;inversion_paired")
br1=  subset(tostore, tostore$type=="conflict;inversion;inversion_partial")
br12= subset(tostore, tostore$type=="conflict;inversion_paired;inversion_partial")
br14= subset(tostore, tostore$type=="conflict;inversion;inversion_paired;inversion_partial")
br4=  subset(tostore, tostore$type=="inversion")
br6=  subset(tostore, tostore$type=="inversion_paired")
br8=  subset(tostore, tostore$type=="inversion_partial")
br10= subset(tostore, tostore$type=="inversion;inversion_paired")
br5=  subset(tostore, tostore$type=="inversion;inversion_partial")
br7=  subset(tostore, tostore$type=="inversion_paired;inversion_partial")
br13= subset(tostore, tostore$type=="inversion;inversion_paired;inversion_partial")
log_info("################################################################################")
log_info("Summary: ")
log_info("################################################################################")
log_info(paste0("Distribution of potential breaks according to assembly"))
log_info("--------------------------------------------------------------------------------")
log_info(paste0("Total: ", nrow(breaks)))
log_info(paste0("Total: ", nrow(unique(breaks[,1:3]))))
table(breaks[,c("type","assembly")])
log_info("--------------------------------------------------------------------------------")
table(breaks[intersect(rownames(unique(breaks[,1:3])), which(breaks$assembly=="hirise")), c("type","assembly")])
log_info("--------------------------------------------------------------------------------")
log_info(paste0("Distribution of (unique) potential breaks from HiRise:"))
log_info("--------------------------------------------------------------------------------")
log_info(paste0("Potential breakpoints (HiRise): ", nrow(tostore)))
log_info(paste0("Potential breakpoints (conflict): ", nrow(brconflicts)))
log_info(paste0("Potential breakpoints (conflict;inversion): ", nrow(br9)))
log_info(paste0("Potential breakpoints (conflict;inversion_paired): ", nrow(br2)))
log_info(paste0("Potential breakpoints (conflict;inversion_partial): ", nrow(br3)))
log_info(paste0("Potential breakpoints (conflict;inversion;inversion_paired): ", nrow(br11)))
log_info(paste0("Potential breakpoints (conflict;inversion;inversion_partial): ", nrow(br1)))
log_info(paste0("Potential breakpoints (conflict;inversion_paired;inversion_partial): ", nrow(br12)))
log_info(paste0("Potential breakpoints (conflict;inversion;inversion_paired;inversion_partial): ", nrow(br14)))
log_info(paste0("Potential breakpoints (inversion): ", nrow(br4)))
log_info(paste0("Potential breakpoints (inversion_paired): ", nrow(br6)))
log_info(paste0("Potential breakpoints (inversion_partial): ", nrow(br8)))
log_info(paste0("Potential breakpoints (inversion;inversion_partial): ", nrow(br5)))
log_info(paste0("Potential breakpoints (inversion;inversion_paired): ", nrow(br10)))
log_info(paste0("Potential breakpoints (inversion_paired;inversion_partial): ", nrow(br7)))
log_info(paste0("Potential breakpoints (inversion;inversion_paired;inversion_partial): ", nrow(br13)))
log_info("################################################################################")

##########################################################################################
tostore=apply(tostore,2,as.character)
brconflicts=apply(brconflicts,2,as.character)
br1=apply(br1,2,as.character)
br2=apply(br2,2,as.character)
br3=apply(br3,2,as.character)
br4=apply(br4,2,as.character)
br5=apply(br5,2,as.character)
br6=apply(br6,2,as.character)
br7=apply(br7,2,as.character)
br8=apply(br8,2,as.character)
br9=apply(br9,2,as.character)
br10=apply(br10,2,as.character)
br11=apply(br11,2,as.character)
br12=apply(br12,2,as.character)
br13=apply(br13,2,as.character)
br14=apply(br14,2,as.character)
##########################################################################################
log_info(paste0("Writing output files..."))
write.table(tostore[,1:4],file = paste0(output, "/",prefix,".breakpoints.bed"), row.names = F, col.names=F, quote = F, sep = "\t")
write.table(breaks[intersect(rownames(unique(breaks[,1:3])), which(breaks$assembly=="hirise")),],file = paste0(output, "/",prefix,".breaks"), row.names = F, col.names=F, quote = F, sep = "\t")
write.table(tostore,file = paste0(output, "/",prefix,".hirise.auto.bed"), row.names = F, col.names=F, quote = F, sep = "\t")
write.table(brconflicts,file = paste0(output, "/",prefix,".hirise.auto.conflicts.bed"), row.names = F,col.names=F, quote = F, sep = "\t")

write.table(br9,  row.names = F,col.names=F, quote = F, sep = "\t", file = paste0(output, "/", prefix,".hirise.auto.conflict_inv.bed"))
write.table(br2,  row.names = F,col.names=F, quote = F, sep = "\t", file = paste0(output, "/", prefix,".hirise.auto.conflicts_invpair.bed"))
write.table(br3,  row.names = F,col.names=F, quote = F, sep = "\t", file = paste0(output, "/", prefix,".hirise.auto.conflicts_invpart.bed"))
write.table(br12, row.names = F,col.names=F, quote = F, sep = "\t", file = paste0(output, "/", prefix,".hirise.auto.conflicts_invpair_invpart.bed"))
write.table(br11, row.names = F,col.names=F, quote = F, sep = "\t", file = paste0(output, "/", prefix,".hirise.auto.conflicts_inv_invpair.bed"))
write.table(br1,  row.names = F,col.names=F, quote = F, sep = "\t", file = paste0(output, "/", prefix,".hirise.auto.conflicts_inv_invpart.bed"))
write.table(br14, row.names = F,col.names=F, quote = F, sep = "\t", file = paste0(output, "/", prefix,".hirise.auto.conflicts_inv_invpair_invpart.bed"))
write.table(br4,  row.names = F,col.names=F, quote = F, sep = "\t", file = paste0(output, "/", prefix,".hirise.auto.inv.bed"))
write.table(br6,  row.names = F,col.names=F, quote = F, sep = "\t", file = paste0(output, "/", prefix,".hirise.auto.invpair.bed"))
write.table(br8,  row.names = F,col.names=F, quote = F, sep = "\t", file = paste0(output, "/", prefix,".hirise.auto.invpart.bed"))
write.table(br5,  row.names = F,col.names=F, quote = F, sep = "\t", file = paste0(output, "/", prefix,".hirise.auto.inv_invpart.bed"))
write.table(br10, row.names = F,col.names=F, quote = F, sep = "\t", file = paste0(output, "/", prefix,".hirise.auto.inv_invpair.bed"))
write.table(br7,  row.names = F,col.names=F, quote = F, sep = "\t", file = paste0(output, "/", prefix,".hirise.auto.invpair_invpart.bed"))
write.table(br13, row.names = F,col.names=F, quote = F, sep = "\t", file = paste0(output, "/", prefix,".hirise.auto.inv_invpair_invpart.bed"))
##########################################################################################
log_info(paste0("End..."))