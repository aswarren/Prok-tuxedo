#!/homes/clarkc/miniconda3/bin/Rscript

args = commandArgs(trailingOnly=TRUE)

numContrasts = length(args) - 4

if (numContrasts < 1) {
    stop("Not enough parameters: plot_subsystem_heatmaps.R <subsystem_file> <counts_file> <metadata_file> <output_prefix> <feature_count> <contrast 1> ... <contrast n>")
}

#check argument length
subsystem.file = args[1]
counts.file = args[2]
metadata.file = args[3]
subsystem.level = args[4]
feature.count = args[5]

#check file extensions
if (grepl("htseq",feature.count)) {
    count_sep = "\t"
} else if (grepl("stringtie",feature.count)) {
    count_sep = ","
} else {
    print("Error in plot_subsystem_heatmaps.R: can't determine counts file delimeter")
    print(counts.file)
    stop()
}

#load libraries quietly
library(ggplot2,quietly=TRUE)
library(reshape)

#TODO: create a heatmap with the log(gene_expression)
counts.mtx <- read.table(counts.file,sep=count_sep,header=T,row.names=1,stringsAsFactors=FALSE)
metadata <- read.table(metadata.file,sep="\t",header=T,stringsAsFactors=FALSE)
subsystem.map <- read.table(subsystem.file,sep="\t",header=T,stringsAsFactors=FALSE)
subsystem.map = subsystem.map[!grepl("NONE",subsystem.map[,2]),]
counts.mtx = counts.mtx[subsystem.map[,1],]

#unique conditions
conditions = unique(metadata$Condition)

#average-log matrix 
subsystem.df = data.frame(Subsystem=subsystem.map[,2])
#Average expression and log transform data
for (c in conditions) {
    print(c)
    #Subset data on current contrast
    curr.metadata = subset(metadata,subset=Condition==c)
    curr.count.mtx = counts.mtx[,curr.metadata$Sample]+1
    log.df = data.frame(Gene=rownames(curr.count.mtx),Log_Expression=log(rowMeans(curr.count.mtx)))
    log.df$Map <- mapply(function(gene) subset(subsystem.map,subset=Patric_ID==gene)[,2],gene=subsystem.map$Patric_ID)
    log.df = log.df[order(match(subsystem.map$Patric_ID,log.df$Gene)),]
    subsystem.df = cbind(subsystem.df,log.df$Log_Expression)
}
matrix_headers = c("Subsystem",conditions)
colnames(subsystem.df) = matrix_headers

#melt data frame
tmp.melt = melt(subsystem.df)
tmp.headers = c("Subsystem","Condition","Log_Expression")
colnames(tmp.melt) = tmp.headers
subsystem.melt = aggregate(tmp.melt$Log_Expression,by=list(tmp.melt$Condition,tmp.melt$Subsystem),FUN=sum)
heatmap.headers = c("Condition","Subsystem","Sum_Log_Expression")
colnames(subsystem.melt) <- heatmap.headers
#create heatmap
png("test_heatmap.png")
ggplot(subsystem.melt,aes(x=Subsystem,y=Condition,fill=Sum_Log_Expression))+geom_tile()+theme(axis.text.x=element_text(angle=90))
dev.off()
