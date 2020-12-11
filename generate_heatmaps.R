#!/homes/clarkc/miniconda3/bin/Rscript
args = commandArgs(trailingOnly=TRUE)

numContrasts = length(args) - 4

if (numContrasts < 1) {
    stop("Not enough parameters: generate_heatmaps.R <counts_file> <metadata_file> <heatmap_genes_file> <output_prefix> <feature_count> <specialty_genes")
}

#TODO: adjust so if specialty genes or other features aren't past in the script will still run
#TODO: as in, make the coloring of genes and such on the heatmap optional
#check argument length
counts.file = args[1]
metadata.file = args[2]
genes.file = args[3]
prefix = args[4]
feature.count = args[5]
specialty_genes.file = args[6]

#check file extensions
if (grepl("htseq",feature.count)) {
    count_sep = "\t"
} else if (grepl("stringtie",feature.count)) {
    count_sep = ","
} else {
    print("Error in generate_heatmaps.R: can't determine counts file delimeter")
    print(counts.file)
    stop()
}

#load libraries quietly
#library(ggplot2,quietly=TRUE)
library(ComplexHeatmap)
library(reshape2)

###TODO: add units to the heatmap (TPM, TPKM, etc)

#TODO: create a heatmap with the log(gene_expression)
counts.mtx <- read.table(counts.file,sep=count_sep,header=T,row.names=1,stringsAsFactors=FALSE)
metadata <- read.table(metadata.file,sep="\t",header=T,stringsAsFactors=FALSE)
genes.list <- read.table(genes.file,stringsAsFactors=FALSE)
colnames(genes.list) <- c("Genes")

###Process specialty genes
#TODO:total 9 different fields for the moment
sp.colors = c("#E41A1C","#377EB8","#4DAF4A","#984EA3","#FF7F00","#FFFF33","#A65628","#F781BF","#999999")
specialty.genes <- read.table(specialty_genes.file,header=T,sep="\t")
print(intersect(genes.list$Genes,specialty.genes$Patric_ID))
specialty.genes = subset(specialty.genes,subset=Patric_ID %in% genes.list$Genes)
print(specialty.genes$Patric_ID)
stop()
uniq_fields = unique(specialty.genes$SP_Field)
sp.colors = sp.colors[1:length(uniq_fields)]
specialty.genes$SP_Color = rep(0,length.out=length(specialty.genes$SP_Field))
for (i in 1:length(uniq_fields)) {
    specialty.genes[which(specialty.genes$SP_Field == uniq_fields[i]),][,3] = sp.colors[i] 
}
###

#conditions
conditions = unique(metadata$Condition)

#Average expression and log transform data
expression.df <- data.frame(Genes=genes.list$Genes)
matrix_headers = c("Genes")
sample_split = c()
for (c in conditions) {
    print(c)
    #Subset data on current contrast
    curr.metadata = subset(metadata,subset=Condition==c)
    sample_split = c(sample_split,rep(c,length(curr.metadata$Sample)))
    curr.count.mtx = log(counts.mtx[genes.list$Genes,curr.metadata$Sample]+1)
    matrix_headers = c(matrix_headers,curr.metadata$Sample)
    #log.df = data.frame(Gene=rownames(curr.count.mtx),Log_Expression=log(rowMeans(curr.count.mtx)))
    #expression.df = cbind(expression.df,log.df$Log_Expression)
    expression.df = cbind(expression.df,curr.count.mtx)
}
colnames(expression.df) = matrix_headers
rownames(expression.df) = expression.df$Genes
expression.mtx = as.matrix(expression.df[-1])
###TODO: setup the colors list
colors.list = rep("black",length.out=length(rownames(expression.mtx)))
names(colors.list) = rownames(expression.df)
colors.list[specialty.genes$Patric_ID] = specialty.genes$SP_Color
###
out_png = paste(prefix,"_heatmap_mqc.png",sep="")
#TODO: calculate image width and height based on number of samples
png(out_png,width=600,heigh=600)
#heatmap(expression.mtx,cexCol=1,Colv=NULL)
Heatmap(expression.mtx,name="log-expression",cluster_columns=FALSE,row_names_gp = gpar(fontsize=8,col=colors.list),column_names_gp = gpar(fontsize=8),column_names_rot = 45, column_split = sample_split, border=TRUE)
dev.off()
