#!/home/cc8dm/miniconda3/bin/Rscript

#parameter format: 
#Rscript RunDESeq2.R <counts_file.txt> <metadata_file.txt> <contrast 1> <contrast 2> ... <contrast n>
#contrasts should be a csv pair and list all the contrasts to make in this dataset
args = commandArgs(trailingOnly=TRUE)

numContrasts = length(args) - 3

if (numContrasts < 1) {
    stop("Not enough parameters: Rscript RunDESeq2.R <counts_file.txt> <metadata_file.txt> <output_prefix> <contrast 1> ... <contrast n>")
}

library(DESeq2)
library(tools)

counts.file = args[1]
metadata.file = args[2]
out_prefix = args[3]

#Check file extensions
if (file_ext(counts.file) == "txt") {
    count_sep = "\t"
}
if (file_ext(counts.file) == "csv") { 
    count_sep = ","
}
if (file_ext(metadata.file) == "txt") {
    meta_sep = "\t"
}
if (file_ext(metadata.file) == "csv") { 
    meta_sep = ","
}
#Load counts table and metadata table and replace invalid characters
count.mtx <- read.table(counts.file,sep=count_sep,header=T,row.names=1,stringsAsFactors=FALSE)
rownames(count.mtx) = gsub("-","___",rownames(count.mtx))
colnames(count.mtx) = gsub("-","___",colnames(count.mtx))
metadata <- read.table(metadata.file,sep=meta_sep,header=T,row.names=1,stringsAsFactors=FALSE)
metadata$Condition = gsub("-","___",metadata$Condition)
rownames(metadata) = gsub("-","___",rownames(metadata))

#Reorder columns of counts table to match row order of metadata 
#Note: May be able to omit this line
count.mtx = count.mtx[,rownames(metadata)]

#create a list, size of number of comparisons, that will hold the diff_exp data
#diff_exp_list <- rep(0,numContrasts)
#contrast_list <- rep(0,numContrasts)
#gmx_frame = data.frame(matrix(ncol=0,nrow=length(rownames(count.mtx))))
#rownames(gmx_frame) <- rownames(count.mtx)
#iterate over contrasts
for (i in 4:length(args)) 
{
    #Subset data on current contrast
    curr_contrast = unlist(strsplit(args[i],","))
    curr_contrast = gsub("-","_",curr_contrast)
    curr.metadata = subset(metadata,(subset=Condition==curr_contrast[1])|(subset=Condition==curr_contrast[2]))
    curr.count.mtx = count.mtx[,rownames(curr.metadata)] 
    #Remove all zero rows and add pseudocount to genes
    curr.count.mtx = curr.count.mtx[rowSums(curr.count.mtx) != 0,]
    curr.count.mtx = curr.count.mtx + 1
    #Run standard DESeq2 pipeline
    dds <- DESeqDataSetFromMatrix(countData = curr.count.mtx, colData = curr.metadata, design = ~Condition)
    dds <- DESeq(dds)
    res <- results(dds,contrast=c("Condition",curr_contrast[1],curr_contrast[2]))
    rownames(res) = gsub("___","-",rownames(res))
    colnames(res) = gsub("___","-",colnames(res))
    #store information for gmx frame
    #res = res[rownames(gmx_frame),]
    #gmx_frame = cbind(gmx_frame,res)
    #contrast_list[i-3] = paste(curr_contrast[1],curr_contrast[2],sep="_vs_")
    #write to output file
    results_file = paste(out_prefix,"_",curr_contrast[1],"_vs_",curr_contrast[2],".txt",sep="")     
    write.table(res,file=results_file,sep="\t",quote=FALSE)
} 

#Write out gmx file: genes as rownames, comparisons as colnames, logfoldchange as values
#gmx_file = "gene_exp.gmx"
#colnames(gmx_frame) <- contrast_list
#write.table(gmx_frame,file=gmx_file,sep="\t",quote=FALSE)
