#!/home/cc8dm/miniconda3/envs/patric_env/bin/Rscript

#parameter format: 
#Rscript RunDESeq2.R <counts_file.txt> <metadata_file.txt> <contrast 1> <contrast 2> ... <contrast n>
#contrasts should be a csv pair and list all the contrasts to make in this dataset
args = commandArgs(trailingOnly=TRUE)

numContrasts = length(args) - 4

if (numContrasts < 1) {
    stop("Not enough parameters: RunDESeq2.R <counts_file> <metadata_file> <output_prefix> <feature_count> <contrast 1> ... <contrast n>")
}
counts.file = args[1]
metadata.file = args[2]
out_prefix = args[3]
feature_count = args[4]

#Check file extensions
if (grepl("htseq",feature_count)) {
    count_sep = "\t"
} else if (grepl("stringtie",feature_count)) {
    count_sep = "," 
} else { 
    print("Error in RunDESeq2.R: can't determine counts file delimeter")
    print(counts.file)
    stop() 
}

if (grepl("txt",metadata.file)) {
    meta_sep = "\t"
}
if (grepl("csv",metadata.file)) { 
    meta_sep = ","
}

#Differential expression library, load quietly so it doesn't completely fill the log output
suppressMessages(library(DESeq2,quietly=TRUE))
library(EnhancedVolcano,quietly=TRUE)

#Load counts table and metadata table and replace invalid characters
count.mtx <- read.table(counts.file,sep=count_sep,header=T,row.names=1,stringsAsFactors=FALSE)
rownames(count.mtx) = gsub("gene-","",rownames(count.mtx))
rownames(count.mtx) = gsub("rna-","",rownames(count.mtx))
colnames(count.mtx) = gsub("\\.","_",colnames(count.mtx))
metadata <- read.table(metadata.file,sep=meta_sep,header=T,row.names=1,stringsAsFactors=FALSE)
metadata$Condition = gsub("\\.","_",metadata$Condition)
rownames(metadata) = gsub("\\.","_",rownames(metadata))
#Reorder columns of counts table to match row order of metadata 
#TODO: May be able to omit this line
count.mtx = count.mtx[,rownames(metadata)]
#iterate over contrasts
#index 4 in args is where the contrasts currently start
for (i in 5:length(args)) 
{
    #Subset data on current contrast
    curr_contrast = unlist(strsplit(args[i],","))
    curr_contrast = gsub("-","_",curr_contrast)
    curr.metadata = subset(metadata,(subset=Condition==curr_contrast[1])|(subset=Condition==curr_contrast[2]))
    curr.count.mtx = count.mtx[,rownames(curr.metadata)] 
    #Remove all zero rows and add pseudocount to genes: Will not work for some samples otherwise
    curr.count.mtx = curr.count.mtx[rowSums(curr.count.mtx) != 0,]
    curr.count.mtx = curr.count.mtx + 1
    #Run standard DESeq2 pipeline
    print("running DESeq")
    dds <- DESeqDataSetFromMatrix(countData = curr.count.mtx, colData = curr.metadata, design = ~Condition)
    dds <- DESeq(dds)
    res <- results(dds,contrast=c("Condition",curr_contrast[1],curr_contrast[2]))
    #Try and remove any "gene-" or "rna-" prefixes
    #rownames(res) = gsub("gene-","",rownames(res))
    #rownames(res) = gsub("rna-","",rownames(res))
    #If invalid characters were present, put them back
    rownames(res) = gsub("___","-",rownames(res))
    colnames(res) = gsub("___","-",colnames(res))
    res = cbind(data.frame(Gene_Name=rownames(res)),res)
    #write to output file
    results_file = paste(curr_contrast[1],"_vs_",curr_contrast[2],".",feature_count,".",out_prefix,".deseq2",sep="")     
    write.table(res,file=results_file,sep="\t",quote=FALSE,row.names=FALSE)
    #Create volcano plot
    ev_image_name = paste(out_prefix,"_",curr_contrast[1],"_vs_",curr_contrast[2],"_mqc.png",sep="")
    png(ev_image_name)
    ev_png <- EnhancedVolcano(res,lab=rownames(res),x='log2FoldChange',y='padj',subtitle="",legendPosition="right")
    print(ev_png)
    dev.off()
} 
