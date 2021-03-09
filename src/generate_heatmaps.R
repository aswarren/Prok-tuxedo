#!/homes/clarkc/miniconda3/bin/Rscript

args = commandArgs(trailingOnly=TRUE)

#check argument length
counts.file = args[1]
metadata.file = args[2]
genes.file = args[3]
prefix = args[4]
feature.count = args[5]
specialty_genes.file = args[6]
subsystem.file = args[7]

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
library(ComplexHeatmap,quietly=TRUE)
library(svglite)

###TODO: add units to the heatmap (TPM, TPKM, etc)

#create a heatmap with normalized expression 
counts.mtx <- read.table(counts.file,sep=count_sep,header=T,row.names=1,stringsAsFactors=FALSE)
metadata <- read.table(metadata.file,sep="\t",header=T,stringsAsFactors=FALSE)
genes.list <- read.table(genes.file,stringsAsFactors=FALSE)
colnames(genes.list) <- c("Genes")
#check for if specialty genes were found for this genome: if not don't include in heatmap
if (specialty_genes.file == "NONE") {
    specialty.genes <- NULL
} else {
    specialty.genes <- read.table(specialty_genes.file,header=T,sep="\t")
}
#check for if subsystem genes were found for this genome: if not don't include in heatmap
if (subsystem.file == "NONE") {
    subsystem.map <- NULL
} else {
    subsystem.map <- read.table(subsystem.file,sep="\t",header=T,stringsAsFactors=FALSE)
}

#normalize: put into the standard normal space 
counts.mtx = scale(counts.mtx) 

###Calculate picture widtth based on the number of samples
#heatmap width and offset are set in the draw() method at the bottom
png_width = nrow(metadata)*100 + 200   
svg_width = nrow(metadata) + ncol(metadata) 

###Process specialty genes if they were found
if (!is.null(specialty.genes)) 
{
    sp.labels = LETTERS 
    specialty.genes = subset(specialty.genes,subset=Patric_ID %in% genes.list$Genes)
    uniq_sp = unique(specialty.genes$SP_Field)
    sp.labels = sp.labels[1:length(uniq_sp)]
    specialty.genes$SP_Label = rep("",length.out=length(specialty.genes$SP_Field))
    for (i in 1:length(uniq_sp)) {
        specialty.genes[which(specialty.genes$SP_Field == uniq_sp[i]),][,3] = paste(" ",sp.labels[i],sep="") 
    }
    specialty.genes$Patric_ID = as.character(specialty.genes$Patric_ID)
}

###Process subsystem genes if they were found
#Filter entries with no label and get any included in the heatmap 
#TODO:total 9 different fields for the moment
if (!is.null(subsystem.map)) 
{
    sub.colors = c("#E41A1C","#377EB8","#4DAF4A","#984EA3","#FF7F00","#FFFF33","#A65628","#F781BF","#999999")
    subsystem.map = subsystem.map[!grepl("NONE",subsystem.map[,2]),]
    subsystem.map = subset(subsystem.map,subset=Patric_ID %in% genes.list$Genes)
    uniq_subsystems = unique(subsystem.map[,2]) 
    sub.colors = sub.colors[1:length(uniq_subsystems)] 
    subsystem.map$Sub_Color = rep(0,length.out=length(subsystem.map[,2]))
    for (i in 1:length(uniq_subsystems)) {
        subsystem.map[which(subsystem.map[,2] == uniq_subsystems[i]),][,3] = sub.colors[i]
    }
    subsystem.map$Patric_ID = as.character(subsystem.map$Patric_ID)
}

#Create Matrix
expression.df <- data.frame(Genes=genes.list$Genes)
matrix_headers = c("Genes")
sample_split = c()
conditions = unique(metadata$Condition)
for (c in conditions) {
    print(c)
    #Subset data on current contrast
    curr.metadata = subset(metadata,subset=Condition==c)
    sample_split = c(sample_split,rep(c,length(curr.metadata$Sample)))
    #curr.count.mtx = log(counts.mtx[genes.list$Genes,curr.metadata$Sample]+1)
    curr.count.mtx = counts.mtx[genes.list$Genes,curr.metadata$Sample]
    matrix_headers = c(matrix_headers,curr.metadata$Sample)
    expression.df = cbind(expression.df,curr.count.mtx)
}
colnames(expression.df) = matrix_headers
rownames(expression.df) = expression.df$Genes
expression.mtx = as.matrix(expression.df[-1])

###setup the subsystem colors list
colors.list = rep("black",length.out=length(rownames(expression.mtx)))
if (!is.null(subsystem.map))
{
    names(colors.list) = rownames(expression.df)
    colors.list[subsystem.map$Patric_ID] = subsystem.map$Sub_Color
    ###Create subsystem genes legend
    sub_legend = Legend(labels = uniq_subsystems, title = "Subsystem Genes", legend_gp = gpar(fill=sub.colors))
}
###setup the specialty genes label additions
if (!is.null(specialty.genes)) 
{
    for (usp in uniq_sp) {
        usp.genes = specialty.genes[which(specialty.genes$SP_Field == usp),]$Patric_ID
        exp.index = which(rownames(expression.mtx) %in% usp.genes)
        rownames(expression.mtx)[exp.index] = paste(rownames(expression.mtx)[exp.index],specialty.genes[which(specialty.genes$SP_Field == usp),]$SP_Label,sep="")
    }
    sp.index = which(rownames(expression.mtx) %in% specialty.genes$Patric_ID)
    ###Create specialty genes legend
    #eval(substitute()) exchanges the value of i for the integer at that step when evaluating the function call in draw()
    #Warning: not using eval(substitute()) causes sp.labels[i] uses the last value of i in the loop, which sets the icons in the legend ot the same value
    sp_list <- lapply(1:length(sp.labels),function(i){
        retval <- function(x,y,w,h) grid.text(eval(substitute(sp.labels[i],list(f=as.name(i)))),x,y)
        retval
    })
    sp_legend = Legend(labels = uniq_sp, title = "Specialty Genes", graphics = sp_list)
}

###create legend list from available heatmap legends
legend.list = list()
num_legends = 0
if (!is.null(subsystem.map))
{
    num_legends = num_legends + 1
    legend.list[[num_legends]] <- sub_legend
}
if (!is.null(specialty.genes)) 
{
    num_legends = num_legends + 1
    legend.list[[num_legends]] <- sp_legend
}

###Testing: take out "fig|" from gene names in heatmap so the legend doesn't overlap with the genes
for (i in 1:length(rownames(expression.mtx))) {
    split_gene = unlist(strsplit(rownames(expression.mtx)[i],"\\."))
    rownames(expression.mtx)[i] = paste("*",split_gene[length(split_gene)],sep=".")
}

###Create heatmap: SVG
#out_svg = paste("Normalized_Top_50_Differentially_Expressed_Genes_mqc.svg",sep="")
out_svg = paste("Normalized_Top_50_Differentially_Expressed_Genes.svg",sep="")
ht = Heatmap(expression.mtx,name="Z-score",cluster_columns=FALSE,row_names_max_width = unit(8, "cm"),show_row_dend = FALSE,row_names_gp = gpar(fontsize=8,col=colors.list),column_names_gp = gpar(fontsize=8),column_names_rot = 45, column_split = sample_split, border=TRUE, column_title = "Normalized Top 50 Differentially Expressed Genes")
#svg(out_svg,width=svg_width)
#svglite(out_svg,width=svg_width)
svglite(out_svg)
draw(ht,heatmap_legend_list=legend.list,padding=unit(c(1,1,1,15),"mm"))
dev.off()

###Name of output png: must end with "_mqc" in order to multiqc to recognize it
#out_png = paste("Normalized_Top_50_Differentially_Expressed_Genes_mqc.png",sep="")
###Create heatmap
#png(out_png,width=png_width,height=600)
#ht = Heatmap(expression.mtx,name="Normalized-Counts",cluster_columns=FALSE,row_names_gp = gpar(fontsize=8,col=colors.list),column_names_gp = gpar(fontsize=8),column_names_rot = 45, column_split = sample_split, border=TRUE)
#draw(ht,heatmap_legend_list=list(sub_legend,sp_legend),padding=unit(c(1,1,1,15),"mm"))
#dev.off()


