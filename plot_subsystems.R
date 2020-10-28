#!/homes/clarkc/miniconda3/bin/Rscript

#parameter format and check parameter inputs
#plot_subsystems.R <subsystem_map.csv> <counts_file.txt|csv> <subsystem_level>
args = commandArgs(trailingOnly=TRUE)

if (length(args) != 5) {
    stop("Incorrect parameters: plot_subsystems.R <subsystem_map.txt> <counts_file.txt|csv> <metadata.txt>  <subsystem_level> <feature_count>")
}

subsystem.file = args[1]
counts.file = args[2]
metadata.file = args[3]
subsystem.level = args[4]
feature.count = args[5]

#check files and extensions
if (grepl("htseq",feature.count)) {
    count_sep = "\t"
} else if (grepl("stringtie",feature.count)) {
    count_sep = ","
} else {
    print("Error in plot_subsystems.R: can't determine counts file delimeter")
    print(counts.file)
    stop()
}

#load libraries quietly
library(ggplot2,quietly=TRUE)

#read tables
counts.mtx <- read.table(counts.file,sep=count_sep,header=T,row.names=1,stringsAsFactors=FALSE)
metadata <- read.table(metadata.file,sep="\t",header=T,stringsAsFactors=FALSE)
subsystem.map <- read.table(subsystem.file,sep="\t",header=T,stringsAsFactors=FALSE)
#Filter entries with no label and get the intersection of patric_ids
subsystem.map = subsystem.map[!grepl("NONE",subsystem.map[,2]),]
counts.mtx = counts.mtx[subsystem.map[,1],]

#Get all unique conditions
conditions = unique(metadata$Condition)

#Create a plot for each condition
for (c in conditions) {
    print(c)
    cond_label = paste(c,"_Expression",sep="")
    curr.mtx = counts.mtx[,subset(metadata,subset=Condition==c)$Sample]
    if (ncol(curr.mtx) > 1) {
        curr.mtx <- cbind(curr.mtx,rowMeans(curr.mtx))
    }
    colnames(curr.mtx)[-1] <- cond_label 
    plot.df <- data.frame(Genes=rownames(curr.mtx),Expression=log(curr.mtx[,c(cond_label)]+1),System=subsystem.map[,c(subsystem.level)])
    #box-whisker plot
    #boxplot <- ggplot(plot.df,aes(x=System,y=Expression))+geom_boxplot()+theme(axis.text.x=element_text(angle=90))
    #png(paste("test_boxplot",cond_label,".png",sep=""))
    #png(paste("test_boxplot",cond_label,".png",sep=""))
    #print(boxplot) 
    #violin plot
    vln_plot <- ggplot(plot.df,aes(x=System,y=Expression))+geom_violin()+theme(axis.text.x=element_text(angle=90))+ggtitle(cond_label)+ylab("Log-Expression")+xlab(subsystem.level)
    png(paste("Violin_",cond_label,"_",subsystem.level,"_mqc.png",sep="")) 
    print(vln_plot)
    dev.off()
}
