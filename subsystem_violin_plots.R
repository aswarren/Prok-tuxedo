#!/homes/clarkc/miniconda3/bin/Rscript

#parameter format and check parameter inputs
#subsystem_violin_plots.R <subsystem_map.csv> <counts_file.txt|csv> <subsystem_level>
args = commandArgs(trailingOnly=TRUE)

if (length(args) != 5) {
    stop("Incorrect parameters: subsystem_violin_plots.R <subsystem_map.txt> <counts_file.txt|csv> <metadata.txt>  <subsystem_level> <feature_count>")
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
    print("Error in subsystem_violin_plots.R: can't determine counts file delimeter")
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
#TODO: increasing the image width and height
#TODO: change color palette
#Create a plot for each condition
for (c in conditions) {
    print(c)
    #cond_label = paste(c,"_Expression",sep="")
    cond_label = c
    curr.mtx = counts.mtx[,subset(metadata,subset=Condition==c)$Sample]
    if (ncol(curr.mtx) > 1) {
        curr.mtx <- cbind(curr.mtx,rowMeans(curr.mtx))
    }
    colnames(curr.mtx)[-1] <- cond_label 
    plot.df <- data.frame(Genes=rownames(curr.mtx),Expression=log(curr.mtx[,c(cond_label)]+1),System=subsystem.map[,c(subsystem.level)])
    #violin plot
    cond_title = paste(cond_label," Log-Expression",sep="")
    vln_plot <- ggplot(plot.df,aes(x=System,y=Expression,fill=System))+geom_violin(trim=FALSE)+ggtitle(cond_title)+ylab("Log-Expression")+xlab("System")
    #add boxplot
    vln_plot = vln_plot + geom_boxplot(width=0.1,fill="white") 
    #remove x-axis text
    #+theme(axis.text.x=element_text(size=0,angle=45,vjust=-0.5))
    vln_plot = vln_plot + theme(axis.text.x = element_blank(),axis.ticks.x = element_blank())
    #plot image and send to png file
    vln_png = paste(cond_label,"_",subsystem.level,"_Expression","_mqc.png",sep="")
    png(vln_png,res=100,width=700) 
    print(vln_plot)
    dev.off()
    #ggsave(file=vln_png)
}
