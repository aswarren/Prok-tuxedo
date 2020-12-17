#!/homes/clarkc/miniconda3/bin/Rscript

###Functions
#https://stackoverflow.com/questions/12041042/how-to-plot-just-the-legends-in-ggplot2
g_legend <- function(a.gplot) {
    tmp <- ggplot_gtable(ggplot_build(a.gplot))
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
    legend <- tmp$grobs[[leg]]
    legend
}

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
library(gridExtra,quietly=TRUE)
library(reshape2,quietly=TRUE)

#read tables
counts.mtx <- read.table(counts.file,sep=count_sep,header=T,row.names=1,stringsAsFactors=FALSE)
metadata <- read.table(metadata.file,sep="\t",header=T,stringsAsFactors=FALSE)
subsystem.map <- read.table(subsystem.file,sep="\t",header=T,stringsAsFactors=FALSE)

#Filter entries with no system label and get the intersection of patric_ids
subsystem.map = subsystem.map[!grepl("NONE",subsystem.map[,2]),]
counts.mtx = counts.mtx[subsystem.map[,1],]

#Testing: min and max values
log_min = log(min(counts.mtx)+1)
log_max = log(max(counts.mtx))

#Get all unique subsystems and conditions
conditions = unique(metadata$Condition)
subsystems = unique(subsystem.map[,2])  

#Calculate image width and height
num_columns <- ceiling(sqrt(length(subsystems)))
num_samples <- ncol(counts.mtx)
num_rows <- ceiling(length(subsystems)/num_columns)
png_width = (num_columns + num_samples)*100
png_height = num_rows*200 

#Get all unique features/subsystems
legend <- NULL 
plot_list = vector("list",length(subsystems)+1)
for (i in 1:length(subsystems)) {
    curr.system = subsystems[i] 
    curr.mtx = counts.mtx[rownames(counts.mtx) %in% subsystem.map[which(subsystem.map[,2] == curr.system),1],] 
    curr.mtx = data.frame(curr.mtx)
    curr.mtx$Genes <- rownames(curr.mtx)
    melt.df = melt(curr.mtx,id.vars=c("Genes"),measure.vars=colnames(curr.mtx)[-c(length(colnames(curr.mtx)))]) 
    colnames(melt.df) <- c("Gene","Sample","Counts")
    melt.df$LogCounts <- log(melt.df$Counts+1)
    melt.df$Condition <- rep(0,length.out=nrow(melt.df))
    for (c in conditions) {
        melt.df[melt.df$Sample %in% subset(metadata,subset=Condition==c)$Sample,]$Condition = c
    }
    vln_plot <- ggplot(melt.df,aes(x=Sample,y=LogCounts,fill=Condition))+geom_violin(trim=FALSE)+ylim(log_min,log_max)+ggtitle(curr.system)+ylab("Log-Counts")+xlab("Sample") 
    vln_plot = vln_plot + geom_boxplot(width=0.1,fill="white")
    if (is.null(legend)) {
        legend <- g_legend(vln_plot)
    }
    vln_plot = vln_plot + theme(axis.text.x = element_text(size=8,angle=45,vjust=0.5), legend.position = "none")
    plot_list[[i]] <- vln_plot
}
plot_list[[length(subsystems)+1]] <- legend
vln_png = paste(subsystem.level,"_Subsystem_Distribution_mqc.png",sep="")
png(vln_png,width=png_width,height=png_height)
do.call("grid.arrange",c(plot_list,ncol=num_columns))
dev.off()
