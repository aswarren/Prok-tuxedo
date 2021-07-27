#!/usr/bin/env python3

import sys,os,subprocess

#Writes introduction text to introduction.pipeline, which is read in by the module
#all variables are assumed to be passed in as strings, so cast to ints when necessary
def write_introduction_pipeline(recipe,feature_count,num_samples,num_conditions,num_comparisons):
    ###First part of the output section
    output_str = "The " + recipe + " recipe was executed with " + num_samples      
    output_str = output_str + " samples"
    ###If differential expression is not turned on, do not add condition and comparison content
    if int(num_comparisons) > 0:
        output_str = output_str + " across " + num_conditions + " conditions and " + num_comparisons + " differential expression comparisons.<br/>"
    else:
        output_str = output_str + ".<br/>"
    ###Paragraph section
    output_str = output_str + get_intro_paragraph() + "<br/>"
    ###QC Report section
    #output_str = output_str + "<b>Sample assessment portion</b><br/>"
    ###Tool order section 
    output_str = output_str + "Structure of the report:<br/>"
    output_str = output_str + "- <b>General Statistics</b>: A table providing summary information about the quantity and quality of alignments for each sample<br/>"
    output_str = output_str + "- <b>Introduction</b>: Provides an introduction to the service and summarizes the structure of the report<br/>"
    output_str = output_str + "- <b>Fastqc</b>: A quality control too for assessing high throught sequence data<br/>"  
    output_str = output_str + "- <b>RSeQC</b>: A package providing a number of modules used to evaluated high throughput RNA-seq data<br/>"
    if feature_count == "htseq":
        output_str = output_str + "- <b>HTSeq</b>: Analysis of high-throughput sequencing data for abundence estimation<br/>"  
    else: #samtools
        output_str = output_str + "- <b>Samtools</b>: A set of utilities used to interact with and post-process short sequence aligments in SAM, BAM, and CRAM formats<br/>"
        
    output_str = output_str + "- <b>Pathways</b>: Violin plots depicting the distributions for PATRIC Subsystems and KEGG Pathways<br/>"  
    if int(num_comparisons) > 0:
        output_str = output_str + "- <b>Differential Gene Expression</b>: Displays a heatmap of the top 50 differentially expressed genes among all comparisons and a set of volcano plots for each condition comparison<br/>"  
    output_str = output_str + "- <b>Samstat</b>: Displays statistics for sequence files from next generation sequencing projects<br/>"
    output_str = output_str + "- <b>References</b>: A list of references for the tools utilized in this RNASeq analysis protocol"
    ###write to output file
    with open("introduction.pipeline","w") as o:
        o.write(output_str)

def get_intro_paragraph():
    output = "The BVBRC transcriptomic service is a multi-step pipeline used to summarize RNA-Seq experimental results. Sample quality is assessed using FastQC, Samtools and Samstat. Condensed FastQC and Samstat reports can be found at the beginning and end of this report, respectively."
    output = output + " Alignment is performed using Bowtie2 (bacterial) and Hisat2 (host). Abundance estimation is performed using HTSeq, Stringtie, or Cufflinks, which depends on the invoked recipe. Pathway distribution visualizations, using PATRIC Subsystems and Kegg Pathways if present for the specified genome, have also been provided. If differential expression is enabled, then the abundance estimation ouput is processed through DESeq2 (HTSeq, Stringtie) or cuffdiff (Cufflinks). Volcano Plots and a heatmap isolating the top 50 differentially expressed genes are provided to visualize this information."
    return output
