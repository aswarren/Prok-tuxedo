#!/usr/bin/env python

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
    ###QC Report section
    #output_str = output_str + "<b>Sample assessment portion</b><br/>"
    ###Tool order section 
    output_str = output_str + "Report Contents Order:<br/>"
    output_str = output_str + "- <b>Fastqc</b>: A quality control too for assessing high throught sequence data.<br/>"  
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
