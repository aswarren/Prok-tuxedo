#!/usr/bin/env python

import sys,os,subprocess

#Writes introduction text to introduction.pipeline, which is read in by the module
#all variables are assumed to be passed in as strings, so cast to ints when necessary
def write_introduction_pipeline(recipe,num_samples,num_conditions,num_comparisons):
    ###First part of the output section
    output_str = "The " + recipe + " recipe was executed with " + num_samples      
    output_str = output_str + " samples"
    ###If differential expression is not turned on, do not add condition and comparison content
    if int(num_conditions) > 0:
        output_str = output_str + " across " + num_conditions + " conditions and " + num_comparisons + " differential expression comparisons.<br/>"
    else:
        output_str = output_str + ".<br/>"
    ###QC Report section
    output_str = output_str + "<b>Sample assessment portion</b><br/>"
    ###Tool order section 
    output_str = output_str + "Report Contents Order:<br/>"
    output_str = output_str + "- Fastqc<br/>"  
    output_str = output_str + "- HTSeq<br/>"  
    output_str = output_str + "- Subsystems<br/>"  
    output_str = output_str + "- Differential Gene Expression"  
    ###write to output file
    with open("introduction.pipeline","w") as o:
        o.write(output_str)
