#!/usr/bin/env python

import os,sys,subprocess

def create_counts_table_host(genome_list,condition_dict,job_data):
    #Remove the last 5 lines from htseq-count output
    omit_list = ["__no_feature","__ambiguous","__too_low_aQual","__not_aligned","__alignment_not_unique"]
    for genome in genome_list:
        genome_file = genome["genome"]
        genome_annotation = genome["annotation"]
        genome_dir = genome["output"]
        feature_count = job_data.get("feature_count","htseq")
        #change to genome directory
        os.chdir(genome_dir)
        gene_set = set()
        transcript_set = set()
        replicate_list = []
        gc_dict = {}
        tc_dict = {}
        for condition in condition_dict:
            for replicate in condition_dict[condition]["replicates"]:
                gc_file = replicate[genome["genome"]]["gene_counts"]
                tc_file = replicate[genome["genome"]]["transcript_counts"]
                replicate_id = os.path.basename(gc_file).replace(".gene.counts","")
                replicate_list.append(replicate_id)
                gc_dict[replicate_id] = {}
                tc_dict[replicate_id] = {}
                with open(gc_file,"r") as gf:
                    for line in gf:
                        feature,count = line.strip().split("\t")
                        if feature not in omit_list:
                            gene_set.add(feature)
                            gc_dict[replicate_id][feature] = count
                with open(tc_file,"r") as tf:
                    for line in tf:
                        feature,count = line.strip().split("\t")
                        if feature not in omit_list:
                            transcript_set.add(feature)
                            tc_dict[replicate_id][feature] = count
        #output counts table
        genome_id = os.path.basename(genome_dir)
        #Delimeter: , (csv files)
        delim = "\t"
        gene_counts_mtx = genome_id+"."+feature_count+".gene_counts"
        transcript_counts_mtx = genome_id+"."+feature_count+".transcript_counts"
        with open(gene_counts_mtx,"w") as gcm:
            #write headers
            gcm.write("Feature")
            for replicate_id in replicate_list:
                gcm.write("%s%s"%(delim,replicate_id))
            gcm.write("\n")
            #write gene info 
            for gene in gene_set:
                gcm.write(gene)
                for replicate_id in replicate_list:
                    if gene in gc_dict[replicate_id]:
                        gcm.write("%s%s"%(delim,gc_dict[replicate_id][gene]))
                    else:
                        gcm.write("%s0"%delim)
                gcm.write("\n")
        genome["gene_matrix"] = os.path.join(genome["output"],gene_counts_mtx)
        with open(transcript_counts_mtx,"w") as tcm:
            #write headers
            tcm.write("Feature")
            for replicate_id in replicate_list:
                tcm.write("%s%s"%(delim,replicate_id))
            tcm.write("\n")
            for transcript in transcript_set:
                tcm.write(transcript)
                for replicate_id in replicate_list:
                    if transcript in tc_dict[replicate_id]:
                        tcm.write("%s%s"%(delim,tc_dict[replicate_id][transcript]))
                    else:
                        tcm.write("%s0"%delim)
                tcm.write("\n")
        genome["transcript_matrix"] = os.path.join(genome["output"],transcript_counts_mtx)

#Merges the counts file generated for each replicate from htseq-count for each genome. Outputs file to genome directory
#Names file according to genome identifier
#Parameters:
# - genome: The current genome dictionary object from genome_list
# - condition_dict: complete condition dictionary object
def create_counts_table(genome_list,condition_dict,job_data):
    #Remove the last 5 lines from htseq-count output
    omit_list = ["__no_feature","__ambiguous","__too_low_aQual","__not_aligned","__alignment_not_unique"]
    for genome in genome_list:
        genome_file = genome["genome"]
        genome_annotation = genome["annotation"]
        genome_dir = genome["output"]
        feature_count = job_data.get("feature_count","htseq")
        #change to genome directory
        #change to genome directory
        os.chdir(genome_dir)
        feature_set = set()
        replicate_list = []
        counts_dict = {}
        #load counts file information into memory
        for condition in condition_dict: 
            for replicate in condition_dict[condition]["replicates"]:
                counts_file = replicate[genome["genome"]]["counts"]
                replicate_id = os.path.basename(counts_file).replace(".counts","")
                replicate_list.append(replicate_id)
                counts_dict[replicate_id] = {}
                with open(counts_file,"r") as cf:
                    for line in cf:
                        feature,count = line.strip().split("\t")
                        #skip the last five lines in counts file
                        if feature not in omit_list:
                            feature_set.add(feature)
                            counts_dict[replicate_id][feature] = count
        #output counts table
        genome_id = os.path.basename(genome_dir)
        #Delimeter: , (csv files)
        delim = "\t"
        genome_counts_mtx = genome_id+"."+feature_count+".gene_counts"
        with open(genome_counts_mtx,"w") as gcm:
            #write headers
            gcm.write("Feature")
            for replicate_id in replicate_list:
                gcm.write("%s%s"%(delim,replicate_id))
            gcm.write("\n")
            #write feature info 
            for feature in feature_set:
                gcm.write(feature)
                for replicate_id in replicate_list:
                    if feature in counts_dict[replicate_id]:
                        gcm.write("%s%s"%(delim,counts_dict[replicate_id][feature]))
                    else:
                        gcm.write("%s0"%delim)
                gcm.write("\n")
        genome["gene_matrix"] = os.path.join(genome["output"],genome_counts_mtx)

#Put a metadata file in each genome directory
#Subsetting the data on current conditions can be done in R
def create_metadata_file(genome_list,condition_dict,output_dir):
    #From genome list, get any genome key to get a replicate identifier  
    genome_key = genome_list[0]["genome"]
    #stores conditions as keys and replicate list as value
    info_dict = {}
    for condition in condition_dict:
        info_dict[condition] = [] 
        for replicate in condition_dict[condition]["replicates"]:
            replicate_id = os.path.basename(replicate[genome_key]["bam"]).replace(".bam","")
            info_dict[condition].append(replicate_id)
    #write metadata file: column 1 == replicate ids and column 2 == condition
    #Print file to top level output directory, only 1 metadat file is needed
    #Get the top level output directory
    metadata_file = os.path.join(output_dir,"Metadata.txt")
    with open(metadata_file,"w") as mf:
        mf.write("Sample\tCondition\n")
        for condition in info_dict:
            for replicate in info_dict[condition]:
                mf.write("%s\t%s\n"%(replicate,condition))
    #Add metadata file to each genome in genome_list
    for genome in genome_list:
        genome["deseq_metadata"] = metadata_file

#Writes the input file for prepDE.py, a prerequisite for running the Stringtie -> DESEq2 pipeline
#Contents: <Sample> <Path/To/GTF>
def write_gtf_list(genome_list,condition_dict):
     for genome in genome_list:
        genome_file = genome["genome"]
        genome_annotation = genome["annotation"]
        genome_dir = genome["output"]
        #change to genome directory
        os.chdir(genome_dir)
        rep_gtf_list = []
        for condition in condition_dict:
            for replicate in condition_dict[condition]["replicates"]:
                replicate_id = os.path.basename(replicate[genome_file]["bam"]).replace(".bam","")
                #TODO: set as based on a run condition
                #Set to false if the stringtie --merge command was run
                skip_merged_annotation = True
                if skip_merged_annotation:
                    rep_gtf_file = replicate[genome_file]["gtf"]
                else:
                    rep_gtf_file = replicate[genome_file]["merged_gtf"]
                rep_gtf_list.append((replicate_id,rep_gtf_file))
        gtf_path_filename = genome_dir+"_GTF_Sample_Paths.txt"                 
        with open(gtf_path_filename,"w") as gpf:
            for rep_gtf in rep_gtf_list:
                gpf.write("%s\n"%("\t".join(rep_gtf)))
        genome["prepDE_input"] = gtf_path_filename

#Calls the prepDE.py script that transforms stringtie output into a format usable by DESeq2
def prep_stringtie_diffexp(genome_list,condition_dict,host_bool,pipeline_log):
    for genome in genome_list:
        avg_length = str(average_read_length_total(condition_dict,genome))
        genome_file = genome["genome"]
        genome_dir = genome["output"]
        genome_id = os.path.basename(genome_dir)
        os.chdir(genome_dir)
        genome_counts_mtx = genome_id+".stringtie.gene_counts" 
        genome["gene_matrix"] = os.path.join(genome["output"],genome_counts_mtx)
        prep_cmd = ["prepDE.py","-i",genome["prepDE_input"],"-l",avg_length,"-g",genome_counts_mtx]
        if host_bool:
            transcript_counts_mtx = genome_id+".stringtie.transcript_counts"
            genome["transcript_matrix"] = os.path.join(genome["output"],transcript_counts_mtx)
            prep_cmd+=["-t",transcript_counts_mtx]
        if not os.path.exists(genome_counts_mtx):
            print(" ".join(prep_cmd))
            pipeline_log.append(" ".join(prep_cmd))
            subprocess.check_call(prep_cmd)

def average_read_length_total(condition_dict,genome):
    num_replicates = 0
    total_length = 0
    for condition in condition_dict:
        for r in condition_dict[condition]["replicates"]:
            num_replicates = num_replicates + 1
            total_length = total_length + int(r[genome["genome"]]["avg_read_length"])
    avg_length = int(total_length/num_replicates)
    return avg_length


