#!/opt/patric-common/runtime/bin/python3

import os,sys,glob

#NOTE: unit tests are not supporting the cufflinks pipeline at this time 
#NOTE: does not check for the existence of files that should not exist

#flag either of the unit tests for RNASeq
#basic: checks if the subdirectory structure in the output is correct: should it be before or after cleanup files???
#other: given a file with a list of N genes, checks the top TPM genes and top differentially expressed genes 
def run_unit_test(test_param,genome_list,condition_dict,output_dir,contrast_list,job_data):
    if test_param == 'basic':
        result = run_basic_test(genome_list,condition_dict,output_dir,contrast_list,job_data)
    else:
        result = run_top_genes_test(test_file,genome_list,condition_dict,output_dir,contrast_list,job_data)
    if result == 0:
        sys.stderr.write("%s unit test passed with no errors\n"% "BASIC" if test_param == "basic" else "TOP_N_GENES")
    else:
        sys.stderr.write("%s unit test finished with an error\n"%"BASIC" if test_param == "basic" else "TOP_N_GENES")

#Go through the actual output directory and confirm everything in the expected files list is found 
#If everything is found then the test passes, otherwise it fails and outputs the names of files not found
def run_basic_test(genome_list,condition_dict,output_dir,contrast_list,job_data):
    expected_file_list = get_expected_files(genome_list,condition_dict,output_dir,contrast_list,job_data)
    for item in expected_file_list[:]: #iterate over a copy of the list
        if os.path.exists(item):
            expected_file_list.remove(item) 
    if len(expected_file_list) == 0:
        return 0
    else:
        sys.stderr.write("BASIC test failed: The following files were not found: \n")
        sys.stderr.write("%s\n"%"\n".join(expected_file_list))
    return 1

#Go through and create the expected filenames for outputs
def get_expected_files(genome_list,condition_dict,output_dir,contrast_list,job_data):
    expected_files = []
    quant_method = "htseq" if job_data.get("feature_count","htseq") == "htseq" else "stringtie" 
    genome_id = job_data.get("reference_genome_id",None)
    if not genome_id:
        return 1 #shouldn't happen at all but put here just in case
    #files at the top output_dir level
    pipeline_file = os.path.join(output_dir,"Pipeline.txt") 
    expected_files.append(pipeline_file)
    #files at the genome level
    for genome in genome_list:
        #abundance estimages
        gene_counts = os.path.join(genome["output"],".".join([genome_id,quant_method,"gene_counts","tsv"]))
        tpm_counts = os.path.join(genome["output"],".".join([genome_id,quant_method,"tpms","tsv"]))
        if job_data.get("recipe","RNA-Rocket") == "Host":
            transcript_counts = os.path.join(genome["output"],".".join([genome_id,quant_method,"transcript_counts","tsv"]))
            expected_files.append(transcript_counts)
        expected_files.append(gene_counts)
        expected_files.append(tpm_counts)
        #strandedness output file
        strand_file = os.path.join(genome["output"],"library_geometry.txt")
        expected_files.append(strand_file)
        #multiqc report
        mqc_report = os.path.join(genome["output"],genome_id + "_report.html")
        expected_files.append(mqc_report)
        #if differential expression, check for contrast files at genome folder level
        if len(contrast_list) > 0:
            #gene matrix output
            gmx_file = os.path.join(genome["output"],"gene_exp.gmx")
            expected_files.append(gmx_file)
            #DESeq2 output
            for c in contrast_list:  
                contrast_file = os.path.join(genome["output"],".".join([c[0]+"_vs_"+c[1],quant_method,"Genes.deseq2.tsv"]))
                expected_files.append(contrast_file)
        #files within a condition/replicate 
        for library in condition_dict:
            for rep in condition_dict[library]["replicates"]:
                rep_name = rep["name"]
                #strandedness file 
                rep_strand_file = os.path.join(rep["target_dir"],rep_name+".infer")
                expected_files.append(rep_strand_file)
                #check for bam file if bam file was not given as input
                #bam file would have been generated in pipeline, so if there was not read file then assume bam file 
                if "read" in rep or "read1" in rep:
                    rep_bam_file = os.path.join(rep["target_dir"],rep["name"]+".bam")
                    rep_bam_idx = os.path.join(rep["target_dir"],rep["name"]+".bam.bai")
                    expected_files.append(rep_bam_file)
                    expected_files.append(rep_bam_idx)
                if job_data.get("feature_count","htseq") == "htseq": 
                    quant_output_file = os.path.join(rep["target_dir"],rep["name"]+".counts")
                else:
                    quant_output_file = os.path.join(rep["target_dir"],rep["name"]+".gtf")
                rep_samtools_stats = os.path.join(rep["target_dir"],rep["name"]+".samtools_stats")
                expected_files.append(quant_output_file)
                expected_files.append(rep_samtools_stats)
    return expected_files

def run_top_genes_test(test_file,genome_list,condition_dict,output_dir,contrast_list,job_data):
    genes_list = []
    with open(test_file,"r") as tf:
        for line in tf:
            genes_list.append(tf.strip())
    
