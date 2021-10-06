#!/opt/patric-common/runtime/bin/python3

import os,sys,glob
import pandas,json

#NOTE: unit tests are not supporting the cufflinks pipeline at this time 
#NOTE: does not check for the existence of files that should not exist

#flag either of the unit tests for RNASeq
#basic: checks if the subdirectory structure in the output is correct: should it be before or after cleanup files???
#other: given a file with a list of N genes, checks the top TPM genes and top differentially expressed genes 
def run_unit_test(test_param,genome_list,condition_dict,output_dir,contrast_list,job_data):
    if test_param == 'basic':
        result = run_basic_test(genome_list,condition_dict,output_dir,contrast_list,job_data)
    else:
        #result = run_top_genes_test(test_param,genome_list,condition_dict,output_dir,contrast_list,job_data) TOP_N_GENES Test: replace with gene set test
        result = run_gene_set_test(test_param,genome_list,condition_dict,output_dir,contrast_list,job_data)
    if result == 0:
        sys.stderr.write("%s unit test passed with no errors\n"% "BASIC" if test_param == "basic" else "GENE_SET")
    else:
        sys.stderr.write("%s unit test finished with an error\n"%"BASIC" if test_param == "basic" else "GENE_SET")
    return result

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
        if job_data.get("feature_count","htseq") == "htseq": 
            gene_counts = os.path.join(genome["output"],".".join([genome_id,quant_method,"gene_counts","tsv"]))
        else:
            gene_counts = os.path.join(genome["output"],".".join([genome_id,quant_method,"gene_counts","csv"]))
        tpm_counts = os.path.join(genome["output"],".".join([genome_id,quant_method,"tpms","tsv"]))
        if job_data.get("recipe","RNA-Rocket") == "Host":
            if job_data.get("feature_count","htseq") == "htseq":
                transcript_counts = os.path.join(genome["output"],".".join([genome_id,quant_method,"transcript_counts","tsv"]))
            else:
                transcript_counts = os.path.join(genome["output"],".".join([genome_id,quant_method,"transcript_counts","csv"]))
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
        check_reps_flag = check_replicates(condition_dict,job_data,contrast_list)
        if len(contrast_list) > 0 and not job_data.get("novel_features",False):
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
                #check for bam file if bam file was not given as input
                #bam file would have been generated in pipeline, so if there was not read file then assume bam file 
                #TODO: currently there is  not strandedness determination for single-end libraries 
                #if "read" in rep or "read1" in rep:(if single-end strandedness is implemented)
                if "read1" in rep:
                    #strandedness file 
                    if "read2" in rep:
                        rep_strand_file = os.path.join(rep["target_dir"],rep_name+".infer")
                        expected_files.append(rep_strand_file)
                    rep_bam_file = os.path.join(rep["target_dir"],rep["name"]+".bam")
                    rep_bam_idx = os.path.join(rep["target_dir"],rep["name"]+".bam.bai")
                    expected_files.append(rep_bam_file)
                    expected_files.append(rep_bam_idx)
                if job_data.get("feature_count","htseq") == "htseq": 
                    if job_data.get("recipe","RNA-Rocket") == "Host":
                        quant_output_file = os.path.join(rep["target_dir"],rep["name"]+".gene.counts")
                        transcript_quant_file = os.path.join(rep["target_dir"],rep["name"]+".transcript.counts")
                        expected_files.append(transcript_quant_file)
                    else:
                        quant_output_file = os.path.join(rep["target_dir"],rep["name"]+".counts")
                else:
                    quant_output_file = os.path.join(rep["target_dir"],"transcripts.gtf")
                rep_samtools_stats = os.path.join(rep["target_dir"],rep["name"]+".samtools_stats")
                expected_files.append(quant_output_file)
                expected_files.append(rep_samtools_stats)
    return expected_files

def run_gene_set_test(test_json,genome_list,condition_dict,output_dir,contrast_list,job_data):
    #with open(json_file,"r") as tj_handle:
    #    test_json = json.loads(tj_handle.read()) 
    gene_list = test_json.get("genes_list",[])
    if len(gene_list) == 0:
        sys.stderr.write("No genes in test genes list, check json file\n")
        sys.exit(-1)
    #rpoB_Exposed_vs_rpoB_Unexposed.htseq.Genes.deseq2
    test_contrasts = test_json.get("contrasts")
    test_conditions = test_json.get("experimental_conditions")
    test_conditions = [x.replace(" ","_") for x in test_conditions]
    contrast_key = ".".join([test_conditions[test_contrasts[0]-1]+"_vs_"+test_conditions[test_contrasts[1]-1],job_data.get("feature_count","htseq"),"Genes.deseq2"])
    for genome in genome_list:
        os.chdir(genome["output"])
        gmx_file = genome["gmx"] 
        gmx_df = pandas.read_csv(gmx_file,sep="\t",header=0)
        top_hundred_genes = gmx_df.sort_values(contrast_key,ascending=False)["Gene_ID"].iloc[0:100].tolist() 
        #if 90% of genes in gene_list are in top_hundred_genes, test passes. Otherwise, it fails
        intersect_length = len(list(set(gene_list)&set(top_hundred_genes)))
        intersect_prop = float(intersect_length)/float(len(gene_list))*100
        print("proportion of shared differentially expressed genes = %s"%str(intersect_prop))
        if intersect_prop >= 90:
            sys.stderr.write("Differential Expression Unit Test: PASSED at {}%\n".format(str(intersect_prop)))
            return 0
        else:
            sys.stderr.write("Differential Expression Unit Test: FAILED at {}%\n".format(str(intersect_prop)))
            return -2

def run_top_genes_test(test_file,genome_list,condition_dict,output_dir,contrast_list,job_data):
    genes_list = []
    quant_method = "htseq" if job_data.get("feature_count","htseq") == "htseq" else "stringtie" 
    with open(test_file,"r") as tf:
        for line in tf:
            genes_list.append(line.strip())
    for genome in genome_list:
        genome_id = job_data.get("reference_genome_id",None)
        tpm_counts = os.path.join(genome["output"],".".join([genome_id,quant_method,"tpms","tsv"]))
        tpm_df = pandas.read_csv(tpm_counts,sep="\t",header=0)
        tpm_headers = list(tpm_df.columns.values) 
        tpm_headers = tpm_headers[1:len(tpm_headers)] #get rid of gene names column for looping over samples
        #Loop through the tpm_headers (samples)
        #get the number of genes at the top of the TPM list that are in the gene list
        #store that statistic and report for each sample
        tpm_stat_dict = {}
        for h in tpm_headers:
            curr_df = tpm_df.sort_values(h,ascending=False)
            top_n_genes = curr_df['Gene'].iloc[0:len(genes_list)].tolist() 
            tpm_stat_dict[h] = len(list(set(top_n_genes)&set(genes_list)))
        tpm_stats_out = [] 
        tpm_stats_out.append("Intersection of "+str(len(genes_list))+" Top Sample Genes with Unit Testing Genes:")
        for h in tpm_headers:
            tpm_stats_out.append(h+"\t"+str(tpm_stat_dict[h])) 
        print("\n".join(tpm_stats_out))
        #check differential expression results
        #see if genes in the list are up/down regulated in the top N genes
        if "gmx" not in genome:
            return #differential expression results do not exist, unless there was an error making the gmx file
        gmx_file = genome["gmx"] 
        gmx_df = pandas.read_csv(gmx_file,sep="\t",header=0)
        gmx_headers = list(gmx_df.columns.values)
        gmx_headers = gmx_headers[1:len(gmx_headers)] 
        #Loop through the gmx_headers (comparisons)
        #get the number of geens at the top of diff_exp list that are in the gene list
        #store that statistic and report for each comparison, upregulated and downregulated
        gmx_stat_dict = {}
        for h in gmx_headers:
            up_df = gmx_df.sort_values(h,ascending=False)
            down_df = gmx_df.sort_values(h,ascending=True)
            up_n_genes = up_df['Gene_ID'].iloc[0:len(genes_list)].tolist()
            down_n_genes = down_df['Gene_ID'].iloc[0:len(genes_list)].tolist()
            gmx_stat_dict[h] = {}
            gmx_stat_dict[h]["up"] = len(list(set(up_n_genes)&set(genes_list)))
            gmx_stat_dict[h]["down"] = len(list(set(down_n_genes)&set(genes_list)))
        gmx_down_stats = []
        gmx_down_stats.append("Intersection of "+str(len(genes_list))+" Top Downregulated Genes with Unit Testing Genes:")
        gmx_up_stats = []
        gmx_up_stats.append("Intersection of "+str(len(genes_list))+" Top Upregulated Genes with Unit Testing Genes:")
        for h in gmx_headers:
            gmx_down_stats.append(h+"\t"+str(gmx_stat_dict[h]["down"]))
            gmx_up_stats.append(h+"\t"+str(gmx_stat_dict[h]["up"]))
        print("\n".join(gmx_up_stats))
        print("\n".join(gmx_down_stats))

def check_replicates(condition_dict,job_data,contrasts):
    exp_dict = {}
    for condition in condition_dict:
        exp_dict[condition] = 0
        for rep in condition_dict[condition]["replicates"]:
            exp_dict[condition] = exp_dict[condition] + 1
    exp_conds = job_data.get("experimental_conditions",[])
    if len(exp_conds) == 0:
        return False
    for pair in contrasts:
        c1 = pair[0]
        c2 = pair[1]
        if exp_dict[c1] < 2 or exp_dict[c2] < 2:
            return False
    return True
