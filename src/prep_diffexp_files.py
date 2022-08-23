#!/usr/bin/env python3

import os,sys,glob,subprocess,concurrent.futures,shutil

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
        gc_dict = {} #gene count dictionary
        tc_dict = {} #transcript count dictionary
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
        genome_id = genome["genome_id"]
        #Delimeter: , (csv files)
        delim = "\t"
        gene_counts_mtx = genome_id+"."+feature_count+".gene_counts.tsv"
        transcript_counts_mtx = genome_id+"."+feature_count+".transcript_counts.tsv"
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
        genome_id = genome["genome_id"]
        #Delimeter: , (csv files)
        delim = "\t"
        genome_counts_mtx = genome_id+"."+feature_count+".gene_counts.tsv"
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
                skip_merged_annotation = False 
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
        genome["prepDE_list"] = rep_gtf_list

#Calls the prepDE.py script that transforms stringtie output into a format usable by DESeq2
def prep_stringtie_diffexp(genome_list,condition_dict,host_bool,pipeline_log):
    for genome in genome_list:
        avg_length = str(average_read_length_total(condition_dict,genome))
        genome_file = genome["genome"]
        genome_dir = genome["output"]
        genome_id = genome["genome_id"]
        os.chdir(genome_dir)
        genome_counts_mtx = genome_id+".stringtie.gene_counts.csv" 
        genome["gene_matrix"] = os.path.join(genome["output"],genome_counts_mtx)
        rep_gtf_list = genome["prepDE_list"]
        for rep, genome_gtf_tmp in rep_gtf_list:
            genome_gtf_tmp2 = genome_gtf_tmp+".tmp"
            with open(genome_gtf_tmp,"r") as tmp_g, open(genome_gtf_tmp2,'w') as tmp_g2:
                for line in tmp_g:
                    if line[0] == "#":
                        tmp_g2.write(line)
                    else:
                        parts = line.strip().split("\t")
                        if parts[2] in ["exon", "transcript"]:
                            parsed_line = [tuple(i.strip().split(" ")[0:2]) for i in parts[-1].replace('"',"").split(";") if i.strip()]
                            ordered_keys = [i[0] for i in parsed_line]
                            features = dict(parsed_line)
                            gene_id = features.get("gene_id","").replace("\"","").replace("gene-","")
                            transcript_id = features.get("transcript_id","").replace("\"","").replace("rna-","")
                            ref_gene_name = features.get("ref_gene_name",features.get("gene_name",""))
                            features["transcript_id"] = transcript_id
                            if (gene_id.startswith("MSTRG.") or gene_id.startswith("STRG.")) and ref_gene_name:
                                ordered_keys.append("stie_id")
                                features["stie_id"] = gene_id
                                features["gene_id"] = ref_gene_name
                            else:
                                features["gene_id"] = gene_id
                            parts[-1] = "; ".join([i+' "'+features[i]+'"' for i in ordered_keys])
                        tmp_g2.write("\t".join(parts)+"\n")
            shutil.move(genome_gtf_tmp2, genome_gtf_tmp) 
                                
        prep_cmd = ["prepDE.py","-i",genome["prepDE_input"],"-l",avg_length,"-g",genome_counts_mtx]
        '''
        if host_bool:
            transcript_counts_mtx = genome_id+".stringtie.transcript_counts.csv"
            genome["transcript_matrix"] = os.path.join(genome["output"],transcript_counts_mtx)
            prep_cmd+=["-t",transcript_counts_mtx]
        '''
        transcript_counts_mtx = genome_id+".stringtie.transcript_counts.csv"
        genome["transcript_matrix"] = os.path.join(genome["output"],transcript_counts_mtx)
        prep_cmd+=["-t",transcript_counts_mtx]
        if not os.path.exists(genome_counts_mtx) or not os.path.exists(transcript_counts_mtx):
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

def create_tpm_matrix_htseq(genome_list,condition_dict,host_flag,num_threads,pipeline_log):
    #TPMCalculator -g ~/Genomes/9606/GCF_000001405.39_GRCh38.p13_genomic.mod.gtf -b Test_Htseq/9606/avirulent/replicate1/SRR10307420.bam -e
    #might need -a flag to output genes/transcripts with 0 counts
    #tpm_calc_cmd = ["TPMCalculator","-g",<genome_gtf>,"-b",<rep_bam>,<transcript_flag>]
    for genome in genome_list:
        os.chdir(genome["output"])
        genome_file=genome["genome"]
        genome_gff = genome["annotation_link"]
        genome_gtf_tmp = genome_gff.replace("gff","temp.gtf")
        genome_gtf = genome_gff.replace("gff","gtf")
        genome["gtf"] = genome_gtf
        gff_to_gtf_cmd = ["gffread",genome_gff,"-T","-o",genome_gtf_tmp]
        print(" ".join(gff_to_gtf_cmd))
        pipeline_log.append(" ".join(gff_to_gtf_cmd))
        if not os.path.exists(genome_gtf_tmp):
            subprocess.check_call(gff_to_gtf_cmd)
        if not os.path.exists(genome_gtf):
            gtf_lines = []
            with open(genome_gtf_tmp,"r") as tmp_g:
                for line in tmp_g:
                    if line[0] == "#":
                        gtf_lines.append(line.strip())
                    else:
                        if "gene_id" not in line and "gene_name" in line:
                            line = line.replace("gene_name","gene_id")
                        #TPMCalculator only works wth exon
                        #skip this if using bacterial annotations
                        if not host_flag:
                            line = line.replace("CDS","exon")
                        gtf_lines.append(line.strip())
            with open(genome_gtf,"w") as o:
                o.write("\n".join(gtf_lines))
        rm_tmp_gff = ["rm",genome_gtf_tmp]
        pipeline_log.append(" ".join(rm_tmp_gff))
        subprocess.check_call(rm_tmp_gff)
        ###Sometimes the gene_id field is not present: switch gene_name to gene_id if not present
        tpm_calc_list = []
        for condition in condition_dict: 
            for replicate in condition_dict[condition]["replicates"]:
                rep_bam = replicate[genome_file]["bam"]
                rep_id = os.path.basename(replicate[genome_file]["bam"]).replace(".bam","")
                rep_outdir = os.path.dirname(rep_bam)
                tpm_calc_cmd = ["TPMCalculator","-g",genome_gtf,"-b",rep_bam]
                if host_flag:
                    tpm_calc_cmd += ["-e"]
                tpm_stdout = os.path.join(os.path.dirname(replicate[genome_file]["bam"]),rep_id+"_tpm_calculator.out")
                tpm_stderr = os.path.join(os.path.dirname(replicate[genome_file]["bam"]),rep_id+"_tpm_calculator.err")
                #TPMCalculator takes too long and is single threaded, try to use the Pool api to run concurrently
                tpm_calc_tuple = (tpm_calc_cmd,tpm_stdout,tpm_stderr)
                if True:
                    tpm_calc_list.append(tpm_calc_tuple)
        for tpm_cmd in tpm_calc_list:
            pipeline_log.append(" ".join(tpm_cmd[0]))
        ###Run using a pool thread
        #calcs the run_tpm_calc function for executing the TPMCalculator program with a single thread
        tpm_dir = os.path.join(genome["output"],"TPMCalculator_Output")
        if not os.path.exists(tpm_dir):
            os.mkdir(tpm_dir)
        os.chdir(tpm_dir)
        with concurrent.futures.ThreadPoolExecutor(max_workers=num_threads) as pool:
            pool.map(run_tpm_calc,tpm_calc_list)
        ###Create the tpm matrix
        #assuming current directory is genome["output"]/TPMCalculator_Output
        tpm_dict = {}
        gene_set = set()
        for calc_file in glob.glob("*_genes.out"): 
            rep_id = calc_file.replace("_genes.out","")
            tpm_dict[rep_id] = {}
            with open(calc_file,"r") as cf:
                next(cf) 
                for line in cf:
                    line = line.strip().split()
                    gene = line[0]
                    gene = gene.replace("gene-","",).replace("rna-","")
                    if "#" in gene:
                        continue
                    gene_set.add(gene)
                    tpm_val = line[6]
                    tpm_dict[rep_id][gene] = tpm_val
        #write out tpm file and add to genome dictionary
        os.chdir(genome["output"])
        counts_matrix = genome["gene_matrix"]
        tpm_file = counts_matrix.replace("gene_counts","tpms")
        genome["genes_tpm_matrix"] = tpm_file
        with open(tpm_file,"w") as o:
            o.write("Gene")
            for rep_id in tpm_dict:
                o.write("\t%s"%rep_id)
            o.write("\n")
            for gene in gene_set:
                o.write("%s"%gene)
                for rep_id in tpm_dict:
                    tpm_val = tpm_dict[rep_id][gene] if gene in tpm_dict[rep_id] else "0"
                    o.write("\t%s"%tpm_val)
                o.write("\n")
            
def run_tpm_calc(job):
    print(" ".join(job[0]))
    with open(job[1],"w") as out, open(job[2],"w") as err:
        subprocess.check_call(job[0],stdout=out,stderr=err)

#either map them correctly or make the identifiers consisten (merge gtfs and rerun)
def create_tpm_matrix_stringtie(genome_list,condition_dict,host_flag):
    for genome in genome_list:
        tpm_dict = {}
        #counts_matrix = genome["gene_matrix"]
        counts_matrix = genome["transcript_matrix"]
        tpm_file = counts_matrix.replace("gene_counts","tpms").replace("transcript_counts","tpms").replace("csv","tsv")
        #feature_field = "gene_id"
        feature_field = "transcript_id"
        rep_order = []
        feature_set = set()
        for condition in condition_dict: 
            for rep in condition_dict[condition]["replicates"]:
                rep_id = os.path.basename(rep[genome["genome"]]["bam"]).replace(".bam","")
                rep_order.append(rep_id)
                tpm_dict[rep_id] = {}
                if "merged_gtf" in rep[genome["genome"]]:
                    rep_gtf = rep[genome["genome"]]["merged_gtf"]
                else:
                    rep_gtf = rep[genome["genome"]]["gtf"] 
                with open(rep_gtf,"r") as gtf_handle:
                    for line in gtf_handle:
                        if line[0] == "#":
                            continue
                        else:
                            parts = line.strip().split("\t")
                            if parts[2] in ["transcript"]:
                                try:
                                    parsed_line = [tuple(i.strip().split(" ")[0:2]) for i in parts[-1].replace('"',"").split(";") if i.strip()]
                                    features = dict(parsed_line)
                                except:
                                    sys.stderr.write("Problem parsing line for tpms file "+rep_gtf+"\n"+line)
                                    continue
                                ordered_keys = [i[0] for i in parsed_line]
                                gene_id = features.get("gene_id","").replace("\"","").replace("gene-","")
                                transcript_id = features.get("transcript_id","").replace("\"","").replace("rna-","")
                                ref_gene_name = features.get("ref_gene_name","")
                                tpm_val = features.get("TPM", 0)
                                if (gene_id.startswith("MSTRG.") or gene_id.startswith("STRG.")) and ref_gene_name:
                                    gene_id = ref_gene_name
                                if feature_field == "gene_id":
                                    feature_set.add(gene_id)
                                    tpm_dict[rep_id][gene_id] = tpm_val
                                else:
                                    feature_set.add(transcript_id)
                                    tpm_dict[rep_id][transcript_id] = tpm_val
        with open(tpm_file,"w") as o:
            o.write("Gene")
            for rep in tpm_dict:
                o.write("\t%s"%rep)
            o.write("\n")
            for gene in feature_set:
                o.write("%s"%gene)
                for rep_id in rep_order:
                    if gene not in tpm_dict[rep_id]:
                        o.write("\t0")
                    else:
                        o.write("\t%s"%(tpm_dict[rep_id][gene]))
                o.write("\n")
        genome["genes_tpm_matrix"] = tpm_file 
