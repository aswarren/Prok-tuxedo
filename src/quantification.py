#!/usr/bin/env python

import os,sys,glob,math,shutil,subprocess

#Run the feature count program specified in json input, or run htseq-count by default
def run_featurecount(genome_list, condition_dict, parameters, output_dir, job_data, pipeline_log):
    #check parameters for stringtie option. Assuming htseq2
    program = job_data.get("feature_count","htseq")
    if program == "htseq":
        run_htseq_count(genome_list, condition_dict, parameters, job_data, output_dir, pipeline_log)
    elif program == "stringtie":
        run_stringtie(genome_list, condition_dict, parameters, job_data, output_dir, pipeline_log)
    else: #not a valid program
        sys.stderr.write("Invalid feature count program: htseq or stringtie only\n")
        os.exit(1)

def run_stringtie(genome_list, condition_dict, parameters, job_data, output_dir, pipeline_log):
    thread_count= parameters.get("stringtie",{}).get("-p",0)
    #defaults to not searching for novel features, adds -e flag
    find_novel_features = job_data.get("novel_features",False) 
    skip_merged_annotation = False 
    if thread_count == 0:
        thread_count=2 #multiprocessing.cpu_count()
    for genome in genome_list:
        genome_file=genome["genome"]
        genome_link = genome["genome_link"]
        #transcriptome assembly
        gtf_list = []
        for library in condition_dict:
            for r in condition_dict[library]["replicates"]:
                cur_dir=os.path.dirname(os.path.realpath(r[genome_file]["bam"]))
                os.chdir(cur_dir)
                r[genome_file]["dir"]=cur_dir
                cuff_gtf=os.path.join(cur_dir,"transcripts.gtf")
                stringtie_cmd = ["stringtie",r[genome_file]["bam"],"-p",str(thread_count),"-A","gene_abund.tab"]
                #check if novel features is turned off: add -e flag
                if not find_novel_features:
                    stringtie_cmd = stringtie_cmd + ["-e"]
                stringtie_cmd = stringtie_cmd + ["-G",genome["annotation"],"-o",cuff_gtf]
                r[genome_file]["gtf"] = cuff_gtf
                gtf_list.append(cuff_gtf)
                print (" ".join(stringtie_cmd))
                if not os.path.exists(cuff_gtf):
                    pipeline_log.append(" ".join(stringtie_cmd))
                    subprocess.check_call(stringtie_cmd)
                else:
                    sys.stderr.write(cuff_gtf+" stringtie file already exists. skipping\n")
        #True: Skip merged annotation pipeline
        #False: Run merged annotation pipeline
        if skip_merged_annotation:
            continue 
        #merge reconstructed transcriptomes
        ##stringtie --merge -G <reference annotation> -o <merged annotation> <gtf list>
        os.chdir(genome["output"])
        #os.mkdir("merged_annotation")
        merge_file = os.path.join(genome["output"],"merged_annotation","merged.gtf")
        if not os.path.exists(merge_file):
            merge_cmd = ["stringtie","--merge","-G",genome["annotation"],"-o",merge_file]+gtf_list
            print(" ".join(merge_cmd))
            pipeline_log.append(" ".join(merge_cmd))
            subprocess.check_call(merge_cmd)
        genome["merged_annotation"] = merge_file 
        #requantify transcriptome results with merged annotation
        for library in condition_dict:
            for r in condition_dict[library]["replicates"]:
                cur_dir=os.path.dirname(os.path.realpath(r[genome_file]["bam"]))
                os.chdir(cur_dir)
                merge_gtf = os.path.join(cur_dir,"transcripts_merged.gtf")
                stringtie_cmd = ["stringtie",r[genome_file]["bam"],"-p",str(thread_count),"-A","merged_abund.tab","-e","-G",genome["merged_annotation"],"-o",merge_gtf]
                r[genome_file]["merged_gtf"] = merge_gtf
                if not os.path.exists(merge_gtf):
                    print (" ".join(stringtie_cmd))
                    pipeline_log.append(" ".join(stringtie_cmd))
                    subprocess.check_call(stringtie_cmd)
                else:
                    sys.stderr.write(merge_gtf+" stringtie file already exists. skipping\n")
    create_stringtie_gene_transcript_dict(genome_list,condition_dict,False if skip_merged_annotation else True)

#Sometimes strintie puts in STRG genes 
#Put a dictionary for each replicate into the condition_dict
#Swap STRG genes with transcript names when appropriate
def create_stringtie_gene_transcript_dict(genome_list,condition_dict,merged_flag):
    gtf_key = "merged_gtf" if merged_flag else "gtf"
    for genome in genome_list:
        genome_file=genome["genome"]
        genome_link = genome["genome_link"]
        for library in condition_dict:
            for rep in condition_dict[library]["replicates"]:
                gtf_file = rep[genome_file][gtf_key]
                gtf_dict = {} 
                gtf_dict["gene_to_transcript"] = {} 
                gtf_dict["transcript_to_gene"] = {} 
                with open(gtf_file,"r") as gtf_handle:
                    for line in gtf_handle:
                        if line[0] == "#":
                            continue
                        line = line.strip().split("\t")
                        if line[2] == "transcript": #transcript line is the only one with TPM values, which is mainly what this dictionary is used for
                            features = line[-1].split(";")
                            gene_id = features[0].replace("gene_id ","").replace("transcript_id ","").replace("\"","").replace("gene-","")
                            transcript_id = features[1].replace("gene_id ","").replace("transcript_id ","").replace("\"","").replace(" ","").replace("rna-","")
                            if gene_id not in gtf_dict["gene_to_transcript"]:
                                gtf_dict["gene_to_transcript"][gene_id] = []
                            gtf_dict["gene_to_transcript"][gene_id].append(transcript_id)  
                            gtf_dict["transcript_to_gene"][transcript_id] = gene_id
                rep[genome_file]["gtf_dict"] = gtf_dict
                        

def split_bam_file(replicate_bam,threads,pipeline_log):
    print("Splitting %s into %s files for parallel htseq-count"%(replicate_bam,str(threads)))
    #count number of reads in bam file
    count_cmd_list = ["samtools","view","-c","--threads",str(threads),replicate_bam]
    print(" ".join(count_cmd_list))
    pipeline_log.append(" ".join(count_cmd_list))
    count_cmd = subprocess.Popen(count_cmd_list,stdout=subprocess.PIPE)
    num_lines,err = count_cmd.communicate()
    num_lines = float(num_lines.strip())
    #number of lines per file
    num_lines_per_file = int(math.floor(num_lines/float(threads)))
    #Get header from bam file: htseq won't work without the header
    header_cmd = ["samtools","view","-H",replicate_bam]
    print(" ".join(header_cmd))
    pipeline_log.append(" ".join(header_cmd))
    header_output = subprocess.Popen(header_cmd,stdout=subprocess.PIPE)
    header,err = header_output.communicate()
    #separate bam file into multiple bam files
    view_cmd = ["samtools","view","--threads",str(threads),replicate_bam]
    split_cmd = ["split","-l",str(num_lines_per_file)]
    if os.path.exists("Split_Bams"):
        #if directory wasn't deleted for some reason, remove before running
        shutil.rmtree("Split_Bams") 
    print(" ".join(view_cmd))
    pipeline_log.append(" ".join(view_cmd))
    bam_var = subprocess.Popen(view_cmd,stdout=subprocess.PIPE)
    os.mkdir("Split_Bams")
    os.chdir("Split_Bams")
    print(" ".join(split_cmd))
    pipeline_log.append(" ".join(split_cmd))
    subprocess.check_call(split_cmd,stdin=bam_var.stdout)
    sys.stdout.flush() #Get warnings if buffer isn't flushed for some reason, buffer overflow? seems unlikely
    for sam in glob.glob("*"): #iterate through all split files and format
        new_sam = sam+".sam"
        with open(new_sam,"w") as ns:
            ns.write(header)
            sam_file = open(sam,"r")
            sam_lines = sam_file.readlines()
            ns.write("".join(sam_lines))
            sam_file.close()
        os.remove(sam)
    os.chdir("..")

def run_htseq_parallel(genome_annotation,replicate_bam,counts_file,strand,feature,feature_type,threads, pipeline_log):
    #run htseq-count
    htseq_cmd = ["htseq-count","-t",feature_type,"-m","intersection-nonempty","--nonunique=all","-f","sam","-r","pos","-s",strand,"-i",feature,"-n",str(threads)]
    #htseq_cmd = ["htseq-count","-t",feature_type,"-f","sam","-r","pos","-s",strand,"-i",feature,"-n",str(threads)]
    for sam in glob.glob("Split_Bams/*"): #iterate through all split files and format
        htseq_cmd.append(sam)
    htseq_cmd.append(genome_annotation)
    print(" ".join(htseq_cmd))
    pipeline_log.append(" ".join(htseq_cmd))
    #htseq_stdout = subprocess.Popen(htseq_cmd,stdout=subprocess.PIPE,shell=True)
    htseq_stdout = subprocess.Popen(htseq_cmd,stdout=subprocess.PIPE)
    print("running htseq-parallel: communicate()")
    htseq_output,htseq_err = htseq_stdout.communicate()
    sys.stdout.flush()
    print("finished communicate(), writing to Output.txt")
    with open("Output.txt","w") as o:
        o.write(htseq_output)
    #    subprocess.check_call(htseq_cmd,stdout=o)
    #Combine output into counts file and delete Split_Bams directory
    #TODO: maybe rewrite using a bash command?
    with open("Output.txt","r") as split_counts:
        with open(counts_file,"w") as cf:
            for line in split_counts:
                line = line.strip().split()
                gene = line[0]
                counts = [int(x) for x in line[1:]]
                cf.write("%s\t%s\n"%(gene,sum(counts))) 
    os.remove("Output.txt")

# -s: (yes,no,reverse) 
# -i: feature to look for in annotation file (final column)
# -t: feature type to be used (3rd column), all others ignored. default = gene
def run_htseq_count(genome_list, condition_dict, parameters, job_data, output_dir, pipeline_log):
    strand = job_data.get("htseq",{}).get("-s","no")
    feature = job_data.get("htseq",{}).get("-i","ID")
    feature_type = job_data.get("htseq",{}).get("-t","gene")
    #recipe to check for host and change parameters
    recipe = job_data.get("recipe","RNA-Rocket")
    threads = parameters.get("htseq",{}).get("-p","1")
    for genome in genome_list:
        genome_file = genome["genome"]
        genome_annotation = genome["annotation"]
        for condition in condition_dict:
            for replicate in condition_dict[condition]["replicates"]:
                cur_dir = os.path.dirname(os.path.realpath(replicate[genome_file]["bam"]))
                os.chdir(cur_dir)
                replicate[genome_file]["dir"] = cur_dir
                counts_file = os.path.basename(replicate[genome_file]["bam"]).replace(".bam",".counts")
                counts_file_path = os.path.join(cur_dir,counts_file) 
                #Add counts file to replicate entry
                if recipe == "Host":
                    replicate[genome_file]["gene_counts"] = counts_file_path.replace(".counts",".gene.counts") 
                    replicate[genome_file]["transcript_counts"] = counts_file_path.replace(".counts",".transcript.counts") 
                else:
                    replicate[genome_file]["counts"] = counts_file_path
                if int(threads) > 1:
                    if recipe == "Host":
                        #Split the bam file
                        if not os.path.exists(replicate[genome_file]["gene_counts"]) or not os.path.exists(replicate[genome_file]["transcript_counts"]):
                            split_bam_file(replicate[genome_file]["bam"],threads,pipeline_log)
                        if os.path.exists(replicate[genome_file]["gene_counts"]):
                            sys.stderr.write("%s exists for genome file %s: skipping htseq-count for genes\n"%(replicate[genome_file]["gene_counts"],genome_file))
                        else:
                            #Run gene_count matrix parameters
                            #Create two counts files and append them together
                            genes_file = os.path.basename(replicate[genome_file]["gene_counts"])
                            genes_file_1 = genes_file+".tmp" 
                            run_htseq_parallel(genome_annotation,replicate[genome_file]["bam"],genes_file,strand,"ID","gene",threads,pipeline_log)
                            run_htseq_parallel(genome_annotation,replicate[genome_file]["bam"],genes_file_1,strand,"ID","pseudogene",threads,pipeline_log)
                            cf_open = open(genes_file,"a")
                            cf1_open = open(genes_file_1,"r")
                            cf1_lines = cf1_open.readlines()
                            cf1_open.close()
                            cf_open.write("".join(cf1_lines))
                            cf_open.close()
                            os.remove(genes_file_1)
                        if os.path.exists(replicate[genome_file]["transcript_counts"]):
                            sys.stderr.write("%s exists for genome file %s: skipping htseq-count for transcripts\n"%(replicate[genome_file]["transcript_counts"],genome_file))
                        else:
                            #Run transcript_count matrix parameters
                            transcript_file = os.path.basename(replicate[genome_file]["transcript_counts"])
                            feature_list = ["mRNA","lnc_RNA","transcript","snRNA","V_gene_segment","snoRNA","enhancer","biological_region","primary_transcript","miRNA","C_gene_segment","rRNA","tRNA"]
                            for f in feature_list:
                                transcript_file_tmp = transcript_file+".tmp" 
                                run_htseq_parallel(genome_annotation,replicate[genome_file]["bam"],transcript_file_tmp,strand,"ID",f,threads,pipeline_log)
                                tf_open = open(transcript_file,"a")
                                tf1_open = open(transcript_file_tmp,"r") 
                                tf1_lines = tf1_open.readlines()
                                tf1_open.close()
                                tf_open.write("".join(tf1_lines))
                                tf_open.close()
                                os.remove(transcript_file_tmp)
                    else: #bacterial pipeline
                        if os.path.exists(counts_file_path):
                            sys.stderr.write("%s exists for genome file %s: skipping htseq-count\n"%(counts_file,genome_file))
                        else:
                            split_bam_file(replicate[genome_file]["bam"],threads, pipeline_log)
                            run_htseq_parallel(genome_annotation,replicate[genome_file]["bam"],counts_file,strand,feature,feature_type,threads,pipeline_log)
                    if os.path.exists("Split_Bams"):
                        shutil.rmtree("Split_Bams")
                else:
                    if recipe == "Host":
                        print("Host htseq pipeline not enabled for 1 thread, will take too long")
                        sys.exit(1)
                    print("running htseq-count and writing to %s"%counts_file)
                    htseq_cmd = ["htseq-count","-t",feature_type,"-f","bam","-r","pos","-s",strand,"-i",feature,replicate[genome_file]["bam"],genome_annotation]
                    print(" ".join(htseq_cmd))
                    #prints to stdout, so redirect output to file
                    with open(counts_file,"w") as cf:
                        subprocess.check_call(htseq_cmd,stdout=cf)


