#!/usr/bin/env python

import os, sys
import argparse
import subprocess
import multiprocessing
import cuffdiff_to_genematrix
import glob
import tarfile, json
import shutil
import math

#take genome data structure and condition_dict and make directory names. processses condition to ensure no special characters, or whitespace
def make_directory_names(genome, condition_dict):
    if genome["dir"].endswith('/'):
        genome["dir"]=genome["dir"][:-1]
    genome["dir"]=os.path.abspath(genome["dir"])
    genome["output"]=os.path.join(os.path.abspath(output_dir),os.path.basename(genome["dir"]))
    for condition in condition_dict:
        rcount=0
        condition_index=condition_dict[condition]["condition_index"]
        for r in condition_dict[condition]["replicates"]:
            cur_cleanup=[]
            rcount+=1
            target_dir=os.path.join(genome["output"], str(condition_index),"replicate"+str(rcount))
            target_dir=os.path.abspath(target_dir)
            r["target_dir"]=target_dir

#hisat2 has problems with spaces in filenames
#prevent spaces in filenames. if one exists link the file to a no-space version.
def link_space(file_path):
    result=file_path
    name=os.path.splitext(os.path.basename(file_path))[0]
    if " " in name:
        clean_name=name.replace(" ","")
        result= file_path.replace(name,clean_name)
        if not os.path.exists(result):
            subprocess.check_call(["ln","-s",file_path,result])
    return result

#If the hisat_index is tared, untar it
#If the files exist and are untared, copy to correct location
#Return hisat index prefix that is passed into hisat command: -x <hisat_index_prefix>
def setup_hisat_indices(genome,output_dir,final_cleanup,job_data):
    if genome["hisat_index"].endswith("ht2.tar"):
        genome_link = genome["genome_link"]
        archive = tarfile.open(genome["hisat_index"])
        indices = [os.path.join(output_dir,os.path.basename(x)) for x in archive.getnames()]
        final_cleanup+=indices
        archive.close()
        subprocess.check_call(["tar","-xvf",genome["hisat_index"],"-C",output_dir])
        if job_data.get("host_ftp",None):
            index_prefix = os.path.join(output_dir,"/".join(genome["hisat_index"].split("/")[-2:]).replace("ht2.tar","")) 
        else:
            index_prefix = os.path.join(output_dir,os.path.basename(genome["hisat_index"]).replace(".ht2.tar",""))
    else:
        for index in glob.glob(os.path.join(genome["dir"],"*.ht2")): #for each hisat index file in the genome directory
            copy_index = os.path.join(output_dir,os.path.basename(index))
            if not os.path.exists(copy_index):#if it doesn't exist where the index files should be
                shutil.copy(index,copy_index) #copy it over
            final_cleanup+=[copy_index]
            index_prefix = os.path.join(output_dir,".".join(os.path.basename(genome["hisat_index"]).split(".")[:-2]))
    return index_prefix

#pretty simple: its for prokaryotes in that parameters will be attuned to give best performance and no tophat
#TODO: rewrite and separate bowtie2 and hisat2 runs into separate functions
#TODO: check for which feature count so samtools sort can be adjusted as needed (-r name for htseq)
def run_alignment(genome_list, condition_dict, parameters, output_dir, job_data): 
    #modifies condition_dict sub replicates to include 'bowtie' dict recording output files
    for genome in genome_list:
        genome_link = genome["genome_link"]
        final_cleanup=[]
        if "hisat_index" in genome and genome["hisat_index"]:
            index_prefix = setup_hisat_indices(genome,output_dir,final_cleanup,job_data) 
            cmd=["hisat2","--dta-cufflinks", "-x", index_prefix] 
            #cmd=["hisat2","--dta-cufflinks", "-x", index_prefix,"--no-spliced-alignment"] 
            thread_count= parameters.get("hisat2",{}).get("-p",0)
        else:
            try:
                subprocess.check_call(["bowtie2-build", genome_link, genome_link])
            except Exception as err:
                sys.stderr.write("bowtie build failed: %s %s\n" % (err, genome_link))
                os.exit(1)
            #cmd=["hisat2","--dta-cufflinks", "-x", genome_link, "--no-spliced-alignment"] 
            cmd=["bowtie2", "-x", genome_link]
            thread_count= parameters.get("bowtie2",{}).get("-p",0)
        if thread_count == 0:
            thread_count=2 #multiprocessing.cpu_count()
        cmd+=["-p",str(thread_count)]
        for condition in condition_dict:
            rcount=0
            for r in condition_dict[condition]["replicates"]:
                cur_cleanup=[]
                rcount+=1
                target_dir=r["target_dir"]
                fastqc_cmd=["fastqc","--outdir",target_dir]
                samstat_cmd=["samstat"]
                cur_cmd=list(cmd)
                if "read2" in r:
                    cur_cmd+=["-1",link_space(r["read1"])," -2",link_space(r["read2"])]
                    name1=os.path.splitext(os.path.basename(r["read1"]))[0].replace(" ","")
                    name2=os.path.splitext(os.path.basename(r["read2"]))[0].replace(" ","")
                    sam_file=os.path.join(target_dir,name1+"_"+name2+".sam")
                    fastqc_cmd+=[r["read1"],r["read2"]]
                else:
                    cur_cmd+=[" -U",link_space(r["read1"])]
                    name1=os.path.splitext(os.path.basename(r["read1"]))[0].replace(" ","")
                    sam_file=os.path.join(target_dir,name1+".sam")
                    fastqc_cmd+=[r["read1"]]
                cur_cleanup.append(sam_file)
                bam_file=sam_file[:-4]+".bam"
                samstat_cmd.append(bam_file)
                r[genome["genome"]]={}
                r[genome["genome"]]["bam"]=bam_file
                r[genome["genome"]]["fastqc"] = bam_file.replace(".bam","_fastqc.html") 
                cur_cmd+=["-S",sam_file]
                print " ".join(fastqc_cmd)
                #TODO: uncomment
                #subprocess.check_call(fastqc_cmd)
                if os.path.exists(bam_file):
                    sys.stderr.write(bam_file+" alignments file already exists. skipping\n")
                else:
                    print cur_cmd
                    #TODO: sent hisat/bowtie2 log output to file so multiqc can read it
                    alignment_log = bam_file.replace("bam","bowtie") 
                    with open(alignment_log,"w") as al:
                        subprocess.check_call(cur_cmd,stdout=al,stderr=al) #call bowtie2 or hisat2
                if not os.path.exists(bam_file):
                    sort_threads = parameters.get("samtools",{}).get("-@",1)
                    #TODO: testing -n option in samtools sort, which is passed into htseq
                    subprocess.check_call("samtools view -Su "+sam_file+" | samtools sort -o - - -@ "+str(sort_threads)+" > "+bam_file, shell=True)#convert to bam
                    subprocess.check_call("samtools index "+bam_file, shell=True)
                    #subprocess.check_call('samtools view -S -b %s > %s' % (sam_file, bam_file+".tmp"), shell=True)
                    #subprocess.check_call('samtools sort %s %s' % (bam_file+".tmp", bam_file), shell=True)
                print " ".join(samstat_cmd)
                #TODO: uncomment
                #subprocess.check_call(samstat_cmd)

                for garbage in cur_cleanup:
                    subprocess.call(["rm", garbage])
        for garbage in final_cleanup:
            subprocess.call(["rm", garbage])

#Run the feature count program specified in json input, or run htseq-count by default
def run_featurecount(genome_list, condition_dict, parameters, output_dir, job_data):
    #check parameters for stringtie option. Assuming htseq2
    program = job_data.get("feature_count","htseq")
    if program == "htseq":
        run_htseq_count(genome_list, condition_dict, parameters, job_data, output_dir)
    elif program == "stringtie":
        run_stringtie(genome_list, condition_dict, parameters, job_data, output_dir)
    else: #not a valid program
        sys.stderr.write("Invalid feature count program: htseq or stringtie only\n")
        os.exit(1)

#Always set default parameter behavior when updating functions
def run_stringtie(genome_list, condition_dict, parameters, job_data, output_dir):
    thread_count= parameters.get("stringtie",{}).get("-p",0)
    #defaults to not searching for novel features, adds -e flag
    find_novel_features = job_data.get("novel_features",False) 
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
                stringtie_cmd = ["stringtie",r[genome_file]["bam"],"-p",str(thread_count),"-A","gene_abund.tab","-M","1.0"]
                #check if novel features is turned off: add -e flag
                if not find_novel_features:
                    stringtie_cmd = stringtie_cmd + ["-e"]
                stringtie_cmd = stringtie_cmd + ["-G",genome["annotation"],"-o",cuff_gtf]
                r[genome_file]["gtf"] = cuff_gtf
                gtf_list.append(cuff_gtf)
                if not os.path.exists(cuff_gtf):
                    print " ".join(stringtie_cmd)
                    subprocess.check_call(stringtie_cmd)
                else:
                    print (" ".join(stringtie_cmd))
                    sys.stderr.write(cuff_gtf+" stringtie file already exists. skipping\n")
        #True: Skip merged annotation pipeline
        #False: Run merged annotation pipeline
        #TODO: set as based on a run condition
        skip_merged_annotation = True
        if skip_merged_annotation:
            return
        #merge reconstructed transcriptomes
        ##stringtie --merge -G <reference annotation> -o <merged annotation> <gtf list>
        os.chdir(genome["output"])
        #os.mkdir("merged_annotation")
        merge_file = os.path.join(genome["output"],"merged_annotation","merged.gtf")
        if not os.path.exists(merge_file):
            merge_cmd = ["stringtie","--merge","-G",genome["annotation"],"-o",merge_file]+gtf_list
            print(" ".join(merge_cmd))
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
                    subprocess.check_call(stringtie_cmd)
                else:
                    sys.stderr.write(merge_gtf+" stringtie file already exists. skipping\n")

#TODO: smallRNA estimation in a future release
def run_cufflinks(genome_list, condition_dict, parameters, output_dir):
    for genome in genome_list:
        if genome.get("host",False):
            continue
        genome_file=genome["genome"]
        genome_link = genome["genome_link"]
        annotation_only =  int(parameters.get("cufflinks",{}).get("annotation_only",0))
        if annotation_only:
            use_annotation="-G"
        else:
            use_annotation="-g"
        
        cmd=["cufflinks","-q","-G",genome["annotation"],"-b",genome_link,"-I","50"]
        thread_count= parameters.get("cufflinks",{}).get("-p",0)
        if thread_count == 0:
            thread_count=2 #multiprocessing.cpu_count()

        cmd+=["-p",str(thread_count)]
        for library in condition_dict:
            for r in condition_dict[library]["replicates"]:
                cur_dir=os.path.dirname(os.path.realpath(r[genome_file]["bam"]))
                os.chdir(cur_dir)
                cur_cmd=list(cmd)
                r[genome_file]["dir"]=cur_dir
                bam_file = r[genome_file]["bam"]
                # Review: Memory mapped location on each system
                # Attempt to copy to /dev/shm. cufflinks seeks a lot in the file.
                # If that fails, try tmp.
                #
                bam_to_use = None
                
                try:
                    bam_tmp = os.tempnam("/dev/shm", "CUFFL")
                    shutil.copy(bam_file, bam_tmp)
                    print "Copy succeeded to %s" % (bam_tmp)
                    bam_to_use = bam_tmp
                    
                except IOError as err:
                    os.unlink(bam_tmp)
                    bam_to_use = None
		    bam_tmp = None
                    sys.stderr.write("Can't copy %s to %s: %s\n" % (bam_file, bam_tmp, err))

                if bam_to_use == None:
                    try:
                        bam_tmp = os.tempnam(None, "CUFFL")
                        shutil.copy(bam_file, bam_tmp)
                        bam_to_use = bam_tmp

                    except IOError as err:
                        os.unlink(bam_tmp)
                        bam_to_use = None
			bam_tmp = None
                        sys.stderr.write("Can't copy %s to %s: %s\n" % (bam_file, bam_tmp, err))
                    
                if bam_to_use == None:
                    sys.stderr.write("Can't copy %s to tmp space\n" % (bam_file))
                    bam_to_use = bam_file
		    bam_tmp = None
                
                cur_cmd += [bam_to_use]
                cuff_gtf=os.path.join(cur_dir,"transcripts.gtf")
                if not os.path.exists(cuff_gtf):
                    print " ".join(cur_cmd)
		    try:
                        sys.stderr.write("Invoke cufflinks: %s\n" % (cur_cmd))
			subprocess.check_call(cur_cmd)
		    except Exception as e:
		    	if bam_tmp != None:
			    sys.stderr.write("remove temp %s in exception handler\n" % (bam_tmp))
			    os.unlink(bam_tmp)
                        sys.stderr.write("Cufflinks error: %s\n" % (e))
			raise
                else:
                    sys.stderr.write(cuff_gtf+" cufflinks file already exists. skipping\n")
		if bam_tmp != None:
		    sys.stderr.write("remove temp %s\n" % (bam_tmp))
		    os.unlink(bam_tmp)

#Differential expression pipeline using the cufflinks protocol 
def run_cuffdiff(genome_list, condition_dict, parameters, output_dir, gene_matrix, contrasts, job_data, map_args, diffexp_json):
    #run cuffquant on every replicate, cuffmerge on all resulting gtf, and cuffdiff on the results. all per genome.
    for genome in genome_list:
        genome_file=genome["genome"]
        genome_link = genome["genome_link"]
        is_host=genome.get("host",False)
        merge_manifest=os.path.join(genome["output"],"gtf_manifest.txt")
        merge_folder=os.path.join(genome["output"],"merged_annotation")
    	subprocess.call(["mkdir","-p",merge_folder])
        merge_file=os.path.join(merge_folder,"merged.gtf")
        if is_host:
            merge_cmd=["stringtie","--merge","-g","0","-G",genome["annotation"],"-o", merge_file]
            thread_count= parameters.get("stringtie",{}).get("-p",0)
        else:
            merge_cmd=["cuffmerge","-g",genome["annotation"],"-o", merge_folder]
            thread_count= parameters.get("cuffmerge",{}).get("-p",0)
        if thread_count == 0:
            thread_count=2 #multiprocessing.cpu_count()
        merge_cmd+=["-p",str(thread_count)]
        with open(merge_manifest, "w") as manifest: 
            for library in condition_dict:
                for r in condition_dict[library]["replicates"]:
                    manifest.write("\n"+os.path.join(r[genome_file]["dir"],"transcripts.gtf"))
        merge_cmd+=[merge_manifest]

        if not os.path.exists(merge_file):
            print " ".join(merge_cmd)
            subprocess.check_call(merge_cmd)
        else:
            sys.stderr.write(merge_file+" cuffmerge file already exists. skipping\n")

        #setup diff command
        cur_dir=genome["output"]
        thread_count= parameters.get("cuffdiff",{}).get("-p",0)
        if thread_count == 0:
            thread_count=2
        diff_cmd=["cuffdiff",merge_file,"-p",str(thread_count),"-b",genome_link,"-L",",".join(condition_dict.keys())]
        os.chdir(cur_dir)
        cds_tracking=os.path.join(cur_dir,"cds.fpkm_tracking")
        contrasts_file = os.path.join(cur_dir, "contrasts.txt")
        with open(contrasts_file,'w') as contrasts_handle:
            contrasts_handle.write("condition_A\tcondition_B\n")
            for c in contrasts:
                #contrasts_handle.write(str(c[0])+"\t"+str(c[1])+"\n")
                contrasts_handle.write(str(c[1])+"\t"+str(c[0])+"\n")
        diff_cmd+=["--contrast-file",contrasts_file,"-o",cur_dir]

        #create quant files and add to diff command
        for library in condition_dict:
            quant_list=[]
            for r in condition_dict[library]["replicates"]:
                #quant_cmd=["cuffquant",genome["annotation"],r[genome_file]["bam"]]
                quant_cmd=["cuffquant",merge_file,r[genome_file]["bam"]]
                cur_dir=r[genome_file]["dir"]#directory for this replicate/genome
                os.chdir(cur_dir)
                quant_file=os.path.join(cur_dir,"abundances.cxb")
                quant_list.append(quant_file)
                if not os.path.exists(quant_file):
                    subprocess.check_call(quant_cmd)
                else:
                    print " ".join(quant_cmd)
                    sys.stderr.write(quant_file+" cuffquant file already exists. skipping\n")
            diff_cmd.append(",".join(quant_list))

        cur_dir=genome["output"]
        os.chdir(cur_dir)
        cds_tracking=os.path.join(cur_dir,"cds.fpkm_tracking")
        if not os.path.exists(cds_tracking):
            print " ".join(diff_cmd)
            subprocess.check_call(diff_cmd)
        else:
            sys.stderr.write(cds_tracking+" cuffdiff file already exists. skipping\n")
        #write gene_gmx file now for cuffdiff pipeline
        de_file=os.path.join(cur_dir,"gene_exp.diff")
        gmx_file=os.path.join(cur_dir,"gene_exp.gmx")
        if os.path.exists(de_file) and not os.path.exists(gmx_file):
            cuffdiff_to_genematrix.main([de_file],gmx_file)

#Runs the differential expression import protocol 
def run_diff_exp_import(genome_list, condition_dict, parameters, output_dir, contrasts, job_data, map_args, diffexp_json):
    for genome in genome_list:
        cur_dir = genome["output"]
        gmx_file=os.path.join(cur_dir,"gene_exp.gmx")
        transform_script = "expression_transform.py"
        if os.path.exists(gmx_file):
            experiment_path=os.path.join(output_dir, map_args.d)
            subprocess.call(["mkdir","-p",experiment_path])
            transform_params = {"output_path":experiment_path, "xfile":gmx_file, "xformat":"tsv",\
                    "xsetup":"gene_matrix", "source_id_type":"patric_id",\
                    "data_type":"Transcriptomics", "experiment_title":"RNA-Seq", "experiment_description":"RNA-Seq",\
                    "organism":job_data.get("reference_genome_id")}
            diffexp_json["parameters"]=transform_params
            params_file=os.path.join(cur_dir, "diff_exp_params.json")
            with open(params_file, 'w') as params_handle:
                params_handle.write(json.dumps(transform_params))
            convert_cmd=[transform_script, "--ufile", params_file, "--sstring", map_args.sstring, "--output_path",experiment_path,"--xfile",gmx_file]
            print " ".join(convert_cmd)
            try:
               subprocess.check_call(convert_cmd)
            except(subprocess.CalledProcessError):
               sys.stderr.write("Running differential expression import failed.\n")
               #subprocess.call(["rm","-rf",experiment_path])
               return
            diffexp_obj_file=os.path.join(output_dir, os.path.basename(map_args.d.lstrip(".")))
            with open(diffexp_obj_file, 'w') as diffexp_job:
                diffexp_job.write(json.dumps(diffexp_json))

def split_bam_file(replicate_bam,threads):
    #print("Splitting %s into %s files for parallel htseq-count"%(replicate_bam,str(threads)))
    #count number of reads in bam file
    count_cmd_list = ["samtools","view","-c","--threads",str(threads),replicate_bam]
    print(" ".join(count_cmd_list))
    count_cmd = subprocess.Popen(count_cmd_list,stdout=subprocess.PIPE)
    num_lines,err = count_cmd.communicate()
    num_lines = float(num_lines.strip())
    #number of lines per file
    num_lines_per_file = int(math.floor(num_lines/float(threads)))
    #Get header from bam file: htseq won't work without the header
    header_cmd = ["samtools","view","-H",replicate_bam]
    print(" ".join(header_cmd))
    header_output = subprocess.Popen(header_cmd,stdout=subprocess.PIPE)
    header,err = header_output.communicate()
    #separate bam file into multiple bam files
    view_cmd = ["samtools","view","--threads",str(threads),replicate_bam]
    split_cmd = ["split","-l",str(num_lines_per_file)]
    if os.path.exists("Split_Bams"):
        #if directory wasn't deleted for some reason, remove before running
        shutil.rmtree("Split_Bams") 
    print(" ".join(view_cmd))
    bam_var = subprocess.Popen(view_cmd,stdout=subprocess.PIPE)
    os.mkdir("Split_Bams")
    os.chdir("Split_Bams")
    print(" ".join(split_cmd))
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

def run_htseq_parallel(genome_annotation,replicate_bam,counts_file,strand,feature,feature_type,threads):
    #run htseq-count
    #TODO: have check to add "-r pos" to command in case of paired end reads
    htseq_cmd = ["htseq-count","-t",feature_type,"-m","intersection-nonempty","--nonunique=all","-f","sam","-s",strand,"-i",feature,"-n",str(threads)]
    for sam in glob.glob("Split_Bams/*"): #iterate through all split files and format
        htseq_cmd.append(sam)
    htseq_cmd.append(genome_annotation)
    print(" ".join(htseq_cmd))
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
def run_htseq_count(genome_list, condition_dict, parameters, job_data, output_dir):
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
                            split_bam_file(replicate[genome_file]["bam"],threads)
                        if os.path.exists(replicate[genome_file]["gene_counts"]):
                            sys.stderr.write("%s exists for genome file %s: skipping htseq-count for genes\n"%(replicate[genome_file]["gene_counts"],genome_file))
                        else:
                            #Run gene_count matrix parameters
                            #Create two counts files and append them together
                            genes_file = os.path.basename(replicate[genome_file]["gene_counts"])
                            genes_file_1 = genes_file+".1" 
                            run_htseq_parallel(genome_annotation,replicate[genome_file]["bam"],genes_file,strand,"ID","gene",threads)
                            run_htseq_parallel(genome_annotation,replicate[genome_file]["bam"],genes_file_1,strand,"ID","pseudogene",threads)
                            cf_open = open(counts_file,"a")
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
                            #First run one feature, then iterate through and append results to output transcripts file 
                            transcript_file = os.path.basename(replicate[genome_file]["transcript_counts"])
                            #run_htseq_parallel(genome_annotation,replicate[genome_file]["bam"],transcript_file,strand,"ID","mRNA",threads)
                            #feature_list = ["mRNA","pseudogene","lnc_RNA","transcript","snRNA","V_gene_segment"]
                            feature_list = ["mRNA","pseudogene","lnc_RNA","transcript","snRNA","V_gene_segment","snoRNA","enhancer","biological_region","primary_transcript","miRNA","C_gene_segment","rRNA","tRNA"]
                            for f in feature_list:
                                transcript_file_tmp = transcript_file+".tmp" 
                                run_htseq_parallel(genome_annotation,replicate[genome_file]["bam"],transcript_file_tmp,strand,"ID",f,threads)
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
                            split_bam_file(replicate[genome_file]["bam"],threads)
                            run_htseq_parallel(genome_annotation,replicate[genome_file]["bam"],counts_file,strand,feature,feature_type,threads)
                    if os.path.exists("Split_Bams"):
                        shutil.rmtree("Split_Bams")
                else:
                    #TODO: add support for host pipeline 
                    if recipe == "Host":
                        print("Host htseq pipeline not enabled for 1 thread, will take too long")
                        sys.exit(1)
                    print("running htseq-count and writing to %s"%counts_file)
                    #TODO: update this command
                    htseq_cmd = ["htseq-count","-t",feature_type,"-f","bam","-s",strand,"-i",feature,replicate[genome_file]["bam"],genome_annotation]
                    print(" ".join(htseq_cmd))
                    #prints to stdout, so redirect output to file
                    with open(counts_file,"w") as cf:
                        subprocess.check_call(htseq_cmd,stdout=cf)
 
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
def create_deseq_metadata(genome_list,condition_dict,output_dir):
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
    metadata_file = os.path.join(output_dir,"DESeq_Metadata.txt")
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
def prep_stringtie_diffexp(genome_list,condition_dict):
    for genome in genome_list:
        genome_file = genome["genome"]
        genome_dir = genome["output"]
        genome_id = os.path.basename(genome_dir)
        os.chdir(genome_dir)
        genome_counts_mtx = genome_id+".stringtie.gene_counts" 
        transcript_counts_mtx = genome_id+".stringtie.transcript_counts"
        genome["gene_matrix"] = os.path.join(genome["output"],genome_counts_mtx)
        genome["transcript_matrix"] = os.path.join(genome["output"],transcript_counts_mtx)
        if not os.path.exists(genome_counts_mtx):
            #TODO: add check for read length
            prep_cmd = ["prepDE.py","-i",genome["prepDE_input"],"-l","101","-g",genome_counts_mtx,"-t",transcript_counts_mtx]
            print(" ".join(prep_cmd))
            subprocess.check_call(prep_cmd)

#Gets the lists of contrasts and runs DESeq2
#Runs DESeq2 once for each genome
def run_deseq2(genome_list,contrasts,job_data):
    #Get list of contrasts to pass into deseq2 R script
    contrast_cmd = []
    for pair in contrasts:
        #remove any commas in names if present
        pair = [x.replace(",","_") for x in pair]
        contrast_cmd.append(",".join(pair))  

    #For each genome, run the deseq2 R script passing in the genome counts matrix, metadata file, and all contrasts
    #changes directory to the top genome directory for each genome
    #invoking RunDESeq2.R: RunDESeq2.R <counts_file.txt> <metadata_file.txt> <output_prefix> <feature_count> <contrast_1> <contrast_2> ... <contrast_n>
    for genome in genome_list:
        os.chdir(genome["output"])
        genome_prefix = os.path.basename(genome["output"])
        diffexp_list = [] #list of files to run differential expression on 
        diffexp_list.append((genome["gene_matrix"],"genes",job_data.get("feature_count","htseq")))
        #transcript_matrix doesn't exist for bacterial pipeline
        if "transcript_matrix" in genome:
            diffexp_list.append((genome["transcript_matrix"],"transcripts",job_data.get("feature_count","htseq")))
        metadata_file = genome["deseq_metadata"] 
        if not os.path.exists(genome["gene_matrix"]):
            print("%s: file doesn't exist"%genome["gene_matrix"])
            continue
        if "transcript_matrix" in genome and not os.path.exists(genome["transcript_matrix"]):
            print("%s: file doesn't exist"%genome["transcript_matrix"])
            continue
        if not os.path.exists(metadata_file):
            print("%s: file doesn't exist"%metadata_file)
            continue
        #diffexp_cmd_gene = ["RunDESeq2.R",genome["gene_matrix"],metadata_file,genome_prefix+".genes"]+contrast_cmd
        #diffexp_cmd_transcript = ["RunDESeq2.R",genome["transcript_matrix"],metadata_file,genome_prefix+".transcripts"]+contrast_cmd
        for diffexp_params in diffexp_list:
            diffexp_cmd = ["RunDESeq2.R",diffexp_params[0],metadata_file,diffexp_params[1],diffexp_params[2]]+contrast_cmd
            print("%s\n"%" ".join(diffexp_cmd))
            subprocess.check_call(diffexp_cmd)
            genome["diff_exp_contrasts"] = []
            for pair in contrasts:
                pair = [x.replace(",","_") for x in pair]
                pair = "_vs_".join(pair) 
                #diffexp_file = genome_prefix + "_" + pair + ".txt"
                diffexp_output = pair+"."+diffexp_params[2]+"."+diffexp_params[1]+".deseq2"
                genome["diff_exp_contrasts"].append(os.path.join(genome["output"],diffexp_output))
         
#Writes the gene_exp.gmx file used in expression_transform.py from DESeq2 output
def write_gmx_file(genome_list):
    gene_set = set()
    gene_count_dict = {}
    contrast_list = []
    for genome in genome_list:
        os.chdir(genome["output"])
        contrast_file_list = genome["diff_exp_contrasts"]
        for contrast_file in contrast_file_list:
            contrast_name = os.path.basename(contrast_file).replace(".txt","")
            contrast_list.append(contrast_name)
            gene_count_dict[contrast_name] = {}
            with open(contrast_file,"r") as cf:
                next(cf)
                for line in cf:
                    gene,baseMean,log2FC,lfcSE,stat,pvalue,padj = line.strip().split("\t")
                    gene_set.add(gene)
                    gene_count_dict[contrast_name][gene] = log2FC

        with open("gene_exp.gmx","w") as o:
            o.write("Gene_ID\t%s\n"%"\t".join(contrast_list))
            for gene in gene_set:
                o.write(gene)
                for contrast in contrast_list:
                    if gene in gene_count_dict[contrast]:
                        o.write("\t%s"%gene_count_dict[contrast][gene])
                    else:
                        o.write("\t0")
                o.write("\n")

def setup(genome_list, condition_dict, parameters, output_dir, job_data):
    for genome in genome_list:
        genome_link=os.path.join(output_dir, os.path.basename(genome["genome"]))
        genome["genome_link"]=genome_link
        genome["host"] = ("hisat_index" in genome and genome["hisat_index"])
        if not os.path.exists(genome_link):
            subprocess.check_call(["ln","-s",genome["genome"],genome_link])
        make_directory_names(genome, condition_dict)
        for condition in condition_dict:
            rcount=0
            for r in condition_dict[condition]["replicates"]:
                target_dir=r["target_dir"]
                rcount+=1
                subprocess.call(["mkdir","-p",target_dir])
                if "srr_accession" in r:
                    srr_id = r["srr_accession"] 
                    meta_file = os.path.join(target_dir,srr_id+"_meta.txt")
                    subprocess.check_call(["p3-sra","--out",target_dir,"--metadata-file", meta_file, "--id",srr_id])
                    with open(meta_file) as f:
                        job_meta = json.load(f)
                        files = job_meta[0].get("files",[])
                        for i,f in enumerate(files):
                            if f.endswith("_2.fastq"):
                                r["read2"]=os.path.join(target_dir, f)
                            elif f.endswith("_1.fastq"):
                                r["read1"]=os.path.join(target_dir, f)
                            elif f.endswith(".fastq"):
                                r["read1"]=os.path.join(target_dir,f)
                            elif f.endswith("fastqc.html"):
                                r["fastqc"].append(os.path.join(target_dir, f))
                    






def main(genome_list, condition_dict, parameters_str, output_dir, gene_matrix=False, contrasts=[], job_data=None, map_args=None, diffexp_json=None):
    #arguments:
    #list of genomes [{"genome":somefile,"annotation":somefile}]
    #dictionary of library dictionaries structured as {libraryname:{library:libraryname, replicates:[{read1:read1file, read2:read2file}]}}
    #parametrs_str is json parameters list keyed as bowtie, cufflinks, cuffdiff.
    output_dir=os.path.abspath(output_dir)
    subprocess.call(["mkdir","-p",output_dir])
    #TODO: Read objects(?) that know their read parameters:
    #TODO: fastq utils to replace certain functionality?
    if parameters_str != None :
        parameters=json.loads(parameters_str)
    else:
        parameters = {}
    setup(genome_list, condition_dict, parameters, output_dir, job_data)

    #TRUE: runs cufflinks then cuffdiff if differential expression is turned on
    #FALSE: runs either htseq-count or stringtie
    run_cuffdiff_pipeline = job_data.get("feature_count","htseq") == "cuffdiff"
    run_alignment(genome_list, condition_dict, parameters, output_dir, job_data)
    if run_cuffdiff_pipeline:
        run_cufflinks(genome_list, condition_dict, parameters, output_dir)
    else:
        run_featurecount(genome_list, condition_dict, parameters, output_dir, job_data)
    #TODO: change novel_features condition when novel-isoform differential expression is implemented
    if len(condition_dict.keys()) > 1 and not job_data.get("novel_features",False):
        #If running cuffdiff pipeline, terminated after running expression import
        if run_cuffdiff_pipeline:
            run_cuffdiff(genome_list, condition_dict, parameters, output_dir, gene_matrix, contrasts, job_data, map_args, diffexp_json)
            run_diff_exp_import(genome_list, condition_dict, parameters, output_dir, contrasts, job_data, map_args, diffexp_json)
            sys.exit(0)
        #Create deseq2 metadata file. Ordering the matrix correctly will be done in R script 
        create_deseq_metadata(genome_list,condition_dict,output_dir)
        #run prepDE.py to create counts files if using stringtie
        if job_data.get("feature_count","htseq") == "stringtie":
            write_gtf_list(genome_list,condition_dict) #function that writes the input for prepDE.py, which is a list of samples and paths to their gtf files. Do this for each genome
            prep_stringtie_diffexp(genome_list,condition_dict)   
        else: #htseq
            if job_data.get("recipe","RNA-Rocket") == "Host":
                create_counts_table_host(genome_list,condition_dict,job_data)
            else:
                create_counts_table(genome_list,condition_dict,job_data)
        run_deseq2(genome_list,contrasts,job_data)
        write_gmx_file(genome_list)
        #differential expression import will fail for host
        #if not job_data.get("recipe","RNA-Rocket") == "Host":
        run_diff_exp_import(genome_list, condition_dict, parameters, output_dir, contrasts, job_data, map_args, diffexp_json)
        ###TODO: Remove
        #with open("Variable_Outputs.txt","w") as o:
            #o.write("genome_list:\n")
            #o.write(genome_list)
            #o.write("\n\ncondition_dict:\n")
            #o.write(condition_dict)
            #o.write("\n\nparameters:\n")
            #o.write(parameters)
            #o.write("\n\ncontrasts:\n")
            #o.write(contrasts)
            #o.write("\n\njob_data:\n")
            #o.write(job_data)
            #o.write("\n\nmap_args:\n")
            #o.write(map_args)
            #o.write("\n\ndiffexp_json:\n")
            #o.write(diffexp_json)
        ###

if __name__ == "__main__":
    #modelling input parameters after rockhopper
    parser = argparse.ArgumentParser()
    #if you want to support multiple genomes for alignment you should make this json payload an nargs+ parameter
    parser.add_argument('--jfile',
            help='json file for job {"reference_genome_id": "1310806.3", "experimental_conditions":\
                    ["c1_control", "c2_treatment"], "output_file": "rnaseq_baumanii_1505311", \
                    "recipe": "RNA-Rocket", "output_path": "/anwarren@patricbrc.org/home/test",\
                    "paired_end_libs": [{"read1": "/anwarren@patricbrc.org/home/rnaseq_test/MHB_R1.fq.gz",\
                    "read2": "/anwarren@patricbrc.org/home/rnaseq_test/MHB_R2.fq.gz", "condition": 1},\
                    {"read1": "/anwarren@patricbrc.org/home/rnaseq_test/MERO_75_R1.fq.gz",\
                    "read2": "/anwarren@patricbrc.org/home/rnaseq_test/MERO_75_R2.fq.gz", "condition": 2}], "contrasts": [[1, 2]]}', required=True)
    parser.add_argument('--sstring', help='json server string specifying api {"data_api":"url"}', required=False, default=None)
    parser.add_argument('-g', help='csv list of directories each containing a genome file *.fna and annotation *.gff', required=True)
    parser.add_argument('--index', help='flag for enabling using HISAT2 indices', action='store_true', required=False)
    #parser.add_argument('-L', help='csv list of library names for comparison', required=False)
    #parser.add_argument('-C', help='csv list of comparisons. comparisons are library names separated by percent. ', required=False)
    parser.add_argument('-p', help='JSON formatted parameter list for tuxedo suite keyed to program', default="{}", required=False)
    parser.add_argument('-o', help='output directory. defaults to current directory.', required=False, default=None)
    parser.add_argument('-d', help='name of the folder for differential expression job folder where files go', required=True) 
    #parser.add_argument('-x', action="store_true", help='run the gene matrix conversion and create a patric expression object', required=False)
    #parser.add_argument('readfiles', nargs='+', help="whitespace sep list of read files. shoudld be \
    #        in corresponding order as library list. ws separates libraries,\
    #        a comma separates replicates, and a percent separates pairs.")
    if len(sys.argv) ==1:
        parser.print_help()
        sys.exit(2)
    map_args = parser.parse_args()
    assert map_args.d.startswith(".") # job object folder name needs a .
    condition_dict={}
    condition_list=[]
    comparison_list=[]
    #if map_args.L:
    #    condition_list=map_args.L.strip().split(',')
    #if args.C:
    #    comparison_list=[i.split("%") for i in args.C.strip().split(',')]
        
    #create library dict
    with open(map_args.jfile, 'r') as job_handle:
        job_data = json.load(job_handle)

    condition_list= job_data.get("experimental_conditions",[])
    got_conditions=False
    if not len(condition_list):
        condition_list.append("results")
    else:
        got_conditions=True
    gene_matrix=True
    if map_args.o == None:
        output_dir="./"
    else:
        output_dir=map_args.o
    for cond in condition_list:
        condition_dict[cond]={"condition":cond}
    count=0
    #add read/replicate structure to library dict
    replicates={}
    contrasts=[]
    for c in job_data.get("contrasts",[]):
        contrasts.append([condition_list[c[0]-1],condition_list[c[1]-1]])
    for read in job_data.get("paired_end_libs",[])+job_data.get("single_end_libs",[])+job_data.get("srr_libs",[]):
        if "read" in read:
            read["read1"] = read.pop("read")
        condition_index = int(read.get("condition", count+1))-1 #if no condition return position so everything is diff condition
        condition_id = condition_list[condition_index] if got_conditions else "results"
        #organize things by condition/replicate
        condition_dict[condition_id].setdefault("replicates",[]).append(read)
        #store the condition index if it doesn't exist
        condition_dict[condition_id].setdefault("condition_index",condition_index)
        count+=1
    genome_dirs=map_args.g.strip().split(',')
    genome_list=[]
    for g in genome_dirs:
        cur_genome={"genome":[],"annotation":[],"dir":g,"hisat_index":[]}
        for f in os.listdir(g):
            if f.endswith(".fna") or f.endswith(".fa") or f.endswith(".fasta"):
                cur_genome["genome"].append(os.path.abspath(os.path.join(g,f)))
            elif f.endswith(".gff"):
                cur_genome["annotation"].append(os.path.abspath(os.path.join(g,f)))
            #Change to check for ht2 files instead of just the tar file
            elif f.endswith(".ht2.tar") and len(cur_genome["hisat_index"]) == 0: 
                cur_genome["hisat_index"].append(os.path.abspath(os.path.join(g,f)))
            elif f.endswith(".1.ht2"):
                cur_genome["hisat_index"].append(os.path.abspath(os.path.join(g,f)))

        if len(cur_genome["genome"]) != 1:
            sys.stderr.write("Too many or too few fasta files present in "+g+"\n")
            sys.exit(2)
        else:
            cur_genome["genome"]=cur_genome["genome"][0]
        if len(cur_genome["annotation"]) != 1:
            sys.stderr.write("Too many or too few gff files present in "+g+"\n")
            sys.exit(2)
        else:
            cur_genome["annotation"]=cur_genome["annotation"][0]
        if map_args.index:
            if len(cur_genome["hisat_index"]) != 1:
                sys.stderr.write("Missing hisat index tar file for "+g+"\n")
                sys.exit(2)
            else:
                cur_genome["hisat_index"]=cur_genome["hisat_index"][0]
        genome_list.append(cur_genome)

    #job template for differential expression object
    diffexp_json = json.loads("""                    {
                        "app": {
                            "description": "Parses and transforms users differential expression data",
                            "id": "DifferentialExpression",
                            "label": "Transform expression data",
                            "parameters": [
                                {
                                    "default": null,
                                    "desc": "Comparison values between samples",
                                    "id": "xfile",
                                    "label": "Experiment Data File",
                                    "required": 1,
                                    "type": "wstype",
                                    "wstype": "ExpList"
                                },
                                {
                                    "default": null,
                                    "desc": "Metadata template filled out by the user",
                                    "id": "mfile",
                                    "label": "Metadata File",
                                    "required": 0,
                                    "type": "wstype",
                                    "wstype": "ExpMetadata"
                                },
                                {
                                    "default": null,
                                    "desc": "User information (JSON string)",
                                    "id": "ustring",
                                    "label": "User string",
                                    "required": 1,
                                    "type": "string"
                                },
                                {
                                    "default": null,
                                    "desc": "Path to which the output will be written. Defaults to the directory containing the input data. ",
                                    "id": "output_path",
                                    "label": "Output Folder",
                                    "required": 0,
                                    "type": "folder"
                                },
                                {
                                    "default": null,
                                    "desc": "Basename for the generated output files. Defaults to the basename of the input data.",
                                    "id": "output_file",
                                    "label": "File Basename",
                                    "required": 0,
                                    "type": "wsid"
                                }
                            ],
                            "script": "App-DifferentialExpression"
                        },
                        "elapsed_time": null,
                        "end_time": null,
                        "hostname": "",
                        "id": "",
                        "is_folder": 0,
                        "job_output": "",
                        "output_files": [
                        ],
                        "parameters": {
                            "mfile": "",
                            "output_file": "",
                            "output_path": "",
                            "ustring": "",
                            "xfile": ""
                        },
                        "start_time": ""
                    }
""")
    main(genome_list,condition_dict,map_args.p,output_dir,gene_matrix,contrasts,job_data,map_args,diffexp_json)


