#!/usr/bin/env python

import os,sys,glob,subprocess
import concurrent.futures
import shutil
import tarfile
import json,gzip

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

#general function for linking a file to a target directory
def link_file(target_dir,filename):
    outfile = os.path.join(target_dir,os.path.basename(filename))
    if not os.path.exists(outfile):
        subprocess.check_call(["ln","-s",filename,outfile])
    return outfile

#pretty simple: its for prokaryotes in that parameters will be attuned to give best performance and no tophat
def run_alignment(genome_list, condition_dict, parameters, output_dir, job_data, pipeline_log): 
    #modifies condition_dict sub replicates to include 'bowtie' dict recording output files
    for genome in genome_list:
        genome_link = genome["genome_link"]
        #samstat dict for putting the samstat reports in the same portion in the multiqc report
        samstat_dict = {}
        samstat_dict["reports"] = []
        final_cleanup=[]
        if "hisat_index" in genome and genome["hisat_index"]:
            if type(genome["hisat_index"]) is list:
                genome["hisat_index"] = genome["hisat_index"][0] 
            archive = tarfile.open(genome["hisat_index"])
            indices = [os.path.join(output_dir,os.path.basename(x)) for x in archive.getnames()]
            final_cleanup+=indices
            num_parts = archive.getnames()[0].count("/") #strip any subdirectory structure from the tar output
            archive.close()
            tar_cmd = ["tar","-xvf",genome["hisat_index"],"-C",output_dir]
            if num_parts > 0:
                tar_cmd += ["--strip-components",str(num_parts)]
            print (" ".join(tar_cmd))
            pipeline_log.append(" ".join(tar_cmd))
            subprocess.check_call(tar_cmd)
            index_prefix = os.path.join(output_dir, os.path.basename(genome["hisat_index"]).replace(".ht2.tar","")) #somewhat fragile convention. tar prefix is underlying index prefix
            print("index_prefix = %s"%index_prefix)
            genome["index_prefix"] = index_prefix
            #for index in glob.glob(index_prefix+"*.ht2"):
            #    print(index)
            if not os.path.exists(index_prefix+".1.ht2"):
                print("hisant indices were not unpacked correctly: %s"%index_prefix)
                #for h_index in glob.glob(os.path.dirname(genome["hisat_index"])+"*ht2"):
                #    os.symlink(h_index,os.path.join(output_dir,os.path.basename(h_index)))
                sys.exit(1)
            cmd=["hisat2","--dta-cufflinks", "-x", index_prefix] 
            #cmd=["hisat2","--dta-cufflinks", "-x", index_prefix,"--no-spliced-alignment"] 
            thread_count= parameters.get("hisat2",{}).get("-p",0)
        else:
            #cmd=["hisat2","--dta-cufflinks", "-x", genome_link, "--no-spliced-alignment"] 
            try:
                bowtie_build_cmd = ["bowtie2-build", genome_link, genome_link]
                print(" ".join(bowtie_build_cmd))
                pipeline_log.append(" ".join(bowtie_build_cmd))
                subprocess.check_call(["bowtie2-build", genome_link, genome_link])
            except Exception as err:
                sys.stderr.write("bowtie build failed: %s %s\n"%(err, genome_link))
                sys.exit(1)
            cmd=["bowtie2", "-x", genome_link]
            thread_count= parameters.get("bowtie2",{}).get("-p",0)
        ###
        # strandedness determination requires a bam file, so inserting here since the above creates the hisat/bowtie indices
        # get strandedness parameter: stores the value in the replicates dictionary. Key does not exist if single end
        # return if there is an issue in strandedness
        for condition in condition_dict:
            for r in condition_dict[condition]["replicates"]:
                if genome["genome"] not in r:
                    r[genome["genome"]] = {}
        ret_val = run_fastqc(genome,condition_dict,parameters,pipeline_log)
        if ret_val != 0:
            print("Error in fastqc")
            return -1
        if not job_data.get("skip_sampling",False): #assigning strandedness requires sampling
            ret_val = run_sample_alignment(genome,condition_dict,parameters,pipeline_log)
            if ret_val != 0:
                print("Error in running sample alignment")
                return -1
            ret_val = check_sample_alignment(genome,condition_dict,parameters)
            if ret_val != 0:
                print("Error from checking sample results")
                return -1
            ret_val = assign_strandedness_parameter(genome,condition_dict,parameters,pipeline_log)
            if ret_val != 0:
                print("Error in assigning strandedness")
                return -1
        ###
        if thread_count == 0:
            thread_count=2 #multiprocessing.cpu_count()
        cmd+=["-p",str(thread_count)]
        scount = 0 #used in modifying labels for the samstat reports
        #loop through bam files present in condition dict  
        for condition in condition_dict:
            for r in condition_dict[condition]["replicates"]:
                if "bam" not in r:
                    continue
                scount+=1
                bam_file = link_file(r["target_dir"],r["bam"])
                r[genome["genome"]]["bam"] = bam_file
                samstat_cmd = ["samstat",bam_file]
                samtools_threads = parameters.get("samtools",{}).get("-@",1)
                stats_cmd = ["samtools","stats","--threads",str(samtools_threads),bam_file]
                stats_outfile = os.path.join(r["target_dir"],os.path.basename(bam_file.replace("bam","samtools_stats")))
                if not os.path.exists(stats_outfile):
                    pipeline_log.append(" ".join(stats_cmd))
                    try:
                        with open(stats_outfile,"w") as o:
                            subprocess.check_call(stats_cmd,stdout=o)
                    except Exception as e:
                        sys.stderr.write("ERROR in samtools stats:\n{0}\n".format(e))
                        return(-3)
                try:
                    r[genome["genome"]]["avg_read_length"] = get_average_read_length_per_file(stats_outfile)
                except Exception as e:
                    sys.stderr.write("ERROR in reading from samtools stats output file\n")  
                    return (-3)
                samstat_file = os.path.join(r["target_dir"],os.path.basename(bam_file)+".samstat.html")
                mod_samstat_file = os.path.join(os.path.dirname(bam_file),"Samstat_"+os.path.basename(bam_file)+".samstat.html")
                print(mod_samstat_file)
                r[genome["genome"]]["samstat"] = mod_samstat_file
                samstat_dict["reports"].append(mod_samstat_file)
                if not os.path.exists(samstat_file):
                    pipeline_log.append(" ".join(samstat_cmd))
                    try:
                        subprocess.check_call(samstat_cmd)
                    except Exception as e:
                        sys.stderr.write("ERROR running samstat on {0}:\n{1}\n".format(bam_file,e))
                        return(-3)
                #if not os.path.exists(mod_samstat_file):
                modify_samstat_for_multiqc(samstat_file,scount)
        #loop through reads
        for condition in condition_dict:
            rcount=0
            for r in condition_dict[condition]["replicates"]:
                #skip this for loop if bam file is present: processed in the previous for loop 
                if "bam" in r:
                    continue
                scount+=1
                cur_cleanup=[]
                rcount+=1
                target_dir=r["target_dir"]
                samstat_cmd=["samstat"]
                cur_cmd=list(cmd)
                if "read2" in r:
                    cur_cmd+=["-1",link_space(r["read1"]),"-2",link_space(r["read2"])]
                    #name1=os.path.splitext(os.path.basename(r["read1"]))[0].replace(" ","")
                    #name2=os.path.splitext(os.path.basename(r["read2"]))[0].replace(" ","")
                    #sam_file=os.path.join(target_dir,name1+"_"+name2+".sam")
                else:
                    cur_cmd+=[" -U",link_space(r["read1"])]
                    #name1=os.path.splitext(os.path.basename(r["read1"]))[0].replace(" ","")
                    #sam_file=os.path.join(target_dir,name1+".sam")
                sam_file=os.path.join(target_dir,r["name"]+".sam")
                cur_cleanup.append(sam_file)
                bam_file=sam_file[:-4]+".bam"
                samstat_cmd.append(bam_file)
                r[genome["genome"]]={}
                r[genome["genome"]]["bam"]=bam_file
                cur_cmd+=["-S",sam_file]
                #add strandedness parameter
                if "strand_param" in r and r["strand_param"] != "undetermined":
                    if cur_cmd[0] == "hisat2":
                        cur_cmd.insert(1,"--rna-strandness")
                        cur_cmd.insert(2,r["strand_param"])
                    else:
                        cur_cmd.insert(1,"--"+r["strand_param"].lower()) #strand_param is uppercase by default, but bowtie requires lowercase
                if os.path.exists(bam_file):
                    sys.stderr.write(bam_file+" alignments file already exists. skipping\n")
                else:
                    print (" ".join(cur_cmd))
                    if job_data.get("recipe","RNA-Rocket") == "Host":
                        alignment_log = bam_file.replace("bam","hisat") 
                    else:
                        alignment_log = bam_file.replace("bam","bowtie") 
                    with open(alignment_log,"w") as al:
                        pipeline_log.append(" ".join(cur_cmd))
                        subprocess.check_call(cur_cmd,stdout=al,stderr=al) #call bowtie2 or hisat2
                samtools_threads = parameters.get("samtools",{}).get("-@",1)
                if not os.path.exists(bam_file):
                    pipeline_log.append("samtools view -Su "+sam_file+" | samtools sort -o - - -@ "+str(samtools_threads)+" > "+bam_file)
                    subprocess.check_call("samtools view -Su "+sam_file+" | samtools sort -o - - -@ "+str(samtools_threads)+" > "+bam_file, shell=True)#convert to bam
                    pipeline_log.append("samtools index "+bam_file)
                    subprocess.check_call("samtools index "+bam_file, shell=True)
                    #subprocess.check_call('samtools view -S -b %s > %s' % (sam_file, bam_file+".tmp"), shell=True)
                    #subprocess.check_call('samtools sort %s %s' % (bam_file+".tmp", bam_file), shell=True)
                print (" ".join(samstat_cmd))
                stats_cmd = ["samtools","stats","--threads",str(samtools_threads),bam_file]
                stats_outfile = bam_file.replace("bam","samtools_stats")
                if not os.path.exists(stats_outfile):
                    pipeline_log.append(" ".join(stats_cmd))
                    with open(stats_outfile,"w") as o:
                        subprocess.check_call(stats_cmd,stdout=o)
                r[genome["genome"]]["avg_read_length"] = get_average_read_length_per_file(stats_outfile)
                samstat_file = bam_file+".samstat.html"
                #mod_samstat_file = os.path.join(os.path.dirname(bam_file),"Samstat_"+os.path.basename(bam_file)+".samstat_mqc.html")
                mod_samstat_file = os.path.join(os.path.dirname(bam_file),"Samstat_"+os.path.basename(bam_file)+".samstat.html")
                print(mod_samstat_file)
                r[genome["genome"]]["samstat"] = mod_samstat_file
                samstat_dict["reports"].append(mod_samstat_file)
                if not os.path.exists(samstat_file):
                    pipeline_log.append(" ".join(samstat_cmd))
                    subprocess.check_call(samstat_cmd)
                #if not os.path.exists(mod_samstat_file):
                try:
                    modify_samstat_for_multiqc(samstat_file,scount)
                except Exception as e:
                    sys.stderr.write("ERROR modifying samstat report for multiqc, {0}:\n{1}\n".format(samstat_file,e))
                    return(-3)
                for garbage in cur_cleanup:
                    subprocess.call(["rm", garbage])
        ###write out samstat json file
        os.chdir(genome["output"])
        with open("samstat.json","w") as so:
            so.write(json.dumps(samstat_dict))
        ###cleanup files
        for garbage in final_cleanup:
            subprocess.call(["rm", garbage])
    return 0

def run_fastqc(genome,condition_dict,parameters,pipeline_log):
    fastqc_cmd_list = []
    for condition in condition_dict:
        for r in condition_dict[condition]["replicates"]:
            target_dir=r["target_dir"]
            fastqc_cmd=["fastqc","--outdir",target_dir]
            if "read2" in r:
                fastqc_cmd+=[r["read1"],r["read2"]]
                r[genome["genome"]]["fastqc"] = r["name"]+"_1_fastqc.html" 
            else:
                fastqc_cmd+=[r["read1"]]
                r[genome["genome"]]["fastqc"] = r["name"]+"_fastqc.html" 
            if not os.path.exists(r[genome["genome"]]["fastqc"]):
                print (" ".join(fastqc_cmd))
                pipeline_log.append(" ".join(fastqc_cmd))
                fastqc_cmd_list.append(fastqc_cmd)
    future_returns = []
    with concurrent.futures.ThreadPoolExecutor(max_workers=8) as pool:
        future_returns = list(pool.map(run_fastqc_pool,fastqc_cmd_list))
    for f in future_returns:
        if f != 0:
            return -1
    return 0

def run_fastqc_pool(job):
    print(" ".join(job))
    try:
        subprocess.check_call(job)
    except Exception as e:
        sys.stderr.write("ERROR in fastqc for sample {0}:\n{1}\n".format(r["name"]," ".join(fastqc_cmd)))
        return -1
    return 0

#Only applies to paired-end data
def assign_strandedness_parameter(genome,condition_dict,parameters,pipeline_log):
    #generate bed files for genome
    ret_val = generate_bed_from_gff(genome)
    if ret_val != 0:
        return(-3)
    remove_list = []
    num_sample = 1000
    for condition in condition_dict:
        for r in condition_dict[condition]["replicates"]:
            #skip if not paired end reads
            if "read2" not in r:
                continue
            #skip if already determined somehow: probably only likely to happend during dualRNASeq, which is not implemented yet
            if "strand_param" in r:
                continue
            #skip if bam file is present
            if "bam" in r:
                continue
            #add sample sam file and alignment statistics to cleanup
            remove_list.append(r[genome["genome"]]["sample_sam_file"])
            remove_list.append(r[genome["genome"]]["sample_alignment_output"])
            ###generate RSEQc strandedness file using infer_experiment.py
            infer_cmd = ["infer_experiment.py","-i",r[genome["genome"]]["sample_sam_file"],"-r",genome["bed"],"-s",str(num_sample)]
            infer_file = os.path.join(r["target_dir"],r["name"]+".infer") 
            print(" ".join(infer_cmd))
            pipeline_log.append(" ".join(infer_cmd))
            with open(infer_file,"w") as o:
                try:
                    subprocess.check_call(infer_cmd,stdout=o)
                except Exception as e:
                    sys.stderr.write("ERROR in infer_experiment for sampled file {0}:\n{1}\n".format(r[genome["genome"]]["sample_sam_file"],e))
                    return(-3)
            ###assess infer_file and set strandedness parameter for replicate
            with open(infer_file,"r") as i:
                for x in range(0,3): #skip first 3 lines
                    next(i)    
                undetermined = float(next(i).split(":")[1].strip()) 
                rf_stranded = float(next(i).split(":")[1].strip())
                fr_stranded = float(next(i).split(":")[1].strip())
                if undetermined > fr_stranded and undetermined > rf_stranded:
                    r["strand_param"] = "undetermined"
                elif abs(fr_stranded - rf_stranded) <= 0.1: #subjective threshold: should work in most cases
                    r["strand_param"] = "undetermined"
                else:
                    r["strand_param"] = "FR" if fr_stranded > rf_stranded else "RF"
    ###remove all files in the remove_list
    for trash in remove_list:
        rm_cmd = ["rm",trash]
        subprocess.check_call(rm_cmd)
    ###output library geometry file for user
    os.chdir(genome["output"])
    with open("library_geometry.txt","w") as geom_handle:
        geom_handle.write("Sample\tCondition\tStrandedness\n")
        for condition in condition_dict:
            for r in condition_dict[condition]["replicates"]:
                if "strand_param" in r: 
                    geom_handle.write("%s\t%s\t%s\n"%(r["name"],condition,r["strand_param"]))
                else:
                    geom_handle.write("%s\t%s\t%s\n"%(r["name"],condition,"NA"))
    return 0
                    
def run_sample_alignment(genome,condition_dict,parameters,pipeline_log):
    remove_list = []
    for condition in condition_dict:
        for r in condition_dict[condition]["replicates"]:
            #if the number of reads is less than 1000, skip file:
            num_sample = 2000
            min_reads_flag = has_min_reads(num_sample,r["read1"])
            if not min_reads_flag:
                continue
            #sample fastq files using seqtk
            #setting defulat sampling to 1000 reads
            #-s parameter sets the seed and ensures the same reads are sampled from each file, assuming the files are properly paired
            ###TODO: replace shell=True for redirecting the output to a file
            ###TODO: for the moment, termination run if any of the alignment steps outright fails
            sample1_file = os.path.join(r["target_dir"],r["name"] + ".s1.fastq")
            sample1_cmd = " ".join(["seqtk","sample","-s","42",r["read1"],str(num_sample),">",sample1_file])
            remove_list.append(sample1_file)
            print(sample1_cmd)
            pipeline_log.append(" ".join(sample1_cmd))
            if not os.path.exists(sample1_file):
                try:
                    subprocess.check_call(sample1_cmd,shell=True)
                except Exception as e:
                    sys.stderr.write("ERROR in seqtk sampling of {0}:\n{1}\n".format(r["read1"],e))
                    return(-3)
            if "read2" in r:
                sample2_file = os.path.join(r["target_dir"],r["name"] + ".s2.fastq")
                sample2_cmd = " ".join(["seqtk","sample","-s","42",r["read2"],str(num_sample),">",sample2_file])
                remove_list.append(sample2_file)
                print(sample2_cmd)
                pipeline_log.append(" ".join(sample2_cmd))
                if not os.path.exists(sample2_file):
                    try:
                        subprocess.check_call(sample2_cmd,shell=True)
                    except Exception as e:
                        sys.stderr.write("ERROR in seqtk sampling of {0}:\n{1}\n".format(r["read2"],e))
                        return(-3)
            ###run bowtie2 or hisat2 on sampled files
            sample_align_cmd = []
            sam_file = os.path.join(r["target_dir"],r["name"] + ".sample.sam")
            r[genome["genome"]]["sample_sam_file"] = sam_file
            ###TODO: replace hisat_index condition implementation with something more robust
            if len(genome["hisat_index"]) > 0:
                r[genome["genome"]]["sample_alignment_output"] = sam_file.replace("sample.sam","sample.hisat")
                sample_align_cmd = ["hisat2","-x",genome["index_prefix"]]
                if "read2" in r:
                    sample_align_cmd+=["-1",sample1_file,"-2",sample2_file]
                else:
                    sample_align_cmd+=["-U",sample1_file]
                sample_align_cmd+=["--mp","1,0","--pen-noncansplice","20","-S",sam_file,"-p",str(parameters.get("hisat2",{}).get("-p","1"))] 
            else: #bowtie
                sample_align_cmd = ["bowtie2","-x",genome["genome_link"]]
                if "read2" in r:
                    sample_align_cmd+=["-1",sample1_file,"-2",sample2_file]
                else:
                    sample_align_cmd+=["-U",sample1_file]
                sample_align_cmd+=["-S",sam_file,"-p",str(parameters.get("bowtie2",{}).get("-p","1"))]
                r[genome["genome"]]["sample_alignment_output"] = sam_file.replace("sample.sam","sample.bowtie")
            if not os.path.exists(sam_file):
                print(" ".join(sample_align_cmd))
                pipeline_log.append(" ".join(sample_align_cmd))
                try:
                    with open(r[genome["genome"]]["sample_alignment_output"],"w") as sample_out: 
                        subprocess.check_call(sample_align_cmd,stderr=sample_out)
                    with open(r[genome["genome"]]["sample_alignment_output"],"r") as sample_out:
                        sys.stdout.write("{}\n".format("".join(sample_out.readlines())))
                except Exception as e:
                    sys.stderr.write("ERROR in sampling alignment {0}:\n{1}".format(" ".join(sample_align_cmd),e))
                    return(-3)
    ###remove all files in the remove_list
    for trash in remove_list:
        rm_cmd = ["rm",trash]
        subprocess.check_call(rm_cmd)
    return(0)

def check_sample_alignment(genome,condition_dict,parameters):
    error_threshold = 5.0 
    error_list = []
    error_alignment = []
    warning_threshold = 50.0
    warning_list = []
    warning_alignment = []
    for condition in condition_dict:
        for r in condition_dict[condition]["replicates"]:
            if not os.path.exists(r[genome["genome"]]["sample_alignment_output"]):
                print("%s does not exist"%r[genome["genome"]]["sample_alignment_output"])
            with open(r[genome["genome"]]["sample_alignment_output"],"r") as sample_out:
                alignment_output = sample_out.readlines()
            ###Alignment output is as a percentage
            #TODO: make sure hisat outputs in the same format as bowtie
            alignment_value = float(alignment_output[-1].split()[0].replace("%","")) 
            if alignment_value <= error_threshold:
                error_list.append(r["name"])
                error_alignment.append("".join(alignment_output))
            elif alignment_value <= warning_threshold:
                warning_list.append(r["name"])
                warning_alignment.append("".join(alignment_output))
    os.chdir(genome["output"])
    if len(error_list) > 0:
        with open("SAMPLE_ALIGNMENT_ERRORS.txt","w") as o:
            o.write("The following samples contained poor alignment to reference genome_id {0}.\n".format(genome["genome"]))
            o.write("SOLUTION: Check the fastqc output for errors or consider choosing a different reference genome:\n") 
            for idx,sample in enumerate(error_list): 
                o.write("{}\n".format("".join(["-"]*10)))
                o.write("Sample {0}:\n{1}\n".format(sample,error_alignment[idx])) 
    if len(warning_list) > 0:
        with open("SAMPLE_ALIGNMENT_WARNINGS.txt","w") as o:
            o.write("The following samples contained moderate alignment to reference genome_id {0}.\n".format(genome["genome"]))
            o.write("SOLUTION: No action is necessary, but it is recommended to review sample quality or choose a different reference genome to improve results:\n")
            for idx,sample in enumerate(warning_list):
                o.write("{}\n".format("".join(["-"]*10)))
                o.write("Sample {0}:\n{1}\n".format(sample,warning_alignment[idx]))
    if len(error_list) > 0:
        return (-3)
    else:
        return (0)

def generate_bed_from_gff(genome):
    bed_file = genome["annotation_link"].replace(".gff",".bed")
    genome["bed"] = bed_file
    bed_cmd = "gff2bed --do-not-sort < " + genome["annotation_link"] + " > " + bed_file
    try:
        if not os.path.exists(bed_file):
            subprocess.check_call(bed_cmd,shell=True)
    except Exception as e:
        sys.stderr.write("ERROR generating bed {0} file from gff file {1}:\n{2}\n".format(bed_file,genome["annotation_link"],e)) 
        return(-3)
    return 0

#Reads the output from samtools stat and grabs the average read length value
def get_average_read_length_per_file(stats_file):
    with open(stats_file,"r") as sf:
        for line in sf:
            if "average length:" in line: 
                line = line.strip().split()
                avg_length = line[-1]
                return avg_length

def has_min_reads(num_reads,reads_file):
    if not reads_file.endswith(".gz"):
        r_count = 0
        with open(reads_file,"r") as r1:
            for line in r1:
                r_count = r_count + 1  
                for i in range(0,3):
                    next(r1)
                if r_count == num_reads:
                    return True
    else:
        r_count = 0
        with gzip.GzipFile(reads_file) as r1:
            for line in r1:
                r_count = r_count + 1  
                for i in range(0,3):
                    next(r1)
                if r_count == num_reads:
                        return True
    return False

#Pass in a samstat html file to remove portion that conflicts with the multiqc report
#Uses sed to remove the block of css code that clashes with multiqc
#Then goes through and removes all content other than the intro statistics and base quality distribution by position
def modify_samstat_for_multiqc(samstat_filename,count):
    sed_cmd = ["sed","-i","/body {/,/}/ d;",samstat_filename]
    subprocess.check_call(sed_cmd)
    #out_samstat = os.path.join(os.path.dirname(samstat_filename),"Samstat_"+os.path.basename(samstat_filename.replace(".html","_mqc.html")))    
    out_samstat = os.path.join(os.path.dirname(samstat_filename),"Samstat_"+os.path.basename(samstat_filename))    
    mqc_id = os.path.basename(out_samstat).split(".")[0]
    with open(samstat_filename,"r") as sf:
        samstat_lines = sf.readlines()
    start_remove = 0
    end_remove = 0
    #the string conditions below are unique in the samstat file
    for index,line in enumerate(samstat_lines): 
        if "Base quality distributions separated by mapping quality thresholds" in line:
            start_remove = index
        if "</footer>" in line:
            end_remove = index
    ###Need to rename canvas element variables since they have the same name across samstat reports and cause problems
    ###in the multiqc report. Number of canvas divs is a function of the sequence length: set to a high threshold (90)
    #canvas_vars = ["canvas1","canvas2","canvas3","canvas4"] 
    canvas_vars = ["canvas"+str(i) for i in range(1,90)] 
    canvas_replace = [x[:len(x)-1]+str(count)+x[len(x)-1] for x in canvas_vars]
    output_lines = samstat_lines[:start_remove] + samstat_lines[end_remove:]
    output_lines = "".join(output_lines) 
    for ci,canvas in enumerate(canvas_vars):
        output_lines = output_lines.replace(canvas,canvas_replace[ci])
    output_lines = output_lines.replace("animation: true","animation: false")
    with open(out_samstat,"w") as o:
        o.write("%s"%output_lines)
