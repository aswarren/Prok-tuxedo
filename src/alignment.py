#!/usr/bin/env python

import os,sys,subprocess
import shutil
import tarfile
import json

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
            archive = tarfile.open(genome["hisat_index"])
            indices = [os.path.join(output_dir,os.path.basename(x)) for x in archive.getnames()]
            final_cleanup+=indices
            num_parts = archive.getnames()[0].count("/") #strip any subdirectory structure from the tar output
            archive.close()
            tar_cmd = ["tar","-xvf",genome["hisat_index"],"-C",output_dir]
            if num_parts > 0:
                tar_cmd += ["--strip-components",str(num_parts)]
            print (" ".join(tar_cmd))
            subprocess.check_call(tar_cmd)
            index_prefix = os.path.join(output_dir, os.path.basename(genome["hisat_index"]).replace(".ht2.tar","")) #somewhat fragile convention. tar prefix is underlying index prefix
            genome["index_prefix"] = index_prefix
            if not os.path.exists(index_prefix+".1.ht2"):
                print("hisant indices were not unpacked correctly: %s"%index_prefix)
                sys.exit()
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
        assign_strandedness_parameter(genome,condition_dict,parameters)
        ###
        if thread_count == 0:
            thread_count=2 #multiprocessing.cpu_count()
        cmd+=["-p",str(thread_count)]
        scount = 0 #used in modifying labels for the samstat reports
        for condition in condition_dict:
            rcount=0
            for r in condition_dict[condition]["replicates"]:
                scount+=1
                cur_cleanup=[]
                rcount+=1
                target_dir=r["target_dir"]
                fastqc_cmd=["fastqc","--outdir",target_dir]
                samstat_cmd=["samstat"]
                cur_cmd=list(cmd)
                if "read2" in r:
                    cur_cmd+=["-1",link_space(r["read1"]),"-2",link_space(r["read2"])]
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
                #add strandedness parameter
                if "strand_param" in r:
                    if cur_cmd[0] == "hisat2":
                        cur_cmd.insert["--rna-strandedness",r["strand_param"]]
                    else:
                        cur_cmd.insert(1,"--"+r["strand_param"].lower()) #strand_param is uppercase by default, but bowtie requires lowercase
                if not os.path.exists(r[genome["genome"]]["fastqc"]):
                    print (" ".join(fastqc_cmd))
                    pipeline_log.append(" ".join(fastqc_cmd))
                    #subprocess.check_call(fastqc_cmd)
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
                print " ".join(samstat_cmd)
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
                modify_samstat_for_multiqc(samstat_file,scount)
                for garbage in cur_cleanup:
                    subprocess.call(["rm", garbage])
        ###write out samstat json file
        os.chdir(genome["output"])
        with open("samstat.json","w") as so:
            so.write(json.dumps(samstat_dict))
        ###cleanup files
        for garbage in final_cleanup:
            subprocess.call(["rm", garbage])

#Only applies to paired-end data
def assign_strandedness_parameter(genome,condition_dict,parameters):
    #generate bed files for genome
    generate_bed_from_gff(genome)
    remove_list = []
    for condition in condition_dict:
        for r in condition_dict[condition]["replicates"]:
            #skip if not paired end reads
            if "read2" not in r:
                continue
            #skip if already determined somehow: probably only likely to happend during dualRNASeq, which is not implemented yet
            if "strand_param" in r:
                continue
            #if the number of reads is less than 1000, skip file:
            num_sample = "1000"
            r_count = 0
            with open(r["read1"],"r") as r1:
                for line in r1:
                    r_count = r_count + 1  
                    if r_count == int(num_sample):
                        break
            if r_count < int(num_sample):
                continue
            #sample fastq files using seqtk
            #TESTING: do I need to make sure the same reads have been sampled??
            sample1_file = os.path.join(r["target_dir"],".".join(os.path.basename(r["read1"]).split(".")[:-1]) + ".s1.fastq")
            sample2_file = os.path.join(r["target_dir"],".".join(os.path.basename(r["read1"]).split(".")[:-1]) + ".s2.fastq")
            remove_list.append(sample1_file)
            remove_list.append(sample2_file)
            #setting defulat sampling to 1000 reads
            sample1_cmd = " ".join(["seqtk","sample","-s","42",r["read1"],num_sample,">",sample1_file])
            sample2_cmd = " ".join(["seqtk","sample","-s","42",r["read2"],num_sample,">",sample2_file])
            ###TODO: replace shell=True for redirecting the output to a file
            print(sample1_cmd)
            if not os.path.exists(sample1_file):
                subprocess.check_call(sample1_cmd,shell=True)
            print(sample2_cmd)
            if not os.path.exists(sample2_file):
                subprocess.check_call(sample2_cmd,shell=True)
            ###run bowtie2 or hisat2 on sampled files
            sample_align_cmd = []
            sam_file = os.path.join(r["target_dir"],".".join(os.path.basename(r["read1"]).split(".")[:-1]) + ".sample.sam")
            remove_list.append(sam_file)
            ###TODO: replace hisat_index implementation with something more robust
            if len(genome["hisat_index"]) > 0:
                sample_align_cmd = ["hisat2","-x",genome["index_prefix"],"-1",sample1_file,"-2",sample2_file,"--mp","1,0","--pen_noncansplice","20","-S",sam_file,"-p",str(parameters.get("hisat2",{}).get("-p","1"))] 
            else: #bowtie
                sample_align_cmd = ["bowtie2","-x",genome["genome_link"],"-1",sample1_file,"-2",sample2_file,"-S",sam_file,"-p",str(parameters.get("bowtie2",{}).get("-p","1"))]
            if not os.path.exists(sam_file) and len(sample_align_cmd) > 0:
                print(sample_align_cmd)
                subprocess.check_call(sample_align_cmd)
            ###generate RSEQc strandedness file using infer_experiment.py
            infer_cmd = ["infer_experiment.py","-i",sam_file,"-r",genome["bed"],"-s",num_sample]
            infer_file = os.path.join(r["target_dir"],os.path.basename(r["read1"]).split(".")[0]+".infer") 
            print(" ".join(infer_cmd))
            with open(infer_file,"w") as o:
                subprocess.check_call(infer_cmd,stdout=o)
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
    ###
    os.chdir(genome["output"])
    with open("library_geometry.txt","w") as geom_handle:
        geom_handle.write("Sample\tCondition\tStrandedness\n")
        for condition in condition_dict:
            for r in condition_dict[condition]["replicates"]:
                if "strand_param" in r: 
                    geom_handle.write("%s\t%s\t%s\n"%(r["name"],condition,r["strand_param"]))
                else:
                    geom_handle.write("%s\t%s\t%s\n"%(r["name"],condition,"NA"))
                    

#genome = one of the genome dictionaries in the genome_list variable
def generate_bed_from_gff(genome):
    bed_file = genome["annotation_link"].replace(".gff",".bed")
    genome["bed"] = bed_file
    bed_cmd = "gff2bed --do-not-sort < " + genome["annotation_link"] + " > " + bed_file
    if not os.path.exists(bed_file):
        subprocess.check_call(bed_cmd,shell=True)

#Reads the output from samtools stat and grabs the average read length value
def get_average_read_length_per_file(stats_file):
    with open(stats_file,"r") as sf:
        for line in sf:
            if "average length:" in line: 
                line = line.strip().split()
                avg_length = line[-1]
                return avg_length

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
    ###in the multiqc report. Number of canvas divs is a function of the sequence length: set to a high threshold 
    #canvas_vars = ["canvas1","canvas2","canvas3","canvas4"] 
    canvas_vars = ["canvas"+str(i) for i in range(1,20)] 
    canvas_replace = [x[:len(x)-1]+str(count)+x[len(x)-1] for x in canvas_vars]
    output_lines = samstat_lines[:start_remove] + samstat_lines[end_remove:]
    output_lines = "".join(output_lines) 
    for ci,canvas in enumerate(canvas_vars):
        output_lines = output_lines.replace(canvas,canvas_replace[ci])
    output_lines = output_lines.replace("animation: true","animation: false")
    with open(out_samstat,"w") as o:
        o.write("%s"%output_lines)
