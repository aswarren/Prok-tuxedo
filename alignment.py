#!/usr/bin/env python

import os,sys,subprocess
import tarfile

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
def run_alignment(genome_list, condition_dict, parameters, output_dir, job_data): 
    #modifies condition_dict sub replicates to include 'bowtie' dict recording output files
    for genome in genome_list:
        genome_link = genome["genome_link"]
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
            if not os.path.exists(index_prefix+".1.ht2"):
                print("hisant indices were not unpacked correctly: %s"%index_prefix)
                sys.exit()
            cmd=["hisat2","--dta-cufflinks", "-x", index_prefix] 
            #cmd=["hisat2","--dta-cufflinks", "-x", index_prefix,"--no-spliced-alignment"] 
            thread_count= parameters.get("hisat2",{}).get("-p",0)
        else:
            #cmd=["hisat2","--dta-cufflinks", "-x", genome_link, "--no-spliced-alignment"] 
            try:
                subprocess.check_call(["bowtie2-build", genome_link, genome_link])
            except Exception as err:
                sys.stderr.write("bowtie build failed: %s %s\n"%(err, genome_link))
                os.exit(1)
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
                if not os.path.exists(r[genome["genome"]]["fastqc"]):
                    subprocess.check_call(fastqc_cmd)
                if os.path.exists(bam_file):
                    sys.stderr.write(bam_file+" alignments file already exists. skipping\n")
                else:
                    print cur_cmd
                    if job_data.get("recipe","RNA-Rocket") == "Host":
                        alignment_log = bam_file.replace("bam","hisat") 
                    else:
                        alignment_log = bam_file.replace("bam","bowtie") 
                    with open(alignment_log,"w") as al:
                        subprocess.check_call(cur_cmd,stdout=al,stderr=al) #call bowtie2 or hisat2
                samtools_threads = parameters.get("samtools",{}).get("-@",1)
                if not os.path.exists(bam_file):
                    subprocess.check_call("samtools view -Su "+sam_file+" | samtools sort -o - - -@ "+str(samtools_threads)+" > "+bam_file, shell=True)#convert to bam
                    subprocess.check_call("samtools index "+bam_file, shell=True)
                    #subprocess.check_call('samtools view -S -b %s > %s' % (sam_file, bam_file+".tmp"), shell=True)
                    #subprocess.check_call('samtools sort %s %s' % (bam_file+".tmp", bam_file), shell=True)
                print " ".join(samstat_cmd)
                stats_cmd = ["samtools","stats","-@",str(samtools_threads),bam_file]
                stats_outfile = bam_file.replace("bam","samtools_stats")
                if not os.path.exists(stats_outfile):
                    with open(stats_outfile,"w") as o:
                        subprocess.check_call(stats_cmd,stdout=o)
                r[genome["genome"]]["avg_read_length"] = get_average_read_length_per_file(stats_outfile)
                samstat_file = bam_file+".samstat.html"
                if not os.path.exists(samstat_file):
                    subprocess.check_call(samstat_cmd)
                for garbage in cur_cleanup:
                    subprocess.call(["rm", garbage])
        for garbage in final_cleanup:
            subprocess.call(["rm", garbage])

#Reads the output from samtools stat and grabs the average read length value
def get_average_read_length_per_file(stats_file):
    with open(stats_file,"r") as sf:
        for line in sf:
            if "average length:" in line: 
                line = line.strip().split()
                avg_length = line[-1]
                return avg_length

