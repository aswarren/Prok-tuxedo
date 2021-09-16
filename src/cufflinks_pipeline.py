#!/usr/bin/env python3

import os,sys,subprocess
import cuffdiff_to_genematrix 

#TODO: smallRNA estimation in a future release
def run_cufflinks(genome_list, condition_dict, parameters, output_dir, pipeline_log):
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
        
        cmd=["cufflinks","--quiet","-G",genome["annotation"],"-b",genome_link,"-I","50"]
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
                    bam_tmp = os.path.join("/dev/shm","CUFFL")
                    shutil.copy(bam_file, bam_tmp)
                    print ("Copy succeeded to %s" % (bam_tmp))
                    bam_to_use = bam_tmp
                except IOError as err:
                    os.unlink(bam_tmp)
                    bam_to_use = None
                    bam_tmp = None
                    sys.stderr.write("Can't copy %s to %s: %s\n" % (bam_file, bam_tmp, err))

                if bam_to_use == None:
                    try:
                        bam_tmp = os.path.join(".","CUFFL")
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
                    print (" ".join(cur_cmd))
                    pipeline_log.append(" ".join(cur_cmd))
                    try:
                        sys.stderr.write("Invoke cufflinks: %s\n" % (cur_cmd))
                        subprocess.check_call(cur_cmd)
                    except Exception as e:
                        if bam_tmp != None:
                            sys.stderr.write("remove temp %s in exception handler\n" % (bam_tmp))
                            os.unlink(bam_tmp)
                            sys.stderr.write("Cufflinks error: %s\n" % (e))
                else:
                    sys.stderr.write(cuff_gtf+" cufflinks file already exists. skipping\n")
                if bam_tmp != None:
                    sys.stderr.write("remove temp %s\n" % (bam_tmp))
                    os.unlink(bam_tmp)

#Differential expression pipeline using the cufflinks protocol 
def run_cuffdiff(genome_list, condition_dict, parameters, output_dir, gene_matrix, contrasts, job_data, map_args, diffexp_json, pipeline_log):
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
            print (" ".join(merge_cmd))
            pipeline_log.append(" ".join(merge_cmd))
            subprocess.check_call(merge_cmd)
        else:
            sys.stderr.write(merge_file+" cuffmerge file already exists. skipping\n")

        #setup diff command
        cur_dir=genome["output"]
        thread_count= parameters.get("cuffdiff",{}).get("-p",0)
        if thread_count == 0:
            thread_count=2
        diff_cmd=["cuffdiff","--quiet",merge_file,"-p",str(thread_count),"-b",genome_link,"-L",",".join(condition_dict.keys())]
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
                quant_cmd=["cuffquant","--quiet",merge_file,r[genome_file]["bam"]]
                cur_dir=r[genome_file]["dir"]#directory for this replicate/genome
                os.chdir(cur_dir)
                quant_file=os.path.join(cur_dir,"abundances.cxb")
                quant_list.append(quant_file)
                if not os.path.exists(quant_file):
                    pipeline_log.append(quant_cmd)
                    subprocess.check_call(quant_cmd)
                else:
                    print (" ".join(quant_cmd))
                    sys.stderr.write(quant_file+" cuffquant file already exists. skipping\n")
            diff_cmd.append(",".join(quant_list))

        cur_dir=genome["output"]
        os.chdir(cur_dir)
        cds_tracking=os.path.join(cur_dir,"cds.fpkm_tracking")
        if not os.path.exists(cds_tracking):
            print (" ".join(diff_cmd))
            pipeline_log.append(" ".join(diff_cmd))
            subprocess.check_call(diff_cmd)
        else:
            sys.stderr.write(cds_tracking+" cuffdiff file already exists. skipping\n")
        #write gene_gmx file now for cuffdiff pipeline
        de_file=os.path.join(cur_dir,"gene_exp.diff")
        gmx_file=os.path.join(cur_dir,"gene_exp.gmx")
        if os.path.exists(de_file) and not os.path.exists(gmx_file):
            cuffdiff_to_genematrix.main([de_file],gmx_file)


