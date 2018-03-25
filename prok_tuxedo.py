#!/usr/bin/env python

import os, sys
import argparse
import subprocess
import multiprocessing
import cuffdiff_to_genematrix
import tarfile, json

#pretty simple: its for prokaryotes in that parameters will be attuned to give best performance and no tophat

def run_alignment(genome_list, condition_dict, parameters, output_dir): 
    #modifies condition_dict sub replicates to include 'bowtie' dict recording output files
    for genome in genome_list:
        genome_link=os.path.join(output_dir, os.path.basename(genome["genome"]))
        final_cleanup=[]
        if not os.path.exists(genome_link):
            subprocess.check_call(["ln","-s",genome["genome"],genome_link])
        if "hisat_index" in genome and genome["hisat_index"]:
            archive = tarfile.open(genome["hisat_index"])
            indices= [os.path.join(output_dir,os.path.basename(x)) for x in archive.getnames()]
            final_cleanup+=indices
            #archive.extractall(path=output_dir)
            archive.close()
            subprocess.check_call(["tar","-xvf", genome["hisat_index"], "-C", output_dir])
            index_prefix = os.path.join(output_dir, os.path.basename(genome["hisat_index"]).replace(".ht2.tar","")) #somewhat fragile convention. tar prefix is underlying index prefix
            cmd=["hisat2","--dta-cufflinks", "-x", index_prefix] 
        else:
            subprocess.check_call(["hisat2-build", genome_link, genome_link])
            cmd=["hisat2","--dta-cufflinks", "-x", genome_link, "--no-spliced-alignment"] 
            #cmd=["bowtie2", "-x", genome_link]
        thread_count=multiprocessing.cpu_count()
        cmd+=["-p",str(thread_count)]
        if genome["dir"].endswith('/'):
            genome["dir"]=genome["dir"][:-1]
        genome["dir"]=os.path.abspath(genome["dir"])
        genome["output"]=os.path.join(output_dir,os.path.basename(genome["dir"]))
        for library in condition_dict:
            rcount=0
            for r in condition_dict[library]["replicates"]:
                cur_cleanup=[]
                rcount+=1
                target_dir=os.path.join(genome["output"],library,"replicate"+str(rcount))
                target_dir=os.path.abspath(target_dir)
                subprocess.call(["mkdir","-p",target_dir])
                cur_cmd=list(cmd)
                if "read2" in r:
                    cur_cmd+=["-1",r["read1"]," -2",r["read2"]]
                    name1=os.path.splitext(os.path.basename(r["read1"]))[0]
                    name2=os.path.splitext(os.path.basename(r["read2"]))[0]
                    sam_file=os.path.join(target_dir,name1+"_"+name2+".sam")
                else:
                    cur_cmd+=[" -U",r["read1"]]
                    name1=os.path.splitext(os.path.basename(r["read1"]))[0]
                    sam_file=os.path.join(target_dir,name1+".sam")
                cur_cleanup.append(sam_file)
                bam_file=sam_file[:-4]+".bam"
                r[genome["genome"]]={}
                r[genome["genome"]]["bam"]=bam_file
                cur_cmd+=["-S",sam_file]
                if os.path.exists(bam_file):
                    sys.stderr.write(bam_file+" alignments file already exists. skipping\n")
                else:
                    print cur_cmd
                    subprocess.check_call(cur_cmd) #call bowtie2
                if not os.path.exists(bam_file):
                    subprocess.check_call("samtools view -Su "+sam_file+" | samtools sort -o - - > "+bam_file, shell=True)#convert to bam
                    subprocess.check_call("samtools index "+bam_file, shell=True)
                    #subprocess.check_call('samtools view -S -b %s > %s' % (sam_file, bam_file+".tmp"), shell=True)
                    #subprocess.check_call('samtools sort %s %s' % (bam_file+".tmp", bam_file), shell=True)
                for garbage in cur_cleanup:
                    subprocess.call(["rm", garbage])
        for garbage in final_cleanup:
            subprocess.call(["rm", garbage])

def run_cufflinks(genome_list, condition_dict, parameters, output_dir):
    for genome in genome_list:
        genome_file=genome["genome"]
        genome_link=os.path.join(output_dir, os.path.basename(genome["genome"]))
        if not os.path.exists(genome_link):
            subprocess.check_call(["ln","-s",genome["genome"],genome_link])
        cmd=["cufflinks","-q","-g",genome["annotation"],"-b",genome_link,"-I","50"]
        thread_count=multiprocessing.cpu_count()
        cmd+=["-p",str(thread_count)]
        for library in condition_dict:
            for r in condition_dict[library]["replicates"]:
                cur_dir=os.path.dirname(os.path.realpath(r[genome_file]["bam"]))
                os.chdir(cur_dir)
                cur_cmd=list(cmd)
                r[genome_file]["dir"]=cur_dir
                cur_cmd+=[r[genome_file]["bam"]]#each replicate has the bam file
                cuff_gtf=os.path.join(cur_dir,"transcripts.gtf")
                if not os.path.exists(cuff_gtf):
                    print " ".join(cur_cmd)
                    subprocess.check_call(cur_cmd)
                else:
                    sys.stderr.write(cuff_gtf+" cufflinks file already exists. skipping\n")

def run_diffexp(genome_list, condition_dict, parameters, output_dir, gene_matrix, contrasts, job_data):
    #run cuffquant on every replicate, cuffmerge on all resulting gtf, and cuffdiff on the results. all per genome.
    for genome in genome_list:
        genome_file=genome["genome"]
        genome_link=os.path.join(output_dir, os.path.basename(genome["genome"]))
        if not os.path.exists(genome_link):
            subprocess.check_call(["ln","-s",genome["genome"],genome_link])
        merge_cmd=["cuffmerge","-g",genome["annotation"]]
        merge_manifest=os.path.join(genome["output"],"gtf_manifest.txt")
        merge_folder=os.path.join(genome["output"],"merged_annotation")
        merge_file=os.path.join(merge_folder,"merged.gtf")
        thread_count=multiprocessing.cpu_count()
        merge_cmd+=["-p",str(thread_count),"-o", merge_folder]
        for library in condition_dict:
            for r in condition_dict[library]["replicates"]:
                with open(merge_manifest, "a") as manifest: manifest.write("\n"+os.path.join(r[genome_file]["dir"],"transcripts.gtf"))
        merge_cmd+=[merge_manifest]
        diff_cmd=["cuffdiff",merge_file,"-p",str(thread_count),"-b",genome_link,"-L",",".join(condition_dict.keys())]
        if not os.path.exists(merge_file):
            print " ".join(merge_cmd)
            subprocess.check_call(merge_cmd)
        else:
            sys.stderr.write(merge_file+" cuffmerge file already exists. skipping\n")
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
        contrasts_file = os.path.join(cur_dir, "contrasts.txt")
        with open(contrasts_file,'w') as contrasts_handle:
            constrasts_handle.write("condition_A\tcondition_B\n")
            for c in contrasts:
                contrasts_handle.write(str(c[0])+"\t"+str(c[1])+"\n")
        diff_cmd+=["--contrast-file",contrasts_file]
        if not os.path.exists(cds_tracking):
            print " ".join(diff_cmd)
            subprocess.check_call(diff_cmd)
        else:
            sys.stderr.write(cds_tracking+" cuffdiff file already exists. skipping\n")
        if gene_matrix:
            de_file=os.path.join(cur_dir,"gene_exp.diff")
            gmx_file=os.path.join(cur_dir,"gene_exp.gmx")
            cuffdiff_to_genematrix.main([de_file],gmx_file)
            trasform_script =os.path.join(os.path.realpath(__file__),"p3diffexp","expression_transform.py")
            if os.path.exists(transform_script) and os.path.exists(gmx_file):
                transform_params = {"output_path":cur_dir, "xfile":gmx_file, "xformat":"tsv",\
                        "xsetup":"gene_matrix", "source_id_type":"feature_id",\
                        "data_type":"Transcriptomics", "title":"RNA-Seq", "description":"RNA-Seq",\
                        "organism":job_data.genome_id}
                params_file=os.path.join(cur_dir, "diff_exp_params.json")
                with open(params_file, 'w') as params_handle:
                    params_handle.write(json.dumps(transform_params))
                convert_cmd=[transform_script, "--ufile", params_file, "sstring", job_data.sstring]
                print " ".join(convert_cmd)
            #convert_cmd+=]
            #subprocess.check_call(convert_cmd)
            

def main(genome_list, condition_dict, parameters_file, output_dir, gene_matrix=False, contrasts=[], job_data=None):
    #arguments:
    #list of genomes [{"genome":somefile,"annotation":somefile}]
    #dictionary of library dictionaries structured as {libraryname:{library:libraryname, replicates:[{read1:read1file, read2:read2file}]}}
    #parametrs_file is json parameters list keyed as bowtie, cufflinks, cuffdiff.
    output_dir=os.path.abspath(output_dir)
    subprocess.call(["mkdir","-p",output_dir])
    if parameters_file and os.path.exists(parameters_file):
    	parameters=json.load(open(parameters_file,'r'))
    else:
        parameters=[]
    run_alignment(genome_list, condition_dict, parameters, output_dir)
    run_cufflinks(genome_list, condition_dict, parameters, output_dir)
    if len(condition_dict.keys()) > 1:
        run_diffexp(genome_list, condition_dict, parameters, output_dir, gene_matrix, contrasts, job_data)


if __name__ == "__main__":
    #modelling input parameters after rockhopper
    parser = argparse.ArgumentParser()
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
    parser.add_argument('-p', help='JSON formatted parameter list for tuxedo suite keyed to program', required=False)
    parser.add_argument('-o', help='output directory. defaults to current directory.', required=False)
    #parser.add_argument('-x', action="store_true", help='run the gene matrix conversion and create a patric expression object', required=False)
    #parser.add_argument('readfiles', nargs='+', help="whitespace sep list of read files. shoudld be \
    #        in corresponding order as library list. ws separates libraries,\
    #        a comma separates replicates, and a percent separates pairs.")
    if len(sys.argv) ==1:
        parser.print_help()
        sys.exit(2)
    map_args = parser.parse_args()
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
    if not len(condition_list): condition_list.append("results")
    gene_matrix=True
    if not map_args.o:
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
    for read in job_data.get("paired_end_libs",[])+job_data.get("single_end_libs",[]):
        if "read" in read:
            read["read1"] = read.pop("read")
        condition_index = int(read.get("condition", count+1))-1 #if no condition return position so everything is diff condition
        rep_store=condition_dict[condition_list[condition_index]].setdefault("replicates",[]).append(read)
        count+=1
    genome_dirs=args.g.strip().split(',')
    genome_list=[]
    for g in genome_dirs:
        cur_genome={"genome":[],"annotation":[],"dir":g,"hisat_index":[]}
        for f in os.listdir(g):
            if f.endswith(".fna") or f.endswith(".fa") or f.endswith(".fasta"):
                cur_genome["genome"].append(os.path.abspath(os.path.join(g,f)))
            elif f.endswith(".gff"):
                cur_genome["annotation"].append(os.path.abspath(os.path.join(g,f)))
            elif f.endswith(".ht2.tar"):
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
        if args.index:
            if len(cur_genome["hisat_index"]) != 1:
                sys.stderr.write("Missing hisat index tar file for "+g+"\n")
                sys.exit(2)
            else:
                cur_genome["hisat_index"]=cur_genome["hisat_index"][0]


        genome_list.append(cur_genome)
    
    main(genome_list,condition_dict,args.p,output_dir,gene_matrix,contrasts, map_args)


