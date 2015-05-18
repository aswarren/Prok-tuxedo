#!/usr/bin/env python

import os, sys
import argparse
import subprocess
import multiprocessing


def run_alignment(genome_list, library_dict, parameters, output_dir): 
    #modifies library_dict sub replicates to include 'bowtie' dict recording output files
    for genome in genome_list:
        genome_link=os.path.join(output_dir, os.path.basename(genome["genome"]))
        if not os.path.exists(genome_link):
            subprocess.check_call(["ln","-s",genome["genome"],genome_link])
        subprocess.check_call(["bowtie2-build", genome_link, genome_link])
        cmd=["bowtie2", "-x", genome_link]
        if genome["dir"].endswith('/') and not genome["dir"]=="./":
            genome["dir"]=genome["dir"][:-1]
        target_dir=os.path.join(output_dir,os.path.basename(genome["dir"]))
        subprocess.call(["mkdir","-p",target_dir])
        for library in library_dict:
            for r in library_dict[library]["replicates"]:
                cur_cmd=list(cmd)
                if "read2" in r:
                    cur_cmd+=["-1",r["read1"]," -2",r["read2"]]
                    name1=os.path.splitext(os.path.basename(r["read1"]))[0]
                    name2=os.path.splitext(os.path.basename(r["read2"]))[0]
                    sam_file=os.path.join(target_dir,name1+"_"+name2+".sam")
                else:
                    cur_cmd+=" -U "+r["read1"]
                    name1=os.path.splitext(os.path.basename(r["read1"]))
                    sam_file=os.path.join(target_dir,name1+".sam")
                bam_file=sam_file[:-4]+".bam"
                r[genome["genome"]]=bam_file
                if os.path.exists(sam_file) or os.path.exists(bam_file):
                    sys.stderr.write(sam_file+" alignments file already exists. skipping\n")
                    continue
                cur_cmd+=["-S",sam_file]
                thread_count=multiprocessing.cpu_count()
                if thread_count < 1: thread_count=1
                cur_cmd+=["-p",str(thread_count)]
                print cur_cmd
                subprocess.check_call(cur_cmd) #call bowtie2
                subprocess.check_call("samtools view -Shu "+sam_file+" | \ samtools sort -o - - > "+bam_file, shell=True)#convert to bam
                #subprocess.call(["rm", sam_file])

def run_cufflinks(genome_list, library_dict, parameters):
    for genome in genome_list:
        genome_link=os.path.join(output_dir, os.path.basename(genome["genome"]))
        if not os.path.exists(genome_link):
            subprocess.check_call(["ln","-s",genome["genome"],genome_link])
        cmd=["cufflinks","-g",genome["annotation"],"-b",genome_link,"-I","50"]
        for library in library_dict:
            for r in library_dict[library]["replicates"]:
                cur_dir=os.path.dirname(os.path.realpath(r[genome["genome"]]))
                os.chdir(cur_dir)
                cur_cmd=list(cmd)
                cur_cmd+=r[genome["genome"]]
                subprocess.check_call(cur_cmd)

def run_cuffdiff(genome_list, library_dict, parameters):
    pass

def main(genome_list, library_dict, parameters_file, output_dir):
    #arguments:
    #list of genomes [{"genome":somefile,"annotation":somefile}]
    #dictionary of library dictionaries structured as {libraryname:{library:libraryname, replicates:[{read1:read1file, read2:read2file}]}}
    #parametrs_file is json parameters list keyed as bowtie, cufflinks, cuffdiff.
    if os.path.exists(parameters_file):
    	parameters=json.load(open(parameters_file,'r'))
    else:
        parameters=[]
    run_alignment(genome_list, library_dict, parameters, output_dir)
    run_cufflinks(genome_list, library_dict, parameters)


if __name__ == "__main__":
    #modelling input parameters after rockhopper
    parser = argparse.ArgumentParser()
    parser.add_argument('-g', help='csv list of directories each containing a genome file *.fna and annotation *.gff', required=True)
    parser.add_argument('-L', help='csv list of library names for comparison', required=False)
    parser.add_argument('-p', help='JSON formatted parameter list for tuxedo suite keyed to program', required=True)
    parser.add_argument('-o', help='output directory. defaults to current directory.', required=False)
    parser.add_argument('readfiles', nargs='+', help="whitespace sep list of read files. shoudld be \
            in corresponding order as library list. ws separates libraries,\
            a comma separates replicates, and a percent separates pairs.")
    if len(sys.argv) ==1:
        parser.print_help()
        sys.exit(2)
    args = parser.parse_args()
    library_dict={}
    library_list=[]
    if args.L:
        library_list=args.L.strip().split(',')
    #create library dict
    if not len(library_list): library_list.append("results")
    if not args.o:
        output_dir="./"
    else:
        output_dir=args.o
    for lib in library_list:
        library_dict[lib]={"library":lib}
    count=0
    #add read/replicate structure to library dict
    for read in args.readfiles:
        replicates=read.split(',')
        rep_store=library_dict[library_list[count]]["replicates"]=[]
        for rep in replicates:
            pair=rep.split('%')
            pair_dict={"read1":pair[0]}
            if len(pair) == 2:
                pair_dict["read2"]=pair[1]
            rep_store.append(pair_dict)
        count+=1
    genome_dirs=args.g.strip().split(',')
    genome_list=[]
    for g in genome_dirs:
        cur_genome={"genome":[],"annotation":[],"dir":g}
        for f in os.listdir(g):
            if f.endswith(".fna") or f.endswith(".fa") or f.endswith(".fasta"):
                cur_genome["genome"].append(os.path.join(g,f))
            elif f.endswith(".gff"):
                cur_genome["annotation"].append(os.path.join(g,f))
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


        genome_list.append(cur_genome)
    
    main(genome_list,library_dict,args.p,output_dir)


