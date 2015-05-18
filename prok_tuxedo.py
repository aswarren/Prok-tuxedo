#!/usr/bin/env python

import os, sys
import argparse
import subprocess

def run_alignment(genome_list, library_dict, parameters, output_dir): 
    #modifies library_dict sub replicates to include 'bowtie' dict recording output files
    for genome in genome_list:
        genome_link=os.path.join(output_dir, os.path.basename(genome["genome"][0]))
        if not os.path.exists(genome_link):
            subprocess.check_call(["ln","-s",genome["genome"][0],genome_link])
        subprocess.check_call(["bowtie2-build", genome_link])
        cmd="bowtie2 -x "+genome_link
        target_dir=os.path.join(output_dir,os.path.basename(genome["dir"]))
        subprocess.call(["mkdir","-p","target_dir"])
        for library in library_dict:
            for r in library_dict[library]["replicates"]:
                cur_cmd=cmd
                if "read2" in r:
                    cur_cmd+=" -1 "+r["read1"]+" -2 "+r["read2"]
                    name1=os.path.splitext(os.path.basename(r["read1"]))
                    name2=os.path.splitext(os.path.basename(r["read2"]))
                    sam_file=os.path.join(target_dir,name1+"_"+name2+".sam")
                else:
                    cur_cmd+=" -U "+r["read1"]
                    name1=os.path.splitext(os.path.basename(r["read1"]))
                    sam_file=os.path.join(target_dir,name1+".sam")
                if os.path.exists(sam_file):
                    sys.stderr.write(sam_file+" alignments file already exists. skipping\n")
                    continue
                bam_file=sam_file[:-4]+".bam"
                cur_cmd+=" -S "+sam_file
                subprocess.check_call(cur_cmd) #call bowtie2
                subprocess.check_call("samtools view -Shu "+sam_file+" | \ samtools sort -o - - > "+bam_file)#convert to bam
                subprocess.call("rm "+sam_file)

def run_cufflinks(genome_list, library_dict, parameters):
    for genome in genome_list:
        subprocess.check_call(["bowtie2-build", "-l"])
        for library in library_dict:
            subprocess.check_call(["bowtie2", "-l"])

def run_cuffdiff(genome_list, library_dict, parameters):
    pass

def main(genome_list, library_dict, parameters_file, output_dir):
    #arguments:
    #list of genomes [{"genome":somefile,"annotation":somefile}]
    #dictionary of library dictionaries structured as {libraryname:{library:libraryname, replicates:[{read1:read1file, read2:read2file}]}}
    #parametrs_file is json parameters list keyed as bowtie, cufflinks, cuffdiff.
    parameters=json.load(open(args.sfile,'r'))
    run_alignment(genome_list, library_dict, parameters, output_dir)


if __name__ == "__main__":
    #modelling input parameters after rockhopper
    parser = argparse.ArgumentParser()
    parser.add_argument('-g', help='csv list of directories each containing a genome file *.fna and annotation *.gff', required=True)
    parser.add_argument('-L', help='csv list of library names for comparison', required=False)
    parser.add_argument('-p', help='JSON formatted parameter list for tuxedo suite keyed to program', required=True)
    parser.add_argument('readfiles', nargs='+', help="whitespace sep list of read files. shoudld be \
            in corresponding order as library list. ws separates libraries,\
            a comma separates replicates, and a percent separates pairs.")
    if len(sys.argv) ==1:
        parser.print_help()
        sys.exit(2)
    library_dict={}
    library_list=args.L.strip().split(',')
    #create library dict
    if not len(library_list): library_list.append("results")
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
    for g in genome_dirs:
        cur_genome={"genome":[],"annotation":[],"dir":g}
        for f in os.listdir("g"):
            if f.endswith(".fna") or f.endswith(".fa") or f.endswith(".fasta"):
                cur_genome["genome"].append(os.path.join(g,f))
            elif f.endswith(".gff"):
                cur_genome["annotation"].append(os.path.join(g,f))
        if len(cur_genome["genome"]) != 1:
            sys.stderr.write("Too many or too few fasta files present in "+g+"\n")
            sys.exit(2)
        if len(cur_genome["annotation"]) != 1:
            sys.stderr.write("Too many or too few gff files present in "+g+"\n")
            sys.exit(2)


        genome_list.append(cur_genome)
    
    main(genome_list,library_dict,args.p)


