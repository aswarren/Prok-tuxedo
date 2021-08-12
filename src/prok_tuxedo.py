#!/usr/bin/env python3

#import standard libraries
import os, sys, subprocess, glob, argparse
import multiprocessing
import tarfile, json
import requests,shutil
import math

#import scripts
import cuffdiff_to_genematrix
import alignment
import quantification
import cufflinks_pipeline
import prep_diffexp_files
import pathways
import multiqc_report
import multiqc_module_output as mmo
import unit_tests


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
            #target_dir=os.path.join(genome["output"], str(condition_index),"replicate"+str(rcount))
            target_dir=os.path.join(genome["output"], str(condition),"replicate"+str(rcount))
            target_dir=os.path.abspath(target_dir)
            r["target_dir"]=target_dir

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
            print (" ".join(convert_cmd))
            try:
               subprocess.check_call(convert_cmd)
            except(subprocess.CalledProcessError):
               sys.stderr.write("Running differential expression import failed.\n")
               #subprocess.call(["rm","-rf",experiment_path])
               return
            diffexp_obj_file=os.path.join(output_dir, os.path.basename(map_args.d.lstrip(".")))
            with open(diffexp_obj_file, 'w') as diffexp_job:
                diffexp_job.write(json.dumps(diffexp_json))


#Gets the lists of contrasts and runs DESeq2
#Runs DESeq2 once for each genome
def run_deseq2(genome_list,contrasts,job_data,dge_dict,output_dir):
    #Get list of contrasts to pass into deseq2 R script
    contrast_cmd = []
    for pair in contrasts:
        #remove any commas in names if present
        pair = [x.replace(",","_") for x in pair]
        contrast_cmd.append(",".join(pair))  

    #For each genome, run the deseq2 R script passing in the genome counts matrix, metadata file, and all contrasts
    #changes directory to the top genome directory for each genome
    #invoking run_deseq2.R: run_deseq2.R <counts_file.txt> <metadata_file.txt> <output_prefix> <feature_count> <contrast_1> <contrast_2> ... <contrast_n>
    for genome in genome_list:
        os.chdir(genome["output"])
        genome_prefix = genome["genome_id"]
        diffexp_list = [] #list of files to run differential expression on 
        diffexp_list.append((genome["gene_matrix"],"Genes",job_data.get("feature_count","htseq")))
        #transcript_matrix doesn't exist for bacterial pipeline
        if "transcript_matrix" in genome:
            diffexp_list.append((genome["transcript_matrix"],"Transcripts",job_data.get("feature_count","htseq")))
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
        genome["diff_exp_contrasts"] = []
        for diffexp_params in diffexp_list:
            #diffexp_cmd = ["run_deseq2.R",diffexp_params[0],metadata_file,diffexp_params[1],diffexp_params[2]]+contrast_cmd
            diffexp_cmd = ["run_deseq2",diffexp_params[0],metadata_file,diffexp_params[1],diffexp_params[2]]+contrast_cmd
            print("%s\n"%" ".join(diffexp_cmd))
            subprocess.check_call(diffexp_cmd)
            dge_dict[genome["genome"]]["volcano"] = wrap_svg_in_html("Volcano_Plots.svg",output_dir)
            for pair in contrasts:
                pair = [x.replace(",","_") for x in pair]
                pair = "_vs_".join(pair) 
                #diffexp_file = genome_prefix + "_" + pair + ".txt"
                diffexp_output = pair+"."+diffexp_params[2]+"."+diffexp_params[1]+".deseq2.tsv"
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
            #when creating gene_exp.gmx, ignore transcripts
            if "Transcript" in contrast_file:
                continue
            contrast_name = os.path.basename(contrast_file).replace(".tsv","")
            contrast_list.append(contrast_name)
            gene_count_dict[contrast_name] = {}
            with open(contrast_file,"r") as cf:
                next(cf)
                for line in cf:
                    gene,baseMean,log2FC,lfcSE,stat,pvalue,padj = line.strip().split("\t")
                    gene_set.add(gene.replace("gene-",""))
                    gene_count_dict[contrast_name][gene] = log2FC
        genome["gmx"] = os.path.join(genome["output"],"gene_exp.gmx")
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

#Strategy: 
#   - go through and keep two sorted lists: upregulated and downregulated
#   - pass top 50(?) into heatmap program 
#Returns: number of genes in the heatmap, skip heatmap if 0
def top_diffexp_genes(genome_list,dge_dict,feature_count): 
    #report number of genes that are significant in the heatmap and in the overall report
    pval_threshold = 0.05
    num_sig_genes = 0
    sig_genes = set()
    for genome in genome_list:
        os.chdir(genome["output"])
        contrast_file_list = genome["diff_exp_contrasts"]
        #List for upregulated genes
        up_genes_list = []
        up_log_list = []
        #List for downregulated genes
        down_genes_list = []
        down_log_list = []
        #grab all signification genes
        for contrast_file in contrast_file_list:
            if "Transcripts" in contrast_file and feature_count == "htseq":
                continue
            if "Genes" in contrast_file and feature_count == "stringtie":
                continue
            with open(contrast_file,"r") as cf:
                next(cf) #skip header
                for line in cf:
                    gene,baseMean,logFC,lfcSE,stat,pvalue,padj = line.strip().split("\t") 
                    #do not include genes with these label
                    if padj == "NA":
                        continue
                    if "STRG" in gene:
                        continue
                    if float(padj) < pval_threshold:
                        num_sig_genes = num_sig_genes + 1   
                        sig_genes.add(gene)
                    if float(logFC) < 0 and gene not in down_genes_list and gene not in up_genes_list:
                        down_genes_list.append(gene)
                        down_log_list.append(float(logFC))
                    elif float(logFC) > 0 and gene not in up_genes_list and gene not in down_genes_list:
                        up_genes_list.append(gene)
                        up_log_list.append(float(logFC))
        #zip and sort lists 
        up_sorted = sorted(tuple(zip(up_genes_list,up_log_list)),key=lambda x:x[1],reverse=True)
        down_sorted = sorted(tuple(zip(down_genes_list,down_log_list)),key=lambda x:x[1])
        #get subset of genes for heatmap
        num_up = 25 if len(up_sorted) > 25 else len(up_sorted)
        num_down = 25 if len(down_sorted) > 25 else len(down_sorted)
        heatmap_genes = list(up_sorted)[:num_up]+list(down_sorted)[:num_down]
        num_sig_heatmap_genes = 0
        for g in heatmap_genes:
            if g in sig_genes:
                num_sig_heatmap_genes = num_sig_heatmap_genes + 1
        if len(heatmap_genes) == 0:
            sys.stderr.write("No significant genes found in differential expression file: no heatmap output\n")
            continue 
        #write genes to file
        heatmap_genes_file = os.path.join(genome["output"],genome["genome_id"]+"_heatmap_genes.txt")    
        genome["heatmap_genes"] = heatmap_genes_file 
        with open(heatmap_genes_file,"w") as o:
            for gene in heatmap_genes:
                o.write("%s\n"%gene[0]) 
        dge_dict[genome["genome"]]["significant_genes"] = str(num_sig_genes)
        dge_dict[genome["genome"]]["heatmap_significant_genes"] = str(num_sig_heatmap_genes)

def generate_heatmaps(genome_list,job_data,dge_dict,output_dir):
    feature_count = "htseq" if job_data.get("feature_count","htseq") == "htseq" else "stringtie"
    for genome in genome_list:
        os.chdir(genome["output"])
        if not "heatmap_genes" in genome:
            continue
        #<heatmap_script.R> <gene_counts.txt> <metaata.txt> <heatmap_genes.txt> <output_prefix> <feature_count> <specialty_genes> <subsystem_genes>
        #Test_Htseq_four_cond/83333.13/Normalized_Top_50_Differentially_Expressed_Genes_mqc.svg
        #heatmap_cmd = ["generate_heatmaps.R",genome["genes_tpm_matrix"],genome["deseq_metadata"],genome["heatmap_genes"],genome["genome_id"],feature_count]
        heatmap_cmd = ["generate_heatmaps",genome["genes_tpm_matrix"],genome["deseq_metadata"],genome["heatmap_genes"],genome["genome_id"],feature_count]
        if "specialty_genes_map" in genome:
            heatmap_cmd += [genome["specialty_genes_map"]]
        else:
            heatmap_cmd += ["NONE"]
        if "superclass_map" in genome:
            heatmap_cmd += [genome["superclass_map"]]
        else:
            heatmap_cmd += ["NONE"]
        genome["heatmap_svg"] = os.path.join(genome["output"],"Normalized_Top_50_Differentially_Expressed_Genes.svg")
        print(" ".join(heatmap_cmd))
        subprocess.check_call(heatmap_cmd)
        dge_dict[genome["genome"]]["heatmap"] = wrap_svg_in_html(genome["heatmap_svg"],output_dir)

#Places <DOCTYPE> and <html></html> around svg code
#Returns html string, replacing svg with html from the svg_file parameter
def wrap_svg_in_html(svg_file,output_dir):
    html_file = svg_file.replace(".svg",".html")
    print("Wrapping %s in html: %s"%(svg_file,html_file))
    with open(svg_file,"r") as sf:
        svg_lines = sf.readlines()
    html_lines = ["<!DOCTYPE html>\n","<html>\n"]+svg_lines+["</html>"]
    with open(html_file,"w") as hf:
        hf.write("%s"%"".join(html_lines))
    move_svg(svg_file,output_dir)
    return move_html(html_file,output_dir)

#calls the list of svg files and moves them to one directory
#necessary for structure in the user workspace
def move_svg(svg,output_dir):
    svg_dir = os.path.join(output_dir,"report_images")
    if not os.path.exists(svg_dir):
        os.mkdir(svg_dir)
    out_svg = os.path.join(svg_dir,os.path.basename(svg))
    os.rename(svg,out_svg) 

def move_html(html,output_dir):
    html_dir = os.path.join(output_dir,"html_images")
    if not os.path.exists(html_dir):
        os.mkdir(html_dir)
    out_html = os.path.join(html_dir,os.path.basename(html))
    os.rename(html,out_html)
    return out_html

#removes html files 
#part of cleanup before sending files to user workspace
def remove_html_dir(output_dir):
    html_dir = os.path.join(output_dir,"html_images")
    if os.path.exists(html_dir):
        subprocess.check_call(["rm","-r",html_dir])

#Cleanup all remaining files not to be submitted to the users workspace
def cleanup_files(genome_list,output_dir):
    os.chdir(output_dir)
    cleanup_list = []
    #cleanup genome files
    for gff_file in glob.glob("*gff"):
        cleanup_list.append(os.path.abspath(gff_file))
    for bowtie_file in glob.glob("*bt2"):
        cleanup_list.append(os.path.abspath(bowtie_file))
    for fasta_file in glob.glob("*fna"):
        cleanup_list.append(os.path.abspath(fasta_file))
    for bed_file in glob.glob("*bed"):
        cleanup_list.append(os.path.abspath(bed_file))
    for hisat_file in glob.glob("*ht2"):
        cleanup_list.append(os.path.abspath(hisat_file))
    #cleanup html images
    if os.path.exists("html_images"):
        cleanup_list.append(os.path.abspath("html_images"))
    #cleanup files in each genome directory
    for genome in genome_list:
        os.chdir(genome["output"])
        #extra R plot pdf is generated: haven't figured out how to remove yet
        if os.path.exists("Rplots.pdf"):
            cleanup_list.append(os.path.abspath("Rplots.pdf"))
        for map_file in glob.glob("*mapping"):
            cleanup_list.append(os.path.abspath(map_file))
        for genes_file in glob.glob("*genes"):
            cleanup_list.append(os.path.abspath(genes_file))
        for json_file in glob.glob("*json"):
            cleanup_list.append(os.path.abspath(json_file))
        for intro_file in glob.glob("*pipeline"):
            cleanup_list.append(os.path.abspath(intro_file))
        for yaml_file in glob.glob("*yaml"): #multiqc config file
            cleanup_list.append(os.path.abspath(yaml_file))
    #remove all files in cleanup_list
    for rm_file in cleanup_list:
        subprocess.check_call(["rm","-r",rm_file])
        

###eventually if dual-RNASeq is enabled in the pipeline, need to adjust target_dir
#since it only references one genome directory, and a few function implementations depend on target_dir
# - assign_stradedness_parameter() in alignment.py
def setup(genome_list, condition_dict, parameters, output_dir, job_data):
    for genome in genome_list:
        genome_link=os.path.join(output_dir, os.path.basename(genome["genome"]))
        annotation_link = os.path.join(output_dir, os.path.basename(genome["annotation"]))
        genome["genome_link"]=genome_link
        genome["annotation_link"] = annotation_link
        genome["host"] = ("hisat_index" in genome and genome["hisat_index"])
        if not os.path.exists(genome_link):
            subprocess.check_call(["ln","-s",genome["genome"],genome_link])
        if not os.path.exists(annotation_link):
            subprocess.check_call(["ln","-s",genome["annotation"],annotation_link])
        make_directory_names(genome, condition_dict)
        rep_names = []
        for condition in condition_dict:
            rcount=0
            for r in condition_dict[condition]["replicates"]:
                ###assign name key in replicate dictionary
                if "bam" in r:
                    r_name = os.path.basename(r["bam"]).replace(".bam","")
                elif "srr_accession" in r:
                    r_name = r["srr_accession"]
                else:
                    r_name = os.path.basename(r["read1"]).split(".")[0]
                    if "read2" in r:
                        r_name = r_name + "_" + os.path.basename(r["read2"]).split(".")[0]
                if r_name in rep_names: #shouldn't happen but just in case
                    sys.stderr.write("Error, duplicate sample names: %s. Check json file or rename a sample\n"%(r_name))
                    sys.exit(-1)
                else:
                    rep_names.append(r_name)
                r["name"] = r_name
                ###
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
                                r["name"] = r["srr_accession"]+"_1_"+r["srr_accession"]+"_2"
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
    
    ###Setup dictionaries for the multiqc modules: use json dump to put these files at the top of each genome directory
    dge_dict = {}
    pathway_dict = {}
    ref_dict = {}
    ref_dict["recipe"] = job_data.get("recipe","RNA-Rocket")
    for genome in genome_list:
        dge_dict[genome["genome"]] = {}
        #dge_dict[genome["genome"]]["genome_id"] = os.path.basename(genome["genome"]).replace(".fna","")
        dge_dict[genome["genome"]]["genome_id"] = genome["genome_id"]
        pathway_dict[genome["genome"]] = {}
        #pathway_dict[genome["genome"]]["genome_id"] = os.path.basename(genome["genome"]).replace(".fna","")
        pathway_dict[genome["genome"]]["genome_id"] = genome["genome_id"]
    
    #pipeline_log holds the commands at each step of the pipeline and prints it to an output file Pipeline.txt
    pipeline_log = []
    #TRUE: runs cufflinks then cuffdiff if differential expression is turned on
    #FALSE: runs either htseq-count or stringtie
    run_cuffdiff_pipeline = job_data.get("feature_count","htseq") == "cuffdiff"
    alignment_value = alignment.run_alignment(genome_list, condition_dict, parameters, output_dir, job_data, pipeline_log)
    if alignment_value != 0:
        sys.stdout.write("Error in alignment step: look at stderr\n")
        cleanup_files(genome_list,output_dir)
        sys.exit(0)
    if run_cuffdiff_pipeline:
        cufflinks_pipeline.run_cufflinks(genome_list, condition_dict, parameters, output_dir)
    else:
        quant_val = quantification.run_featurecount(genome_list, condition_dict, parameters, output_dir, job_data, pipeline_log)
        if quant_val != 0:
            sys.stdout.write("Error in quantification step for %s: look at stderr\n"%job_data.get("feature_count","host"))
            #cleanup_files(genome_list,output_dir)
            sys.exit(0)
    prep_diffexp_files.create_metadata_file(genome_list,condition_dict,output_dir)
    if not run_cuffdiff_pipeline and job_data.get("feature_count","htseq") == "stringtie":
        prep_diffexp_files.write_gtf_list(genome_list,condition_dict) #function that writes the input for prepDE.py, which is a list of samples and paths to their gtf files. Do this for each genome
        prep_diffexp_files.prep_stringtie_diffexp(genome_list,condition_dict,job_data.get("recipe","RNA-Rocket") == "Host",pipeline_log)   
        prep_diffexp_files.create_tpm_matrix_stringtie(genome_list,condition_dict,host_flag=job_data.get("recipe","RNA-Rocket") == "Host")
    elif not run_cuffdiff_pipeline and job_data.get("feature_count","htseq") == "htseq": #htseq
        if job_data.get("recipe","RNA-Rocket") == "Host":
            prep_diffexp_files.create_counts_table_host(genome_list,condition_dict,job_data)
            prep_diffexp_files.create_tpm_matrix_htseq(genome_list,condition_dict,True,8) #true for host_flag
        else:
            prep_diffexp_files.create_counts_table(genome_list,condition_dict,job_data)
            prep_diffexp_files.create_tpm_matrix_htseq(genome_list,condition_dict,False,8) #false for host_flag
    if not run_cuffdiff_pipeline and not job_data.get("recipe","RNA-Rocket") == "Host":
        pathways.run_subsystem_analysis(genome_list,job_data,pathway_dict,output_dir)
        pathways.run_kegg_analysis(genome_list,job_data,pathway_dict,output_dir)
    
    dge_flag = False #used for multiqc_report
    if len(condition_dict.keys()) > 1 and not job_data.get("novel_features",False):
        dge_flag = True
        #If running cuffdiff pipeline, terminated after running expression import
        if run_cuffdiff_pipeline:
            cufflinks_pipeline.run_cuffdiff(genome_list, condition_dict, parameters, output_dir, gene_matrix, contrasts, job_data, map_args, diffexp_json)
            run_diff_exp_import(genome_list, condition_dict, parameters, output_dir, contrasts, job_data, map_args, diffexp_json)
            sys.exit(0)
        #volcano plots are generated in the same script that runs deseq2
        run_deseq2(genome_list,contrasts,job_data,dge_dict,output_dir)
        write_gmx_file(genome_list)
        try:
            if not job_data.get("recipe","RNA-Rocket") == "Host":
                run_diff_exp_import(genome_list, condition_dict, parameters, output_dir, contrasts, job_data, map_args, diffexp_json)
        except:
            print("Expression import failed")
        #get amr and specialty genes for labeling the heatmap
        if not job_data.get("recipe","RNA-Rocket") == "Host":
            pathways.setup_specialty_genes(genome_list)
        #generate heatmaps 
        top_diffexp_genes(genome_list,dge_dict,job_data.get("feature_count","htseq")) 
        generate_heatmaps(genome_list,job_data,dge_dict,output_dir)
    ###write the introduction to the multiqc report based on the recipe, number of samples, conditions, and contrasts
    for genome in genome_list:
        os.chdir(genome["output"])
        num_samples = len(job_data.get("single_end_libs",[])) + len(job_data.get("paired_end_libs",[])) + len(job_data.get("srr_libs",[]))
        mmo.write_introduction_pipeline(job_data.get("recipe","RNA-Rocket"),job_data.get("feature_count","htseq"),str(num_samples),str(len(job_data.get("experimental_conditions","0"))),str(len(job_data.get("contrasts","0"))))
    #write out any of the json files to the top level of the genomes
    for genome in genome_list:
        os.chdir(genome["output"])
        with open("differential_expression.json","w") as dge_handle:
            dge_handle.write(json.dumps(dge_dict[genome["genome"]]))
        with open("pathways.json","w") as pathway_handle:
            pathway_handle.write(json.dumps(pathway_dict[genome["genome"]]))
        with open("references.json","w") as ref_handle:
            ref_handle.write(json.dumps(ref_dict))
            
    multiqc_report.run_multiqc(genome_list,condition_dict,job_data,dge_flag=dge_flag)
    ###output pipeline txt file
    os.chdir(output_dir)
    with open("Pipeline.txt","w") as o:
        o.write("\n".join(pipeline_log))
    ###cleanup files not to be submitted to the user workspace
    cleanup_files(genome_list,output_dir)
    ###Run 'basic' test to check if all expected files exist
    #run top_n_genes test if parameter and file is included
    basic_result = unit_tests.run_unit_test("basic",genome_list,condition_dict,output_dir,contrasts,job_data)
    if "unit_json" in map_args: #unit test json is loaded in the pre-main block
        unit_tests.run_unit_test(map_args.unit_json,genome_list,condition_dict,output_dir,contrasts,job_data)
    if basic_result == 0:
        sys.exit(0) 
    else:
        sys.exit(2) #exit code specific to basic test 
        

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
    parser.add_argument('-u','--unit_test',
            help='follow a filename:\
                    \t- filename: file containing N genes is checked against the genes with the greatest abundance and differential expression (if flagged)',required=False)
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
    ###add bam files to library dictionary
    bam_count = 0
    for bam in job_data.get("bam_libs",[]):
        condition_index = int(bam.get("condition", count+1))-1 #copying condition position from read loop above
        condition_id = condition_list[condition_index] if got_conditions else "results"
        if "replicates" not in condition_dict[condition_id]:
            condition_dict[condition_id].setdefault("replicates",[])
        condition_dict[condition_id]["replicates"].append(bam)
        if "condition_index" not in condition_dict[condition_id]:
            condition_dict[condition_id].setdefault("condition_index",condition_index)
        bam_count = bam_count + 1
    genome_dirs=map_args.g.strip().split(',')
    genome_list=[]
    for g in genome_dirs:
        cur_genome={"genome":[],"annotation":[],"dir":g,"hisat_index":[]}
        cur_genome["genome_id"] = job_data.get("reference_genome_id",None)
        if not cur_genome["genome_id"]:
            sys.stderr("Cannot find reference genome ID in job submission json: aborting\n")
            sys.exit(-1)
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
        if map_args.index:
            if len(cur_genome["hisat_index"]) != 1:
                sys.stderr.write("Missing hisat index tar file for "+g+"\n")
                sys.exit(2)
            else:
                cur_genome["hisat_index"]=cur_genome["hisat_index"][0]
        genome_list.append(cur_genome)

    #check if unit test json
    if map_args.unit_test:
        with open(map_args.unit_test,"r") as unit_handle:
            map_args.unit_json = json.loads(unit_handle.read())  
    
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

