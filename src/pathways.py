#!/usr/bin/env python3

import sys,os,subprocess
import requests
import json
from prok_tuxedo import wrap_svg_in_html
from authenticate import authenticateByEnv

def run_subsystem_analysis(genome_list,job_data,pathway_dict,output_dir):
    for genome in genome_list:
        #retrieve subsystem information
        subsystem_dict = get_subsystem_mapping(genome)
        if not subsystem_dict:
            sys.stderr.write("Error in subsystem analysis for genome_id %s\n"%genome["genome"])
            return
        else:
            genome["subsystem_dict"] = subsystem_dict
    #write mapping file
    write_subsystem_mapping_files(genome_list) 
    #Run subsytem plotting R script
    #grid_violin_plots.R <subsystem_map.txt> <counts_file.txt|csv> <metadata.txt>  <subsystem_level> <feature_count>
    feature_count = "htseq" if job_data.get("feature_count","htseq") == "htseq" else "stringtie"
    subsystem_levels = ["Superclass"]
    subsystem_map = ["superclass_map"]
    for genome in genome_list:
        os.chdir(genome["output"])
        if not "superclass_map" in genome:
            continue
        for i,level in enumerate(subsystem_levels):
            #subsystem_plot_cmd = ["subsystem_violin_plots.R",genome[subsystem_map[i]],genome["gene_matrix"],genome["deseq_metadata"],level,feature_count]    
            #subsystem_plot_cmd = ["grid_violin_plots.R",genome[subsystem_map[i]],genome["genes_tpm_matrix"],genome["deseq_metadata"],level,feature_count]    
            subsystem_plot_cmd = ["grid_violin_plots",genome[subsystem_map[i]],genome["genes_tpm_matrix"],genome["deseq_metadata"],level,feature_count]    
            #output_grid_file = level + "_Pathway_Distribution_mqc.svg"
            output_grid_file = os.path.join(genome["output"],level + "_Pathway_Distribution.svg")
            print(" ".join(subsystem_plot_cmd))
            try:
                subprocess.check_call(subsystem_plot_cmd)
                pathway_dict[genome["genome"]]["subsystem_grid"] = wrap_svg_in_html(output_grid_file,output_dir)
            except Exception as e:
                sys.stderr.write("ERROR generating subsystem plots,skipping heatmap generation:\n{0}\n{1}\n".format(" ".join(subsystem_plot_cmd,e))
            #subsystem_violin_plot(subsystem_dict,genome["gene_matrix"],genome["deseq_metadata"],level,feature_count)
            
def run_kegg_analysis(genome_list,job_data,pathway_dict,output_dir):
    for genome in genome_list:
        kegg_dict = get_kegg_genes_mapping(genome) 
        if not kegg_dict:
            sys.stderr.write("Error in kegg analysis for genome_id %s\n"%genome["genome"])
            return
        else:
            genome["kegg_dict"] = kegg_dict
    #write mapping file
    write_kegg_mapping_files(genome_list)
    #Run kegg plotting R script
    feature_count = "htseq" if job_data.get("feature_count","htseq") == "htseq" else "stringtie"
    #Currently only one kegg level is processed, but the implementation can be extended to support more. Just add to the lists
    kegg_levels = ["Pathway_Class"]
    kegg_map = ["kegg_map"]
    for genome in genome_list:
        os.chdir(genome["output"])
        if not "kegg_map" in genome:
            continue
        for i,level in enumerate(kegg_levels):
            #kegg_plot_cmd = ["grid_violin_plots.R",genome[kegg_map[i]],genome["genes_tpm_matrix"],genome["deseq_metadata"],level,feature_count]
            kegg_plot_cmd = ["grid_violin_plots",genome[kegg_map[i]],genome["genes_tpm_matrix"],genome["deseq_metadata"],level,feature_count]
            output_kegg_grid_file = os.path.join(genome["output"],level + "_Pathway_Distribution.svg") 
            print(" ".join(kegg_plot_cmd))
            try:
                subprocess.check_call(kegg_plot_cmd)
                pathway_dict[genome["genome"]]["kegg_grid"] = wrap_svg_in_html(output_kegg_grid_file,output_dir)
            except Exception as e:
                sys.stderr.write("ERROR generating kegg plots, skipping heatmap generation:\n{0}\n{1}\n".format(" ".join(kegg_plot_cmd),e))

def get_subsystem_mapping(genome):
    genome_url = "https://patricbrc.org/api/subsystem/?eq(genome_id,"+genome["genome_id"]+")&limit(10000000)&http_accept=application/solr+json"
    req = requests.Request('GET',genome_url)
    authenticateByEnv(req)
    prepared = req.prepare()
    s = requests.Session()
    response = s.send(prepared)
    subsystem_dict = {}
    superclass_set = set()
    print("Retrieving subsystem ids for genome_id %s"%(genome["genome_id"]))
    print(genome_url)
    if not response.ok:
        sys.stderr.write("Failed to retrieve subsystem ids for genomd_id %s"%genome["genome_id"])
        return None
    ###Check if there are any entries for this genome id
    if len(response.json()['response']['docs']) == 0:
        sys.stderr.write("No subsystem gene ids found for genome_id %s\n"%genome["genome_id"])
        return None
    for entry in response.json()['response']['docs']: 
        subsystem_dict[entry['patric_id']] = {}
        subsystem_dict[entry['patric_id']]['Superclass'] = entry['superclass'] if len(entry['superclass']) > 0 else 'NONE'
        superclass_set.add(entry['superclass'])
    return subsystem_dict

def write_subsystem_mapping_files(genome_list):
    for genome in genome_list:
        if "subsystem_dict" not in genome:
            continue 
        subsystem_dict = genome["subsystem_dict"] 
        #put path into genome dictionary before writing the output file
        os.chdir(genome["output"])     
        superclass_map = genome["genome_id"]+".superclass_mapping"
        superclass_path = os.path.join(genome["output"],superclass_map)
        genome["superclass_map"] = superclass_path
        #write superclass file
        with open(superclass_map,"w") as sm:
            #write headers
            sm.write("Patric_ID\tSuperclass\n")
            for patric_id in subsystem_dict:
                sm.write("%s\t%s\n"%(patric_id,subsystem_dict[patric_id]["Superclass"])) 

def write_kegg_mapping_files(genome_list):
    for genome in genome_list:
        if "kegg_dict" not in genome:
            continue
        kegg_dict = genome["kegg_dict"]
        #put path into genome dictionary before writing the output file
        os.chdir(genome["output"])
        kegg_map = genome["genome_id"]+".kegg_mapping"
        kegg_path = os.path.join(genome["output"],kegg_map)
        genome["kegg_map"] = kegg_path
        with open(kegg_map,"w") as km: 
            #write headers
            km.write("Patric_ID\tPathway_Class\n")
            for patric_id in kegg_dict:
                km.write("%s\t%s\n"%(patric_id,kegg_dict[patric_id]["pathway_class"]))

def setup_specialty_genes(genome_list):
    for genome in genome_list:
        genome["specialty_genes_dict"] = get_specialty_genes_mapping(genome)
    write_specialty_genes_mapping_files(genome_list)

def write_specialty_genes_mapping_files(genome_list):
    for genome in genome_list:
        if "specialty_genes_dict" not in genome:
            continue 
        sp_dict = genome["specialty_genes_dict"]
        os.chdir(genome["output"])
        sp_map = genome["genome_id"]+".specialty_genes"
        sp_path = os.path.join(genome["output"],sp_map)
        genome["specialty_genes_map"] = sp_path
        with open(sp_map,"w") as o:
            o.write("Patric_ID\tSP_Field\n")
            for p_id in sp_dict:
                o.write("%s\t%s\n"%(p_id,sp_dict[p_id]['property']))

#https://patricbrc.org/api/sp_gene/?in(genome_id,(242231.10))&limit(8000)&select(property,patric_id)&http_accept=application/solr+json
def get_specialty_genes_mapping(genome):
    prefix_url = "https://patricbrc.org/api/sp_gene/?in(genome_id,("
    suffix_url = "))&limit(8000)&select(property,patric_id)&http_accept=application/solr+json"
    sp_gene_url = prefix_url + genome["genome_id"] + suffix_url
    print("Retrieving specialty gene ids for genome_id %s\n"%(genome["genome_id"]))
    print(sp_gene_url)
    req = requests.Request('GET',sp_gene_url)
    authenticateByEnv(req)
    prepared = req.prepare()
    s = requests.Session()
    response = s.send(prepared)
    sp_dict = {}
    if not response.ok:
        sys.stderr.write("Failed to retrieve specialty_gene ids for genome_id %s"%genome["genome_id"])
        return None
    for entry in response.json()['response']['docs']:
        sp_dict[entry['patric_id']] = {}
        sp_dict[entry['patric_id']]['property'] = entry['property']
    return sp_dict

def get_kegg_genes_mapping(genome):
    prefix_url = "https://patricbrc.org/api/pathway/?eq(genome_id,"
    suffix_url = ")&limit(10000000)&http_accept=application/solr+json"
    pathway_url = prefix_url + genome["genome_id"] + suffix_url
    print("Retrieving pathway mapping for genome_id %s\n"%(genome["genome_id"]))
    print(pathway_url) 
    req = requests.Request('GET',pathway_url)
    authenticateByEnv(req)
    prepared = req.prepare()
    s = requests.Session()
    response = s.send(prepared)
    kegg_dict = {}
    if not response.ok:
        sys.stderr.write("Failed to retrieve kegg gene ids for genome_id %s\n"%genome["genome_id"])
        return None 
    ###Check if there are any entries for this genome id
    if len(response.json()['response']['docs']) == 0:
        sys.stderr.write("No kegg gene ids found for genome_id %s\n"%genome["genome_id"])
        return None
    ###Two options for Kegg categories: pathway_name and pathway_class
    #quick distribution check for 208964.12 shows >100 pathway_name categories and around 15 for pathway_class
    #using pathway_class
    for entry in response.json()['response']['docs']:
        #some kegg pathway genes do not have corresponding PATRIC IDs: skipping those genes for now
        if "patric_id" not in entry:
            continue
        kegg_dict[entry["patric_id"]] = {}
        kegg_dict[entry["patric_id"]]["pathway_class"] = entry["pathway_class"]
    return kegg_dict

