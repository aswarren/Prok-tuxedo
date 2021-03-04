#!/usr/bin/env python

import sys,os,subprocess
import requests
import json
from prok_tuxedo import wrap_svg_in_html
from authenticate import authenticateByEnv

def run_subsystem_analysis(genome_list,job_data,pathway_json):
    for genome in genome_list:
        #retrieve subsystem information
        subsystem_dict = get_subsystem_mapping(genome)
        if not subsystem_dict:
            sys.stderr.write("Error in subsystem analysis for genome_id %s"%genome["genome"])
        else:
            genome["subsystem_dict"] = subsystem_dict
    #write mapping file
    write_subsystem_mapping_files(genome_list) 
    #Run subsytem plotting R script
    #grid_violin_plots.R <subsystem_map.txt> <counts_file.txt|csv> <metadata.txt>  <subsystem_level> <feature_count>
    feature_count = "htseq" if job_data.get("feature_count","htseq") == "htseq" else "stringtie"
    add_class = False
    if add_class:
        subsystem_levels = ["Superclass","Class"]
        subsystem_map = ["superclass_map","class_map"]
    else:
        subsystem_levels = ["Superclass"]
        subsystem_map = ["superclass_map"]
    for genome in genome_list:
        os.chdir(genome["output"])
        if not "superclass_map" in genome:
            continue
        for i,level in enumerate(subsystem_levels):
            #subsystem_plot_cmd = ["subsystem_violin_plots.R",genome[subsystem_map[i]],genome["gene_matrix"],genome["deseq_metadata"],level,feature_count]    
            subsystem_plot_cmd = ["grid_violin_plots.R",genome[subsystem_map[i]],genome["gene_matrix"],genome["deseq_metadata"],level,feature_count]    
            #output_grid_file = level + "_Pathway_Distribution_mqc.svg"
            output_grid_file = level + "_Pathway_Distribution.svg"
            print(" ".join(subsystem_plot_cmd))
            subprocess.check_call(subsystem_plot_cmd)
            pathway_json["subsystem_grid"] = wrap_svg_in_html(output_grid_file)
            #subsystem_violin_plot(subsystem_dict,genome["gene_matrix"],genome["deseq_metadata"],level,feature_count)
            
def run_kegg_analysis(genome_list,job_data,pathway_json):
    for genome in genome_list:
        kegg_dict = get_kegg_genes_mapping(genome) 
        if not kegg_dict:
            sys.stderr.write("Error in kegg analysis for genome_id %s"%genome["genome"])
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
            kegg_plot_cmd = ["grid_violin_plots.R",genome[kegg_map[i]],genome["gene_matrix"],genome["deseq_metadata"],level,feature_count]
            output_kegg_grid_file = level + "_Pathway_Distribution.svg" 
            print(" ".join(kegg_plot_cmd))
            subprocess.check_call(kegg_plot_cmd)
            pathway_json["kegg_grid"] = wrap_svg_in_html(output_kegg_grid_file)

def get_subsystem_mapping(genome):
    genome_url = "https://patricbrc.org/api/subsystem/?eq(genome_id,"+os.path.basename(genome["output"])+")&limit(10000000)&http_accept=application/solr+json"
    req = requests.Request('GET',genome_url)
    authenticateByEnv(req)
    prepared = req.prepare()
    s = requests.Session()
    response = s.send(prepared)
    subsystem_dict = {}
    superclass_set = set()
    class_set = set()
    print("Retrieving subsystem ids for genome_id %s"%(os.path.basename(genome["output"])))
    print(genome_url)
    if not response.ok:
        sys.stderr.write("Failed to retrieve subsystem ids for genomd_id %s"%os.path.basename(genome["output"]))
        return None
    ###Check if there are any entries for this genome id
    if len(response.json()['response']['docs']) == 0:
        sys.stderr.write("No subsystem gene ids found for genome_id %s\n"%os.path.basename(genome["output"]))
        return None
    for entry in response.json()['response']['docs']: 
        subsystem_dict[entry['patric_id']] = {}
        subsystem_dict[entry['patric_id']]['Superclass'] = entry['superclass'] if len(entry['superclass']) > 0 else 'NONE'
        superclass_set.add(entry['superclass'])
        subsystem_dict[entry['patric_id']]['Class'] = entry['class'] if len(entry['class']) > 0 else 'NONE'
        class_set.add(entry['class'])
    #subsystem_dict["superclass_set"] = superclass_set
    #subsystem_dict["class_set"] = class_set
    return subsystem_dict

def write_subsystem_mapping_files(genome_list):
    for genome in genome_list:
        if "subsystem_dict" not in genome:
            continue 
        subsystem_dict = genome["subsystem_dict"] 
        #put path into genome dictionary before writing the output file
        os.chdir(genome["output"])     
        superclass_map = os.path.basename(genome["output"])+".superclass_mapping"
        superclass_path = os.path.join(genome["output"],superclass_map)
        genome["superclass_map"] = superclass_path
        #TODO: remove class map
        class_map = os.path.basename(genome["output"])+".class_mapping"
        class_path = os.path.join(genome["output"],class_map)
        genome["class_map"] = class_path
        #write superclass file
        with open(superclass_map,"w") as sm, open(class_map,"w") as cm:
            #write headers
            sm.write("Patric_ID\tSuperclass\n")
            if "class_map" in genome:
                cm.write("Patric_ID\tClass\n")
            for patric_id in subsystem_dict:
                sm.write("%s\t%s\n"%(patric_id,subsystem_dict[patric_id]["Superclass"])) 
                if "class_map" in genome:
                    cm.write("%s\t%s\n"%(patric_id,subsystem_dict[patric_id]["Class"])) 

def write_kegg_mapping_files(genome_list):
    for genome in genome_list:
        if "kegg_dict" not in genome:
            continue
        kegg_dict = genome["kegg_dict"]
        #put path into genome dictionary before writing the output file
        os.chdir(genome["output"])
        kegg_map = os.path.basename(genome["output"])+".kegg_mapping"
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
        sp_map = os.path.basename(genome["output"])+".specialty_genes"
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
    sp_gene_url = prefix_url + os.path.basename(genome["output"]) + suffix_url
    print("Retrieving specialty gene ids for genome_id %s\n"%(os.path.basename(genome["output"])))
    print(sp_gene_url)
    req = requests.Request('GET',sp_gene_url)
    authenticateByEnv(req)
    prepared = req.prepare()
    s = requests.Session()
    response = s.send(prepared)
    sp_dict = {}
    if not response.ok:
        sys.stderr.write("Failed to retrieve specialty_gene ids for genome_id %s"%os.path.basename(genome["output"]))
        return None
    for entry in response.json()['response']['docs']:
        sp_dict[entry['patric_id']] = {}
        sp_dict[entry['patric_id']]['property'] = entry['property']
    return sp_dict

def get_kegg_genes_mapping(genome):
    prefix_url = "https://patricbrc.org/api/pathway/?eq(genome_id,"
    suffix_url = ")&limit(10000000)&http_accept=application/solr+json"
    pathway_url = prefix_url + os.path.basename(genome["output"]) + suffix_url
    print("Retrieving pathway mapping for genome_id %s\n"%(os.path.basename(genome["output"])))
    print(pathway_url) 
    req = requests.Request('GET',pathway_url)
    authenticateByEnv(req)
    prepared = req.prepare()
    s = requests.Session()
    response = s.send(prepared)
    kegg_dict = {}
    if not response.ok:
        sys.stderr.write("Failed to retrieve kegg gene ids for genome_id %s\n"%os.path.basename(genome["output"]))
        return None 
    ###Check if there are any entries for this genome id
    if len(response.json()['response']['docs']) == 0:
        sys.stderr.write("No kegg gene ids found for genome_id %s\n"%os.path.basename(genome["output"]))
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

#TODO: write this function
#TODO: url does not result in a list of amr genes, just an empty list
#https://patricbrc.org/api/genome_amr/?in(genome_id,(242231.10))&in(resistant_phenotype,(Resistant,Susceptible,Intermediate))&limit(1)&facet((pivot,(antibiotic,resistant_phenotype,genome_id)),(mincount,1),(limit,-1))&json(nl,map)
def get_amr_mapping(genome):
    prefix_url = "https://patricbrc.org/api/genome_amr/?in(genome_id,("
    suffix_url = "))&in(resistant_phenotype,(Resistant,Susceptible,Intermediate))&limit(1)&facet((pivot,(antibiotic,resistant_phenotype,genome_id)),(mincount,1),(limit,-1))&json(nl,map)"
    amr_url = prefix_url + os.path.basename(genome["output"]) + suffix_url
    print("Retrieving amr ids for genome_id %s\n"%(os.path.basename(genome["output"])))
    print(amr_url)
    req = requests.Request('GET',amr_url)
    authenticateByEnv(req)
    prepared = req.prepare()
    s = requests.Session()
    response = s.send(prepared)
    print(response)
    amr_dict = {}
    if not response.ok:
        sys.stderr.write("Failed to retrieve amr ids for genome_id %s"%os.path.basename(genome["output"]))
        return None
    for entry in response.json()['response']['docs']:
        print(entry)
        return None


