#!/homes/clarkc/miniconda3/bin/python3

import sys,os
sys.path.append('/homes/clarkc/miniconda3/lib/python3.7/site-packages/')
import seaborn 
import pandas as pd

def subsystem_violin_plot(subsystem_dict,counts_mtx_file,metadata_file,level,feature_count):
    counts_sep = "," if feature_count == "stringtie" else "\t"
    count_mtx = pd.read_csv(counts_mtx_file,sep=counts_sep)
    counts_mtx.head()
    #for patric_id in subsystem_dict: 
         

def run_subsystem_analysis(genome_list,job_data):
    for genome in genome_list:
        #retrieve subsystem information
        subsystem_dict = get_subsystem_mapping(genome)
        if not subsystem_dict:
            sys.stderr.write("Error in subsystem analysis for genome_id %s"%genome["genome"])
        else:
            genome["subsystem_dict"] = subsystem_dict
    ###Testing: Don't write mapping file, keep as dictionary
    #write mapping file
    #write_subsystem_mapping_files(genome_list) 
    #Run subsytem plotting R script
    #subsystem_violin_plots.R <subsystem_map.txt> <counts_file.txt|csv> <metadata.txt>  <subsystem_level> <feature_count>
    feature_count = "htseq" if job_data.get("feature_count","htseq") == "htseq" else "stringtie"
    subsystem_levels = ["Superclass","Class"]
    subsystem_map = ["superclass_map","class_map"]
    for genome in genome_list:
        os.chdir(genome["output"])
        for i,level in enumerate(subsystem_levels):
            #subsystem_plot_cmd = ["subsystem_violin_plots.R",genome[subsystem_map[i]],genome["gene_matrix"],genome["deseq_metadata"],level,feature_count]    
            #print(" ".join(subsystem_plot_cmd))
            #subprocess.check_call(subsystem_plot_cmd)
            subsystem_violin_plot(subsystem_dict,genome["gene_matrix"],genome["deseq_metadata"],level,feature_count)
            
def write_subsystem_mapping_files(genome_list):
    for genome in genome_list:
        if "subsystem_dict" not in genome:
            continue 
        subsystem_dict = genome["subsystem_dict"] 
        os.chdir(genome["output"])     
        superclass_map = os.path.basename(genome["output"])+".superclass_mapping"
        class_map = os.path.basename(genome["output"])+".class_mapping"
        superclass_path = os.path.join(genome["output"],superclass_map)
        class_path = os.path.join(genome["output"],class_map)
        genome["superclass_map"] = superclass_path
        genome["class_map"] = class_path
        #write superclass file
        with open(superclass_map,"w") as sm, open(class_map,"w") as cm:
            #write headers
            sm.write("Patric_ID\tSuperclass\n")
            cm.write("Patric_ID\tClass\n")
            for sub_id in subsystem_dict:
                #if sub_id == "superclass_set" or sub_id == "class_set":
                #    continue
                sm.write("%s\t%s\n"%(sub_id,subsystem_dict[sub_id]["superclass"])) 
                cm.write("%s\t%s\n"%(sub_id,subsystem_dict[sub_id]["class"])) 

#TODO: incorporate KB_Auth token check
def get_subsystem_mapping(genome):
    genome_url = "https://patricbrc.org/api/subsystem/?eq(genome_id,"+os.path.basename(genome["output"])+")&limit(10000000)&http_accept=application/solr+json"
    req = requests.Request('GET',genome_url)
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
    for entry in response.json()['response']['docs']: 
        subsystem_dict[entry['patric_id']] = {}
        subsystem_dict[entry['patric_id']]['Superclass'] = entry['superclass'] if len(entry['superclass']) > 0 else 'NONE'
        superclass_set.add(entry['superclass'])
        subsystem_dict[entry['patric_id']]['Class'] = entry['class'] if len(entry['class']) > 0 else 'NONE'
        class_set.add(entry['class'])
    #subsystem_dict["superclass_set"] = superclass_set
    #subsystem_dict["class_set"] = class_set
    return subsystem_dict

