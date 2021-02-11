#!/usr/bin/env python

import sys,os,subprocess

SPACE = "    " #multiqc does not like the tab character, using a 4-space macro

# The easiest way to get the report in to the correct order
# is to write a config file and setup the sections
# based on each run
# - function for creating config file
# - base section names on the files: example, sections for each of the samstat reports then setup the order
def setup_multiqc_configs(genome_list,condition_dict):
    section_order = ["title","logo","sp","module_order","section_comments","remove_sections","custom_content"]
    config_dict = setup_shared_config_sections() 
    config_list = []
    for genome in genome_list:
        os.chdir(genome["output"])
        config_path = os.path.join(genome["output"],os.path.basename(genome["output"])+"_multiqc_config.yaml")
        config_list.append(config_path)
        ###write config file
        curr_config_list = []
        for section in section_order:
            if section not in config_dict:
                sys.stderr.write("%s not in section_comments section in the multiqc configuration file\n"%section)
                continue
            curr_config_list = curr_config_list + config_dict[section] 
        config_str = "\n".join(curr_config_list)
        with open(config_path,"w") as o:
            o.write(config_str)
    return config_list

###Adds <entry> to the <section> in the config_dict at position <index>
def add_to_config_section(config_dict,section,index,entry):
    if not section in config_dict:
        sys.stderr.write("%s is not a valid section in the config_dict\n"%section)
        return config_dict
    if index == 0:
        print("cannot add entries at index 0, overwrites section header")
        return config_dict
    if section == "custom_content":
        if index <= 1:
            print("for custom_content section, index must be greater than 1")
            return config_dict
        config_dict[section].insert(index,entry)
    elif section == "sp":
        print("haven't included support for section \"sp\" yet")
    else:
        config_dict[section].insert(index,entry)  
    return config_dict

def run_multiqc(genome_list,condition_dict):
    #config_path = "/homes/clarkc/RNASeq_Pipeline/Prok-tuxedo/Multiqc/multiqc_config.yaml"
    config_path_list = setup_multiqc_configs(genome_list,condition_dict) 
    debug_multiqc = False
    remove_data_dir = False
    force_overwrite = True
    for index,genome in enumerate(genome_list):
        os.chdir(genome["output"])
        report_name = os.path.basename(genome["output"])+"_report.html"
        multiqc_cmd = ["multiqc","--flat","-o",".","-n",report_name,"-t","simple",".","-c",config_path_list[index]]
        if remove_data_dir:
            multiqc_cmd += ["--no-data-dir"]
        if force_overwrite:
            multiqc_cmd += ["-f"]
        if debug_multiqc:
            multiqc_cmd += ["--lint"]
        print(" ".join(multiqc_cmd))
        subprocess.check_call(multiqc_cmd)

#returns a dictionary with strings for the sections contained in this function
def setup_shared_config_sections():
    config_dict = {}
    title_list = [
                "title: \"BVBRC Transcriptomic Service\"",
                "subtitle: \"RNASeq Analysis\"" 
                ]
    config_dict["title"] = title_list
    logo_list = [
                "custom_logo: \'/homes/clarkc/RNASeq_Pipeline/Prok-tuxedo/Multiqc/BV_BRC.png\'",
                "custom_logo_title: \'BV-BRC\'"
                ]
    config_dict["logo"] = logo_list
    sp_list = [
                "sp:",
                SPACE+"general_stats:",
                SPACE+SPACE+"fn: \'*.bam\'",
                SPACE+"hisat2:",
                SPACE+SPACE+"fn: \'*.hisat\'",
                SPACE+"bowtie2:",
                SPACE+SPACE+"fn: \'*.bowtie\'",
                SPACE+"samtools/stats:",
                SPACE+SPACE+"fn: \'*.samtools_stats\'",
                SPACE+"htseq:",
                SPACE+SPACE+"fn: \'*.counts\'"
                ]
    config_dict["sp"] = sp_list
    module_list = [
                "module_order:",
                SPACE+"- introduction",
                SPACE+"- fastqc",
                SPACE+"- bowtie2",
                SPACE+"- hisat2",
                SPACE+"- subsystems",
                SPACE+"- htseq",
                SPACE+"- samtools",
                SPACE+"- differential_expression",
                SPACE+"- custom_content",
                SPACE+"- samstat",
                SPACE+"- references"
                ]
    config_dict["module_order"] = module_list
    comments_list = [
        "section_comments:",
        SPACE+"general_stats: \""+get_general_stats_intro()+"\""
    ]
    config_dict["section_comments"] = comments_list
    remove_list = [
                "remove_sections:",
                SPACE+"- fastqc_status_checks",
                SPACE+"- fastqc_per_base_sequence_content",
                SPACE+"- samtools-stats"
                ]
    config_dict["remove_sections"] = remove_list
    custom_content_list = [
                "custom_content:",
                SPACE+"order:",
                #SPACE+SPACE+"- Superclass_Subsystem_Distribution",
                SPACE+SPACE+"- Volcano_Plots",
                SPACE+SPACE+"- Normalized_Top_50_Differentially_Expressed_Genes"
                ]
    config_dict["custom_content"] = custom_content_list
    return config_dict

#Returns a string used in the introduction of the general_stats section of the multiqc report
#strings in this section are read using a markup interpreter 
def get_general_stats_intro():
    gs_string = "<p style=font-size:60%;>Column Descriptors (hover over headers for detail):<br />"    
    gs_string = gs_string + "<b>%Assigned</b>: Percentage of sample reads HTSeq attributes to genomic feature abundance"
    gs_string = gs_string + " <b>M Assigned</b>: Number of reads (millions) HTSeq attributes to genomic feature abundance"
    gs_string = gs_string + " <b>Error rate</b>: mismatches (NM) / bases mapped (CIGAR), calculated by Samtools;"
    gs_string = gs_string + " <b>M Reads Mapped</b>: Number of reads (millions) mapped in the bam file, calculated by Samtools;"
    gs_string = gs_string + " <b>% Mapped</b>: Percentage of sample reads mapped in the bam file, calculated by Samtools;"
    gs_string = gs_string + " <b>M Total seqs</b>: Total sequences in the bam file (millions), calculated by Samtools"
    gs_string = gs_string + "</p>"
    return gs_string
