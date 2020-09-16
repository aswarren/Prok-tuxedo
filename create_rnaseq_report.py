#!/usr/bin/env python

import sys
import os
import argparse
from bs4 import BeautifulSoup  

#Dictionaries with info used in the report
image_dict = {}

#Return whether condition_dict object contains all necessary fields for report
def validate_condition_dict(condition_dict):
    return True

#Return whether genome_list object contains all necessary fields for report
def validate_genome_list(genome_list):
    return True

#Return whether job_data object contains all necessary fields for report
def validate_job_data(job_data):
    return True

#function to parse out the images from the fastqc file
def get_images_from_fastqc(fastqc_file):
    global image_dict
    #TODO: modify to accept data instead
    image_list = []
    with open(fastqc_file,"r") as ff:
        fastqc_contents = ff.read()
        soup = BeautifulSoup(fastqc_contents,'lxml') 
    images = soup.find_all("img")
    #Put images into dictionary based on their 'alt' attribute name
    #Take the first "Per base quality graph" and skip the second
    alt_names_list = ["Per base quality graph","Per Sequence quality graph","Per base sequence content","Per sequence GC content graph","N content graph","Sequence length distribution","Duplication level graph","Adapter graph"]
    for i in images:
        if len(i.attrs) < 4:
            continue
        if i.attrs["alt"] in alt_names_list and i.attrs["alt"] not in image_dict:
            image_dict[i.attrs["alt"]] = i.attrs["src"]
    
def output_images():
    global image_dict
    header = "<html>\n"
    footer = "</html>"
    out_lines = []
    out_lines.append(header)
    for alt in image_dict:
        out_lines.append("<h2>")
        out_lines.append(alt)
        out_lines.append("</h2>\n")
        out_lines.append("<div>\n<img src=\"")
        out_lines.append(image_dict[alt])    
        out_lines.append("\"/>\n\n</div>\n")
    out_lines.append(footer)
    with open("Test_output.html","w") as o:
        o.write("".join(out_lines))

def setup(genome_list,condition_dict,job_data):
    global image_dict

def main(genome_list,condition_dict,job_data,report_file):
    #feature_count: htseq-count or stringtie, or cufflinks
    #recipe: Host or RNA-Rocket. Rockhopper depricated??
    #reference genome id
    get_images_from_fastqc()


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('-g','--genome_list',required=True)
    parser.add_argument('-c','--condition_dict',required=True)
    parser.add_argument('-j','--job_data',required=True)
    parser.add_argument('-o','--output_file',required=False,default="RNASeq_Report.txt")

    args = parser.parse_args()
    print(type(args))

    #check for correct arguments
    #TODO: check that error message prints object properly when validation is false
    if not validate_job_data(args.job_data):
        sys.stderr.write("job_data object does not contain all necessary fields:\n%s"%args.job_data)
        sys.exit(-1)
    if not validate_genome_list(args.genome_list):
        sys.stderr.write("genome_list object does not contain all necessary fields:\n%s"%args.genome_list)
        sys.exit(-1)
    if not validate_condition_dict(args.condition_dict):
        sys.stderr.write("condition_dict object does not contain all necessary fields:\n%s"%args.condition_dict)
        sys.exit(-1)
    main(args.genome_list,args.condition_dict,args.job_data,args.output_file)
