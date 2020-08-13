#!/usr/bin/env python

import sys
import os
import argparse

#Return whether condition_dict object contains all necessary fields for report
def validate_condition_dict(condition_dict):
    return False

#Return whether genome_list object contains all necessary fields for report
def validate_genome_list(genome_list):
    return False

#Return whether job_data object contains all necessary fields for report
def validate_job_data(job_data):
    return False

def main(genome_list,condition_dict,job_data,report_file):
    #feature_count: htseq-count or stringtie, or cufflinks
    #recipe: Host or RNA-Rocket. Rockhopper depricated??
    #reference genome id
    print("nothing implemented yet")

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
