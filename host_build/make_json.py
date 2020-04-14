#!/usr/bin/env python

import os, sys
import argparse
import json



def main(map_args):
    #arguments:
    #file tsv of ncbi host genome information
    #ftp prefix for patric
    ftp_prefix=map_args.get('ftp_prefix')
    result={'genomes':[]}
    genomes=result.get('genomes')
    warned=False
    with open(map_args.get('host_table')) as fh:
        header=[]
        for line in fh:
            if line.startswith('#'):
                if "assembly_accession" in line:
                    header=[i.strip() for i in line.lstrip('#').strip().split('\t') if i.strip()]
                continue
            parts=line.strip().split("\t")
            diff = -1
            if len(parts)< len(header):
                if not warned:
                    sys.stderr.write("warning header more fields than table\n")
                    warned=True
                diff = len(parts)-len(header)-1
            x={v:parts[i] for i,v in enumerate(header[:diff])}
            file_prefix=x.get('assembly_accession')+"_"+x.get('asm_name')
            x.setdefault("patric_ftp","/".join([ftp_prefix,x.get('taxid'),x.get('asm_name'),file_prefix]))
            genomes.append(x)
    with sys.stdout as fh:
        fh.write(json.dumps(result))

if __name__ == "__main__":
    #modelling input parameters after rockhopper
    parser = argparse.ArgumentParser()
    #if you want to support multiple genomes for alignment you should make this json payload an nargs+ parameter
    parser.add_argument('--host_table', help='host table tsv as built by build_host.sh', required=True)
    parser.add_argument('--ftp_prefix', help='patric url of host genome locations', required=False, default='ftp://ftp.patricbrc.org/host_genomes')
    if len(sys.argv) ==1:
        parser.print_help()
        sys.exit(2)
    pargs = parser.parse_args()
    main(vars(pargs))


