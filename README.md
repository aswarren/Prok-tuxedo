<h2>Prok-Tuxedo</h2>
<p> 
An RNA-Seq analysis tool that incorporates a number of bioinformatics programs in a streamlined fashion and displays summary results in an html report. <br/>
This pipeline is used for running the PATRIC Transcriptomic Service pipeline, which allows users to submit their RNA-Seq processing requests and provides a number of additional features.<br/>
Visit PATRIC for more information: <br/>
<a href="https://patricbrc.org/app/Rnaseq" >PATRIC Transcriptomic Service </a> 
</p>

usage: prok_tuxedo.py [-h] --jfile JFILE [--sstring SSTRING] -g G [--index]
                      [-p P] [-o O] -d D

optional arguments:
  -h, --help         show this help message and exit
  --jfile JFILE      json file for job {"reference_genome_id": "1310806.3",
                     "experimental_conditions": ["c1_control",
                     "c2_treatment"], "output_file":
                     "rnaseq_baumanii_1505311", "recipe": "RNA-Rocket",
                     "output_path": "/anwarren@patricbrc.org/home/test",
                     "paired_end_libs": [{"read1":
                     "/anwarren@patricbrc.org/home/rnaseq_test/MHB_R1.fq.gz",
                     "read2":
                     "/anwarren@patricbrc.org/home/rnaseq_test/MHB_R2.fq.gz",
                     "condition": 1}, {"read1": "/anwarren@patricbrc.org/home/
                     rnaseq_test/MERO_75_R1.fq.gz", "read2": "/anwarren@patric
                     brc.org/home/rnaseq_test/MERO_75_R2.fq.gz", "condition":
                     2}], "contrasts": [[1, 2]]}
  --sstring SSTRING  json server string specifying api {"data_api":"url"}
  -g G               csv list of directories each containing a genome file
                     *.fna and annotation *.gff
  --index            flag for enabling using HISAT2 indices
  -p P               JSON formatted parameter list for tuxedo suite keyed to
                     program
  -o O               output directory. defaults to current directory.
  -d D               name of the folder for differential expression job folder
                     where files go

An example run (small pair against itself):
python prok_tuxedo.py -o ./rnaseq_test/ -g ./test/baumanii_1505311/ -d .rnaseq_baumanii_1505311_diffexp --jfile ./test/baumanii_1505311/2cond_1comp_local.json --sstring {"data_api":"url_base_data_api"}
