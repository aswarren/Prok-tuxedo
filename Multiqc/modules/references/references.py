from multiqc.modules.base_module import BaseMultiqcModule

class MultiqcModule(BaseMultiqcModule):
    def __init__(self):
        # Initialise the parent object
        super(MultiqcModule, self).__init__(name='References', anchor='ref_mod',
        href="",
        info="")
        self.mod_data = dict()
        ###example of finding files specified in the search_patterns.yaml file
        ref_json_file = None
        for f in self.find_log_files('references/json'): 
            ref_json_file = f['fn']
        
        ###add section text content
        self.add_section(
            content = get_references_html(ref_json_file)
        )


def get_references_html(ref_json):
    ref_dict = get_references_dictionary()
    content_str = "<table>\n"
    for idx,program in enumerate(ref_dict):
        content_str = content_str + "<tr><td style=\"padding:15px\">" + str(idx+1) + ".</td><td style=\"vertial-align:bottom\">" + ref_dict[program] + "</td></tr>"
    content_str = content_str + "</table>"
    return content_str
         

def get_references_dictionary():
    ref_dict = {
            "fastqc" : "Andrews, S. (2010). FastQC:  A Quality Control Tool for High Throughput Sequence Data [Online]. Available online at: http://www.bioinformatics.babraham.ac.uk/projects/fastqc/",
            "htseq" : "Simon Anders, Paul Theodor Pyl, Wolfgang Huber HTSeq — A Python framework to work with high-throughput sequencing data Bioinformatics (2014), in print, online at doi:10.1093/bioinformatics/btu638",
            "stringtie" : "Kovaka S, Zimin AV, Pertea GM, Razaghi R, Salzberg SL, Pertea M Transcriptome assembly from long-read RNA-seq alignments with StringTie2, Genome Biology 20, 278 (2019), doi:10.1186/s13059-019-1910-1",
            "samstat" : "Lassmann et al. (2010) \"SAMStat: monitoring biases in next generation sequencing data.\" Bioinformatics doi:10.1093/bioinformatics/btq614 [PMID: 21088025]",
            "samtools" : "Li H, Handsaker B, Wysoker A, Fennell T, Ruan J, Homer N, Marth G, Abecasis G, Durbin R; 1000 Genome Project Data Processing Subgroup. The Sequence Alignment/Map format and SAMtools. Bioinformatics. 2009 Aug 15;25(16):2078-9. doi: 10.1093/bioinformatics/btp352. Epub 2009 Jun 8. PMID: 19505943; PMCID: PMC2723002.",
            "multiqc" : "MultiQC: Summarize analysis results for multiple tools and samples in a single report Philip Ewels, Måns Magnusson, Sverker Lundin and Max Käller Bioinformatics (2016) doi: 10.1093/bioinformatics/btw354 PMID: 27312411",
            "deseq" : "Love, M.I., Huber, W., Anders, S. Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2 Genome Biology 15(12):550 (2014)",
            "enhanced_volcano" : "Kevin Blighe, Sharmila Rana and Myles Lewis (2019). EnhancedVolcano: Publication-ready volcano plots with enhanced colouring and labeling. R package version 1.4.0. https://github.com/kevinblighe/EnhancedVolcano",
            "complex_heatmap" : "Gu, Z. (2016) Complex heatmaps reveal patterns and correlations in multidimensional genomic data. Bioinformatics.",
            "hisat" : "Kim D, Langmead B and Salzberg SL. HISAT: a fast spliced aligner with low memory requirements. Nature Methods 2015",
            "bowtie" : "Langmead, B., Trapnell, C., Pop, M. et al. Ultrafast and memory-efficient alignment of short DNA sequences to the human genome. Genome Biol 10, R25 (2009). https://doi.org/10.1186/gb-2009-10-3-r25"
        }
    return ref_dict

