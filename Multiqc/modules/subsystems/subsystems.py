from multiqc.modules.base_module import BaseMultiqcModule
import json

class MultiqcModule(BaseMultiqcModule):
    def __init__(self):
        # Initialise the parent object
        super(MultiqcModule, self).__init__(name='Subsystems', anchor='mymod',
        href="",
        info="")
        self.mod_data = dict()
        ###example of finding files specified in the search_patterns.yaml file
        for f in self.find_log_files('subsystems/vln_grid'): 
            path_file = f['fn']
    
        img_file = None
        with open(path_file,"r") as pf:
            path_json = json.load(pf)
            img_file = path_json["subsystem_grid"]
        #self.add_section(plot=None,content="/homes/clarkc/RNASeq_Pipeline/SRP220530_Ecoli_Single/Test_Htseq_four_cond_new_structure/83333.113/Superclass_Subsystem_Distribution_subsystems.html")  
        if img_file:
            img_html = None
            with open(img_file,"r") as i:
                img_html = i.read()
            self.add_section(content=img_html)
