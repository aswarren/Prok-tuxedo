from multiqc.modules.base_module import BaseMultiqcModule
import json

class MultiqcModule(BaseMultiqcModule):
    def __init__(self):
        # Initialise the parent object
        super(MultiqcModule, self).__init__(name='Pathways', anchor='mymod',
        href="",
        info="")
        self.mod_data = dict()
        ###example of finding files specified in the search_patterns.yaml file
        path_file = None
        for f in self.find_log_files('subsystems/vln_grid'): 
            path_file = f['fn']
    
        if not path_file:
            return
        
        superclass_grid = None
        kegg_grid = None
        with open(path_file,"r") as pf:
            path_json = json.load(pf)
            superclass_grid = path_json["subsystem_grid"]
            kegg_grid = path_json["kegg_grid"]
            
        #self.add_section(plot=None,content="/homes/clarkc/RNASeq_Pipeline/SRP220530_Ecoli_Single/Test_Htseq_four_cond_new_structure/83333.113/Superclass_Subsystem_Distribution_subsystems.html")  
        img_list = []
        if superclass_grid:
            with open(superclass_grid,"r") as i:
                img_list.append("<h2>BV-BRC Superclass Distribution</h2>")
                img_list.append(i.read())
        if kegg_grid:
            with open(kegg_grid,"r") as i:
                img_list.append("<h2>Kegg Pathway Class Distribution</h2>")
                img_list.append(i.read())
        if len(img_list) > 0:
            img_html = "\n".join(img_list)
            self.add_section(content=img_html)
