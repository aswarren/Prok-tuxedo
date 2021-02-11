from multiqc.modules.base_module import BaseMultiqcModule

class MultiqcModule(BaseMultiqcModule):
    def __init__(self):
        # Initialise the parent object
        super(MultiqcModule, self).__init__(name='Differential Gene Expression', anchor='dge_module',
        href="",
        info="")
        self.mod_data = dict()
        ###example of finding files specified in the search_patterns.yaml file
        for f in self.find_log_files('example'): 
            print(f)
            print(f['fn'])
        
        ###add section text content
        self.add_section(
            content = '<p>Insert report section</p>'
        )


