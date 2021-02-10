from multiqc.modules.base_module import BaseMultiqcModule

class MultiqcModule(BaseMultiqcModule):
    def __init__(self):
        # Initialise the parent object
        super(MultiqcModule, self).__init__(name='Introduction', anchor='introduction',
        href="",
        info="")
        self.mod_data = dict()
        
        intro_str = None        
        for f in self.find_log_files('introduction/pipeline'):
            print(f['fn'])
            with open(f['fn'],"r") as i:
                intro_str = i.read()

        ###add section text content
        if intro_str:
            intro_str = "<p>" + intro_str + "</p>"
            self.add_section(
                content = intro_str 
            )


