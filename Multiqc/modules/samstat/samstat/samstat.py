from multiqc.modules.base_module import BaseMultiqcModule
import json,os

class MultiqcModule(BaseMultiqcModule):
    def __init__(self):
        # Initialise the parent object
        super(MultiqcModule, self).__init__(name='Samstat Reports', anchor='mymod',
        href="",
        info="")
        self.mod_data = dict()
        ###example of finding files specified in the search_patterns.yaml file
        for f in self.find_log_files('samstat/report'): 
            path_file = f['fn']
        
        path_json = None
        with open(path_file,"r") as pf:
            path_json = json.load(pf)

        if path_json:
            '''
            self.css = {
                'assets/css/multiqc_samstat.css' :
                os.path.join(os.path.dirname(__file__), 'assets', 'css', 'multiqc_samstat.css')
            }
            self.js = {
                'assets/js/multiqc_samstat.js' :
                os.path.join(os.path.dirname(__file__), 'assets', 'js', 'multiqc_samstat.js')
            }  
            '''
            img_list = []
            report_list = path_json["reports"] 
            for report in report_list:
                with open(report,"r") as r:
                    img_list.append(r.read())
            ###Add in collapsible section
            #collapse_prefix = "<button type=\"button\" onclick=\"click_collapse_samstat()\">Open Collapsible</button>\n<div class=\"samstat_block\">"
            #collapse_suffix = "</div>"
            #img_list = [collapse_prefix] + img_list + [collapse_suffix] 
            ###
            img_html = "\n".join(img_list)
            self.add_section(content=img_html)

        ###add section text content
        #self.add_section(
        #    content = '<p>Insert report section</p>'
        #)


