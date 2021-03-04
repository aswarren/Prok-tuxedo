<h2>Prok-Tuxedo</h2>
<p> 
An RNA-Seq analysis tool that incorporates a number of bioinformatics programs in a streamlined fashion and displays summary results in an html report. <br/>
This pipeline is used for running the <a href="https://patricbrc.org/app/Rnaseq" >PATRIC Transcriptomic Service </a>
 pipeline, which provides an interface for users to submit and customise their RNA-Seq experiment analyses along with a number of additional features.<br/>
Visit <a href="https://patricbrc.org/">PATRIC</a> for more information <br/>
</p>



usage: prok_tuxedo.py [-h] --jfile JFILE [--sstring SSTRING] -g G [--index]
                      [-p P] [-o O] -d D

An example run (small pair against itself): <br/>
python prok_tuxedo.py -o ./rnaseq_test/ -g ./test/baumanii_1505311/ -d .rnaseq_baumanii_1505311_diffexp --jfile ./test/baumanii_1505311/2cond_1comp_local.json --sstring {"data_api":"url_base_data_api"}
