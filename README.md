usage: prok_tuxedo.py [-h] -g G [-L L] [-p P] [-o O] readfiles [readfiles ...]

positional arguments:
  readfiles   whitespace sep list of read files. shoudld be in corresponding
              order as library list. ws separates libraries, a comma separates
              replicates, and a percent separates pairs.

optional arguments:
  -h, --help  show this help message and exit
  -g G        csv list of directories each containing a genome file *.fna and
              annotation *.gff
  -L L        csv list of library names for comparison
  -p P        JSON formatted parameter list for tuxedo suite keyed to program
  -o O        output directory. defaults to current directory.


An example run (small pair against itself):

python ./prok_tuxedo.py -L lib1,lib2 -g /home/anwarren/mcclelland/reference/14028s_test -p somefile /home/anwarren/mcclelland/test/READ1_SHORT.fastq%/home/anwarren/mcclelland/test/READ2_SHORT.fastq /home/anwarren/mcclelland/test/READ1_SHORT.fastq%/home/anwarren/mcclelland/test/READ2_SHORT.fastq
