# d-chimer_v1

Disentangle chimeric (d-chimer) sequences in de novo assembled viromes/metagenomes is a BLAST based pipeline conceived for taxonomic assignments by taking into account that :
- Contigs can be only partly covered in a single BLAST search 
- Contigs can be covered by different genes from different organisms


## 1. **Usage**
 
Once, d-chimer installed and the parameters configured (through the `d-chimer_config.yaml`, file), it can be used as follow :

    python3 /path/to/d-chimer/dchimer.py  -p blastn -L True -f query_file.fasta *
    
\* We recommand users to execute d-chimer in the directory where they input data are. 
    

By default, d-chimer runs in a recursive manner. The number of cycles can be limited by using the -m option. For example,  to limit d-chimer to two cycles :

      python3 /path/to/d-chimer/dchimer.py  -p blastn -L True -f query_file.fasta -m 2


For full view of d-chimer options, type :


      python3 /path/to/d-chimer/dchimer.py -h

      run a blast program and filter its outputs

      optional arguments:
        -h, --help            show this help message and exit
        -p blastprogram, --program blastprogram
                              blastProgram (required) : blastn or blastx
        -L local, --local local
                              use local machine BLAST programs (required) :
                              True/False (default : True)
        -f fastafile, --fasta fastafile
                              input fasta_file (required)
        -m max_loops, --max_recursive_loops max_loops
                              maximum number of times to to process uncovered zones
                              fasta      

      usage : dchimer [-h] -p blastProgram -L True -f fastafile [-m] max_loops

#### Usage exemple with the provided input file "testSequences.fasta"
- run d-chimer with the input file :
      
           python3 /path/to/d-chimer/dchimer.py -p blastn -L True -f testSequences.fasta

for each cycle, will be written :
- the blastn output : *testSequences.1.bn.csv* (cycle 1 blastn output).
It is a 17 column tabultated BLASTn output file : the 12 firts are default BLAST outputs + five columns standing respectively for *taxid*, *kingdom*, *query length*, *query sequence alignment* and *subject sequence alignment* 
- the uncovered_zones : *testSequences.1.bn.fas*
- the blastn negative sequences : *testSequences.1.bn_neg.fasta*
- a file to ensure that there no mismatches during taxonomy joining to filtered outputs (must be empty, if correct): *testSequences.1.bn.filtered.notfound_taxo.tsv*.
- A unique folder containing named : *testSequences_bn_out/*
  - the filtered outputs : *testSequences.1.bn.filtered.sorted.tsv*
  
    It is a ten column file which contains contains : *contig identifiant*, *taxonomic identifiant*, *subject accession number*, *identity percent*, *e-value*, *alignment length*, *query alignment start coordinate*, *query alignment end coordinate*, *subject alignment start coordinate*, *subject alignment end coordinate*

        NODE_r2_2_length_2963_cov_4.73413       2202954 MH618085        70.058  4.68e-46        521     954     1461    832     338
 
 - taxomony associated filered outputs : *testSequences.1.bn.filtered.taxo*
   
   It is an eleven column file which contains contains : *taxonomic identifiant*, *contig identifiant*, *subject accession number*, *identity percent*, *e-value*, *alignment length*, *query alignment start coordinate*, *query alignment end coordinate*, *subject alignment start coordinate*, *subject alignment end coordinate* and *full taxomic path*
 
         2202954 NODE_r2_2_length_2963_cov_4.73413       MH552562        65.545  6.95e-25        624     1163    1765    2514    1919    Viruses;unclassifiedviruses;            Circulargeneticelementsp.
 
 
 - query/subjects alignments : *testSequences.1.bn.filtered.aln*
 Its a three lines item with a description followed by two lines each representing the BLAST aligned sequence for contig and subject.

All these files and the folder are produced where d-chimer were executed. The same set of file are produced with blastx... everywhere "bn" is replaced by "bx".

## 2. **Installation**
Users can execute the `INSTALL.sh` bash script once this repo is cloned, to get the d-chimer and python libraries and custom data.
### 2.1 **d-chimer code and python libraries:** 
- Clone or download the d-chimer repository.   

             git clone https://github.com/stirera/d-chimer_v1


d-chimer depends on several python3 libraries and ncbi-BLAST and databases.
- Python libraries
  - Biopython and PyYAML

   - For example, using python pip command as below:
   

            pip install biopython

            pip install PyYAML

### 2.2 **d-chimer custom data :**
The taxonomic information is need in *****step 2.4 above.*****
It is available at https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz

- Using unix wget command (for example) to download:  

      wget "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz"
      
 The file named "new_taxdump.tar.gz", takes few seconds to be downloded.
    
 - Decompress it using, for example: 
    
       tar xzf new_taxdump.tar.gz
       
 - Create a tab delimited taxo file:
 
        cat fullnamelineage.dmp|sed "s/\s\+\|\s\+//g"|sed 's/\|//2g'|awk 'BEGIN{FS="|"}{print $1,$3,$4,$2}'|sed "s/\s/\t/g"|sort -k1,1 > fullnamelineage_taxid_sorted.dmp


#### - SOME CAUTIONS:

##### /!\ DOUBLE BLASTX-call : 
we used a custom database for the BLASTx version. 

In fact, using the ncbi nr database directly were too long for a BLASTx search. 

For viromes, a database of only viruses were created and deduplicated using cd-hit (http://weizhong-lab.ucsd.edu/cd-hit/ or https://github.com/weizhongli/cdhit). 
Using the d-chimer BLASTx version provided here requires a two ways BLASTx. 

The solutions are :
- Solution 1: users can create such database or a similar one according to their needs.
   - download proteins from ncbi and dereplicate using cd-hit :
 
        cd-hit -i input_viruses.fna -o viruses_nr_prots.100p.faa -c 1.00 -t 1
      
      where :
      
      -i : input file
      
      -o : output file
      
      -c : %identity `number of identical amino acids in alignment divided by the full length of the shorter sequence`
      
      -t : tolerance for redundance
   

- or Solution 2: change the code (*call_blastx_and_filter* method in *d-chimer_methods.py* file). The corresponding “yaml” configuration file have to be changed accordingly (see *****section 3 ***** ).

- for review need, we provide the BLAST formatted viral proteins database used in the d-chimer paper it is available here : https://drive.google.com/file/d/1Pmv4Rt6bFf5pgO-9k_HGHw2qEUSxsxTz/view?usp=sharing.

Once the zip file downloaded and decompressed, users must specify its path in the *d-chimer_config.yaml* file.  

##### /!\ Biopython based blast+ is used for ease of installation.
In fact, if Biopython is installed along all databases, d-chimer can be used after configuring the "yaml" file. 
The biopython blast calls may fail with blast v5 databases but work find with V4 ones. When using biopython embedded BLAST+ programs, please use the option "-L False".
If using v5 databases, then one have to install local BLAST+ programs, configure accordingly the yaml file and use the "-L True" option when running d-chimer.


##### /!\ BLAST V5 database :

users must provide the path to the local installation of BLAST in the configuration file d-chimer_config.yaml (see exemple below).


## 3. How to configure d-chimer (d-chimer_config.yaml) ? 
- Fill the different fields in the *d-chimer_config.yaml* file :

        blastn_parameters :
          dbpath_nt : /home/user/db/ncbi/nt/nt
          nb_threads_bn : 28
          evalue_nt : 0.01

        filter_blastn_parameters (-d is the shift parameter, -l is the minimum alignment length allowed -I is and identifiant over which cycles are incremented):
          d : 10
          l : 50
          I : 1

        blastx_parameters :
          dbpath_vrl : /home/user/db/ncbi/virnr/viruses_nr_prots.100p.faa
          dbpath_nr : /home/user/db/ncbi/nr/nr
          nb_threads_bx : 28
          evalue_vir : 0.1
          evalue_nr : 0.01

        filter_blastx_parameters (-d is the shift parameter, -l is the minimum alignment length allowed -I is and identifiant over which cycles are incremented):
          d : 10
          l : 17
          I : 1

        add_taxo_parameters :
          tax_lineages_file : /home/user/db/fullnamelineage_taxid_sorted.dmp
         
          blast_path : /home/user/tools/ncbi-blast-2.12.0+/bin

## 4. **How does d-chimer proceed ? (d-chimer overview scheme)**
d-chimer handles the chimeric sequences by using:
### 4.1 ***BLAST:*** 
BLAST is used to make homology search of input sequences against a reference database; either nucleotides or proteins database.
### 4.2 ***A filter:*** 
The d-chimer filter analyses coordinates of subjects aligned by BLAST on queries to build stacks of subjects aligning at the same regions of the query. The top scoring subject for each stack is kept.
The regions without any subject (uncovered zones) are cut and saved in a new fasta file. They will be re-submitted automatically to BLAST.

<img src="../main_d/d-chimer_filter.png?raw=true" class="left" >

### 4.3 ***Recursive execution of BLAST and the filter:***
Uncovered zones produced are taken and submitted to BLAST and filtered. The process ends when no uncovered zones are found or no BLAST hit is produced.

<img src="../main_d/d-chimer_recursion.PNG?raw=true" class="left" >

### 4.4 ***Taxonomic Information add to filtered BLAST outputs:***
After this process completes, filtered outputs (from the filter) are joined to the taxonomic information (using `fullnamelineage_taxid_sorted.dmp` file ; see section 3).

### 5. Appendix : INSTALL BLAST+ programs and databases
Users can either use custon databases and configure conveniently database paths in yaml files and taxonomic information as indicated in section 3.

BLAST+ can be installed, by downloading binaries here : https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/

download with wget : wget https://ftp.ncbi.nih.gov/blast/executables/blast+/

BLAST databases can be downloaded from this link : https://ftp.ncbi.nlm.nih.gov/blast/db/


If blast+ programs are installed, ncbi nucleotide and protein databases can be downloaded by using the following commands :
- nucleotide database

      update_blastdb.pl nt
     
- protein database

      update_blastdb.pl nr
