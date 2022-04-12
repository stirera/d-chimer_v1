# d-chimer_v1
## 1. **What is d-chimer for ?**
Disentangle chimeric (d-chimer) sequences in de novo assembled viromes/metagenomes is a BLAST based pipeline conceived for taxonomic assignments by taking into account that :
- Contigs can be only partly covered in a single BLAST search 
- Contigs can be covered by different genes from different organisms

## 2. **How does d-chimer proceed ?**
d-chimer handles the chimeric sequences by using:
### 2.1 ***BLAST:*** 
BLAST is used to make homology search of input sequences against a reference database; either nucleotides or proteins database.
### 2.2 ***A filter:*** 
The d-chimer filter analyses coordinates of subjects aligned by BLAST on queries to build stacks of subjects aligning at the same regions of the query. The top scoring subject for each stack is kept.
The regions without any subject (uncovered zones) are cut and saved in a new fasta file. They will be re-submitted automatically to BLAST.

### 2.3 ***Recursive execution of BLAST and the filter:***
Uncovered zones produced are taken and submitted to BLAST and filtered. The process ends when no uncovered zones are found or no BLAST hit is produced.

*****[precisions par rapport à l'échange de mail de la semaine dernière.... la recursivité est codée dans la fonction call_blast_and_filter() dans d-chimer_methodes.py]*****

### 2.4 ***Taxonomic Information add to filtered BLAST outputs:***
After this process completes, filtered outputs (from the filter) are joined to the taxonomic information.




d-chimer is provided here with to execute BLASTn or BLASTx (using protein database) using the -p option.

## 3. **How to install d-chimer on your own system ?**
Clone or download the d-chimer repository into your system. d-chimer depends on several python3 libraries and ncbi databases.
### 3.1 Python libraries
- Biopython

   - For example, using python pip command as below:
   

            pip install biopython


- Yaml

   - For example, using python pip command as below:


            pip install PyYAML


### 3.2 BLAST+ programs
BLAST+ can be installed, by downloading binaries here : https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/

download with wget : wget https://ftp.ncbi.nih.gov/blast/executables/blast+/

### 3.3 BLAST databases
BLAST databases can be downloaded from this link : https://ftp.ncbi.nlm.nih.gov/blast/db/

If blast+ programs are installed, ncbi nucleotide and protein databases can be downloaded by using the following commands :
- nucleotide database

      update_blastdb.pl nt
     
- protein database

      update_blastdb.pl nr

The taxonomic information is need in *****step 2.4 above.*****
It is available at https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz

- Using unix wget command (for example) to download:  

      wget "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz"
      
 The file named "new_taxdump.tar.gz", takes few seconds to be downloded into you system.
    
 - Decompress it using, for example: 
    
       tar xzf new_taxdump.tar.gz
       
 - Create a tab delimited taxo file:
 
        cat fullnamelineage.dmp|sed "s/\s\+\|\s\+//g"|sed 's/\|//2g'|awk 'BEGIN{FS="|"}{print $1,$3,$4,$2}'|sed "s/\s/\t/g"|sort -k1,1 > fullnamelineage_taxid_sorted.dmp

#### - Optional:
taxdb.bti and taxdb.btd are available from https://ftp.ncbi.nlm.nih.gov/blast/db/taxdb.tar.gz.

These files add-up the "Kingdom" (i.e. Eukaryota, Bacteria, Viruses) information included in the BLAST command-line used here. 

These two files have to be put into the directory d-chimer BLAST is executed.


#### - SOME CAUTIONS:

##### /!\ DOUBLE BLASTX-call : 
we used a custom database for the BLASTx version. 

In fact, using the ncbi nr database directly were too long for a BLASTx search. 

For viromes, a database of only viruses were created and deduplicated using cd-hit (http://weizhong-lab.ucsd.edu/cd-hit/ or https://github.com/weizhongli/cdhit). 
Using the d-chimer BLASTx version provided here requires a two ways BLASTx. 

The solutions are :
- users can create such database or a similar one according to their needs.

- or change the code (*call_blast_and_filter* method in *d-chimer_methods.py* file). The corresponding “yaml” configuration file have to be changed accordingly (see *****section 3.3***** ).

- for review need, we provide the BLAST formatted viral proteins database used in the d-chimer paper it is available here : https://drive.google.com/file/d/1Pmv4Rt6bFf5pgO-9k_HGHw2qEUSxsxTz/view?usp=sharing.

Once the zip file downloaded and decompressed, users must specify its path in the *d-chimer_config.yaml* file.  

##### /!\ Biopython based blast+ is used for ease of installation.
In fact, if Biopython is installed and all databases along all databases, d-chimer can be used after configuring the "yaml" file. 
The biopython blast calls may fail with blast v5 databases but work find with V4 ones. When using biopython embedded BLAST+ programs, please use the option "-L False".
If using v5 databases, then one have to install local BLAST+ programs, configure accordingly the yaml file and use the "-L True" option when running dchimer.


##### /!\ BLAST V5 database :

users must provide the path to the local installation of BLAST in the configuration file d-chimer_config.yaml (see exemple below).


## 3.3 How to configure d-chimer on your own system (d-chimer_config.yaml) ? 
- Fill the different fields in the *d-chimer_config.yaml* file :

        blastn_parameters :
          dbpath_nt : /home/user/db/ncbi/nt/nt
          nb_threads_bn : 28
          evalue_nt : 0.01

        filter_blastn_parameters :
          d : 10
          l : 50
          I : 1

        blastx_parameters :
          dbpath_vrl : /home/user/db/ncbi/virnr/viruses_nr_prots.100p.faa
          dbpath_nr : /home/user/db/ncbi/nr/nr
          nb_threads_bx : 28
          evalue_vir : 0.1
          evalue_nr : 0.01

        filter_blastx_parameters :
          d : 10
          l : 17
          I : 1

        add_taxo_parameters :
          tax_lineages_file : /home/user/db/fullnamelineage_taxid_sorted.dmp
         
        blast_path : /home/user/tools/ncbi-blast-2.12.0+/bin



            

## 4. **How to run d-chimer on your own system  ?**
 
Once the parameters are configured (through the `d-chimer_config.yaml`, file), d-chimer can be used as follow :

    python3 dchimer.py  -p blastn -L True -f query_file.fasta

By default, d-chimer runs in a recursive manner. The number of cycles can be limited by using the -m option. For example,  to limit d-chimer to two cycles :

      python3 dchimer.py -p blastn -L True -f query_file.fasta -m 2


For full view of d-chimer options, type :


      python3 dchimer.py -h

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
 ############
 # d-chimer_v1
## 1. **What is d-chimer for ?**
Disentangle chimeric (d-chimer) sequences in de novo assembled viromes/metagenomes is a BLAST based pipeline conceived for taxonomic assignments by taking into account that :
- Contigs can be only partly covered in a single BLAST search 
- Contigs can be covered by different genes from different organisms

## 2. **How does d-chimer proceed ?**
d-chimer handles the chimeric sequences by using:
### 2.1 ***BLAST:*** 
BLAST is used to make homology search of input sequences against a reference database; either nucleotides or proteins database.
### 2.2 ***A filter:*** 
The d-chimer filter analyses coordinates of subjects aligned by BLAST on queries to build stacks of subjects aligning at the same regions of the query. The top scoring subject for each stack is kept.
The regions without any subject (uncovered zones) are cut and saved in a new fasta file. They will be re-submitted automatically to BLAST.

### 2.3 ***Recursive execution of BLAST and the filter:***
Uncovered zones produced are taken and submitted to BLAST and filtered. The process ends when no uncovered zones are found or no BLAST hit is produced.

*****[precisions par rapport à l'échange de mail de la semaine dernière.... la recursivité est codée dans la fonction call_blast_and_filter() dans d-chimer_methodes.py]*****

### 2.4 ***Taxonomic Information add to filtered BLAST outputs:***
After this process completes, filtered outputs (from the filter) are joined to the taxonomic information.




d-chimer is provided here with to execute BLASTn or BLASTx (using protein database) using the -p option.

## 3. **How to install d-chimer on your own system ?**
- Clone or download the d-chimer repository into your system.   

             git clone https://github.com/stirera/d-chimer_v1


d-chimer depends on several python3 libraries and ncbi databases.
- Python libraries
  - Biopython and PyYAML

   - For example, using python pip command as below:
   

            pip install biopython

            pip install PyYAML


The taxonomic information is need in *****step 2.4 above.*****
It is available at https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz

- Using unix wget command (for example) to download:  

      wget "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz"
      
 The file named "new_taxdump.tar.gz", takes few seconds to be downloded into you system.
    
 - Decompress it using, for example: 
    
       tar xzf new_taxdump.tar.gz
       
 - Create a tab delimited taxo file:
 
        cat fullnamelineage.dmp|sed "s/\s\+\|\s\+//g"|sed 's/\|//2g'|awk 'BEGIN{FS="|"}{print $1,$3,$4,$2}'|sed "s/\s/\t/g"|sort -k1,1 > fullnamelineage_taxid_sorted.dmp

#### - Optional:
taxdb.bti and taxdb.btd are available from https://ftp.ncbi.nlm.nih.gov/blast/db/taxdb.tar.gz.

These files add-up the "Kingdom" (i.e. Eukaryota, Bacteria, Viruses) information included in the BLAST command-line used here. 

These two files have to be put into the directory d-chimer BLAST is executed.


#### - SOME CAUTIONS:

##### /!\ DOUBLE BLASTX-call : 
we used a custom database for the BLASTx version. 

In fact, using the ncbi nr database directly were too long for a BLASTx search. 

For viromes, a database of only viruses were created and deduplicated using cd-hit (http://weizhong-lab.ucsd.edu/cd-hit/ or https://github.com/weizhongli/cdhit). 
Using the d-chimer BLASTx version provided here requires a two ways BLASTx. 

The solutions are :
- users can create such database or a similar one according to their needs.

- or change the code (*call_blast_and_filter* method in *d-chimer_methods.py* file). The corresponding “yaml” configuration file have to be changed accordingly (see *****section 3.3***** ).

- for review need, we provide the BLAST formatted viral proteins database used in the d-chimer paper it is available here : https://drive.google.com/file/d/1Pmv4Rt6bFf5pgO-9k_HGHw2qEUSxsxTz/view?usp=sharing.

Once the zip file downloaded and decompressed, users must specify its path in the *d-chimer_config.yaml* file.  

##### /!\ Biopython based blast+ is used for ease of installation.
In fact, if Biopython is installed and all databases along all databases, d-chimer can be used after configuring the "yaml" file. 
The biopython blast calls may fail with blast v5 databases but work find with V4 ones. When using biopython embedded BLAST+ programs, please use the option "-L False".
If using v5 databases, then one have to install local BLAST+ programs, configure accordingly the yaml file and use the "-L True" option when running dchimer.


##### /!\ BLAST V5 database :

users must provide the path to the local installation of BLAST in the configuration file d-chimer_config.yaml (see exemple below).


## 4. How to configure d-chimer on your own system (d-chimer_config.yaml) ? 
- Fill the different fields in the *d-chimer_config.yaml* file :

        blastn_parameters :
          dbpath_nt : /home/user/db/ncbi/nt/nt
          nb_threads_bn : 28
          evalue_nt : 0.01

        filter_blastn_parameters :
          d : 10
          l : 50
          I : 1

        blastx_parameters :
          dbpath_vrl : /home/user/db/ncbi/virnr/viruses_nr_prots.100p.faa
          dbpath_nr : /home/user/db/ncbi/nr/nr
          nb_threads_bx : 28
          evalue_vir : 0.1
          evalue_nr : 0.01

        filter_blastx_parameters :
          d : 10
          l : 17
          I : 1

        add_taxo_parameters :
          tax_lineages_file : /home/user/db/fullnamelineage_taxid_sorted.dmp
         
        blast_path : /home/user/tools/ncbi-blast-2.12.0+/bin



            

## 5. **How to run d-chimer on your own system  ?**
 
Once the parameters are configured (through the `d-chimer_config.yaml`, file), d-chimer can be used as follow :

    python3 dchimer.py  -p blastn -L True -f query_file.fasta

By default, d-chimer runs in a recursive manner. The number of cycles can be limited by using the -m option. For example,  to limit d-chimer to two cycles :

      python3 dchimer.py -p blastn -L True -f query_file.fasta -m 2


For full view of d-chimer options, type :


      python3 dchimer.py -h

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



### 6. Appendix : INSTALL BLAST+ programs and databases
Users can either use custon databases and configure conveniently database paths in yaml files and taxonomic information as indicated in section 3.

BLAST+ can be installed, by downloading binaries here : https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/

download with wget : wget https://ftp.ncbi.nih.gov/blast/executables/blast+/

### 3.3 BLAST databases
BLAST databases can be downloaded from this link : https://ftp.ncbi.nlm.nih.gov/blast/db/


If blast+ programs are installed, ncbi nucleotide and protein databases can be downloaded by using the following commands :
- nucleotide database

      update_blastdb.pl nt
     
- protein database

      update_blastdb.pl nr

 
 
