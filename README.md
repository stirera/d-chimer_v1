# d-chimer_review
## 1. **What is d-chimer for ?**
Disentangle Chimeric (d-chimer) sequences in de novo assembled viromes/metagenomes is a BLAST based pipeline conceived for taxonomic assignments by taking into account that :
- Contigs can be only partly covered in a single BLAST search 
- Contigs can be covered by different genes from different organisms

## 2. **How does d-chimer proceed ?**
d-chimer handles the chimeric sequences by using:
### 2.1 BLAST 
BLAST is used to make homology search of input sequences against a reference database; either nucleotides or proteins database.
### 2.2 ***A filter:*** 
The d-chimer filter analyses coordinates of subjects aligned by BLAST on queries to build stacks of subjects aligning at the same regions of the query. The longuest and top scoring subjects for each stack are kept.
The regions without any subject (uncovered zones) are cut and saved in a new fasta file. They will be submitted automatically to BLAST.

### 2.3 ***Recursive execution of BLAST and the filter:***
Uncovered zones produced are taken and submitted to BLAST and filtered. The process ends when no uncovered zones are found or no BLAST hit is produced.

*****[precisions par rapport à l'échange de mail de la semaine dernière.... la recursivité est codée dans la fonction call_blast_and_filter() dans d-chimer_methodes.py]*****

### 2.4 ***Taxonomic Information add to filtered BLAST outputs:***
After this process completes, filtered outputs (from the filter) are joined to the taxonomic information.




d-chimer is provided here with a BLASTn version (using nucleotide reference database) and a BLASTx version (using protein database)

## 3. **How to install d-chimer on your own system ?**
Clone or download the d-chimer repository into your system. d-chimer depends on several python3 libraries and ncbi databases.
### 3.1 Python libraries
- Biopython

   - For example, using python pip command as below:
   

            pip install biopython


- Yaml

   - For example, using python pip command as below:


            pip install PyYAML



### 3.2 BLAST databases
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

For viromes, a database of only viruses were created and deduplicated using cd-hit (http://weizhong-lab.ucsd.edu/cd-hit/). 
Using the d-chimer BLASTx version provided here requires a two ways BLASTx. 

The solutions are :
- ousers can create such database or a similar one according to their needs.

- or change the code (*call_blast_and_filter* method in *d-chimerr_methods.py* file). The corresponding “yaml” configuration file have to be changed accordingly (see *****section 3.3***** ).

- for review need, we provide the BLAST formatted viral proteins database used in the d-chimer paper it is available here : https://drive.google.com/file/d/1Pmv4Rt6bFf5pgO-9k_HGHw2qEUSxsxTz/view?usp=sharing.

Once the zip file downloaded and decompressed, users must specify its path in the *blastx_config.yaml* file.  

##### /!\ Biopython based blast+ is used for ease of installation.
In fact, if Biopython is installed and databases + taxonomic information available (*fullnamelineage_taxid_sorted.dmp* file), d-chimer can be used after configuring the "yaml" files.

## 3.3 How to configure d-chimer on your own system ?
- blastn configuration (*blastn_config.yaml* file):

      blastn_parameters :
        - dbpath=/path_to/db/nt
        - nb_threads=56
        - evalue=0.01

      filter_blastn_parameters :
        - d=10
        - l=50
        - I=1

      add_taxo_parameters :
        - tax_lineages_file=/path_to/db/fullnamelineage_taxid_sorted.dmp


- blastx configuration (*blastx_config.yaml* file):

      blastx_parameters :
         - dbpath_vrl=/path_to/db/virnr/viruses_nr_prots.100p.faa
         - dbpath_nr=/path_to/db/nr/nr
         - nb_threads=28
         - evalue=0.1
         - evalue1=0.01

      filter_blastx_parameters :
         - d=10
         - l=17
         - I=1

      add_taxo_parameters :
         - tax_lineages_file=/path_to/db/fullnamelineage_taxid_sorted.dmp


## 4. **How to run d-chimer on your own system  ?**
 
Once the parameters are configured, the d-chimer programs can be used as follow :

    python3.7 d-chimer/blastn/d-chimer_main.py -f query_file.fasta

or

    python3.7 d-chimer/blastx/d-chimer_main_bx.py -f query_file.fasta 

for blastx.


By default, d-chimer runs in a recursive manner. The number of cycles can be limited by using the -m option. For example,  to limit d-chimer to two cycles :

      python3.7 d-chimer/blastn/d-chimer_main.py -f query_file.fasta -m 2
