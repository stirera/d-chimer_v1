# d-chimer_v1

The aim of d-chimer pipeline is to detect and assign fragments of contigs included in chimeric sequences during taxonomic assignment with BLAST. It first executes BLAST and filters its results and recycles unassaigned contig fragments.


d-chimer can execute a BLASTn or a BLASTx independently when searching against reference databases. In the BLASTn version, a cycle is made of i) a search against a nucleotide database (ncbi nt), ii) filtering of the outputs where recycling is done and ends with iii) taxonomic identification by adding full lineage path and names.

The BLASTx version runs by i) using a custom viral protein database ii) the matching sequences  are extracted iii) the matching sequences undergo a second BLASTx against the whole protein database (ncbi nr), iV) the filter is executed as in BLAST (step ii) and taxonomic identification done as in BLASTn (step iv).

In the the d-chimer publication, datasets were first submitted to d-chimer BLASTn and the non-matching sequences were then re-submitted to d-chimer BLASTx (see figure below).


  <img src="../img/global_d-chimer.png?raw=true" class="left">

###### figure 1 : d-chimer usage with BLASTn version followed by BLASTx.


## 1. Usage
 
Once d-chimer is installed (see section 3 for installation) and the parameters configured (through the `d-chimer_config.yaml`file, section 4), it can be used as follow *:

    python3 /path/to/d-chimer/dchimer.py  -p blastn -L True -f query_file.fasta
or

    python3 /path/to/d-chimer/dchimer.py  -p blastx -L True -f query_file.fasta
    
\* We recommand users to execute d-chimer in the directory where they input data are. 
    

By default, d-chimer runs in a recursive manner. The number of cycles can be limited by using the -m option. For example,  to limit d-chimer to two cycles :

      python3 /path/to/d-chimer/dchimer.py  -p blastn -L True -f query_file.fasta -m 2


For full view of d-chimer options, type :



      python3 /path/to/d-chimer/dchimer.py -h


      usage: dchimer.py [-h] -p blastprogram -L True | False -f fastafile
                  [-m max_loops]

      run a blast program and filter its outputs

      optional arguments:
        -h, --help            show this help message and exit
        -p blastprogram, --program blastprogram
                              blastProgram (required) : blastn or blastx
        -L True | False, --local True | False
                              use local machine BLAST+ programs or Biopython
                              embedded ones (required) : True = use local AND False
                              = use Biopython BLAST+
        -f fastafile, --fasta fastafile
                              input fasta_file (required)
        -m max_loops, --max_recursive_loops max_loops
                              maximum number of times to process uncovered zones
                              fasta


## 2.  Usage exemple with the provided input file "testSequences.fasta"
### 2.1 run d-chimer blastn with the input file :
      
           python3 /path/to/d-chimer/dchimer.py -p blastn -L True -f testSequences.fasta

  for each cycle, will be written :
 - The blastn hits : *testSequences.1.bn.csv* (cycle 1 blastn hits).
  It is a 17 column tabultated BLASTn output file : the 12 firts are default BLAST outputs + five columns standing respectively for *taxid*, *kingdom*, *query length*,   *query sequence alignment* and *subject sequence alignment* 
 - The uncovered_zones : *testSequences.1.bn.fas*
 - The blastn negative sequences : *testSequences.1.bn_neg.fasta*
 - A file to ensure that there no mismatches during taxonomy joining to filtered outputs (must be empty, if correct): *testSequences.1.bn.filtered.notfound_taxo.tsv*.
 - A folder containing the final outputs, named : *testSequences_bn_out/*
 - The filtered outputs : *testSequences.1.bn.filtered.sorted.tsv*
  
    It is a ten column file which contains contains : *contig identifiant* (line 1), *taxonomic identifiant* (line 2), *subject accession number* (line 3), *identity percent* (line 4), *e-value* (line 5),          *alignment length* (line 6), *query alignment start coordinate* (line 7), *query alignment end coordinate* (line 7), *subject alignment start coordinate* (line 8), *subject alignment end coordinate* (line 8).

 <img src="../img/d_chimer_test_bn_output2.PNG?raw=true" class="left" >
 
 ###### figure 2 : description of filtered outputs of BLAST before taxonomic paths are added.
        
 
  - Taxomony associated filered output (see the image below) : *testSequences.1.bn.filtered.taxo*
   
   It is an eleven column file which contains contains : *taxonomic identifiant* (line 1), *contig identifiant* (line 2), *subject accession number* (line 3), *identity percent* (line 4), *e-value* (line 5), *alignment length* (line 6), *query alignment start coordinate* (line 7), *query alignment end coordinate* (line 7), *subject alignment start coordinate* (line 8), *subject alignment end coordinate* (line 8) and *full taxomic path* (line 9)
   
  <img src="../img/d_chimer_test_bn_output.PNG?raw=true" class="left" >

###### figure 3 : description of filtered outputs of BLAST with taxonomic paths.

   - query/subjects alignments : *testSequences.1.bn.filtered.aln*
   Its a three lines item with a description followed by two lines each representing the BLAST aligned sequence for contig and subject.

### 2.2 run d-chimer blastx with the input file :

  A similar set of file are produced with blastx... everywhere "bn" is replaced by "bx". Additional intermediate files will be produced :

   - BLASTx against viral database output : *testSequences.bx.vir.csv*  
   - List of matching sequences against the viral database : *test_1.1.bx.vir.list*
   - Fasta file of sequneces matching the viral database : *testSequences.bx.vir.fas*

d-chimer produces all the outputs in the repository where it was launched.


## 3. Installation
Users can execute the `INSTALL.sh` bash script once this repo is cloned, to get d-chimer source code, the python libraries and custom data.
### 3.1 d-chimer code and python libraries : 
- Clone or download the d-chimer repository.   

             git clone https://github.com/stirera/d-chimer_v1


d-chimer depends on several python3 libraries and ncbi-BLAST and databases.
- Python libraries
  - Biopython and PyYAML

   - For example, using python pip command as below:
   

            pip install biopython

            pip install PyYAML

### 3.2 **d-chimer custom data :**
#### 3.2.1 taxonomic lineages :
The taxonomic full lineage paths are needed for d-chimer to complete its run correctly ***** see also section 5.4 below fo details.*****
It is available at https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz

- Using unix wget command (for example) to download:  

      wget "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz"
      
 The file named "new_taxdump.tar.gz", takes few seconds to be downloaded.
    
 - Decompress it using, for example: 
    
       tar xzf new_taxdump.tar.gz
       
 - Create a tab delimited taxo file :
 
 The command used here reformats the original NCBI `fullnamelineage.dmp` file into a `taxid sorted` tabulated file, named `fullnamelineage_taxid_sorted.dmp`.
 
        cat fullnamelineage.dmp|sed "s/\s\+\|\s\+//g"|sed 's/\|//2g'|awk 'BEGIN{FS="|"}{print $1,$3,$4,$2}'|sed "s/\s/\t/g"|sort -k1,1 > fullnamelineage_taxid_sorted.dmp

#### 3.2.2 custom viral database.
Using the ncbi nr database directly were too long for BLASTx searche, we used a custom database for the BLASTx version. 
It is made by downloading all proteing sequences in ncbi nr protein database.
The produced fasta file was deduplicated using cd-hit (http://weizhong-lab.ucsd.edu/cd-hit/ or https://github.com/weizhongli/cdhit). 
Using the d-chimer BLASTx version provided here requires a two ways BLASTx.
 
        cd-hit -i input_viruses.fna -o viruses_nr_prots.100p.faa -c 1.00 -t 1
      
      where :
      
      -i : input file
      
      -o : output file
      
      -c : %identity `number of identical amino acids in alignment divided by the full length of the shorter sequence`
      
      -t : tolerance for redundance
   

-  We provide the BLAST formatted viral proteins database used in the d-chimer paper it is available here : https://drive.google.com/file/d/1Pmv4Rt6bFf5pgO-9k_HGHw2qEUSxsxTz/view?usp=sharing.

Once the zip file downloaded and decompressed, users must specify its path in the *d-chimer_config.yaml* file (section 4, below).  


## 4. How to configure d-chimer (d-chimer_config.yaml) ? 
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

## 5. **How does d-chimer proceed ? (d-chimer overview scheme)**
d-chimer handles the chimeric sequences by using:
### 5.1 ***BLAST:*** 
BLAST is used to make homology search of input sequences against a reference database; either nucleotides or proteins database.
### 5.2 ***A filter:*** 
The d-chimer filter analyses coordinates of subjects aligned by BLAST on queries to build stacks of subjects aligning at the same regions of the query. The top scoring subject for each stack is kept.
The regions without any subject (uncovered zones) are cut and saved in a new fasta file. They will be re-submitted automatically to BLAST.

<img src="../img/d-chimer_filter.png?raw=true" class="left" >

###### figure 4 : d-chimer filter description

### 5.3 ***Recursive execution of BLAST and the filter:***
Uncovered zones produced are taken and re-submitted to BLAST and filtered, in a new cycle. The process ends when no uncovered zones are found or no BLAST hit is produced or maximum number cycles is reached.

<img src="../img/d-chimer_recursion.PNG?raw=true" class="left" >

###### figure 5 : d-chimer recursive execution 

### 5.4 ***Taxonomic Information add to filtered BLAST outputs:***
After this process completes, filtered outputs (from the filter) are joined to the taxonomic information (using `fullnamelineage_taxid_sorted.dmp` file ; see section 3.2.1).

## 6. Appendix : INSTALL BLAST+ programs and databases
Users can either use custon databases and configure conveniently database paths in yaml files and taxonomic information as indicated in section 3.

BLAST+ can be installed, by downloading binaries here : https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/

download with wget : wget https://ftp.ncbi.nih.gov/blast/executables/blast+/

BLAST databases can be downloaded from this link : https://ftp.ncbi.nlm.nih.gov/blast/db/


If blast+ programs are installed, ncbi nucleotide and protein databases can be downloaded by using the following commands :
- nucleotide database

      update_blastdb.pl nt
     
- protein database

      update_blastdb.pl nr

##### /!\ Biopython based BLAST+ is used for ease of installation but cannot handle the ncbi V5 databases.
The biopython blast calls fail with blast v5 databases but works fine with V4 ones. When using biopython embedded BLAST+ programs, please use the option "-L False".
To use v5 databases, one have to install local BLAST+ programs and databases, configure accordingly the `yaml` file by providing absolute path to the local installation of BLAST+ binaries repository and databases ;  the "-L True" option must be used when running d-chimer.

## 7. best-hit filter
A best-hit filter is available `besthit_filter.py`. It chooses the best-scoring match subject for each query sequence in the csv file. For that a fasta file (queries used with BLASTn/x), the resulting output csv file and produces the same outputs as d-chimer at a given cycle (cf `section 2` of this document). 
It shares the  `-p, -f -c -l -I` options with the same significances as d-chimer.

      python3 /path/to/d-chimer/besthit_filter.py -h
      usage: besthit_filter.py [-h] -p blastprogram -f fastafile -c csvfile [-l minBpLength] -I idrecursif

      Postprocessing of blast csv files using best scoring match

      optional arguments:
        -h, --help            show this help message and exit
        -p blastprogram, --program blastprogram
                              blastProgram (required) : blastn or blastx
        -f fastafile, --fasta fastafile
                              fasta file used as BLAST input
        -c csvfile, --csv csvfile
                              csv file, BLAST output
        -l minBpLength, --length minBpLength
                              length (default 50): i.e minimum length of an alignment
        -I idrecursif, --Id idrecursif
                              A number to identify files

      usage : besthit_filter [-h] -p blastprogram -f fastafile -c input csvfile -l length -I id
