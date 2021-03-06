#!/usr/bin/python3
# -*-coding:Latin-1 -*


"""
Licence header
--------------
dchimer : python3 main program of d-chimer pipeline.
Copyright (C) 2022, INSTITUT PASTEUR DE LA GUYANE
This program is free software: you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by the Free Software Foundation,
either version 3 of the License, or (at your option) any later version.
This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
See the GNU General Public License for more details.
You should have received a copy of the GNU General Public License along with this program.
If not, see <https://www.gnu.org/licenses/>.
date: 19.04.2022
--
Author : Sourakhata Tirera <stirera@pasteur-cayenne.fr>
Author : Jean-Marc Frigerio <jean-marc.frigerio@inrae.fr>
Author : Alix de Thoisy <alixdet@protonmail.com>
"""

import os
import sys

import datetime
import shutil
import subprocess
import yaml
import re

from Bio import SeqIO
from Bio.Blast.Applications import NcbiblastxCommandline,\
                                   NcbiblastnCommandline
import CsvIO

# End of imports

path = os.path.dirname(os.path.abspath(__file__))
config = yaml.safe_load(open(path+"/dchimer_config.yaml"))

# BLASTn and nt database parameters
ntdbpath = config['blastn_parameters']['dbpath_nt']
nb_threads_bn = config['blastn_parameters']['nb_threads_bn']
evalue_nt = config['blastn_parameters']['evalue_nt']

# BLASTn filter parameters
filter_d_bn = config['filter_blastn_parameters']['d']
filter_l_bn = config['filter_blastn_parameters']['l']

# BLASTx parameters
vdbpath = config['blastx_parameters']['dbpath_vrl']
nrdbpath = config['blastx_parameters']['dbpath_nr']
nb_threads_bx = config['blastx_parameters']['nb_threads_bx']
evalue_vir = config['blastx_parameters']['evalue_vir']
evalue_nr = config['blastx_parameters']['evalue_nr']

# BLASTx filter parameters
filter_d_bx = config['filter_blastx_parameters']['d']
filter_l_bx = config['filter_blastx_parameters']['l']

# taxa lineages file : taxids with fullpaths
taxlineages = config['add_taxo_parameters']['tax_lineages_file']

blast_path = config['blast_path']


def runblast(program, local, qfile, db, evalue, outcsv):
    """This method takes a query fasta file, and parameters
    and runs a blast(n/x) over it
    """
    fmtbn = "6 qseqid sacc pident length mismatch gaps qstart qend " + \
            "sstart send evalue bitscore staxid sskingdom qlen qseq sseq"

    fmtbx = "6 qseqid sacc pident length mismatch gaps qstart qend " + \
            "sstart send evalue bitscore staxid sskingdom sframe " + \
            "qframe qlen qseq sseq"

    # Here are the BLAST commands
    if program == 'blastn' :
        if local is True :

          cmd = [str(blast_path)+"/" + program,
              "-task", program,
              "-query", str(qfile),
              "-db", str(db),
              "-num_threads", str(nb_threads_bn),
              "-evalue", str(evalue),
              "-culling_limit", str(10),
              "-num_alignments", str(10),
              "-outfmt", str(fmtbn),
              "-out", str(outcsv)]

          subprocess.run(cmd)

        else :
            cline=(NcbiblastnCommandline(cmd='blastn',
                                         out=outcsv,
                                         outfmt=fmtbn,
                                         query=qfile,
                                         db=db,
                                         evalue=evalue,
                                         num_threads=nb_threads_bn,
                                         culling_limit=10 ,
                                         num_alignments=10 ,
                                         task='blastn'))
            cline()

    elif program == 'blastx' :
      # this have to be modified so that users can choose to
      # run directly over nr.
        if local is True :

          cmd = [str(blast_path)+'/'+ program,
              "-query", str(qfile),
              "-db",    str(db),
              "-num_threads", str(nb_threads_bx),
              "-evalue", str(evalue),
              "-culling_limit", str(10),
              "-num_alignments", str(10),
              "-outfmt", str(fmtbx),
              "-out", str(outcsv)]

          subprocess.run(cmd)

        else :
            cline=(NcbiblastxCommandline(cmd=program,
                                     out=outcsv,
                                     outfmt=fmtbx,
                                     query=qfile,
                                     db=db,
                                     evalue=evalue,
                                     num_threads=nb_threads_bx,
                                     culling_limit=10,
                                     num_alignments=10))

            cline()


def extract_idLists_fromFasta(seqlist, fasta, listfasta):
    """Takes a list of sequence ids and a fasta file
    and creates an output fasta with only the sequence ids
    contained in the list.
    """
    with open(seqlist, "r") as listfh, \
         open(fasta, "r") as fastafh, \
         open(listfasta, "w") as outfastah :

        fastadict = SeqIO.to_dict(SeqIO.parse(fastafh, "fasta"))

        for line in listfh:
            k = re.sub('\n', '', line)
            k2 = re.sub('^', '>', fastadict[k].id)
            outfastah.write(k2+"\n")
            outfastah.write(str(fastadict[k].seq)+"\n")


def gather_matchingSequences(outvircsv, qfile, vposqfile):
    """get the list of contigs matching viral database and
    create a fasta file with those squences.
    """
    vposlist=".".join(outvircsv.split(".")[0:-1]) + ".list"
    tmp=open("tmplist", "w+")
    list=open(vposlist, "w")
    p1=subprocess.run(['cut', '-f', '1', outvircsv], stdout=subprocess.PIPE)
    tmp.write(p1.stdout.decode("utf-8")); tmp.close()
    p2=subprocess.run(['sort', '-u', './tmplist'], stdout=subprocess.PIPE)
    list.write(p2.stdout.decode("utf-8")) ; list.close()
    subprocess.run(['rm', './tmplist'])
    extract_idLists_fromFasta(vposlist, qfile, vposqfile)


def taxo_postprocess(program, root, cpt, filt_type):
    """
    This function adds taxonomy to the filtered files
    using a bash script. It also manages files :
    empty files are removed and results to the output repository
    This function is called at every stage where
    the the blast-filter recursive process can end.
    """
    if program == 'blastn' :
        blasttag = 'bn'
    elif program == 'blastx' :
        blasttag = 'bx'
    
    # filt_type can be recursive (=1)
    if (filt_type == 1) :
      c = 1
      repertoire = root[0] + "_" + blasttag + "_out"
      os.makedirs(repertoire)
      while c <= cpt :
          tsv = root[0] + "." + str(c) + "." + blasttag + ".filtered.tsv"
          tax = root[0] + "." + str(c) + "." + blasttag + ".filtered.taxo"
          aln = root[0] + "." + str(c) + "." + blasttag + ".aln"

          if os.stat(tsv).st_size == 0:
              os.remove(tsv)

          else:
              subprocess.run(['/usr/bin/bash', path+'/add_taxo.sh',
                             tsv, taxlineages, blasttag], check=True)
              shutil.move(tsv, repertoire)
              shutil.move(tax, repertoire)

          if os.stat(aln).st_size == 0:
              os.remove(aln)

          else:
              shutil.move(aln, repertoire)
              c += 1
              timeend = datetime.datetime.now()
          print("end of the whole " + program + " with sample :",
                root[0], " at :\n", timeend)

    if (filt_type == 0) :
      repertoire = root[0] + "_bh_" + blasttag + "_out"
      os.makedirs(repertoire)
      tsv = root[0] + "." + str(cpt) + ".bh." + blasttag + ".filtered.tsv"
      tax = root[0] + "." + str(cpt) + ".bh." + blasttag + ".filtered.taxo"
      aln = root[0] + "." + str(cpt) + ".bh." + blasttag + ".aln"

      if os.stat(tsv).st_size == 0:
        os.remove(tsv)

      else:
        subprocess.run(['/usr/bin/bash', path+'/add_taxo.sh',
                        tsv, taxlineages, blasttag], check=True)
        shutil.move(tsv, repertoire)
        shutil.move(tax, repertoire)

      if os.stat(aln).st_size == 0:
       os.remove(aln)

      else:
        shutil.move(aln, repertoire)
        timeend = datetime.datetime.now()
        print("end of the whole " + program + " with sample :",
              root[0], " at :\n", timeend)

def call_blastx_and_filter(program, local, qfile, cpt, cpt_max):
    """
    Recursive executions function of BLAST and filter.
    It uses all the fuctions above and runs cycles :
      BLASTx over viral database (proteins)
      Extract viral database matching queries
      Creates a fasta file from the matching queries
      Runs BLASTx over the nr database, using the precedent
      step file output file.
      Runs filter over the output of the precedent step
      If there is uncovered zones found during the
      filtering process, a new cycle will begin
      If there is no uncovered zones the recursive process
      will stop
      At the end of the recursive (blastx and filter), the
      next step starts auomatically. At this step,
      taxonmic information (taxpath are added to the
      BLASTx filtered outputs)
    """
    root = qfile.split(".")
    cpt += 1
    outvircsv = root[0]+"."+str(cpt)+".bx.vir.csv"
    vposlist = root[0]+"."+str(cpt)+".bx.vir.list"
    vposqfile = root[0]+"."+str(cpt)+".bx.vir.fas"
    outcsv = root[0] + "." + str(cpt) + ".bx.csv"

    print("BLASTx is running... on viral db")

    # Here is the BLASTx command on viral proteins database
    runblast(program, local, qfile, vdbpath, evalue_vir, outvircsv)
    print("BLASTx finished running... on viral db")

    if os.stat(outvircsv).st_size == 0:
        taxo_postprocess(program, root, cpt, 1)
        sys.exit("empty file : " + outvircsv + "\nNoting to do !")

    # get contigs that matched the viral protein database
    gather_matchingSequences(outvircsv, qfile, vposqfile)

    print("BLASTx is running... on NR db")

    # Here is the BLASTx command on nr database
    runblast(program, local, vposqfile, nrdbpath, evalue_nr, outcsv)
    print("BLASTx finished running... on NR db")

    print("Filtering the BLASTx outputs ...")

    filter_bx_output(vposqfile, outcsv, int(filter_d_bx),
                     int(filter_l_bx), str(cpt), program)
    print("cycle number ", cpt, " completed")

    if ((os.stat(root[0] + "." + str(cpt) + ".bx.fas").st_size == 0) or
       (cpt == cpt_max)):
        taxo_postprocess(program, root, cpt, 1)
        sys.exit("empty file or max number of cycles :" +
                 "\nNothing (more) to do !")

    call_blastx_and_filter(program, local, root[0] +
                           "." + str(cpt) +
                           ".bx.fas", cpt, cpt_max)


def call_blastn_and_filter(program, local, qfile, cpt, cpt_max):
    """
    Recursive executions function of BLAST and filter.
    It uses all the fuctions above and runs cycles :
      Run BLASTn over the nt database
      Run filter over the output of the precedent step
      If there is uncovered zones found during the
      filtering process, a new cycle will begin
      If there is no uncovered zones the process will stop
      and the subsequent step starts auomatically. At this
      step, taxonmic information (taxpath are added to the
      BLASTn filtered outputs)
    """
    root = qfile.split(".")
    cpt += 1
    outcsv = root[0] + "." + str(cpt) + ".bn.csv"

    print("BLASTn is running...")
    runblast(program, local, qfile, ntdbpath, evalue_nt, outcsv)
    print("BLASTn finished running...")

    print("Filtering the BLASTn outputs ...")
    filter_bn_output(qfile, outcsv, int(filter_d_bn), int(filter_l_bn), cpt, program)
    print("cycle number ", cpt, " completed")

    if ((os.stat(root[0] + "." + str(cpt) + ".bn.fas").st_size == 0) or
       (cpt == cpt_max)):
        taxo_postprocess(program, root, cpt, 1)
        sys.exit("empty file or max number of cycles :" +
                 "\nNothing (more) to do !")

    call_blastn_and_filter(program, local, root[0] +
                           "." + str(cpt) +
                           ".bn.fas", cpt, cpt_max)

def filter_bx_output(fastafile, csvfile, minDec, minBpLength, loopId, program) :
    '''
    After BLASTx is run over the nr protein database, the d-chimer filter
    (i.e., this function), takes the output csv and the BLASTx input fasta.
    The ouputs of this function are :
    - filtered csv file
    - uncovered zones (fasta)
    - BLASTx negative (fasta)
    - local alignments (aln)
    # --------- algorithm
    0 Get args
    1 open files

        1.1 read fasta and csv files :
            1.1.1 if contig is in  csv :
                1.1.1.1 get contigSubjects (filtered):
                1.1.1.2 split contig uncovered zones and recycle

                next
            1.1.1 else :
                    contig -> neg.fasta

            next


    2 End exit
    '''
    # Step 0 : Get arguments and create process id along start time
    
    root = fastafile.split(".")
    
    p_id = os.getpid()
    print ("filtering started with proces ID: ", p_id)

    timestart = datetime.datetime.now()

    # Step 1 : opening files with  root string shared between all files
    with open(fastafile, 'r') as fasta, \
         open(csvfile, 'r') as csv, \
         open(root[0]+"."+str(loopId)+".bx_neg.fasta", 'w') \
        as neg_fasta,\
        open(root[0]+"."+str(loopId)+".bx.fas", 'w') as uncovered, \
        open(root[0]+"."+str(loopId)+".bx.aln", 'w') as alnfile, \
        open(root[0]+"."+str(loopId)+".bx.filtered.tsv", 'w') \
            as filteredcsv:

        seqio = SeqIO.parse(fasta, "fasta")
        csvio = CsvIO.CsvIO(csv, int(minDec), loopId, program)

        for sequence in seqio:
            contig = csvio.next(sequence.id)

            if contig:
                print ("filtering query sequence : ", contig.id)
                F = 0
                for sub in contig.subject:
                    if int (sub.length) >= int(minBpLength) :
                        F=1
                        
                if F==0 :
                    SeqIO.write(sequence, neg_fasta, "fasta")
                    
                contig.filter_subjects()
                for sub in contig.selectedSubjects:
                    x = len(contig.selectedSubjects)
                    print ("found ", x, "subjects")

                    subjectdata = contig.id + "\t" + \
                        sub.staxid + "\t" + \
                        sub.sacc + "\t" + \
                        sub.pident + "\t" + \
                        sub.evalue + "\t" + \
                        sub.length + "\t" + \
                        sub.qstart + "\t" + \
                        sub.qend + "\t" + \
                        sub.sstart + "\t" + \
                        sub.send + "\n"

                    filteredcsv.write(subjectdata)

                    h = "aln : query=" + contig.id + \
                        " [ length =" + sub.qlen + "]" + \
                        " [frame =" + sub.qframe + "]" + \
                        " [" + sub.qstart + " " + \
                        sub.qend + "]" + \
                        " : subject=" + sub.sacc + \
                        " [" + sub.sstart + \
                        " " + sub.send + \
                        " " + str(sub.sens) + \
                        "]" + "\t" + sub.sskingdom + \
                        "\n" + sub.qseq + \
                        "\n" + sub.sseq + \
                        "\n----\n"

                    alnfile.write(h)

                contig.get_uncovered_zones()

                if contig.departs:
                    print ("recycling :")
                    rec_pool = contig.recycleseq(sequence.seq)
                    
                    for r in rec_pool:
                        uncovered.write(r[0]+"\n")
                        uncovered.write(str(r[1])+"\n")

            else:
                print ("query sequence not in csv, negative: ", sequence.id)
                SeqIO.write(sequence, neg_fasta, "fasta")

    print ("fin de process, ID =", p_id)
    timeend = datetime.datetime.now()
    print ("start_time= ", timestart, "\n", "end_time=   ", timeend)


def filter_bn_output(fastafile, csvfile, minDec, minBpLength, recursId, program):
    """
    After BLASTx is run over the nt nucleotides database, the d-chimer filter
    (i.e., this function), takes the output csv and the BLASTn input fasta.
    The ouputs of this function are :
    - filtered csv file
    - uncovered zones (fasta)
    - BLASTx negative (fasta)
    - local alignments (aln)
    # --------- algorithm
    0 Get args
    1 open files

        1.1 read fasta and csv files :
            1.1.1 if contig is in  csv :
                1.1.1.1 get contigSubjects (filtered):
                1.1.1.2 split contig uncovered zones and recycle

                next
            1.1.1 else :
                    contig -> neg.fasta

            next


    2 End exit
    """
    # Step 0 : Get arguments and create process id along start time
    root = fastafile.split(".")

    p_id = os.getpid()
    print("filtering started with proces ID: ", p_id)

    timestart = datetime.datetime.now()

    with open(fastafile, 'r') as fasta, \
        open(csvfile, 'r') as csv, \
        open(root[0]+"."+str(recursId)+".bn_neg.fasta", 'w') \
        as neg_fasta,\
        open(root[0]+"."+str(recursId)+".bn.fas", 'w') as uncovered, \
        open(root[0]+"."+str(recursId)+".bn.aln", 'w') as alnfile, \
        open(root[0]+"."+str(recursId)+".bn.filtered.tsv", 'w') \
            as filteredcsv:

        seqio = SeqIO.parse(fasta, "fasta")
        csvio = CsvIO.CsvIO(csv, minDec, minBpLength, program)

        for sequence in seqio:
            contig = csvio.next(sequence.id)

            if contig:
                print("filtering query sequence : ", contig.id)
                F = 0
                for sub in contig.subject:
                    if int(sub.length) >= int(minBpLength) :
                        F = 1

                if F == 0 :
                    SeqIO.write(sequence, neg_fasta, "fasta")

                contig.filter_subjects()
                for sub in contig.selectedSubjects:
                    x = len(contig.selectedSubjects)
                    print("found ", x, "subjects")

                    subjectdata = contig.id + "\t" + \
                        sub.staxid + "\t" + \
                        sub.sacc + "\t" + \
                        sub.pident + "\t" + \
                        sub.evalue + "\t" + \
                        sub.length + "\t" + \
                        sub.qstart + "\t" + \
                        sub.qend + "\t" + \
                        sub.sstart + "\t" + \
                        sub.send + "\n"

                    filteredcsv.write(subjectdata)

                    h = "aln : query=" + contig.id + \
                        " [ length =" + sub.qlen + "]" + \
                        " [" + sub.qstart + " " + \
                        sub.qend + "]" + \
                        " : subject=" + sub.sacc + \
                        " [" + sub.sstart + \
                        " " + sub.send + \
                        " " + str(sub.sens) + \
                        "]" + "\t" + sub.sskingdom + \
                        "\n" + sub.qseq + \
                        "\n" + sub.sseq + \
                        "\n----\n"

                    alnfile.write(h)

                contig.get_uncovered_zones()

                if contig.departs:

                    print("recycling :")
                    rec_pool = contig.recycleseq(sequence.seq)

                    for r in rec_pool:
                        uncovered.write(r[0]+"\n")
                        uncovered.write(str(r[1])+"\n")

            else:
                print("query sequence not in csv, negative: ", sequence.id)
                SeqIO.write(sequence, neg_fasta, "fasta")

    print("fin de process, ID =", p_id)
    timeend = datetime.datetime.now()
    print("start_time= ", timestart, "\n", "end_time=   ", timeend)


def filter_besthit_blastn(fastafile, csvfile, recursId, minDec, minBpLength, program):
    """
    0 Get args
    1 open files

        1.1 read fasta and csv files :
            1.1.1 if contig is in  csv :
                1.1.1.1 get contigSubjects (filtered):
                1.1.1.2 split contig uncovered zones and recycle

                next
            1.1.1 else :
                    contig -> neg.fasta

            next


    2 End exit
    
    # Step 0 : Get arguments and create process id along start time
    program = args.blastprogram
    args = opts_and_args()
    root = args.fastafile.split(".")
    """
    p_id = os.getpid()
    print ("filtering on best hit started with proces ID: ", p_id)

    timestart = datetime.datetime.now()

    # Step 1 : opening files with  root string shared between all files
    root = fastafile.split(".")

    with open(fastafile, 'r') as fasta, \
        open(csvfile, 'r') as csv, \
        open(root[0]+"."+str(recursId)+".bh.bn_neg.fasta", 'w') \
        as neg_fasta,\
        open(root[0]+"."+str(recursId)+".bh.bn.aln", 'w') as alnfile, \
        open(root[0]+"."+str(recursId)+".bh.bn.filtered.tsv", 'w') \
            as filteredcsv:

        seqio = SeqIO.parse(fasta, "fasta")
        csvio = CsvIO.CsvIO(csv, minDec, minBpLength, program)

        for s in seqio:
            contig = csvio.next(s.id)

            if contig:
                l=len(contig.subject)
                print ("filtering query sequence : ", contig.id)
                print (contig.id, "contained : ", l, "subjects")
                #contig.filter_subjects()
                
                sub=contig.select_besthit(contig.subject)


                subjectdata = contig.id + "\t" + \
                    sub.staxid + "\t" + \
                    sub.sacc + "\t" + \
                    sub.pident + "\t" + \
                    sub.evalue + "\t" + \
                    sub.length + "\t" + \
                    sub.qstart + "\t" + \
                    sub.qend + "\t" + \
                    sub.sstart + "\t" + \
                    sub.send + "\n"

                filteredcsv.write(subjectdata)

                h = "aln : query=" + contig.id + \
                    " [ length =" + sub.qlen + "]" + \
                    " [" + sub.qstart + " " + \
                    sub.qend + "]" + \
                    " : subject=" + sub.sacc + \
                    " [" + sub.sstart + \
                    " " + sub.send + \
                    " " + str(sub.sens) + \
                    "]" + "\t" + sub.sskingdom + \
                    "\n" + sub.qseq + \
                    "\n" + sub.sseq + \
                    "\n----\n"

                alnfile.write(h)

                print(subjectdata)
                print(h)
                   

            else:
                print ("query sequence not in csv, negative: ", s.id)
                SeqIO.write(s, neg_fasta, "fasta")

    taxo_postprocess(program, root, recursId, 0)

    print ("fin de process, ID =", p_id)
    timeend = datetime.datetime.now()
    print ("start_time= ", timestart, "\n", "end_time=   ", timeend)


def filter_besthit_blastx(fastafile, csvfile, recursId, minDec, minBpLength, program):
    """
    0 Get args
    1 open files

        1.1 read fasta and csv files :
            1.1.1 if contig is in  csv :
                1.1.1.1 get contigSubjects (filtered):
                1.1.1.2 split contig uncovered zones and recycle

                next
            1.1.1 else :
                    contig -> neg.fasta

            next


    2 End exit
    
    # Step 0 : Get arguments and create process id along start time
    args = opts_and_args()
    root = args.fastafile.split(".")
    """
    p_id = os.getpid()
    print ("filtering started with proces ID: ", p_id)

    timestart = datetime.datetime.now()

    # Step 1 : opening files with  root string shared between all files
    root = fastafile.split(".")
    with open(fastafile, 'r') as fasta, \
        open(csvfile, 'r') as csv, \
        open(root[0]+"."+str(recursId)+".bh.bx_neg.fasta", 'w') \
        as neg_fasta,\
        open(root[0]+"."+str(recursId)+".bh.bx.aln", 'w') as alnfile, \
        open(root[0]+"."+str(recursId)+".bh.bx.filtered.tsv", 'w') \
            as filteredcsv:

        seqio = SeqIO.parse(fasta, "fasta")
        csvio = CsvIO.CsvIO(csv, minDec, minBpLength, program)

        for s in seqio:
            contig = csvio.next(s.id)

            if contig:
                l=len(contig.subject)
                print ("filtering query sequence : ", contig.id)
                print (contig.id, "contained : ", l, "subjects")
                #contig.filter_subjects()
                
                sub=contig.select_besthit(contig.subject)


                subjectdata = contig.id + "\t" + \
                    sub.staxid + "\t" + \
                    sub.sacc + "\t" + \
                    sub.pident + "\t" + \
                    sub.evalue + "\t" + \
                    sub.length + "\t" + \
                    sub.qstart + "\t" + \
                    sub.qend + "\t" + \
                    sub.sstart + "\t" + \
                    sub.send + "\n"

                filteredcsv.write(subjectdata)

                h = "aln : query=" + contig.id + \
                    " [ length =" + sub.qlen + "]" + \
                    " [" + sub.qstart + " " + \
                    sub.qend + "]" + \
                    " : subject=" + sub.sacc + \
                    " [" + sub.sstart + \
                    " " + sub.send + \
                    " " + str(sub.sens) + \
                    "]" + "\t" + sub.sskingdom + \
                    "\n" + sub.qseq + \
                    "\n" + sub.sseq + \
                    "\n----\n"

                alnfile.write(h)

                print(subjectdata)
                print(h)
                   
            else:
                print ("query sequence not in csv, negative: ", s.id)
                SeqIO.write(s, neg_fasta, "fasta")

    taxo_postprocess(program, root, recursId, 0)

    print ("fin de process, ID =", p_id)
    timeend = datetime.datetime.now()
    print ("start_time= ", timestart, "\n", "end_time=   ", timeend)
