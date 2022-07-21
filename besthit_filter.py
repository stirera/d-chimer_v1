#!/usr/bin/python3
# -*-coding:Latin-1 -*

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


import argparse
import datetime
import os
import sys

from Bio import SeqIO, SeqRecord
from Bio.Seq import Seq

import Contig
import CsvIO
import Subject

from dchimer_methods import filter_besthit_blastn, filter_besthit_blastx

def opts_and_args():
    description = 'Postprocessing of blast csv files using best scoring match'
    usage = """usage : besthit_filter [-h] -p blastprogram -f fastafile -c input
            csvfile -l length -I id"""
    parser = argparse.ArgumentParser(description=description,
                                     epilog=usage)
    parser.add_argument('-p',
                        '--program',
                        type=str,
                        help='blastProgram (required) : blastn or blastx',
                        dest='blastprogram',
                        metavar='blastprogram',
                        required=True,
                        default='blastn')
    
    parser.add_argument('-f',
                        '--fasta',
                        type=str,
                        help='fasta file used as BLAST input',
                        dest='fastafile',
                        metavar='fastafile',
                        required=True,
                        default='')
    parser.add_argument('-c',
                        '--csv',
                        type=str,
                        help='csv file, BLAST output',
                        dest='csvfile',
                        metavar='csvfile',
                        required=True,
                        default='')
    parser.add_argument('-l',
                        '--length',
                        type=int,
                        help="""length (default 50): i.e minimum length of an alignment
                             """,
                        dest='minBpLength',
                        metavar='minBpLength',
                        required=False,
                        default='50')
    parser.add_argument('-I',
                        '--Id',
                        type=int,
                        help='A number to identify files',
                        dest='recursId',
                        metavar='idrecursif',
                        required=True,
                        default=0)
    args = parser.parse_args()
    return(args)


def main():
    args = opts_and_args()
    root = args.fastafile.split(".")
    minDec = 0 # set for functionnality as not applicable to besthit
    if args.blastprogram == 'blastn': 
        repertoire = root[0] + "_bh_bn_out"

        if os. path. isdir(repertoire):

            sys.exit('ERROR : BLASTn output directory exists : ' +
                     repertoire +
                     '\n\n' + 'Please remove it before'+'\n')

        elif args.minBpLength < 50 :

            input('WARNING !!!! : CHECK THE -l OPTION VALUE : '
                     '\nOR\n' + 'PRESS ENTER TO CONTINUE'+'\n')    

        filter_besthit_blastn(args.fastafile, args.csvfile, args.recursId, minDec, args.minBpLength, args.blastprogram)

    if args.blastprogram == 'blastx': 
        repertoire = root[0] + "_bh_bx_out"

        if os. path. isdir(repertoire):

            sys.exit('ERROR : BLASTx output directory exists : ' +
                     repertoire +
                     '\n\n' + 'Please remove it before'+'\n')

        elif args.minBpLength < 17 :

            input('WARNING !!!! : CHECK THE -l OPTION VALUE : '
                     '\nOR\n' + 'PRESS ENTER TO CONTINUE'+'\n')  
        filter_besthit_blastx(args.fastafile, args.csvfile, args.recursId, minDec, args.minBpLength, args.blastprogram)

if __name__ == '__main__':
    main()
