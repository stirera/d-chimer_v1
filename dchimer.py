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
"""

import os
import sys

from dchimer_methods import call_blastx_and_filter, call_blastn_and_filter

import argparse


def opts_and_args():
    """d-chimer arguments and options"""
    description = 'run a blast program and filter its outputs'
    usage = 'usage : dchimer [-h] -p blastProgram -L True -f fastafile [-m] max_loops'
    parser = argparse.ArgumentParser(description=description, epilog=usage)

    parser.add_argument('-p',
                        '--program',
                        type=str,
                        help='blastProgram (required) : blastn or blastx',
                        dest='blastprogram',
                        metavar='blastprogram',
                        required=True,
                        default='blastn')

    parser.add_argument('-L',
                        '--local',
                        type=bool,
                        help="""use local machine BLAST+ programs or Biopython
                                embedded ones (required) :
                                True = use local AND FALSE = use Biopython BLAST+ (default : True)""",
                        dest='local',
                        metavar='Tue | False',
                        required=True,
                        default=True)

    parser.add_argument('-f',
                        '--fasta',
                        type=str,
                        help='input fasta_file (required)',
                        dest='fastafile',
                        metavar='fastafile',
                        required=True,
                        default='')

    parser.add_argument('-m',
                        '--max_recursive_loops',
                        type=int,
                        help="""maximum number of times to to
                                process uncovered zones fasta""",
                        dest='max_loops',
                        metavar='max_loops',
                        required=False,
                        default=10000)
    args = parser.parse_args()
    return args


def main():
    """Recursive process (implemented in call_blast(n/x)_and_filter
    functions) which calls a blast program, filters its
    outputs until max_loops number of cycles reached or process ends
    """
    args = opts_and_args()

    program = args.blastprogram
    local = args.local
    qfile = args.fastafile

    cpt = 0
    cpt_max = args.max_loops

    root = args.fastafile.split(".")

    if program == 'blastn':
        repertoire = root[0] + "_bn_out"

        if os. path. isdir(repertoire):

            sys.exit('ERROR : BLASTn output directory exists : ' +
                     repertoire +
                     '\n\n' + 'Please remove it before'+'\n')

        call_blastn_and_filter(program, local, qfile, cpt, cpt_max)
    elif program == 'blastx':
        repertoire = root[0] + "_bx_out"

        if os. path. isdir(repertoire):

            sys.exit('ERROR : BLASTx output directory exists : ' +
                     repertoire +
                     '\n\n' + 'Please remove it before'+'\n')

        call_blastx_and_filter(program, local, qfile, cpt, cpt_max)
    else :
        sys.exit('ERROR : unknown blast program ! ' +
                 'allowed values are : \n\tblastn or blastx')


if __name__ == '__main__':
    main()
