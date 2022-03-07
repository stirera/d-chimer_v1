#!/usr/bin/python3
# -*-coding:Latin-1 -*

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
                        help="""use local machine BLAST programs (required) :
                                True/False (default : True)""",
                        dest='local',
                        metavar='local',
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
