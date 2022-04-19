#!/usr/bin/python3
# -*-coding:Latin-1 -*

"""
Licence header
--------------
Subject class : python3 Class for handling pairs of query/subjct in d-chimer pipeline.
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


class Subject():
    """This Class juste registers the full line of
    BLASTx output, except query/contig
    identifiant, the first coulumn (bl[0])
    """

    def __init__(self, bl, program):
        """It is initialized by a BLASTx output line : bl[]
        Here it is callend by CsvIO.next() method
        """
        if program == 'blastx':

            if int(bl[15]) < 0:
                self.sacc, self.pident, self.length, self.mismatch, \
                    self.gaps, self.qend, self.qstart, self.sstart, \
                    self.send, self.evalue, self.bitscore, self.staxid, \
                    self.sskingdom, self.sframe, self.qframe, self.qlen, \
                    self.qseq, self.sseq = bl[1:19]
                self.sens = False
            else:
                self.sacc, self.pident, self.length, self.mismatch, \
                    self.gaps, self.qstart, self.qend, self.send, \
                    self.sstart, self.evalue, self.bitscore, self.staxid, \
                    self.sskingdom, self.sframe, self.qframe, self.qlen, \
                    self.qseq, self.sseq = bl[1:19]

                self.sens = True

        elif program == 'blastn' :
            self.sacc, self.pident, self.length, \
                self.mismatch, self.gaps, self.qstart, \
                self.qend, self.sstart, self.send, \
                self.evalue, self.bitscore, self.staxid, \
                self.sskingdom, self.qlen, self.qseq, self.sseq = bl[1:17]
            self.sens = True
            if int(bl[9]) < int(bl[8]):
                self.sens = False
