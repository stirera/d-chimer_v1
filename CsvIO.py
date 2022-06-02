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

import Contig
import Subject


class CsvIO():
    """
    This class takes input csv filehandle and returns contigs that
    matched in blast search.
    It returns "null" in case contig did not match or end of file.
    """

    def __init__(self, handle, minDec, minBpLength, program):
        self.handle = handle
        self.minDec = int(minDec)
        self.minBpLength = int(minBpLength)
        self.program = program

    def next(self, s_id):
        """For a given contig search its id in BLAST csv output
        Load it with subject class if exists... else return null
        """
        pos = self.handle.tell()
        bl = self.handle.readline().split()

        # No BLAST output line (end of file)
        if not bl:
            return None

        c = Contig.Contig(bl[0], self.minDec, self.minBpLength)  # init Contig
        # if Current contig (from fasta), does not appear in (output) csv file
        if c.id != s_id:
            self.handle.seek(pos)  # retour a position initiale
            return None

        c.subject.append(Subject.Subject(bl, self.program))

        # takes subsequent lines with the same query contig id
        while c.id == s_id:
            pos = self.handle.tell()
            bl = self.handle.readline().split()
            if not bl:
                return c  # end of file
            c.id = bl[0]
            c.subject.append(Subject.Subject(bl, self.program))

#       else:  # end of this contig
        c.id = s_id
        c.subject = c.subject[:-1]
        self.handle.seek(pos)

        return c
