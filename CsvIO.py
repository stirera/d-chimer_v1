#!/usr/bin/python3
# -*-coding:Latin-1 -*

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
