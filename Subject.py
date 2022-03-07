#!/usr/bin/python3
# -*-coding:Latin-1 -*


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
