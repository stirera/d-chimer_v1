#!/usr/bin/python3
# -*-coding:Latin-1 -*


class Contig():
    """This class corresponds to query contig object at the the end of blast.
    It is associated with subjects as a table (see subject class)
    It contains variables for subjects' filtering along the filtering variables
    minDec and minBpLength.
    Its most important methods are Contig.filter_subjects()
    and Contig.get_uncovered_zones().
    """

    def __init__(self, id, minDec, minBpLength):
        self.id = id  # contig or query id
        self.subject = []  # array of subject(s) (class subjects)
        self.minBpLength = int(minBpLength)
        '''staight from the user defined variable of the filter_main -l option
        = minimum alignment length'''
        self.minDec = int(minDec)
        '''staight from the user defined variable of the filter_main -d option'''
        self.selectedSubjects = []  # array of the subjects after filtration
        self.departs = []  # start coordinates of subjects
        self.arrivees = []  # end of alignments of subjects

    def filter_subjects(self):
        """This function produces the final set of subject retained
        after filtering. It uses contig.build_stacks,
        it builds a matrix of subjects with stacks,
        overlapping, or single subjects.
        Using the contig.select_besthit(), it chosees the best-hit
        for each stack. It is called in in filter_blastn_lib.
        It calls : Contig.sort_subjects_on_starting(), to first
        sort subjects  on their starting coordiantes befiore starting
        stacking procedure.
        """
        new_subjects = self.sort_subjects_on_starting(self.subject)

        classes = []
        idtab = []
        smatrix = []
        for s, _element in enumerate(new_subjects):
            idtab = [new_subjects[s].sacc,
                     new_subjects[s].qstart,
                     new_subjects[s].qend]
            if idtab not in classes:
                classes.append(idtab)  # seen ?
                if len(smatrix) >= 1:
                    t = self.build_stacks(new_subjects[s], smatrix)
                    if t[0] == 1:
                        smatrix[t[1]].append(new_subjects[s])
                    else:
                        smatrix.append([new_subjects[s]])
                else:
                    smatrix.append([new_subjects[s]])

        # Finish with selecting best-hits
        for a, _element in enumerate(smatrix):
            if len(smatrix[a]) > 0:
                s = self.select_besthit(smatrix[a])
                self.selectedSubjects.append(s)
            else:
                self.selectedSubjects.append(smatrix[a])

    def sort_subjects_on_starting(self, subjects):
        """Given a contig object with all its subjects:
        filter on length (min = minBpLength)
        Order according to alignment starting coordinates.
        """
        deps = []
        new_subjects = []
        for s, _element in enumerate(subjects):
            if int(subjects[s].length) >= int(self.minBpLength):
                deps.append([int(subjects[s].qstart), subjects[s]])

        sorted_deps = sorted(deps, key=lambda d: int(d[0]))
        for d in sorted_deps:
            new_subjects.append(d[1])

        return new_subjects

    def build_stacks(self, subject, matrix):
        """Given a subject and a matrix of subjects :
        Compare coordinates and determine if its alignment
        coordinates are covered by subject(s) in the matrix (+/- minDec)
        fusion is a flag variable set to 0, If covered flag = 1
        Returns fusion, coordinates in the matrix and subject as a table
        """
        fusion = x = y = 0
        tabfus = []
        for x, _element in enumerate(matrix):
            for y, _element in enumerate(matrix[x]):
                dep = int(subject.qstart)
                arriv = int(subject.qend)
                # cell of the matrix => a subject
                borders = [int(matrix[x][y].qstart), int(matrix[x][y].qend)]

                # matrix element have a starting and ending coordinates
                # framed (by > minDec) by subject
                if (((dep - self.minDec) > borders[0]) and
                   ((arriv + self.minDec) < borders[1])):
                    fusion = 1
                    tabfus = [fusion, x, y, subject]
                    return tabfus

                m = 0
                while m <= self.minDec:
                    # from 0 to minDec (Shift parameter):
                    # search if border of (qstart/qend) fit in,
                    # in other words searching for overlaps <= minDec value
                    # print ("build_stacks loop number 3 while", y)
                    if ((int(dep-m) == borders[0]) or
                        (int(dep+m) == borders[0]) or
                       (int(arriv-m) == borders[1]) or
                       (int(arriv+m) == borders[1])):
                        fusion = 1
                        # print "valeur de fusion if fusion :", fusion
                        tabfus = [fusion, x, y, subject]
                        return tabfus

                    m += 1

        if fusion == 0:
            tabfus = [fusion, x, y, subject]
            return tabfus

    def select_besthit(self, tab):
        """Selects the best bitscoring alignment in a stack."""
        bitscores = []
        bestbitscore = 0.00
        for s, _element in enumerate(tab):
            bitscores.append(float(tab[s].bitscore))
        bestbitscore = max(bitscores)
        i = bitscores.index(bestbitscore)

        return tab[i]

    def get_uncovered_zones(self):
        """Given final retained alignments:
        Get coordinates of retained alignments and
        cut where there is no alignments (minBplength applies)
        It fullfill the two arrays of contig class : Contig.departs
        and Contig.arrivees.
        """
        new_selectedSubjects = \
            self.sort_subjects_on_starting(self.selectedSubjects)
        if len(new_selectedSubjects) == 1:
            s = new_selectedSubjects[0]
            if int(s.qstart) > self.minBpLength:
                splitdep = 1
                splitarr = int(s.qstart) - 1
                self.departs.append(splitdep)
                self.arrivees.append(splitarr)
            if int(s.qend) + self.minBpLength < int(s.qlen):
                splitdep = int(s.qend)+1
                splitarr = int(s.qlen)
                self.departs.append(splitdep)
                self.arrivees.append(splitarr)
        else:

            for s1 in (range(len(new_selectedSubjects) - 1)):
                s = new_selectedSubjects[s1]

                if s1 == 0 and int(s.qstart) > self.minBpLength:
                    splitdep = 1
                    splitarr = int(s.qstart)-1
                    self.departs.append(splitdep)
                    self.arrivees.append(splitarr)
                    if int(new_selectedSubjects[s1].qend) +\
                       self.minBpLength <\
                       int(new_selectedSubjects[s1+1].qstart):
                        splitarr = int(new_selectedSubjects[s1+1].qstart)-1
                        splitdep = int(s.qend)+1
                        self.departs.append(splitdep)
                        self.arrivees.append(splitarr)
                else:
                    if int(new_selectedSubjects[s1].qend) +\
                       self.minBpLength <\
                       int(new_selectedSubjects[s1+1].qstart):
                        splitarr = int(new_selectedSubjects[s1+1].qstart)-1
                        splitdep = int(s.qend)+1
                        self.departs.append(splitdep)
                        self.arrivees.append(splitarr)

            if len(new_selectedSubjects) > 1:
                s = new_selectedSubjects.pop()
                if int(s.qend)+self.minBpLength < int(s.qlen):
                    splitdep = int(s.qend) + 1
                    splitarr = int(s.qlen)
                    self.departs.append(splitdep)
                    self.arrivees.append(splitarr)

            print("get_uncovered_zones method completed with :\n d=",
                  self.departs,
                  "\n a=",
                  self.arrivees)

    def recycleseq(self, sequence):
        """Given a sequence split, according to coordinates"""
        pool = []
        for i, _element in enumerate(self.departs):
            depart = int(self.departs[i]) - 1
            arrivee = int(self.arrivees[i]) - 1
            idseq = ">" + self.id + \
                    "_sect:" + str(self.departs[i]) + \
                    "_" + str(self.arrivees[i])
            section = sequence[depart:arrivee]
            pool.append([idseq, section])

        return pool
