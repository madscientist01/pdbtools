#!/usr/bin/python
# -*- coding: utf-8 -*-

#
# Wrapper of Superpose (CCP4)
# Written by MadScientist (http://madscientist.wordpress.com)
# Require pdbextract.py and Biopython
#

import sys
import subprocess
import re
import os
import glob
import argparse
import jinja2
from pdbextract import *


class Aligned:

    #
    # Storage class for the structural alignment in Superpose output file
    #

    def __init__(self, **kwargs):
        self.pdb = kwargs.get('pdb')
        self.chain = kwargs.get('chain')
        self.number = kwargs.get('number')
        self.aa = kwargs.get('aa')
        self.distance = kwargs.get('distance')
        self.aligned = kwargs.get('aligned')
        self.secondary = kwargs.get('secondary')
        if self.chain and self.number and self.aa:
            self.id = self.chain + '_' + str(self.number) + '_' \
                + self.aa
        else:
            self.id = None
        self.alignid = None
        self.gap = kwargs.get('gap')
        self.outfile = ''


class Superpose:

    def __init__(self, **kwargs):
        self.queryPDB = kwargs.get('queryPDB')
        self.subjectPDB = kwargs.get('subjectPDB')
        self.superposedPDB = kwargs.get('superposedPDB')
        self.outfile = ''
        self.rmsd = 0
        self.aligned = 0
        self.qscore = 0
        self.queryAlign = []
        self.subjectAlign = []
        self.querySeq = ''
        self.subjectSeq = ''
        self.upper = 0
        self.lower = 0
        self.querySecondary = ''
        self.alignSecondary = ''
        self.alignDistance = []

    def fastaSplit(self, fasta, width):

        #
        # Split sequence in fasta string into predefined width of string lists
        #

        cursor = 0
        end = 0
        buffer = []
        while cursor < len(fasta):
            if len(fasta) - cursor > width:
                end += width
            else:
                end = len(fasta)
            buffer.append(fasta[cursor:end] + '\n')
            cursor += width
        return buffer

    def JSONoutput(self):

        [header, chains] = PDBParse(self.queryPDB, None, None)
        querypdbid = extractHeader(header, '^DBREF  (\S{4}) ', False)
        if querypdbid:
            queryHeader = querypdbid
        else:
            queryHeader = self.queryPDB[:len(self.queryPDB) - 4]
        [header, chains] = PDBParse(self.subjectPDB, None, None)
        subjectpdbid = extractHeader(header, '^DBREF  (\S{4}) ', False)
        if subjectpdbid:
            subjectHeader = subjectpdbid
        else:
            subjectHeader = self.subjectPDB[:len(self.subjectPDB) - 4]
        querySeq = self.querySeq
        subjectSeq = self.subjectSeq
        querySecondary = self.querySecondary
        subjectSecondary = self.subjectSecondary
        distance = self.distanceOutput()
        queryStart = self.queryAlign[self.upper].number
        subjectStart = self.subjectAlign[self.upper].number
        jsons = \
            """
queryid:"{0}",
query:"{1}", 
subjectid:"{2}",
subject:"{3}",
distance:[{4}],
querysecondary:"{5}",
subjectsecondary:"{6}",
querystart:{7},
subjectstart:{8}
		"""
        output = '{' + jsons.format(
            queryHeader,
            querySeq,
            subjectHeader,
            subjectSeq,
            distance,
            querySecondary,
            subjectSecondary,
            queryStart,
            subjectStart,
            ) + '}'
        return output

    def FASTAoutput(self):

        #
        # Generate structural alignment as FASTA formats
        #

        fasta = []
        if len(self.queryAlign) > 0:

            queryHeader = self.queryAlign[self.upper].chain + ':' \
                + str(self.queryAlign[self.upper].number) + '-' \
                + self.queryAlign[self.lower].chain + ':' \
                + str(self.queryAlign[self.lower].number)

            fasta.append('>{0}:{1}\n'.format(self.queryPDB[:len(self.queryPDB)
                         - 4], queryHeader))
            fasta = fasta + self.fastaSplit(self.querySeq, 60)
            subjectHeader = self.subjectAlign[self.upper].chain + ':' \
                + str(self.subjectAlign[self.upper].number) + '-' \
                + self.subjectAlign[self.lower].chain + ':' \
                + str(self.subjectAlign[self.lower].number)

            fasta.append('>{0}:{1}\n'.format(self.subjectPDB[:len(self.subjectPDB)
                         - 4], subjectHeader))
            fasta = fasta + self.fastaSplit(self.subjectSeq, 60)
        return fasta

    def distanceOutput(self):

        distance = ','.join(self.alignedDistance)
        return distance

    def generateAlignment(self):

        oneLetter = {
            'VAL': 'V',
            'ILE': 'I',
            'LEU': 'L',
            'GLU': 'E',
            'GLN': 'Q',
            'ASP': 'D',
            'ASN': 'N',
            'HIS': 'H',
            'TRP': 'W',
            'PHE': 'F',
            'TYR': 'Y',
            'ARG': 'R',
            'LYS': 'K',
            'SER': 'S',
            'THR': 'T',
            'MET': 'M',
            'ALA': 'A',
            'GLY': 'G',
            'PRO': 'P',
            'CYS': 'C',
            'MSE': 'M',
            '-': '-',
            }
        sequence = ''
        aligned = ''
        distance = []
        for i in range(len(self.queryAlign)):
            residue = self.queryAlign[i]

            if residue.aa in oneLetter:
                sequence = sequence + oneLetter[residue.aa]
            else:
                sequence = sequence + 'X'
            if residue.alignid.aa in oneLetter:
                aligned = aligned + oneLetter[residue.alignid.aa]
            else:
                aligned = aligned + 'X'

        uTrim = 0
        lTrim = len(sequence)
        for i in range(len(sequence)):
            if self.queryAlign[i].distance:
                uTrim = i
                break

        for i in reversed(range(len(sequence))):
            if self.queryAlign[i].distance:
                lTrim = i
                break

        self.querySeq = sequence[uTrim:lTrim + 1]
        self.subjectSeq = aligned[uTrim:lTrim + 1]
        self.upper = uTrim
        self.lower = lTrim
        qSecondary = ''
        sSecondary = ''
        for i in range(uTrim, lTrim + 1):
            residue = self.queryAlign[i]
            if residue.distance:
                distance.append(str(residue.distance))
            else:
                distance.append('99.9')
            if residue.secondary:
                qSecondary = qSecondary + residue.secondary
            else:
                qSecondary = qSecondary + ' '

            if residue.alignid.secondary:
                sSecondary = sSecondary + residue.alignid.secondary
            else:
                sSecondary = sSecondary + ' '

        self.alignedDistance = distance
        self.querySecondary = qSecondary
        self.subjectSecondary = sSecondary
        return

    def extractSubjectPDB(self):

        if os.path.exists(self.subjectPDB) and len(self.subjectAlign) \
            > 0:
            startchain = self.subjectAlign[self.upper].chain
            start = self.subjectAlign[self.upper].number
            endchain = self.subjectAlign[self.lower].chain
            end = self.subjectAlign[self.lower].number
            [header, chains] = PDBParse(self.subjectPDB, None, None)

            if startchain == endchain:
                extractRegion = '{0}:{1}_{2}'.format(startchain, start,
                        end)
            else:
                extractRegion = '{0}:{1}_{2} {3}:{4}_{5}'.format(
                    startchain,
                    start,
                    9999,
                    endchain,
                    1,
                    end,
                    )
            pdbExtract = PDBExtract(extractregion=extractRegion,
                                    header=True,overwrite=True)
            newFileName = pdbExtract.extractRegions(header, chains,
                    self.subjectPDB)
            return newFileName
        else:
            print 'problem'
            return filename

    def extractSuperposedPDB(self):

        if os.path.exists(self.superposedPDB) and len(self.queryAlign) \
            > 0:
            startchain = self.queryAlign[self.upper].chain
            start = self.queryAlign[self.upper].number
            endchain = self.queryAlign[self.lower].chain
            end = self.queryAlign[self.lower].number
            [header, chains] = PDBParse(self.superposedPDB, None, None)

            if startchain == endchain:
                extractRegion = '{0}:{1}_{2}'.format(startchain, start,
                        end)
            else:
                extractRegion = '{0}:{1}_{2} {3}:{4}_{5}'.format(
                    startchain,
                    start,
                    9999,
                    endchain,
                    1,
                    end,
                    )
            pdbExtract = PDBExtract(extractregion=extractRegion,
                                    header=True)
            newFileName = pdbExtract.extractRegions(header, chains,
                    self.superposedPDB)
            return newFileName
        else:
            print 'problem'
            return filename

    def run(self):

        #
        # Run superpose with two PDB file (self.queryPDB and self.subjectPDB)
        #

        if os.path.exists(self.queryPDB) \
            and os.path.exists(self.subjectPDB):
            p = subprocess.Popen(['superpose', self.queryPDB,
                                 self.subjectPDB, self.superposedPDB],
                                 stdout=subprocess.PIPE)

        #
        # stdout from superpose will be saved into p_stdout
        #

            p_stdout = p.stdout.read()
            self.outfile = self.superposedPDB[:len(self.superposedPDB)
                - 4] + '.out'
            fw = open(self.outfile, 'w')
            fw.write(p_stdout)
            fw.close()

        #
        # Parsing output file. RMSD, Q-score, # of aligned residues are parsed using regex
        #

            regex = \
                re.compile(" at RMSD =\s+([\d.]+),\s+Q=\s+([\d.]+)\s+and alignment length\s+([\d.]+)"
                           )
            match = regex.search(p_stdout)

            if match:
                self.rmsd = float(match.group(1))
                self.qscore = float(match.group(2))
                self.aligned = float(match.group(3))
                success = True
            else:
                success = False
                return success

        #
        # Parsing sequence alignment
        # ............

            data = p_stdout.split('\n')
            capture = False
            regex = re.compile('^\|(.*)\|(.*)\|(.*)\|')

            for line in data:
                if line == '|-------------+------------+-------------|':
                    capture = True

                if line \
                    == "`-------------\'------------\'-------------\'":
                    capture = False

                if capture and line \
                    != '|-------------+------------+-------------|':
                    mat = regex.match(line)
                    if mat:
                        query = mat.group(1)
                        distance = mat.group(2)
                        subject = mat.group(3)
                    querySecondary = query[0:1].strip()
                    queryChain = query[3:4].strip()
                    queryAmino = query[5:8].strip()
                    queryAminoNum = query[8:12].strip()
                    queryDistance = distance[4:8].strip()
                    subjectSecondary = subject[0:1].strip()
                    subjectChain = subject[3:4].strip()
                    subjectAmino = subject[5:9].strip()
                    subjectAminoNum = subject[8:12].strip()

                    # print querySecondary, queryChain, queryAmino, queryAminoNum, queryDistance
                    # print subjectSecondary, subjectChain, subjectAmino, subjectAminoNum

                    if len(querySecondary) + len(queryChain) \
                        + len(queryAmino) + len(queryAminoNum) > 0:
                        queryresidue = Aligned(
                            pdb=self.queryPDB,
                            chain=queryChain,
                            aa=queryAmino,
                            number=int(queryAminoNum),
                            gap=False,
                            secondary=querySecondary,
                            )
                    else:
                        queryresidue = Aligned(pdb=self.queryPDB, aa='-'
                                , gap=True)

                    if len(subjectSecondary) + len(subjectChain) \
                        + len(subjectAmino) + len(subjectAminoNum) > 0:
                        subjectresidue = Aligned(
                            pdb=self.subjectPDB,
                            chain=subjectChain,
                            aa=subjectAmino,
                            number=int(subjectAminoNum),
                            gap=False,
                            secondary=subjectSecondary,
                            )
                    else:
                        subjectresidue = Aligned(pdb=self.subjectPDB,
                                aa='-', gap=True)

                    if len(queryDistance):
                        queryresidue.distance = float(queryDistance)
                        subjectresidue.distance = float(queryDistance)

                    queryresidue.alignid = subjectresidue
                    subjectresidue.alignid = queryresidue
                    self.queryAlign.append(queryresidue)
                    self.subjectAlign.append(subjectresidue)

            self.generateAlignment()
            return success


def htmlout(table, argument, pdb):
    #
    # Generate render.html which display rendered png files
    #
    #
    jinja_environment = \
    jinja2.Environment(loader=jinja2.FileSystemLoader(os.path.dirname(__file__)))
    template = jinja_environment.get_template('template.html')
    filename = pdb[:len(pdb) - 4] + '.html'
    scriptname = pdb[:len(pdb) - 4] + '.js'
            # header for the hits information table
    headers = [
        'Query',
        'PDB ID',
        'Title',
        'RMSD',
        'Q-Score',
        'Aligned Amino Acid',
        'Aligned PDB'
        ]

    [header, chains] = PDBParse(pdb, None, None)
    queryTitle = extractHeader(header, '^TITLE    [ \d](.*)$', True)
    queryPDBId = extractHeader(header, '^DBREF  (\S{4}) ', False)
    scriptdataList = []
    tableData = []
    for (alignpdb, align) in table.items():
        [header, chains] = PDBParse(align.superposedPDB, None, None)
        title = extractHeader(header, '^TITLE    [ \d](.*)$', True)
        pdbid = extractHeader(header, '^DBREF  (\S{4}) ', False)
        canvasHeight = int(100 + (len(align.querySeq) / 40 + 1) * 70)
        singleRecord = {
                'query' : queryPDBId, 
                'queryTitle' : queryTitle,
                'pdbid' : pdbid,
                'title' : title,
                'rmsd' : align.rmsd,
                'qscore' : align.qscore,
                'aligned' : probability,
                'alignedpdb' : alignpdb,
                'canvas' : 'canvas_'+alignpdb[:len(alignpdb)-4],
                'canvasHeight':canvasHeight
            }
        tableData.append(singleRecord)
        data = {
        'pdb': alignpdb[:len(alignpdb) - 4], 
        'json': align.JSONoutput()
        }
        scriptdataList.append(data)

    t = template.render(scriptname=scriptname, tableData= tableData, headers=headers, scriptdata=scriptdataList)
    f = open(hhblits.path + filename, 'w')
    f.write(t)
    f.close()

   

def pymolRendering(table, argument, pdb):

    #
    # Pymol Rendering of Stepwise comparison
    #
    #

    PyMOLPath = '/Applications/MacPyMOL.app/Contents/MacOS/MacPyMol'
    pmlFile = pdb[:len(pdb) - 4] + '.pml'
    pymolFile = open(pmlFile, 'w')

    if argument.view:
        if os.path.exists(argument.view):
            f = open(argument.view)
            view = f.readlines()
            if re.match('^set_view', view[0]):
                for line in view:
                    pymolFile.write(line)
            f.close()

    pymolFile.write('bg_color white\n')
    pymolFile.write('load {0}/{1}\n'.format(os.getcwd(), pdb))
    pymolFile.write('hide all\n')

    for (superposed, align) in table.items():

        if argument.stepwise:
            pymolFile.write('load {0}/{1}\n'.format(os.getcwd(),
                            align.subjectPDB))
            pymolFile.write('load {0}/{1}\n'.format(os.getcwd(),
                            align.superposedPDB))
            pymolFile.write('hide all\n')
            pymolFile.write('show cartoon, {0}\n'.format(align.superposedPDB[:-4]))
            pymolFile.write('show cartoon, {0}\n'.format(align.subjectPDB[:-4]))
            pymolFile.write('util.cbc {0}\n'.format(align.subjectPDB[:-4]))
            if not argument.view:
                pymolFile.write('orient {0}\n'.format(pdb[:-4]))
            pymolFile.write('ray\n'.format(align.superposedPDB[:-4]))
            pymolFile.write('png {0}.png\n'.format(align.superposedPDB[:-4]))
        else:
            pymolFile.write('load {0}/{1}\n'.format(os.getcwd(),
                            align.superposedPDB))

    pymolFile.write('hide all\n')
    pymolFile.write('show cartoon\n')
    pymolFile.close()

    if argument.render:
        print 'Running PyMol...'
        p = subprocess.Popen([PyMOLPath, '-c', pmlFile],
                             stdout=subprocess.PIPE)
        p_stdout = p.stdout.read()
    if argument.stepwise:
        htmlout(table, argument, pdb)


def main(argument):

    pdbList = argument.files

    if len(pdbList) > 0:
        for pdb in pdbList:
            if os.path.exists(pdb):

                subjectList = glob.glob('*.pdb')
                pdbLoadList = []

                if len(subjectList) > 0:
                    table = {}
                    for onePdb in subjectList:
                        if onePdb != pdb:
                            print 'Processing {0}'.format(onePdb)
                            superposedFilename = pdb[:len(pdb) - 4] \
                                + '_' + onePdb
                            sup = Superpose(queryPDB=onePdb,
                                    subjectPDB=pdb,
                                    superposedPDB=superposedFilename)
                            if sup.run():
                                if argument.extract:
                                    superposedFilename = sup.extractSuperposedPDB()
                                    sup.superposedPDB = superposedFilename
                                    subjectFilename = sup.extractSubjectPDB()
                                    sup.subjectPDB = subjectFilename

                                write = True
                                if argument.score and sup.qscore >= float(argument.score):
                                    print '{0} is skipped. Q Score:{1}'.format(sup.queryPDB,sup.qscore)
                                    write = False

                                if argument.rmsd and sup.rmsd >= float(argument.rmsd):
                                    print '{0} is skipped. RMSD:{1}'.format(sup.queryPDB, sup.rmsd)
                                    write = False

                                if argument.align and sup.aligned >= float(argument.align):
                                    print '{0} is skipped. aligned:{1}'.format(sup.queryPDB, sup.aligned)
                                    write = False

                                if write:

                                    table[superposedFilename] = sup
                                    pdbLoadList.append(superposedFilename)
                                    if argument.stalign:
                                        aligned = sup.FASTAoutput()
                                        if aligned:
                                            f = open(superposedFilename[:len(superposedFilename) - 4] 
                                                + '.fasta', 'w')
                                            f.writelines(aligned)
                                            f.close()

                    pymolRendering(table, argument, pdb)
                else:
                    print 'File is not exist!'
                    sys.exit()


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    
    parser.add_argument(
        '-r',
        '--reference_pdb',
        nargs='+',
        dest='files',
        default=[],
        help='Files to process',
        )
    parser.add_argument(
        '-e',
        '--extract_superimpose',
        action='store_true',
        dest='extract',
        default=False,
        help='Extract superposed region',
        )
    parser.add_argument(
        '-d',
        '--rmsd_cutoff',
        action='store',
        dest='rmsd',
        type=float,
        help='Set RMSD cutoff',
        )
    parser.add_argument(
        '-s',
        '--score_cutoff',
        action='store',
        dest='score',
        type=float,
        help='Set Q Score cutoff',
        )
    parser.add_argument(
        '-a',
        '--align_cutoff',
        action='store',
        dest='align',
        type=float,
        help='Set alignment length cutoff',
        )
    parser.add_argument(
        '-t',
        '--stepwise_compare_rendering',
        action='store_true',
        dest='stepwise',
        default=False,
        help='Render 1:1 comparison',
        )
    parser.add_argument(
        '-p',
        '--pymol',
        action='store_true',
        dest='render',
        default=False,
        help='Execute PyMOL',
        )
    parser.add_argument(
        '-n',
        '--align',
        action='store_true',
        dest='stalign',
        default=False,
        help='Save structure based sequence alignment',
        )
    parser.add_argument('-v', '--fixed_view', action='store',
                        dest='view',
                        help='Filename which contains set_view command of PyMOL'
                        )

    results = parser.parse_args()
    main(results)
