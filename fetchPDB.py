#!/usr/bin/python
# -*- coding: utf-8 -*-

import urllib2
import os
import argparse

class PDBFetch(object):

    '''
    
    '''
    def __init__(self, **kwargs):
        self.pdblist = kwargs.get('pdblist')
        self.biologicalUnit = kwargs.get('biologicalUnit', True)
        self.verbose = kwargs.get('verbose', True)
        self.path = kwargs.get('path','')
        if len(self.path) > 0 and not os.path.isdir(self.path):
            os.mkdir(self.path)
        if len(self.path) > 0 and self.path[len(self.path) - 1] != '/':
            self.path = self.path + '/'

    def download(self):

        if self.biologicalUnit:
            pdburl = \
                'ftp://ftp.wwpdb.org/pub/pdb/data/biounit/coordinates/all/{0}.pdb.gz'
        else:
            pdburl = 'http://www.rcsb.org/pdb/files/{0}.pdb.gz '

        for pdb in self.pdblist:
            n = 1
            filename = pdb + '.pdb'
            if os.path.exists(pdb + '.pdb'):
                filename = pdb + '.' + str(n) + '.pdb'
                while os.path.exists(filename):
                    n += 1
                    filename = pdb + '.' + str(n) + '.pdb'

            if self.verbose:
                print 'download {0}'.format(pdb)
            try:
                f = urllib2.urlopen(pdburl.format(pdb.lower()))
                data = f.read()
                with open(filename + '.gz', 'wb') as code:
                    code.write(data)
                os.system('gunzip {0}.gz'.format(filename))
            except urllib2.URLError:
                print '{0} : Error in pdb code. Skipped.'.format(pdb)


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-l',
        '--list',
        dest='list',
        default=[],
        nargs='+',
        help='pdblist to search',
        )
    parser.add_argument(
        '-b',
        '--biological_unit',
        dest='biologicalUnit',
        action='store_true',
        default=False,
        help='download biological unit of PDB',
        )
    parser.add_argument(
        '-v',
        '--verbose',
        dest='verbose',
        action='store_true',
        default=False,
        help='verbose',
        )

    parser.add_argument('-p', '--path', dest='path', default='',
                        help='generate subdirectory')
    results = parser.parse_args()
    fetch = PDBFetch(pdblist=results.list, biologicalUnit=results.biologicalUnit, verbose=results.verbose, path=results.path)
    fetch.download()

