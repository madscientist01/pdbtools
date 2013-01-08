#!/usr/bin/python
# -*- coding: utf-8 -*-
# fetchPDB
# Script and Class for the easy pdb fetch from RCSB PDB
#
#
# The MIT License
#
# Copyright (c) 2012 Suk Namgoong (suk.namgoong@gmail.com)
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.
import urllib2
import os
import argparse

class fetchPDB(object):

    '''
    PDB Fetch utility and class
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
    fetch = fetchPDB(pdblist=results.list, biologicalUnit=results.biologicalUnit, verbose=results.verbose, path=results.path)
    fetch.download()

