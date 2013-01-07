#!/usr/bin/python

from setuptools import setup
from setuptools.extension import Extension

exts = [Extension('cpairwise2',
                ['pdbtools/cpairwise2module.c'],
                )]
setup(name='pdbtools',
      version='0.1',
      description='PDBTools : General Purpose PDB Manipulation Utility',
      author='Suk Namgoong',
      ext_modules = exts,
      author_email='suk.namgoong@gmail.com',
      url='https://github.com/madscientist01/pdbtools',
      platforms='any',
      install_requires=['biopython'],
      packages=['pdbtools'],
      )
