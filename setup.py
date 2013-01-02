#!/usr/bin/python

from setuptools import setup

setup(name='pdbtools',
      version='0.1',
      description='PDBTools : General Purpose PDB Manipulation Utility',
      author='Suk Namgoong',
      author_email='suk.namgoong@gmail.com',
      url='https://github.com/madscientist01/pdbtools',
      platforms='any',
      install_requires=['biopython'],
      packages=['pdbtools'],
      )
