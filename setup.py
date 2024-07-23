import os
import sys
from setuptools import setup

setup(name='lsaBGC-Pan',
      version='1.0.0',
      description='Suite for comparative genomic, population genetic and evolutionary analysis of micro-evolutionary novelty in BGCs all in the context of a single species or genus.',
      url='http://github.com/Kalan-Lab/lsaBGC-Pan/',
      author='Rauf Salamzade',
      author_email='salamzader@gmail.com',
      license='BSD-3',
      packages=['lsaBGC'],
      scripts=['docker/LSABGC',
               'scripts/processAndReformatUserProvidedGenbanks.py',
               'scripts/genbankToProkkaGFF.py',
               'scripts/popstrat',
               'scripts/phylogenate',
               'scripts/GSeeF',
               'bin/lsaBGC-Cluster',
               'bin/lsaBGC-See',
               'bin/lsaBGC-ComprehenSeeIve',
               'bin/lsaBGC-MIBiGMapper',
               'bin/lsaBGC-Reconcile',
               'bin/lsaBGC-Sociate',
               'workflows/lsaBGC-Pan'],
      zip_safe=False)

