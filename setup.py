#!/usr/bin/python

from distutils.core import setup

def readme():
    with open('README.md') as f: return f.read()

setup(
    name='PSGfinder',
    description='Finding signals of positive selection in pairwise alignments',
    url='https://genev.unige.ch/research/laboratory/Juan-Montoya',
    author='Joel Tuberosa',
    author_email='joel.tuberosa@unige.ch',
    license='GNU',
    scripts=['script/psgfinder.py'],
    packages=['psgfindertools'],
    platforms=['any'],
    version='2.5',
    long_description=readme(),
      )
