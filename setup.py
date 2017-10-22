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
    version='2.6.b1',
    long_description=readme(),
    classifiers=[
        'Development Status :: 4 - Beta',
        'Environment :: Console',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Natural Language :: English',
        'Operating System :: MacOS :: MacOS X',
        'Operating System :: Microsoft :: Windows',
        'Operating System :: POSIX :: Linux',
        'Programming Language :: Python :: 2.7',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        ]
      )
