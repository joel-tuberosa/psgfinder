#!/usr/bin/python

import unittest
import subprocess
import os
import pickle
import gzip

class PSGfinderTest(unittest.TestCase):
    
    # import numpy
    def test_import_numpy(self):
        '''Try to import the numpy module (needed for the fisher module)'''
        try: import numpy
        except ImportError:
            self.assertTrue(False)
            
    # import fisher module
    def test_import_numpy(self):
        '''Try to import the fisher module'''
        try: from fisher import pvalue
        except ImportError:
            self.assertTrue(False)    
    
    # try to run yn00
    def test_yn00(self):
        '''Try to run yn00'''
        try: r = subprocess.call(['yn00'], 
                stdout=open(os.devnull, 'w'),
                stderr=subprocess.STDOUT)
        except OSError as e:
            self.asserFalse(e.errno == os.errno.ENOENT)
        self.assertEqual(r, 127)
    
    # search for windows and calculate dN/dS
    def test_parse(self):
        '''Define windows and estimate rates of substitutions'''
        
        from psgfindertools.dna import handle_fasta
        from psgfindertools.dwin import parse
        with open("cd4-human_chimp.fas") as f:
            a = handle_fasta(f, aligned=True)
        estimations = parse(a, fname="cd4-human_chimp.fas")
        estimations[0]["nucleotides"] = a
        with gzip.open("cd4-estimations") as f:
            original_estimations = pickle.load(f)
        self.assertEqual(estimations, original_estimations)
    
    # Find windows with a dN/dS > 1 and compute significance
    def test_dnds_stat(self):
        '''Find windows with a dN/dS > 1 and compute Fisher's p-value'''
        
        from psgfindertools.dwin import dnds_stat
        with gzip.open("cd4-estimations") as f:
            original_estimations = pickle.load(f)
        filtered_estimations = dnds_stat(original_estimations)
        with gzip.open("cd4-filtered_estimations") as f:
            original_filtered_estimations = pickle.load(f)
        self.assertEqual(filtered_estimations, original_filtered_estimations)
    
    # Calculate test overlap
    def test_fdr(self):
        '''Estimate windows overlapping level'''
        
        from psgfindertools.fdr import fdr
        with gzip.open("cd4-filtered_estimations") as f:
            original_filtered_estimations = pickle.load(f)
        self.assertEqual(
            fdr(original_filtered_estimations, 1), 0.41304347826086957)

if __name__ == "__main__":
    unittest.main()
