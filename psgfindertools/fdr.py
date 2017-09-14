#!/usr/bin/python

'''Functions to estimate the number of independant tests out of a set of
overlaping windows'''

from psgfindertools import PSGParam
from psgfindertools.dna import aadiff

# FUNCTIONS - FDR
def rho(windows):
    '''returns the rate of overlaping (rho)'''
    
    if len(windows) == 0: return 1
    windows.sort()
    coordinates, n_coordinates = set(), 0
    for w in windows: 
        coordinates.update(set(w))
        n_coordinates += len(w)
    coverage = len(coordinates)
    return n_coordinates / float(coverage)

def fdr(estimations, alpha, stats=None):
    '''calculates the degree of independance; returns a p-value threshold'''
        
    independent_tests, subset = 0, []
    fname = estimations[0]['file name']
    tr_al = estimations[0]['nucleotides'].get_translated()
    for x in estimations:
    
        # Check overlapping windows in within the same alignment
        if fname == x['file name']:
            a, b = x['window'].split('-')
            a, b = int(a)/3, int(b)/3
            subset.append([ i+a for i in aadiff(tr_al.get_slice(a, b)) ])
        
        # Sum the estimated number of independant tests and initialize
        # for the next alignment.
        else:
            independent_tests += (len(subset) / rho(subset))
            fname = x['file name']
            tr_al = x['nucleotides'].get_translated()
            subset = [[ c for c in aadiff(tr_al) ]]
    
    # Add the last estimation.
    independent_tests += (len(subset) / rho(subset))
    
    # Calculate the corrected p-value threshold
    threshold = (alpha/independent_tests) if independent_tests > 0 else alpha
    if stats is not None: stats['RIT'] = independent_tests
    return threshold
