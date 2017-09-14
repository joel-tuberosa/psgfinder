#!/usr/bin/python

'''Functions to map putative windows of divergent positive selection and
estimate their dN/dS with yn00'''

import re
import sys
import os
import subprocess
import tempfile
from fisher import pvalue
from functools import partial
from StringIO import StringIO

from psgfindertools import PSGParam
from psgfindertools.misc import recursive_rm, is_number
from psgfindertools.dna import phylip, map_gap, aadiff

# FUNCTION - YN00
def read_yn00output(f):
    c = 0
    for line in f:
        if c == 2: 
            line = line.replace("-", "").replace("+", "").split()
            yield [ x if is_number(x) else 'nan' for x in line ]
            c = 0
        elif c == 1:
            c += 1
        elif line.startswith("seq. seq."):
            c = 1

def yn00_pairwise(*pairwise_als, **kwargs):
    '''run yn00 in a subprocess, pairwise alignment only'''
    
    # defaults
    yn00ctl = '''      seqfile = al.nuc * sequence data file name
      outfile = al.nuc.yn           * main result file
      verbose = 0  * 1: detailed output (list sequences), 0: concise output

        icode = 0  * 0:universal code; 1:mammalian mt; 2-10:see below

    weighting = 0  * weighting pathways between codons (0/1)?
   commonf3x4 = 0  * use one set of codon freqs for all pairs (0/1)? 
        ndata = 1


* Genetic codes: 0:universal, 1:mammalian mt., 2:yeast mt., 3:mold mt.,
* 4: invertebrate mt., 5: ciliate nuclear, 6: echinoderm mt., 
* 7: euplotid mt., 8: alternative yeast nu. 9: ascidian mt., 
* 10: blepharisma nu.
* These codes correspond to transl_table 1 to 11 of GENEBANK.
'''

    ndata = len(pairwise_als)
    main_output = []
    flush = True
    current_dir = os.getcwd()
    estimf = None
    base_header = ['sequences', 'size', '%gap', 'S', 'N', 't', 'kappa',
        'omega', 'dN', 'dN SE', 'dS', 'dS SE']
    estimf_values = list(base_header) ### discociated instances!
    estimf_values_to_add = []
    wdir = tempfile.mkdtemp(prefix='yn00_')
    extra_val = {}
    
    # parse function's keyword arguments
    for k in kwargs:
        
        # make a tempdir or use an existing path if wdir is specified
        if k == 'wdir' and kwargs['wdir'] is not None: 
            wdir = kwargs['wdir']
            flush = False
    
        # values are returned by the function
        elif k == 'values':
            values = kwargs['values']
        
        # store estimf variable, which tells the function to directly
        # write estimations in the temp file instead of storing a huge
        # amount of data in the RAM. estimf is an open file object
        elif k == 'estimf':
            estimf = kwargs['estimf']
        
        # write those values to estimf, otherwise write everything
        elif k == 'estimf_values':
            estimf_values = kwargs['estimf_values']
        
        # extra values can be added to the output via kwargs
        else:
            extra_val[k] = kwargs[k]
            estimf_values_to_add.append(k)
            
        # they are added to estimf_values if they are not already
        # stated
        for k in estimf_values_to_add:
            if k not in estimf_values: estimf_values.append(k)
        
    # cd in it
    if os.path.isdir(wdir): os.chdir(wdir)
    else: raise IOError('{} does not exist or is not a directory'.format(wdir))
    
    try:
        sizes = []
        with open('al.nuc', 'w') as pamlinput:
            als_count = 0
            for pairwise_al in pairwise_als:
                pamlinput.write(phylip(pairwise_al, l=-1, s=-1) + '\n')
                als_count += 1
                assert pairwise_al.length != 0, pairwise_al.get_names()[0]
                sizes.append(pairwise_al.length)
            
        with open('yn00.ctl', 'w') as ctlf:
            ctlf.write(yn00ctl.replace('ndata = 1', 'ndata = {}'.format(
                als_count)))
        
        # run yn00 and raise error if something's wrong
        returncode = subprocess.call(['yn00'], stdout=open(os.devnull, 'w'))
        
        # if an error occured, it will be prompt in stderr by yn00
        if returncode != 0: raise ValueError()
        with open('al.nuc.yn', 'r') as yn00output:
            i = 0
            for m in read_yn00output(yn00output):
                
                # match the window size and the sequence names
                size = sizes[i]
                sequences_names = '||'.join(pairwise_als[i].get_names())
                
                # 11 items must be found for each yn00 estimations
                if len(m) != 11:
                    sys.stderr.write(
                        'Warning: cannot retrieve all estimations from ' +
                        '{}||{}\n'.format(*pairwise_al.names[:2]))
                    i += 1
                    continue                
                
                m = m[2:]
                pcgap = map_gap(pairwise_als[i]).count('-')/float(size)
                estim = dict(zip(
                    base_header, [sequences_names, size, pcgap] + m))
                
                # extra values
                for k in extra_val:
                    estim[k] = extra_val[k][i] if           \
                        isinstance(extra_val[k], list) else \
                        extra_val[k]
                
                # write everything in estimf if defined
                if estimf is not None:
                    estimf.write('\t'.join(
                        ( str(estim[key]) for key in estimf_values )) + 
                        '\n')

                # the kwarg "values" gives the values to be returned
                main_output.append(dict(
                    ( (k, estim[k]) 
                      for k in values + sorted(extra_val.keys()) )))
                    
                # increment i
                i += 1
    finally:
        os.chdir(current_dir)
        if flush: recursive_rm(wdir)
    return main_output

# FUNCTIONS - WINDOWS PROCESSING
def map_win(al, msize=3, mmut=3):
    '''returns every possible windows given a list of coordinates'''

    if mmut > msize: raise ValueError('msize must be greater or equal to mmut')
    win_map = []
    c = aadiff(al)
    for x in range(len(c)-(mmut-1)):
        for y in range(mmut-1+x, len(c)):
            if (c[y]+1)-c[x] >= msize: win_map.append((c[x], c[y]+1))
            
    # coordinates are for amino acids
    return win_map

def parse(al, values=None, fname='-', wdir=None, estimf=None, config=None,
          stats=None):
    '''Searches for windows with putative dN/dS > 1'''
    
    # runs with default parameters
    if config is None: config = PSGParam()
    
    if values is None:
        values = ['sequences', 'size', 'N', 'S', '%gap', 't', 'kappa', 'dN', 
            'dN SE', 'dS', 'dS SE', 'omega']
    
    # everything is written in estimf as long as this variable points
    # toward an open file object.
    
    # whole gene estimations
    whole_gene_estimations = yn00_pairwise(al, 
        values=values+['dN', 'dS'], 
        wdir=wdir, 
        estimf=estimf, 
        estimf_values=config.estimf_header,
        **{'window': '0-{:d}'.format(al.length), 
           'whole gene': '1', 
           'file name': fname})[0]
    
    # add whole gene's dN and dS to the stats
    if stats is not None:
        stats['dN'].append(float(whole_gene_estimations['dN']))
        stats['dS'].append(float(whole_gene_estimations['dS']))

    # get window coordinates (in nucleotides)
    windows_coordinates = [ (x*3, y*3) for x, y in 
       map_win(al, msize=config.msize, mmut=config.mmut) ]
    
    # maximum by 5000 windows to avoid yn00 crash!
    windows_estimations = []
    for i in xrange(0, len(windows_coordinates), 5000):
        step = i
        windows = [ al.get_slice(start=x, stop=y)
            for x, y in windows_coordinates[i:i+5000] ]
        
        # substitution rates estimations for each window
        windows_estimations += yn00_pairwise(*windows, 
            values=values+['dN', 'dS'], 
            wdir=wdir, 
            estimf=estimf,
            estimf_values=config.estimf_header,
            **{'window': 
                map(lambda xy: '{:d}-{:d}'.format(*xy), windows_coordinates), 
               'whole gene': '0', 
               'file name': fname})
        
    # add the number of analyzed windows to the stats
    if stats is not None:
        stats['n_win'] += len(windows_coordinates)
    
    return [whole_gene_estimations] + windows_estimations
    
def dnds_stat(estimations):
    '''return estimations of windows with dN/dS > 1'''
    
    filtered_estimations = []
    fname = estimations[0]['file name']
    genedS = float(estimations[0]['dS'])
    for i in range(len(estimations)):
        if 'nan' in estimations[i].values():
            continue
        name = estimations[i]['file name']
        if name != fname:
            fname = estimations[i]['file name']
            genedS = float(estimations[i]['dS'])
        if genedS != 0:
            estimations[i]['dN/dS(whole gene)'] = float(
                estimations[i]['dN'])/genedS
        else: continue
        
        # process numbers with fisher module
        if estimations[i]['whole gene'] == '1' or \
           estimations[i]['dN/dS(whole gene)'] > 1:
            n = round(float(estimations[i]['dN']) * float(estimations[i]['N']))
            N = round(float(estimations[i]['N'])) - n
            s = round(float(estimations[i]['dS']) * float(estimations[i]['S']))
            S = round(float(estimations[i]['S'])) - s
            mat = [[n,N],[s,S]]
            p = pvalue(n,N,s,S)
            estimations[i]['p-value'] = p.two_tail
            filtered_estimations.append(estimations[i])
    return filtered_estimations