#!/usr/bin/python

'''Functions for alignment cleaning.'''

import re
from psgfindertools.dna import map_gap, aadiff, map_gap_coordinates,           \
                               map_al_coordinates
from psgfindertools.dwin import yn00_pairwise

# FUNCTIONS - ALIGNMENTS CLEANING
def pairwise_cleaning(al, n=30, dsmax=1.0, k=2, d=4, q=6):
    '''cleans pairwise alignment regarding gaps, dS and aa differences'''
    
    if len(al) == 2:
        al = gap_cleaning(al, n=n)
        if map_gap(al).count('-') == al.length: return al
        if dsmax > 0: 
            al = ds_cleaning(al, dsmax=dsmax)
            if map_gap(al).count('-') == al.length: return al
        al = aadiff_cleaning(al, k=k, d=d, q=q)
        return al
    else: raise ValueError('need pairwise alignment (only 2 sequences)')
    
def gap_cleaning(al, n=30):
    '''masks uninterrupted regions shorter than n nucleotides'''
    
    if al.__class__.__name__ != 'Alignment':
        raise ValueError('need aligned data')
    gap_map = map_gap(al)
    if '-' not in gap_map:
        if al.length < n: al.mask(0, al.length)
        return al
    pattern = re.compile('^X{l}1,{n}{r}-|-X{l}1,{n}{r}-|-X{l}1,{n}{r}$'.format(
        l='{', n=n, r='}'))
    matches = pattern.findall(gap_map)
    if matches:
        if matches[0][0] == 'X':
            gap_map = re.sub('^'+matches[0], '-'*len(matches[0]), gap_map)
            del matches[0]
        if matches and matches[-1][-1] == 'X':
            gap_map = re.sub(matches[-1]+'$', '-'*len(matches[-1]), gap_map)
            del matches[-1]
        for m in matches:
            gap_map = gap_map.replace(m, '-'*len(m))
        assert len(gap_map) == al.length ### debug...
        
        for x, y in map_gap_coordinates(gap_map):
            al.mask(x, y)
    return al
    
def ds_cleaning(al, dsmax=1.0):
    '''
    masks uninterrupted regions whose ds is lower than a threshold uses
    yn00 for pairwise alignment only
    '''
    
    if len(al) != 2: raise ValueError('need a pairwise alignment')
    al_map = map_al_coordinates(map_gap(al))
    if not al_map: return al
    results = yn00_pairwise(*[ al.get_slice(x, y) for x, y in al_map ],
        values=['dS'])
    for i in range(len(results)):
        ds = results[i]['dS']
        if ds == 'nan': continue
        if float(ds) > dsmax: al.mask(*tuple(al_map[i]))
    return al
    
def aadiff_cleaning(al, coordinates=None, k=2, d=4, q=6):
    '''masks regions with high concentration of amino acid differences'''
    
    # translate alignment
    tr_al = al.get_translated()
    gap_map = map_gap(tr_al)
    
        
    # Find clusters
    if coordinates is None: coordinates = aadiff(tr_al)
    if k <= 1: clusters = [ [c] for c in coordinates ]    
    else:
        clusters, cluster = [], []
        for i in range(len(coordinates)-1):
            i_ngaps = gap_map[coordinates[i]:coordinates[i+1]].count('-')
            if coordinates[i+1] == coordinates[i] + 1 + i_ngaps:
                if coordinates[i] not in cluster: 
                    cluster.append(coordinates[i])
                cluster.append(coordinates[i+1])
            elif cluster:
                clusters.append(cluster)
                cluster = []
        if cluster: clusters.append(cluster)
    
    # Discard clusters if they are smaller than k
    clusters = [ cluster for cluster in clusters if len(cluster) >= k ]
    if not clusters: return al
    
    # Concatenate clusters if they are less distant than d
    if len(clusters) > 1:
        i, N = 0, len(clusters)
        while i <= N - 2:
            i_ngaps = gap_map[clusters[i][-1]:clusters[i+1][-1]].count('-')
            if clusters[i+1][0] - clusters[i][-1] - i_ngaps <= d:
                clusters[i+1] = clusters[i] + clusters[i+1]
                del clusters[i]
                N -= 1
            else: i += 1
    
    # Add the amino acid changes that are within the clusters but not
    # accounted in the previous procedures
    for i in range(len(clusters)):
        clusters[i] += [ c for c in coordinates 
            if (c > clusters[i][0] 
            and c < clusters[i][-1]) 
            and c not in clusters[i]]

    # Discard clusters if they have less than q amino acid changes
    clusters = [ cluster for cluster in clusters if len(cluster) > q ]
    if not clusters: return al
    
    # Mask bad alignment regions with gaps
    for cluster in clusters:
        al.mask(cluster[0]*3, (cluster[-1]+1)*3)
    return al
