#!/usr/bin/python

'''Miscellaneous functions used in psgfinder modules'''

# EXTERNAL MODULES
import os, shutil, tempfile
from math import floor, exp, factorial

# CLASS - TEMP DIR
class TemporaryDirectory(object):
    '''context manager for tempdir - usable with "with" statement'''
    
    def __enter__(self):
        self.name = tempfile.mkdtemp()
        
        # record origin directory
        self.origin = os.getcwd()
        
        # cd in the temp directory
        os.chdir(self.name)
        return self.name

    def __exit__(self, exc_type, exc_value, traceback):
        
        # delete temp directory
        shutil.rmtree(self.name)
        
        # cd in the origin directory
        os.chdir(self.origin)
        
# FUNCTIONS - MISCELLANEOUS
def readrange(s):
    s = s.strip()
    if s.find('-') == -1: return slice(int(s)-1, int(s))
    else:
        if s.count('-') == 1:
            x, y = s.split('-')
            if s[0] == '-': return slice(int(y))
            elif s[-1] == '-': return slice(int(x), None)
            else: return slice(int(x), int(y))
        else: raise ValueError('range must contain at most one "-"')

def is_number(x):
    '''check if a variable can be converted to a float'''
    try: x = float(x)
    except ValueError: return False
    return True
        
def plot_cat(_f,_c):
    '''adds values to a frequency distribution; (for positive integer)'''
    if _c < 0: raise ValueError('positive integer needed, negative found')
    if _c+1 > len(_f): _f.extend([0]*(_c+1-len(_f)))
    _f[_c] = _f[_c] + 1
    return _f

def parsing_order(l):
    '''returns a reordered list for efficient parsing'''
    dividers = [ i+1 for i in range(l) if l % (i+1) == 0 ]
    new_order = []
    for d in dividers:
        x = l/d
        while x <= l:
            if x-1 not in new_order: new_order.append(x-1)
            x += l/d
    return new_order    

def trim_gap(al):
    '''trims gaps at the beginning and the end of an alignment'''
    gap = set('-')
    
    # leftward
    i = 0
    while i < al.length and set(al.get_site(i)) == gap:
        i += 1
    new_al = al.get_slice(start=i)
    
    # rightward
    j = new_al.length-1
    while j >= 0 and set(new_al.get_site(j)) == gap:
        j -= 1
    new_al = new_al.get_slice(stop=j+1)
    
    return new_al
    
def poisson_threshold(syn_stats, s_cat=36, lambda_p=0, sfthresh=.99):
    '''returns dS threshold given a Poisson distribution'''
    
    # Only retrieve dS data if lambda has not been set up or if the
    # entered value is negative.
    if lambda_p <= 0:
        n = 0
        
        # Cumulative distribution of the number of synonymous mutations
        # in small windows
        selected_s = ( S*dS for S, dS in syn_stats
            if round(S) == s_cat )        
        f, step = [0], 1
        for s in selected_s:
            n += 1
            categorie = int(floor(s / step))
            f = plot_cat(f, categorie)
        if f == [0]: return None

        # Discard 0 values and values greater than 1, then find the mode, 
        # which will estimate the Poisson distribution's lambda
        # parameter.
        f = dict(zip(f[1:s_cat], range(1, s_cat)))
        if f == {}: return None
        lambda_p = f[sorted(f.keys())[-1]] + 1

        # force the 0 categorie to be higher in the Poisson distribution
        # if lambda is lower than 4
        while (lambda_p < 4) \
        and (f.keys()[0] > exp(-lambda_p) * n):
            lambda_p -= 1

    # calculate the Poisson distribution threshold
    P, k = 0, 0
    while P < sfthresh:
        P += ((lambda_p**k) * (exp(-lambda_p))) / factorial(k)
        k += 1
    else:
        
        # At this stage, the upper limit of acceptable values is k - 2, 
        # so k - 1 is the downer value of unacceptable values. The final 
        # calcul give the maximum value of dS for windows.
        return (k - 1) / float(s_cat)  
