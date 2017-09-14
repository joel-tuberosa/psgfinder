#!/usr/bin/python

'''Classes and functions to handle fasta formatted files and work with 
coding sequence and amino acid sequence alignments.'''

import re
from StringIO import StringIO

# CLASSES
class eukariot_gc(dict):
    '''translates codon into amino acid according eukariotic genetic code'''
    
    def __init__(self, value={}, mol_type='rna'):
        dict.__init__(value)
        codons = []
        if mol_type == 'dna':
            for first in ('T','C','A','G') :
                for second in ('T','C','A','G') :
                    for third in ('T','C','A','G'):
                        codons += [first + second + third]
        elif mol_type == 'rna':
            for first in ('U','C','A','G') :
                for second in ('U','C','A','G') :
                    for third in ('U','C','A','G'):
                        codons += [first + second + third]
        AA = ['F','F','L','L','S','S','S','S',
              'Y','Y','*','*','C','C','*','W',
              'L','L','L','L','P','P','P','P',
              'H','H','Q','Q','R','R','R','R',
              'I','I','I','M','T','T','T','T',
              'N','N','K','K','S','S','R','R',
              'V','V','V','V','A','A','A','A',
              'D','D','E','E','G','G','G','G']
        self.update(dict(zip(codons, AA)))
        self.update(value)

    def __getitem__(self, codon):
        if 'N' in codon or '?' in codon: return 'X'
        elif codon == '---': return '-'
        elif '-' in codon or len(codon) < 3:
            raise KeyError('incomplete codon: "{}"'.format(codon))
        return dict.__getitem__(self, codon)
        
class Sequences(dict):
    '''handles an unaligned sequences collection'''
        
    def baseinit(self, value={}, mol_type='dna'):
        
        # init from sequences
        if value.__class__.__name__ == self.__class__.__name__:
            dict.__init__(self, value)
            self.names = value.names
            self.mol_type = value.mol_type
            
        # init from dict
        elif type(value) == dict:
            dict.__init__(self, value)
            self.names = value.keys()
            self.mol_type = mol_type
        
        # init from iterable
        else:
            self.names = [] 
            for name, sequence in value:
                if name in self.names:
                    raise ValueError(
                        'two sequences have the same name: "{}"'.format(name))
                self[name] = sequence
            self.mol_type = mol_type
        
        # compute parameters
        if self.mol_type not in ('dna', 'rna', 'prot'):
            raise ValueError('unkown value for mol_type')
        
        # sequences are always upper case
        for name in self.names: self[name] = self[name].upper()
            
    def __init__(self, value={}, mol_type='dna'):
        self.baseinit(value, mol_type)
    
    def __iter__(self): ### allows editing
        return ( (name, self[name]) for name in self.get_names() )
    
    def __setitem__(self, key, value):
        dict.__setitem__(self, key, value)
        if key not in self.names: self.names.append(key)
    
    def __delitem__(self, key):
        dict.__delitem__(self, key)
        self.names.remove(key)
    
    def update(self, arg=None, **kwargs):
        '''updates sequences collection with an other object'''
        
        if arg is None: return
        elif (self.__class__.__name__ == arg.__class__.__name__
            or type(arg) is dict):
            for name, sequence in arg: self[name] = sequence
        else:
            raise TypeError(
                'cannot update {} object with {} object'.format(
                    (self.__class__.__name__, arg.__class__.__name__)))
        for name in kwargs: self[name] = kwargs[name]
    
    def __repr__(self):
        return '<{}.{} instance at {}>'.format(
            __name__, self.__class__.__name__, hex(id(self)))  
    
    def __str__(self):
        return '\n'.join(
            [ (name + '\n' + seq) for name, seq in self ]) + '\n'
    
    def get_names(self): return [ name for name in self.names ]
    
    def get_translated(self, gc=None, ignore_length=True):
        '''translates sequences and return a new object'''
        new = self.__class__(self, mol_type=self.mol_type)
        new.translate(gc, ignore_length)
        return new
   
    def translate(self, gc=None, ignore_length=True):
        '''translates sequences within the object'''
        
        # check the molecule type
        if self.mol_type == 'prot':
            raise AttributeError(
                'protein sequences do not have translate method')
        gc = dna_gc if self.mol_type == 'dna' else rna_gc
        
        # no sequence
        if not self:
            self.mol_type = 'prot'
            return
        
        # translate
        translated = []
        for name, sequence in self:
            
            # if length is ignored and len(sequences)%3 != 0, the last
            # codon (which is incomplete) will not be translated.
            if not ignore_length and len(sequence)%3:
                raise ValueError(
                    'sequence "{}" is not multiple of 3'.format(name))
            translated.append((name, ''.join(( gc[sequence[i:i+3]]
                for i in xrange(0, len(sequence)-(len(sequence)%3), 3) ))))
        self.__init__(translated, mol_type='prot')
        
        # update attributes
        self.mol_type = 'prot'
        if self.__class__.__name__ == 'alignment':
            self.length = len(self.values()[0])
        
class Alignment(Sequences):
    '''handles a sequences alignment'''
    
    def __init__(self, value={}, mol_type='dna'):
        self.baseinit(value, mol_type)
        self.length = self.get_length()
                
    def get_length(self):
        all_len = set(( len(s) for s in self.values() ))
        if len(all_len) > 1:
            raise ValueError('aligned sequences must be the same length')
        return list(all_len).pop() if all_len else None
    
    def get_site(self, index):
        '''returns sites at position "index"'''
        
        return [ self[k][index] for k in self.names ]
    
    def get_slice(self, start=None, stop=None, step=None):
        '''returns a slice of the original alignment'''
        
        s = slice(start, stop, step)
        new = ( (name, sequence[s]) for name, sequence in self )
        return Alignment(new, mol_type=self.mol_type)
     
    def __add__(self, value):
        '''concatenates alignments'''
        
        if not isinstance(value, Alignment):
            raise ValueError('you can only append alignments')
        if set(value.names) != set(self.names):
            raise ValueError('sequences names must be the same')
        return Alignment([ (name, sequence + value[name]) 
            for name, sequence in self ], mol_type=self.mol_type)
    
    def mask(self, start=None, stop=None):
        if start is None:
            if stop is None: length = self.length
            else: length = stop
        else:
            if stop is None: length = self.length-start
            else: length = stop-start
        sl = slice(start, stop, None)
        for name in self.names:
            self[name] = list(self[name])
            self[name][sl] = [ '-' for i in xrange(length) ]
            self[name] = ''.join(self[name])
    
    def update(self, arg=None, **kwargs):
        Sequences.update(self, arg, **kwargs)
        self.length = self.get_length()
        
# GLOBALS - GENETIC CODE
dna_gc = eukariot_gc(mol_type='dna')
rna_gc = eukariot_gc()

# FUNCTIONS - SEQUENCE FILES
def fasta_iter(f):
    '''Yields tuples of names and sequences out of a fasta file'''
    
    p = re.compile('[^A-z?*.-]')
    name, sequence = '', ''
    for line in f:
        if line[0] == '>':
            if name: yield name, sequence
            name, sequence = line.rstrip()[1:], ''
        else: sequence += p.sub('', line)
    if (name, sequence) == ('', ''): return
    yield name, sequence

def handle_fasta(f, mol_type='dna', aligned=False):
    '''retrieves alignment from a FASTA formatted file'''
    
    if isinstance(f, str): f = StringIO(f)
    if aligned: return Alignment(fasta_iter(f), mol_type=mol_type)
    else: return Sequences(fasta_iter(f), mol_type=mol_type)
    
def phylip(al, options='', l=60, s=10):
    '''returns a string containing the alignment in PHYLIP format'''
    
    for xname, x in [('l', l), ('s', s)]:
        if not isinstance(x, int):
            raise TypeError('{} must be an integer, {} found'.format( 
                (xname, type(x).__name__)))
        if x == 0 or x < -1: 
            raise ValueError('{} must be positive or equal to -1'.format(xname))
    
    a_string = ''    
    if options:
        a_string += '  {:d} {:d} {}\n'.format(len(al), al.length, options)
    else: a_string += '  {:d} {:d}\n'.format(len(al), al.length)

    if l == -1 or l >= al.length:
        if s == -1: return a_string + str(al)
        else:
            for name, seq in al:
                seq = list(seq)
                seq = [ seq[j] + ' ' 
                    if j%s == s-1 else seq[j]
                    for j in range(len(seq)) ]
                a_string += '{}\n{}\n\n'.format(name, ''.join(seq))
            return a_string
            
    if 'I' in options:
        for name, seq in al:
            a_string += name + '\n'
        for i in range(0, al.length, l):
            a_string += '{}\n'.format((i + 1))
            for name, seq in al:
                seqslice = seq[i:i+l]
                if s != -1:
                    seqslice = list(seqslice)
                    seqslice = [ seqslice[j] + ' ' 
                        if j%s == s-1 else seqslice[j]
                        for j in range(len(seqslice)) ]
                    seqslice = ''.join(seqslice)
                a_string += seqslice + '\n'
        a_string += '\n'
    
    else:
        for name, seq in al:
            a_string += name + '\n'
            for i in range(0, al.length, l):
                seqslice = seq[i:i+l]
                if s != -1:
                    seqslice = list(seqslice)
                    seqslice = [ seqslice[j] + ' ' 
                        if j%s == s-1 else seqslice[j]
                        for j in range(len(seqslice)) ]
                    seqslice = ''.join(seqslice)
                a_string += seqslice + '\n'
            a_string += '\n'
    
    return a_string

def fasta(*args):
    '''writes sequences to FASTA format'''
    if isinstance(args[0], str):
        name, sequence = args
        seqlines = ( sequence[i:i+60] for i in xrange(0, len(sequence), 60) )
        return '>{}\n{}\n'.format(name, '\n'.join(seqlines))
    else:
        return ''.join(( fasta(name, sequence) for name, sequence in args[0] ))

def stop_codon_in_site(al, i):
    '''returns True when a stop codon is found'''
    
    if i < 0: i = al.length - i
    site = al.get_slice(i, i+3).values()
    for stop in ('TGA', 'TAG', 'TAA'):
        if stop in site: return True
    return False
    
def mask_end_stop_codon(al):
    '''replaces the end stop codon by 3 gap characters'''
    for name in al.get_names():
        j = len(al[name])
        while j >= 3 and al[name][j-3:j] == '---':
            j -= 3

        # check for a stop codon, replace it with gaps
        if al[name][j-3:j] in ('TGA', 'TAG', 'TAA'):
            al[name] = list(al[name])
            al[name][j-3:j] = ['-', '-', '-']
            al[name] = ''.join(al[name])

        # check if there is an incomplete codon
        elif '-' in al[name][j-3]:
            raise ValueError(
                'Codons are not properly aligned, or the sequence length is' + 
                'not multiple of 3.')
    return al

def stop_codon_in_alignment(al):
    '''checks if any sequence of an alignment contains a stop codon'''
    for i in xrange(0, al.length, 3):
        if stop_codon_in_site(al, i): return True
    return False

def map_gap(al):
    '''returns a string that map gap on a sequences of X'''
    return ''.join([ '-' 
        if '-' in al.get_site(i) else 'X'
        for i in range(al.length) ])

def map_gap_coordinates(gap_map):
    '''returns a list of gaps coordinates'''
    c = []
    switch = 0
    for i in range(len(gap_map)):
        if gap_map[i] != '-':
            if switch == 0:
                continue
            else:
                if i != 0: c[-1] += [i]
                switch = 0
        else:
            if switch == 0:
                c.append([i])
                switch = 1
            else:
                continue
    if len(c[-1]) == 1: c[-1].append(len(gap_map))
    return c
        
def map_al_coordinates(gap_map):
    '''returns a list of continuous alignments coordinates'''
    
    c = []
    switch = 0
    for i in range(len(gap_map)):
        if gap_map[i] == '-':
            if switch == 0:
                continue
            else:
                if i != 0: c[-1] += [i]
                switch = 0
        else:
            if switch == 0:
                c.append([i])
                switch = 1
            else:
                continue
    if len(c[-1]) == 1: c[-1].append(len(gap_map))
    return c

def aadiff(al):
    '''returns a list of position of amino acid differences'''
    if al.mol_type in ('dna', 'rna'): tr = al.get_translated()
    elif al.mol_type == 'prot': tr = al
    else: raise ValueError('incompatible mol. type: "{}"'.format(al.mol_type))
    return [ i 
        for i in xrange(tr.length)
        if len(set(tr.get_site(i)) - set(['X', '-', '.'])) > 1 ]
