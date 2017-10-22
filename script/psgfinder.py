#!/usr/bin/python

'''
USAGE
    psgfinder.py [OPTION] [CONTROL FILE]

DESCRIPTION
    PSGfinder analyses pairwise alignments of orthologous CDSs in the
    search of regions under positive selection. It uses the program yn00 
    from the package PAML 4.6-9 (Ziheng Yang 2007) to estimate synonymous
    and non-synonymous substitution rates (dS and dS, respectively) of 
    candidate regions previously defined by a parametrable algorithm. A  
    Fisher's Exact test allows to discriminate positive selection  
    (dN/dS>1) from neutral evolution (dN/dS=1) and an adapted method of
    correction for multiple testing can be applied (option -f) to prevent 
    false positives. The program returns a report to the standard output.

    Additionnaly, psgfinder implements a method of alignment cleaning 
    (option -c),masking regions of unsure orthology.
    
    In order to manipulate parameter values of the various algorithms
    (window definition, alignment cleaning, p-value threshold), PSGfinder
    will read psgfinder.ctl if found in the working directory. This
    control file is a text file with the following lines (ignoring lines
    starting with '%'): 
%------------------------------------------------------------------------ 
    # alignment cleaning parameters 
                              n = 30 
                     maximum dS = 1 
                              k = 2 
                              d = 4 
                              q = 6 

    # windows parsing and filterng parameters 
                      min. size = 4 
        min. number of aa diff. = 3 
                         dS max = -1 
                     S category = 36 
                         lambda = -1 
                          alpha = .05 
%------------------------------------------------------------------------
 
OPTIONS
    -a, --alignments=DIR        Specify alignments directory
                                 (default=alignments).
    -c, --clean-up              Clean alignments before analysis.
        --save-cleaned=PREFIX   Save the cleaned alignments in the 
                                 alignment directory with the given
                                 PREFIX.
    -e, --estimations=FILE      Locate an estimations FILE and analyze it
                                 (1).
    -f, --fdr                   Apply a correction for multiple testing
                                 (2).
    -l, --log-file=FILE         Write parameters and various stats to
                                 FILE.
    -r, --data-range=X-Y        Define a range of files to be analysed
                                 (from X to Y, where X and Y are
                                 respectively the Xth and the Yth file
                                 found in the alignments and the Yth file
                                 found in the alignments directory. You 
                                 can also specify one number which will 
                                 correspond to one file.
    --sliding-windows=X:Y       Run a sliding windows analysis, with
                                 window length X and step Y
    -x, --no-test               Do not test positive selection, only
                                 estimate substitution rates and print
                                 results to standard output.
    -h, --help                  Display this message.
    
    
    (1) If the estimation file contains a header, use option --header.
    (2) When using --fdr option combined with --estimations, you must
         provide the alignment folder because the program needs to parse 
         the alignment to calculate the relative number of independant
         tests.
    
'''

import sys
import os
import tempfile
import getopt

from psgfindertools.dna import *
from psgfindertools.dwin import *
from psgfindertools.cleaning import *
from psgfindertools.fdr import *
from psgfindertools.misc import *
from psgfindertools.parameters import *

# OPTION HANDLING
class Options(object):

    def __init__(self, argv):

        '''read options and arguments with getopt and edit config'''

        # set default
        self.set_default()
        
        # handle options with getopt
        try:
            options_names = ['alignments=', 'clean-up',  'save-cleaned=', 
            'estimations=', 'fdr', 'log-file=', 'data-range=', 'no-test',
            'sliding-windows=', 'help']
            opts, args = getopt.getopt(argv[1:], 'a:ce:fl:r:xh', options_names)
        except getopt.GetoptError, e:
            sys.stderr.write(str(e) + '\n\n' + __doc__)
            sys.exit(1)

        for o, a in opts:
            if o in ('-a', '--alignments'):
                self.alignments_path = os.path.abspath(a) + '/'
                p, d = os.path.split(a.rstrip('/'))
                self.cleaned_alignments_path = os.path.abspath(
                    os.path.join(p, 'cleaned-' + d)) + '/'
            elif o in ('-c', '--clean-up'):
                self.apply_cleaning = True
            elif o == '--save-cleaned':
                self.save_cleaned = a
            elif o in ('-e', '--estimations'):
                self.estimations_fname = a
            elif o in ('-f', '--fdr'):
                self.apply_fdr = True
            elif o in ('-l', '--log-file'):
                self.log_file = a
            elif o in ('-r', '--data-range'):
                self.data_range = readrange(a)
            elif o == ('--sliding-windows'):
                self.method = "sliding_windows"
                self.wsize, self.wstep = self.get_sliding_windows_parameters(a)
            elif o in ('-x', '--no-test'):
                self.no_test = True
            elif o in ('-h', '--help'):
                sys.stdout.write(__doc__)
                sys.exit(0)
            else:
                assert False, 'unhandled option: {}'.format(o)

        self.args = args
    
    def get_sliding_windows_parameters(self, a):
        wsize, wstep = a.strip().split(":")
        wsize, wstep = int(wsize), int(wstep)
        for name, x in (("wsize", wsize), ("wstep", wstep)):
            if x < 1: raise ValueError(
                "{} must be greater than 1".format(name))
        return (wsize, wstep)
                
    def set_default(self):
    
        cwd = os.getcwd()
        
        # default value for options
        self.alignments_path = cwd + '/alignments/'
        self.cleaned_alignments_path = cwd + '/cleaned_alignments/'
        self.estimations_fname = None
        self.log_file = None
        self.apply_cleaning = False
        self.save_cleaned = None
        self.save_estimations = False
        self.apply_fdr = False
        self.no_test = False
        self.data_range = slice(None)
        self.method = 'dwin'
        self.wsize = None
        self.wstep = None

# GLOBALS - STATS FOR LOG
stats = {'n_al': 0,
    'n_win': 0,
    'dN': [],
    'dS': [],
    'RIT': None,
    'alpha': None,
    'skipped': []}
        
# FUNCTIONS - LOG UPDATE AND WRITING
def write_log_file(fname):
    global stats
    
    with open(fname, 'w') as f:
        f.write(('''### PSGfinder logfile
    aligments analyzed          {}
    windows analyzed            {}
    average dN (whole gene)     {}
    average dS (whole gene)     {}
    RIT                         {}
    p-value threshold           {}\n'''.format(
        stats['n_al'], 
        stats['n_win'],
        stats['dS'], 
        stats['dN'], 
        stats['RIT'], 
        stats['alpha'])
        ).replace('None', 'nan'))
        
        f.write('''
    skipped alignment files
        {}\n'''.format('\n        '.join(stats['skipped'])))

# FUNCTIONS - RANGE EXPRESSION INTERPRETER
def readrange(s):
    '''inteprets a range expression such as "X-Y"'''
    s = s.strip()
    if s.find('-') == -1: return slice(int(s)-1, int(s))
    else:
        if s.count('-') == 1:
            x, y = s.split('-')
            if s[0] == '-': return slice(int(y))
            elif s[-1] == '-': return slice(int(x), None)
            else: return slice(int(x), int(y))
        else: raise ValueError('range must contain at most one "-"')

# FUNCTIONS - MAIN        
def main(argv=sys.argv):

    global stats
    
    # Get parameter values if psgfinder.ctl is present in the current
    # directory, otherwise use default; get options from the command 
    # line.
    if os.path.isfile('psgfinder.ctl'):
        config = PSGParam(ctl='psgfinder.ctl')
    else: config = PSGParam()
    options = Options(argv)
    
    # set up method, DWin or sliding windows
    config.method = options.method
    config.wsize = options.wsize
    config.wstep = options.wstep

    # The file object estimf is used to store values estimated with yn00
    # for every window. The default procedure is 1) to open it as a
    # temporary file; 2) read alignment files, parse windows and write
    # all yn00 results in that file; 3) recatch dS values to calculate a
    # filtering criterion; 4) filter and compute the final dN/dS value
    # and the associated p-value.
    
    # With option -x, the script write the values calculated with yn00 in
    # the standard output instead and does not compute the final dN/dS
    # value and the associated p-value. In this case, the temporary file
    # stays empty.
    
    # With option -e, the script does not look at alignments but recatch
    # the content of a file produced with psgfinder using the option -x.
    # In this case, estimf is not a temporary file, but its location is
    # defined by the option's value.
    with (tempfile.TemporaryFile(prefix='psgfinder_') 
        
          ### ... write values to a temp file
          if options.estimations_fname is None
        
          ### ... or open the file define with option -e
          else open(options.estimations_fname, 'U')) as estimf: 
        
        # (option -x) Print a header for values estimated with yn00.
        if options.no_test:
            sys.stdout.write('\t'.join(config.estimf_header) + '\n')
            sys.stdout.flush()
        estimations = []
            
        # If the option -e is not defined, get all files names to 
        # work with.
        if options.estimations_fname is None: 
            al_files_names = os.listdir(options.alignments_path)
            al_files_names = al_files_names[options.data_range]
            stats['n_al'] = len(al_files_names)
            
            # Make a directory (if not already existing) for cleaned 
            # alignments.
            if options.apply_cleaning and \
               options.save_cleaned is not None:
                if not os.path.isdir(options.cleaned_alignments_path):
                    os.mkdir(options.cleaned_alignments_path)
        
        # ... otherwise make an empty list to jump this step later
        # on.
        else: al_files_names = []
        
        # If al_files_names is empty: no loop; otherwise, process 
        # alignments one by one (option -x or default procedure)
        for fname in al_files_names:
            with open(options.alignments_path + fname, 'r') as f:
                al = handle_fasta(f, aligned=True)
            if len(al) != 2:
                raise TypeError(
                    '{} does not contain a pairwise alignment'. format(name))
            
            # Mask end stop codons
            try: 
                al = mask_end_stop_codon(al)
            except ValueError, e:
                sys.stderr.write(e + 
                    '\n\nWarning: {} will be skipped.\n'.format(fname))
                stats['n_al'] -= 1
                stats['skipped'].append(fname)
                continue                    

            # check for stop codons
            if stop_codon_in_alignment(al):
                sys.stderr.write(
                    'Warning: {} will be skipped because it contains an in' +
                    ' frame stop codon\n'.format(fname))
                stats['n_al'] -= 1
                stats['skipped'].append(fname)
                continue
            
            # check if the alignment length is multiple of 3.
            if al.length%3 != 0:
                sys.stderr.write(
                    'Warning: {} will be skipped because the alignment length' +
                    ' is not multiple of 3.\n'.format(fname))
                stats['n_al'] -= 1
                stats['skipped'].append(fname)
                continue
                
            # clean alignments if needed
            if options.apply_cleaning:
                try:
                    al = pairwise_cleaning(al, config.cleaning_n,
                        config.cleaning_dsmax, config.cleaning_k,
                        config.cleaning_d, config.cleaning_q)
                except:
                    sys.stderr.write(
                        'problem with {}'.format(fname) +
                        'when cleaning alignment.\n')
                    raise
                
                # save the cleaned alignment to a file with the given
                # prefix.
                if options.save_cleaned is not None:
                    cleaned_al_fname = options.cleaned_alignments_path + \
                        '{}{}'.format(options.save_cleaned, fname)
                    with open(cleaned_al_fname, 'w') as fout:
                        fout.write(fasta(al))
                
            # estimation on candidate windows, write estimations in
            # stdout if options.no_test is True (option -x)
            try:
                parse(al, 
                    values=[],
                    fname=fname,
                    estimf=estimf if not options.no_test else sys.stdout,
                    config=config,
                    stats=stats)
                
            except ValueError, e:
                sys.stderr.write(
                    'problem with {}'.format(fname) + 
                    ' when estimating windows substitution rates\n' + str(e) + 
                    ' \n')
                raise

        if options.no_test: return 0
      
        # when finish, recatch estimf content
        estimf.seek(0)
        
        # ignore the header line
        if options.estimations_fname:
            header = estimf.readline()
        
        if config.w_dsmax is None:
            Si = config.estimf_header.index('S')
            dSi = config.estimf_header.index('dS')
            syn_stats = []
            for line in estimf:
                line = line.rstrip().split('\t')
                syn_stats.append((float(line[Si]), float(line[dSi])))
            config.w_dsmax = poisson_threshold(syn_stats)
            if config.w_dsmax is None:
                sys.stderr.write(
                    'Warning: psgfinder could not estimate windows maximum' +
                    ' dS from dataset, the threshold will be 1\n')
                config.w_dsmax = 1.0
            estimf.seek(0)
            if options.estimations_fname:
                header = estimf.readline()
        
        estimations = [ 
            dict(zip(config.estimf_header, line.rstrip().split('\t')))
            for line in estimf ]
                
        estimations = [ x for x in estimations 
            if x['whole gene'] == '1' or float(x['dS']) < config.w_dsmax ]
                
    stats['dN'] = sum(stats['dN'])/len(stats['dN']) if \
        len(stats['dN']) else float('nan')
    stats['dS'] = sum(stats['dS'])/len(stats['dS']) if \
        len(stats['dS']) else float('nan')
    
    # calculate dNdS and keep windows with dN(window) > dS(gene), 
    # calculate a p-value for H0(dN = dS)
    if estimations:
        estimations = dnds_stat(estimations)
        number_of_windows = sum(( 1 for x in estimations 
            if x['whole gene'] == '0' ))
        if number_of_windows:
            sys.stderr.write(
                ' ... Number of windows with dN/dS > 1: {}\n'.format(
                number_of_windows))
        else:
            sys.stderr.write(' ... No window found with dN/dS > 1.\n')
            if options.log_file is not None: write_log_file(options.log_file)
            return 0            
    
    # otherwise, write log file and exit
    else:
        sys.stderr.write(' ... No window found with dN/dS > 1.\n')
        if options.log_file is not None: write_log_file(options.log_file)
        return 0
    
    # keep genes that have at least one candidate window
    genes_to_keep = { x['file name'] for x in estimations
        if x['whole gene'] == '0' }
    estimations = [ x for x in estimations 
        if x['file name'] in genes_to_keep ]
    
    # estimate relative number of independant tests and correct p-value
    if options.apply_fdr:
        
        # attach alignments to whole genes
        for i in xrange(len(estimations)):
            if estimations[i]['whole gene'] == '1':
                fname = options.alignments_path + estimations[i]['file name']
                with open(fname, 'r') as f:
                    s = handle_fasta(f, aligned=True)
                estimations[i]['nucleotides'] = s
        
        # apply false discovery rate correction to p-value
        config.alpha = fdr(estimations, config.alpha, stats=stats)
    
    # add alpha to stat
    stats['alpha'] = config.alpha

    # filter significant results
    estimations = [ x for x in estimations 
        if x['whole gene'] == '1' or x['p-value'] < config.alpha ]
    number_of_windows = sum(( 1 for x in estimations 
            if x['whole gene'] == '0' ))
    if number_of_windows:
        sys.stderr.write(
            ' ... Number of windows with significant dN/dS > 1: {}\n'.format(
            number_of_windows))
    else:
        sys.stderr.write(' ... No window found with significant dN/dS > 1.\n')
        if options.log_file is not None: write_log_file(options.log_file)
        return 0
        
    # (again...) keep genes that have at least one candidate window
    genes_to_keep = { x['file name'] for x in estimations
        if x['whole gene'] == '0' }
    estimations = ( x for x in estimations if x['file name'] in genes_to_keep)
    
    # write log file
    if options.log_file is not None: write_log_file(options.log_file)
    
    # print results to standard ouput
    header = [ 'file name', 'sequences', 'window', 'size', 'N', 'S', '%gap',
        'dN', 'dS', 'dN/dS(whole gene)', 'p-value', 'whole gene' ]
    sys.stdout.write('\t'.join(header) + '\n')
    for x in estimations:
        sys.stdout.write('\t'.join(( str(x[e]) for e in header)) + '\n')
        try: sys.stdout.flush()
        except IOError, e: ### handle SIGTERM
            if str(e) == '[Errno 32] Broken pipe': return 0
            else: raise IOError, e
    return 0
    
# RUN MAIN FUNCTION
if __name__ == '__main__':
    sys.exit(main())
