#!/usr/bin/python

import re

# CLASSES - PARAMETER
class PSGParam:
    '''pack PSGfinder parameters and stats'''

    def __init__(self, **parameters):

        # default parameters for alignment cleaning
        self.cleaning_n = 30
        self.cleaning_dsmax = 1
        self.cleaning_k = 2
        self.cleaning_d = 4
        self.cleaning_q = 6
        
        # default parameters for windows parsing
        self.method = "dwin"
        self.wsize = None # sliding windows' length
        self.wstep = None # sliding windows' step
        self.stand_by = False
        self.no_ds_filter = False
        self.msize = 4  
        self.mmut = 3
        self.w_dsmax = None
        self.s_cat = 36
        self.sfthresh = .99
        self.alpha = 0.05
        self.lambda_p = None
        
        # header of the estimation file
        self.estimf_header = ['file name', 'sequences', 'window', 'size', 'N', 
            'S', '%gap', 't', 'kappa', 'dN', 'dN SE', 'dS', 'dS SE', 'omega',
            'whole gene']
        
        # set up...
        if 'ctl' in parameters.keys():
            with open(parameters['ctl'], 'r') as ctlf: self.read_ctl(ctlf)
        self.config(**parameters)

    def read_ctl(self, f):
        '''read a ctl file containing parameters for the whole procedure'''
        p = {}
        for line in f.readlines():
            line = re.sub('#.*', '', line.strip())
            m = re.match('.+?=.+', line) 
            if m is not None:
                name, value = [ x.strip() for x in m.group().split('=') ]
                p[name] = float(value)
        
        # add values to config...
        self.config(**p)
    
    def config(self, **parameters):

        # parse parameters, using either variable names or names in ctl 
        # file
        for key in parameters:
            value = parameters[key]
            
            # alignment cleaning parameters
            if key in ('cleaning_n', 'n'):
                self.cleaning_n = int(value)
            elif key in ('cleaning_dsmax', 'maximum dS'):
                self.cleaning_dsmax = int(value)
            elif key in ('cleaning_k', 'k'):
                self.cleaning_k = int(value)
            elif key in ('cleaning_d', 'd'):
                self.cleaning_d = int(value)
            elif key in ('cleaning_q', 'q'):
                self.cleaning_q = int(value)
            
            # windows parsing parameters
            elif key in ('msize', 'min. size'): 
                self.msize = int(value)
            elif key in ('mmut', 'min. number of aa diff.'): 
                self.mmut = int(value)  
            elif key in ('w_dsmax', 'dS max'):
                if value is None or value == -1: self.w_dsmax = None
                else: self.w_dsmax = float(value)
            elif key in ('s_cat', 'S category'):
                self.s_cat = int(value)
            elif key in ('lambda_p', 'lambda'):
                if value is None or value == -1: self.lambda_p = None
                else: self.lambda_p = float(value)
            elif key in ('sfthresh', 'static filter threshold'):
                self.sfthresh = float(value)
            elif key in ('significance_threshold', 'alpha'):
                self.alpha = float(value)
            
            # flow options
            elif key == 'stand_by':
                self.stand_by = bool(value)
            elif key == 'no_test':
                self.no_test = bool(value)
            elif key == 'no_ds_filter':
                self.no_ds_filter = bool(value)
            
            # unknown parameter
            elif key != 'ctl':
                raise TypeError('unknown parameter: "{}"'.format(key))

