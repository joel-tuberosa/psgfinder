#!/usr/bin/python

'''
Tools for PSGfinder

    cleaning        Methods for nucleotide alignment cleaning.
                        
                        functions:
                        - pairwise_cleaning
                        - aadiff_cleaning
                        - ds_cleaning
                        - gap_cleaning
    
    dna             Classes and methods for FASTA file processing, nucleotide
                     alignment handling and coding sequence translation.
                        
                        classes
                        - Alignment
                        - eukariot_gc
                        - Sequences
                        
                        functions
                        - aadiff
                        - fasta
                        - fasta_iter
                        - handle_fasta
                        - map_al_coordinates
                        - map_gap
                        - map_gap_coordinates
                        - mask_end_stop_codon
                        - phylip
                        - stop_codon_in_alignment
                        - stop_codon_in_site

    dwin            Implementation of the Dwin (windows defined between amino
                     acid differences) and the sliding windows methods.
                     
                        functions
                        - dnds_stat
                        - map_win
                        - parse
                        - read_yn00output
                        - setup_method
                        - sliding_windows
                        - window_slider
                        - yn00_pairwise
                     
    fdr             Methods for the estimation of the degree of independance 
                     within a set of partially overlapping tests, intended for
                     multiple comparisons correction.
                        
                        functions
                        - fdr
                        - rho
                    
    misc            Misceallenous operations (temporary file handling, poisson
                     distribution, parsing order...)
                     
                        functions
                        - is_number
                        - parsing_order
                        - plot_cat
                        - poisson_threshold
                        - trim_gap
                        - 
    
    parameters      Object (PSGParam) that contains all parameter values for
                     PSGfinder.
                     
                        class
                        - PSGParam
'''
