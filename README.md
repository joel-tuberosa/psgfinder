# PSGfinder
Finding signals of divergent positive selection in pairwise alignments.

PSGfinder analyses pairwise alignments of orthologous CDSs in the search of
regions under positive selection. It uses the program yn00 from the package PAML
4.6-9 (Ziheng Yang 2007) to estimate synonymous and non-synonymous substitution
rates (dS and dS, respectively) of candidate regions previously defined by a
parametrable algorithm. A Fisher's Exact test allows to discriminate positive
selection (dN/dS>1) from neutral evolution (dN/dS=1) and an adapted method of 
correction for multiple testing can be applied (option -f) to prevent false
positives. The program returns a report to the standard output.

## Getting started

### Prerequisites
- [PAML](http://abacus.gene.ucl.ac.uk/software/paml.html)
- [Numpy](http://www.numpy.org/)
- [python-fisher](https://pypi.python.org/pypi/fisher/)

### Installing
1. Install prerequisites. PAML's yn00 executable must be in a directory known
  by the PATH enviromnent variable. You can test whether Numpy and python-fisher
  are correctly installed by importing them in Python.
 2. Execute the setup.py script:
  ```
   python setup.py install
  ```
### Testing
You can start by printing the help message to see if the psgfinder.py script is
correctly installed:
  ```
  psgfinder.py --help
  ```
Then, go in the example directory and run the following commands to match the 
three example outputs:

#### No alignment cleaning, no correction for multiple testing
  ```
  psgfinder.py -l test01-log.txt >test01.txt
  ```

#### Alignment cleaning, no correction for multiple testing
  ```
  psgfinder.py -c -l test02-log.txt >test02.txt
  ```

#### Alignment cleaning, correction for multiple testing
  ```
  psgfinder.py -c -f -l test03-log.txt >test03.txt
  ```

#### What are these options for?
 - -l or --log allows you to save a log file with a few statistics on your run
 - -c or --cleaning will make psgfinder.py clean the alignments prior to their
   analysis
 - -f or --fdr is correcting the p-value for multiple testing by accounting for
   an estimated number of independant tests (RIT in the log file).

## Content                             
This package contains the script psgfinder.py, allowing to screen pairwise
alignments for signal of divergent positive selection as described above. It
also install the module psgfindertools, which contains functions and classes
used in the script.

