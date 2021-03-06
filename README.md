# PSGfinder
Finding signals of divergent positive selection in pairwise alignments.

PSGfinder analyses pairwise alignments of orthologous CDSs in the search of
regions under positive selection. It uses the program yn00 from the package PAML
4.6-9 (Ziheng Yang 2007) to estimate synonymous and non-synonymous substitution
rates (dS and dN, respectively) of candidate regions previously defined by a
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
  by the PATH environment variable. You can test whether Numpy and python-fisher
  are correctly installed by importing them in Python;
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
Each of these analyses should take about 20 seconds. Information about the number
of analyzed windows and the number of windows with significant signal of positive
selection should be written in the standard error stream. You can compare the
results of these tests (the log files and the result files) with the files in
example/results.

These alignments were kindly provided by Alexandra Weber. They were screened for
signals of divergent positive selection in the study: "Positive selection on
sperm ion channels in a brooding brittle star: consequence of life-history
traits evolution." (Weber et al. 2017).

#### What are these options for?
 - -l or --log allows you to save a log file with a few statistics on your run;
 - -c or --cleaning will make psgfinder.py clean the alignments prior to their
   analysis;
 - -f or --fdr is correcting the p-value for multiple testing by accounting for
   an estimated number of independant tests (RIT in the log file).

#### How did PSGfinder find the input?
By default, the script will look for a directory named 'alignments' in the
current working directory. If you want to set up an other directory for 
screening, use the -a or --alignments option with the name of the directory
in argument.
  ```
  psgfinder.py -a /path/to/alignments/ ...
  ```
## Content                             
This package contains the script psgfinder.py, allowing to screen pairwise
alignments for signal of divergent positive selection as described above. It
also install the module psgfindertools, which contains functions and classes
used in the script.

