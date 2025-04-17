#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""A collection of useful functions to report statistics of for fasta files.

Works as a script if it is called by its name.

@author: Costas Bouyioukos
@organization: The Sainsbury Laboratory
@since: January 2011
@copyright: The program is coming as it is. You have the right to redistribute,
transform and change the source code presuming the appropriate reference and
the license is kept free.
@license: GNU GPL3 or newer.
@contact: U{Costas Bouyioukos<mailto:cbouyio@gmail.com>}
@version: 0.2.1"""


import sys
import re

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt


# A small high order function
def comma_me(amount):
  """Commafying recipe.
  Taken from: http://code.activestate.com/recipes/146461-commafying-an-integer/
  The locale module has a similar functionality too locale.format("%.2f", num, 1)
  """
  orig = amount
  new = re.sub(r"^(-?\d+)(\d{3})", r"\g<1>,\g<2>", amount)
  if orig == new:
    return new
  else:
    return comma_me(new)


# Start the FASTA statistics functions
def parse_fasta(fastaFile, typ = "fasta") :
  """Low level function to return an iterator of sequence record objects.

  This is the primitive function that all the rest of the functions in that
  module should call.
  """
  if typ not in ("fasta", "fastq") :
    raise StandardError("Only fasta or fastq file formats are supported.")
  seqList = []
  with open(fastaFile.name) as hf:
    for seqRec in SeqIO.parse(hf, typ):
      seqList.append(seqRec)
  return seqList

def count_fasta(fastaFile) :
  """Return the number of sequences a fasta file contains.
  """
  return len(parse_fasta(fastaFile))

def count_lengths(fastaFile, sort = False) :
  """Return a list with the lengths of the sequences.
  """
  lengths = []
  for seqRec in parse_fasta(fastaFile) :
    lengths.append(len(seqRec))
  if sort :
    lengths.sort(reverse = True)
  lens = np.array(lengths)
  return lens

def print_lengths(fastaFile, sort = False) :
  """Print a list with the lengths of the sequences.
  
  Elsewhere, pandas.cut() is a convenient way to bin values into arbitrary intervals.
  Let’s say you have some data on ages of individuals and want to bucket them sensibly:
  """
  # TODO I do not know if we need this function.
  print('length')
  for seqLen in count_lengths(fastaFile, sort) :
    print(str(seqLen))

def count_nucleotides(fastaFile) :
  """Return the number of bases the file contains.
  """
  return sum(count_lengths(fastaFile))

def length_histogram(fastaFile, bins) :
  """Return a tuple of frequencies and bin edges.
  """
  lengthsList = []
  for seqRec in parse_fasta(fastaFile) :
    lengthsList.append(len(seqRec))
  return np.histogram(lengthsList, bins)

def min_max_lengths(lengths) :
  """Return the min and max of the sequence lengths.
  """
  return (lengths.max(), lengths.min())

def calculate_N50(fastaFile) :
  """Return the N50 of a fasta file.
  """
  #TODO Include a reference file as an option for the calculation of the N50.
  halfLenTotal = count_nucleotides(fastaFile) / 2.0
  currentLen = 0
  n50 = 0
  for length in count_lengths(fastaFile, sort = True) :
    currentLen = currentLen + length
    if currentLen >= halfLenTotal :
      n50 = length
      break
  return n50

def length_statistics(lengths) :
  """Return some descriptive statistics of the sequence lengths.

  Mean, SD, Median, Coefficient of variation.
  """
  mean = lengths.mean()
  sd = lengths.std()
  median = np.median(lengths)
  coefVar = sd/mean
  iqr = np.quantile(lengths, 0.75) - np.quantile(lengths, 0.25)
  return (mean, sd, coefVar, median, iqr)

def mode(lengths) :
  """Return the mode (the most frequent value of the data.)
  """
  md = stats.mode(lengths, axis=None)
  return md

def print_histogram(lengths, bins) :
  """Return a "pretty" print out of the histogram.
  """
  # Compute the histogram
  hist, bin_edges = np.histogram(lengths, bins = bins)
  # Print the histogram
  for i in range(len(bin_edges) - 1):
    bar = '.' * int(hist[i])
    print(f'{comma_me(str(int(bin_edges[i])))} - {comma_me(str(int(bin_edges[i+1])))}: {bar} ({hist[i]})')

def plot_histogram(lengths, bins):
  """Return a matplotlib figure of the histogram.
  """
  # Compute the histogram
  hist, bin_edges = np.histogram(lengths, bins = bins)
  # Plot the histogram
  plt.figure(figsize=(8, 6))
  plt.stairs(hist, bin_edges, fill=True)
  plt.title('Histogram of contig length')
  plt.xlabel('nts')
  plt.ylabel('Count')
  plt.show()


if __name__ == "__main__" :
  import argparse

  # Command line arguments
  parser = argparse.ArgumentParser(description = 'Python script (and module) to calculate a series of statistics related to sequences contained in a Fasta/q file.')
  parser.add_argument('fasta', nargs = '?', default = '-', type = argparse.FileType('r'), metavar = 'path_fasta/q', help = 'A fasta/q input file path or STDIN. (Default: STDIN)')
  parser.add_argument('out', nargs = '?', default = '-', type = argparse.FileType('w'), metavar = 'path_output', help = 'The output file path or STDOUT. (Default: STDOUT)')
  parser.add_argument('-b', '--bins', metavar = 'no_of_bins', help = "The number of bins for the calculation of the length distribution. (Default: 'sqrt')", default = 'sqrt', dest = 'bins')
  parser.add_argument('-t', '--type', type = str, metavar = 'type_of_file', help = "Designate the type of the Fasta/q file. (Default: fasta)", default = 'fasta', dest = 'tp')

  ns = parser.parse_args()
  bins = ns.bins
  tp = ns.tp
  fastaFile = ns.fasta
  lengthsArray = count_lengths(fastaFile)

  print('File "' + fastaFile.name + '" has the following statistics:')
  print(f'Number of sequences          : {comma_me(str(count_fasta(fastaFile)))}')
  print(f'Number of nucleotides        : {comma_me(str(count_nucleotides(fastaFile)))}')
  mn, mx = min_max_lengths(lengthsArray)
  print(f'Longest Sequence             : {comma_me(str(mx))}')
  print(f'Shortest Sequence            : {comma_me(str(mn))}')
  mean, sd, cv, median, iqr = length_statistics(lengthsArray)
  mean = f'{mean:.2f}'
  sd = f'{sd:.2f}'
  cv = f'{cv:.3f}'
  print(f'Mean sequence length         : {comma_me(str(mean))}')
  print(f'STDV sequence length         : {sd}')
  print(f'Variation of sequence length : {cv}')
  print(f'Median sequence length       : {comma_me(str(median))}')
  print(f'IQR sequence length          : {str(iqr)}')
  print(f'Mode of sequence length      : {str(mode(lengthsArray))}')  # TODO: needs to be computed on a histogram from numpy with bins etc.
  print(f'N50 of the current file      : {comma_me(str(calculate_N50(fastaFile)))}')
  print(f'-Sequence length histogram:')  # TODO: take the proper numpy implementation from here https://numpy.org/doc/stable/reference/routines.statistics.html
  print_histogram(lengthsArray, bins)
  plot_histogram(lengthsArray, bins)
