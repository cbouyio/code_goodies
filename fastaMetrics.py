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


def parse_fasta(fastaFile, typ = "fasta"):
  """Low level function to return an iterator of sequence record objects.

  This is the primitive function that all the rest of the functions in that
  module should call."""
  if typ not in ("fasta", "fastq") :
    raise StandardError("Only fasta or fastq file formats are supported.")
  seqList = []
  with open(fastaFile.name) as hf:
    seqList.extend(iter(SeqIO.parse(hf, typ)))
  return seqList

def count_fasta(fastaFile) :
  """Return the number of sequences a fasta file contains."""
  return len(parse_fasta(fastaFile))

def count_lengths(fastaFile, sort = False):
  """Return a list with the lengths of the sequences."""
  dist = [len(seqRec) for seqRec in parse_fasta(fastaFile)]
  if sort :
    dist.sort(reverse = True)
  return dist

def print_lengths(fastaFile, sort = False):
  """Print a list with the lengths of the sequences. Elsewhere, pandas.cut() is a convenient way to bin values into arbitrary intervals. Let’s say you have some data on ages of individuals and want to bucket them sensibly:"""
  # TODO I do not know if we need this function.
  print('length')
  for seqLen in count_lengths(fastaFile, sort):
    print(seqLen)

def count_nucleotides(fastaFile) :
  """Return the number of bases the file contains."""
  return sum(count_lengths(fastaFile))

def length_histogram(fastaFile, bins):
  """Return a tuple of frequencies and bin edges."""
  lengthsList = [len(seqRec) for seqRec in parse_fasta(fastaFile)]
  return np.histogram(lengthsList, bins)

def print_histogram(fastaFile, bin = 'auto'):
  """Return a "pretty" print out of the histogram."""
  #Compute the bin centres
  hist, bins = length_histogram(fastaFile, bin)
  # binRanges = []
  # for i in xrange(bins) :
  #   binRange = [int(round(hist[1] + i * hist[2], 0)), int(round(hist[1] + (i + 1) * hist[2], 0))]
  #   if binRange[0] < 0.0 :
  #     binRange[0] = 0
  #   binRanges.append(binRange)
  # freqs = hist[0]
  maxDigits1 = 2 * len(comma_me(str(bins[-1]))) + 1
  maxDigits2 = len(comma_me(str(hist[-1]))) + 1
  print('Length'.center(maxDigits1), ':', 'Counts (freq)'.center(maxDigits2), ':', 'Hist'.center(maxDigits1 + maxDigits2))
  for i in range(len(hist)):
    binStr = f"{comma_me(str(int(bins[i])))}-{comma_me(str(int(bins[i + 1]) - 1))}"
    freqStr = comma_me(str(hist[i])) + " ({0:.3f})".format(hist[i]/float(count_fasta(fastaFile)))
    histStr = f"{'+' * hist[i]}"
    print(binStr.rjust(maxDigits1), ':', freqStr.rjust(maxDigits2), ':', histStr.ljust(maxDigits1 + maxDigits2))

def min_max_lengths(fastaFile) :
  """Return the min and max of the sequence lengths."""
  lengths = count_lengths(fastaFile, sort = True)
  return (lengths[-1], lengths[0])

def calculate_N50(fastaFile) :
  """Return the N50 of a fasta file."""
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

def mean_sd_median_lengths(fastaFile):
  """Return some descriptive statistics of the sequence lengths.

  Mean, SD, Median, Coefficient of variation.
  """
  lengths = count_lengths(fastaFile)

def length_statistics(fastaFile) :
  """Return some descriptive statistics of the sequence lengths.

  Mean, SD, Median, Coefficient of variation.
  """
  lengths = count_lengths(fastaFile)
  mean = np.mean(lengths)
  sd = np.std(lengths)
  median = np.median(lengths)
  coefVar = sd/mean
  iqr = np.quantile(lengths, 0.75) - np.quantile(lengths, 0.25)
  return (mean, sd, coefVar, median, iqr)

def mode(fastaFile):
  """Return the mode (the most frequent value of the data.)"""
  lengths = count_lengths(fastaFile)
  return stats.mode(lengths, axis=None)

def comma_me(amount):
  """Commafying recipe.
  Taken from: http://code.activestate.com/recipes/146461-commafying-an-integer/
  The locale module has a similar functionality too   locale.format("%.2f", num, 1)"""
  orig = amount
  new = re.sub("^(-?\d+)(\d{3})", '\g<1>,\g<2>', amount)
  return new if orig == new else comma_me(new)


if __name__ == "__main__":
  import argparse

  parser = argparse.ArgumentParser(description = 'Python script (and module) to calculate a series of statistics related to sequences contained in a Fasta/q file.')
  parser.add_argument('fasta', nargs = '?', default = '-', type = argparse.FileType('r'), metavar = 'path_fasta/q', help = 'A fasta/q input file path or STDIN. (Default: STDIN)')
  parser.add_argument('out', nargs = '?', default = '-', type = argparse.FileType('w'), metavar = 'path_output', help = 'The output file path or STDOUT. (Default: STDOUT)')
  parser.add_argument('-b', '--bins', type = int, metavar = 'no_of_bins', help = "The number of bins for the calculation of the length distribution. (Default: 'auto')", default = 10, dest = 'bins')
  parser.add_argument('-t', '--type', type = str, metavar = 'type_of_file', help = "Designate the type of the Fasta/q file. (Default: fasta)", default = 'fasta', dest = 'tp')

  ns = parser.parse_args()
  bins = ns.bins
  tp = ns.tp
  fastaFile = ns.fasta

  print(f'File "{fastaFile.name}" has the following statistics:')
  print(f'Number of sequences          : {comma_me(str(count_fasta(fastaFile)))}')
  print(f'Number of nucleotides        : {comma_me(str(count_nucleotides(fastaFile)))}')
  mn, mx = min_max_lengths(fastaFile)
  print(f'Longest Sequence             : {comma_me(str(mx))}')
  print(f'Shortest Sequence            : {comma_me(str(mn))}')
  mean, sd, cv, median, iqr = length_statistics(fastaFile)
  sd = f'{sd:.1f}'
  cv = f'{cv:.4f}'
  print(f'Mean sequence length         : {comma_me(str(mean))}')
  print(f'STDV sequence length         : {comma_me(sd)}')
  print(f'Variation of sequence length : {comma_me(cv)}')
  print(f'Median sequence length       : {comma_me(str(median))}')
  print(f'IQR sequence length          : {comma_me(str(iqr))}')
  #print('Mode of sequence length      : %s' % comma_me(str(mode(fastaFile))))  # TODO: needs to be computed on a histogram from numpy with bins etc.
  print(f'N50 of the current file      : {comma_me(str(calculate_N50(fastaFile)))}')
  print('-Sequence length histogram:')
  print_histogram(fastaFile, bins)
