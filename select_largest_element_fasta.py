#!/usr/bin/python3
"""select_largest_element_fasta.py is a python script that selects and print in fasta the longest genomic element (like a UTR, coding region etc.) from the same gene out of a multi fasta file.
"""


__version__ = "0.0.1"

# import sys
import argparse
from Bio import SeqIO


parser = argparse.ArgumentParser(prog='select_largest_element_fasta.py', description='Select and print in fasta the longest genomic element (like a UTR, coding region etc.) from the same gene out of a multi fasta file.', epilog="Author: Costas Bouyioukos, 2020, Universite de Paris and UMR7216.")
parser.add_argument('infile', nargs='?', default='-', type=argparse.FileType('r'), metavar="input_file", help='Path to a multi fasta. (or STDIN).')
parser.add_argument("outfile", nargs='?', default='-', type=argparse.FileType('w'), metavar='output_file', help="Path to output FASTA file. (or STDOUT).")
parser.add_argument("-l", "--length", help="Write only sequences longer than 'LEN' (Default=8).", type=int, default=8,
    dest="length", metavar="LEN",)

# Parse the command line arguments.
optArgs = parser.parse_args()

genesElement = {}  # The genes -> elements dictionary
for record in (rec.upper() for rec in SeqIO.parse(optArgs.infile, "fasta")):  # This generator renders all sequences uppercase.
    names = record.id.split()[0].split("|")  # The first field until the first whitespace (FASTA conversion) is selected as ID and then split with | for ENSEMPLE ad so. conversions.
    gene = names[0]  # This is the gene name.
    # Populate the dictionary.
    if gene in genesElement:
        # This condition checks for the longest sequence.
        if len(record.seq) >= len(genesElement[gene].seq):
            genesElement[gene] = record

    elif len(record.seq) >= optArgs.length:
        genesElement[gene] = record
# Write to output
SeqIO.write(list(genesElement.values()), optArgs.outfile, "fasta")
