#!/usr/bin/env python3
"""A python tool module that exends with some functionalities the CodonUsage module of BioPython and provides interface to calculate CAI (Codon Adaptation Index) from a fasta file.

Runs as a command line application, or can be imported as a module.

Author: Costas Bouyioukos, @UMR7216, 2019
"""

# TODO develope it as a biopython extension to calculate and update the CAI indexices from direct connection to ENSEMBL.! An interesting addition.

import argparse
import sys
#import datetime
import pandas as pd
from Bio import SeqIO
from Bio.SeqUtils import GC
from Bio.SeqUtils import CodonUsage

# CONSTANTS
# Here we store (HARDCODED) the dictionaries of codon counts for various organsims.
CODON_COUNT_HUMAN = {"TTT": 699720, "TTC": 794184, "TTA": 320356, "TTG": 535679, "CTT": 576736, "CTC": 785652, "CTA": 294189, "CTG": 1596864, "ATT": 654002, "ATC": 829032, "ATA": 307388, "ATG": 913885, "GTT": 461894, "GTC": 602589, "GTA": 303612, "GTG": 1133771, "TAT": 494637, "TAC": 598949, "TAA": 22853, "TAG": 18231, "CAT": 462703, "CAC": 624310, "CAA": 529770, "CAG": 1449133, "AAT": 707587, "AAC": 771805, "AAA": 1035431, "AAG": 1324003, "GAT": 939387, "GAC": 1053184, "GAA": 1273418, "GAG": 1674275, "TCT": 652998, "TCC": 766054, "TCA": 554595, "TCG": 182321, "CCT": 794065, "CCC": 827321, "CCA": 740822, "CCG": 283264, "ACT": 563608, "ACC": 806532, "ACA": 664892, "ACG": 248351, "GCT": 775658, "GCC": 1137563, "GCA": 688457, "GCG": 287751, "TGT": 426306, "TGC": 489888, "TGA": 45761, "TGG": 513091, "CGT": 193933, "CGC": 408392, "CGA": 265448, "CGG": 476812, "AGT": 529539, "AGC": 821479, "AGA": 511837, "AGG": 495809, "GGT": 460219, "GGC": 884845, "GGA": 701790, "GGG": 676009}  # Calculated from Homo_sapiens.GRCh38.cds.all.fa

CODON_COUNT_MOUSE = {'TTT': 487247, 'TTC': 583751, 'TTA': 206411, 'TTG': 390039, 'CTT': 391179, 'CTC': 555903, 'CTA': 231967, 'CTG': 1090334, 'ATT': 443375, 'ATC': 610734, 'ATA': 217617, 'ATG': 650659, 'GTT': 314637, 'GTC': 426050, 'GTA': 217110, 'GTG': 781639, 'TAT': 338951, 'TAC': 431874, 'TAA': 15337, 'TAG': 12576, 'CAT': 318108, 'CAC': 430402, 'CAA': 360443, 'CAG': 1001063, 'AAT': 458748, 'AAC': 568132, 'AAA': 664865, 'AAG': 961342, 'GAT': 625070, 'GAC': 740154, 'GAA': 830718, 'GAG': 1146504, 'TCT': 483182, 'TCC': 514643, 'TCA': 362254, 'TCG': 120745, 'CCT': 543418, 'CCC': 509721, 'CCA': 517421, 'CCG': 172796, 'ACT': 398121, 'ACC': 518374, 'ACA': 471543, 'ACG': 159765, 'GCT': 574099, 'GCC': 722716, 'GCA': 468071, 'GCG': 176080, 'TGT': 319358, 'TGC': 330180, 'TGA': 28664, 'TGG': 339642, 'CGT': 132247, 'CGC': 252169, 'CGA': 192601, 'CGG': 293090, 'AGT': 382348, 'AGC': 567082, 'AGA': 365636, 'AGG': 356903, 'GGT': 318192, 'GGC': 573158, 'GGA': 480814, 'GGG': 423631}  # Calculated from Mus_musculus.GRCm38.cds.all.fa

CODON_COUNT_YEAST = {"TTT": 79148, "TTC": 53944, "TTA": 77583, "TTG": 78560, "CTT": 36968, "CTC": 16801, "CTA": 39392, "CTG": 31703, "ATT": 88446, "ATC": 49093, "ATA": 53840, "ATG": 61056, "GTT": 63153, "GTC": 32925, "GTA": 35748, "GTG": 32143, "TAT": 55653, "TAC": 42041, "TAA": 3138, "TAG": 1511, "CAT": 40077, "CAC": 22270, "CAA": 77278, "CAG": 35917, "AAT": 105623, "AAC": 71012, "AAA": 123449, "AAG": 88653, "GAT": 109757, "GAC": 58078, "GAA": 132048, "GAG": 56930, "TCT": 68479, "TCC": 41295, "TCA": 55198, "TCG": 25761, "CCT": 38941, "CCC": 20258, "CCA": 51045, "CCG": 15538, "ACT": 58292, "ACC": 36147, "ACA": 51798, "ACG": 23982, "GCT": 58801, "GCC": 35734, "GCA": 47400, "GCG": 18449, "TGT": 24113, "TGC": 14726, "TGA": 2079, "TGG": 30566, "CGT": 18305, "CGC": 7918, "CGA": 9151, "CGG": 5641, "AGT": 42680, "AGC": 29642, "AGA": 61208, "AGG": 28019, "GGT": 65720, "GGC": 28880, "GGA": 32778, "GGG": 18168}  # Calculated from Saccharomyces_cerevisiae.R64-1-1.cds.all.fa

CurrentCodonIndex = {}


def generate_codon_dict_fasta(fasta_file):
    """Calculate a codon counts dictionary from a fasta file.

    Return the dict.
    """
    codon_counts = {'TTT': 0, 'TTC': 0, 'TTA': 0, 'TTG': 0,
                    'CTT': 0, 'CTC': 0, 'CTA': 0, 'CTG': 0,
                    'ATT': 0, 'ATC': 0, 'ATA': 0, 'ATG': 0,
                    'GTT': 0, 'GTC': 0, 'GTA': 0, 'GTG': 0,
                    'TAT': 0, 'TAC': 0, 'TAA': 0, 'TAG': 0,
                    'CAT': 0, 'CAC': 0, 'CAA': 0, 'CAG': 0,
                    'AAT': 0, 'AAC': 0, 'AAA': 0, 'AAG': 0,
                    'GAT': 0, 'GAC': 0, 'GAA': 0, 'GAG': 0,
                    'TCT': 0, 'TCC': 0, 'TCA': 0, 'TCG': 0,
                    'CCT': 0, 'CCC': 0, 'CCA': 0, 'CCG': 0,
                    'ACT': 0, 'ACC': 0, 'ACA': 0, 'ACG': 0,
                    'GCT': 0, 'GCC': 0, 'GCA': 0, 'GCG': 0,
                    'TGT': 0, 'TGC': 0, 'TGA': 0, 'TGG': 0,
                    'CGT': 0, 'CGC': 0, 'CGA': 0, 'CGG': 0,
                    'AGT': 0, 'AGC': 0, 'AGA': 0, 'AGG': 0,
                    'GGT': 0, 'GGC': 0, 'GGA': 0, 'GGG': 0}
    with fasta_file as fh:
        for rec in SeqIO.parse(fh, "fasta"):
            # Do not take into account very short ORFs.
            if len(rec.seq) < 45:
                continue;
            # Make all UPPERcase.
            dna_seq = str(rec.seq).upper() if str(rec.seq).islower() else str(rec.seq)
            for i in range(0, len(dna_seq), 3):
                codon = dna_seq[i:i+3]
                if codon in codon_counts:
                    codon_counts[codon] += 1
                elif len(codon) >= 3 and "N" not in codon:
                    raise TypeError(f"Illegal codon {codon} in gene: {rec.id}")
    #with open("codon_counts_" + datetime.datetime.now().strftime('%Y%m%d_%H%M%S') + ".json", "w") as js:
    #    json.dump(codon_counts, js)
    return codon_counts


def generate_index_dict(organism, dict_ext = None):
    """Generate a codon usage index from a codon usage dictionary.

    Calculate the codon index as the CurrentCodonIndex module constant.

    RCSU values.
    """
    # TODO ABSOLUTELY INCORPORATE ALL THIS IN BIOPYTHON.
    #To calculate the index we first need to sum the number of times synonymous codons were used all together.
    codon_dict = {}
    if organism == "human":
        codon_dict = CODON_COUNT_HUMAN
    elif organism == "mouse":
        codon_dict = CODON_COUNT_MOUSE
    elif organism == "yeast":
        codon_dict = CODON_COUNT_YEAST
    elif dict_ext:
        codon_dict = dict_ext
    else:
        raise ValueError(f"Organism {organism} is not supported")
    for aa in CodonUsage.SynonymousCodons:
        total = 0.0
        # RCSU values are CodonCount/((1/num of synonymous codons) * sum of all synonymous codons)

        codons = CodonUsage.SynonymousCodons[aa]
        for codon in codons:
            total += codon_dict[codon]
        # calculate the RSCU value for each of the codons
        rcsu = []
        for codon in codons:
            denominator = float(total) / len(codons)
            rcsu.append(codon_dict[codon] / denominator)
        # now generate the index W=RCSUi/RCSUmax:
        rcsu_max = max(rcsu)
        for codon_index, codon in enumerate(codons):
            CurrentCodonIndex[codon] = rcsu[codon_index] / rcsu_max


def calculate_CAI(file, cci=CurrentCodonIndex):
    """Calculate the Codon Adaptation Index from a FASTA file.

    Return a data frame with it.
    """
    k =0
    caidf = pd.DataFrame(columns=["CAI"])
    SeqCai = CodonUsage.CodonAdaptationIndex()
    SeqCai.set_cai_index(cci)
    for seq in SeqIO.parse(file, "fasta"):
        k+=1
        cai = SeqCai.cai_for_gene(str(seq.seq))
        idt = seq.id
        caidf.loc[idt] = float(cai)
    print(f"Total seqs: {k}")
    print(f"Total number of wrong codons: {SeqCai.j}")
    return caidf


if __name__ == '__main__':
    parser = argparse.ArgumentParser(prog='biopython_CAI.py', description='Fetch transcript features (sequence and data) from ENSEMBL, for each gene ID in a list and output a multi FASTA file.', epilog="Author: Costas Bouyioukos, 2019, Paris UMR7216.")
    parser.add_argument('infile', nargs='?', default='-', type=argparse.FileType('r'), metavar="input_file", help='Path to the input FASTA file. (or STDIN).')
    parser.add_argument("outfile", nargs='?', default='-', type=argparse.FileType('w'), metavar='output_file', help="Path to output table file. (or STDOUT).")
    codon_counts = parser.add_mutually_exclusive_group(required=True)
    codon_counts.add_argument('-f', '--orfs', nargs = "?", help="Path to a FASTA file contianing all ORFs that we need to calculate CAI for a new organism.  In that case the program produces in the STDOUT a dictionary of codon counts. (Default=None).", type=argparse.FileType('r'), default=None, dest="orfs", metavar="ORFsFile")
    codon_counts.add_argument('-o', '--organism', nargs = "?", choices=['human', 'mouse', 'yeast'], help="Selection of the organism of study. In that case the program will calculate the CAI based on the hardcoded codon counts dictionaries. (Default=human, accepted values one of: ['human', 'mouse', 'yeast']).", type=str, default='human', dest="org", metavar="Organism")

    # Parse the command line arguments.
    optArgs = parser.parse_args()

    # Load or generate the codon counts.
    if optArgs.orfs:
        sys.stdout.write(str(generate_codon_dict_fasta(optArgs.orfs)) + "\n")
        sys.exit("Calculation of codon counts dictionary finished!")

    # Generate the CAIndex dictionary.
    generate_index_dict(optArgs.org)

    #print(CurrentCodonIndex)

    # Calculate the CAI pre gene.
    pdCai = calculate_CAI(optArgs.infile)

    # Print it to outfile.
    pdCai.to_csv(optArgs.outfile, sep = "\t")
