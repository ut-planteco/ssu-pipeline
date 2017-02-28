#!/usr/bin/env python
from __future__ import division

import os
import argparse
import console
import sys

"""
	Prepare arguments to be parsed for reverse complemeting sequences with 
	wrong direction.
"""
parser = argparse.ArgumentParser(description = """ 
	Reverse complements FASTA sequences that BLAST+ has identified with 
	+/- strand against MaarjAM database.
	""")
parser.add_argument(
	'-f', metavar = 'FASTA_FILE', required = True, type = open, help = """
	FASTA file, where sequences are stored
	""")
parser.add_argument(
	'-b', metavar = 'BLAST_FILE', required = True, type = open, help = """
	BLAST file, where strand is stored
	""")

args = parser.parse_args()

def reverse_complement(seq):   
	complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'} 
	bases = list(seq) 
	bases = reversed([complement.get(base,base) for base in bases])
	return ''.join(bases)

"""
	Read strand info into memory
"""
lookup = {}
for r in args.b:
	col = r.strip().split("\t")
	if col[8] is "+/-":
		lookup[col[0]] = True
"""
	Read FASTA file and reverse complement sequences with wrong strand
	Read FASTA and QUALITY file simultaneously and write output
"""
header = ""
sequence = ""

for r in args.f:
	r = r.strip()
	if r.startswith(">"):
		if len(sequence) > 0:
			if header in lookup:
				sequence = reverse_complement(sequence)
			sys.stdout.write("%s\n%s\n" % (header, sequence))
		header = r
	else:
		sequence += r

if len(sequence) > 0:
	if header in lookup:
		sequence = reverse_complement(sequence)
	sys.stdout.write("%s\n%s\n" % (header, sequence))
