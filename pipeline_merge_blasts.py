#!/usr/bin/env python
from __future__ import division

import os
import argparse
import console
import fileinput

"""
	When running BLAST with GenBank (INSDC) data using seperate partitions, this script
	helps to merge different BLAST results together by sorting best hits using BLAST 
	score. Input is taken from STDIN, where you can define list of files and output is
	written to STDOUT.
"""

parser = argparse.ArgumentParser(description = """
	Merges BLAST outputs into one BLAST file by selecting best BLAST score hits only for 
	each sequence. Input is taken from STDIN and output is written to STDOUT.
	""")
args = parser.parse_args()

hits = {}

for line in fileinput.input():
	col = line.strip().split("\t")
	if len(col) < 15:
		continue
	seq = col[0]
	score = col[15]
	if seq in hits:
		comp_score = hits[seq][15]
		if score > comp_score:
			hits[seq] = col
	else:
		hits[seq] = col
		
for key in hits:
	print("\t".join(hits[key]))
