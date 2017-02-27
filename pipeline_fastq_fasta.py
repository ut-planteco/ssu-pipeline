#!/usr/bin/env python
from __future__ import division

import os
import argparse
import console
import fileinput

"""
	Converts FASTQ file to FASTA. Input is taken from STDIN and outpt is written to STDOUT
"""

parser = argparse.ArgumentParser(description = """
	Converts FASTQ file to FASTA. Input is taken from STDIN and output is written to STDOUT.
	""")
args = parser.parse_args()

i = 0
for line in fileinput.input():
	i += 1
	if i % 4 == 1:
		print(line.replace("@", ">").strip())
	elif i % 4 == 2:
		print(line.strip())
    
