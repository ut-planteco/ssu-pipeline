#!/usr/bin/env python
from __future__ import division

import os
import argparse
import gzip
import console
import sys

"""
	Prepare arguments to be parsed for cleaning Illumina sequences
"""

parser = argparse.ArgumentParser(description = """
	Cleaning of Illumina paired-end sequences by quality filtering and primer checking 
	both reads and outputing interleaved FASTQ file for FLASh program for combing the reads.
	""")
parser.add_argument(
	'-folder', metavar = 'FOLDER', required = True, type = None, help = """
	define FOLDER where FASTQ or FASTQ.tar.gz files are stored
	""")
parser.add_argument(
	'-fprimer', metavar = 'SEQUENCE', required = False, type = None, help = """
	define forward read primer
	""")
parser.add_argument(
	'-rprimer', metavar = 'SEQUENCE', required = False, type = None, help = """
	define reverse read primer
	""")
parser.add_argument(
	'-fadapter', metavar = 'SEQUENCE', required = False, type = None, help = """
	define adapter for forward read
	""")
parser.add_argument(
	'-radapter', metavar = 'SEQUENCE', required = False, type = None, help = """
	define adapter for reverse read
	""")
parser.add_argument(
	'-quality', metavar = 'QUALITY', default = 30, type = int, help = """
	average quality of sequence to be accepted
	""")
parser.add_argument(
	'-phred', metavar = 'QUALITY', default = 33, type = int, help = """
	FASTQ file phred quality score (33)
	""")
args = parser.parse_args()

def avgQuality(qual, phred):
	total = 0
	for i in qual:
		total += (ord(i) - phred)
	return float(total) / len(qual)

def startsWith(full, sub):
	if len(full) >= len(sub) and sub == full[0:len(sub)]:
		return True
	else:
		return False

files = [f for f in os.listdir(args.folder) if os.path.isfile(os.path.join(args.folder, f)) and "_R1_" in f]

delfiles = []

if args.folder[-1:] != "/":
	args.folder += "/"

seq_i = 0
seq_sel = 0

if args.fprimer is None:
	args.fprimer = ""
if args.rprimer is None:
	args.rprimer = ""

j = 0

for f in files:
	sample = f.split("_")[0]
	f = args.folder + f
	console.log("Parsing file %s\n" % (f))
	r = f.replace("_R1_", "_R2_")
	if not os.path.isfile(r):
		r = None
		console.log("Skipping file %s as it does not have other pair" % (f))
		continue
	fh1 = None
	fh2 = None
	if f[-3:] == ".gz":
		console.unpack(f, args.folder)
		f = f.replace(".tar.gz", "")
		delfiles.append(f)
		if r is not None:
			console.unpack(r, args.folder)
			r = r.replace(".tar.gz", "")
			delfiles.append(r)
	fh1 = open(f, "r+")
	if r is not None:
		fh2 = open(r, "r+")
	if fh1 is not None:
		c1 = []
		c2 = []
		i = 0
		with fh1 as r1, fh2 as r2:
			for l1, l2 in zip(r1, r2):
				i += 1
				j += 1
				c1.append(l1.strip())
				c2.append(l2.strip())
				if i % 4 == 0:
					seq_i += 1
					if(seq_i % 1000 == 0):
						console.log("%d/%d sequences selected/parsed\r" % (seq_sel, seq_i))
					if startsWith(c1[1], args.fprimer) and startsWith(c2[1], args.rprimer) and avgQuality(c1[3], args.phred) >= args.quality and avgQuality(c2[3], args.phred) >= args.quality:
						if args.fadapter is not None and args.fadapter in c1[1]:
							continue
						if args.radapter is not None and args.radapter in c2[1]:
							continue
						seq_sel += 1
						c1[0] = "@%s-%d" % (sample, j)
						c2[0] = "@%s-%d" % (sample, j)
						sys.stdout.write("%s\n%s\n" % ("\n".join(c1), "\n".join(c2)))
					c1 = []
					c2 = []

console.log("%d/%d sequences selected/parsed\n" % (seq_sel, seq_i))
console.log("Cleaning up filesystem\n")
for f in delfiles:
	os.remove(f)
