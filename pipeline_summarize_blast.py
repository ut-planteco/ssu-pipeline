#!/usr/bin/env python
from __future__ import division

import os
import argparse
import console
import sys
import operator

"""
	Prepare arguments to be parsed for summarizing BLAST
"""
parser = argparse.ArgumentParser(description = """ 
	Summarizing BLAST results with user define identity and alignment length
	thresholds. Providing FASTA file will result nohit file that can be used 
	to run additional BLAST with another database.
	""")
parser.add_argument(
	'-b', metavar = 'BLAST_FILE', required = True, type = file, help = """
	BLAST tabulated output that was generated with pipeline_parseblast.py
	""")
parser.add_argument(
	'-f', metavar = 'FASTA_FILE', type = file, help = """
	FASTA file to be used to output list of bad hits that did not match
	thresholds
	""")
parser.add_argument(
	'-i', metavar = 'IDENTITY[0-100]', required = True, type = int, help = """
	hit identity in percentage to be accepted as a hit, recommended 97
	""")
parser.add_argument(
	'-l', metavar = 'ALIGNMENT[0-100]', required = True, type = int, help = """
	hit aliginment length in percentage to be accepted a hit, recommended 95
	""")
parser.add_argument(
	'-vs', metavar = 'VARIABLE_START', required = False, type = int, help = """
	reference sequence variable region start
	""")
parser.add_argument(
	'-ve', metavar = 'VARIABLE_END', required = False, type = int, help = """
	reference sequence variable region end
	""")
parser.add_argument(
	'-t', metavar = 'BLAST_TYPE[0-2]', required = True, type = int, help = """
	defines which section of the BLAST to be used to summarize results.
	0 - suitable for MaarjAM, only last portion of hit description is used,
	1 - all hit description is used,
	2 - hit identificator is used
	""")
args = parser.parse_args()

nohits = []
rows = {}
cols = {}
hits = {}
total = 0
# for fasta sequences
sequences = {}
samples = {}
sequences_total = 0

i = 0
console.log("Parsing BLAST\n")
for f in args.b:
	f = f.strip()
	col = f.split("\t")
	mlen = min(int(col[13]), int(col[14]))
	alen = int(col[6])
	i += 1
	if i % 10000 == 0:
		console.log("%d/%d hits found/parsed\r" % (total, i))
	if float(col[4]) >= args.i and alen >= mlen * (args.l / 100.0):
		if args.vs is not None and args.ve is not None and args.vs < args.ve:
			if col[11] > col[12]:
				col[11], col[12] = col[12], col[11]
			if args.vs < int(col[11]) or args.ve > int(col[12]):
				continue
		sample = col[0].split("-")[0]
		if args.t == 0:
			hit = col[2].split(" ")[-1]
		elif args.t == 1:
			hit = col[2].split("| ")[-1]
		else:
			hit = col[1]
		if hit in rows:
			rows[hit] += 1
		else:
			rows[hit] = 1
		if sample in cols:
			cols[sample] += 1
		else:
			cols[sample] = 1
		index = "_".join([hit, sample])
		if index in hits:
			hits[index] += 1
		else:
			hits[index] = 1
		total += 1
		sequences[col[0]] = True
	
console.log("%d/%d hits found/parsed\n" % (total, i))	
i = 0
if args.f is not None:
	console.log
	fh = open("%s.i%d.a%d.nohits.fasta" % (args.b.name, args.i, args.l), "w+")
	write = False
	for f in args.f:
		if len(f) > 0 and f[0] == ">":
			sequences_total += 1
			if sequences_total % 10000 == 0:
				console.log("%d/%d nohits sequences found/parsed\r" % (i, sequences_total))
			index = f[1:].split("-")[0]
			if index in samples:
				samples[index] += 1
			else:
				samples[index] = 1
			if f[1:] not in sequences:
				i += 1
				write = True
			else:
				write = False
		if write:
			fh.write(f)
	if fh is not None:
		fh.close()
	console.log("%d/%d nohits sequences found/parsed\n" % (i, sequences_total))	

sorted_rows = sorted(rows.items(), key=operator.itemgetter(1), reverse=True)
sorted_cols = sorted(cols.items(), key=operator.itemgetter(1), reverse=True)

console.log("Writing out pivot table with results\n")
fh = open("%s.i%d.a%d.tsv" % (args.b.name, args.i, args.l), "w+")
### print header of the tsv file
fh.write("Hit/Sample\tSamples\tTotal\t")
### go through samples as columns
for k in sorted_cols:
	fh.write("%s\t" % (k[0]))
fh.write("\nTaxa\t\t%d\t" % (len(rows)))
for k in sorted_cols:
	cnt = 0
	for _k in sorted_rows:
		index = "_".join([_k[0], k[0]])
		if index in hits:
			cnt += 1
	fh.write("%d\t" % (cnt))
# only if fasta file is defined
if sequences_total > 0:
	fh.write("\nCleaned reads\t\t%d\t" % (sequences_total))
	for k in sorted_cols:
		if k[0] in samples:
			fh.write("%d\t" % (samples[k[0]]))
		else:
			fh.write("\t")
fh.write("\nTotal\t%d\t%d\t" % (len(cols), total))
for k in sorted_cols:
	fh.write("%d\t" % (k[1]))
for k in sorted_rows:
	cnt = 0
	for _k in sorted_cols:
		index = "_".join([k[0], _k[0]])
		if index in hits:
			cnt += 1
	fh.write("\n%s\t%d\t%d\t" % (k[0], cnt, k[1]))
	for _k in sorted_cols:
		index = "_".join([k[0], _k[0]])
		if index in hits:
			fh.write("%d\t" % (hits[index]))
		else:
			fh.write("\t")
fh.write("\n")
if fh is not None:
	fh.close()
