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
	Summarizing INSDC BLAST results with user define identity and alignment length
	thresholds. 
	""")
parser.add_argument(
	'-b', metavar = 'BLAST_FILE', required = True, type = open, help = """
	BLAST tabulated output that was generated with pipeline_parseblast.py
	""")
parser.add_argument(
	'-ti', metavar = 'ID_FILE', required = True, type = open, help = """
	Taxonomy file, where for each GenBank ID node ID is specified
	""")
parser.add_argument(
	'-tt', metavar = 'TAXONOMY_FILE', required = True, type = open, help = """
	Taxonomy file, where for each node ID scientific name is provided
	""")
parser.add_argument(
	'-tn', metavar = 'NODE_FILE', required = True, type = open, help = """
	Taxonomy file, where full tree node connections of node IDs are provided to 
	build full taxonomy tree 
	""")
parser.add_argument(
	'-i', metavar = 'IDENTITY[0-100]', required = True, type = int, help = """
	hit identity in percentage to be accepted as a hit, recommended 90
	""")
parser.add_argument(
	'-l', metavar = 'ALIGNMENT[0-100]', required = True, type = int, help = """
	hit aliginment length in percentage to be accepted a hit, recommended 90
	""")
args = parser.parse_args()

rows = {}
cols = {}
hits = {}
ids = {}
nodes = {}
names = {}
taxonomy_tree = {}

total = 0
# for fasta sequences
sequences = {}
samples = {}
sequences_total = 0

def buildTree(tax_id):
	global taxonomy_tree
	if tax_id in taxonomy_tree:
		return taxonomy_tree[tax_id]
	elif tax_id in ids:
		convert_id = ids[tax_id]
		l = []
		while convert_id != "1":
			if convert_id in nodes:
				if convert_id in names:
					l.append(names[convert_id])
				else:
					l.append("root")
				convert_id = nodes[convert_id]
				if convert_id == "1":
					l.append("root")
			else:
				convert_id = "1"
		tmp = "; ".join(l[::-1])
		if len(tmp) == 0:
			tmp = tax_id
		taxonomy_tree[tax_id] = tmp
		return tmp
	else:
		taxonomy_tree[tax_id] = tax_id
		return taxonomy_tree[tax_id]

console.log("Fetching list of GenBank IDs from BLAST\n")
for f in args.b:
	f = f.strip()
	col = f.split("\t")
	ids[col[1]] = True

console.log("Fetching conversion of GenBank IDs to nodes\n")
for f in args.ti:
	f = f.strip()
	col = f.split("\t")
	if col[0] in ids:
		ids[col[0]] = col[1]

console.log("Building taxonomy tree lookup\n")
for f in args.tn:
	f = f.strip()
	col = f.split("\t|\t")
	nodes[col[0]] = col[1]

console.log("Building taxonomy name lookup\n")
for f in args.tt:
	f = f.replace("\t|\n", "\n").strip()
	col = f.split("\t|\t")
	if col[3] == "scientific name":
		names[col[0]] = col[1]

i = 0
console.log("Parsing BLAST\n")
with open(args.b.name) as infile:
	for f in infile:
		f = f.strip()
		col = f.split("\t")
		mlen = min(int(col[13]), int(col[14]))
		alen = int(col[6])
		i += 1
		if i % 10000 == 0:
			console.log("%d/%d hits found/parsed\r" % (total, i))
		if float(col[4]) >= args.i and alen >= mlen * (args.l / 100.0):
			sample = col[0].split("-")[0]
			hit = buildTree(col[1])
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