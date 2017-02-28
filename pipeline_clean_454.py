#!/usr/bin/env python
from __future__ import division

import os
import argparse
import console

"""
	Prepare arguments to be parsed for cleaning 454 sequences
"""
parser = argparse.ArgumentParser(description = """ 
	Demultiplex 454 sequences by barcode and primer, filtering out low quality 
	sequences using average quality threshold.
	""")
parser.add_argument(
	'-f', metavar = 'FASTA_FILE', required = True, type = open, help = """
	FASTA file, where sequences are stored
	""")
parser.add_argument(
	'-qf', metavar = 'QUALITY_FILE', type = open, help = """
	QUALITY file, where sequence qualities are stored for each nucleotide. This
	file headers should match FASTA file headers.
	""")
parser.add_argument(
	'-b', metavar = 'BARCODE_FILE', required = True, type = open, help = """
	BARCODE file, where sample ID, barcodes and primers are stored in tabular 
	file format
	""")
parser.add_argument(
	'-bs', metavar = 'SAMPLE_COLUMN', required = True, type = int, help = """
	sample column in BARCODE file (first column is 1)
	""")
parser.add_argument(
	'-bb', metavar = 'BARCODE_COLUMN', required = True, type = int, help = """
	barcode column in BARCODE file
	""")
parser.add_argument(
	'-bp', metavar = 'PRIMER_COLUMN', type = int, help = """
	primer column in BARCODE file
	""")
parser.add_argument(
	'-q', metavar = 'AVERAGE_QUALITY', type = float, default = 25, help = """
	lower limit of average quality of the sequence to be filtered out 
	(%(default)s)
	""")
parser.add_argument(
	'-trimq', metavar = 'TRIM_QUALITY', type = float, help = """
	lower limit of average quality that is trimmed away when it drops below the 
	threshold (recommended = 20)
	""")
parser.add_argument(
	'-trimw', metavar = 'TRIM_WINDOW', type = int, help = """
	window size to calculate average quality for trimming sequence end (recommended = 50)
	""")
parser.add_argument(
	'-ml', metavar = 'MINIMUM_LENGTH', type = int, default = 170, help = """
	minimum allowed length of the sequences after trimming (%(default)s)
	""")
parser.add_argument(
	'-tl', metavar = 'TRUNCATE_LENGTH', type = int, help = """
	truncate sequences longer than provided length to remove reverse primer
	""")
args = parser.parse_args()

def calculateAverageQuality(quality):
	if len(quality) > 0:
		return sum(quality) / len(quality)
	else:
		return 0

def calculateSlidingWindow(quality, threshold, window):
	total = 0
	for i in range(len(quality) - window):
		avg = sum(quality[i:(i + window)]) / float(window)
		if avg < threshold:
			return i
	return len(quality)

def parseSequence(quality, args, lookup, sequence, f):
	avg_quality = 40
	if len(quality) > 0:
		if args.trimq is not None and args.trimw is not None:
			trim_len = calculateSlidingWindow(quality, args.trimq, args.trimw)
			if trim_len < len(quality):
				quality = quality[:trim_len]
				sequence = sequence[:trim_len]
		avg_quality = calculateAverageQuality(quality)
	for k in lookup:
		_length = len(sequence)
		if sequence.startswith(k) and avg_quality >= args.q \
				and len(sequence) >= args.ml + len(k):
			sequence = sequence[len(k):]
			_tmp = lookup[k]
			if args.tl is not None:
				sequence[0:args.tl]
			f.write(">%s-%s\t%s\t%s\t%s\t%.2f\t%d\t%d\n%s\n" % (_tmp[0], header[1:], _tmp[0], _tmp[1], _tmp[2], avg_quality, _length, len(sequence), sequence))
			return True
	return False
"""
	Read BARCODES into memory
"""
lookup = {}
for r in args.b:
	col = r.strip().split("\t")
	_sample = col[args.bs - 1].replace(" ", "-")
	_barcode = col[args.bb - 1]
	_primer = ""
	if args.bp is not None:
		_primer = col[args.bp - 1]
	index = _barcode + _primer
	lookup[index] = [_sample, _barcode, _primer]

"""
	Read FASTA and QUALITY file simultaneously and write output
"""
f = open(args.f.name.replace(".fna", "").replace(".fasta", "") + '.cleaned.fasta', 'w+')	# write results into file

header = ""
sequence = ""
quality = []
seq_i = 0
seq_c = 0
for r1 in args.f:
	r1 = r1.strip()
	r2 = None
	if args.qf is not None:
		r2 = args.qf.readline().strip()

	if r1.startswith(">"):
		seq_i += 1
		if seq_i % 1000 == 0:
			console.log("%d/%d sequences cleaned/parsed\r" % (seq_c, seq_i))
			
		if r1.split(" ")[0] != r2.split(" ")[0] and r2 is not None:
			raise ValueError('Fasta and quality file headers do not match')
			break
		if parseSequence(quality, args, lookup, sequence, f):
			seq_c += 1
		header = r1.split(" ")[0].split("\t")[0]
		sequence = ""
		quality = []
	else:
		sequence += r1
		if r2 is not None:
			# add together average quality
			col = r2.split(" ")
			for i in col:
				quality.append(int(i))
# add final sequence also
if parseSequence(quality, args, lookup, sequence, f):
	seq_c += 1
console.log("%d/%d sequences cleaned/parsed\n" % (seq_c, seq_i))
if f is not None:
	f.close()
