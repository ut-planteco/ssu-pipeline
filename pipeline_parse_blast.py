#!/usr/bin/env python
from __future__ import division

import os
import argparse
import console
import sys

"""
	Parsing XML BLAST output for tab delimited file for easier parsing with samples,
	hit descripttions, hit identities and hit aligment length values. Input is taken
	from STDIN, use it in pipe command.
"""

parser = argparse.ArgumentParser(description = """
	Parsing XML BLAST output for tab delimited file for easier parsing with samples,
	hit descripttions, hit identities and hit aligment length values. Input is taken
	from STDIN, use it in pipe command.
	""")
args = parser.parse_args()

i = 0

for line in sys.stdin:
	tmp = line.strip().split("<")
	if len(tmp) > 1:
		tmp2 = tmp[1].split(">")
		tag = tmp2[0]
		if len(tmp2) > 1:
			value = tmp2[1]
		else:
			value = None
		if tag == "Iteration_query-def":
			i += 1
			if i % 100 == 0:
				console.log("%d hits parsed\r" % (i))
			qry = {}
			qry['qry-id'] = value
		if tag == "Iteration_query-len":
			qry['qry-len'] = value
		if tag == "Hit_num":
			qry['hit'] = {}
			qry['hit']['nr'] = value
		if tag == "Hit_id":
			if value[:2] == "gi":
				qry['hit']['id'] = value.split("|")[1]
			else:
				qry['hit']['id'] = ""
		if tag == "Hit_accession" and qry['hit']['id'] == "":
			qry['hit']['id'] = value
		if tag == "Hit_def":
			qry['hit']['def'] = value.split("&gt;")[0]
		if tag == "Hit_len":
			qry['hit']['len'] = value
		if tag == "Hsp_num":
			qry['hit']['hsp'] = {}
			qry['hit']['hsp']['nr'] = value
		if tag == "Hsp_bit-score":
			qry['hit']['hsp']['score'] = value
		if tag == "Hsp_evalue":
			qry['hit']['hsp']['evalue'] = value
		if tag == "Hsp_query-from":
			qry['hit']['hsp']['qfrom'] = value
		if tag == "Hsp_query-to":
			qry['hit']['hsp']['qto'] = value
		if tag == "Hsp_hit-from":
			qry['hit']['hsp']['rfrom'] = value
		if tag == "Hsp_hit-to":
			qry['hit']['hsp']['rto'] = value
		if tag == "Hsp_identity":
			qry['hit']['hsp']['identity'] = value
		if tag == "Hsp_align-len":
			qry['hit']['hsp']['alen'] = value
		if tag == "Hsp_hit-frame":
			if value == "1":
				value = "+/+"
			else:
				value = "+/-"
			qry['hit']['hsp']['frame'] = value

		if tag == "Hsp_midline":
			# print our result
			identity = float(qry['hit']['hsp']['identity']) / float(qry['hit']['hsp']['alen']) * 100
			tmp = qry['qry-id'].split("-")
			if len(tmp) > 1:
				sample = tmp[0]
			else:
				sample = "NA"
			mlen = min(int(qry['qry-len']), int(qry['hit']['len']))
			alen = float(qry['hit']['hsp']['alen']) / float(mlen) * 100
			if alen > 100:
				alen = 100
			sys.stdout.write("\t".join([qry['qry-id'], qry['hit']['id'], qry['hit']['def'], qry['hit']['hsp']['evalue'], 
				"{0:.2f}".format(identity), qry['hit']['hsp']['identity'], 
				qry['hit']['hsp']['alen'], qry['hit']['hsp']['nr'], qry['hit']['hsp']['frame'],
				qry['hit']['hsp']['qfrom'], qry['hit']['hsp']['qto'], qry['hit']['hsp']['rfrom'],
				qry['hit']['hsp']['rto'], qry['qry-len'], qry['hit']['len'], qry['hit']['hsp']['score'], "{0:.2f}".format(alen), "\n"]))
console.log("%d hits parsed\n" % (i))
