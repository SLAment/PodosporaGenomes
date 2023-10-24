#!/usr/bin/env python
# encoding: utf-8

import os # For the input name
import sys # To exit the script 
import re # Regular expressions

gtf = "data/Podans_v2016/genome_annotation_PODANS_SCV.gff"

functionre = re.compile("(ID=[\w:\_\d\.\s;=]*)(;function [\w\-_; #]*)")

for line in open(gtf, 'r'):
	cols = line.rstrip("\n").split("\t")		# break the line into columns defined by the tab
	attributes = cols[8]
	matchy = functionre.search(attributes)
	if matchy:
		cleanattributes = matchy.group(1)
		cleanattributes
		# print(cleanattributes)
	else:
		cleanattributes = attributes

	newcols = cols
	newcols[8] = cleanattributes.replace(" ;", ";").rstrip(" ;")
	newline = '\t'.join(newcols) + ";" + "\n"
	sys.stdout.write(newline)
