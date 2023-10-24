#!/usr/bin/env python
# encoding: utf-8

# ================== addtRNA2tbl =================

# Funannotate tool gff2tbl removes the tRNA genes from my annotation so this
# script puts them back into the tbl file.

# ==================================================
# Sandra Lorena Ament Velasquez
# 2023/07/02
# +++++++++++++++++++++++++++++++++++++++++++++++++

# ------------------------------------------------------
import os # For the input name
import sys # For reading the input
import argparse # For the fancy options
# ------------------------------------------------------

version = 1.0
versiondisplay = "{0:.2f}".format(version)

# ============================
# Check input file
# ============================
# Make a nice menu for the user
parser = argparse.ArgumentParser(description="* Append extra information to Podospora gff annotations *", epilog="BLAST must be locally installed.") # Create the object using class argparse

# Add options
parser.add_argument('tbl', help="tbl file")
parser.add_argument('gff', help="Gff3 file with the missing tRNA genes")

# extras
parser.add_argument('--version', '-v', action='version', version='%(prog)s ' + versiondisplay)

try:
	# ArgumentParser parses arguments through the parse_args() method You can
	# parse the command line by passing a sequence of argument strings to
	# parse_args(). By default, the arguments are taken from sys.argv[1:]
	args = parser.parse_args()
	tblopen = open(args.tbl, 'r')
	gffopen = open(args.gff, 'r')
except IOError as msg:  # Check that the file exists
    parser.error(str(msg)) 
    parser.print_help()
# ============================

# Make a function to disect the gene ID
def getGeneID(attributes):
	attris = attributes.split(";")
	for a in attris:
		if "Parent=" in a: # Assuming all gene models have an ID
			geneID = a.split("=")[1]
	return(geneID)

def getProduct(attributes):
	attris = attributes.split(";")
	for a in attris:
		if "product=" in a or "Product=" in a:
			product = a.split("=")[1]
	return(product)

# def getNote(attributes):
# 	attris = attributes.split(";")
# 	for a in attris:
# 		if "Note=" in a or "note=" in a:
# 			note = a.split("=")[1]
# 	return(note)	

# Make a dictionary for all contigs, each containing another dictionary with coordinates
# The big assumption here is that there are no genes with shared starting coordinates
gffDic = {}
for line in gffopen: # Without saving in memory
	if line.startswith('#') or (line == '\n'):
		sys.stdout.write(line)
	else:
		cols = line.rstrip("\n").split("\t")		# break the line into columns defined by the tab

		# My tRNA annotation has everything I need in the tRNA feature
		if cols[2] == "tRNA":
			contig = cols[0]
			geneID = getGeneID(cols[8])
			product = getProduct(cols[8])
			start = int(cols[3])
			end = int(cols[4])

			if 'Pseudo_Undet' in product:
				# notetrna = getNote(cols[8]).lstrip("probable Pseudo_Undet tRNA: anti-codon ") # Then the information of the anticodon is not in the product
				newgene = f"{start}\t{end}\tgene\n\t\t\tlocus_tag\t{geneID}\n\t\t\tpseudo\n{start}\t{end}\ttRNA\n\t\t\tproduct\ttRNA-Xxx\n\t\t\tpseudo\n"
			else:
				newgene = f"{start}\t{end}\tgene\n\t\t\tlocus_tag\t{geneID}\n{start}\t{end}\ttRNA\n\t\t\tproduct\t{product}\n"

			if contig not in gffDic.keys():
				gffDic[contig] = {start: (end, newgene)} # assuming there are no missing tRNA genes that start with the same coordinates
			else:
				gffDic[contig][start] = (end, newgene)

# Now read the tbl file
notthisgene = True
for line in tblopen:
	if '>Feature' in line:
		cols = line.rstrip("\n").split(" ")
		contig = cols[1]
		currentstart = 1
		currentend = 1
		sys.stdout.write(line)
	elif 'gene\n' in line:
		notthisgene = True
		cols = line.rstrip("\n").split("\t")		# break the line into columns defined by the tab

		# The tbl format expresses sense with the coordinates
		left = int(cols[0].lstrip("<"))
		right = int(cols[1].lstrip(">"))

		if left < right:
			startgene, endgene = left, right
		else:
			startgene, endgene = right, left

		# Not very efficient :/
		if contig in gffDic.keys(): # Some contigs won't have tRNAs
			for item in gffDic[contig].keys():
				if item > currentstart and item <= startgene: # The tRNA has to go before this gene
					if item == startgene and gffDic[contig][item][0] == endgene:
						notthisgene = False # This is the same tRNA gene but it's badly formatted so raplace it
					
					if "QC763_0108930" not in gffDic[contig][item][1]: # this gene overlaps Pa_7_9730 in CBS411.78m and NCBI doesn't like it, so just ignore it.
						sys.stdout.write(gffDic[contig][item][1])

		if notthisgene:
			sys.stdout.write(line) # It's not the same gene, so print the line

		# Assign this gene as the last point of reference
		currentstart = startgene
		currentend = endgene			
	elif notthisgene:
		sys.stdout.write(line)


