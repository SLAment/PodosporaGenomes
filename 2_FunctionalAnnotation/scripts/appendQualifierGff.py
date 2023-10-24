#!/usr/bin/env python
# encoding: utf-8

# ================== appendQualifierGff =================

# From a list of genes, append a qualifier to genes. The input is a list of
# gene IDs and the desired qualifiers. For example,
# the "pseudogene=unprocessed" qualifier, as preferred by NCBI

# https://www.ncbi.nlm.nih.gov/genbank/genomes_gff/
# https://www.insdc.org/submitting-standards/pseudogene-qualifier-vocabulary/

# unprocessed: the pseudogene has arisen from a copy of the parent gene by
# duplication followed by accumulation of random mutation. The changes,
# compared to their functional homolog, include insertions, deletions,
# premature stop codons, frameshifts and a higher proportion of
# non-synonymous versus synonymous substitutions.

# ==================================================
# Sandra Lorena Ament Velasquez
# 2023/06/30
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
parser.add_argument('gff3', help="GFF3 with names from GffgenesIDFix.py")
parser.add_argument('genelist', help="List of GeneIDs and their desired qualifiers (separated by a tab)")

# extras
parser.add_argument('--version', '-v', action='version', version='%(prog)s ' + versiondisplay)

try:
	# ArgumentParser parses arguments through the parse_args() method You can
	# parse the command line by passing a sequence of argument strings to
	# parse_args(). By default, the arguments are taken from sys.argv[1:]
	args = parser.parse_args()
	gffopen = open(args.gff3, 'r')
	listopen = open(args.genelist, 'r')
except IOError as msg:  # Check that the file exists
    parser.error(str(msg)) 
    parser.print_help()
# ============================

# Make dictionary with the genes and their desired qualifiers
qualifiers =  {}
for line in listopen:
	cols = line.rstrip("\n").split("\t")
	qualifiers[cols[0]] = cols[1]

# Make a function to disect the gene ID
def getGeneID(attributes):
	attris = attributes.split(";")
	for a in attris:
		if "ID=" in a: # Assuming all gene models have an ID
			geneID = a.split("=")[1]
			return(geneID)

# Append qualifier to the right genes
for line in gffopen: # Without saving in memory
	if line.startswith('#') or (line == '\n'):
		sys.stdout.write(line)
	else:
		cols = line.rstrip("\n").split("\t")		# break the line into columns defined by the tab

		if cols[2] == "gene":
			geneID = getGeneID(cols[8])
			if geneID in qualifiers.keys():
				newline = line.rstrip("\n").rstrip(";") + ";" + qualifiers[geneID] + ";\n" # Assuming it should end in ;
				sys.stdout.write(newline)
			else:
				sys.stdout.write(line)
		else:
			sys.stdout.write(line)
