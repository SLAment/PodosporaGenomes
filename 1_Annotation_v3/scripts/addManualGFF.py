#!/usr/bin/env python
# encoding: utf-8

# ================== addManualGFF =================
# Rudimentary script to replace gene models in a focal gff overlapping with
# those present in a manually curated gff.
# Working only for genes at the moment!
# ==================================================
# Sandra Lorena Ament Velasquez
# 2023/06/02
# +++++++++++++++++++++++++++++++++++++++++++++++++

# ------------------------------------------------------
import sys # To exit the script 
import os # For the input name
import argparse  # For the fancy options
import gffutils
import re
# ------------------------------------------------------

version = 1.00
versiondisplay = "{0:.2f}".format(version)

# ============================
# Make a nice menu for the user
# ============================
parser = argparse.ArgumentParser(description="* Script to replace gene models in a focal GFF3 overlapping with those present in a manually curated GFF3*")  # Create the object using class argparse

# Add options
parser.add_argument('GFFfocal', help="Focal GFF3 file")
parser.add_argument('GFFmanual', help="Manually curated GFF3 file")

# extras
parser.add_argument('--version', '-v', action='version', version='%(prog)s ' + versiondisplay)

try:
	# ArgumentParser parses arguments through the parse_args() method You can
	# parse the command line by passing a sequence of argument strings to
	# parse_args(). By default, the arguments are taken from sys.argv[1:]
	args = parser.parse_args()
except IOError as msg:  # Check that the file exists
	parser.error(str(msg)) 
	parser.print_help()

# ---------------------------------
# Make databases
# ---------------------------------
# t0 = time.time()
# This will parse the file, infer the relationships among the features in the file, and store the features and relationships
# See https://pythonhosted.org/gffutils/autodocs/gffutils.create_db.html
# I expect a MAKER gff, where CDS have no unique ID

dbfnchoice = ':memory:'

# http://daler.github.io/gffutils/database-ids.html
id_spec={"gene": ["ID", "Name"], 
	"mRNA": ["ID", "transcript_id"], 
	"repeat": ["repeat", "simple_repeat", "similarity"]
	}

first_db = gffutils.create_db(data = args.GFFfocal, 
	keep_order = True,
	dbfn = dbfnchoice,
	# force = True, # force=True overwrite any existing databases.
	id_spec = id_spec, 
	merge_strategy = "create_unique") # Add an underscore an integer at the end for each consecutive occurrence of the same ID 
	# verbose = True,) 

second_db = gffutils.create_db(data = args.GFFmanual, 
	keep_order = True,
	dbfn = dbfnchoice,
	# force = True, # force=True overwrite any existing databases.
	id_spec = id_spec, 
	merge_strategy = "create_unique") # Add an underscore an integer at the end for each consecutive occurrence of the same ID 
	# verbose = True,) 

# t1 = time.time()
# db_results = inspect.inspect(db) # Report
# print("\n\nIt took {0:.1f}s to create database".format(t1 - t0))
# ---------------------------------
# ---------------------------------
# Define functions
# ---------------------------------
def printWholeGene(gene, db):
	print(gene)
	for child in list(db.children(gene, order_by='start', level = 1)): 
		print(child)
	for grandchild in list(db.children(gene, order_by='start', level = 2)):
		print(grandchild)

# ---------------------------------
# Process GFFs
# ---------------------------------
## Get iterator of all genes
genes_1 = [gene for gene in first_db.features_of_type("gene")]
genes_2 = [gene for gene in second_db.features_of_type("gene")]

# Presumably the manually curated gff is smaller so faster to handle it
# Make a dictionary containing the ranges of the protein coding part of manually curated genes
# To try to speed up things if the second gff is too big, I'll make a larger dictionary per contig to contain the ranges of each gene
contigdic = {}
ranges_manual = {}
for gene in genes_2:
	geneID = gene['ID'][0]
	# The overlap is only a problem if the CDS overlap
	posranges = []
	for child in list(second_db.children(gene, order_by='start', featuretype='CDS')): 
		posranges.append(child.start)
		posranges.append(child.end)
	# ranges_manual[geneID] = (min(posranges), max(posranges))

	if gene.seqid not in contigdic.keys():
		contigdic[gene.seqid] = {geneID: (min(posranges), max(posranges))}
	else:
		# print(gene, posranges)
		contigdic[gene.seqid].update({geneID: (min(posranges), max(posranges))})

# Print a header
print("##gff-version 3")

# Find the genes overlapping in the two GFFs
for gene in genes_1:
	if gene.seqid not in contigdic.keys(): # There are no annotated genes in this contig, so just print them as they are
		printWholeGene(gene, first_db)
	else:
		overlapgenes = []
		# Loop through the manually curated genes
		for manualgene in contigdic[gene.seqid].keys():
			start_focal = gene.start
			end_focal = gene.end

			start_manual = contigdic[gene.seqid][manualgene][0]
			end_manual = contigdic[gene.seqid][manualgene][1]

			# Case 1: the focal gene is upstream and it doesn't overlap at all with the manually-curated gene
			# Case 2: the focal gene is downstream and it doesn't overlap at all with the manually-curated gene
			if not (end_focal < start_manual or start_focal > end_manual):
				overlapgenes.append(gene['ID'][0])

			# Acutally I don't have to print the other cases, I can just print the
			# manual curated genes at the end and then do the sorting externally.

			# # Case 3: the focal gene overlaps with the N' terminus of the manually-curated gene
			# # Case 6: the focal gene contains the whole manually-curated gene
			# elif end_focal > start_manual and start_focal <= start_manual:
			# 	print(gene)
			# # Case 4: the focal gene overlaps with the C' terminus of the manually-curated gene
			# elif start_focal >= start_manual and start_focal <= end_manual:
			# 	print(gene)
			# # Case 5: the focal gene is contained within the manually-curated gene
			# elif start_focal >= start_manual and end_focal <= end_manual:
			# 	print(gene)

		if len(overlapgenes) == 0:
			printWholeGene(gene, first_db)

# Print the manually-curated genes at the end
for gene in genes_2:
	printWholeGene(gene, second_db)








