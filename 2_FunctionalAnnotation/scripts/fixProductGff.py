#!/usr/bin/env python
# encoding: utf-8

# ================== fixProductGff =================

# The product names of Funannotate are not always acceptable or correct. Here
# I fix a few of them.
# The sister script for tbl files is `tbl4ncbi.py`

# ==================================================
# Sandra Lorena Ament Velasquez
# 2023/07/04
# +++++++++++++++++++++++++++++++++++++++++++++++++

# ------------------------------------------------------
import os # For the input name
import sys # For reading the input
import argparse # For the fancy options
import gffutils
# ------------------------------------------------------

version = 1.0
versiondisplay = "{0:.2f}".format(version)

# ============================
# Check input file
# ============================
# Make a nice menu for the user
parser = argparse.ArgumentParser(description="* Append extra information to Podospora gff annotations *", epilog="BLAST must be locally installed.") # Create the object using class argparse

# Add options
parser.add_argument('gff', help="Gff3 file")
parser.add_argument('--genelist', '-g', help="List of GeneIDs and their desired qualifiers with description (separated by a tab)")

# extras
parser.add_argument('--version', '-v', action='version', version='%(prog)s ' + versiondisplay)

try:
	# ArgumentParser parses arguments through the parse_args() method You can
	# parse the command line by passing a sequence of argument strings to
	# parse_args(). By default, the arguments are taken from sys.argv[1:]
	args = parser.parse_args()
	gffopen = open(args.gff, 'r')
except IOError as msg:  # Check that the file exists
    parser.error(str(msg)) 
    parser.print_help()
# ============================


# ---------------------------------
# Make database
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

db = gffutils.create_db(data = args.gff, 
	keep_order = True,
	dbfn = dbfnchoice,
	force = True, # force=True overwrite any existing databases.
	id_spec = id_spec, 
	merge_strategy = "create_unique") # Add an underscore an integer at the end for each consecutive occurrence of the same ID 
	# verbose = True,) 

# t1 = time.time()
# db_results = inspect.inspect(db) # Report
# print("\n\nIt took {0:.1f}s to create database".format(t1 - t0))
# ---------------------------------

## Get iterator of all genes
genes = [gene for gene in db.features_of_type("gene")]

# ---------------------------
# Make a dictionary with the pseudogenes list
qualifiers =  {}
if args.genelist != None:
	listopen = open(args.genelist, 'r')
	for line in listopen:
		cols = line.rstrip("\n").split("\t")
		qualifiers[cols[0]] = (cols[1], cols[2])

# ---------------------------
### Fix gff
# ---------------------------
# Print headers
for line in gffopen:
	if line.startswith( '#' ):
		sys.stdout.write(line)
gffopen.close()

for gene in genes:
	thisgeneID = gene['ID'][0]

	# Use the information in the dictionary to add pseudo and pseudogene qualifiers
	if thisgeneID in qualifiers.keys(): 
		qtype, quality = qualifiers[thisgeneID][0].split('=')
		qnote = qualifiers[thisgeneID][1]

		# The qualifier about the type of pesudogene might already be there and will be replaced
		gene[qtype] = quality 

		# Add the note too do not replace what was there before.
		if 'note' in gene.attributes:
			gene['note'].append(qnote)
		else:
			gene['note'] = qnote

		# Print whole gene
		print(gene)
		for child in list(db.children(gene, order_by='start', level = 1)): 
			print(child)
		for grandchild in list(db.children(gene, order_by='start', level = 2)):
			print(grandchild)

	# Fix products with typos
	else:
		print(gene)
		for child in list(db.children(gene, order_by='start', level = 1)): # The products should be in the mRNAs (funannotate)
			if 'product' in child.attributes:
				# Full thing replacement
				if child['product'][0] == 'conserved fungal protein': child['product'][0] = 'hypothetical protein'
				elif child['product'][0] == 'Belongs to the peptidase S8': child['product'][0] = 'Mycosin-4'
				elif child['product'][0] == 'proteinral transcription repressor': child['product'][0] = 'General transcriptional corepressor TUP1' # https://www.uniprot.org/uniprotkb/P16649/entry
				elif child['product'][0] == 'proteinral amino acid permease agp2': child['product'][0] = 'General amino acid permease AGP2' # https://www.uniprot.org/uniprotkb/P38090/entry
				elif child['product'][0] == 'proteinral negative regulator of transcription subunit 5': child['product'][0] = 'General negative regulator of transcription subunit 5' # https://www.uniprot.org/uniprotkb/Q12514/entry
				elif child['product'][0] == 'transcriptional repressor proteinral negative regulator of transcription subunit 4': child['product'][0] = 'General negative regulator of transcription subunit 4' # https://www.uniprot.org/uniprotkb/Q12514/entry
				# Replace words
				else:
					if 'bioproteinsis' in child['product'][0]: child['product'][0] = child['product'][0].replace("bioproteinsis", "biogenesis")
					if 'diacyglycerol' in child['product'][0]: child['product'][0] = child['product'][0].replace("diacyglycerol", "diacylglycerol")
					if 'morphoproteinsis' in child['product'][0]: child['product'][0] = child['product'][0].replace("morphoproteinsis", "morphogenesis")
					if 'phopho' in child['product'][0]: child['product'][0] = child['product'][0].replace("phopho", "phospho")
					if 'redicted' in child['product'][0]: child['product'][0] = child['product'][0].replace("redicted", "utative")
					if 'robable' in child['product'][0]: child['product'][0] = child['product'][0].replace("robable", "utative")
					if 'human' in child['product'][0]: child['product'][0] = child['product'][0].replace("human", "")
			print(child)

		for grandchild in list(db.children(gene, order_by='start', level = 2)):
			print(grandchild)

# https://www.ncbi.nlm.nih.gov/genbank/genomes_gff/

# [2] pseudogenes should be flagged with pseudogene=<TYPE> qualifier in column
# 9 on the gene feature and optionally on any child features. Further details
# about the TYPE values allowed for the pseudogene qualifier are available
# at: http://www.insdc.org/documents/pseudogene-qualifier-vocabulary .

# [3] annotate with pseudo=true any genes that are 'broken' but are not
# thought to be pseudogenes. These are genes that do not encode the expected
# translation, for example because of internal stop codons. These are often
# caused by problems with the sequence and/or assembly.

