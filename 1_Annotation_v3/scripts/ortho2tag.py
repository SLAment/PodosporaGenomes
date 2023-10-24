#!/usr/bin/env python
# encoding: utf-8

# ================== ortho2tag =================

# I want to retain information of orthology to the P. anserina strain S+
# annotation (Podans_v2016) in the feauture IDs of the annotation of other
# strains from P. anserina species complex. The idea is to fuse it with the
# locus_tag provided by NCBI for some of the assemblies. This is what was 
# done for the PODCO annotation in NCBI.

# ==================================================
# Sandra Lorena Ament Velasquez
# 2023/05/05
# +++++++++++++++++++++++++++++++++++++++++++++++++

# ------------------------------------------------------
import os # For the input name
import sys # For reading the input
import argparse # For the fancy options
import gffutils
import re # Regular expressions
# ------------------------------------------------------
version = 1.0
versiondisplay = "{0:.2f}".format(version)

# ============================
# Check input file
# ============================
# Make a nice menu for the user
parser = argparse.ArgumentParser(description="* Append the names of a reference database of genes to a new annotation in another genome. *", epilog="BLAST must be locally installed.") # Create the object using class argparse

# Add options
parser.add_argument('gff3', help="GFF3 with names from GffgenesIDFix.py")
# parser.add_argument('ortholist', help="List of one-to-one orthologs produced by GffgenesIDFix.py")
parser.add_argument('--version', "-v", action='version', version='%(prog)s ' + versiondisplay)

try:
	# ArgumentParser parses arguments through the parse_args() method You can
	# parse the command line by passing a sequence of argument strings to
	# parse_args(). By default, the arguments are taken from sys.argv[1:]
	args = parser.parse_args()
	gffopen = open(args.gff3, 'r')
	# onetoonesopen = open(args.onetoones, 'r')
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

db = gffutils.create_db(data = args.gff3, 
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

# ---------------------------
### Read the onetoones file
# ---------------------------
# Get a list of the reference gene names and order 
# orthodict = {line.rstrip("\n").split("\t")[0]: line.rstrip("\n").split("\t")[1] for line in open(args.ortholist, 'r')}

def printWholeGene(gene):
	print(gene)
	for child in list(db.children(gene, order_by='start', level = 1)): 
		print(child)
	for grandchild in list(db.children(gene, order_by='start', level = 2)):
		print(grandchild)

# ---------------------------
### Change IDs
# ---------------------------
# Print headers
for line in gffopen:
	if line.startswith( '#' ):
		sys.stdout.write(line)
gffopen.close()

## Get iterator of all genes
genes = [gene for gene in db.features_of_type("gene")]
genenames = [feature.id for feature in db.all_features() if feature.featuretype in ["gene"]]

PaGenePattern = re.compile('(P[a-z]+)_([\d])_([\d]+)([\.\-\w\_\d]*)')
for gene in genes:
	print()
	name = gene.attributes['Name'][0] # Assuming there is only one
	matchy = PaGenePattern.search(name)
	if matchy:
		species = matchy.group(1)
		chromosome = matchy.group(2)
		genenumber = matchy.group(3)
		extras = matchy.group(4)

		if not extras: # This is a one-to-one ortholog
			tag, genenum = gene['ID'][0].split('_') # Get locus tag
			newidgene = f"{tag}_{chromosome}{'0'*(5-len(genenumber))}{genenumber}" # change to tag_XNNNN, where X is the chromosome number. E.g., PODCO_311360
			if newidgene in genenames: # Some of the original Podan2 genes were not placed, so had a Pa_0_NNNN notation
				newidgene = f"{tag}_9{'0'*(5-len(genenumber))}{genenumber}" # To avoid overlaps, there is no chromosome 9
			
			gene['ID'] = newidgene

			print(gene)

			# Fix all children
			mRNA = 1
			for child in list(db.children(gene)):
				if child.featuretype == "tRNA":
					child['Parent'] = newidgene
				
					# make a new ID for the mRNA
					newidrna = newidgene + "-tRNA"
					child['ID'] = newidrna
					child['Name'] = newidrna
					
					print(child)


				elif child.featuretype == "mRNA":
					child['Parent'] = newidgene
				
					# make a new ID for the mRNA
					newidrna = newidgene + "-T" + str(mRNA)
					mRNA += 1
					child['ID'] = newidrna
					child['Name'] = newidrna
					
					print(child)

					# Update all the features that depend on the mRNA, assuming each feature has a SINGLE parent
					typeids = {'exon':1, 'CDS':1, 'five_prime_UTR':1, 'three_prime_UTR':1, 'intron':1, 'repeat':1, 'sequence_feature': 1 }
					for grandchild in list(db.children(gene, level = 2)):
						grandchild['Parent'] = newidrna # TODO Add something here, a conditional, in case there are multiple parents

						typefea = grandchild.featuretype # What type are we dealing with?
						if typefea == 'five_prime_UTR': typename = '5pUTR'
						elif typefea == 'three_prime_UTR': typename = '3pUTR'
						else: typename = typefea
						
						newidchild = newidrna + '-' + typename + "{0:01d}".format(typeids[typefea])
						grandchild['ID'] = newidchild

						typeids[typefea] += 1 # Increase the count of the corresponding type for this gene				
						print(grandchild)

		else: # Not a one-to-one ortholog
			printWholeGene(gene)
	else:
		printWholeGene(gene)

