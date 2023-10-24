#!/usr/bin/env python
# encoding: utf-8

# ================== decoratePodoGff3 =================

# I want to retain information of orthology to the P. anserina strain S+
# annotation (Podans_v2016) in the feauture IDs of the annotation of other
# strains from P. anserina species complex. The idea is to fuse it with the
# locus_tag provided by NCBI for some of the assemblies. This is what was 
# done for the PODCO annotation in NCBI.

# ==================================================
# Sandra Lorena Ament Velasquez
# 2023/05/14
# +++++++++++++++++++++++++++++++++++++++++++++++++

# ------------------------------------------------------
import os # For the input name
import sys # For reading the input
# print(sys.version)
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
parser = argparse.ArgumentParser(description="* Append extra information to Podospora gff annotations *", epilog="BLAST must be locally installed.") # Create the object using class argparse

# Add options
parser.add_argument('gff3', help="GFF3 with names from GffgenesIDFix.py")
parser.add_argument('ortholist', help="List of one-to-one orthologs produced by GffgenesIDFix.py")
parser.add_argument('--sppcode', '-s', help="Species code. Default 'Pa'.", default="Pa")
parser.add_argument('--namedgenes', '-g', help="Table of genes codes with named gene products from Podospora anserina", type=str) 
parser.add_argument('--version', "-v", action='version', version='%(prog)s ' + versiondisplay)

try:
	# ArgumentParser parses arguments through the parse_args() method You can
	# parse the command line by passing a sequence of argument strings to
	# parse_args(). By default, the arguments are taken from sys.argv[1:]
	args = parser.parse_args()
	gffopen = open(args.gff3, 'r')
except IOError as msg:  # Check that the file exists
    parser.error(str(msg)) 
    parser.print_help()
# ============================

if args.namedgenes:
	products = open(args.namedgenes, 'r')
	proddic = {}
	for line in products:
		cols = line.rstrip("\n").split("\t")
		proddic[cols[0]] = (cols[1], cols[2])

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

# ortho2tag.py already change the gene IDs so I can't use that information...
# orthodict = {}
# for line in open(args.ortholist, 'r'):
# 	geneid, genealias = line.rstrip("\n").split("\t")
# 	orthodict[geneid] = genealias

# So get a list of the orthologs that got assigned  
orthos = [line.rstrip("\n").split("\t")[1] for line in open(args.ortholist, 'r')]

## Get iterator of all genes
genes = [gene for gene in db.features_of_type("gene")]

# # Get ids of genes in database
# allIDs = [feature.id for feature in db.all_features() if feature.featuretype in ["gene"]]

def printWholeGene(gene):
	print(gene)
	for child in list(db.children(gene, order_by='start', level = 1)): 
		print(child)

	cdscount = 1
	for grandchild in list(db.children(gene, order_by='start', level = 2)):

		# Only necessary because funannotate makes the cds IDs not-unique
		if grandchild.featuretype == "CDS":
			grandchild['ID'] = grandchild['ID'][0] + str(cdscount)
			cdscount += 1
		
		print(grandchild)

# ---------------------------
### Change IDs
# ---------------------------
# Print headers
for line in gffopen:
	if line.startswith( '#' ):
		sys.stdout.write(line)
gffopen.close()

for gene in genes:
	thisgeneID = gene['ID'][0]
	# Infer the name of the one-to-one ortholog (appended by ortho2tag.py)
	tagpattern = re.compile("([\w]*)_(\d)(\d*)")
	matchy = tagpattern.search(thisgeneID)
	if matchy: # It should always be the case
		tag = matchy.group(1)
		chromosome = int(matchy.group(2))
		genecode = matchy.group(3)

		if chromosome == 9: 
			chromosome = 0 # I had to do that trick for a few orthologs in ortho2tag.py
		
		orthocode = f"{args.sppcode}_{chromosome}_{genecode.lstrip('0')}"
		thisproduct = f"Pa_{chromosome}_{genecode.lstrip('0')}"

		# Genes with one-to-one orthologs to P. anserina
		if (orthocode in orthos and len(genecode) + 1 == 6) or (thisproduct in proddic.keys()): 
			# Verify it's in the orthos list but also, if the length of the genecode+1 is 7, then it's just a look-alike, but actually not an ortholog
			# This is a consequence of my choice of length for the IDs with GFFnumerator.py and then by the renaming with ortho2tag.py
			# Finally, if it didn't make it to the onetoone file, I might have added it later during manual curation

			# gene['Alias'] = orthocode
			gene['Note'] = orthocode
			gene['color'] = "#000000" # Add a color for visualization in IGV

			# Add the known products
			samegene = False
			if args.namedgenes:	
				if thisproduct in proddic.keys():
					gene['color'] = "#3399ff" # this is a special gene

					# Multiple names can be specified by providing the primary
					# name first, and additional names as a comma-separated
					# list. https://www.ncbi.nlm.nih.gov/genbank/genomes_gff/#attributes
					# print(gene)
					if 'Name' in gene.attributes:
						for name in gene['Name']:
							if name.lower() == proddic[thisproduct][0].lower().replace('-', ''): 
								samegene = True # e.g. HETQ1 and Het-Q1
							elif 'HETE1' in name: 
								samegene = True # Ignore it later on and use my name
						if samegene:
							gene['Name'] = proddic[thisproduct][0] # Use mine
						else: 
							gene['Name'].append(proddic[thisproduct][0])

					else:
						gene['Name'] = proddic[thisproduct][0]

			# Print gene and children
			print(gene)
			for child in list(db.children(gene, order_by='start', level = 1)): 
				child['color'] = "#000000"
				if args.namedgenes and thisproduct in proddic.keys():
					child['color'] = "#3399ff"
					if not samegene and 'Name' in child.attributes: 
						# print(gene.attributes, file=sys.stderr)
						child['Name'].append(proddic[thisproduct][0])
						child['product'].append(proddic[thisproduct][1])
					else:
						child['Name'] = proddic[thisproduct][0]
						child['product'] = proddic[thisproduct][1]
				print(child)

			cdscount = 1
			for grandchild in list(db.children(gene, order_by='start', level = 2)):
				grandchild['color'] = "#000000"

				# Only necessary because funannotate makes the cds IDs not-unique
				if grandchild.featuretype == "CDS":
					grandchild['ID'] = grandchild['ID'][0] + str(cdscount)
					cdscount += 1

				print(grandchild)

		# Some genes I manually curated
		elif thisgeneID in proddic.keys(): 
			gene['Name'] = proddic[thisgeneID][0] # I don't care what funannotate had to say
			gene['color'] = "#279583" # Add a color for visualization in IGV

			print(gene)

			for child in list(db.children(gene, order_by='start', level = 1)): 
				gene['Name'] = proddic[thisgeneID][0]
				child['product'] = proddic[thisgeneID][1] # I don't care what funannotate had to say
				child['color'] = "#279583"
				print(child)

			cdscount = 1
			for grandchild in list(db.children(gene, order_by='start', level = 2)):
				grandchild['color'] = "#279583"

				# Only necessary because funannotate makes the cds IDs not-unique
				if grandchild.featuretype == "CDS":
					grandchild['ID'] = grandchild['ID'][0] + str(cdscount)
					cdscount += 1
				print(grandchild)

		else: # Genes that match the Pa_X_XXX codes but are not onetoones and are not in the manually-curated set: print as they are
			if 'Name' in gene.attributes and 'HETE1' in gene['Name'][0]: # It's an HNWD paralog, probably a pseudogene I didn't annotate and funannotate is giving it a misleading name
				gene['Name'] = 'HNWD'
			printWholeGene(gene)

	else: # This should never happen after using ortho2tag.py, but just in case
		print(f"Gene with ID {gene['ID']} is not following the expected format of TAG_NNNNN")
		sys.exit(1)



