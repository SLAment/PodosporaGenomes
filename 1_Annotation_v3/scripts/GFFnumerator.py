#!/usr/bin/env python
# encoding: utf-8

# ================== GFFnumerator =================

# Script to re-name the IDs of the genes and features of a GFF3 file. It
# prints to the standard output, and it's designed for python3. It depends on
# the library gffutils.

# Version 4: Changed naming scheme to match the NCBI specs better:
# - The IDs should have a "locus_tag", so the IDs look like "tag_NNNNN", where
#   NNNNN stands for numbers. The tag itself should not contain `-_*` characters 
# - Instead of "mRNA" for transcripts, the ID has have "T"
# - Added option --step to change the numbering in IDs. They now start at 0 instead of 1
# - Creates a tRNA feature for the tRNAscan annotation

# Version 3: Added argument --namestoo

# Sources
# https://pythonhosted.org/gffutils/autodocs/gffutils.create_db.html
# http://daler.github.io/gffutils/database-ids.html#merge-strategy
# https://github.com/daler/gffutils/blob/master/gffutils/test/test.py
# https://github.com/daler/gffutils/issues/61
# https://pythonhosted.org/gffutils/database-schema.html
# ==================================================
# Sandra Lorena Ament Velasquez
# Johannesson Lab, Evolutionary Biology Center, Uppsala University, Sweden
# 2019/03/26
# +++++++++++++++++++++++++++++++++++++++++++++++++

# ------------------------------------------------------
import sys # To exit the script 
import os # For the input name
import argparse  # For the fancy options
# import time
import datetime
import gffutils
import gffutils.inspect as inspect
import re
# ------------------------------------------------------

version = 4.04
versiondisplay = "{0:.2f}".format(version)

# ============================
# Make a nice menu for the user
# ============================
parser = argparse.ArgumentParser(description="* Script to re-name the IDs of the genes and features of a GFF3 file *")  # Create the object using class argparse

# Add options
parser.add_argument('GFF', help="Sorted GFF3 file")
parser.add_argument('--sample', '-s', help="String representing sample that gets appended into the gene IDs (or the locus_tag of NCBI). Default: 'FUN'", default='FUN')
parser.add_argument('--namestoo', '-n', help="Change the names of the genes too, not only the IDs", default=False, action='store_true')
parser.add_argument("--printdb", "-p", help="Print the database into a file instead of saving it in memory", default=False, action='store_true')
parser.add_argument('--version', '-v', action='version', version='%(prog)s ' + versiondisplay)

# Decoration
parser.add_argument('--step', '-t', help="Interval between gene IDs. NCBI suggests a step of 10 in case future annotation updates find new genes in between existing ones. Default: 1", default='1', type=int)


try:
	# ArgumentParser parses arguments through the parse_args() method You can
	# parse the command line by passing a sequence of argument strings to
	# parse_args(). By default, the arguments are taken from sys.argv[1:]
	args = parser.parse_args()
	GFFopen = open(args.GFF, 'r')
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

if args.printdb:
	input_base = os.path.splitext(args.GFF)[0] # Taking out the prefix of the file
	input_name = os.path.basename(input_base) # Remove the path
	dbfnchoice = input_name + '.db'
else:
	dbfnchoice = ':memory:'

# http://daler.github.io/gffutils/database-ids.html
id_spec={"gene": ["ID", "Name"], 
	"mRNA": ["ID", "transcript_id"], 
	"repeat": ["repeat", "simple_repeat", "similarity"]
	} 


db = gffutils.create_db(data = args.GFF, 
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
nogenes = len(genes)

## Get iterator of all repeats
repeats = [r for r in db.features_of_type("repeat")]
norepeats = len(repeats)

# Get ids of repeats in database
allIDs = [feature.id for feature in db.all_features() if feature.featuretype in ["gene", "repeat"]]

## Make new IDs for the genes and the repeats
newIDs = [args.sample + "_" + "{0:07d}".format(n) for n in range(args.step, (nogenes + norepeats)*args.step + 1, args.step)]

# ---------------------------
### Change IDs
# ---------------------------

def getnewID(id):
	indexfeature = allIDs.index(id)
	newid = newIDs[indexfeature]
	return(newid)

def gen():
	## Get the original headers
	headcount = 1
	for line in GFFopen:
		if (headcount == 1) and ("##gff-version" not in line): sys.stdout.write("##gff-version 3\n")
		if line.startswith( '#' ):
			sys.stdout.write(line)
		headcount += 1
	GFFopen.close()

	# Add a line to mark the file with this script
	now = datetime.datetime.now()
	newhead = '# Original file ' + os.path.basename(args.GFF) + ' modified with GFFnummerator.py v. ' + str(versiondisplay) + ' on ' + str(now) +  '\n'
	sys.stdout.write(newhead)

	# Genes and repeats
	focalfeatures = [feature for feature in db.all_features() if feature.featuretype in ["gene", "repeat"]]

	# Actual GFF
	for thing in focalfeatures:
		# The annotation of trnascan is in the name
		trnastatus = False
		if 'trnascan' in thing['ID'][0]: # Produced by tRNAscan, so a specific format is expected	
			trnascanpattern = re.compile('trnascan-[\w\.]*-noncoding-([\w]*)_([\w]{3})-[\w\d\.-]*')
			matchy = trnascanpattern.search(thing['ID'][0])
			if matchy:
				trnastatus = True
				thing.attributes['Note'] = f"tRNA-{matchy.group(1)}_anti-codon_{matchy.group(2)}"

		newidthing = getnewID(thing.id)
		thing['ID'] = newidthing

		if args.namestoo: thing['Name'] = newidthing

		print() # I like to have space to see the individual items
		print(thing)

		if thing.featuretype == "gene":
			# The tRNAscan annotation has only the gene feature, but not children
			if trnastatus:
				newchild = thing
				newchild.source = "tRNAscan-SE"
				newchild.featuretype = "tRNA"
				newchild['ID'] = newidthing + '-tRNA'
				newchild['Parent'] = newidthing
				
				if matchy.group(2) == "NNN":
					newchild.attributes['product'] = f"tRNA-Xxx" # https://www.ncbi.nlm.nih.gov/genbank/eukaryotic_genome_submission_annotation/
					newchild.attributes['Note'] = f"probable {matchy.group(1)} tRNA"
				else:
					newchild.attributes['product'] = f"tRNA-{matchy.group(1)}"
					newchild.attributes['Note'] = f"probable {matchy.group(1)} tRNA: anti-codon {matchy.group(2)}"
				print(newchild)

			# Normal coding genes
			mRNA = 1
			for child in list(db.children(thing)):	
				if child.featuretype == "mRNA":
					child['Parent'] = newidthing

					# make a new ID for the mRNA
					newidrna = newidthing + "-T" + str(mRNA)
					mRNA += 1
					child['ID'] = newidrna
					if args.namestoo: child['Name'] = newidrna
		
					print(child)

					# Update all the features that depend on the mRNA, assuming each feature has a SINGLE parent
					typeids = {'gene':1, 'mRNA':1, 'exon':1, 'CDS':1, 'five_prime_UTR':1, 'three_prime_UTR':1, 'intron':1, 'repeat':1, 'sequence_feature': 1 }
					for grandchild in list(db.children(thing, level = 2)):
						grandchild['Parent'] = newidrna # TODO Add something here, a conditional, in case there are multiple parents

						typefea = grandchild.featuretype # What type are we dealing with?
						if typefea == 'five_prime_UTR': typename = '5pUTR'
						elif typefea == 'three_prime_UTR': typename = '3pUTR'
						else: typename = typefea

						newidchild = newidrna + '-' + typename + "{0:01d}".format(typeids[typefea])
						grandchild['ID'] = newidchild
						# if args.namestoo: grandchild['Name'] = newidchild # They didn't have a name to beggining with in MAKER

						typeids[typefea] += 1 # Increase the count of the corresponding type for this gene
						print(grandchild)
gen()


# ---------------------------------
# Learning how to use gffutils 
# ---------------------------------
## What features does it have?
# print(db_results)
# print(list(db_results['featuretype'].keys()))
# allfeats = list(db.all_features()) # Iterator of all features
# allIDsallfeats = [f.id for f in db.all_features()]
# print(allIDsallfeats)


## Get number of elements of a feature type
# nogenes = db.count_features_of_type("gene")

## Loop
# for gene in db.features_of_type('gene', order_by='start'):
# 	print(gene['ID'])

## Print the whole thing
# for gene in db.features_of_type('gene'):
# 	print(gene)
# 	for child in list(db.children(gene)):
# 		print(child)
# --------------------------------------------

# def newline(cols, newattributes):
# 	newcols = cols 
# 	newcols[8] = newattributes
# 	newline = '\t'.join(newcols) + '\n' # Stitch it together as a line
# 	return(newline)

# ---------------------------------
# Name of output
# ---------------------------------
# if not args.outputname:
# 	input_base = os.path.splitext(args.GFF)[0] # Taking out the prefix of the file
# 	input_name = os.path.basename(input_base) # Remove the path
# 	database_filename = input_name + '.db'
# 	# newgff = open(input_name + '_newID.gff', "w")
# else:
# 	database_filename = args.outputname + '.db'
# 	# newgff =  open(args.outputname, "w")
# dir(gene) --> 'astuple', 'attributes', 'bin', 'calc_bin', 'chrom', 'dialect', 'end', 'extra', 'featuretype', 'file_order', 'frame', 'id', 'keep_order', 'score', 'seqid', 'sequence', 'sort_attribute_values', 'source', 'start', 'stop', 'strand'
# ---------------------------------

# for line in GFFopen:
# 	if '##gff-version' in line:
# 		sys.stdout.write(line)

# 		# Add a line to mark the file with this script
# 		now = datetime.datetime.now()
# 		newhead = '# IDs of ' + os.path.basename(args.GFF) + ' renamed with GFFnummerator.py v. ' + str(versiondisplay) + ' on ' + str(now) +  '\n'
# 		sys.stdout.write(newhead)

# 	elif '#' in line: # Print headers as it is
# 		sys.stdout.write(line)

# 	elif line not in ['\n', '\r\n']: # Ignore empty lines
# 		cols = line.rstrip("\n").split("\t")
# 		print(cols)


