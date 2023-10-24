#!/usr/bin/env python
# encoding: utf-8

# ================== tbl4ncbi =================

# NCBI is very particular about the submitted tbl files, so try to get closer
# to what they want. The original tbl file was produced by `funannotate util
# gff2tbl` from a gff3 file. This is part of the `gff2function.smk` pipeline.

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
parser.add_argument('--genelist', '-g', help="List of GeneIDs and their desired qualifiers with description (separated by a tab)")

# extras
parser.add_argument('--version', '-v', action='version', version='%(prog)s ' + versiondisplay)

try:
	# ArgumentParser parses arguments through the parse_args() method You can
	# parse the command line by passing a sequence of argument strings to
	# parse_args(). By default, the arguments are taken from sys.argv[1:]
	args = parser.parse_args()
	tblopen = open(args.tbl, 'r')
except IOError as msg:  # Check that the file exists
    parser.error(str(msg)) 
    parser.print_help()
# ============================

qualifiers =  {}
if args.genelist != None:
	listopen = open(args.genelist, 'r')
	for line in listopen:
		cols = line.rstrip("\n").split("\t")
		qualifiers[cols[0]] = (cols[1], cols[2])


for line in tblopen:
	cols = line.rstrip("\n").split("\t")		# break the line into columns defined by the tab
	if 'locus_tag' in line: # this indicates a new gene
		unprocessed = False
		pseudo = False
		locus_tag = cols[4]
		sys.stdout.write(line)

		# in the funannotate output the locus_tag is always the last thing in the gene part, so use that to print the pseudos qualifiers if required
		if locus_tag in qualifiers.keys():
			# https://www.ncbi.nlm.nih.gov/genbank/eukaryotic_genome_submission_annotation/
			if qualifiers[locus_tag][0] == "pseudogene=unprocessed":
				newline = f"\t\t\tpseudogene\tunprocessed"
				unprocessed = True
			elif qualifiers[locus_tag][0] == "pseudo=true":
				newline = f"\t\t\tpseudo"
				pseudo = True
			print(newline)
			print(f"\t\t\tnote\t{qualifiers[locus_tag][1]}")
				
	elif 'protein_id' in line: # funannotate util gff2tbl used the name of the transcript for the protein_id, but I don't like this
		# Here we are asumming a structure like \t\t\tprotein_id\tgnl|ncbi|QC762_102900-T1
		protein_id = cols[4].rstrip('-T1') # Remove the T1 in the end that comes from the pipeline's mRNA ID
		newline = f"\t\t\tprotein_id\t{protein_id}"
		print(newline)
		if pseudo: print(f"\t\t\tpseudo")
		if unprocessed: print(f"\t\t\tpseudogene\tunprocessed")
	# Fixed now with `fixProductGff.py`
	# elif 'product\t' in line: # Here is where we try to fix all the typos
	# 	# Full thing replacement
	# 	if '\t\t\tproduct\tconserved fungal protein\n' == line: line = "\t\t\tproduct\thypothetical protein\n"
	# 	elif '\t\t\tproduct\tBelongs to the peptidase S8\n' == line: line = "\t\t\tproduct\tMycosin-4\n" # https://www.uniprot.org/uniprotkb/I6YC58/entry
	# 	elif '\t\t\tproduct\tproteinral transcription repressor\n' == line: line = "\t\t\tproduct\tGeneral transcriptional corepressor TUP1\n" # https://www.uniprot.org/uniprotkb/P16649/entry
	# 	elif '\t\t\tproduct\tproteinral amino acid permease agp2\n' == line: line = "\t\t\tproduct\tGeneral amino acid permease AGP2\n" # https://www.uniprot.org/uniprotkb/P38090/entry
	# 	elif '\t\t\tproduct\tproteinral negative regulator of transcription subunit 5\n' == line: line = "\t\t\tproduct\tGeneral negative regulator of transcription subunit 5\n" # https://www.uniprot.org/uniprotkb/Q12514/entry
	# 	elif '\t\t\tproduct\ttranscriptional repressor proteinral negative regulator of transcription subunit 4\n' == line: line = "\t\t\tproduct\tGeneral negative regulator of transcription subunit 4\n" # https://www.uniprot.org/uniprotkb/Q12514/entry
	# 	else:
	# 		if 'bioproteinsis' in line: line = line.replace("bioproteinsis", "biogenesis")
	# 		if 'diacyglycerol' in line: line = line.replace("diacyglycerol", "diacylglycerol")
	# 		if 'morphoproteinsis' in line: line = line.replace("morphoproteinsis", "morphogenesis")
	# 		if 'phopho' in line: line = line.replace("phopho", "phospho")
	# 		if 'redicted' in line: line = line.replace("redicted", "utative")
	# 		if 'robable' in line: line = line.replace("robable", "utative")
	# 		if 'human' in line: line = line.replace("human ", "")
			
	# 	sys.stdout.write(line)
	elif 'transcript_id' in line:
		transcript_id = cols[4].rstrip('-T1_mrna') # Remove the T1 in the end that comes from the pipeline's mRNA ID
		newline = f"\t\t\ttranscript_id\t{transcript_id}_rna"
		print(newline)
	elif 'CFMR\t12345' in line: # I think this is a placeholder? (CFMR: The Reference Culture Collection at the Center for Forest Mycology Research)
		pass #remove it
	else:
		sys.stdout.write(line)


