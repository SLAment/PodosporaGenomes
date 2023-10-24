#!/usr/bin/env python
# encoding: utf-8
from __future__ import print_function

# ================== GffgenesIDFix =================
# Script to append the names of a reference database of genes to a new annotation in another genome.
# BLAST must be locally installed.

# Version 2.5: 
# Updated to remove the use of Bio.Alphabet
# Made for python3, but adapted back to python2. From version 2.3 I went back
# to python3 however.

# Whenever there is a reference genome split into several new target models,
# the script attaches a subindex in the form originalgene.1, originalgene.2,
# etc.

# There is a stochastic element related to the order in which the dictionaries
# keys are read, I suspect.

# Unfortunately, for hard cases sometimes names might be overlapping (same
# name for more than one gene). Example: maker-chromosome_1-augustus-gene-13.5

# ==================================================
# Sandra Lorena Ament Velasquez
# 2017/06/16 - 2023/05/05
# +++++++++++++++++++++++++++++++++++++++++++++++++

# ------------------------------------------------------
import os # For the input name
import sys # For reading the input
# import re
import argparse # For the fancy options
from Bio import SeqIO # biopython
from Bio.Blast.Applications import * # NCBI BLASTX wrapper from the Bio.Blast.Applications module
# from Bio.Alphabet import generic_dna
import subprocess # For the database
from collections import defaultdict #This tells python that the dictionary contains a list so you can freely append things to the value of each key
from shutil import rmtree
import re # For the sorting function
# ------------------------------------------------------
version = 2.53
versiondisplay = "{0:.2f}".format(version)

# ============================
# Check input file
# ============================
# Make a nice menu for the user
parser = argparse.ArgumentParser(description="* Append the names of a reference database of genes to a new annotation in another genome. *", epilog="BLAST must be locally installed.") # Create the object using class argparse

# Add options
parser.add_argument('refgenes', help="Multifasta file with the reference gene sequences with the desired names IN THE SAME ORDER as in the reference genome")
parser.add_argument('targetgenes', help="Multifasta file with the target gene sequences from the new annotation IN ORDER")
parser.add_argument('targetgff', help="Gff3 file (sorted) of the target genome containing the genemodels to be re-named")

parser.add_argument("--dbtype", "-d", help="Type of database for local BLAST. Default: nucl", default='nucl', choices = ["nucl", "prot"])
parser.add_argument("--identity", "-i", help="Minimum identity percentage to consider a gene homologous. Default: 95", type=float, default=95.0)
parser.add_argument("--threads", "-t", help="Number of threads for BLAST. Default: 1", default="1", type=int) # nargs='+' All, and at least one, argument
parser.add_argument('--clean', '-c', help="Erase the created directories and files when finish", default=False, action='store_true')
parser.add_argument('--noblast', '-n', help="Do not make the BLAST again (useful for testing)", default=False, action='store_true')

parser.add_argument('--output', '-o', help="Name of output file (default is to append '_ID.gff3' to the input name)", default=None)
parser.add_argument('--temp', '-u', help="Path and directory where temporary files are written in style path/to/dir/. Default: working directory", default='./')
parser.add_argument('--verbose', '-b', help="Print a little bit of extra information", default=False, action='store_true')
parser.add_argument('--superverbose', '-B', help="Print a lot more of extra information, including report files", default=False, action='store_true')
parser.add_argument('--version', "-v", action='version', version='%(prog)s ' + versiondisplay)

try:
	# ArgumentParser parses arguments through the parse_args() method You can
	# parse the command line by passing a sequence of argument strings to
	# parse_args(). By default, the arguments are taken from sys.argv[1:]
	args = parser.parse_args()
	refgenesopen = open(args.refgenes, 'r')
	targetgenesopen = open(args.targetgenes, 'r')
	targetgffopen = open(args.targetgff, 'r')
	# targetfastaopen = open(args.targetfasta, 'r')
except IOError as msg:  # Check that the file exists
    parser.error(str(msg)) 
    parser.print_help()
# ============================

if args.superverbose:
	args.verbose = True

# ------------------------------------------------------
# Define functions
# ------------------------------------------------------
# Get a clean name for a file
def namebase(file):
	input_name = os.path.splitext(file)[0] # Remove the extension
	input_name = os.path.basename(input_name) # Remove the path
	return(input_name)

# Make a local database with the given fasta file
def makeBLASTdb(fasta, databasename, dbtype):
	# databasename = name + '_db/' + name + '_db' # Name of the database
	createdb = "makeblastdb -in " + fasta + " -out " + databasename + " -dbtype " + dbtype + " -parse_seqids" # BLAST command
	if args.superverbose: # Print it if necessary
		print("Blast command:", createdb)
	process = subprocess.Popen(createdb.split(), stdout=subprocess.PIPE) # pipe the command to the shell
	stdout, stderr = process.communicate() # rn it

def thebesthit(query, hitslist): # Assuming the best alignment is the longest and there is only one subject being hit
	alignlens = [int(hit[2]) for hit in hitslist] # make a list of the alignment lengths of all the hits
	# hitnames = [hit[0] for hit in hitslist]
	# peridents = [hit[1] for hit in hitslist]
	indexbestone = alignlens.index(max(alignlens)) # What is the index of the longest (assuming there are never two alignments of the same length)
	bestsubject = hitslist[indexbestone][0] # Get the name of the subject in that hit
	return([bestsubject, hitslist[indexbestone]]) # Return a list of the name of the subject and the entire hit

# Using the position of the genes in the 
def genNeighbours(tabcollection, key):
	keyindex = tabcollection.index(key)

	if keyindex > 0: # There is a left neighbor
		left = tabcollection[keyindex - 1]
	else:
		left = None

	if keyindex < len(tabcollection) - 1:
		right = tabcollection[keyindex + 1]
	else:
		right = None

	return([left, right])

# Got these functions from http://stackoverflow.com/questions/5967500/how-to-correctly-sort-a-string-with-a-number-inside by unutbu user
def atoi(text):
    return int(text) if text.isdigit() else text

def natural_keys(text):
    '''
    alist.sort(key=natural_keys) sorts in human order
    http://nedbatchelder.com/blog/200712/human_sorting.html
    (See Toothy's implementation in the comments)
    '''
    return [ atoi(c) for c in re.split('(\d+)', text) ]

# Print a dictionary into a table
def print_dic( dicty, output_name, header = '' ):
	ofile = open(output_name, 'w')
	if header != '':
		ofile.write((header+'\n'))

	sortedkeys = sorted(list(dicty.keys()), key = natural_keys)

	for key in sortedkeys:
		line = str(key) + '\t' + str(dicty[key]) + '\n' # Assume a single value for each key
		ofile.write(line)
	ofile.close()


# ------------------------------------------------------
# Output names
# ------------------------------------------------------

nameref = namebase(args.refgenes)
nametar = namebase(args.targetgenes)
namegff = namebase(args.targetgff)

# ------------------------------------------------------
# Deal with the input files
# ------------------------------------------------------
# Get a list of the reference gene names and order 
genes = [seq_record.id for seq_record in SeqIO.parse(refgenesopen, "fasta")]

# Get a list of the target gene names and order 
targenes = [seq_record.id for seq_record in SeqIO.parse(targetgenesopen, "fasta")]

if not args.noblast:
	# ------------------------------------------------------
	# # Create a local BLAST database
	# ------------------------------------------------------
	# Define the names of the databases
	databasename_ref = args.temp + nameref + '_db/' + nameref + '_db'
	databasename_tar = args.temp + nametar + '_db/' + nametar + '_db'

	# Make the local BLAST databases
	makeBLASTdb(args.refgenes, databasename_ref, args.dbtype)
	makeBLASTdb(args.targetgenes, databasename_tar, args.dbtype)

	# ------------------------------------------------------
	# BLAST
	# ------------------------------------------------------
	# # Define the blast command based on the database type
	# if args.dbtype == "nucl":
	# 	blast_cline = NcbiblastnCommandline(query=args.targetfasta, db=databasename, evalue=0.001, outfmt=6, out="BLAST_hits.tab", num_threads=args.threads)
	# elif args.dbtype == "prot":
	# 	blast_cline = NcbiblastxCommandline(query=args.targetfasta, db=databasename, evalue=0.001, outfmt=6, out="BLAST_hits.tab", num_threads=args.threads)


	# Target genes as query
	if args.dbtype == 'nucl':
		blast_cline = NcbiblastnCommandline(query=args.targetgenes, db=str(databasename_ref), evalue=0.001, outfmt=6, out= args.temp + "TargetVsRef_hits.tab", num_threads=args.threads)
	elif args.dbtype == 'prot':
		blast_cline = NcbiblastpCommandline(query=args.targetgenes, db=str(databasename_ref), evalue=0.001, outfmt=6, out= args.temp + "TargetVsRef_hits.tab", num_threads=args.threads)
	if args.superverbose:
		print("Blast command:", blast_cline)
	stdout, stderr = blast_cline() # Do the actual blasting

	# Reference genes as query
	if args.dbtype == 'nucl':
		blast_cline = NcbiblastnCommandline(query=args.refgenes, db=str(databasename_tar), evalue=0.001, outfmt=6, out= args.temp + "RefVsTarget_hits.tab", num_threads=args.threads)
	elif args.dbtype == 'prot':
		blast_cline = NcbiblastpCommandline(query=args.refgenes, db=str(databasename_tar), evalue=0.001, outfmt=6, out= args.temp + "RefVsTarget_hits.tab", num_threads=args.threads)
	if args.superverbose:
		print("Blast command:", blast_cline)
	stdout, stderr = blast_cline()
else:
	print("Skipping the BLASTing ...")

# ------------------------------------------------------
# Read the results back
# ------------------------------------------------------
TargetVsRef_tabs = [line.rstrip("\n").split("\t") for line in open(args.temp + "TargetVsRef_hits.tab", 'r')] 			# Read tab file into a list
RefVsTarget_tabs = [line.rstrip("\n").split("\t") for line in open(args.temp + "RefVsTarget_hits.tab", 'r')] 			# Read tab file into a list

# ------------------------------------------------------
# Output report files
# ------------------------------------------------------
if args.superverbose: fusedlocifile = open(args.temp + "fusedloci.txt", "w")

# ------------------------------------------------------
# Make a dictionary with the tab file for each query
# ------------------------------------------------------
print("Loading BLAST results ...")

# heads = "query_id	subject_id	percent_identity	alignment_length	N_mismatches	N_gaps	query_start	query_end	subject_start	subject_end	evalue	bit_score"

tarhits = defaultdict(list) #This tells python that the dictionary contains a list so you can freely append things to the value of each key
for line in TargetVsRef_tabs:
	tarhits[line[0]].append(line[1:])
	# if line[2] >= 50: tarhits[line[0]].append(line[1:])

refhits = defaultdict(list) #This tells python that the dictionary contains a list so you can freely append things to the value of each key
for line in RefVsTarget_tabs:
	refhits[line[0]].append(line[1:])
	# if line[2] >= 50: refhits[line[0]].append(line[1:])

""" First check the one-to-one's """ 
tarcount = 0 # Count the single hits
fusedloci = 0 # How many times where gene models from the reference fused into the new annotation?

onetoones = {} # Record the good one-to-one orthologs
splitgenes = {}
newnamelist = [] # To make sure there are no repeated IDs

pararef = {}
paratar = {}

print("Processing hits...")
for hit in tarhits.keys():
	# print(hit, tarhits[hit])
	# if hit == "maker-chromosome_5-snap-gene-9.147": print("hola")

	# How many different subject (reference genes) does this query (target gene) have?
	refhitbythistarget = list(set(gene[0] for gene in tarhits[hit]))

	if len(refhitbythistarget) == 1: # Only one hit for that target query
		tarcount += 1
		theonesubject = thebesthit(hit, tarhits[hit])[0] # when len(tarhits[hit]) == 1, the same as tarhits[hit][0][0]

		for ref in refhits.keys():
			# Make sure that the reference hits are not broken (there might be several hits, but actually only one subject)
			refhitnames = list(set([tarhit[0] for tarhit in refhits[ref]])) # Get list of the target genes for this given reference gene

			if (len(refhitnames) == 1) and (ref == theonesubject): # One-to-one hits
				# print("\t", refhits[ref], thebesthit(ref, refhits[ref]))

				# What is the identity percentage?
				perident = float(thebesthit(ref, refhits[ref])[1][1]) # when ref not broken, same as float(refhits[ref][0][1])

				if perident >= args.identity: # It passes, so it's consider a one-to-one ortholog
					onetoones[hit] = ref
					
				elif args.superverbose:
					print("Gene might be one-to-one to %s, but too low identity (%.2f): %s" % (ref, perident, hit))

			elif (len(refhitnames) != 1) and (ref == theonesubject): # The ref and the target hit each other, but the ref also hits other things
				# print(ref, refhits[ref])
				""" It's possible that there is only one reference gene that got split into more than one target gene """
				# --- Are the split genes continuous? (and therefore, the assignment to this reference gene is correct?)

				# Make a list of the indexes of the target genes hit by the query reference
				geneindextar = [targenes.index(this_hit) for this_hit in refhitnames]

				# Are they continuous (syntenic) genes?
				for number in range(0, len(geneindextar) - 1):
					if (geneindextar[number + 1] == geneindextar[number] + 1) or (geneindextar[number - 1] == geneindextar[number] - 1): # They are continuous 
						# # Trying to catch the issue with Pa_1_5500, where splitted genes also match overlapping genes
						# # --------
						# if hit in list(splitgenes.keys()): # This target gene is already in the collection of splitgenes
						# 	print(hit)
						# 	if ref == 'Pa_1_5500': print("Pa_1_5500")
						
						# 	# mydict = {'george':16,'amber':19}
						# 	# print(list(mydict.keys())[list(mydict.values()).index(16)])	

						# 	# # The hit is already in the dictionary, but with the same ref?

						# 	# Very difficult cases
						# 	if ref in list(pararef.keys()): # Is this reference also in the list of those with paralogs?
						# 		print( "para", hit)
						# 		print(pararef[hit])
						# 	if ref in list(onetoones.keys()): # Is this reference also in the list of those with paralogs?
						# 		print( "onetoones", hit)
						# 		print(pararef[hit])
						# # --------
						splitgenes[hit] = ref  # The reference gene was split into new genes in the targets	

					else: # It's actually some paralogs detected by the reference vs. target
						pararef[ref] = refhits[ref] # Let's deal with them later
						
	else: # The target hits many things
		# Make a list of the (unique) indexes of the reference genes hit by the query target
		geneindex = list(set([genes.index(this_hit[0]) for this_hit in tarhits[hit] if float(this_hit[1]) >= args.identity])) # Filter by identity to filter out at least obvious paralogs

		# Are they continuous (collinear) genes?
		continuous = False
		# if hit == "maker-chromosome_5-snap-gene-9.147": print("hola")

		for number in range(0, len(geneindex) - 1): # This might not work well if there are more than 3 hits, or if 1 of the 3 is not continuous
			if (geneindex[number + 1] == geneindex[number] + 1) or (geneindex[number - 1] == geneindex[number] - 1): # They are continuous 
				continuous = True
			else: # There is some inconsistency in the full list of hits, but then some real cases will be lost
				continuous = False

		if continuous:
			# --- Make a super name
			newgene = ""
			for ind in geneindex:
				newgene += genes[ind] + '-'
			newgene = newgene.rstrip("-")

			# # --- or take the first gene's name:
			# newgene = genes[min(geneindex)]

			# Does this gene name already exist?
			if newgene in newnamelist: 
				# Maybe this happened before, so get all names that match this gene
				maxnewgene = 1
				for genename in newnamelist:
					if newgene in genename: # check the whole list of genenames
						
						partsofgenename = genename.split('.') # I made the merged genes in this format "Ppc_7_3550-Ppc_5_8050-Ppc_5_8060.2"
						if len(partsofgenename) > 1: # It's not the first time we see this merged gene name
							maxnewgene = max(maxnewgene, int(partsofgenename[1]))
				
				newgene = newgene + '.' + str(maxnewgene + 1)

			newnamelist.append(newgene)

			# Append it to the one-to-one list
			onetoones[hit] = newgene
			fusedloci += 1

			# Write the gene in a report 
			if args.superverbose: fusedlocifile.write(newgene + '\n') 

			""" Notice that in this case I'm assuming synteny is enough evidence, as opposed to percentage of identity, which I'm not checking """
		else: # They are not continuous, most likely
			paratar[hit] = tarhits[hit] # Deal with them later
		


# ------------------------------------------------------
# Deal with the paralogs
# ------------------------------------------------------
print("Dealing with paralogs...")

onetooneskeys = list(onetoones.keys())
withparalogs = {} 

for ref in pararef.keys():
	# print(ref, pararef[ref])
	refleft, refright = genNeighbours(genes, ref) # Get neighbors to the sides, if possible for the current reference gene

	paralogs = list(set([r[0] for r in pararef[ref]])) # How many paralogs?

	# -- Check by synteny
	matching = False
	for par in paralogs:	
		parleft, parright = genNeighbours(targenes, par) # Neighbors for the current paralog

		if parleft in onetooneskeys: # Is already present in the one-to-one cases
			# print("left", parleft, onetoones[parleft], refleft)
			if onetoones[parleft] == refleft:
				matching = True

		if parright in onetooneskeys: # Is already present in the one-to-one cases
			# print("right", parright, onetoones[parright], refright)
			if onetoones[parright] == refright:
				matching = True

		if matching:
			# Has this gene been named already?
			if par in splitgenes.keys():
				if args.superverbose:
					print("... Re-naming annotation of %s for %s" % (par, splitgenes[par]))
				splitgenes.pop(par, None) # Remove that gene from the splitgenes dictionary

			if par not in onetoones.keys(): # Those have priority
				withparalogs[par] = ref
			
			matching = False

# ------------------------------------------------------
# Fix names of split reference genes:
# ------------------------------------------------------
splitrefgenes = set(splitgenes.values())

for locus in splitrefgenes: # For every reference gene in these collection
	subgenes = 1 
	for target in splitgenes.keys(): # Check all the target genes that hit it
		if splitgenes[target] == locus:		
			splitgenes[target] = locus + '.' + str(subgenes) # Rename it by appending a subfix
			subgenes += 1

# ------------------------------------------------------
# Report some numbers
# ------------------------------------------------------
if args.verbose:
	# How many times was a reference gene broken in new genes?
	# splitgenescount = len(set([i.split(".")[0] for i in splitgenes.values()]))

	print()
	print("Threshold for percentage of identity used:", args.identity)
	print("Number of target genes:", len(tarhits.keys()))
	print("Target genes with single hits:", tarcount)
	print("* Target genes with a one-to-one reference hit:", len(onetoones))
	print("\tNo. of reference genes fused into one target model:", fusedloci)
	print("* Target genes hitting several reference paralogs, but only one chosen:", len(withparalogs))
	print("* Target genes with split reference gene names:", len(splitgenes))
	print("\tNo. of reference genes split into several target models:", len(splitrefgenes))
	print("Total number of genes re-named:", len(onetoones) + len(withparalogs) + len(splitgenes))
	
	# Print actual gene lists
	print_dic( onetoones, args.temp + "onetoones.txt", header = 'model\tgene' )
	print_dic( splitgenes, args.temp + "splitgenes.txt", header = 'model\tgene' )
	fusedlocifile.close() # Close it since it's being written on the fly



# Merge the dictionaries (only works in python 3.5 and up) # https://stackoverflow.com/questions/38987/how-to-merge-two-python-dictionaries-in-a-single-expression
acceptednames = {**onetoones, **withparalogs, **splitgenes}

## Python 2
# acceptednames = onetoones.copy()
# acceptednames.update(withparalogs) # which returns None since it mutates acceptednames
# acceptednames.update(splitgenes)

# ------------------------------------------------------
# Read the GFF
# ------------------------------------------------------
# ---- Names for input and output ----
if args.output: # User defined
	newgff = open(args.output, "w") # Name of output
	# output_handle = open(args.output + '.fas', "w")
else:
	newgff = open(namegff + "_ID.gff3", "w")

# Read the GFF file
for line in targetgffopen:
	if ('\tgene\t' in line): # Look for the genes definition
		cols = line.rstrip("\n").split("\t")		# break the line into columns defined by the tab
		attributes = cols[8].split(";")	

		for element in attributes:
			if "ID=" in element: 
				geneID = element.split('=')[1] # Extract the current gene ID
			# elif ("Name=" in element) and (geneID in acceptednames.keys()):
			# 	# Look for the new name
			# 	newname = "Name=" + acceptednames[geneID]
			elif 'Name=' in element:
				if geneID in acceptednames.keys(): #assuming ID always comes before Name
					newname = "Name=" + acceptednames[geneID] # Look for the new name
				else:
					newname = element
				IndexName = attributes.index(element)
			elif 'ote=' in element: # Note= or note=
				pass
			else: # there is no name?			
				newname = element
			# IndexName = attributes.index(element)

		if newname in newnamelist: # This can happen in horrible situations where split genes also overlap other genes (eg. Pa_1_5500)
			newname = newname + '.2' # Assuming this doesn't happen very often ...

		attributes[IndexName] = newname # replace the old name with the new
		newattributes = ';'.join(attributes) # Re-assemble the attributes

		newcols = cols[0:8] # Constrct a new line with the first columns of the original one
		newcols.append(newattributes) # append the new attributes
		newline = '\t'.join(newcols) + '\n' # Stitch it together as a line

		# Add it to the new gff file
		newgff.write(newline)

		newnamelist.append(newname)
	else:
		newgff.write(line)

# ------------------------------------------------------
# Clean
# ------------------------------------------------------
if args.clean:
	if args.verbose:
		print("\nRemoving intermediate files ...")
	# Remove BLAST databases	
	rmtree(args.temp + nameref + '_db')
	rmtree(args.temp + nametar + '_db')
	# Remove the BLAST results
	os.remove(args.temp + "TargetVsRef_hits.tab")
	os.remove(args.temp + "RefVsTarget_hits.tab")


