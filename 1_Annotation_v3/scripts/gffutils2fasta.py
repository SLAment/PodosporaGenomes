#!/usr/bin/env python
# encoding: utf-8

# import time
# start_time = time.time()

# ================== gffutils2fasta =================
# Script to extract fasta subsequences out of an input fasta file using a
# corresponding GFF3 file.

# SOURCES
# https://pythonhosted.org/gffutils/autodocs/gffutils.create_db.html
# http://daler.github.io/gffutils/database-ids.html#merge-strategy
# https://github.com/daler/gffutils/blob/master/gffutils/test/test.py
# https://github.com/daler/gffutils/issues/61
# https://pythonhosted.org/gffutils/database-schema.html
# https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md

# CAVEATS:
# - The script assumes that the gene is the higher level in the attributes (i.e. polycistronic transcripts are not compatible)
# - The script CAN'T deal with genes with multiple isoforms (mRNAs). If forced, it might fuse exons were it shouldn't.
# - It can't deal with circular features


## VERSION notes
# version 1.7 - Made it more tolerable of CDS without IDs and added an option to extract pseudogenes
# version 1.5 - "from Bio.Alphabet import generic_dna" was removed from BioPython >1.76. See https://biopython.org/wiki/Alphabet
# ==================================================
# Sandra Lorena Ament Velasquez
# Johannesson Lab, Evolutionary Biology Center, Uppsala University, Sweden
# 2019/03/27
# +++++++++++++++++++++++++++++++++++++++++++++++++

# ------------------------------------------------------
import argparse  # For the fancy options
from Bio import SeqIO
from Bio.Seq import Seq
import sys  # To exit the script
import os  # For the input name
import re
import gffutils
# ------------------------------------------------------

version = 1.70
versiondisplay = "{0:.2f}".format(version)
supportedtypes = ["gene", "CDS", "cds", "exon", "noutrs", "similarity", "expressed_sequence_match", "repeat", "pseudogene"] # Unlike the CDS, Exons may contain the UTRs; noutrs is from start to stop codon without introns in nuleotides

# ============================
# Make a nice menu for the user
# ============================
parser = argparse.ArgumentParser(description="* Extract GFF3 features to fasta file *", epilog="The name in the first column of the GFF3 has to be the same as the sequence in the fasta file.")  # Create the object using class argparse

# Add options
parser.add_argument('fastafile', help="Fasta file with the genomic secuences")
parser.add_argument('GFF', help="GFF3 file")
parser.add_argument("--type", "-t", help="The feature you want to extract. Options are: "+str(supportedtypes)+", Default: gene. Notice that 'exon' would include the UTRs if present (eg. MAKER output). The option 'noutrs' includes all the body of the gene, except UTRs.", default='gene', choices = supportedtypes)
parser.add_argument('--proteinon', '-p', help="Return protein sequences (only useful for CDS)", default=False, action='store_true')
parser.add_argument('--join', '-j', help="Join together the features that belong to a single gene (only useful for CDS and exon)", default=False, action='store_true')
parser.add_argument('--specificgene', '-g', help="Extract only this specific gene or list of genes separated by commas and no spaces (give ID or Name of gene)") 
parser.add_argument("--extrabp", "-e", help="Extra base pairs on the side (only works with -t gene, default 0)", type=int, default=0)

parser.add_argument("--onlyids", "-n", help="Only keep the name of the gene ID in the output", default=False, action='store_true')
parser.add_argument("--onlynames", "-N", help="Only keep the name of the gene in the output", default=False, action='store_true')

# # parser.add_argument("--appendtoname", "-a", help="Append string to the sequence names", default='') # Maybe one day
# # parser.add_argument("--othertype", "-t", help="The feature you want to extract, with your own name", default=NULL) # One day

# # extras
parser.add_argument('--output', '-o', help="Name of output file (default is to append 'ed' to the input name)", default=None)
parser.add_argument('--version', '-v', action='version', version='%(prog)s ' + versiondisplay)
parser.add_argument('--verbose', '-b', help="Print a little bit of extra information", default=False, action='store_true')

try:
	# ArgumentParser parses arguments through the parse_args() method You can
	# parse the command line by passing a sequence of argument strings to
	# parse_args(). By default, the arguments are taken from sys.argv[1:]
	args = parser.parse_args()
	GFFopen = open(args.GFF, 'r')
	fastaopen = open(args.fastafile, 'r')
except IOError as msg:  # Check that the file exists
	parser.error(str(msg)) 
	parser.print_help()

# ------------------------------------------------------

# I often write it in small case, so to catch that
if args.type == "cds": args.type = "CDS" # it's necessary in capitals for later

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
	"repeat": ["ID", "Name"], 
	"similarity": ["Target"], 
	"pseudogene": ["ID", "Name"], 
	"pseudogenic_transcript": ["ID", "Name"], 
	"expressed_sequence_match": ["ID", "Name"]} 

db = gffutils.create_db(data = args.GFF, 
	keep_order = True,
	dbfn = dbfnchoice,
	# force = True, # force=True overwrite any existing databases.
	id_spec = id_spec, 
	verbose = True,
	merge_strategy = "create_unique") # Add an underscore an integer at the end for each consecutive occurrence of the same ID 
print() # if verbose = True

# t1 = time.time()
from gffutils import inspect; db_results = inspect.inspect(db) # Report
# print("\n\nIt took {0:.1f}s to create database".format(t1 - t0))
# ---------------------------------


# ---------------------------
### Read the fasta file
# ---------------------------
# This stores in memory
records_dict = SeqIO.to_dict(SeqIO.parse(fastaopen, "fasta"))

# ---------------------------
# Names for input and output
# ---------------------------
if args.output: # User defined
	output_handle = open(args.output, "w")
else:
	input_name = os.path.splitext(args.GFF)[0]
	input_name = os.path.basename(input_name) # Put the output in the working directory, not in the input dir
	if args.extrabp > 0:
		output_handle = open(input_name + '_' + str(args.type) + '_p' + str(args.extrabp) +'.fas', "w")
	else:
		output_handle = open(input_name + '_' + str(args.type) + '.fas', "w")

# ---------------------------
### Define some useful functions
# ---------------------------

def getseqbasic(dbobject, seq_record, extrabp = args.extrabp):
	""" Get the sequence of a simple feature """
	start = dbobject.start - 1 - extrabp # The minus one is because in the GFF3 the count starts at 1, but in python it starts at 0s
	stop = dbobject.end + extrabp # Because I use string slice, the last element won't be included
	dbobjectseq = seq_record[start:stop] # Precisely because GFF3 is based 1, so no - 1 is needed
	return(dbobjectseq)

def seqnamer(geneID, genename, geneseq, typeseq = args.type): #, onlyids = args.onlyids, onlynames = args.onlynames):
	""" Rename sequence so it's not the chromosome name """
	if args.onlyids: # Name of the output
		geneseq.id = geneID
		geneseq.description = ''
	elif args.onlynames: # Name of the output
		geneseq.id = genename
		geneseq.description = ''
	else:
		geneseq.id = geneID + '|' + genename + '|' + typeseq + '|'	
	return(geneseq)

# ---------------------------
## Only look for some specific gene
# ---------------------------
if args.specificgene: # Only one specific gene, or string of genes, is required
	focalgenes = [gene for gene in args.specificgene.split(",")]

# ---------------------------
## Retrieve the desire features
# ---------------------------
# For pseudogenes
if args.type == "pseudogene": 
	for gene in db.features_of_type('pseudogene'):
		geneID = gene['ID'][0]
		seq_record = records_dict[gene.chrom] # The chromosome sequence
		# Get DNA sequence
		geneseq = getseqbasic(gene, seq_record)
		# Update the naming of the sequence
		geneseq = seqnamer(geneID, geneID, geneseq, typeseq = args.type)

		# Print the sequence 
		SeqIO.write(geneseq, output_handle, "fasta")


### Typical format from a RepeatMasker output
if args.type == "similarity": 
	for gene in db.features_of_type('similarity'):
		geneID = gene['Target'][0]
		seq_record = records_dict[gene.chrom] # The chromosome sequence
		# Get DNA sequence
		geneseq = getseqbasic(gene, seq_record)
		# Update the naming of the sequence
		geneseq = seqnamer(geneID, geneID, geneseq, typeseq = 'repeat')

		# Print the sequence 
		SeqIO.write(geneseq, output_handle, "fasta")

### blastn alignments or the output of `gtfRM2gff.py`
elif (args.type == "expressed_sequence_match") or (args.type == "repeat"):
	for gene in db.features_of_type(args.type):
		geneID = gene['ID'][0]
		genename = gene['Name'][0] # same as gene.attributes['Name'][0]
		seq_record = records_dict[gene.chrom] # The chromosome sequence
		# Get DNA sequence
		geneseq = getseqbasic(gene, seq_record)
		# Update the naming of the sequence
		geneseq = seqnamer(geneID, genename, geneseq, typeseq = args.type)

		# Print the sequence 
		SeqIO.write(geneseq, output_handle, "fasta")

### More canonical genes
else:
	childwarning = False
	for gene in db.features_of_type('gene'):
		# dir(gene) --> 'astuple', 'attributes', 'bin', 'calc_bin', 'chrom', 'dialect', 'end', 'extra', 'featuretype', 'file_order', 'frame', 'id', 'keep_order', 'score', 'seqid', 'sequence', 'sort_attribute_values', 'source', 'start', 'stop', 'strand'
		geneID = gene['ID'][0]

		## -- Some gffs don't have a name feature (eg. funannotate output)
		try:
			genename = gene['Name'][0] # same as gene.attributes['Name'][0]
		except:
			genename = "" # An empty string
		# --
		
		seq_record = records_dict[gene.chrom] # The chromosome sequence

		if args.specificgene and ((geneID not in focalgenes) and (genename not in focalgenes)):
			continue # If the gene is not in the user's list, then skip the rest of the code and go to the next gene

		if args.type == 'gene':
			geneseq = getseqbasic(gene, seq_record) # Get the sequence for this gene
			geneseq = seqnamer(geneID, genename, geneseq, typeseq = args.type) # Update the naming of the sequence

			# Print the sequence 
			SeqIO.write(geneseq, output_handle, "fasta")

		elif args.type == 'noutrs':
			# We want from the start to the stop codon, including the introns, but excluding the UTRs
			allchildren = [child for child in db.children(gene, featuretype='CDS', order_by='start')]
			
			if len(allchildren) >= 1: # Some genes might not have CDS
				start = allchildren[0].start - 1 - args.extrabp # Start of gene excluding UTRs
				stop = allchildren[len(allchildren) - 1].end + args.extrabp # End of gene excluding UTRs
				geneseq = seq_record[start:stop] # Precisely because GFF3 is based 1, so no - 1 is needed
				geneseq = seqnamer(geneID, genename, geneseq, typeseq = 'gene_noutrs') # Update the naming of the sequence

				# Print the sequence 
				SeqIO.write(geneseq, output_handle, "fasta")

		# Getting exons	
		elif args.type == 'exon':
			child_counter = 1
			child_concat = Seq('')

			for child in db.children(gene, featuretype=args.type, order_by='start'): # if the exon has children features like CDS and UTRs.
				parents = db.parents(child, featuretype='mRNA') # Assuming there are mRNA features in the gtf/gff file
				childID = child['ID'][0]
				exonparent = list(parents)[0]['ID'][0]
					
				start = child.start - 1
				stop = child.end

				if args.join:
					# Get DNA seq
					child_concat += seq_record[start:stop] # Last letter will be excluded
				else:
					# Get DNA seq
					exonseq = seq_record[start:stop] # Last letter will be excluded

					# The naming is a bit different than the normal, so I cannot use the function seqnamer()
					if args.onlyids: # Name of the output
						exonseq.id = childID
						exonseq.description = ''
					elif args.onlynames: # Name of the output
						exonseq.id = genename + '_exon.' + str(child_counter)
						exonseq.description = ''
					else:
						# Print them into the file (and append a number for each CDS)
						exonseq.id = childID + '|' + exonparent + '|exon.' + str(child_counter) + '|'
					
					SeqIO.write(exonseq, output_handle, "fasta")
				

				child_counter += 1

			if args.join:
				child_concat = seqnamer(geneID, genename, child_concat, typeseq = 'concat_exons') # Update the naming of the sequence
				SeqIO.write(child_concat, output_handle, "fasta")

		# Getting CDS	
		elif args.type == 'CDS':
			child_counter = 1
			child_concat = Seq('')

			for child in db.children(gene, featuretype=args.type, order_by='start'):

				## --- Get the cdsparent name from the mRNA ---
				# parents = db.parents(child, featuretype='mRNA') # Might be useful for different transcripts, but I don't encounter that problem
				# cdsparent = list(parents)[0]['ID'][0]
				## --------------------------------------------
				
				## -- Some gffs don't have an ID for their CDS
				try:
					childID = child['ID'][0]
				except:
					childID = ""
					childwarning = True # To print later

				# --

				cdsparent = geneID
				strand = child.strand
				phase = int(child.frame)

				if (strand != '+') and (strand != '-'):
					print("There is no information of strand direction!")
					sys.exit(1)

				# Produce a protein sequence for the CDS in each exon (taking care of the phases)
				if not args.join:
					if strand == '+':
						start = child.start - 1 + phase
						raw_stop = child.end

						# Trim the extra bases at the end that are lost because the CDS are not joint
						howmanycodons = (raw_stop - start) // 3
						stop = (howmanycodons * 3) + start 

						# Get DNA seq and translate
						cdsseq = seq_record[start:stop] # Last letter will be excluded

						# If protein sequences are required
						if args.proteinon:
							cdsseq.seq = cdsseq.seq.translate()

					elif strand == '-': 
						raw_start = child.start - 1
						stop = child.end - phase 

						# Trim the extra bases at the end that are lost because the CDS are not joint
						howmanycodons = (stop - raw_start) // 3
						start = stop - (howmanycodons * 3)

						if args.proteinon: # Get DNA seq, reverse complement, and translate
							cdsseq = seq_record[start:stop].reverse_complement()
							cdsseq.id = seq_record.id
							cdsseq.description = seq_record.description
							cdsseq.seq = cdsseq.seq.translate()
						else:
							cdsseq = seq_record[start:stop]

					if args.onlyids: # Name of the output
						cdsseq.id = childID
						cdsseq.description = ''
					elif args.onlynames: # Name of the output
						cdsseq.id = genename + '_CDS.' + str(child_counter)
						cdsseq.description = ''
					else:
						# Print them into the file (and append a number for each CDS)
						cdsseq.id = childID + '|' + cdsparent + '|CDS.' + str(child_counter) + '|'
				
					SeqIO.write(cdsseq, output_handle, "fasta")
					
					child_counter += 1
				# Produce a single CDS for the entire gene, linked to **hey**
				else:
					# Get phase for the very first exon
					# phase = int(feature[7])

					if strand == '+':
						if child_counter == 1: # Consider it only on the beginning of the protein
							start = child.start - 1 + phase
							child_counter += 1
						else: # All the other cds should be in frame
							start = child.start - 1

					elif strand == '-':
						start = child.start - 1

					stop = child.end
					# Get DNA seq
					child_concat += seq_record[start:stop] # Last letter will be excluded

					# SeqIO.write(child_concat, output_handle, "fasta")


			if args.join and args.type == 'CDS': # **hey** Print into the file if the output is to be joined
				try: # Sometimes genes won't have a CDS feature, in which case we better move on...
					if args.proteinon:

						if strand == '+':	
							cdsseq_concat = child_concat.seq.translate()
							child_concat.seq = cdsseq_concat
						elif strand == '-':
							# Use last phase to correct the reading frame in the minus strand for this gene (and reverse complement it)
							child_concat = child_concat[:len(child_concat) - phase].reverse_complement()
							# Translate to protein
							cdsseq_concat = child_concat.seq.translate()
							child_concat.seq = cdsseq_concat
							# Edit Seq object to have the same names again
							child_concat.id = seq_record.id
							child_concat.description = seq_record.description

					child_concat = seqnamer(geneID, genename, child_concat, typeseq = args.type) # Update the naming of the sequence

				except:
					if args.verbose: print("The gene " + gene['ID'][0] + " " + gene['Name'][0] + " is weird. Maybe it doesn't have a CDS? Skipped.")

				# Finally write the output sequence
				SeqIO.write(child_concat, output_handle, "fasta")
	if childwarning and args.onlyids and not args.join: print("WARNING: The CDS have no IDs so the sequences will have no names!")


# Close all opened files
output_handle.close()
fastaopen.close()
GFFopen.close()

"""
gene --> A region (or regions) that includes all of the sequence elements necessary to encode a functional transcript.
exons --> A region of the transcript sequence within a gene which is not removed from the primary RNA transcript by RNA splicing.
CDS --> A contiguous sequence which begins with, and includes, a start codon and ends with, and includes, a stop codon.

A GFF3 file has the following fields

##gff-version 3
0	ref_seq	Thamnoliak127a3n4s400-scf_the-scaffold
1	source	SNAP
2	type	gene	( This is constrained to be either: (a) a term from the "lite" sequence ontology, SOFA; or (b) a SOFA accession number. The latter alternative is distinguished using the syntax SO:000000.)	
3	start_position	835
4	stop_position	1002
5	score	16.007
6	strand	+
7	phase_codon	.
8	attributes	ID=snap.gene.000001; Name=Thamnoliak127a3n4s400-scf_the-scaffold-snap.1
"""

# ---------------------------------
# Learning how to use gffutils 
# ---------------------------------
# What features does it have?
# print(db_results)
# print(list(db_results['featuretype'].keys()))
# allfeats = list(db.all_features()) # Iterator of all features
# allIDsallfeats = [f.id for f in db.all_features()]
# print(allIDsallfeats)

# # Print the whole thing
# for gene in db.features_of_type('gene'):
# 	print(gene)
# 	for child in list(db.children(gene)):
# 		print(child)
# --------------------------------------------

# print("--- %s seconds ---" % (time.time() - start_time))


