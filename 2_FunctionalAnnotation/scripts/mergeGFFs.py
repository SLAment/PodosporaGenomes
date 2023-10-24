#!/usr/bin/env python
# encoding: utf-8

# ================== mergeGFFs =================
# Append missing features in one gff3 file from another
# ==================================================
# Sandra Lorena Ament Velasquez
# 2023/05/17
# +++++++++++++++++++++++++++++++++++++++++++++++++

# ------------------------------------------------------
import os # For the input name
import sys # For reading the input
import argparse # For the fancy options
import gffutils
import re # Regular expressions
# ------------------------------------------------------
version = 1.1
versiondisplay = "{0:.2f}".format(version)

# ============================
# Check input file
# ============================
# Make a nice menu for the user
parser = argparse.ArgumentParser(description="* Append missing features in one gff3 file from another. Output is unsorted! *", epilog="BLAST must be locally installed.") # Create the object using class argparse

# Add options
parser.add_argument('firstgff', help="Base GFF3 file")
parser.add_argument('secondgff', help="GFF3 file with additional features missing in the first one")
parser.add_argument('--missing', '-m', help="Name of a file to write the missing items alone")
parser.add_argument('--version', "-v", action='version', version='%(prog)s ' + versiondisplay)

try:
	# ArgumentParser parses arguments through the parse_args() method You can
	# parse the command line by passing a sequence of argument strings to
	# parse_args(). By default, the arguments are taken from sys.argv[1:]
	args = parser.parse_args()
except IOError as msg:  # Check that the file exists
    parser.error(str(msg)) 
    parser.print_help()
# ============================

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

first_db = gffutils.create_db(data = args.firstgff, 
	keep_order = True,
	dbfn = dbfnchoice,
	# force = True, # force=True overwrite any existing databases.
	id_spec = id_spec, 
	merge_strategy = "create_unique") # Add an underscore an integer at the end for each consecutive occurrence of the same ID 
	# verbose = True,) 

second_db = gffutils.create_db(data = args.secondgff, 
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

# Get ids of repeats in database
allIDs_1 = [feature.id for feature in first_db.all_features() if feature.featuretype in ["gene", "repeat"]]
allIDs_2 = [feature.id for feature in second_db.all_features() if feature.featuretype in ["gene", "repeat"]]

missingthings = [thing for thing in allIDs_2 if thing not in allIDs_1]

# Print first file, which hopefully has a header
for line in open(args.firstgff, 'r'):
	sys.stdout.write(line)

# Write an additional file for the missing items alone
if args.missing != None: missingout = open(args.missing, 'w')

# Print the missing items
for thing in missingthings:
	print(second_db[thing])
	if args.missing != None: missingout.write(f"{second_db[thing]}\n")

	for child in second_db.children(thing, level = 1):
		print(child)
		if args.missing != None: missingout.write(f"{child}\n")
	for grandchild in second_db.children(thing, level = 2):
		print(grandchild)
		if args.missing != None: missingout.write(f"{grandchild}\n")

