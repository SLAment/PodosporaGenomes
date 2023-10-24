# -*- snakemake -*-

from ete3 import Tree
from Bio import SeqIO
from Bio.Seq import Seq
from os.path import exists
import random

### PodoPhylogeny: Making a phylogeny of the Podospora complex
#############################################################################

# A pipeline to make a phylogeny of *Podospora* complex. 

#############################################################################
# ===========================================================================
# Sandra Lorena Ament Velasquez
# 2023/08/15 - 2023/10/24
# ++++++++++++++++++++++++++++++++++++++++++++++
# Version 1

# -------------------------------------------------
# DATA from the configuration file
configfile: "config/PodoPhylogeny_config.yaml"

AllSamples = config["SampleIDs"]

# Output results from the OrthoTreesUppmax.smk pipeline
treesfile = config["treesfile"]
orthoreport = config["orthoreport"]
OrthoTreesPath = config["OrthoTreesPath"]

# Number of sample orthologs
SAMPLEsize = config["SAMPLEsize"]

# Support filter for collapsing branches
MINsupport = config["MINsupport"]

# Local installations
astral = config["astral"]
quartetscores = config["quartetscores"]
# -------------------------------------------------

# ----------
# Rules not submitted to a job
localrules: concatenation, renametrees, collapselowsupportbraches, scoring_astral, polytomy_test_astral, report_genes_astral, InternodeCertainty
# ----------

rule all:
	input:
		f"astral/astral.bb{MINsupport}.poly.tre",
		"concatenated/onetoones.treefile",
		"results/OneToOne_astral.txt",
		"results/RAxML_IC_Score_BranchLabels.onetoones.tre",

		# QuartetSupport values on the concatenated and astral phylogenies
		f"results/AstralAll.bb{MINsupport}_QuartetScores.tre",
		f"results/Concatenated1000.bb{MINsupport}_QuartetScores.tre"

# ------ Concatenated tree ------

# Make a dictionary of all Pa genes and the alignment number
orthodict = {}
tabs = [line.rstrip("\n").split("\t") for line in open(orthoreport, 'r')]

# import random
for tab in tabs:
	subtabs = tab[0].split(" ") # I accidentally used a white space instead of a tab there...
	
	orthogroup = tab[1]
	PODANref = tab[2]
	num = int(subtabs[0])
	alignfile = subtabs[1]
	
	orthodict[num] = [alignfile, orthogroup, PODANref]

# Select random numbers
if not exists("concatenated/genes_selected.txt"):
	randomnumbers = random.sample(list(orthodict.keys()), SAMPLEsize)
	# randomnumbers = random.sample(range(1, len(orthodict.keys())), SAMPLEsize)
else: #Read the existing file
	randomnumbers = [int(line.rstrip("\n").split("\t")[0]) for line in open("concatenated/genes_selected.txt", 'r') if 'Number' not in line]
randomnumbers = [1,2,3]
samporthos = [] # Sampled orthologs
for randnum in randomnumbers:
	samporthos.append(orthodict[randnum])

randomfastas = [f'{OrthoTreesPath}/{ortho[0]}' for ortho in samporthos]

rule concatenation:
	input:
		fastas = randomfastas,
	output:
		fasta = "concatenated/onetoones.fa",
		random = "concatenated/genes_selected.txt"
	log: "logs/concatenation.log"
	run:# similar to script fastaconcat.py
		# Make a dictionary of the sequences, starting with empty ones
		masternames = {} # A list of all names in all files
		for sample in AllSamples:
			masternames[sample] = Seq("")

		logfile = open(log[0], 'w')
		ofile = open(output.random, 'w')
		ofile.write(f"Number\tOrthogroup\tREFgene\tAlignment\n")
		for randnum in randomnumbers:
			ofile.write(f"{randnum}\t{orthodict[randnum][1]}\t{orthodict[randnum][2]}\t{orthodict[randnum][0]}\n")

		for fastafile in input.fastas:
			# Start the concatenated alignment
			thisfasta = [seq_record for seq_record in SeqIO.parse(fastafile, "fasta")]
			thislen = len(thisfasta[0].seq) # len of alignment in case I need to add missing data
			
			for name in masternames.keys():
				thisnameseq = "-"*thislen # assume it's not there

				for seq in thisfasta:
					if name in seq.id: # it is there!
						thisnameseq = seq.seq
						break # you found it, so stop looping
				# Concat the sequence
				masternames[name] = masternames[name] + thisnameseq	

				if thisnameseq == "-"*thislen:
					logfile.write(f"Missing data in {fastafile} in sample {name}.\n")
		
		# Print concatenated sequence
		result = open(output.fasta, 'w')
		for sample in masternames.keys():
			result.write(">" + sample + '\n')
			result.write( str(masternames[sample]) + '\n')

rule IQTreeConcat:
	""" Make a tree of the concatenated alignment """
	input:
		"concatenated/onetoones.fa"
	output:
		"concatenated/onetoones.treefile"
	params:
		bootstraps = 1000, # UFBoot
	threads: 10
	resources:
		threads = lambda wildcards, threads: threads,
		mem_mb = lambda wildcards, threads: 6800 * threads, # each core gets at most 6.8 GB of RAM in Rackham
		time =  "2-00:00:00",
	shell:
		"iqtree -s {input} -m MFP -seed 1234 -bb {params.bootstraps} -nt {threads} -bnni -pre 'concatenated/onetoones'"

# ------ Trees for ASTRAL -------

rule renametrees:
	""" Rename the leaves in the trees to be all the same """
	input:
		trees = treesfile
	output:
		trees = "trees/alltrees.nom.tre"
	log: "logs/renametrees.log"
	run:
		newtrees = []
		logfile = open(log[0], 'w')

		count = 1
		for line in open(input.trees, 'r'):
			t = Tree(line)

			# Rename the leaves
			present_samples = []
			for node in t.iter_leaves():
				for sample in AllSamples:
					if sample in node.name:
						node.name = sample
						present_samples.append(sample)
						break
			# Save
			if set(present_samples) == set(AllSamples):
				newtrees += [t.write()]
			else:
				logfile.write(f"Tree in line {count} has an incomplete set of taxa.\n")

			count += 1

		# Print trees in a new file
		namedtrees = open(output.trees, 'w')
		for tree in newtrees:
			namedtrees.write(tree + "\n")

rule report_genes_astral:
	""" Report what genes were used for the ASTRAL analysis """
	input:
		badgenes = "logs/renametrees.log",
		allgenes = orthoreport, 
	output:
		report = "results/OneToOne_astral.txt"
	run:
		# Start the output file
		ofile = open(output.report, 'w')
		ofile.write(f"Number\tOrthogroup\tREFgene\tAlignment\n") # header

		# Get the line number of the trees that didn't have all strains
		bad_numbers = [int(line.rstrip("\n").split(" ")[3]) for line in open(input.badgenes, 'r')]
		
		for num in orthodict.keys():
			if num not in bad_numbers:
				ofile.write(f"{num}\t{orthodict[num][1]}\t{orthodict[num][2]}\t{orthodict[num][0]}\n")


rule collapselowsupportbraches:
	""" Remove nodes in the tree that have low support using newick utilities """
	input: 
		"trees/alltrees.nom.tre"
	output:
		"astral/alltrees.nom.bb{MINsupport}.tre"
	shell:
		"nw_ed {input} 'i & b<={wildcards.MINsupport}' o > {output}"
		# i matches internal nodes
		# b > 75 matches nodes whose label has a numerical value of 75 or more (if the label is numeric)
		# o   (splice Out) splice out node, and attach children to parent,
		# 		preserving branch lengths. This is useful for "opening" poorly
		#		supported nodes.

rule astral:
	""" Run ASTRAL on the collapsed trees """
	input:
		"astral/alltrees.nom.bb{MINsupport}.tre"
	output:
		"astral/astral.bb{MINsupport}.tre"
	wildcard_constraints:
		MINsupport="\d+"
	log: "astral/astral.bb{MINsupport}.log",
	params:
		astral = astral,
		# mappingfile = mappingfile, # what sample belongs to what species
	shell:
		"java -jar {params.astral} -i {input} -o {output} 2> {log}"
		# "java -jar {params.astral} -i {input} -a {params.mappingfile} -o {output} 2> {log}"

rule scoring_astral:
	""" Produce extra branch support metrics for the Species tree of ASTRAL """
	input:
		genetrees = "astral/alltrees.nom.bb{MINsupport}.tre",
		spptree = "astral/astral.bb{MINsupport}.tre"
	output:
		"astral/astral.bb{MINsupport}.scored.tre"
	wildcard_constraints:
		MINsupport="\d+"
	log: "astral/astral.bb{MINsupport}.scored.log"
	params:
		astral = astral,
	threads: 1
	resources:
		threads = lambda wildcards, threads: threads,
		mem_mb = lambda wildcards, threads: 6800 * threads, # each core gets at most 6.8 GB of RAM in Rackham
		time =  "2:00:00",
	shell:
		"java -jar {params.astral} -q {input.spptree} -i {input.genetrees} -o {output} -t 8 2> {log}"
		# -t 2 Full annotation (see https://github.com/smirarab/ASTRAL/blob/master/astral-tutorial.md#extensive-branch-annotations)
		# -t 8 Alternative quartet topologies Outputs q1, q2, q3; these three values show quartet support (as defined in the description of -t 1) for the main topology, the first alternative, and the second alternative, respectively.
		#		Main topology: RL|SO, First alternative: RS|LO, and Second alternative: RO|LS

rule polytomy_test_astral:
	""" Test the null hypothesis of polytomy (see doi:10.3390/genes9030132) """
	input:
		genetrees = "astral/alltrees.nom.bb{MINsupport}.tre",
		spptree = "astral/astral.bb{MINsupport}.tre"
	output:
		"astral/astral.bb{MINsupport}.poly.tre"
	wildcard_constraints:
		MINsupport="\d+"
	log: "astral/astral.bb{MINsupport}.poly.log",
	params:
		astral = astral,
	shell:
		"java -jar {params.astral} -q {input.spptree} -i {input.genetrees} -o {output} -t 10 2> {log}"
		# -t 2 Full annotation (see https://github.com/smirarab/ASTRAL/blob/master/astral-tutorial.md#extensive-branch-annotations)

# ------ Various metrics of certainty -------

rule InternodeCertainty:
	""" Calculate the internode Certainty of the ML tree """
	input:
		concat = "concatenated/onetoones.treefile", # RAxML expects to read a strictly bifurcating tree
		trees = "trees/alltrees.nom.tre"
	output:
		"results/RAxML_info.onetoones.tre",
		"results/RAxML_IC_Score_BranchLabels.onetoones.tre"
	shell:
		"raxmlHPC -f i -t {input.concat} -z {input.trees} -m GTRCAT -w $PWD/results -n onetoones.tre"

rule QuartetScores_concat:
	""" Calculate quartet internode certainty scores (Zhou et al. 2020 Systematic Biology) """
	input:
		reftree = "concatenated/onetoones.treefile",
		evaltrees = "astral/alltrees.nom.bb{MINsupport}.tre"
	output:
		newick = "results/Concatenated1000.bb{MINsupport}_QuartetScores.tre"
	params:
		quartetscores = quartetscores
	log:
		"logs/QuartetScores/Concatenated1000.bb{MINsupport}_QuartetScores.log"
	threads: 3
	resources:
		threads = lambda wildcards, threads: threads,
		mem_mb = lambda wildcards, threads: 6800 * threads, # each core gets at most 6.8 GB of RAM in Rackham
		time =  "30:00",
	shell:
		"{params.quartetscores} -v -t {threads} -r {input.reftree} -e {input.evaltrees} -o {output.newick} > {log}"

# -r <file_path>, --ref <file_path>: (required) Path to the reference tree
# -e <file_path>, --eval <file_path>: (required) Path to the evaluation trees
# -o <file_path>, --output <file_path>: (required) Path to the annotated newick output file (for LQIC/QPIC/EQPIC scores)
# -q <file_path>, --qic <file_path>: (optional) Path to the file where to write the raw QIC scores for each quartet
# -s, --savemem: Consume less memory, but with the cost of increased runtime (~50% more)
# -v, --verbose: Verbose mode
# -t <number>, --threads <number>: Maximum number of threads to use
# --version: Displays version information and exits.
# -h, --help: Displays usage information and exits.

rule unroot_astral:
	""" To make it comparable with the Concatenated tree """
	input:
		tree = "astral/astral.bb{MINsupport}.tre",
	output:
		tree = "astral/astral.bb{MINsupport}.unroot.tre",
	run:
		treefile = open(input.tree, 'r').readlines()
		
		# Unroot
		t = Tree(treefile[0], format = 0)
		t.unroot()

		# Print tree in a new file
		rootedtrees = open(output.tree, 'w')
		rootedtrees.write(t.write() + "\n")

rule QuartetScores_astral:
	""" Calculate quartet internode certainty scores (Zhou et al. 2020 Systematic Biology) """
	input:
		reftree = "astral/astral.bb{MINsupport}.unroot.tre",
		evaltrees = "astral/alltrees.nom.bb{MINsupport}.tre"
	output:
		newick = "results/AstralAll.bb{MINsupport}_QuartetScores.tre"
	params:
		quartetscores = quartetscores
	log:
		"logs/QuartetScores/AstralAll.bb{MINsupport}_QuartetScores.log"
	threads: 3
	resources:
		threads = lambda wildcards, threads: threads,
		mem_mb = lambda wildcards, threads: 6800 * threads, # each core gets at most 6.8 GB of RAM in Rackham
		time =  "30:00",
	shell:
		"{params.quartetscores} -v -t {threads} -r {input.reftree} -e {input.evaltrees} -o {output.newick} > {log}"


