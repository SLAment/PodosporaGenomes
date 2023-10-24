# -*- snakemake -*-

from glob import glob

### OrthoTreesUppmax: Getting orthologs of the Podospora complex - Uppmax version
#############################################################################

# The original pipeline was designed in the Uppsala Johannesson Server back in
# 2019, but that one no longer exist. Now I'll migrate the whole analysis to
# Uppmax (in 2023).

# A general pipeline to get orthologous genes from the *Podospora* complex.
# Due to the unpredictable number of ortholog groups, the pipeline relies on a
# [checkpoint function]. I ran into a bug and I had to do some very awkward
# work around. 

# The second version was an attempt to make it faster but it somehow gets
# stuck re-calculating the graph or something and it takes FOREVER, even more
# than the original strategy where I was extracting/BLASTing all genes from
# all samples in a single rule. In the second version I was trying a single
# rule for each blast of each strain. It sounds good for parallelizing, but it
# didn't work out and it turned out to be even slower.

# In the third version I tried to reach a compromise: a rule per sample, for
# all genes. Hopefully it's faster than the first version, without so much of
# the graph problem of the second.

#############################################################################
# ===========================================================================
# Sandra Lorena Ament Velasquez
# 2019/12/19 - 2023/08/15
# ++++++++++++++++++++++++++++++++++++++++++++++
# Version 5

# -------------------------------------------------
# DATA from the configuration file
configfile: "config/OrthoTreesUppmax_config.yaml"

samples = config["SampleIDs"]
SpeciesTreeSamples = config["SpeciesTreeSamples"]

assembliespath = config["assembliespath"]
annotationpath = config["annotationpath"]

# References
PODAN = config["PODAN"]
PODANgff = config["PODANgff"]
PODANgenes = config["PODANgenes"]
PODCO = config["PODCO"]
PODCOgff = config["PODCOgff"]

AllSamples = samples + ["PODAN", "PODCO"]
# -------------------------------------------------

# ----------
# Rules not submitted to a job
localrules: get_gffutils2fasta, get_orthogrs_parser, get_query2haplotype, get_fastaconcat, getgenomes, getrefgenomes, getprotsamples, getprotpodan, getprotcomata, parseorthogroups1n, orthofolders, getPODANgenes, aggregate_trees, make_report
# ----------

rule all:
	input:
		"results/SingleGeneTrees.tre",
		"results/OneToOne_equivalences.txt",
		expand("filtering/{orthodir}/PODAN_1n.clean.txt", orthodir = ["AllSamples", "OnePerSpecies"]) # just so I get the AllSamples


# ------- PREPARE SCRIPTS --------

rule get_gffutils2fasta:
	""" Download script from my GitHub """
	output:
		"scripts/gffutils2fasta.py"
	shell:
		"wget -O {output} https://raw.githubusercontent.com/SLAment/Genomics/master/GenomeAnnotation/gffutils2fasta.py"

rule get_orthogrs_parser:
	""" Download script from my GitHub """
	output:
		"scripts/orthogrs_parser.py"
	shell:
		"wget -O {output} https://raw.githubusercontent.com/SLAment/Genomics/master/Phylogenetics/orthogrs_parser.py"

rule get_query2haplotype:
	""" Download script from my GitHub """
	output:
		"scripts/query2haplotype.py"
	shell:
		"wget -O {output} https://raw.githubusercontent.com/SLAment/Genomics/master/BLAST/query2haplotype.py"

rule get_fastaconcat:
	""" Download script from my GitHub """
	output:
		"scripts/fastaconcat.py"
	shell:
		"wget -O {output} https://raw.githubusercontent.com/SLAment/Genomics/master/FastaManipulation/fastaconcat.py"

# ------- PREPARE ALL DATA --------
## Make symlink of the genomes to work more easily
rule getgenomes:
	""" Make links to the assemblies """
	input: 
		assembliespath + "/{sample}.nice.fa",
	output:
		"genomes/{sample}.fa"
	shell:
		"ln -sf {input} {output}"

rule getrefgenomes:
	""" Make links to the assemblies """
	input: 
		PODAN = PODAN,
		PODCO = PODCO,
	output:
		PODAN = "genomes/PODAN.fa",
		PODCO = "genomes/PODCO.fa",
	shell:
		"cat {input.PODAN} | sed 's;1plus;1;' | sed 's;Chr;chromosome_;' | sed 's;>;>PODAN_;' > {output.PODAN};"
		"cat {input.PODCO} | sed 's;Chromosome;chromosome;' | sed 's;>;>PODCO_;' > {output.PODCO};"

# ----------------------------------
rule getprotsamples:
	""" Get CDS sequences """ 
	input:
		genome = assembliespath + "/{sample}.nice.fa", # This should have been "genomes/{sample}.fa"
		gff = lambda wildcards: glob(annotationpath + "/{sample}.nice*.gff3".format(sample = wildcards.sample)), # Dirty trick to expand * like in bash
		gffutils2fasta = "scripts/gffutils2fasta.py"
	output:
		prots = "proteins/{proteinsdir}/{sample}.fas",
	conda:
		"envs/gffutils.yaml"
	shell:
		"""
		# Get CDS translated
		python {input.gffutils2fasta} {input.genome} {input.gff} --output {output.prots} -t CDS -p -j --onlyids
	
		# Put the sample ID in the name of the sequences
		sed -i 's;>\\(.*\\);>\\1_{wildcards.sample};' {output.prots} 
		"""

rule getprotpodan:
	""" Get CDS sequences of reference genomes""" 
	input:
		genome = PODAN,
		gff = PODANgff,
		gffutils2fasta = "scripts/gffutils2fasta.py"
	output:
		prots = "proteins/{proteinsdir}/PODAN.fas",
	conda:
		"envs/gffutils.yaml"
	shell:
		"""
		# Get CDS translated
		python {input.gffutils2fasta} {input.genome} {input.gff} --output {output.prots} -t CDS -p -j --onlyids
	
		# Put the sample ID in the name of the sequences
		sed -i 's;>gene:\\(.*\\).G;>\\1;' {output.prots}
		"""

rule getprotcomata:
	""" Get CDS sequences of reference genomes""" 
	input:
		genome = PODCO,
		gff = PODCOgff,
		gffutils2fasta = "scripts/gffutils2fasta.py"
	output:
		prots = "proteins/{proteinsdir}/PODCO.fas",
	conda:
		"envs/gffutils.yaml"
	shell:
		"""
		# Get CDS translated
		python {input.gffutils2fasta} {input.genome} {input.gff} --output {output.prots} -t CDS -p -j --onlynames
	
		# Put the sample ID in the name of the sequences
		sed -i 's;>\\(.*\\);>\\1;' {output.prots} 
		"""

# ------- Run OrthoFinder --------
# Do the OrthoFinder analysis but only with one strain per species
rule orthofinder_representatives:
	""" Run OrthoFinder with just representatives """
	input:
		expand("proteins/OnePerSpecies/{sample}.fas", sample = SpeciesTreeSamples)
	output:
		"OrthoFinder/OnePerSpecies/Orthogroups/Orthogroups.tsv"
	threads: 20
	resources:
		threads = lambda wildcards, threads: threads,
		mem_mb = lambda wildcards, threads: 6800 * threads, # each core gets at most 6.8 GB of RAM in Rackham
		time =  "2:00:00",
	shell:
		"module load bioinfo-tools OrthoFinder/2.5.2; " # If I load it in the main environment it breaks the small conda envs
		"rm -r OrthoFinder/OnePerSpecies; " # OrthoFinder fails if the directory exists
		"orthofinder.py -f proteins/OnePerSpecies -t {threads} -o OrthoFinder/OnePerSpecies; "
		"mv OrthoFinder/OnePerSpecies/Results_*/* OrthoFinder/OnePerSpecies/ && rm -r OrthoFinder/OnePerSpecies/Results_*"

rule orthofinder_all:
	""" Run OrthoFinder for all samples """
	input:
		"OrthoFinder/OnePerSpecies/Orthogroups/Orthogroups.tsv",
		expand("proteins/AllSamples/{sample}.fas", sample = AllSamples)
		# expand("proteins/TheRest/{sample}.fas", sample = [x for x in AllSamples if x not in SpeciesTreeSamples])
	output:
		"OrthoFinder/AllSamples/Orthogroups/Orthogroups.tsv"
	threads: 20
	resources:
		threads = lambda wildcards, threads: threads,
		mem_mb = lambda wildcards, threads: 6800 * threads, # each core gets at most 6.8 GB of RAM in Rackham
		time =  "10:00:00",
	shell:
		"module load bioinfo-tools OrthoFinder/2.5.2; "
		"rm -r OrthoFinder/AllSamples; " # OrthoFinder fails if the directory exists
		"orthofinder.py -f proteins/AllSamples -t {threads} -o OrthoFinder/AllSamples; "
		"mv OrthoFinder/AllSamples/Results_*/* OrthoFinder/AllSamples/ && rm -r OrthoFinder/AllSamples/Results_*"

rule parseorthogroups1n:
	""" Parse the output of orthogroups """
	input:
		orthos = "OrthoFinder/{orthodir}/Orthogroups/Orthogroups.tsv",
		orthogrs_parser = "scripts/orthogrs_parser.py"
	log:
		"filtering/{orthodir}/orthogrs_parser.log"
	output:
		podan1n = "filtering/{orthodir}/PODAN_1n.txt",
		podan1nClean = "filtering/{orthodir}/PODAN_1n.clean.txt",
		podanorthos = "filtering/{orthodir}/orthogroups_1n.txt"
	shell:
		"""
		# Filter the Orthogroups.csv for groups of one-to-one orthologs present in all samples
		python {input.orthogrs_parser} {input.orthos} -n1 -b -o filtering/{wildcards.orthodir} --ref PODAN > {log}

		# Remove genes that are not starting with "Pa_"
		grep '^Pa' {output.podan1n} > {output.podan1nClean}
		"""

checkpoint orthofolders:
	""" Make dummy files for every orthogroup """
	input:
		orthologs = "filtering/OnePerSpecies/PODAN_1n.clean.txt",
	output:
		directory("fastahits")
	run:
		shell("mkdir -p fastahits")
		shell("mkdir -p fastahits/dummies")
		num_lines = sum(1 for line in open(input.orthologs)) # How many lines in the file?

		for n in range(1, num_lines + 1):
			number = "{0:04d}".format(n)
			tabopen = open(f"fastahits/dummies/orthologs{number}.dummy", 'w')

rule getPODANgenes:
	""" We need the nucleotide sequence of the reference genes (including introns) """
	input:
		genome = PODAN,
		gff = PODANgff,
		gffutils2fasta = "scripts/gffutils2fasta.py"
	output:
		genes = PODANgenes,
	conda:
		"envs/gffutils.yaml"
	shell:
		"""
		# Get CDS translated
		python {input.gffutils2fasta} {input.genome} {input.gff} --output {output.genes} -t gene --onlyids
	
		# Put the sample ID in the name of the sequences
		sed -i 's;>gene:\\(.*\\).G;>\\1;' {output.genes}
		"""

rule query2haplotype:
	""" Get fasta files of each ortholog and for each sample"""
	input:
		"fastahits",
		orthologs = "filtering/OnePerSpecies/PODAN_1n.clean.txt",
		query2haplo = "scripts/query2haplotype.py",
		refgenes = PODANgenes,
		genomes = expand("genomes/{sample}.fa", sample = AllSamples)
	output:
		"fastahits/orthologs{i}.fas"
	threads: 3
	resources:
		threads = lambda wildcards, threads: threads,
		mem_mb = lambda wildcards, threads: 6800 * threads, # each core gets at most 6.8 GB of RAM in Rackham
		time =  "2:00:00",
	run:
		# Read the ortholog groups file
		tabopen = open(input.orthologs, 'r')
		tabs = [line.rstrip("\n") for line in tabopen] 			# Read tab file into a list
		
		# What orthogroup is this?
		ortho = tabs[int(wildcards.i) - 1]

		# Make a temporary folder
		shell(f"mkdir -p fastahits/orthologs{{wildcards.i}}")

		for strain in AllSamples:
			# Get fasta
			cmd1 = f"python {input.query2haplo} genomes/{strain}.fa {input.refgenes} --seqid {ortho} -t {threads} --tophit --extrabp 0 --temp fastahits/orthologs{wildcards.i}/{strain}/ >> fastahits/orthologs{wildcards.i}/{strain}.temp.fa"
			shell(cmd1)

		# Put them together
		shell(f"cat fastahits/orthologs{{wildcards.i}}/*.temp.fa > {{output}}")

		# Clean
		shell(f"rm -r fastahits/orthologs{{wildcards.i}}")

def individualOrthos(wildcards):
	# From Dima's tutorial
	# https://evodify.com/snakemake-checkpoint-tutorial/
	'''
	Aggregate the file names of the random number of files
	generated at the paralogs step
	'''
	checkpoint_output = checkpoints.orthofolders.get(**wildcards).output[0] # The output folder's name
	return expand("fastahits/{sample}/{sample}_ortholog{i}.fas",
		i = glob_wildcards(os.path.join(checkpoint_output, 'dummies/orthologs{i}.dummy')).i, sample = wildcards.sample)

rule mafft:
	""" Align ortholog groups with MAFFT """
	input:
		ortholog = "fastahits/orthologs{i}.fas",
	output:
		"alignments/orthologs{i}.fas",
	threads: 15 # one long gene (Pa_4_4640) is asking for a lot (6 were fine for all the other genes)
	resources:
		threads = lambda wildcards, threads: threads,
		mem_mb = lambda wildcards, threads: 6800 * threads, # each core gets at most 6.8 GB of RAM in Rackham
		time =  "3:00:00",
	shell:
		"mafft --thread {threads} --threadit 0 --adjustdirection --anysymbol --maxiterate 1000 --retree 1 --localpair {input.ortholog} > {output}"

rule IQTreePerGene: # Careful, some trees don't have all the samples
	""" Run IQTree for each ortholog group """
	input:
		"alignments/orthologs{i}.fas"
	output:
		"iqtree/orthologs{i}/orthologs{i}.treefile"
	threads: 2
	resources:
		threads = lambda wildcards, threads: threads,
		mem_mb = lambda wildcards, threads: 6800 * threads, # each core gets at most 6.8 GB of RAM in Rackham
		time =  "3:00:00",
	params:
		bootstraps = 1000, # UFBoot
	conda:
		"envs/iqtree.yaml"
	shell:
		"""
		iqtree -s {input} -m MFP -seed 1234 -bb {params.bootstraps} -nt {threads} -bnni -pre "iqtree/orthologs{wildcards.i}/orthologs{wildcards.i}" --keep-ident
		"""
		# -bnni to reduce the risk of overestimating branch supports with UFBoot due to severe model violations
		# https://github.com/Cibiv/IQ-TREE/issues/196

def iqtreeoutput(wildcards):
	# From Dima's tutorial
	# https://evodify.com/snakemake-checkpoint-tutorial/
	'''
	Aggregate the file names of the random number of files
	generated at the paralogs step
	'''
	checkpoint_output = checkpoints.orthofolders.get(**wildcards).output[0] # The output folder's name
	return expand("iqtree/orthologs{i}/orthologs{i}.treefile",
		i = glob_wildcards(os.path.join(checkpoint_output, 'dummies/orthologs{i}.dummy')).i)

rule aggregate_trees:
	""" Collect all trees in a single file """ 
	# With this rule I trigger the formation of all other checkpoint files
	input:
		iqtreeoutput,
	output:
		"results/SingleGeneTrees.tre"
	# shell:
	# 	"cat {input} > {output}"
	shell:
		"cat iqtree/orthologs*/*.treefile > {output}"
	# run: # I had to do it like this because the line becomes so long that i get an error ([Errno 7] Argument list too long)
	# 	# Save all trees in a single file
	# 	alltrees = []
	# 	for tree in input:
	# 		# if 'orthologs2437' not in tree: # This is a stupid fix for that ortholog being perfectly conserved in the complex
	# 		t = Tree(tree)
	# 		alltrees += [t.write()] # Save the newick line as is

	# 	# Save all trees in a single file
	# 	result = open(output[0], 'w')
	# 	for tree in alltrees:
	# 		result.write(tree + "\n")


def mafftout(wildcards):
	# From Dima's tutorial
	# https://evodify.com/snakemake-checkpoint-tutorial/
	'''
	Aggregate the file names of the random number of files
	generated at the paralogs step
	'''
	checkpoint_output = checkpoints.orthofolders.get(**wildcards).output[0] # The output folder's name
	return expand("alignments/orthologs{i}.fas",
		i = glob_wildcards(os.path.join(checkpoint_output, 'dummies/orthologs{i}.dummy')).i)

rule make_report:
	""" Make a table with the names equivalence """
	input:
		alignments = mafftout,
		podan1nClean = "filtering/OnePerSpecies/PODAN_1n.clean.txt",
		podanorthos = "filtering/OnePerSpecies/orthogroups_1n.txt"
	output:
		"results/OneToOne_equivalences.txt"
	shell:
		"ls -1 {input.alignments} > results/alignments_names.txt; "
		"paste results/alignments_names.txt {input.podanorthos} {input.podan1nClean} | awk '{{print NR, $0}}' > {output}; "
		"rm results/alignments_names.txt"


