# -*- snakemake -*-

### gff2function.smk - A pipeline to do functional annotation of modified MAKER gffs
#############################################################################

# This pipeline is meant to be ran after PaAnnotation.smk. Unfortunately,
# funannotate scripts remove away a lot of the information I worked hard to
# introduce in PaAnnotation.smk. Here I attempt to recover it, but not
# completely :(.

# https://www.ncbi.nlm.nih.gov/genbank/eukaryotic_genome_submission_annotation/
#############################################################################
# ===========================================================================
# Sandra Lorena Ament Velasquez
# 2023/05/02
# ++++++++++++++++++++++++++++++++++++++++++++++
# Version 1

# -------------------------------------------------
###  CONFIGURATIONS
# -------------------------------------------------
# DATA from the configuration file
configfile: "config/config.yaml"

# Paths and files
path2assemblies = config["path2assemblies"]
path2samples = config["path2samples"]
path2orthos = config["path2orthos"]
knownproducts = config["knownproducts"]
pseudogenes = config["pseudogenes"]

# Scripts
decoratePodoGff3 = config["decoratePodoGff3"]
mergeGFFs = config["mergeGFFs"]
appendQualifierGff = config["appendQualifierGff"]
fixProductGff = config["fixProductGff"]
tbl4ncbi = config["tbl4ncbi"]
addtRNA2tbl = config["addtRNA2tbl"]

annotationversion = config["annotationversion"]

# Samples
longassemblies = config["longassemblies"]
shortassemlies = config["shortassemlies"]

samples = longassemblies + shortassemlies
samples = "Chr4"
# ----------
# Rules not submitted to a job
localrules: getPreniceGFF, newnames_gff, removesmallctgs, fix_product_gff, gff2tbl_final, get_gffutils2fasta, decoratePodoGff3, recoverlostgenes, appendqualifiers, sortgff, newnames, get_oldnames, getproteins, tbl4ncbi, addtRNA_to_tbl, single_reports, full_report
# ----------

# -------------------------------------------------
from Bio import SeqIO
import os.path

# ---------------------------
### Useful functions
# ---------------------------

def assignLocusTag(sample):
	if sample == "CBS237.71m":
		return("EYR66_")
	elif sample == "CBS124.78p":
		return("QC764_")
	elif sample == "CBS411.78m":
		return("QC763_")
	elif sample == "CBS415.72m":
		return("QC762_")
	elif sample == "CBS112042p":
		return("QC761_")
	else:
		return("FUN_") # something generic since I won't put this in NCBI

def assignSpp(sample):
	if "Pa" in sample or (sample in ["CBS433.50p", "CBS455.64m", "Chr4"]):
		return("Podospora anserina")
	elif "Pc" in sample:
		return("Podospora comata")
	elif sample in ["CBS333.63p", "CBS237.71m", "CBS451.62p"]:
		return("Podospora pauciseta")
	elif sample in ["CBS124.78p", "CBS253.71p"]:
		return("Podospora pseudoanserina")
	elif sample in ["CBS411.78m"]:
		return("Podospora pseudopauciseta")
	elif sample in ["CBS415.72m"]:
		return("Podospora pseudocomata")
	elif sample in ["CBS112042p"]:
		return("Podospora bellae-mahoneyi")

def gettemplate(wildcards):
	if wildcards.sample in ["CBS124.78p", "CBS411.78m", "CBS415.72m", "CBS112042p"]:
		return(f"data/GenBankSubTemplate/{wildcards.sample}.sbt")
	else:
		return("data/GenBankSubTemplate/dummy.sbt") # A fake file for all the other samples


# ---------------------------

rule all:
	input:
		"reports/AllSamples.nice-" + annotationversion + ".report.txt",
		expand("results/{sample}/{sample}.nice-" + annotationversion + ".tbl", sample = samples),
		expand("results/{sample}/{sample}.nice-" + annotationversion + ".gff3", sample = samples),
		expand("InterProScan/{sample}.xml", sample = samples),

# ---------------------------

# ---------------------------
# Things break with very long names
# ---------------------------
rule newnames:
	""" Make a new fasta temporary names """
	input:
		genome = path2assemblies + "/{sample}.fa",
	output:
		tempnames = "temp/genomes/{sample}_tempnames.txt", # To keep the original names somewhere
		newgenome = "temp/genomes/{sample}_tempnames.fa"
	run:
		count = 1
		countmt = 1
		dictnames = {} # Setup a dictionary
		output_handle = open(output.newgenome, "w")

		for seq_record in SeqIO.parse(input.genome, "fasta"):
			if seq_record.id in dictnames.keys(): # Just in case
				print(f"ERROR: The sequence {seq_record.id} is more than once in the fasta file {input.genome}")
				raise ValueError(f"ERROR: The sequence {seq_record.id} is more than once in the fasta file {input.genome}")
			else:
				if wildcards.sample in shortassemlies:
					newname = "seq{0:06d}".format(count)
					count += 1 # Increase the count for the next sequence
				elif 'chromosome' in seq_record.id:
					newname = seq_record.id.replace("chromosome_", "")
				elif len(seq_record.id) > 16: # NCBI limits the number of characters in a FASTA header for submission to 16 characters and Augustus also has problems with longer contig/scaffold names.
					if '_mt' in seq_record.id:
						newname = wildcards.sample + "mt_{0:01d}".format(countmt) # There shouldn't be more than 1 but just in case
						countmt += 1
					else:
						newname = wildcards.sample + "_u{0:02d}".format(count) # I don't expect many unplaced contigs
						count += 1
				else: # leave as is
					newname = seq_record.id

				dictnames[seq_record.id] = newname # record the original name somewhere
				seq_record.id = newname
				seq_record.description = '' #Annoying extra info
				SeqIO.write(seq_record, output_handle, "fasta")

		# Keep the names somewhere
		with open(output.tempnames, 'w') as result:
			for key, value in dictnames.items():
				result.write(f"{key}\t{value}\n")

def getPrenice(wildcards):
	highfile = f"{path2samples}/{wildcards.sample}/{wildcards.sample}.prenice-{float(annotationversion) + 0.01}.gff3"
	if os.path.isfile(highfile): # Check if it exists
		return(highfile)
	else:
		return(f"{path2samples}/{wildcards.sample}/{wildcards.sample}.prenice-{annotationversion}.gff3")

rule getPreniceGFF:
	""" Make a symlink to the data """ 
	# This is also a trick to force different versions of the prenice gff to be used, not just 3.00
	input:
		getPrenice
	output:
		temp("data/raw_data/{sample}-prenice.gff3")
	shell:
		"ln -s {input} {output}"

rule newnames_gff:
	""" Get new names into the gff too """
	input:
		tempnames = "temp/genomes/{sample}_tempnames.txt",
		gff = "data/raw_data/{sample}-prenice.gff3"
		# gff = path2samples + "/{sample}/{sample}.prenice-3.00.gff3",
	output:
		gff = "data/raw_data/{sample}.gff3"
	run:
		# Make a dictionary with the contig names
		ctgdic = {line.rstrip("\n").split("\t")[0]: line.rstrip("\n").split("\t")[1] for line in open(input.tempnames, 'r')}

		# Change gff
		with open(output.gff, 'w') as result:
			# Replace this names in the gff file
			for line in open(input.gff, 'r'):
				if '#' in line:
					result.write(line)
				elif line == '\n':
					pass
				else:
					cols = line.rstrip("\n").split("\t")
					contig = cols[0]
					newline = ctgdic[contig] + '\t' + '\t'.join(cols[1:]) + '\n'
					result.write(newline)

# ---------------------------
# antiSMASH
# ---------------------------

rule removesmallctgs:
	""" antiSMASH complains with small contigs that have no gene models """
	input:
		genome = "temp/genomes/{sample}_tempnames.fa",
		gff = "data/raw_data/{sample}.gff3"
	output:
		genome = temp("temp/genomes/{sample}_tempnames_big.fa")
	run:
		# Read gff and make sure there are models for the contigs
		ctglist = []
		for line in open(input.gff):
			cols = line.rstrip("\n").split("\t")

			# This is a SUPER LAZY catch of an error with the PaYm genome that
			# only has tRNAs predicted in its mitochondria; because of that,
			# gff2tbl doesn't print them and that contig ends up without
			# annotation, leading to an error with antiSMASH
			if cols[0] != 'PaYm_mt': 
				ctglist.append(cols[0])
		ctglist = set(ctglist)

		# Read fasta
		output_handle = open(output.genome, "w")
		for seq_record in SeqIO.parse(input.genome, "fasta"):
			if seq_record.id in ctglist:
				SeqIO.write(seq_record, output_handle, "fasta")

rule gff2tbl:
	input:
		genome = "temp/genomes/{sample}_tempnames_big.fa",
		gff = "data/raw_data/{sample}.gff3"
	output:
		"data/raw_data/{sample}.tbl"
	shell:
		"funannotate util gff2tbl --gff3 {input.gff} --fasta {input.genome} > {output}"

rule tbl2gbk:
	input:
		tbl = "data/raw_data/{sample}.tbl",
		genome = "temp/genomes/{sample}_tempnames_big.fa",
	output:
		"data/raw_data/{sample}.gbk"
	run:
		shell(f"funannotate util tbl2gbk --tbl {input.tbl} --fasta {input.genome} --species '{assignSpp(wildcards.sample)}' --strain {wildcards.sample} --output data/raw_data/{wildcards.sample}")

rule set_antiSMASH_env: 
	""" Create a conda enviroment and prepare the databases """
	# https://docs.antismash.secondarymetabolites.org/install/
	output:
		"logs/antishmash_success.txt"
	conda:
		"envs/antismash.yaml"
	resources:
		time = "2:00:00", # Took FOREVER, so it's better in a job
	shell:
		"download-antismash-databases && touch {output}"

rule antiSMASH:
	input:
		env = "logs/antishmash_success.txt",
		gbk = "data/raw_data/{sample}.gbk",
	output:
		"antiSMASH/{sample}/{sample}.gbk"
	resources:
		threads = 1,
		mem_mb = int(1 * 6800), # each core gets at most 6.8 GB of RAM in Rackham
		time = "5:00:00",	
	conda:
		"envs/antismash.yaml"	
	shell:
		"""
		if [ -d antiSMASH/{wildcards.sample} ]; then
			rm -r antiSMASH/{wildcards.sample} # Otherwise it fails
		fi

		antismash {input.gbk} --taxon fungi --output-dir antiSMASH/{wildcards.sample} --genefinding-tool none
		# --genefinding-tool none -- without this, it will raise an error if contigs without genes are found.
		"""

# ---------------------------
# InterProScan
# ---------------------------

rule get_gffutils2fasta:
	""" Get the scripts needed for the rest of the pipeline """
	output:
		"scripts/gffutils2fasta.py",
	shell:
		"wget https://raw.githubusercontent.com/SLAment/Genomics/master/GenomeAnnotation/gffutils2fasta.py -P scripts"

rule getproteins:
	""" Extract proteins from the final gene models """
	input:
		genome = "temp/genomes/{sample}_tempnames.fa",
		gff = "data/raw_data/{sample}.gff3",
		gff3tofasta = "scripts/gffutils2fasta.py"
	output:
		fasta = "temp/proteins/{sample}.fas",
	conda:
		"envs/gffutils.yaml"
	shell:
		"python {input.gff3tofasta} {input.genome} {input.gff} --type CDS --join --proteinon --onlyids --mRNAids --output {output.fasta}; "
		"sed -i 's/*//' {output.fasta}" # InterProScan doesn't like the * and it's quiet about it if run within funannotate


rule InterProScan:
	""" Run the InterProScan module in Uppmax """ # module load bioinfo-tools InterProScan/5.52-86.0
	input:
		"temp/proteins/{sample}.fas",
	output:
		"InterProScan/{sample}.xml"
	threads: 15
	resources:
		threads = lambda wildcards, threads: threads,
		mem_mb = lambda wildcards, threads: 6800 * threads, # each core gets at most 6.8 GB of RAM in Rackham
		time = "1-00:00:00"
	shell:
		"module load bioinfo-tools InterProScan/5.62-94.0 SignalP/6.0g; "
		"funannotate iprscan -i {input} -m local --iprscan_path /sw/bioinfo/InterProScan/5.62-94.0/rackham/interproscan.sh -c {resources.threads} -o {output}"

# ---------------------------
# Funannotate annotate
# ---------------------------

rule funannotate_annotate:
	input:
		genome = "temp/genomes/{sample}_tempnames.fa",
		gff = "data/raw_data/{sample}.gff3",
		antismash = "antiSMASH/{sample}/{sample}.gbk",
		iprscan = "InterProScan/{sample}.xml",
		template = gettemplate
	output:
		"funannotate/{sample}/annotate_results/{species}_{sample}.gbk",
		"funannotate/{sample}/annotate_results/{species}_{sample}.gff3",
	threads: 10
	resources:
		threads = lambda wildcards, threads: threads,
		mem_mb = lambda wildcards, threads: 6800 * threads, # each core gets at most 6.8 GB of RAM in Rackham
		time = "5:00:00",
	run:
		# shell("module load bioinfo-tools SignalP/6.0g")
		if wildcards.sample in ["CBS124.78p", "CBS411.78m", "CBS415.72m", "CBS112042p"]:
			shell(f"funannotate annotate \
				--gff {input.gff} \
				--fasta {input.genome} \
				--out funannotate/{wildcards.sample} \
				--species '{assignSpp(wildcards.sample)}' \
				--strain {wildcards.sample} \
				--antismash {input.antismash} \
				--iprscan {input.iprscan} \
				--cpus {resources.threads} \
				--sbt {input.template}") # In the end I didn't really use it because I added my own annotation afterwards
		else:
			shell(f"funannotate annotate \
				--gff {input.gff} \
				--fasta {input.genome} \
				--out funannotate/{wildcards.sample} \
				--species '{assignSpp(wildcards.sample)}' \
				--strain {wildcards.sample} \
				--antismash {input.antismash} \
				--iprscan {input.iprscan} \
				--cpus {resources.threads}")

# ---------------------------
# Fine tunning of annotation
# ---------------------------

def getFunGFF(wildcards):
	return(f"funannotate/{wildcards.sample}/annotate_results/{assignSpp(wildcards.sample)}_{wildcards.sample}.gff3".replace(" ", "_"))

rule decoratePodoGff3:
	""" Add the name and known products of Podans genes """
	input:
		decoratePodoGff3 = decoratePodoGff3,
		onetoone = path2orthos + "/{sample}/onetoones.txt",
		annotation = knownproducts,
		gff = getFunGFF
		# gff = "funannotate/{sample}/annotate_results/{sample}_OGnames.gff3",
	output:
		"finetuning/{sample}/{sample}_prods.gff3"
	conda:
		"envs/gffutils.yaml"
	shell:
		"python {input.decoratePodoGff3} {input.gff} {input.onetoone} -g {input.annotation} > {output}"

# Multiple note qualifiers can be included and will be concatenated by
# table2asn or Genome Workbench into a single note with semi-colons as
# separators.

rule recoverlostgenes:
	""" Some genes were lost during the annotation (tRNAs)""" # Currently not adding a product 
	# this will add the genes with the old contig names, which is fine for the final nice file, but not for the names.gff3
	# however, I couldn't be bothered to fix it
	input:
		gff1 = "finetuning/{sample}/{sample}_prods.gff3",
		gff2 = "data/raw_data/{sample}.gff3",
		mergeGFFs = mergeGFFs
	output:
		# temp("finetuning/{sample}/{sample}_all_unsort.gff3")
		full = "finetuning/{sample}/{sample}_all_unsort.gff3",
		missing = "finetuning/{sample}/{sample}_missing.gff3",
	conda:
		"envs/gffutils.yaml"
	shell:
		"python {input.mergeGFFs} {input.gff1} {input.gff2} -m {output.missing} > {output.full}"

rule appendqualifiers:
	""" Some genes are known to be pseudogenes """
	input:
		gff = "finetuning/{sample}/{sample}_all_unsort.gff3",
		pseudos = pseudogenes,
		appendQualifierGff = appendQualifierGff
	output:
		temp("finetuning/{sample}/{sample}_all_unsort_pseudos.gff3")
	shell:
		"python {input.appendQualifierGff} {input.gff} {input.pseudos} > {output}"

# https://www.ncbi.nlm.nih.gov/genbank/genomes_gff/
# [2] pseudogenes should be flagged with pseudogene=<TYPE> qualifier in column
# 9 on the gene feature and optionally on any child features. Further details
# about the TYPE values allowed for the pseudogene qualifier are available
# at: http://www.insdc.org/documents/pseudogene-qualifier-vocabulary .

# [3] annotate with pseudo=true any genes that are 'broken' but are not
# thought to be pseudogenes. These are genes that do not encode the expected
# translation, for example because of internal stop codons. These are often
# caused by problems with the sequence and/or assembly.


rule sortgff:
	""" Sort the final gene models """
	input:
		"finetuning/{sample}/{sample}_all_unsort_pseudos.gff3"
	output:
		temp("finetuning/{sample}/{sample}_all_unsort_pseudos_sorted.gff3")
	conda:
		"envs/igvtools.yaml"
	shell:
		"""
		igvtools sort {input} {output} 
		"""

rule fix_product_gff:
	""" Fix badly formatted Funannotate product names """
	input:
		gff = "finetuning/{sample}/{sample}_all_unsort_pseudos_sorted.gff3",
		pseudos = pseudogenes,
		fixProductGff = fixProductGff
	output:
		"results/{sample}/{sample}.nice-" + annotationversion + ".names.gff3"
	conda:
		"envs/gffutils.yaml"
	shell:
		"python {input.fixProductGff} {input.gff} -g {input.pseudos} > {output}"

rule gff2tbl_final:
	input:
		genome = "temp/genomes/{sample}_tempnames.fa",
		gff = "results/{sample}/{sample}.nice-" + annotationversion + ".names.gff3",
	output:
		"results/{sample}/{sample}.nice-" + annotationversion + "_raw.tbl",
		# temp("results/{sample}/{sample}.nice-" + annotationversion + "_raw.tbl"),
	shell:
		"funannotate util gff2tbl --gff3 {input.gff} --fasta {input.genome} > {output}"

rule tbl4ncbi:
	""" The output of gff2tbl is not perfect so try to improve it """
	input:
		tbl = "results/{sample}/{sample}.nice-" + annotationversion + "_raw.tbl",
		pseudos = pseudogenes,
		tbl4ncbi = tbl4ncbi
	output:
		temp("results/{sample}/{sample}.nice-" + annotationversion + "_raw_ncbi.tbl"),
	shell:
		"python {input.tbl4ncbi} {input.tbl} -g {input.pseudos} > {output}"

rule addtRNA_to_tbl:
	""" The output of gff2tbl also removes the tRNAs so put them back into the tbl """
	input:
		tbl = "results/{sample}/{sample}.nice-" + annotationversion + "_raw_ncbi.tbl",
		gff = "finetuning/{sample}/{sample}_missing.gff3",
		addtRNA2tbl = addtRNA2tbl
	output:
		"results/{sample}/{sample}.nice-" + annotationversion + ".tbl",
	shell:
		"python {input.addtRNA2tbl} {input.tbl} {input.gff} > {output}"

# ---------------------------
# Recover the old names
# ---------------------------

rule get_oldnames:
	""" Recover the old names in the latest gff """
	input:
		tempnames = "temp/genomes/{sample}_tempnames.txt",
		gff = "results/{sample}/{sample}.nice-" + annotationversion + ".names.gff3",
	output:
		gff = "results/{sample}/{sample}.nice-" + annotationversion + ".gff3",
		# gff = "funannotate/{sample}/annotate_results/{sample}_OGnames.gff3"
	run:
		# Make a dictionary with the contig names (but inverted compared to before)
		ctgdic = {line.rstrip("\n").split("\t")[1]: line.rstrip("\n").split("\t")[0] for line in open(input.tempnames, 'r')}

		# Change gff
		with open(output.gff, 'w') as result:
			# Replace this names in the gff file
			for line in open(input.gff, 'r'):
				if '##' in line:
					result.write(line)
				elif line == '\n':
					pass
				else:
					cols = line.rstrip("\n").split("\t")
					contig = cols[0]
					if contig in ctgdic.keys(): # in case features were added with the old names in the recoverlostgenes rule 
						newline = ctgdic[contig] + '\t' + '\t'.join(cols[1:]) + '\n'
						result.write(newline)
					else:
						# print(contig)
						result.write(line)

# ---------------------------
# Reports
# ---------------------------

rule single_reports:
	""" Produce a simple report about the results """
	input:
		"results/{sample}/{sample}.nice-" + annotationversion + ".gff3",
	output:
		temp("reports/{sample}.nice-" + annotationversion + ".report.txt")
	shell:
		"""
		protsn=$(grep -cP "\\tmRNA\\t" {input})
		tRNAs=$(grep -cP "\\ttRNA\\t" {input})
		# rRNAs=$(grep -cP "\\trRNA\\t" {input}) # grep will return an error if the output is 0
		echo "{wildcards.sample}\t${{protsn}}\t${{tRNAs}}" > {output}
		"""

rule full_report:
	input:
		expand("reports/{sample}.nice-" + annotationversion + ".report.txt", sample = samples),
	output:
		"reports/AllSamples.nice-" + annotationversion + ".report.txt"
	shell:
		"echo 'Sample\tProteins\ttRNAs' > {output}; "
		"cat {input} >> {output}"

