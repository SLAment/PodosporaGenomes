# -*- snakemake -*-

### CircosSppComplex: Make Circos plots for the *Podospora anserina* species complex. 
#############################################################################

# A pipeline to produce Circos plots for all seven chromosomes of the
# *Podospora anserina* species complex. The annotation files come from running
# `PaAnnotation.smk`, specifically the repeat content ("xxxxx.repeatmasker.gff3"), 
# where "xxxxx" is the sample's name.

#############################################################################
# ==================================================
# Sandra Lorena Ament Velasquez
# Johannesson Lab, Evolutionary Biology Center, Uppsala University, Sweden
# 2020/07/17 - 2020/07/21
# +++++++++++++++++++++++++++++++++++++++++++++++++

from glob import glob
from Bio import SeqIO
from Bio.Alphabet import generic_dna

# -------------------------------------------------
samples = config["SampleIDs"]
reference = config["reference"]
path2assemblies = config["path2assemblies"]
path2repeatmasking = config["path2repeatmasking"]
circosconfpath = config["circosconfpath"]

DotPrep = config["DotPrep"]

NoChrs = 7 # The number of chromosomes in the Podospora species
samplesNoRef = [sample for sample in samples if reference not in sample]
# -------------------------------------------------

rule all:
	input:
		expand("circos_{reference}-vs-{sample}/circos_{reference}-vs-{sample}.png", reference = reference, sample = samplesNoRef),

# ------------ Get data ------------

## Make a dictionary of the samples (key) and their assemblies (values)
assemblies = [(glob(path2assemblies + "/{sample}*.fa*".format(sample=sample))) for sample in samples]
ASSEMBLIESDIC = dict(zip(samples, assemblies))

rule prepareddata:
	""" Get the assemblies so we don't disturb the original files """
	input: 
		reference = lambda wildcards: ASSEMBLIESDIC[wildcards.reference][0],
		assembly = lambda wildcards: ASSEMBLIESDIC[wildcards.sample][0]
	output:
		temp("fastas/{reference}-vs-{sample}.fa")
	wildcard_constraints:
		sample="[a-zA-Z0-9.]+"
	run:
		refseqs = [] # Setup an empty list

		# First the reference
		for seq_record in SeqIO.parse(open(input.reference, "r"), "fasta"):
			if 'chr' in seq_record.id: # Make sure it's a chromosome 
				# seq_record.id = "Podan2_" + seq_record.id.replace("chromosome_", "") # I wanted to change the name of the chromosome here, but I changed my mind
				if wildcards.reference not in seq_record.id:
					seq_record.id = wildcards.reference + "_" + seq_record.id
					seq_record.description = ''
				refseqs.append(seq_record)		

		newfasta = open(output[0], "w") # Make the new fasta
		SeqIO.write(refseqs, newfasta, "fasta")

		## Now the sample
		sampleseqsdic = {} # Setup a dictionary
		contignames = []

		for seq_record in SeqIO.parse(open(input.assembly, "r"), "fasta"):
			if 'chr' in seq_record.id: # Make sure it's a chromosome 
				if wildcards.sample not in seq_record.id:
					seq_record.id = wildcards.sample + "_" + seq_record.id
					seq_record.description = ''
				sampleseqsdic[seq_record.id] = seq_record

		# # Read into memory
		# records_dict = SeqIO.to_dict(SeqIO.parse(open(input.assembly, "r"), "fasta", generic_dna))
		# contignames = [seq for seq in records_dict.keys() if "chr" in seq] # Make a list of the sequences
		# contignames.reverse() # Invert them so they face each other in the plot
		
		contignames = list(sampleseqsdic.keys())
		contignames.reverse() # Invert them so they face each other in the plot

		for seq in contignames:
			SeqIO.write(sampleseqsdic[seq], newfasta, "fasta") # Add them to the output fasta

# ------------ Karyotype ------------


rule indexchr:
	""" Make an index of each chromosome collection """
	input:
		"fastas/{reference}-vs-{sample}.fa"
	output:
		"fastas/{reference}-vs-{sample}.fa.fai"
	shell:
		"samtools faidx {input}"

rule makekaryotype:
	""" Produce a karyotype file for circos """
	input:
		"fastas/{reference}-vs-{sample}.fa.fai"
	output:
		"circos_{reference}-vs-{sample}/karyotype.txt"
	params:
		threads = 1,
	shell:
		"""
		## Change the label to make it shorter and color the contigs differently for the reference and query
		grep {wildcards.reference} {input} | awk '{{a=$1; c=gensub(/(.*)chromosome_([\.0-9]*)/, "\\\\2", "g", a); print "chr","-",$1,c,1,$2,"vdgrey"}}' > {output} 
		grep {wildcards.sample} {input} | awk '{{a=$1; c=gensub(/(.*)chromosome_([\.0-9]*)/, "\\\\2", "g", a); print "chr","-",$1,c,1,$2,"grey"}}' >> {output} 

		# chr - ID LABEL START END COLOR
		# g in gensub stands for global
	 """	

# ------------ Create links ------------

rule mummer:
	input:
		query = lambda wildcards: ASSEMBLIESDIC[wildcards.sample][0],
		reference = lambda wildcards: ASSEMBLIESDIC[wildcards.reference][0],
	output:
		delta = "mummer/{reference}-vs-{sample}.delta",
		deltafilter = "mummer/{reference}-vs-{sample}.filter",
		coords = "mummer/{reference}-vs-{sample}.coords",
		coordsfilter = "mummer/{reference}-vs-{sample}.filter.coords",	
	params:
		threads = 8,
	shell:
		"""
		echo
		echo "MUMmer alignment ..."
		nucmer -b 2000 -c 2000 --maxmatch -p mummer/{wildcards.reference}-vs-{wildcards.sample} {input.reference} {input.query} -t {params.threads}

		# Filter the delta
		delta-filter -q {output.delta} > {output.deltafilter}

		# To view a summary of all the alignments produced by NUCmer
		echo "Running show-coords"
		# For Ribbon http://genomeribbon.com/
		echo "...for Ribbon"
		show-coords -r -lTH {output.delta} > {output.coords}
		show-coords -r -lTH {output.deltafilter} > {output.coordsfilter}

		"""
		# --mum  Use anchor matches that are unique in both the reference and query
		# --mumreference  Use anchor matches that are unique in in the reference
        #           but not necessarily unique in the query (default behavior)
        # -c|mincluster   Sets the minimum length of a cluster of matches (default 65)
        # -b|breaklen     Set the distance an alignment extension will attempt to extend poor scoring regions before giving up (default 200)


rule makelinks:
	""" Prepare a links with the MUMmer alignments """
	input:
		"mummer/{reference}-vs-{sample}.filter.coords",
	output:
		"links/{reference}-vs-{sample}.txt",
	shell:
		"""
		## Prepare a links file assigning the colors based on the chromosome names (and remove the mitochondrial contigs)
		if [ {wildcards.reference} = "Podan2" ]; then # If reference is Podan2, assume it doesn't have the name of the sample in the scaffolds
			grep -v "PaMt_NC_001329.3\\|_mt" {input} | awk '{{print "Podan2_"$10,$1,$2,$11,$3,$4,"color="}}' | awk '{{a=$1; b=$4; c=gensub(/(.*)hromosome_([0-9])([\.0-9]*)/, "chr\\\\2", "g", a); d=gensub(/(.*)hromosome_([0-9])([\.0-9]*)/, "chr\\\\2", "g", b); if( c==d ) {{print $0 c "_a5"}} else {{print $0 c}}}}' > {output}
		elif [ {wildcards.sample} = "Podan2" ]; then
			grep -v "PaMt_NC_001329.3\\|_mt" {input} | awk '{{print $10,$1,$2,"Podan2_"$11,$3,$4,"color="}}' | awk '{{a=$1; b=$4; c=gensub(/(.*)hromosome_([0-9])([\.0-9]*)/, "chr\\\\2", "g", a); d=gensub(/(.*)hromosome_([0-9])([\.0-9]*)/, "chr\\\\2", "g", b); if( c==d ) {{print $0 c "_a5"}} else {{print $0 c}}}}' > {output}		
		else
			grep -v "PaMt_NC_001329.3\\|_mt" {input} | awk '{{print $10,$1,$2,$11,$3,$4,"color="}}' | awk '{{a=$1; b=$4; c=gensub(/(.*)hromosome_([0-9])([\.0-9]*)/, "chr\\\\2", "g", a); d=gensub(/(.*)hromosome_([0-9])([\.0-9]*)/, "chr\\\\2", "g", b); if( c==d ) {{print $0 c "_a5"}} else {{print $0 c}}}}' > {output}
		fi
		"""

# ------------ Distribution of TEs ------------

rule index4bed:
	""" Create and index-like bed file for the sample comparisons """
	input:
		"fastas/{reference}-vs-{sample}.fa.fai"
	output:
		temp("tracks/{reference}-vs-{sample}.bed")
	params:
		threads = 1,
	shell:
		"""
		# It has to be tabs or BEDtools won't like it
		# grep -v {wildcards.sample} {input} | awk {{'print "Podan2_"$1"\\t"1"\\t"$2'}} > {output} 
		# grep {wildcards.sample} {input} | awk {{'print $1"\\t"1"\\t"$2'}} >> {output} 

		cat {input} | awk {{'print $1"\\t"1"\\t"$2'}} >> {output} 
		"""

rule makewindows:
	""" Use the BEDtools makewindows to get windows of the host genome """
	input:
		"tracks/{reference}-vs-{sample}.bed"
	output:
		temp("tracks/{reference}-vs-{sample}_hostwins.bed")
	params:
		threads = 1,
	shell:
		"""
		bedtools makewindows -b {input} -w 50000 -s 10000 > {output}
		"""	

rule TEsbed: 
	""" Prepare data for TE bed from the reference and the compared sample """
	input:
		ref = path2repeatmasking + "/{reference}.repeatmasker.gff3",
		query = path2repeatmasking + "/{sample}.repeatmasker.gff3",
	output:
		temp("tracks/{reference}-vs-{sample}_TEs.bed")
	params:
		threads = 1,
	shell:
		"""
		if [ {wildcards.reference} = "Podan2" ]; then # If reference is Podan2, assume it doesn't have the name of the sample in the scaffolds
			grep 'RepeatMasker' {input.ref} | awk {{'print "Podan2_"$1"\\t"$4"\\t"$5'}} > {output}
			grep 'RepeatMasker' {input.query} | awk {{'print $1"\\t"$4"\\t"$5'}} >> {output}

		elif [ {wildcards.sample} = "Podan2" ]; then
			grep 'RepeatMasker' {input.query} | awk {{'print "Podan2_"$1"\\t"$4"\\t"$5'}} > {output}
			grep 'RepeatMasker' {input.ref} | awk {{'print $1"\\t"$4"\\t"$5'}} >> {output}
		else
			cat {input} | grep 'RepeatMasker' | awk {{'print $1"\\t"$4"\\t"$5'}} >> {output}
		fi

		# grep 'RepeatMasker' {input.ref} | awk {{'print "Podan2_"$1"\\t"$4"\\t"$5'}} > {output}
		# grep 'RepeatMasker' {input.query} | awk {{'print $1"\\t"$4"\\t"$5'}} >> {output}
		"""

rule BEDtools_TEs:
	""" Use BEDtools coverage to produce a distribution of TEs along the host genome """
	# https://bedtools.readthedocs.io/en/latest/content/tools/coverage.html
	input:
		chrs = "tracks/{reference}-vs-{sample}_hostwins.bed",
		tes = "tracks/{reference}-vs-{sample}_TEs.bed",
	output:
		"tracks/{reference}-vs-{sample}_TEdistribution.txt"
	wildcard_constraints:
		sample="[a-zA-Z0-9.]+"
	params:
		threads = 1,
	shell:
		"""
		bedtools coverage -a {input.chrs} -b {input.tes} | awk '{{ print $1,$2,$3,$7 }}' > {output}
		"""

# ------------ Circos ------------


rule makecircosconfig:
	""" Modify a base circos configuration file for each sample """
	input:
		links = "links/{reference}-vs-{sample}.txt",
		trackTEs = "tracks/{reference}-vs-{sample}_TEdistribution.txt",

		# Other files provided by the user
		config = circosconfpath + "/circos.conf",
		colors = circosconfpath + "/colors.conf",
		colors_fonts = circosconfpath + "/colors_fonts_patterns.conf",
		fonts = circosconfpath + "/fonts.conf",
		housekeeping = circosconfpath + "/housekeeping.conf",
		ideogram = circosconfpath + "/ideogram.conf",
		image = circosconfpath + "/image.conf",
		patterns = circosconfpath + "/patterns.conf",
		ticks = circosconfpath + "/ticks.conf",
	output:
		"circos_{reference}-vs-{sample}/etc/circos.conf"
	shell:
		"cat {input.config} | sed 's;links/mummer.txt;../{input.links};' | sed 's;tracks/TEdistribution.txt;../{input.trackTEs};' > {output}"

	# run:
		# shell(f"cp -r {circosconfpath} circos_{wildcards.reference}-vs-{wildcards.sample}/")
	

rule circos:
	""" Run Circos to plot the alignment """
	input:
		"circos_{reference}-vs-{sample}/karyotype.txt",
		"circos_{reference}-vs-{sample}/etc/circos.conf",
		"links/{reference}-vs-{sample}.txt",
		circosconfpath + "/ideogram.conf", # Just so it checks if something changed
	output:
		"circos_{reference}-vs-{sample}/circos_{reference}-vs-{sample}.png"
	shell:
		"cd circos_{wildcards.reference}-vs-{wildcards.sample}; circos ;"
		"mv circos.png circos_{wildcards.reference}-vs-{wildcards.sample}.png"

