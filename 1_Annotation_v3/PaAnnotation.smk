# -*- snakemake -*-

### A pipeline to annotate genomes from the Podospora anserina species complex 
#############################################################################

# Third version of the annotation, prompted by the improved annotation of the S+ genome:

# Lelandais et al. (2022) New insights into genome annotation in Podospora
# anserina through re-exploiting multiple RNA-seq data. BMC Genomics 23, 859.
# https://doi.org/10.1186/s12864-022-09085-4

# Version 3:
# Here, I completely removed the call for abinitio prediction outside of
# MAKER, and deleted the MAKER lines using Augustus because they were not
# even used in the previous version anyway.
# The runTransDecoder.sh script was replaced by a rule

# https://bioinformaticsworkbook.org/dataAnalysis/GenomeAnnotation/Intro_To_Maker.html

# TODO:
# - Transform AbInitio.sh into Snakemake rules to improve transparency and adaptability

# OLD TODO:
# - change the evm.ctl to give more weight to genemark
# - Check that MAKERcpus works (not tested yet)

#############################################################################
# ===========================================================================
# Sandra Lorena Ament Velasquez
# 2023/04/13
# ++++++++++++++++++++++++++++++++++++++++++++++
# Version 3

# -------------------------------------------------
from Bio import SeqIO
# from Bio.Alphabet import generic_dna # depricated
import os.path

# -------------------------------------------------
# DATA from the configuration file
configfile: "config/PaAnnotation_config.yaml"
# -------------------------------------------------

# ------- Set configuration variables --------
annotationversion = "{0:.2f}".format(config["annotationversion"])

## Samples fixed paths
nicegenomespath = config["nicegenomespath"]
illupaths = config["illupaths"]
Podan2 = config["Podan2"]
PODCO = config["PODCO"]

# Samples IDs
samplesnp = config["samplesnp"]
samplespb = config["samplespb"]
samplesillu = config["samplesillu"]
samplessingleend = config["samplessingleend"]

## RNAseq data
RNAseqpath = config["RNAseqpath"]
samplesRNA = config["samplesRNA"]
parentals = config["parentals"]

## TEs
TElib = config["PodoTElib"]

## Scripts
gff2gff3 = config["gff2gff3"]
GffgenesIDFix = config["GffgenesIDFix"]
gff3addproduct = config["gff3addproduct"]
ortho2tag = config["ortho2tag"]
addManualGFF = config["addManualGFF"]

## Samples called
manualcurations = config["manualcurations"]
AllSamples = samplesnp + samplespb + samplesillu + ["Podan2", "PODCO"]

## Training files
snapHMM = config["snapHMM"]
GeneMarkMod = config["GeneMarkMod"]

## Reference genome
refpodan = config["refpodan"]
refgff = config["refgff"]
refgff_scv = config["refgff_scv"]

## Other evidence
comataprot = config["comataprot"]
curatedprots = config["curatedprots"]

## Known gene products
geneproducts = config["geneproducts"]

## Exonerate chunks
exochunks = list(range(1, 21))
nuchunks = len(exochunks)

## Filtering
MINSIZE = config["MINSIZE"] # Minimum size of a contig to get annotated

# -------------------------------------------------
# Hard-coded parameters
# -------------------------------------------------
MAKERcpus = 2 # Jobs are not very efficient, so it's not worth it to give that many threads.

# -------------------------------------------------
# Helper functions
# -------------------------------------------------
distancecodesdic = {"CBS112042p": "Pb", 
				"CBS237.71m": "Pp", 
				"CBS411.78m": "Ppp", 
				"CBS415.72m": "Ppc", 
				"PaTgp": "Pa", 
				"PaWa137m": "Pa", 
				"PaYp": "Pa", 
				"PcWa139m": "Pc", 
				"CBS124.78p": "Ppa", 
				"PaWa100p": "Pa", 
				"PaWa21m": "Pa", 
				"PaWa28m": "Pa", 
				"PaWa46p": "Pa", 
				"PaWa53m": "Pa", 
				"PaWa58m": "Pa", 
				"PaWa63p": "Pa", 
				"PaWa87p": "Pa", 
				"PaYm": "Pa", 
				"PaZp": "Pa",
				"CBS433.50p": "Pa",
				"CBS455.64m": "Pa",
				"CBS253.71p": "Ppa", 
				"CBS307.81m": "Cs", 
				"CBS333.63p": "Pp", 
				"CBS451.62p": "Pp", 
				"PcWa131m": "Pc", 
				"PcWa132p": "Pc", 
				"PcWa133m": "Pc", 
				"Podan2": "Pa", 
				"PODCO": "Pc",
				"Chr4": "Pa"}

def getperind(wildcards):
	samplecode = distancecodesdic[wildcards.sample]
	if samplecode == "Pa": 
		thissampleperiden = 98
	elif samplecode in ["Pc", "Pp"]: #(samplecode == "Pc") or ():
		thissampleperiden = 93
	elif samplecode in ["Pb", "Ppa", "Ppp", "Ppc"]:
		thissampleperiden = 90
	elif samplecode == "Cs":
		thissampleperiden = 77
	return(thissampleperiden)

def distancecodes(wildcards):
	return(distancecodesdic[wildcards.sample])

# The same thing but different names
def typesample(wildcards):
	""" Function to make a distinction between short and long-read assemblies """
	if wildcards.sample in samplesillu:
		genome = f"temp/genomes/{wildcards.sample}_tempnames.fa"
	else:
		genome = f"data/genomes/{wildcards.sample}.fa"
	return(genome)

def typesample_rm(wildcards):
	""" Function to make a distinction between short and long-read assemblies """
	if wildcards.sample in samplesillu:
		genome = f"RepeatMasker/{wildcards.sample}/{wildcards.sample}.repeatmasker.illu.gff"
	else:
		genome = f"RepeatMasker/{wildcards.sample}/{wildcards.sample}.fa.out.gff"
	return(genome)

def typesample_maker(wildcards):
	""" Function to make a distinction between short and long-read assemblies """
	if wildcards.sample in samplesillu:
		genome = f"MAKER/{wildcards.sample}/MAKER_output/{wildcards.sample}.all.illu.gff"
	else:
		genome = f"MAKER/{wildcards.sample}/MAKER_output/{wildcards.sample}.all.gff"
	return(genome)

# For the cluster
# https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html#snakefiles-resources
# def get_mem_mb(wildcards, threads):
#     return threads * 6800 # each core gets at most 6.8 GB of RAM in Rackham

# -------------------------------------------------

# ----------
# Rules not submitted to a job
localrules: assemblies, rnaseqraw, gtf2gff, gff2gff3, referenceprots, referencetrans, referencegens, makerconfig, otherprots, collectmakerresults, collectRMresults, getfastamaker, get_gtfRM2gff, sppcodenicegff, locus_tag_names, newnames, renameMAKERresults, get_gffutils2fasta, get_gffnummerator, sortgff, newcodesgff, migratenames, gtfRM2gff, add_manual_curations, sort_finalgff
# ----------

rule all:
	input:
		#### MAKER
		expand("results/{sample}/{sample}.all.gff3", sample = AllSamples), # MAKER
		expand("results/{sample}/{sample}.prenice-{version}.gff3", sample = AllSamples, version = annotationversion), # Renamed, final models
		expand("results/{sample}/{sample}.prenice-{version}.gff3", sample = manualcurations, version = float(annotationversion) + 0.01), # Renamed, final models

		### RNAseq (produced for my convenience, but only a handful are used for MAKER)
		# Stats
		expand("STAR/{genome}/{rnasample}-to-{genome}_flagstat.txt", zip, rnasample = samplesRNA, genome = parentals),
		# Index for IGV
		expand("STAR/{genome}/{rnasample}-to-{genome}_Aligned.sortedByCoord.out.bam.bai", zip, rnasample = samplesRNA, genome = parentals),
		# Transcripts
		expand("Cufflinks/{genome}/{rnasample}/{rnasample}_transcripts.gff3", zip, rnasample = samplesRNA, genome = parentals),
		# TransDecoder
		expand("TransDecoder/{genome}/{rnasample}/{rnasample}.genome.transdecoder.gff3", zip, rnasample = samplesRNA, genome = parentals),

		### RepeatMasker
		expand("results/{sample}/{sample}.repeatmasker.gff3", sample = AllSamples ), # RepeatMasker

# ---------------------------
### Get additional scripts
# ---------------------------

rule get_gffutils2fasta:
	""" Get the scripts needed for the rest of the pipeline """
	output:
		"scripts/gffutils2fasta.py",
	shell:
		"wget https://raw.githubusercontent.com/SLAment/Genomics/master/GenomeAnnotation/gffutils2fasta.py -P scripts"

rule get_gffnummerator:
	""" Get the scripts needed for the rest of the pipeline """
	output:
		"scripts/GFFnumerator.py",
	shell:
		"wget https://raw.githubusercontent.com/SLAment/Genomics/master/GenomeAnnotation/GFFnumerator.py -P scripts"


rule get_gtfRM2gff:
	""" Get the scripts needed for the rest of the pipeline """
	output:
		"scripts/gtfRM2gff.py",
	shell:
		"wget https://raw.githubusercontent.com/SLAment/Genomics/master/GenomeAnnotation/gtfRM2gff.py -P scripts"


# ------- PREPARE ALL DATA --------

rule assemblies:
	""" Prepare a folder with the assemblies """
	output:
		genome = "data/genomes/{sample}.fa",
	params:
		nicegenomespath = nicegenomespath,
		illupaths = illupaths,
		Podan2 = Podan2,
		PODCO = PODCO,
	wildcard_constraints:
		sample="[a-zA-Z0-9.]+" # to break ambiguity with the tempnames files
	run:
		if (wildcards.sample in samplesnp) or (wildcards.sample in samplespb):
			shell("ln -sf {params.nicegenomespath}/{wildcards.sample}.nice.fa {output.genome}")
		elif wildcards.sample in samplesillu: # Remove very small contigs
			fullassembly = params.illupaths + "/" + wildcards.sample + ".spades_scf.sn.fa" # The name of the SPAdes sense-corrected scaffolds
			filteredseqs = [] # Setup an empty list
			thisminsize = MINSIZE
			for seq_record in SeqIO.parse(fullassembly, "fasta"):
				# To deal with C. samala specifically
				if wildcards.sample in samplessingleend:
					thisminsize = 8000 # To get the small contigs

				if len(seq_record) >= thisminsize:
					filteredseqs.append(seq_record)

			# Write to a new file
			SeqIO.write(filteredseqs, output.genome, "fasta")

			# Or just make a symlink of the whole thing
			## shell("ln -sf {params.illupaths}/{wildcards.sample}.spades_scf.sn.fa {output.genome}")
		elif wildcards.sample == "Podan2":
			shell("ln -s {params.Podan2} {output.genome}")
		elif wildcards.sample == "PODCO":
			shell("ln -s {params.PODCO} {output.genome}")

rule rnaseqraw:
	""" Prepare a folder with the RNAseq data """
	output:
		read1 = "data/rnaseq/{sample}_postQC.1.fq.gz",
		read2 = "data/rnaseq/{sample}_postQC.2.fq.gz",
	params:
		RNAseqpath = RNAseqpath,
		samplesRNA = samplesRNA,
	run:
		shell("ln -sf {params.RNAseqpath}/{wildcards.sample}_postQC.1.fq.gz {output.read1}")
		shell("ln -sf {params.RNAseqpath}/{wildcards.sample}_postQC.2.fq.gz {output.read2}")


rule referenceprots:
	""" Get the protein sequences of the reference genome """
	input:
		genome = refpodan,
		gff = refgff,
		gff3tofasta = "scripts/gffutils2fasta.py"
	output:
		fasta = "data/reffastas/Podan_aa.fas",
	conda:
		"envs/general.yaml"
	shell:
		"python {input.gff3tofasta} {input.genome} {input.gff} --type CDS --join --proteinon --onlyids --output {output.fasta}; "
		"sed -i 's/gene:\\([a-zA-Z0-9_]*\\).G/\\1/' {output.fasta}" # The Podans_v2016 IDs are named gene:Pa_x_xxx.G

rule referencetrans:
	""" Get the transcript sequences of the reference genome """
	input:
		genome = refpodan,
		gff = refgff_scv,
		gff3tofasta = "scripts/gffutils2fasta.py"
	output:
		fasta = "data/reffastas/Podan_mRNA.fas",
	conda:
		"envs/general.yaml"
	shell:
		"python {input.gff3tofasta} {input.genome} {input.gff} --type exon --join --onlyids --output {output.fasta}; "
		"sed -i 's/gene:\\([a-zA-Z0-9_]*\\).G/\\1/' {output.fasta}" # The Podans_v2016 IDs are named gene:Pa_x_xxx.G

rule referencegens:
	""" Get the transcript sequences of the reference genome """
	input:
		genome = refpodan,
		gff = refgff,
		gff3tofasta = "scripts/gffutils2fasta.py"
	output:
		plus = temp("data/reffastas/Podan_genes_plus.fas"),
		plusandminus = "data/reffastas/Podan_genes.fas",
	conda:
		"envs/general.yaml"
	shell:
		"python {input.gff3tofasta} {input.genome} {input.gff} --type gene --onlyids --output {output.plus}; "
		"sed -i 's/gene:\\([a-zA-Z0-9_]*\\).G/\\1/' {output.plus}; " # The Podans_v2016 IDs are named gene:Pa_x_xxx.G

rule otherprots:
	""" Get the protein sequences of other species """
	output:
		comata = "data/OtherSpp/Pcomata_aa.fas",
		curated = "data/OtherSpp/curated_aa.fas",
	params:
		comataprot = comataprot,
		curatedprots = curatedprots,
	shell:
		"""
		cp {params.comataprot} {output.comata}
		sed -i -E 's:([A-Za-z0-9]+)(.1)(.*):\\1:' {output.comata}

		ln -s {params.curatedprots} {output.curated}
		"""

# ---------------------------------

# ---------------------------
### RNAseq
# ---------------------------
rule starindex:
	""" Make an index of the genome for STAR """
	input:
		genome = "data/genomes/{sample}.fa",
	output:
		"STAR/GenomeStarIndex/{sample}_GenomeDir/SAindex"
	threads: 3
	resources:
		threads = 3,
		mem_mb = int(3 * 6800), # each core gets at most 6.8 GB of RAM in Rackham
		time = "30:00",
	params:
		indexdir = "STAR/GenomeStarIndex/{sample}_GenomeDir",
		ram = int(3 * 6.8 * 1000000000) # A Rackham node contains 128 GB of RAM and 20 compute cores (each core gets at most 6.8 GB).
	conda:
		"envs/general.yaml"
	shell:
		"""
		mkdir -p temp/STAR/

		STAR --runMode genomeGenerate --genomeDir {params.indexdir} --genomeFastaFiles {input.genome} --runThreadN {resources.threads} --limitGenomeGenerateRAM {params.ram} --genomeLoad NoSharedMemory --outTmpDir temp/STAR/{wildcards.sample} # --genomeSAindexNbases 3
		
		# genomeLoad=NoSharedMemory, shared memory is not used. This option is recommended if the shared memory is not configured properly on your server.
		# --genomeSAIndexNbases 4 or 5 for small genomes (formula min(14, log2(GenomeLength)/2 - 1)) # log2(37000000)/2 -1 = 11.57051
		"""

rule star:
	""" Map the RNAseq reads to a genome using STAR """
	input:
		genome = "data/genomes/{genome}.fa",
		index = "STAR/GenomeStarIndex/{genome}_GenomeDir/SAindex",
		read1 = "data/rnaseq/{rnasample}_postQC.1.fq.gz",
		read2 = "data/rnaseq/{rnasample}_postQC.2.fq.gz",
	output:
		output = "STAR/{genome}/{rnasample}-to-{genome}_Aligned.sortedByCoord.out.bam"
	resources:
		threads = 3,
		mem_mb = int(3 * 6800), # each core gets at most 6.8 GB of RAM in Rackham
		time = "10:00:00",
	params:
		indexdir = "STAR/GenomeStarIndex/{genome}_GenomeDir",
	conda:
		"envs/general.yaml"
	shell:
		"""
		STAR --genomeDir {params.indexdir} --readFilesIn {input.read1} {input.read2} \
		--runThreadN {resources.threads} \
		--alignIntronMax 1000 \
		--readFilesCommand zcat \
		--outFileNamePrefix STAR/{wildcards.genome}/{wildcards.rnasample}-to-{wildcards.genome}_ \
		--outSAMstrandField intronMotif \
		--outSAMtype BAM SortedByCoordinate

		# --outSAMattributes NH HI AS nM XS 

		# outSAMattributes To make it compatible with SAMtools
		# outSAMstrandField intronMotif to make it compatible with Cufflinks/Cuffdiff

		# For unstranded RNA-seq data, Cufflinks/Cuffdiff require spliced alignments
		# with XS strand attribute, which STAR will generate with --outSAMstrandField
		# intronMotif option. As required, the XS strand attribute will be generated for
		# all alignments that contain splice junctions. The spliced alignments that have
		# undefined strand (i.e. containing only non-canonical unannotated junctions)
		# will be suppressed.

		# If you have stranded RNA-seq data, you do not need to use any specific STAR
		# options. Instead, you need to run Cufflinks with the library option --library-
		# type options. For example, cufflinks ... --library-type fr-firststrand should
		# be used for the standard dUTP protocol, including Illumina’s stranded Tru-Seq.
		# This option has to be used only for Cufflinks runs and not for STAR runs.

		# --outSAMtype BAM SortedByCoordinate output sorted by coordinate Aligned.sortedByCoord.out.bam le, similar to samtools sort command.
		"""

rule bamstats:
	""" Get some statistics of the STAR BAM file """
	input:
		"STAR/{genome}/{rnasample}-to-{genome}_Aligned.sortedByCoord.out.bam"
	output:
		"STAR/{genome}/{rnasample}-to-{genome}_flagstat.txt"
	resources:
		threads = 1,
		mem_mb = 6800, # each core gets at most 6.8 GB of RAM in Rackham
		time = "30:00",
	conda:
		"envs/general.yaml"
	shell:
		"samtools flagstat {input} > {output}"

rule bamindex:
	""" Make an index for the BAM file to read it in IGV """
	input:
		"STAR/{genome}/{rnasample}-to-{genome}_Aligned.sortedByCoord.out.bam"
	output:
		"STAR/{genome}/{rnasample}-to-{genome}_Aligned.sortedByCoord.out.bam.bai"
	resources:
		threads = 1,
		mem_mb = 6800, # each core gets at most 6.8 GB of RAM in Rackham
		time = "30:00",
	conda:
		"envs/general.yaml"
	shell:
		"samtools index {input}"

rule cufflinks:
	""" Produce transcript models with Cufflinks """
	input:
		"STAR/{genome}/{rnasample}-to-{genome}_Aligned.sortedByCoord.out.bam",
	output:
		"Cufflinks/{genome}/{rnasample}/transcripts.gtf"
	resources:
		threads = 10,
		mem_mb = 10 * 6800, # each core gets at most 6.8 GB of RAM in Rackham
		time = "05:00:00",
	conda:
		"envs/cufflinks.yaml"
	shell:
		"""
		# Stranded RNAseq 
		cufflinks {input} --library-type fr-firststrand -L {wildcards.rnasample}_CUFF --output-dir Cufflinks/{wildcards.genome}/{wildcards.rnasample} --num-threads {resources.threads} 2>&1  # Redirect the stderr to stdout
		# fr-firststrand --> dUTP libraries

		# mv Cufflinks/{wildcards.genome}/transcripts.gtf Cufflinks/{wildcards.genome}/{wildcards.rnasample}-to-{wildcards.genome}_transcripts.gtf 
		"""

rule gtf2gff:
	""" Transform the gtf Cufflinks file into gff3 """
	input:
		"Cufflinks/{genome}/{rnasample}/transcripts.gtf"
	output:
		temp("Cufflinks/{genome}/{rnasample}/{rnasample}_transcripts.gff")
	conda:
		"envs/cufflinks.yaml"
	shell:
		"""
		gffread -E {input} -o- > {output} # add | tail -n +3 to take out the headers
		"""

rule gff2gff3:
	""" Transform the gtf Cufflinks file into gff3 """
	input:
		"Cufflinks/{genome}/{rnasample}/{rnasample}_transcripts.gff"
	output:
		"Cufflinks/{genome}/{rnasample}/{rnasample}_transcripts.gff3"
	params:
		gff2gff3 = gff2gff3
	conda:
		"envs/cufflinks.yaml"
	shell:
		"""
		{params.gff2gff3} {input} --prefix {wildcards.rnasample}_CUFF > {output}
		"""

# ---------------------------
### TransDecoder
# ---------------------------

rule transdecoder:
	""" Run TransDecoder to get ORFs from the transcripts """
	# https://github.com/TransDecoder/TransDecoder/wiki
	input:
		transcripts = "Cufflinks/{genome}/{rnasample}/transcripts.gtf",
		genome = "data/genomes/{genome}.fa",
	output:
		fasta = "TransDecoder/{genome}/{rnasample}/{rnasample}.fa", # Transcripts
		gff3temp = "TransDecoder/{genome}/{rnasample}/{rnasample}.gff3",
		gff3 = "TransDecoder/{genome}/{rnasample}/{rnasample}.genome.transdecoder.gff3",
	resources:
		threads = 3,
		mem_mb = (3 * 6800), # each core gets at most 6.8 GB of RAM in Rackham
		time = "01:00:00",
	# conda: 
	# 	"envs.yaml"
	# I couldn't use a conda environment because I was getting an error: `Can't locate URI/Escape.pm in @INC`
	# even tho it worked fine if I manually loaded the environment and run the commands. Somehow the 
	# environment might be getting ignored for some perl paths? It works if loaded externally as an Uppmax module.
	# module load bioinfo-tools TransDecoder/5.7.0
	shell:
		"""
		echo " ... Construct the transcript fasta file using the genome and the transcripts.gtf file ... "
		gtf_genome_to_cdna_fasta.pl {input.transcripts} {input.genome} > {output.fasta}

		echo " ... Convert the transcript structure GTF file to an alignment-GFF3 formatted file ... "
		# perl -e 'use URI::Escape';
		gtf_to_alignment_gff3.pl {input.transcripts} > {output.gff3temp}

		# the final outputs are reported in your current working directory, so that's annoying
		cd TransDecoder/{wildcards.genome}/{wildcards.rnasample}

		echo " ... Generate the best candidate ORF predictions ... "
		TransDecoder.LongOrfs -S -t ../../../{output.fasta} # --output_dir TransDecoder/{wildcards.genome}/{wildcards.rnasample}

		# By default, TransDecoder.LongOrfs will identify ORFs that are at least 100
		# amino acids long. You can lower this via the -m parameter, but know that the
		# rate of false positive ORF predictions increases drastically with shorter
		# minimum length criteria.

		# If the transcripts are oriented according to the sense strand, then include
		# the -S flag to examine only the top strand.

		echo " ... Predict the likely coding regions ... " # This takes the longest time!
		TransDecoder.Predict -t ../../../{output.fasta} --cpu {resources.threads} # --output_dir TransDecoder/{wildcards.genome}/{wildcards.rnasample}
		
		cd ../../..
		
		echo " ... Generate a genome-based coding region annotation file ... "
		cdna_alignment_orf_to_genome_orf.pl {output.fasta}.transdecoder.gff3 {output.gff3temp} {output.fasta} > {output.gff3}
		"""

# ---------------------------
### Repeatmasking
# ---------------------------

# Illumina assemblies from SPAdes have names that are too long.
rule newnames:
	""" Make a new fasta temporary names for the scaffolds for Illumina assemblies """
	input:
		genome = "data/genomes/{sample}.fa"
	output:
		tempnames = "temp/genomes/{sample}_tempnames.txt", # To keep the original names somewhere
		newgenome = "temp/genomes/{sample}_tempnames.fa"
	# wildcard_constraints: # to avoid ambiguity with tempnames
	# 	sample="\d+"
	run:
		count = 1
		dictnames = {} # Setup a dictionary
		output_handle = open(output.newgenome, "w")

		for seq_record in SeqIO.parse(input.genome, "fasta"):
			if seq_record.id in dictnames.keys(): # Just in case
				print(f"ERROR: The sequence {seq_record.id} is more than once in the fasta file {input.genome}")
				raise ValueError(f"ERROR: The sequence {seq_record.id} is more than once in the fasta file {input.genome}")
			else:
				newname = "seq{0:06d}".format(count)
				dictnames[seq_record.id] = newname # record the original name somewhere
				count += 1 # Increase the count for the next sequence
				seq_record.id = newname
				seq_record.description = '' #Annoying extra info
				SeqIO.write(seq_record, output_handle, "fasta")

		# Keep the names somewhere
		with open(output.tempnames, 'w') as result:
			for key, value in dictnames.items():
				result.write(f"{key}\t{value}\n")

rule repeatmasker:
	""" Run RepeatMasker on the input genomes """
	input:
		genome = typesample
	output:
		"RepeatMasker/{sample}/{sample}.fa.out.gff"
	resources:
		threads = 10,
		mem_mb = int(10 * 6800), # each core gets at most 6.8 GB of RAM in Rackham
		time =  "5:00:00",
	params:
		TElib = TElib, # Custom library
	# conda: # The repeatmasker in the environment conflicts with the one loaded for MAKER
	# 	"envs/general.yaml"
	shell: # Conda environments are not allowed for the run directive
		"""
		RepeatMasker -pa {resources.threads} -a -xsmall -gccalc -gff -excln -lib {params.TElib} -dir RepeatMasker/{wildcards.sample} {input.genome}

		if [ -f RepeatMasker/{wildcards.sample}/{wildcards.sample}_tempnames.fa.out.gff ]; then 
			mv RepeatMasker/{wildcards.sample}/{wildcards.sample}_tempnames.fa.out.gff {output}
		fi
		"""
	# 	"RepeatMasker -pa {resources.threads} -a -xsmall -gccalc -gff -excln -lib {params.TElib} -dir RepeatMasker/{wildcards.sample} {input.genome}; "
	# 	"mv RepeatMasker/{wildcards.sample}/{wildcards.sample}*.fa.out.gff {output}" # Rename to deal with tempname
	# run: 
	# 	shell("RepeatMasker -pa {resources.threads} -a -xsmall -gccalc -gff -excln -lib {params.TElib} -dir RepeatMasker/{wildcards.sample} {input.genome}")
	# 	if os.path.exists("RepeatMasker/{wildcards.sample}/{wildcards.sample}_tempnames.fa.out.gff"):
	# 		shell("mv RepeatMasker/{wildcards.sample}/{wildcards.sample}_tempnames.fa.out.gff {output}") # Rename to deal with tempname


# ---------------------------
### MAKER
# ---------------------------

rule makerconfig:
	""" Produce MAKER configuration files for the sample """
	# This depends on Uppmax module
	# module load bioinfo-tools maker/3.01.04
	# See http://weatherby.genetics.utah.edu/MAKER/wiki/index.php/The_MAKER_control_files_explained
	input:
		genome = typesample, # If it's illumina, you need the simple names

		# EST evidence
		transcWa63 = "TransDecoder/PaWa63p/PaWa63m_RNA/PaWa63m_RNA.fa",
		transcWa58BC = "TransDecoder/PaWa58m/PaWa58BCvsS/PaWa58BCvsS.fa",
		transcWa131 = "TransDecoder/PcWa139m/PcWa131m_RNA/PcWa131m_RNA.fa",
		Podans_v2016_trans = "data/reffastas/Podan_mRNA.fas",
		# transcTd = "TransDecoder/PcWa139m/PcTdp_RNA/PcTdp_RNA.fa",

		# Other evidence
		podan2aa = "data/reffastas/Podan_aa.fas",
		comata = "data/OtherSpp/Pcomata_aa.fas",
		curated = "data/OtherSpp/curated_aa.fas",
		# repeatmaskergff = "RepeatMasker/{sample}/{sample}.fa.out.gff", # I wanted this to be used, but I suspect the gff won't be good
	output:
		"MAKER/{sample}/{sample}_bopts.ctl",
		"MAKER/{sample}/{sample}_exe.ctl",
		"MAKER/{sample}/{sample}_opts.ctl",
	params:
		TElib = TElib,
		snapHMM = snapHMM,
		GeneMarkMod = GeneMarkMod,
		repeatmaskergff = "", # I wanted this to be used, but I suspect the gff won't be good
		cpus = MAKERcpus,
	shell:
		"""	
		cd MAKER/{wildcards.sample}
	
		# Produce an empty configuration file
		maker -CTL

		# Rename them to match the sample
		mv maker_bopts.ctl {wildcards.sample}_bopts.ctl
		mv maker_exe.ctl {wildcards.sample}_exe.ctl
		mv maker_opts.ctl {wildcards.sample}_opts.ctl
		mv maker_evm.ctl {wildcards.sample}_evm.ctl

		### Fill the config files with the sample info
		# Set the current genome
		sed -i "s|^genome=|genome=../../{input.genome}|g" {wildcards.sample}_opts.ctl # Add that to the opts file
		
		# #-----Set evidence
		sed -i "s|^est=|est=../../{input.transcWa63},../../{input.transcWa58BC},../../{input.transcWa131},../../{input.Podans_v2016_trans}|g" {wildcards.sample}_opts.ctl 	# The RNAseq STAR+Cufflinks transcripts + Podans_v2016 transcripts
		sed -i "s|^protein=|protein=../../{input.podan2aa},../../{input.comata},../../{input.curated}|g" {wildcards.sample}_opts.ctl # Add that to the opts file

		# #-----Repeat Masking
		sed -i "s|^rmlib=|rmlib={params.TElib}|g" {wildcards.sample}_opts.ctl  #pre-identified repeat elements in a fasta file	
		# sed -i "s|^rm_gff=|rm_gff=../../{params.repeatmaskergff}|g" {wildcards.sample}_opts.ctl  #pre-identified repeat elements from an external GFF3 file

		# Notice I'm leaving it with softmask=1

		# #-----Gene Prediction
		# SNAP HMM file
		sed -i "s|^snaphmm=|snaphmm={params.snapHMM}|g" {wildcards.sample}_opts.ctl  
		# GeneMark HMM file
		sed -i "s|^gmhmm=|gmhmm={params.GeneMarkMod}|g" {wildcards.sample}_opts.ctl  

		# #infer predictions from protein homology, 1 = yes, 0 = no # Otherwise many small genes from Podan2 are lost, even tho some of them might not be real in fact (but some are for sure base on RNAseq)
		sed -i "s|^protein2genome=0|protein2genome=1|g" {wildcards.sample}_opts.ctl  
		# find tRNAs with tRNAscan, 1 = yes, 0 = no
		sed -i "s|^trna=0|trna=1|g" {wildcards.sample}_opts.ctl  
		# also run ab-initio prediction programs on unmasked sequence, 1 = yes, 0 = no
		sed -i "s|^unmask=0|unmask=1|g" {wildcards.sample}_opts.ctl  
				
		# #-----External Application Behavior Options
		sed -i "s|^cpus=1|cpus={params.cpus}|g" {wildcards.sample}_opts.ctl  #max number of cpus to use in BLAST and RepeatMasker (not for MPI, leave 1 when using MPI)
		
		# #-----MAKER Behavior Options
		# skip genome contigs below this length (under 10kb are often useless)
		sed -i "s|^min_contig=1|min_contig=10000|g" {wildcards.sample}_opts.ctl 
		# extra steps to force start and stop codons, 1 = yes, 0 = no
		sed -i "s|^always_complete=0|always_complete=1|g" {wildcards.sample}_opts.ctl # Jesper set this to 1, but I'm actually not sure.
		# length for the splitting of hits (expected max intron size for evidence alignments)
		sed -i "s|^split_hit=10000|split_hit=1000|g" {wildcards.sample}_opts.ctl # Jesper set this to 1000
		# minimum intron length (used for alignment polishing)
		sed -i "s|^min_intron=20|min_intron=2|g" {wildcards.sample}_opts.ctl # Some Podan2 genes have exons as short as 2 bp, although they are rare
		# #consider single exon EST evidence when generating annotations, 1 = yes, 0 = no
		sed -i "s|^single_exon=0|single_exon=1|g" {wildcards.sample}_opts.ctl # It's dangerous but I hope to recover genes/things expressed that are not in Podan2
		# #min length required for single exon ESTs if 'single_exon is enabled'
		sed -i "s|^single_length=250|single_length=500|g" {wildcards.sample}_opts.ctl # to compensate a bit
		
		# #limits use of ESTs in annotation to avoid fusion genes (This results in the loss of three prime UTR on the five prime gene and 
		# # loss of five prime UTR on the three prime gene, but it is better than a merged gene.)
		sed -i "s|^correct_est_fusion=0|correct_est_fusion=1|g" {wildcards.sample}_opts.ctl
		
		## Danger
		# #removes theVoid directory with individual analysis files, 1 = yes, 0 = no
		sed -i "s|^clean_up=0|clean_up=1|g" {wildcards.sample}_opts.ctl
		
		# This option will help save disk space by deleting individual results
		# files (such as blast, exonerate, and gene predictor outputs) once
		# they are no longer needed. If you have the disk space it is usually
		# best to keep this set to 0. Having those files around will make
		# rerunning MAKER much faster if necessary.

		"""

rule maker:
	""" Run MAKER """
	# This depends on Uppmax module
	# module load bioinfo-tools maker/3.01.04
	# See http://weatherby.genetics.utah.edu/MAKER/wiki/index.php/The_MAKER_control_files_explained
	input:
		"MAKER/{sample}/{sample}_bopts.ctl",
		"MAKER/{sample}/{sample}_exe.ctl",
		"MAKER/{sample}/{sample}_opts.ctl",
	output:
		"MAKER/{sample}/{sample}.dummy",
	resources:
		threads = MAKERcpus,
		mem_mb = int(MAKERcpus * 6800), # each core gets at most 6.8 GB of RAM in Rackham
		time = lambda wildcards: "10-00:00:00" if wildcards.sample == "CBS307.81m" else "8-00:00:00",
	shell:
		"""	
		cd MAKER/{wildcards.sample}
		maker -fix_nucleotides {wildcards.sample}_opts.ctl {wildcards.sample}_bopts.ctl {wildcards.sample}_exe.ctl {wildcards.sample}_opts.ctl && touch {wildcards.sample}.dummy
		"""

rule makermerge:
	""" Recover output of MAKER """
	# This depends on Uppmax module
	# module load bioinfo-tools maker/3.01.04
	input:
		"MAKER/{sample}/{sample}.dummy",
	output:
		"MAKER/{sample}/MAKER_output/{sample}.all.gff",
	resources:
		threads = 1,
		mem_mb = int(1 * 6800), # each core gets at most 6.8 GB of RAM in Rackham
		time =  "15:00",
	run:
		sample = wildcards.sample
		if wildcards.sample in samplesillu: # To deal with illumina samples
			sample += "_tempnames"

		shell("""
			cd MAKER/{wildcards.sample}/MAKER_output
			fasta_merge -d ../{sample}.maker.output/{sample}_master_datastore_index.log
			gff3_merge -d ../{sample}.maker.output/{sample}_master_datastore_index.log
			""")

		if wildcards.sample in samplesillu: # Rename the output so it works with the rest of the rulegraph
			shell("mv MAKER/{wildcards.sample}/MAKER_output/{sample}.all.gff {output}")


# ---------------------------
### Renaming
# ---------------------------

rule renameMAKERresults:
	input:
		gff = "MAKER/{sample}/MAKER_output/{sample}.all.gff",
		tempnames = "temp/genomes/{sample}_tempnames.txt"
	output:
		realnames = "MAKER/{sample}/MAKER_output/{sample}.all.illu.gff",
	run:
		# Make a dictionary
		names = {}

		with open(input.tempnames, 'r') as temps:
			for line in temps:
				cleanline = line.rstrip("\n").split("\t")
				names[cleanline[1]] = cleanline[0]

		# Keep the names somewhere
		ofile = open(output.realnames, 'w')

		with open(input.gff, 'r') as gff:
			for line in gff:
				if "#" in line:
					ofile.write(line)
				else:
					cleanline = line.split("\t")
					newname = cleanline[0].rstrip("\n").strip('>')

					if newname in names.keys():
						newline = line.replace(newname, names[newname])
						ofile.write(newline)
					else:
						ofile.write(line)

rule collectmakerresults:
	input:
		maker = typesample_maker,
		# maker = "MAKER/{sample}/MAKER_output/{sample}.all.gff",
	output:
		maker = "results/{sample}/{sample}.all.gff3",
		models = "renaming/{sample}/{sample}.models.gff3",
	shell:
		"""
		# Maker output
		cat {input.maker} | grep -v -P ".\\tcontig\\t" > {output.maker} # I dislike having the full contig when visualizing

		grep -P ".\\tmaker\\t" {output.maker} > {output.models}
		"""

rule sortgff:
	""" Sort the MAKER gene models """
	input:
		namedgff = "renaming/{sample}/{sample}.models.gff3",
	output:
		sortednamedgff = "renaming/{sample}/{sample}.models.sort.gff3",
	conda:
		"envs/general.yaml"
	shell:
		"""
		igvtools sort {input.namedgff} {output.sortednamedgff} 
		"""

def assignLocusTag(wildcards):
	if wildcards.sample == "CBS237.71m":
		return("EYR66")
	elif wildcards.sample == "CBS124.78p":
		return("QC764")
	elif wildcards.sample == "CBS411.78m":
		return("QC763")
	elif wildcards.sample == "CBS415.72m":
		return("QC762")
	elif wildcards.sample == "CBS112042p":
		return("QC761")
	else:
		return(wildcards.sample) # something generic since I won't put this in NCBI

rule newcodesgff:
	""" Give new codes to the MAKER gene models to make them cleaner """
	input:
		sortednamedgff = "renaming/{sample}/{sample}.models.sort.gff3",
		GFFnumerator = "scripts/GFFnumerator.py" 
	output:
		# panomgff = "renaming/{sample}/{sample}.models.sort.PaNom.gff3"
		panomgff = temp("renaming/{sample}/{sample}.models.sort.PaNom.gff3")
	# resources:
	# 	threads = 1,
	# 	mem_mb = int(1 * 6800), # each core gets at most 6.8 GB of RAM in Rackham
	# 	time =  "10:00",
	params:
		locus_tag = assignLocusTag
	conda:
		"envs/general.yaml"
	shell:
		"python {input.GFFnumerator} {input.sortednamedgff} --sample {params.locus_tag} --namestoo --step 10 > {output.panomgff}"

rule getfastamaker:
	""" Get the CDS sequences of the MAKER models """
	input:
		gff = "renaming/{sample}/{sample}.models.sort.PaNom.gff3",
		genome = "data/genomes/{sample}.fa",
		gff3tofasta = "scripts/gffutils2fasta.py",
	log: "renaming/{sample}/{sample}.models.sort.PaNom_weirdgenes.log"
	output:
		# cds = "renaming/{sample}/{sample}.models.sort.PaNom.fas",
		cds = temp("renaming/{sample}/{sample}.models.sort.PaNom.fas"),
	conda:
		"envs/general.yaml"
	shell:
		"python {input.gff3tofasta} {input.genome} {input.gff} --type noutrs --onlyids --output {output.cds} > {log}" # gffutils2fasta.py > v. 1.61
		# "python {input.gff3tofasta} {input.genome} {input.gff} --type noutrs --onlyids --output renaming/{wildcards.sample}/{wildcards.sample}.models.sort.PaNom > {input.gff}_weirdgenes.log" # gffutils2fasta.py v. 1.22

rule migratenames:
	""" Migrate the Podan2 names to the new annotation """
	input:
		makergenes = "renaming/{sample}/{sample}.models.sort.PaNom.fas",
		models = "renaming/{sample}/{sample}.models.sort.PaNom.gff3", # Use the unsorted GFF3 because IGV freaks out with the sorted one, and my script can handle it
		podan_genes = "data/reffastas/Podan_genes.fas", # The podan genes in the right order (fasta file)
		GffgenesIDFix = GffgenesIDFix,
	output:
		namedgff = "renaming/{sample}/{sample}.models.sort.PaNom.ID.gff3",
		onetoones = "renaming/{sample}/onetoones.txt",
		fusedloci = "renaming/{sample}/fusedloci.txt",
		splitgenes = "renaming/{sample}/splitgenes.txt",
	log: "renaming/{sample}/{sample}.models.sort.PaNom.gff3_ID.log",
	threads: 2
	resources:
		threads = lambda wildcards, threads: threads,
		mem_mb = lambda wildcards, threads: 6800 * threads,
		# threads = 2,
		# mem_mb = int(2 * 6800), # each core gets at most 6.8 GB of RAM in Rackham
		time =  "40:00",
	params:
		periden = getperind,
		tempdir = "renaming/{sample}/" # This is to prevent colliding with other samples
	conda:
		"envs/general.yaml"
	shell:
		"python {input.GffgenesIDFix} {input.podan_genes} {input.makergenes} {input.models} --identity {params.periden} --superverbose --clean --threads {resources.threads} --output {output.namedgff} --temp {params.tempdir} > {log}"

# rule addgeneproducts:
# 	""" Add gene products to genes that were identified """
# 	input:
# 		gff = "renaming/{sample}/{sample}.models.sort.PaNom.ID.gff3",
# 		geneproducts = geneproducts
# 	output:
# 		gff = temp("renaming/{sample}/{sample}.models.sort.PaNom.ID.prods.gff3"),
# 	params:
# 		gff3addproduct = gff3addproduct
# 	conda:
# 		"envs/general.yaml"
# 	shell:
# 		"{params.gff3addproduct} {input.gff} {input.geneproducts} > {output.gff}"

rule locus_tag_names:
	""" Change the IDs of genes to match the TAG_NNNN style of NCBI but retaining ortholog information """
	input:
		gff = "renaming/{sample}/{sample}.models.sort.PaNom.ID.gff3",
		ortho2tag = ortho2tag
	output:
		"renaming/{sample}/{sample}.models.sort.PaNom.ID.tag.gff3"
	conda:
		"envs/general.yaml"
	shell:
		"python {input.ortho2tag} {input.gff} > {output}"
	# The locus_tag prefix should be 3-12 alphanumeric characters and the first character may not be a digit. (case-sensitive)


rule sppcodenicegff:
	""" Change the species codes of the gene names """
	input:
		panomgff = "renaming/{sample}/{sample}.models.sort.PaNom.ID.tag.gff3"
		# panomgff = "renaming/{sample}/{sample}.models.sort.PaNom.ID.prods.gff3"
	output:
		nicegff = "results/{sample}/{sample}.prenice-" + annotationversion + ".gff3"
	params:
		sppcode = distancecodes,
	shell:
		"""
		sed "s;Pa_;{params.sppcode}_;g" {input.panomgff} > {output.nicegff} 
		"""

# ---------------------------
### Add manual curations
# ---------------------------

rule add_manual_curations:
	""" I manually curated gene models of the HNWD family in some important samples """
	# Assuming that any gene that overlaps with the annotated features must be removed
	input:
		gff1 = "results/{sample}/{sample}.prenice-" + annotationversion + ".gff3",
		gff2 = "data/manualcuration/{sample}_manual.gff3",
		addManualGFF = addManualGFF
	output:
		temp("results/{sample}/{sample}.prenice-" + str(float(annotationversion) + 0.01) + "-unsorted.gff3")
	conda:
		"envs/general.yaml"
	shell:
		"python {input.addManualGFF} {input.gff1} {input.gff2} > {output}"

rule sort_finalgff:
	""" The output of addManualGFF.py is not sorted so fix it """
	input:
		"results/{sample}/{sample}.prenice-" + str(float(annotationversion) + 0.01) + "-unsorted.gff3",
	output:
		"results/{sample}/{sample}.prenice-" + str(float(annotationversion) + 0.01) + ".gff3",
	conda:
		"envs/general.yaml"
	shell:
		"""
		igvtools sort {input} {output} 
		"""

# ---------------------------
### Collect results
# ---------------------------

rule renameRMresults:
	input:
		repeatmasker = "RepeatMasker/{sample}/{sample}.fa.out.gff",
		tempnames = "temp/genomes/{sample}_tempnames.txt"
	output:
		realnames = "RepeatMasker/{sample}/{sample}.repeatmasker.illu.gff",
	run:
		# Make a dictionary
		names = {}

		with open(input.tempnames, 'r') as temps:
			for line in temps:
				cleanline = line.rstrip("\n").split("\t")
				names[cleanline[1]] = cleanline[0]
		
		# Keep the names somewhere
		ofile = open(output.realnames, 'w')

		with open(input.repeatmasker, 'r') as gff:
			for line in gff:
				if "#" in line:
					ofile.write(line)
				else:
					cleanline = line.rstrip("\n").split("\t")
					cleanline[0] = names[cleanline[0]] # Replace the temporary name with the original
					newline = '\t'.join(cleanline) + '\n'
					ofile.write(newline)

rule gtfRM2gff:
	input:
		gtfRM2gff = "scripts/gtfRM2gff.py",
		repeatmasker = typesample_rm,
	output:
		"results/{sample}/{sample}.repeatmasker.gff3"
	shell:
		"python {input.gtfRM2gff} {input.repeatmasker} > {output}"

rule collectRMresults:
	input:
		repeatmasker = typesample_rm,
	output:
		repeatmasker = "results/{sample}/{sample}.repeatmasker.gff",
	shell:
		"cp {input.repeatmasker} {output.repeatmasker}"
