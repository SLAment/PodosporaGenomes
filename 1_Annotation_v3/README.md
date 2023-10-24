# PaAnnotation: A pipeline to annotate *Podospora* genomes

This is a modified [Snakemake](https://snakemake.github.io) pipeline similar to what I used in [Vogan et al. 2019](https://genome.cshlp.org/content/31/5/789) to annotate *Podospora* genomes. Genome annotation is a painful process with many, many tools, so I ended up using a mixture of modules from the [UPPMAX](https://www.uppmax.uu.se/) cluster in Sweden, and conda environments. This is, of course, not ideal for reproducibility, but hopefully peolple at least get an idea of what to do *after* you managed to install everything.

The goal of this pipeline is to produce a basic annotation with MAKER, while adding a few manually curated gene models in the end. The produced gff file is later functionally annotated with Funannotate in another pipeline.

## The input files

The multitude input files for this pipeline are specified in a configuration (`config/PaAnnotation_config.yaml`), which **has to be modified to match your local paths**!! These are:

### Genome annotation of the reference strain S+ 

- **Podans_v2016_trans**: The transcript DNA sequences from the annotation produced by [Lelandais et al. (2022)](https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-022-09085-4). I made this file by modifying their gff file (`refgff_scv` below) to be more of a gff3 file (with `GTFfunctionremover.py`), followed by my script [`gffutils2fasta.py`](https://github.com/SLAment/Genomics/blob/master/GenomeAnnotation/gffutils2fasta.py). 
- **Podans_v2016_aa**: The protein sequences from the same cleaned gff3 of Lelandais et al. (2022) as above, also produced with `gffutils2fasta.py`.
- **refpodan**: The genome assembly of the strain S+ from Lelandais et al. (2022).
- **refgff**: The annotation file for all genes from Lelandais et al. (2022) for the strain S+.
- **refgff_scv**: The modified annotation file for the transcripts that were curated to have proper transcriptional starting and/or ending sites, from Lelandais et al. (2022).

### Assemblies paths

- **nicegenomespath**: long-read assemblies path
- **illupaths**: I happen to have a different path for the assemblies produced from short-read data
- **Podan2**: For historical reasons I also used the previous genome assembly of the strain S+, available in JGI as Podan2. This fasta file is used simply to annotate it and to compare the old annotation with the one produced from my pipeline. It's also the one I used to make the Circos plots. This assembly is available [here](https://github.com/johannessonlab/SpokBlockPaper/blob/master/GettingTElibrary/data/genomes/Podan2.fa).
- **PODCO**: Genome assembly of the reference *P. comata* genome assembly ("PODCO"), from [Silar et al. (2018)](https://link.springer.com/article/10.1007/s00438-018-1497-3).

### Samples

To keep track of all the different assemblies I have, I divided them by source, although this has no other use.

- **samplesnp**: ["CBS112042p", "CBS237.71m", "CBS411.78m", "CBS415.72m", "PaTgp", "PaWa137m", "PaYp", "PcWa139m", "CBS124.78p"] # Produced with Nanopore
- **samplespb**: ["PaWa100p", "PaWa21m", "PaWa28m", "PaWa46p", "PaWa53m", "PaWa58m", "PaWa63p", "PaWa87p", "PaYm",] # Produced with PacBio
- **samplesillu**: ["CBS253.71p", "CBS307.81m", "CBS333.63p", "CBS451.62p", "PcWa131m", "PcWa132p", "PcWa133m", "PaZp", "CBS433.50p", "CBS455.64m"] # Illumina

- **manualcurations**: Samples were I did manual curation, so they will have the final (3.01) annotation file. Based on this list, the pipeline will look for gff files with the corresponding annotations in `data/manualcuration`, which are provided in the repo.

- **samplessingleend**: the one sample were the reads can only be used as single-end because of problems with the library. This assembly is super fragmented so I had to filter it extra to remove contigs smaller than 8 kb.

### Transposable elements library

- **PodoTElib**: I used our library [PodoTE-1.00.lib](https://github.com/johannessonlab/SpokBlockPaper/blob/master/Annotation/data/PodoTE-1.00.lib).

### RNAseq data

Here I have the path to our RNAseq data after QC, from [Vogan et al. (2019)](https://elifesciences.org/articles/46454) and [Vogan et al. (2021)](https://genome.cshlp.org/content/31/5/789).

- **RNAseqpath**: path to the clean reads

The next two lists match each other. In the first list I have the RNAseq reads name, and in the second I have the name of the strain where they came from, or the closest.

- **samplesRNA**: ["PaWa53BCvsS", "PaWa58BCvsS", "PaWa28BCvsS", "PaYBCvsS", "PaYBCvsS", "PcWa131m_RNA", "PaWa53BCvsS", "PaWa58BCvsS", "PaWa28BCvsS", "PaYBCvsS", "PaWa63m_RNA", "PcTdp_RNA", "PaWa53BCvsS", "PaWa58BCvsS", "PaWa28BCvsS", "PaYBCvsS"]
- **parentals**: ["PaWa53m", "PaWa58m", "PaWa28m", "PaYm", "PaYp", "PcWa139m", "Podan2", "Podan2", "Podan2", "Podan2", "PaWa63p", "PcWa139m", "PaWa137m", "PaWa137m", "PaWa137m", "PaWa137m"]

For example, PaWa53BCvsS came from a backcross of between PaWa53m (its parental) and the S+ genome. I chose PaWa53m because it contains the large Enterprise element, which I want to keep transcripts for.

### Scripts

- **gffread2EVM.py**: available [here](https://github.com/johannessonlab/SpokBlockPaper/blob/285ecc1729f8f4f70e9fd3b41939743237c56818/Annotation/scripts/gffread2EVM.py)
- **GffgenesIDFix**: provided here in the folder `scripts`
- **gff3addproduct.py**: available [here](https://github.com/johannessonlab/SpokBlockPaper/blob/285ecc1729f8f4f70e9fd3b41939743237c56818/Annotation/scripts/gff3addproduct.py#L6)
- **ortho2tag.py**: provided here in the folder `scripts`
- **addManualGFF.py**: provided here in the folder `scripts`

### Training files

- **snapHMM**: The SNAP training file `PaWa28m_PacBioChrPilon.hmm` available [here](https://github.com/johannessonlab/SpokBlockPaper/tree/285ecc1729f8f4f70e9fd3b41939743237c56818/Annotation/hmms).
- **GeneMarkMod**: The GeneMark training file `gmhmm.mod` available [here](https://github.com/johannessonlab/SpokBlockPaper/tree/285ecc1729f8f4f70e9fd3b41939743237c56818/Annotation/hmms).

### Other evidence

- **comataprot**: the protein sequences of the PODCO genome assembly. The pipeline will process them to produce `data/OtherSpp/Pcomata_aa.fas`, provided in the repo.
- **curatedprots**: A set of manually curated proteins available [here](https://github.com/johannessonlab/SpokBlockPaper/blob/285ecc1729f8f4f70e9fd3b41939743237c56818/Annotation/data/OtherSpp/curated_aa.fas), but also provided in the repo as `data/OtherSpp/curated_aa.fas`.

### Extras

- **MINSIZE**: a variable controlling the minimum size of a contig to be annotated (50000 bp)
- **geneproducts**: a file with two columns; the first one has the Pa_X_XXX gene code, and the second has the name of the gene. This is provided in `data/KnownProducts.txt`

## Building the environment

First I loaded the UPPMAX modules.

    $ module load bioinfo-tools snakemake/7.25.0 maker/3.01.04 TransDecoder/5.7.0
    $ export LD_PRELOAD=$MPI_ROOT/lib/libmpi.so   # Probably not used, but to avoid problems

GeneMark requires a key to be put in your home directory.

    $ cp -vf /sw/bioinfo/GeneMark/keyfile/gm_key $HOME/.gm_key
    ‘/sw/bioinfo/GeneMark/keyfile/gm_key’ -> ‘/home/sandral/.gm_key’

This implies the following programs are made available in the environment:

    $ module list

    Currently Loaded Modules:
      1) python3/3.9.5      5) bioinfo-tools              9) perl/5.26.2          13) exonerate/2.4.0    17) blast_databases/latest
      2) snakemake/7.25.0   6) GeneMark/4.38-es          10) perl_modules/5.26.2  14) tRNAscan-SE/1.3.1  18) blast/2.13.0+
      3) gcc/9.3.0          7) BioPerl/1.7.2_Perl5.26.2  11) augustus/3.4.0       15) maker/3.01.04      19) hmmer/3.2.1
      4) openmpi/4.0.3      8) RepeatMasker/4.1.0        12) snap/2013-11-29      16) uppmax             20) TransDecoder/5.7.0


For most Snakemake rules I made little conda environment that separates them from the Uppmax modules. It looks like so:

```yaml
channels:
  - bioconda
  - defaults
  - conda-forge
dependencies:
  - repeatmasker=4.1.2.p1
  - star=2.7.10b
  - samtools=1.17
  - transdecoder=5.7.0
  - igvtools=2.14.1
  - gffutils=0.11.1
```

For some reason I get an issue when using Transdecoder with the local environment when submited as a job in Uppmax. I ended up using the module instead. Same with RepeatMasker, because it conflicts with the MAKER module dependency (see above). Misteriously, if I don't include RepeatMasker and Transcoder in the environment, the other packages are not compatible anymore!!! So I left it as is... 

For Cufflinks I also used an individual environment:

```yaml
channels:
  - bioconda
  - defaults
  - conda-forge
dependencies:
  - cufflinks=2.2.1
```

For this pipeline I also decided to use a "profile", which gives general instructions to submit jobs into the slurm cluster UPPMAX. 

The profile config file looks like so:

```yaml
snakefile: PaAnnotation.smk
cluster:
  mkdir -p logs/{rule} &&
  sbatch
    --account={resources.account}
    --partition={resources.partition}
    --cpus-per-task={resources.threads}
    --mem={resources.mem_mb}
    --job-name={rule}
    --error=logs/{rule}/{rule}-{wildcards}-%j.err
    --output=logs/{rule}/{rule}-{wildcards}-%j.out
    --time={resources.time}
    --mail-type=ALL
    --parsable
default-resources:
  - account="XXXXXXXXX"
  - partition="core"
  - time="1:00:00"
  - threads=1
  - mem_mb=6800

restart-times: 1
max-jobs-per-second: 10
max-status-checks-per-second: 1
jobs: 100
keep-going: True
rerun-incomplete: True
printshellcmds: True
scheduler: greedy
use-conda: True
cluster-cancel: scancel
cluster-cancel-nargs: 50
```
Where XXXXXXXXX is your UPPMAX account.

## Run pipeline in UPPMAX

Get into the folder, for example:

    $ cd /proj/sllstore2017101/b2015200/SnakePipelines/5d_Annotation_v3

Activate the modules:

    $ module load bioinfo-tools snakemake/7.25.0 maker/3.01.04 TransDecoder/5.7.0

First, to get an idea of how the pipeline looks like we can make a rulegraph:

    $ snakemake --profile profile --rulegraph | dot -Tpng > rulegraph.png

To test that everything seems in order:

    $ snakemake --profile profile -pn

Run the pipeline:

    $ screen -R PaAnnotation
    $ module load bioinfo-tools snakemake/7.25.0 maker/3.01.04 TransDecoder/5.7.0
    $ snakemake --profile profile &> PaAnnotation.log &
    [1] 2249


## Results

The ultimate goal is to produce a gff3 called "prenice" that contains the MAKER annotation. For some strains I also did manual curation of a few gene models, mainly NOD-like receptors (NLRs) from the HNWD family [(Paoletti et al. 2007)](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0000283).

For each strain, in the end there will be a folder with the strain name within the strain `results` containing:

- `{sample}.all.gff3` -- The raw result from MAKER, including repeats
- `{sample}.prenice-3.00.gff3` -- A processed MAKER gff of the gene models where I attempted to assign one-to-one orthologs names corresponding to the Podan2 gene codes. For example, the gene Pa_5_10 should now be called NNNN_500010, where NNNN is a locus_tag assigned by GenBank to the genomes, or just the strain name.
- `{sample}.repeatmasker.gff3"` -- A modified gff3 file of the genome annotation using the repeat library [PodoTE-1.00.lib](https://github.com/johannessonlab/SpokBlockPaper/blob/master/Annotation/data/PodoTE-1.00.lib).

For those species with manual annotation, there will also be the ultimate file:

- `{sample}.prenice-3.01.gff3` -- same as the 3.00 version but with manual curation


## Sources of information

- http://weatherby.genetics.utah.edu/MAKER/wiki/index.php/MAKER_Tutorial_for_GMOD_Online_Training_2014#MAKER.27s_Output
- https://isugenomics.github.io/bioinformatics-workbook/dataAnalysis/GenomeAnnotation/Intro_To_Maker.html
- http://weatherby.genetics.utah.edu/MAKER/wiki/index.php/The_MAKER_control_files_explained
- https://www.ebi.ac.uk/about/vertebrate-genomics/software/exonerate-user-guide
- https://www.ncbi.nlm.nih.gov/genbank/genomes_gff/

