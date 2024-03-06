# PodosporaGenomes

Here you'll find the code associated with the paper:

Ament-Velásquez et al. (2024) High-quality genome assemblies of four members of the *Podospora anserina* species complex, *Genome Biology and Evolution*: evae034, [https://doi.org/10.1093/gbe/evae034](https://doi.org/10.1093/gbe/evae034)

The [Snakemake](https://snakemake.readthedocs.io/en/stable/) pipelines follow this orde: 

- 1_Annotation_v3 -- the de novo genome annotation with MAKER
- 2_FunctionalAnnotation -- functional annotation with Funannotate
- 3_CircosSppComplex -- Fig. 1
- 4_OrthoTrees -- Getting the single-copy orthologs
- 5_PodoPhylogeny -- the phylogenomic analyses for Fig. 2

I ran the pipelines in Uppsala University's supercomputer [UPPMAX](https://uppmax.uu.se/), which has a CentOS Linux operating system with a slurm scheduler. However, they should work fine also in other unix environments.

----

Disclaimer: These scripts and files are provided "as is" and without any express or implied warranties, including, without limitation, the implied warranties of merchantability and fitness for a particular purpose.
