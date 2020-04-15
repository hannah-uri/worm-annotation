# worm-annotation
This is a Snakemake pipeline for annotating the genomes of wild C. elegans (and comparing their tRNA transcription rates).

Genomes and RNA-seq read files must be supplied in a directory structure as follows:
  
    genome:   `input/strains/{strain}/pilonN2n.fa`
    RNA-seq:  `input/strains/{strain}/r{replicate number}/rna/{read file}.fastq.gz`
    5'-dep small RNA-seq:  `input/strains/{strain}/r{replicate_number}/small-rna/5-dep/{read file}.fastq.gz`
    5'-indep small RNA-seq:  `input/strains/{strain}/r{replicate_number}/small-rna/5-indep/{read file}.fastq.gz`

To run, install [Conda](https://docs.conda.io/en/latest/) and [Snakemake](https://snakemake.readthedocs.io/en/stable/), then run:

    snakemake --use-conda
