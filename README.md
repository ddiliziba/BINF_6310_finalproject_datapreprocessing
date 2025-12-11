### Reproducibility Study: Improved Haplotype Inference from RNA-seq Data

This repository contains our reproduction of key components from the study:
Mehta et al. (2020). Improved haplotype inference by exploiting long-range linking and allelic imbalance in RNA-seq datasets. Nature Communications.
We aimed to replicate the computational workflow used to process RNA-seq data and prepare inputs for haplotype inference, validating the reproducibility of the original study’s methodology.

## Project Overview

The original study proposed a method that integrates:

Long-range RNA-seq read linking

Allelic imbalance patterns

Splice-aware alignment

Phasing across distal variants

Our project focuses on reproducing the data preprocessing and RNA-seq alignment pipeline, and downstream haplotype reconstruction using tools such as HapTree-X.

## Objectives

Reproduce the data acquisition, preprocessing, and alignment steps used in the study.

Validate the logic, accuracy, and reproducibility of the pipeline.

Identify discrepancies, limitations, or improvements compared to the original workflow.

Document all steps and provide a fully reproducible workflow.

## Workflow Summary
1. Data Acquisition

Downloaded GRCh38 primary assembly FASTA

Downloaded GENCODE v49 GTF annotations

Retrieved publicly available RNA-seq FASTQ files for testing

2. STAR Genome Indexing

We generated a splice-aware STAR index with GENCODE annotations:

STAR --runThreadN 8 \
     --runMode genomeGenerate \
     --genomeDir reference/STAR_GRCh38 \
     --genomeFastaFiles reference/GRCh38.primary_assembly.genome.fa \
     --sjdbGTFfile reference/gencode.v49.annotation.gtf \
     --sjdbOverhang 100

3. RNA-seq Read Alignment

Reads were aligned to the indexed genome:

STAR --runThreadN 8 \
     --genomeDir reference/STAR_GRCh38 \
     --readFilesIn reads/sample_R1.fastq.gz reads/sample_R2.fastq.gz \
     --readFilesCommand zcat \
     --outSAMtype BAM SortedByCoordinate \
     --outFileNamePrefix results/sample_


QC included:

Alignment rate

Uniquely mapped reads

Number of spliced reads

Junction usage statistics

5. Comparison With Original Study

We compared:

Mapping statistics

Junction detection

Coverage patterns

with values reported by Mehta et al. (2020).

## Tools & Versions
Tool	Version
STAR	2.7.11b
samtools	1.18
GENCODE	v49
GRCh38	Primary Assembly
Python	3.10
Linux / Apptainer	HPC environment

## Key Findings

The STAR index and alignment results were reproducible with the original study’s parameters.

Minor differences in mapping rates likely stem from updated genome/annotation versions.

Splice junction detection remained consistent, supporting the credibility of the original methodology.

## Limitations Identified

Computational resource differences between modern HPC and 2020 systems

Annotation version drift (v49 vs. original version) may alter junction counts

Some scripts referenced in the paper were not publicly available

Lack of explicit phasing evaluation metrics in the original methods

## Reference

Berger, E., Yorukoglu, D., Zhang, L. et al. Improved haplotype inference by exploiting long-range linking and allelic imbalance in RNA-seq datasets. Nat Commun 11, 4662 (2020). https://doi.org/10.1038/s41467-020-18320-z

## Acknowledgments

This project was completed as part of the BINF 6310 — Introduction to Computational Methods in Bioinformatics final group project at the Northeastern University.
