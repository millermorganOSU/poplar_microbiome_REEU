---
editor_options: 
  markdown: 
    wrap: 72
---

# Poplar Microbiome Project (WIP)

REEU Big Data Internship

Morgan Miller

Objective of the project:

"...characterize the fungal wood microbiome of *P. trichocarpa* using
ITS metabarcoding and quantitative molecular approaches..."

(from the NSF RAPID proposal to PBI document)

General Workflow:

Raw Data:

Files containing reads from 528 poplar living-wood samples. The ITS2
region of the DNA fragments in the original samples was amplified with
the ITS4 and ITS3_KY01 primers (Sanger) before being sequenced (Illumina
high-throughput).

Quality Assessment using Fast-QC:

Remove Primers Using Cutadapt.

Read Processing using DADA2:

-   Filter and Trim Reads

-   

What is ITS metabarcoding and how can it be used to characterize a
microbiome?

The ITS (Internal Transcribed Spacer) region is a non-coding region
within the structural-ribosomal RNA gene.

Benefits of the ITS regions for fungal metabarcoding:

ITS region has high Inter-specific variability (mutations in this region
are not selected apron)

Flanked by highly conserved regions (structral-ribosomal RNA genes)
useful for primers and amplification
