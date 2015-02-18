---
title       : Operational Protein Families
subtitle    :
author      : Kathryn Iverson
job         : Bioinformagician
framework   : io2012       # {io2012, html5slides, shower, dzslides, ...}
highlighter : highlight.js  # {highlight.js, prettify, highlight}
hitheme     : tomorrow      #
widgets     : []            # {mathjax, quiz, bootstrap}
mode        : selfcontained # {standalone, draft, selfcontained}
knit        : slidify::knit2slides
---


## OPFs

Operational protein families are clusters of similar proteins at the aa sequence level

Creating OPFs is a database free method of grouping sequences

--- .class #id

## Quickly:

1. Assemble
1. Gene prediction
1. Get counts
1. All v all BLAST
1. Cluster

--- &twocol

## Assemble

*** =left
### Iterative velvet

1. Assemble with a high k value
1. Enable read tracking and save unused reads
1. Assemble unused reads at a lower k value
1. Repeat


*** =right
### Megahit

* Runs an iterative assembly by default, small k to large
* Much faster than iterative velvet

---

## Gene prediction

### Metagene annotator (MGA)

* Gene prediction models for Bacteria, Archaea and Phage
* Uses ribosomal binding site patterns as predictors for domain
* Optimised for short fragments
    - uses 700bp fragment as example
    - similar to contig sizes seen after assembly

---

## Estimate gene counts

Map reads with bowtie (and normalize)

---

## Counts file

* Use bowtie to map reads to genes
* Normalize by sequence length 

```
ceiling(num_reads_mapped * 100 / len(gene))
````

---

## BLAST

All v all

Blast all genes against each other

---

## Cluster

mg-cluster with mothur

---

## Clustering algos

* neighbor algorithms
* k-means

---

## Annotation

BLAST genes to a database of your choice, KEGG for example, and annotate OPF by majority vote.

---

# HMP dataset

Dataset was too large to cluster as a whole so OPFs were clustered within the KEGG category they 

---

## eggNOG

---

## KEGG

---

# Pregnancy study

---

## OPFs vs MG-RAST


---

## OPFs vs eggNOG


---

## OPFs vs KEGG


---
