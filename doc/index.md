---
title       : Operational Protein Families
subtitle    :
author      : Kathryn Iverson
job         :
framework   : io2012        # {io2012, html5slides, shower, dzslides, ...}
highlighter : highlight.js  # {highlight.js, prettify, highlight}
hitheme     : tomorrow      #
widgets     : []            # {mathjax, quiz, bootstrap}
mode        : standalone # {standalone, draft, selfcontained}
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

---

## Assemble



### Iterative velvet

1. Assemble with a high k value
1. Enable read tracking and save unused reads
1. Assemble unused reads at a lower k value
1. Repeat



### Megahit

* Runs an iterative assembly by default, small k to large
* Much faster than iterative velvet

---

## Gene prediction

Metagene annotator (MGA)

---

## Get counts

Map reads with bowtie (and normalize)

---

## BLAST

All v all

Blast all genes against each other

---

## Cluster

mg-cluster with mothur

---

## Counts

* Use bowtie to map reads to genes
* Normalize by sequence length

---

## Clustering algos

* neighbor algorithms
* k-means

---

## Annotation

Majority vote

---

# HMP dataset

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
