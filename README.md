# RNA-seq-Model-Benchmarking
This is the comparison of the frequentist and bayesian approach on RNA-seq data

Here’s a clean, **GitHub README-style summary** of your project—high-level, clear, and not too detailed:

---

# RNA-seq Modeling Simulation (R)

## Overview

This project investigates how different statistical modeling choices impact the analysis of RNA-seq data, particularly in detecting **differentially expressed (DE) genes**. We focus on evaluating how assumptions about data distribution and independence affect key performance metrics.

The goal is to build a simulation framework in **R** that mimics realistic RNA-seq data and allows systematic comparison of commonly used models.

---

## Key Idea

RNA-seq data consists of **count measurements of gene expression** across samples. These counts are:

* Overdispersed (variance > mean)
* Potentially correlated across genes
* Influenced by biological differences (e.g., case vs control)

We simulate such data and test how well different models recover true biological signals.

---

## What We Do

### 1. Simulate RNA-seq Data

* 100 samples (50 case, 50 control)
* 1000 genes
* Counts generated using the **Negative Binomial distribution** (realistic for RNA-seq)
* 10% of genes are truly **differentially expressed (ground truth)** 

### 2. Two Simulation Settings

* **Independent genes**: standard assumption in many models
* **Dependent genes**: using a **Gaussian copula** to introduce realistic correlation structure 

### 3. Model Comparison

We implement and compare:

* **Gaussian models** (commonly used in practice, e.g., limma)
* **Negative Binomial models** (more appropriate for count data)

Both are implemented from scratch in **base R**.

### 4. Inference Frameworks

* **Frequentist**: GLMs for each gene
* **Bayesian**: same models with priors for parameters

---

## Objective

For each gene, we test whether expression differs between case and control groups.

We evaluate models based on:

* Power (ability to detect true DE genes)
* False discovery rate (FDR)
* Bias and statistical properties of estimators 

---

## Motivation

Many RNA-seq pipelines rely on simplifying assumptions such as:

* Gaussian approximations
* Independence between genes

This project explores:
How much these assumptions affect real conclusions
 Whether more realistic modeling (e.g., NB + dependency) improves inference

---

##  Implementation

* Language: **R**
* Approach: fully custom implementation (no high-level packages)
* Simulation repeated multiple times for robust evaluation

---

## Big Picture

This project bridges:

* **Statistics** → modeling assumptions and inference
* **Bioinformatics** → gene expression analysis
* **Computation** → simulation and algorithm implementation

Ultimately, we aim to better understand:

*Which models are reliable for detecting true biological signals in RNA-seq data*


