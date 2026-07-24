# Genomic impact of the second plague pandemic on three human populations

This repository contains custom scripts used for the analyses described in
Liu et al., "Genomic impact of the second plague pandemic on three human
populations."

The study analyses whole-genome data from 529 ancient individuals from Lund,
Trondheim and Vilnius, together with present-day reference datasets.

## Repository contents

- `scripts/selectionScan/`: scripts for genome-wide association analyses, meta-analyses, and generation of LocusZoom plots
- `scripts/power/`: statistical power simulations for the selection scans
- `scripts/mahalanobisD/`: scripts for calculating PCA-based Mahalanobis distances
- `scripts/ibd/`: scripts for calling IBD fragments between ancient and present-day individuals using an algorithm similar to IBIS, as well as scripts for fragment filtering and summarization
- `scripts/ancestry_classification/`: ancestry classification based on Mahalanobis distances and IBD sharing with present-day reference populations
- `scripts/plotting/`: scripts used to generate the main manuscript figures
- `demo/`: a small simulated dataset and demonstration workflow


## 1. System requirements

### Operating systems

The analyses were performed and tested on:

- Red Hat Enterprise Linux V8.10
- R 4.2.2
- Python 3.9.9

The code has not been tested on Windows or macOS unless otherwise stated.

### Software dependencies

List each program and version, for example:

- R 4.2.2
- Python 3.9.9
- GLIMPSE2 
- PLINK 1.9/2.0
- bcftools 1.16
- GEMMA 0.98.5
- METAL version released on 2011-03-25

R and Python package dependencies are listed in the individual script documentation.

### Hardware requirements

The demonstration can be run on a standard Linux desktop computer.

The complete analyses used a high-performance computing environment because
of the size of the whole-genome and reference-panel datasets. Memory and CPU
requirements vary by analysis and are described in the relevant directories.

## 2. Installation guide

Clone the repository:

```bash
git clone https://github.com/ivanliu3/second_plague_pandemic_genomics.git
cd second_plague_pandemic_genomics
```

### Typical installation time

Cloning the repository takes less than 1 minute on a standard desktop computer
with a typical internet connection. Installing all required third-party software
and R/Python dependencies typically takes approximately 30–60 minutes, depending
on the operating system, internet speed, and whether some dependencies are already
installed.

## 3. Demo

The `demo/` directory contains example PCA data for ancient test individuals and present-day reference populations. The demo calculates PCA-based Mahalanobis distances using the first seven principal components.

### Instructions to run

From the top-level repository directory, run:

```bash
bash demo/run_demo.sh
```

Alternatively, from the `demo/` directory, run:

```bash
./run_demo.sh
```

The demo requires Python 3, NumPy, and SciPy.

### Expected output

The demo generates:

* `data.mahadist.csv`: Mahalanobis distances between ancient test individuals and reference populations
* `data.pvalue.csv`: p-values associated with the Mahalanobis-distance results

### Expected runtime

The demo completed in 3.85  seconds on a standard macOS desktop using Python 3.12. It is expected to complete in less than 1 minute on a normal desktop computer and requires no specialised hardware.

