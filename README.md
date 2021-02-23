# RNAseq-Based Compound Perturbation Analysis


This repository contains the code used to analyse the compound perturbation RNAseq data. The processing starts with the sequencing files and then uses salmon and a pre-build index to align the sequences to the transcriptome (see salmon commands [here](https://github.com/IDG-Kinase/kinase_perturbation/tree/master/raw_data)). From there, these files are fed into DESeq2 for quantification (see a sample [script](https://github.com/IDG-Kinase/kinase_perturbation/tree/master/analysis/SUM159_DESeq_run) for the SUM159 cells).
