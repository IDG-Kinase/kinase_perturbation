#!/bin/bash

#DMSO
salmon quant --validateMappings -i ../indexes/Homo_sapiens.GRCh38.cdna.all_index.release-95 -l A -r Chinmaya_HCC1806_batch_03/HCC1806_DMSO_R1_S6_R1_001.fastq.gz -o ../salmon_alignments/HCC1806/batch_03/rep_01/DMSO

salmon quant --validateMappings -i ../indexes/Homo_sapiens.GRCh38.cdna.all_index.release-95 -l A -r Chinmaya_HCC1806_batch_03/HCC1806_DMSO_R2_S12_R1_001.fastq.gz -o ../salmon_alignments/HCC1806/batch_03/rep_02/DMSO

#Alisertib
salmon quant --validateMappings -i ../indexes/Homo_sapiens.GRCh38.cdna.all_index.release-95 -l A -r Chinmaya_HCC1806_batch_03/HCC1806_Alisertib_R1_S1_R1_001.fastq.gz -o ../salmon_alignments/HCC1806/batch_03/rep_01/Alisertib

salmon quant --validateMappings -i ../indexes/Homo_sapiens.GRCh38.cdna.all_index.release-95 -l A -r Chinmaya_HCC1806_batch_03/HCC1806_Alisertib_R2_S7_R1_001.fastq.gz -o ../salmon_alignments/HCC1806/batch_03/rep_02/Alisertib

#Dinaciclib
salmon quant --validateMappings -i ../indexes/Homo_sapiens.GRCh38.cdna.all_index.release-95 -l A -r Chinmaya_HCC1806_batch_03/HCC1806_Dinaciclib_R1_S3_R1_001.fastq.gz -o ../salmon_alignments/HCC1806/batch_03/rep_01/Dinaciclib

salmon quant --validateMappings -i ../indexes/Homo_sapiens.GRCh38.cdna.all_index.release-95 -l A -r Chinmaya_HCC1806_batch_03/HCC1806_Dinaciclib_R2_S9_R1_001.fastq.gz -o ../salmon_alignments/HCC1806/batch_03/rep_02/Dinaciclib

#JIB-04
salmon quant --validateMappings -i ../indexes/Homo_sapiens.GRCh38.cdna.all_index.release-95 -l A -r Chinmaya_HCC1806_batch_03/HCC1806_JIB04_R1_S4_R1_001.fastq.gz -o ../salmon_alignments/HCC1806/batch_03/rep_01/JIB04

salmon quant --validateMappings -i ../indexes/Homo_sapiens.GRCh38.cdna.all_index.release-95 -l A -r Chinmaya_HCC1806_batch_03/HCC1806_JIB04_R2_S10_R1_001.fastq.gz -o ../salmon_alignments/HCC1806/batch_03/rep_02/JIB04

#Paclitaxel
salmon quant --validateMappings -i ../indexes/Homo_sapiens.GRCh38.cdna.all_index.release-95 -l A -r Chinmaya_HCC1806_batch_03/HCC1806_Paclitaxel_R1_S2_R1_001.fastq.gz -o ../salmon_alignments/HCC1806/batch_03/rep_01/Paclitaxel

salmon quant --validateMappings -i ../indexes/Homo_sapiens.GRCh38.cdna.all_index.release-95 -l A -r Chinmaya_HCC1806_batch_03/HCC1806_Paclitaxel_R2_S8_R1_001.fastq.gz -o ../salmon_alignments/HCC1806/batch_03/rep_02/Paclitaxel

#Ruxolitinib
salmon quant --validateMappings -i ../indexes/Homo_sapiens.GRCh38.cdna.all_index.release-95 -l A -r Chinmaya_HCC1806_batch_03/HCC1806_Ruxolitinib_R1_S5_R1_001.fastq.gz -o ../salmon_alignments/HCC1806/batch_03/rep_01/Ruxolitinib

salmon quant --validateMappings -i ../indexes/Homo_sapiens.GRCh38.cdna.all_index.release-95 -l A -r Chinmaya_HCC1806_batch_03/HCC1806_Ruxolitinib_R2_S11_R1_001.fastq.gz -o ../salmon_alignments/HCC1806/batch_03/rep_02/Ruxolitinib
