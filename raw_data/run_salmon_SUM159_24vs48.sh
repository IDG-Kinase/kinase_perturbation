#!/bin/bash

#DMSO
salmon quant --validateMappings -i ../indexes/Homo_sapiens.GRCh38.cdna.all_index.release-95 -l A -r Mike_48_hr/DMSO_1_R1_S1_R1_001.fastq.gz -o ../salmon_alignments/24hr_vs_48hr/24_hr/rep_01/DMSO

salmon quant --validateMappings -i ../indexes/Homo_sapiens.GRCh38.cdna.all_index.release-95 -l A -r Mike_48_hr/DMSO_1_R2_S7_R1_001.fastq.gz -o ../salmon_alignments/24hr_vs_48hr/24_hr/rep_02/DMSO

salmon quant --validateMappings -i ../indexes/Homo_sapiens.GRCh38.cdna.all_index.release-95 -l A -r Mike_48_hr/DMSO_2_R1_S2_R1_001.fastq.gz -o ../salmon_alignments/24hr_vs_48hr/48_hr/rep_01/DMSO

salmon quant --validateMappings -i ../indexes/Homo_sapiens.GRCh38.cdna.all_index.release-95 -l A -r Mike_48_hr/DMSO_2_R2_S8_R1_001.fastq.gz -o ../salmon_alignments/24hr_vs_48hr/48_hr/rep_02/DMSO

#Alis
salmon quant --validateMappings -i ../indexes/Homo_sapiens.GRCh38.cdna.all_index.release-95 -l A -r Mike_48_hr/Alis_24_R1_S3_R1_001.fastq.gz -o ../salmon_alignments/24hr_vs_48hr/24_hr/rep_01/Alis

salmon quant --validateMappings -i ../indexes/Homo_sapiens.GRCh38.cdna.all_index.release-95 -l A -r Mike_48_hr/Alis_24_R2_S9_R1_001.fastq.gz -o ../salmon_alignments/24hr_vs_48hr/24_hr/rep_02/Alis

salmon quant --validateMappings -i ../indexes/Homo_sapiens.GRCh38.cdna.all_index.release-95 -l A -r Mike_48_hr/Alis_48_R1_S4_R1_001.fastq.gz -o ../salmon_alignments/24hr_vs_48hr/48_hr/rep_01/Alis

salmon quant --validateMappings -i ../indexes/Homo_sapiens.GRCh38.cdna.all_index.release-95 -l A -r Mike_48_hr/Alis_48_R2_S10_R1_001.fastq.gz -o ../salmon_alignments/24hr_vs_48hr/48_hr/rep_02/Alis

#Palbo
salmon quant --validateMappings -i ../indexes/Homo_sapiens.GRCh38.cdna.all_index.release-95 -l A -r Mike_48_hr/Palbo_24_R1_S5_R1_001.fastq.gz -o ../salmon_alignments/24hr_vs_48hr/24_hr/rep_01/Palbo

salmon quant --validateMappings -i ../indexes/Homo_sapiens.GRCh38.cdna.all_index.release-95 -l A -r Mike_48_hr/Palbo_24_R2_S11_R1_001.fastq.gz -o ../salmon_alignments/24hr_vs_48hr/24_hr/rep_02/Palbo

salmon quant --validateMappings -i ../indexes/Homo_sapiens.GRCh38.cdna.all_index.release-95 -l A -r Mike_48_hr/Palbo_48_R1_S6_R1_001.fastq.gz -o ../salmon_alignments/24hr_vs_48hr/48_hr/rep_01/Palbo

salmon quant --validateMappings -i ../indexes/Homo_sapiens.GRCh38.cdna.all_index.release-95 -l A -r Mike_48_hr/Palbo_48_R2_S12_R1_001.fastq.gz -o ../salmon_alignments/24hr_vs_48hr/48_hr/rep_02/Palbo


