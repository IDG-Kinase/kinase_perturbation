#!/bin/bash

salmon quant --validateMappings -i ../indexes/Homo_sapiens.GRCh38.cdna.all_index.release-95 -l A -r MDAMB231/DMSO_R2_S12_R1_001.fastq.gz -o ../salmon_alignments/MDAMB231/rep_01/DMSO

salmon quant --validateMappings -i ../indexes/Homo_sapiens.GRCh38.cdna.all_index.release-95 -l A -r MDAMB231/DMSO_S6_R1_001.fastq.gz -o ../salmon_alignments/MDAMB231/rep_02/DMSO


salmon quant --validateMappings -i ../indexes/Homo_sapiens.GRCh38.cdna.all_index.release-95 -l A -r MDAMB231/Bort_R2_S11_R1_001.fastq.gz -o ../salmon_alignments/MDAMB231/rep_01/Bort

salmon quant --validateMappings -i ../indexes/Homo_sapiens.GRCh38.cdna.all_index.release-95 -l A -r MDAMB231/Bort_S5_R1_001.fastq.gz -o ../salmon_alignments/MDAMB231/rep_02/Bort


salmon quant --validateMappings -i ../indexes/Homo_sapiens.GRCh38.cdna.all_index.release-95 -l A -r MDAMB231/JIB04_R2_S9_R1_001.fastq.gz -o ../salmon_alignments/MDAMB231/rep_01/JIB04

salmon quant --validateMappings -i ../indexes/Homo_sapiens.GRCh38.cdna.all_index.release-95 -l A -r MDAMB231/JIB04_S3_R1_001.fastq.gz -o ../salmon_alignments/MDAMB231/rep_02/JIB04


salmon quant --validateMappings -i ../indexes/Homo_sapiens.GRCh38.cdna.all_index.release-95 -l A -r MDAMB231/MK_R2_S8_R1_001.fastq.gz -o ../salmon_alignments/MDAMB231/rep_01/MK2206

salmon quant --validateMappings -i ../indexes/Homo_sapiens.GRCh38.cdna.all_index.release-95 -l A -r MDAMB231/MK_S2_R1_001.fastq.gz -o ../salmon_alignments/MDAMB231/rep_02/MK2206


salmon quant --validateMappings -i ../indexes/Homo_sapiens.GRCh38.cdna.all_index.release-95 -l A -r MDAMB231/Palbo_R2_S10_R1_001.fastq.gz -o ../salmon_alignments/MDAMB231/rep_01/Palbo

salmon quant --validateMappings -i ../indexes/Homo_sapiens.GRCh38.cdna.all_index.release-95 -l A -r MDAMB231/Palbo_S4_R1_001.fastq.gz -o ../salmon_alignments/MDAMB231/rep_02/Palbo


salmon quant --validateMappings -i ../indexes/Homo_sapiens.GRCh38.cdna.all_index.release-95 -l A -r MDAMB231/Tram_R2_S7_R1_001.fastq.gz -o ../salmon_alignments/MDAMB231/rep_01/Tram

salmon quant --validateMappings -i ../indexes/Homo_sapiens.GRCh38.cdna.all_index.release-95 -l A -r MDAMB231/Tram_S1_R1_001.fastq.gz -o ../salmon_alignments/MDAMB231/rep_02/Tram
