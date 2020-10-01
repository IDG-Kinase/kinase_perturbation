#!/bin/bash

#DMSO
salmon quant --validateMappings -i ../indexes/Homo_sapiens.GRCh38.cdna.all_index.release-95 -l A -r Chinmaya_HCC1806_batch_02/HCC1806_DMSO_R1_S6_R1_001.fastq.gz -o ../salmon_alignments/HCC1806/batch_02/rep_01/DMSO

salmon quant --validateMappings -i ../indexes/Homo_sapiens.GRCh38.cdna.all_index.release-95 -l A -r Chinmaya_HCC1806_batch_02/HCC1806_DMSO_R2_S12_R1_001.fastq.gz -o ../salmon_alignments/HCC1806/batch_02/rep_02/DMSO

#CoCl2
salmon quant --validateMappings -i ../indexes/Homo_sapiens.GRCh38.cdna.all_index.release-95 -l A -r Chinmaya_HCC1806_batch_02/HCC1806_CoCl2_R1_S1_R1_001.fastq.gz -o ../salmon_alignments/HCC1806/batch_02/rep_01/CoCl2

salmon quant --validateMappings -i ../indexes/Homo_sapiens.GRCh38.cdna.all_index.release-95 -l A -r Chinmaya_HCC1806_batch_02/HCC1806_CoCl2_R2_S7_R1_001.fastq.gz -o ../salmon_alignments/HCC1806/batch_02/rep_02/CoCl2

#Palbo
salmon quant --validateMappings -i ../indexes/Homo_sapiens.GRCh38.cdna.all_index.release-95 -l A -r Chinmaya_HCC1806_batch_02/HCC1806_Palbo_R1_S2_R1_001.fastq.gz -o ../salmon_alignments/HCC1806/batch_02/rep_01/Palbo

salmon quant --validateMappings -i ../indexes/Homo_sapiens.GRCh38.cdna.all_index.release-95 -l A -r Chinmaya_HCC1806_batch_02/HCC1806_Palbo_R2_S8_R1_001.fastq.gz -o ../salmon_alignments/HCC1806/batch_02/rep_02/Palbo

#SGC
salmon quant --validateMappings -i ../indexes/Homo_sapiens.GRCh38.cdna.all_index.release-95 -l A -r Chinmaya_HCC1806_batch_02/HCC1806_SGC_R1_S5_R1_001.fastq.gz -o ../salmon_alignments/HCC1806/batch_02/rep_01/SGC

salmon quant --validateMappings -i ../indexes/Homo_sapiens.GRCh38.cdna.all_index.release-95 -l A -r Chinmaya_HCC1806_batch_02/HCC1806_SGC_R2_S11_R1_001.fastq.gz -o ../salmon_alignments/HCC1806/batch_02/rep_02/SGC

#Tram_100
salmon quant --validateMappings -i ../indexes/Homo_sapiens.GRCh38.cdna.all_index.release-95 -l A -r Chinmaya_HCC1806_batch_02/HCC1806_Tram100_R1_S3_R1_001.fastq.gz -o ../salmon_alignments/HCC1806/batch_02/rep_01/Tram_100

salmon quant --validateMappings -i ../indexes/Homo_sapiens.GRCh38.cdna.all_index.release-95 -l A -r  Chinmaya_HCC1806_batch_02/HCC1806_Tram100_R2_S9_R1_001.fastq.gz -o ../salmon_alignments/HCC1806/batch_02/rep_02/Tram_100

#Tram_500
salmon quant --validateMappings -i ../indexes/Homo_sapiens.GRCh38.cdna.all_index.release-95 -l A -r Chinmaya_HCC1806_batch_02/HCC1806_Tram500_R1_S4_R1_001.fastq.gz -o ../salmon_alignments/HCC1806/batch_02/rep_01/Tram_500

salmon quant --validateMappings -i ../indexes/Homo_sapiens.GRCh38.cdna.all_index.release-95 -l A -r Chinmaya_HCC1806_batch_02/HCC1806_Tram500_R2_S10_R1_001.fastq.gz -o ../salmon_alignments/HCC1806/batch_02/rep_02/Tram_500
