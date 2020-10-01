#!/bin/bash

#DMSO
salmon quant --validateMappings -i ../indexes/Homo_sapiens.GRCh38.cdna.all_index.release-95 -l A -r Chinmaya_HCC1806/HCC1806_DMSO_R1_S6_R1_001.fastq.gz -o ../salmon_alignments/HCC1806/batch_01/rep_01/DMSO

salmon quant --validateMappings -i ../indexes/Homo_sapiens.GRCh38.cdna.all_index.release-95 -l A -r Chinmaya_HCC1806/HCC1806_DMSO_R2_S12_R1_001.fastq.gz -o ../salmon_alignments/HCC1806/batch_01/rep_02/DMSO

#Bort
salmon quant --validateMappings -i ../indexes/Homo_sapiens.GRCh38.cdna.all_index.release-95 -l A -r Chinmaya_HCC1806/HCC1806_Bort_R1_S2_R1_001.fastq.gz -o ../salmon_alignments/HCC1806/batch_01/rep_01/Bort

salmon quant --validateMappings -i ../indexes/Homo_sapiens.GRCh38.cdna.all_index.release-95 -l A -r Chinmaya_HCC1806/HCC1806_Bort_R2_S8_R1_001.fastq.gz -o ../salmon_alignments/HCC1806/batch_01/rep_02/Bort

#Ent
salmon quant --validateMappings -i ../indexes/Homo_sapiens.GRCh38.cdna.all_index.release-95 -l A -r Chinmaya_HCC1806/HCC1806_Ent_R1_S3_R1_001.fastq.gz -o ../salmon_alignments/HCC1806/batch_01/rep_01/Ent

salmon quant --validateMappings -i ../indexes/Homo_sapiens.GRCh38.cdna.all_index.release-95 -l A -r  Chinmaya_HCC1806/HCC1806_Ent_R2_S9_R1_001.fastq.gz -o ../salmon_alignments/HCC1806/batch_01/rep_02/Ent

#MK2206
salmon quant --validateMappings -i ../indexes/Homo_sapiens.GRCh38.cdna.all_index.release-95 -l A -r Chinmaya_HCC1806/HCC1806_MK2206_R1_S4_R1_001.fastq.gz -o ../salmon_alignments/HCC1806/batch_01/rep_01/MK2206

salmon quant --validateMappings -i ../indexes/Homo_sapiens.GRCh38.cdna.all_index.release-95 -l A -r Chinmaya_HCC1806/HCC1806_MK2206_R2_S10_R1_001.fastq.gz -o ../salmon_alignments/HCC1806/batch_01/rep_02/MK2206

#Navit
salmon quant --validateMappings -i ../indexes/Homo_sapiens.GRCh38.cdna.all_index.release-95 -l A -r Chinmaya_HCC1806/HCC1806_Navit_R1_S5_R1_001.fastq.gz -o ../salmon_alignments/HCC1806/batch_01/rep_01/Navit

salmon quant --validateMappings -i ../indexes/Homo_sapiens.GRCh38.cdna.all_index.release-95 -l A -r Chinmaya_HCC1806/HCC1806_Navit_R2_S11_R1_001.fastq.gz -o ../salmon_alignments/HCC1806/batch_01/rep_02/Navit

#Tram
salmon quant --validateMappings -i ../indexes/Homo_sapiens.GRCh38.cdna.all_index.release-95 -l A -r Chinmaya_HCC1806/HCC1806_Tram_R1_S1_R1_001.fastq.gz -o ../salmon_alignments/HCC1806/batch_01/rep_01/Tram

salmon quant --validateMappings -i ../indexes/Homo_sapiens.GRCh38.cdna.all_index.release-95 -l A -r Chinmaya_HCC1806/HCC1806_Tram_R2_S7_R1_001.fastq.gz -o ../salmon_alignments/HCC1806/batch_01/rep_02/Tram

