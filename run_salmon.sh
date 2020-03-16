#!/bin/bash

salmon quant --validateMappings -i indexes/Homo_sapiens.GRCh38.cdna.all_index -l A -r raw_data/11-14-2019/SUM159_1_DMSO_R1_S1_R1_001.fastq.gz -o salmon_alignments/batch_01/rep_01/DMSO

salmon quant --validateMappings -i indexes/Homo_sapiens.GRCh38.cdna.all_index -l A -r raw_data/11-14-2019/SUM159_2_Tram_R1_S2_R1_001.fastq.gz -o salmon_alignments/batch_01/rep_01/Tram

salmon quant --validateMappings -i indexes/Homo_sapiens.GRCh38.cdna.all_index -l A -r raw_data/11-14-2019/SUM159_3_MK2206_R1_S3_R1_001.fastq.gz -o salmon_alignments/batch_01/rep_01/MK2206

salmon quant --validateMappings -i indexes/Homo_sapiens.GRCh38.cdna.all_index -l A -r raw_data/11-14-2019/SUM159_4_Dinac_R1_S4_R1_001.fastq.gz -o salmon_alignments/batch_01/rep_01/Dinac

salmon quant --validateMappings -i indexes/Homo_sapiens.GRCh38.cdna.all_index -l A -r raw_data/11-14-2019/SUM159_5_Ent_R1_S5_R1_001.fastq.gz -o salmon_alignments/batch_01/rep_01/Ent

salmon quant --validateMappings -i indexes/Homo_sapiens.GRCh38.cdna.all_index -l A -r raw_data/11-14-2019/SUM159_6_JIB04_R1_S6_R1_001.fastq.gz -o salmon_alignments/batch_01/rep_01/JIB04


salmon quant --validateMappings -i indexes/Homo_sapiens.GRCh38.cdna.all_index -l A -r raw_data/11-14-2019/SUM159_Dinac_R2_S10_R1_001.fastq.gz -o salmon_alignments/batch_01/rep_02/Dinac

salmon quant --validateMappings -i indexes/Homo_sapiens.GRCh38.cdna.all_index -l A -r raw_data/11-14-2019/SUM159_DMSO_R2_S7_R1_001.fastq.gz -o salmon_alignments/batch_01/rep_02/DMSO

salmon quant --validateMappings -i indexes/Homo_sapiens.GRCh38.cdna.all_index -l A -r raw_data/11-14-2019/SUM159_Ent_R2_S11_R1_001.fastq.gz -o salmon_alignments/batch_01/rep_02/Ent

salmon quant --validateMappings -i indexes/Homo_sapiens.GRCh38.cdna.all_index -l A -r raw_data/11-14-2019/SUM159_JIB04_R2_S12_R1_001.fastq.gz -o salmon_alignments/batch_01/rep_02/JIB04

salmon quant --validateMappings -i indexes/Homo_sapiens.GRCh38.cdna.all_index -l A -r raw_data/11-14-2019/SUM159_MK2206_R2_S9_R1_001.fastq.gz -o salmon_alignments/batch_01/rep_02/MK2206

salmon quant --validateMappings -i indexes/Homo_sapiens.GRCh38.cdna.all_index -l A -r raw_data/11-14-2019/SUM159_Tram_R2_S8_R1_001.fastq.gz -o salmon_alignments/batch_01/rep_02/Tram

salmon quant --validateMappings -i indexes/Homo_sapiens.GRCh38.cdna.all_index -l A -r raw_data/1-28-2020/SUM159_DMSO_R1_S6_R1_001.fastq.gz -o salmon_alignments/batch_02/rep_01/DMSO

salmon quant --validateMappings -i indexes/Homo_sapiens.GRCh38.cdna.all_index -l A -r raw_data/1-28-2020/SUM159_JIB_R1_S1_R1_001.fastq.gz -o salmon_alignments/batch_02/rep_01/JIB04

salmon quant --validateMappings -i indexes/Homo_sapiens.GRCh38.cdna.all_index -l A -r raw_data/1-28-2020/SUM159_Navi_R1_S5_R1_001.fastq.gz -o salmon_alignments/batch_02/rep_01/Navi

salmon quant --validateMappings -i indexes/Homo_sapiens.GRCh38.cdna.all_index -l A -r raw_data/1-28-2020/SUM159_Paclit_R1_S4_R1_001.fastq.gz -o salmon_alignments/batch_02/rep_01/Paclit

salmon quant --validateMappings -i indexes/Homo_sapiens.GRCh38.cdna.all_index -l A -r raw_data/1-28-2020/SUM159_Rux_R1_S2_R1_001.fastq.gz -o salmon_alignments/batch_02/rep_01/Rux

salmon quant --validateMappings -i indexes/Homo_sapiens.GRCh38.cdna.all_index -l A -r raw_data/1-28-2020/SUM159_SGC_R1_S3_R1_001.fastq.gz -o salmon_alignments/batch_02/rep_01/SGC

salmon quant --validateMappings -i indexes/Homo_sapiens.GRCh38.cdna.all_index -l A -r raw_data/1-28-2020/SUM159_Tram_R1_S13_R1_001.fastq.gz -o salmon_alignments/batch_02/rep_01/Tram


salmon quant --validateMappings -i indexes/Homo_sapiens.GRCh38.cdna.all_index -l A -r raw_data/1-28-2020/SUM159_DMSO_R2_S12_R1_001.fastq.gz -o salmon_alignments/batch_02/rep_02/DMSO

salmon quant --validateMappings -i indexes/Homo_sapiens.GRCh38.cdna.all_index -l A -r raw_data/1-28-2020/SUM159_Jib_R2_S8_R1_001.fastq.gz -o salmon_alignments/batch_02/rep_02/JIB04

salmon quant --validateMappings -i indexes/Homo_sapiens.GRCh38.cdna.all_index -l A -r raw_data/1-28-2020/SUM159_Navi_R2_S7_R1_001.fastq.gz -o salmon_alignments/batch_02/rep_02/Navi

salmon quant --validateMappings -i indexes/Homo_sapiens.GRCh38.cdna.all_index -l A -r raw_data/1-28-2020/SUM159_Pacli_R2_S9_R1_001.fastq.gz -o salmon_alignments/batch_02/rep_02/Paclit

salmon quant --validateMappings -i indexes/Homo_sapiens.GRCh38.cdna.all_index -l A -r raw_data/1-28-2020/SUM159_Rux_R2_S11_R1_001.fastq.gz -o salmon_alignments/batch_02/rep_02/Rux

salmon quant --validateMappings -i indexes/Homo_sapiens.GRCh38.cdna.all_index -l A -r raw_data/1-28-2020/SUM159_SGC_R2_S10_R1_001.fastq.gz -o salmon_alignments/batch_02/rep_02/SGC

salmon quant --validateMappings -i indexes/Homo_sapiens.GRCh38.cdna.all_index -l A -r  raw_data/1-22-2020/SUM159_DMSO_R1_S1_R1_001.fastq.gz -o salmon_alignments/batch_03/rep_01/DMSO

salmon quant --validateMappings -i indexes/Homo_sapiens.GRCh38.cdna.all_index -l A -r  raw_data/1-22-2020/SUM159_Alis_R1_S2_R1_001.fastq.gz -o salmon_alignments/batch_03/rep_01/Alis

salmon quant --validateMappings -i indexes/Homo_sapiens.GRCh38.cdna.all_index -l A -r  raw_data/1-22-2020/SUM159_Bort_R1_S3_R1_001.fastq.gz -o salmon_alignments/batch_03/rep_01/Bort

salmon quant --validateMappings -i indexes/Homo_sapiens.GRCh38.cdna.all_index -l A -r  raw_data/1-22-2020/SUM159_CoCl2_R1_S4_R1_001.fastq.gz -o salmon_alignments/batch_03/rep_01/CoCl2

salmon quant --validateMappings -i indexes/Homo_sapiens.GRCh38.cdna.all_index -l A -r  raw_data/1-22-2020/SUM159_FCCP_R1_S5_R1_001.fastq.gz -o salmon_alignments/batch_03/rep_01/FCCP

salmon quant --validateMappings -i indexes/Homo_sapiens.GRCh38.cdna.all_index -l A -r  raw_data/1-22-2020/SUM159_Palbo_R1_S6_R1_001.fastq.gz -o salmon_alignments/batch_03/rep_01/Palbo


salmon quant --validateMappings -i indexes/Homo_sapiens.GRCh38.cdna.all_index -l A -r  raw_data/1-22-2020/SUM159_DMSO_R2_S7_R1_001.fastq.gz -o salmon_alignments/batch_03/rep_02/DMSO

salmon quant --validateMappings -i indexes/Homo_sapiens.GRCh38.cdna.all_index -l A -r  raw_data/1-22-2020/SUM159_Alis_R2_S8_R1_001.fastq.gz -o salmon_alignments/batch_03/rep_02/Alis

salmon quant --validateMappings -i indexes/Homo_sapiens.GRCh38.cdna.all_index -l A -r  raw_data/1-22-2020/SUM159_Bort_R2_S9_R1_001.fastq.gz -o salmon_alignments/batch_03/rep_02/Bort

salmon quant --validateMappings -i indexes/Homo_sapiens.GRCh38.cdna.all_index -l A -r  raw_data/1-22-2020/SUM159_CoCl2_R2_S10_R1_001.fastq.gz -o salmon_alignments/batch_03/rep_02/CoCl2

salmon quant --validateMappings -i indexes/Homo_sapiens.GRCh38.cdna.all_index -l A -r  raw_data/1-22-2020/SUM159_FCCP_R2_S11_R1_001.fastq.gz -o salmon_alignments/batch_03/rep_02/FCCP

salmon quant --validateMappings -i indexes/Homo_sapiens.GRCh38.cdna.all_index -l A -r  raw_data/1-22-2020/SUM159_Palbo_R2_S12_R1_001.fastq.gz -o salmon_alignments/batch_03/rep_02/Palbo

salmon quant --validateMappings -i indexes/Homo_sapiens.GRCh38.cdna.all_index -l A -r  raw_data/1-22-2020/SUM159_Tram_R2_S13_R1_001.fastq.gz -o salmon_alignments/batch_03/rep_02/Tram

