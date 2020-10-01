#!/bin/bash

#DMSO
salmon quant --validateMappings -i ../indexes/Homo_sapiens.GRCh38.cdna.all_index.release-95 -l A -r Zawistowski_bulk_HCC1806/HCC1806_DMSO_24h_replicate1.fastq.gz -o ../salmon_alignments/Zawistowski_HCC1806/batch_01/rep_01/DMSO

salmon quant --validateMappings -i ../indexes/Homo_sapiens.GRCh38.cdna.all_index.release-95 -l A -r Zawistowski_bulk_HCC1806/HCC1806_DMSO_24h_replicate2.fastq.gz -o ../salmon_alignments/Zawistowski_HCC1806/batch_01/rep_02/DMSO

#Tramet
salmon quant --validateMappings -i ../indexes/Homo_sapiens.GRCh38.cdna.all_index.release-95 -l A -r Zawistowski_bulk_HCC1806/HCC1806_500nMtrametinib_24h_replicate1.fastq.gz -o ../salmon_alignments/Zawistowski_HCC1806/batch_01/rep_01/Tramet

salmon quant --validateMappings -i ../indexes/Homo_sapiens.GRCh38.cdna.all_index.release-95 -l A -r Zawistowski_bulk_HCC1806/HCC1806_500nMtrametinib_24h_replicate2.fastq.gz -o ../salmon_alignments/Zawistowski_HCC1806/batch_01/rep_02/Tramet

