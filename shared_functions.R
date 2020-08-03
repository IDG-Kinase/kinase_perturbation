run_DESeq2 <- function(drug_perturb_files) {
	drug_perturb_files = drug_perturb_files %>%
		mutate(names = paste0(treatment,"_",batch,"_",rep))
	
	drug_perturb_files$batch = as.factor(drug_perturb_files$batch)
	drug_perturb_files$rep = as.factor(drug_perturb_files$rep)
	
	drug_perturb_exp = summarizeToGene(tximeta(drug_perturb_files))
	dds <- DESeqDataSet(drug_perturb_exp, design = ~treatment)
	keep <- rowSums(counts(dds)) > 1
	dds <- dds[keep,]
	
	dds_analysis <- DESeq(dds)
	dds_results = results(dds_analysis)
	
	dds_results_filtered = as.data.frame(dds_results)
	dds_results_filtered$ensembl_gene_id = rownames(dds_results)
	dds_results_filtered = dds_results_filtered %>%
		filter(abs(log2FoldChange) >= 1, padj <= 0.05) %>%
		left_join(gene_to_hgnc) %>%
		filter(!is.na(hgnc_symbol))
	
	dds_results_filtered_kinase = dds_results_filtered %>%
		filter(hgnc_symbol %in% all_kinases$symbol)
	
	dds_results_filtered_dark_kinase = dds_results_filtered %>%
		filter(hgnc_symbol %in% dark_kinases$symbol)
	
	return(data.frame(
		number_2fold = dim(dds_results_filtered)[1],
		num_kinases = dim(dds_results_filtered_kinase)[1],
		num_dark_kinases = dim(dds_results_filtered_dark_kinase)[1],
		dk_diff_str = ifelse(dim(dds_results_filtered_dark_kinase)[1] == 0,
												 "",
												 glue_collapse(sort(dds_results_filtered_dark_kinase$hgnc_symbol), sep=" | "))
	))
}