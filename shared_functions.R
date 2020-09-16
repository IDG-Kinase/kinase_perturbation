collect_all_DESeq <- function(exp_info, ...) {
	results = list()
	raw_DESeq_results = get_full_DESeq_results(exp_info, ...)
	results$DESeq_full_results = raw_DESeq_results %>% 
		as.data.frame() %>% 
		rownames_to_column(var = "ensembl_gene_id") %>% 
		left_join(gene_to_hgnc) %>% 
		select("hgnc_symbol",everything())
	
	results$DESeq_summary = run_DESeq2(raw_DESeq_results, ...)
	
	return(results)
}

run_DESeq2 <- function(dds_results, ...) {
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

get_full_DESeq_results <- function(drug_perturb_files, comparison_treatment = "DMSO") {
	drug_perturb_files = drug_perturb_files %>%
		mutate(names = paste0(treatment,"_",batch,"_",rep)) %>%
		#ensure that DMSO is the first level, forcing the fold change values to be
		#calculated as drug_treatment/DMSO
		mutate(treatment = relevel(as.factor(treatment), comparison_treatment),
					 batch = as.factor(batch),
					 rep = as.factor(rep)) %>%
		identity()
	
	drug_perturb_files$batch = as.factor(drug_perturb_files$batch)
	drug_perturb_files$rep = as.factor(drug_perturb_files$rep)
	
	drug_perturb_exp = summarizeToGene(tximeta(drug_perturb_files))
	dds <- DESeqDataSet(drug_perturb_exp, design = ~treatment)
	
	#Code to deal with collapsing technical replicates into single columns in
	#DESeq, will only be used if a "run" column is defined, which should be a
	#unique name for each run, which collapses together in the names field.
	if ("run" %in% colnames(drug_perturb_files)) {
		dds <- collapseReplicates(dds, dds$names, dds$run)
	}
	keep <- rowSums(counts(dds)) > 1
	dds <- dds[keep,]
	
	dds_analysis <- DESeq(dds)
	dds_results = results(dds_analysis)
	dds_results = lfcShrink(dds_analysis, coef=resultsNames(dds_analysis)[2], type="apeglm")
	
	return(dds_results)
}

