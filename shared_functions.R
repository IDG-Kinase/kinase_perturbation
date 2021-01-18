collect_all_DESeq <- function(exp_info, ...) {
	results = list()
	raw_DESeq_results = get_full_DESeq_results(exp_info, ...)
	results$DESeq_full_results = raw_DESeq_results %>% 
		as.data.frame() %>% 
		rownames_to_column(var = "ensembl_gene_id") %>% 
		left_join(gene_to_hgnc) %>% 
		select("hgnc_symbol",everything())
	
	filtered_results = filter_DESeq(results$DESeq_full_results)
	
	results$DESeq_summary = filter_DESeq(filtered_results, ...)
	
	return(results)
}

filter_DESeq <- function(dds_results) {
	dds_results_filtered = dds_results %>%
		filter(abs(log2FoldChange) >= 1, padj <= 0.05) %>%
		filter(!is.na(hgnc_symbol))
	
	return(dds_results_filtered)
}

summarize_DESeq_results <- function(dds_results_filtered) {
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

collect_full_DESeq_and_summary <- function(dds_analysis) {
	full_DESeq = tribble(~hgnc_symbol,~ensembl_gene_id, ~baseMean, ~log2FoldChange, ~lfcSE, ~stat, ~pvalue, ~padj, ~treatment)
	summary_stats = tribble(~treatment,~number_2fold,~num_kinases,~num_dark_kinases,~dk_diff_str)
	
	for (this_comparison in resultsNames(dds_analysis)[-1]) {
		dds_results = lfcShrink(dds_analysis, coef=this_comparison, type="apeglm")
		
		full_DESeq = rbind(
			full_DESeq,
			dds_results %>%
				as.data.frame() %>%
				rownames_to_column(var = "hgnc_symbol") %>%
				mutate(treatment = this_comparison) %>%
				dplyr::select("hgnc_symbol","treatment",everything())
		)
		
		dds_results_filtered = dds_results %>%
			as.data.frame() %>%
			rownames_to_column(var = "hgnc_symbol") %>%
			filter(abs(log2FoldChange) >= 1, padj <= 0.05) %>%
			mutate("treatment" = this_comparison)
		
		dds_results_filtered_kinase = dds_results_filtered %>%
			filter(hgnc_symbol %in% all_kinases$symbol)
		
		dds_results_filtered_dark_kinase = dds_results_filtered %>%
			filter(hgnc_symbol %in% dark_kinases$symbol)
		
		summary_stats = rbind(
			summary_stats,
			data.frame(
				treatment = this_comparison,
				number_2fold = dim(dds_results_filtered)[1],
				num_kinases = dim(dds_results_filtered_kinase)[1],
				num_dark_kinases = dim(dds_results_filtered_dark_kinase)[1],
				dk_diff_str = ifelse(dim(dds_results_filtered_dark_kinase)[1] == 0,
														 "",
														 glue_collapse(sort(dds_results_filtered_dark_kinase$hgnc_symbol), sep=" | "))
			)
		)
	}
	
	return(list(full = full_DESeq, summary = summary_stats))
}

get_DESeq_analysis <- function(drug_perturb_files, 
															 comparison_treatment = "DMSO",
															 DESeq_design = ~treatment) {
	
	drug_perturb_files = exp_info %>%
		mutate(names = paste0(treatment,"_",batch,"_",rep)) %>%
		#ensure that DMSO is the first level, forcing the fold change values to be
		#calculated as drug_treatment/DMSO
		mutate(treatment = relevel(as.factor(treatment), comparison_treatment),
					 batch = as.factor(batch),
					 rep = as.factor(rep)) %>%
		identity()
	
	drug_perturb_files$batch = as.factor(drug_perturb_files$batch)
	drug_perturb_files$rep = as.factor(drug_perturb_files$rep)
	
	transcript_to_hgnc = get_transcript_to_hgnc()
	
	mytx2gene = transcript_to_hgnc %>%
		select(ensembl_transcript_id_version, hgnc_symbol)
	
	imported_tx <- tximport(exp_info$files, 
													type = "salmon",
													tx2gene = mytx2gene)
	
	dds <- DESeqDataSetFromTximport(imported_tx,
																	colData = drug_perturb_files %>% select(-files),
																	design = DESeq_design)
	
	#Toss out the rownames where all the ENSTs without HGNCs are grouped:
	dds = dds[rownames(dds) != "",]
	
	keep <- rowSums(counts(dds)) > 1
	dds <- dds[keep,]
	
	dds_analysis <- DESeq(dds)
}


get_transcript_to_hgnc <- function() {
	if (! file.exists(here('transcript_to_hgnc.csv.gz'))) {
		ensembl = biomaRt::useEnsembl("genes",dataset="hsapiens_gene_ensembl", version=95)
		transcript_to_hgnc = biomaRt::getBM(attributes = c('ensembl_gene_id', 
																											 'ensembl_transcript_id',
																											 'ensembl_transcript_id_version',
																											 'hgnc_symbol'), 
																				mart = ensembl)
		
		transcript_to_hgnc = transcript_to_hgnc %>% 
			#There are a small number of transcripts with multiple HGNC hits, I'm going to
			#select only one
			group_by(ensembl_transcript_id_version) %>% 
			filter(row_number()==1) %>% 
			ungroup()
		
		write_csv(transcript_to_hgnc, here('transcript_to_hgnc.csv.gz'))
	} else {
		transcript_to_hgnc = readr::read_csv(here('transcript_to_hgnc.csv.gz'), 
																				 col_types = cols(
																				 	ensembl_gene_id = col_character(),
																				 	ensembl_transcript_id = col_character(),
																				 	ensembl_transcript_id_version = col_character(),
																				 	hgnc_symbol = col_character()
																				 )) %>%
			#The DESeq transcript to gene level summarization doesn't like NA values in
			#the gene level targets, so I'm going to convert all the NA to "", so they
			#can be grouped and tossed out later
			mutate(hgnc_symbol = ifelse(is.na(hgnc_symbol), "", hgnc_symbol))
		
	}
}

gather_gene_TPM_values <- function(exp_info) {
	transcript_to_hgnc = get_transcript_to_hgnc()
	
	full_TPM_data = data.frame(hgnc_symbol = sort(unique(transcript_to_hgnc$hgnc_symbol))) %>%
		filter(hgnc_symbol != "")
	
	for (i in 1:dim(exp_info)[1]) {
		data_set_name = paste0(exp_info$treatment[i], "-batch_",exp_info$batch[i],"-rep_",exp_info$rep[i])
		this_salmon_data = read_tsv(here(exp_info$files[i]), 
																col_types = cols(
																	Name = col_character(),
																	Length = col_double(),
																	EffectiveLength = col_double(),
																	TPM = col_double(),
																	NumReads = col_double()
																))
		
		this_salmon_data = this_salmon_data %>% 
			left_join(transcript_to_hgnc %>% select(ensembl_transcript_id_version,hgnc_symbol), 
								by = c('Name' = 'ensembl_transcript_id_version')) %>%
			filter(hgnc_symbol != "") %>%
			group_by(hgnc_symbol) %>%
			summarise(TPM = sum(TPM)) %>%
			rename(!!data_set_name := TPM)
		
		full_TPM_data = full_TPM_data %>% left_join(this_salmon_data, by='hgnc_symbol')
	}
	
	return(full_TPM_data)
}