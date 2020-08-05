library(here)

for (number_of_clusters in c(15,20,25,30)) {
	dir.create(here('analysis/DESeq_cluster_genome_display/hallmark'), showWarnings = F)
	output_filename = sprintf('genome_cluster_pathway_enrich_hallmark_cluster%d.html',number_of_clusters)

	rmarkdown::render(here('analysis/DESeq_cluster_genome_display/genome_cluster_pathway_enrich_hallmark.Rmd'),
										output_file = here('analysis/DESeq_cluster_genome_display/hallmark',output_filename),
										params = list(number_of_clusters = number_of_clusters))

	dir.create(here('analysis/DESeq_cluster_genome_display/KEGG'), showWarnings = F)
	output_filename = sprintf('genome_cluster_pathway_enrich_KEGG_cluster%d.html',number_of_clusters)

	rmarkdown::render(here('analysis/DESeq_cluster_genome_display/genome_cluster_pathway_enrich_KEGG.Rmd'),
										output_file = here('analysis/DESeq_cluster_genome_display/KEGG',output_filename),
										params = list(number_of_clusters = number_of_clusters))

	dir.create(here('analysis/DESeq_cluster_genome_display/TF'), showWarnings = F)
	output_filename = sprintf('genome_cluster_pathway_enrich_TFT_cluster%d.html',number_of_clusters)

	rmarkdown::render(here('analysis/DESeq_cluster_genome_display/genome_cluster_pathway_enrich_TFT.Rmd'),
										output_file = here('analysis/DESeq_cluster_genome_display/TF',output_filename),
										params = list(number_of_clusters = number_of_clusters))
}
