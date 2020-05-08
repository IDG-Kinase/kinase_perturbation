Genome Clustering
================
Matthew Berginski
4/7/2020

# Read in DESeq2 Normalized Log Changes

## Data Organizing/Filtering

## Experiment Correlations

``` r
tic()
genome_correlations_experiments = amap::Dist(drug_perturb_ratio_wide_experiments %>%
                                                                                            column_to_rownames(var = "treatment"),
                                                                                         method="correlation")
toc();
```

    ## 0.407 sec elapsed

``` r
correlation_vals = data.frame(
    cor_vals = as.vector((genome_correlations_experiments - 1)*-1)
)

ggplot(correlation_vals, aes(x=cor_vals)) + 
    geom_histogram() +
    xlab('Pairwise Treatment Correlations') +
    ylab('') +
    geom_vline(aes(xintercept = mean(cor_vals)), color='red') +
    theme_berginski()
```

    ## `stat_bin()` using `bins = 30`. Pick better value with `binwidth`.

![](DESeq2_genome_clustering_files/figure-gfm/correlation%20histogram-1.png)<!-- -->

``` r
ggsave(here('figures/genome_clustering/experiment_correlation_vals.png'),height=4,width=4)
```

    ## `stat_bin()` using `bins = 30`. Pick better value with `binwidth`.

``` r
ggplot(segment(ddata_x)) +
    geom_segment(aes(x=x, y=y, xend=xend, yend=yend)) +
    geom_text(data = dendro_labels,
                        aes(label = label, x=x,y=0), angle=90, hjust=1.1) +
    xlab('') +
    theme_berginski() +
    theme(axis.title.x=element_blank(),
                axis.text.x=element_blank(),
                axis.ticks.x=element_blank()) +
    # scale_y_discrete("Correlation Distance", waiver(), waiver(), c(0,2)) +
    scale_y_continuous("Genome Correlation Distance", breaks = c(0,0.5,1), limits = c(-0.8,NA)) +
    NULL
```

![](DESeq2_genome_clustering_files/figure-gfm/experiment%20cluster%20plot-1.png)<!-- -->

``` r
ggsave(here('figures/genome_clustering/experiment_cluster_tree.png'),height=4,width=4)
trimImage(here('figures/genome_clustering/experiment_cluster_tree.png'))
```

## Gene Correlations

``` r
tic()
genome_correlations_ratio = amap::Dist(drug_perturb_ratio_wide_genes %>%
                                                                                column_to_rownames(var = "hgnc_symbol"),
                                                                             method="correlation")
toc();
```

    ## 39.786 sec elapsed

``` r
genome_cluster_ratio = hclust(genome_correlations_ratio)

dd.row <- as.dendrogram(genome_cluster_ratio)
ddata_x <- dendro_data(dd.row)

dendro_labels = label(ddata_x) %>%
    left_join(all_kinases %>% select(symbol,class), by=c('label'='symbol')) %>%
    mutate(class = case_when(
        is.na(class) ~ "Non-kinase",
        TRUE ~ class
    ))
```

    ## Warning: Column `label`/`symbol` joining factors with different levels, coercing
    ## to character vector

``` r
end_points = ddata_x$segments %>% filter(yend == 0) %>%
    left_join(dendro_labels, by=c('x'='x'))

ggplot(segment(ddata_x)) +
    geom_segment(aes(x=x, y=y, xend=xend, yend=yend)) +
    # geom_text(data = dendro_labels,
    #                   aes(label = label, x=x,y=0, color = class), angle=90, hjust=1) +
    # ylim(c(-0.5,2)) +
    # geom_hline(aes(yintercept = 1.2), color='red') +
    xlab('') +
    # geom_point(data = end_points, mapping = aes(x=xend,y=yend-0.02, color=class),shape='\u007C') +
    geom_rect(data = end_points, mapping = aes(xmin=xend-0.5,xmax=xend+0.5,ymin=yend-0.25,ymax=yend-0.01,fill=class)) +
    theme_berginski() +
    theme(axis.title.x=element_blank(),
                axis.text.x=element_blank(),
                axis.ticks.x=element_blank()) +
    scale_y_continuous(name = "Correlation Distance",breaks=seq(0,2,by=0.5)) +
    scale_x_continuous(name = "",expand=c(0,0)) +
    scale_fill_manual(name = "Kinase Type",
                                        values = c("Dark" = "red", "Light" = "blue", "Non-kinase" = rgb(0,0,0,0)))
```

![](DESeq2_genome_clustering_files/figure-gfm/full%20genome%20tree-1.png)<!-- -->

``` r
ggsave(here('figures/genome_clustering/genome_clustering_DESeq_correlation.png'),width=20,height=3*0.9,dpi = 300)
trimImage(here('figures/genome_clustering/genome_clustering_DESeq_correlation.png'))
```

## Associated Genes

``` r
# genome_hier_all_clusters = cutree(genome_cluster_ratio,k = 2:dim(drug_perturb_ratio_wide_genes)[1]) %>%
#   as.data.frame() %>% 
#   rownames_to_column(var = "symbol") %>% 
#   pivot_longer(-symbol, 
#                            names_to = "num_total_clusters", 
#                            values_to = "cluster_id") %>% 
#   mutate(num_total_clusters = as.factor(num_total_clusters)) %>% 
#   left_join(all_kinases %>% select(class,symbol)) %>% 
#   arrange(num_total_clusters,symbol)
# 
# cluster_sizes = genome_hier_all_clusters %>% 
#   group_by(num_total_clusters,cluster_id) %>% 
#   count(name = "num_in_cluster") 
# 
# genome_hier_all_clusters = genome_hier_all_clusters %>% 
#   left_join(cluster_sizes)
```

``` r
# dir.create(here('figures/genome_clustering/associated_genes'), showWarnings = F)
# kinase_in_clusters = genome_hier_all_clusters %>% 
#   filter(symbol %in% all_kinases$symbol) %>%
#   pull(symbol) %>%
#   unique()
# 
# for (this_gene in kinase_in_clusters) {
#   gene_set = get_associated_genes(cluster_sets = genome_hier_all_clusters, 
#                                                                   gene_name = this_gene, 
#                                                                   max_clustered_genes = 20)
#   
#   gene_table = build_gene_associate_table_text(gene_set)
#   p<-tableGrob(gene_table %>%
#                               rename("Associated Genes" = "associated_genes"),
#                            rows = NULL)
#   grid.arrange(p)
#   ggsave(here('figures/genome_clustering/associated_genes',paste0(this_gene,'.png')),p)
#   trimImage(here('figures/genome_clustering/associated_genes',paste0(this_gene,'.png')))
# }
```
