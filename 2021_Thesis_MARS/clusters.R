# Theme loading
# chrdir = TRUE to allow the included source to work
source("lib/preambule.R", chdir = TRUE)

# Out directory if it does not exist
dir.create("out/clusters/", showWarnings = FALSE, recursive = TRUE)

genes_in_cluster <-
    read_tsv("data/clusters/genes_in_cluster.tsv")

cluster_with_lncRNA <-
    genes_in_cluster %>%
    filter(gene_type != "protein_coding") %>%
    pull(cluster_name) %>%
    unique

# The marneral cluster ID based on the presence of MRN1 gene
marneral_cluster <- 
    genes_in_cluster %>%
    filter(gene_id == "AT5G42600") %>%
    pull(cluster_name)

# how many lncRNAs in the clusters ?

graph <-
    genes_in_cluster %>%
    mutate(coding = if_else(gene_type == "protein_coding",
                            "coding", "noncoding")) %>%
    count(cluster_name, cluster_type, coding) %>%
    pivot_wider(names_from="coding", values_from="n") %>%
    mutate(noncoding=ifelse(is.na(noncoding), 0, noncoding)) %>%
    ggplot(aes(x=noncoding, fill=cluster_type)) +
    geom_bar(position = "dodge") +
    scale_fill_brewer(palette="Set1", "cluster_type") +
    xlab("number of lncRNA in the cluster") + ylab("number of cluster")
ggsave("./out/clusters/nb_lncRNAs_in_cluster.pdf", plot=graph, width=7, height=5)

# doeas the cluster have a lncRNA ?

graph <-
    genes_in_cluster %>%
    mutate(coding = if_else(gene_type == "protein_coding",
                            "coding", "noncoding")) %>%
    count(cluster_name, cluster_type, coding) %>%
    pivot_wider(names_from="coding", values_from="n") %>%
    # if NA it is because while counting we did not find any lncRNAs
    mutate(noncoding=!is.na(noncoding)) %>%
    ggplot(aes(x=cluster_type, fill=noncoding)) +
    geom_bar(position = "dodge") +
    scale_fill_brewer(palette="Set1", "lncRNA in cluster") +
    xlab("") + ylab("number of cluster")
ggsave("./out/clusters/lncRNAs_in_cluster.pdf", plot=graph, width=5, height=3.5)

## Correlation of the cluster according to the presence of lncRNAs
# all correlations

# load the correlations
corr_info <-
    read_tsv("data/clusters/corr_info.tsv")

graph <-
    corr_info %>%
    dplyr::select(cluster_type, cluster_name,
                  cluster_max_corr, cluster_median_corr) %>%
    unique %>%
    filter(!is.na(cluster_median_corr)) %>%
    mutate(lncRNA = cluster_name %in% cluster_with_lncRNA,
           cluster_shortname = if_else(cluster_name %in% marneral_cluster, "marneral", NA_character_)) %>%
    ggplot(aes(x = lncRNA, y = cluster_median_corr)) +
    geom_violin(fill = "gray", alpha = .2) +
    geom_boxplot(width=.1, outlier.shape = NA) +
    geom_jitter(height = 0, width = .1, size = 1.5,
                aes(color=cluster_shortname, shape = cluster_shortname)) +
    facet_wrap(~cluster_type) +
    xlab("lncRNA in the cluster") +
    ylab("in between cluster coding genes correlation") +
    scale_colour_manual(values = c("red"), na.translate = TRUE, na.value = "black",  guide = 'none') +
    scale_shape_manual(values = c(19), na.translate = TRUE, na.value = 1,  guide = 'none') +
    ggtitle("Correlation median")
ggsave("out/clusters/cluster_median_corr.pdf", graph, width=5, height=6)

graph <-
    corr_info %>%
    dplyr::select(cluster_type, cluster_name, cluster_max_corr) %>%
    unique %>%
    filter(!is.infinite(cluster_max_corr)) %>%
    mutate(lncRNA = cluster_name %in% cluster_with_lncRNA,
           cluster_shortname = if_else(cluster_name %in% marneral_cluster, "marneral", NA_character_)) %>%
    ggplot(aes(x = lncRNA, y = cluster_max_corr)) +
    geom_violin(fill = "gray", alpha = .2) +
    geom_boxplot(width=.1, outlier.shape = NA) +
    geom_jitter(height = 0, width = .1, size = 1.5,
                aes(color=cluster_shortname, shape = cluster_shortname)) +
    facet_wrap(~cluster_type) +
    xlab("lncRNA in the cluster") +
    ylab("in between cluster coding genes correlation") +
    scale_colour_manual(values = c("red"), na.translate = TRUE, na.value = "black",  guide = 'none') +
    scale_shape_manual(values = c(19), na.translate = TRUE, na.value = 1,  guide = 'none') +
    ggtitle("Correlation maximum")
ggsave("out/clusters/cluster_max_corr.pdf", graph, width=5, height=6)


# only significant correlations
# load the correlations
sig_corr_info <-
    read_tsv("data/clusters/sig_corr_info.tsv")

graph <-
    sig_corr_info %>%
    dplyr::select(cluster_type, cluster_name, cluster_max_corr) %>%
    unique %>%
    filter(!is.infinite(cluster_max_corr)) %>%
    mutate(lncRNA = cluster_name %in% cluster_with_lncRNA,
           cluster_shortname = if_else(cluster_name %in% marneral_cluster, "marneral", NA_character_)) %>%
    ggplot(aes(x = lncRNA, y = cluster_max_corr)) +
    geom_violin(fill = "gray", alpha = .2) +
    geom_boxplot(width=.1, outlier.shape = NA) +
    geom_jitter(height = 0, width = .1, size = 1.5,
                aes(color=cluster_shortname, shape = cluster_shortname)) +
    facet_wrap(~cluster_type) +
    xlab("lncRNA in the cluster") +
    ylab("in between cluster coding genes correlation") +
    scale_colour_manual(values = c("red"), na.translate = TRUE, na.value = "black",  guide = 'none') +
    scale_shape_manual(values = c(19), na.translate = TRUE, na.value = 1,  guide = 'none') +
    ggtitle("Correlation maximum (only significant ones)")
ggsave("out/clusters/cluster_max_sig_corr.pdf", graph, width=5, height=6)

graph <-
    sig_corr_info %>%
    dplyr::select(cluster_type, cluster_name, cluster_median_corr) %>%
    unique %>%
    filter(!is.na(cluster_median_corr)) %>%
    mutate(lncRNA = cluster_name %in% cluster_with_lncRNA,
           cluster_shortname = if_else(cluster_name %in% marneral_cluster, "marneral", NA_character_)) %>%
    ggplot(aes(x = lncRNA, y = cluster_median_corr)) +
    geom_violin(fill = "gray", alpha = .2) +
    geom_boxplot(width=.1, outlier.shape = NA) +
    geom_jitter(height = 0, width = .1, size = 1.5,
                aes(color=cluster_shortname, shape = cluster_shortname)) +
    facet_wrap(~cluster_type) +
    xlab("lncRNA in the cluster") +
    ylab("in between cluster coding genes correlation") +
    scale_colour_manual(values = c("red"), na.translate = TRUE, na.value = "black",  guide = 'none') +
    scale_shape_manual(values = c(19), na.translate = TRUE, na.value = 1,  guide = 'none') +
    ggtitle("Correlation median (only significant ones)")
ggsave("out/clusters/cluster_median_sig_corr.pdf", graph, width=5, height=6)


# expression correlation of lncRNA and the one in the clusters
graph <-
    sig_corr_info %>%
    filter(!is.na(query)) %>%
    ggplot(aes(x = cluster_max_corr, y = max_corr)) +
    geom_point() +
    facet_wrap(~cluster_type) +
    xlab("in between cluster coding genes correlation") +
    ylab("lncRNA with coding gene correlation") +
    ggtitle("Correlation maximum (only significant ones)")
ggsave("out/clusters/sig_corr_max.pdf", graph, width=8, height=5)

graph <-
    sig_corr_info %>%
    filter(!is.na(query)) %>%
    ggplot(aes(x = cluster_median_corr, y = median_corr)) +
    geom_linerange(aes(xmin = cluster_median_corr - cluster_sd_corr,
                       xmax = cluster_median_corr + cluster_sd_corr )) +
    geom_linerange(aes(ymin = median_corr - sd_corr,
                       ymax = median_corr + sd_corr )) +
    geom_point() +
    facet_wrap(~cluster_type) +
    xlab("in between cluster coding genes correlation") +
    ylab("lncRNA with coding gene correlation") +
    ggtitle("Correlation median (only significant ones)")
ggsave("out/clusters/sig_corr_median.pdf", graph, width=8, height=5)

graph <-
    corr_info %>%
    dplyr::select(cluster_type, cluster_name, cluster_median_corr) %>%
    unique %>%
    filter(!is.na(cluster_median_corr)) %>%
    mutate(lncRNA = cluster_name %in% cluster_with_lncRNA,
           sig_corr = cluster_name %in% sig_corr_info$cluster_name) %>%
    ggplot(aes(x=lncRNA, fill=sig_corr)) +
    geom_bar() +
    facet_wrap(~cluster_type) +
    scale_fill_brewer(palette= "Accent", name = "significant correlation") +
    xlab("lncRNA in the cluster") + 
    ylab("number of clusters")
ggsave("out/clusters/lncRNA_sig_clusters.pdf",
       plot = graph, width=4, height=2.25)


sig_corr <- corr_info %>%
    dplyr::select(cluster_type, cluster_name, cluster_median_corr) %>%
    unique %>%
    filter(!is.na(cluster_median_corr)) %>%
    mutate(lncRNA = cluster_name %in% cluster_with_lncRNA,
           sig_corr = cluster_name %in% sig_corr_info$cluster_name) 

write.csv(sig_corr, file="sig_corr")

graph <-
    corr_info %>%
    dplyr::select(cluster_type, cluster_name, cluster_median_corr) %>%
    unique %>%
    filter(!is.na(cluster_median_corr)) %>%
    mutate(lncRNA = cluster_name %in% cluster_with_lncRNA,
           sig_corr = cluster_name %in% sig_corr_info$cluster_name) %>%
    ggplot(aes(x=lncRNA, y=cluster_median_corr, color=sig_corr, fill = sig_corr)) +
    geom_violin(position = position_dodge(.9), alpha = .2, scale="count") +
    facet_wrap(~cluster_type) +
    scale_fill_brewer(palette= "Accent", name = "significant correlation") +
    scale_color_brewer(palette= "Accent", name = "significant correlation") +
    xlab("lncRNA in the cluster") + 
    ylab("median correlation")
ggsave("out/clusters/violin_median_correlation.pdf",
       plot = graph, width = 6, height = 3)


graph <-
    corr_info %>%
    dplyr::select(cluster_type, cluster_name, cluster_median_corr) %>%
    unique %>%
    filter(!is.na(cluster_median_corr)) %>%
    mutate(lncRNA = cluster_name %in% cluster_with_lncRNA,
           sig_corr = cluster_name %in% sig_corr_info$cluster_name) %>%
    ggplot(aes(x=lncRNA, y=cluster_median_corr, color=sig_corr, fill = sig_corr)) +
    geom_boxplot(position = position_dodge(.9, preserve = "single"), alpha = .2, guide = FALSE) +
    facet_wrap(~cluster_type) +
    scale_fill_brewer(palette= "Accent", name = "significant correlation") +
    scale_color_brewer(palette= "Accent", name = "significant correlation") +
    xlab("lncRNA in the cluster") + 
    ylab("median correlation")
ggsave("out/clusters/boxplot_median_correlation.pdf",
       plot = graph, width = 6, height = 3)


lncRNA_corr <- read_tsv("data/clusters/lncRNA_max_corr.tsv")



graph <-
    ggplot(lncRNA_corr, aes(x=cluster_type, y=corr)) +
    geom_violin() +
    geom_boxplot(width = .2) +
    scale_fill_brewer(palette= "Accent", name = "correlation type") +
    scale_color_brewer(palette= "Accent", name = "correlation type") +
    xlab("") + 
    ylab("maximum lncRNA correlation")
ggsave("out/clusters/violin_max_lncRNA_correlation.pdf",
       plot = graph, width=4, height=3.5)



#stat

lncRNA_corr %>%
    ggbarplot(., x = "cluster_type", y = "corr",add = "mean_se") +
    stat_compare_means(method="t.test", 
                       comparisons = list(c("coexp","plantiSMASH") ), 
                       label = "p.label", bracket.size= 0.5) +
    theme_bw()


