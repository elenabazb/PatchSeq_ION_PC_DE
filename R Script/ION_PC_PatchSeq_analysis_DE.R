library(dplyr)
library(clipr)
library(tidyverse)
library(tidyr)
library(data.table)
library(ggplot2)
library(gghalves)
library(ggpubr)
library(ggrepel)
library(ggrastr)
library(rstatix)
library(Seurat)
library(DESeq2)
library(pheatmap)
library(viridis)
library(ragg)
library(Cairo)


# Provided count file should have gene names already translated ####

# Session info ####
# sessioninfo <- sessionInfo()
# writeLines(capture.output(sessionInfo()), "sessionInfo.txt")

# 1) Load raw data and metadata ####
## ION ###
ionraw <- readRDS("ionraw_counts.RDS")
metadata_ion <- readRDS("metadata_ion_export.RDS")

## PC ###
pcraw <- readRDS("pcraw_counts.RDS")
metadata_pc <- readRDS("metadata_pc_export.RDS")





# 2) QC ####
## ION ####
ionraw <- ionraw[,metadata_ion$Library]

# Total counts
metadata_ion$n_counts <- colSums(apply(ionraw, 2, as.numeric))

# Total genes
for (i in 1:ncol(ionraw)){
  metadata_ion[i,"n_genes"] = table(ionraw[,i] > 0)[2] %>% as.numeric()
}
metadata_ion["n_genes"][is.na(metadata_ion["n_genes"])] <- 0


## PC ####
pcraw <- pcraw[,metadata_pc$Library]

# Total counts
metadata_pc$n_counts <- colSums(apply(pcraw, 2, as.numeric))

# Total genes
for (i in 1:ncol(pcraw)){
  metadata_pc[i,"n_genes"] = table(pcraw[,i] > 0)[2] %>% as.numeric()
}
metadata_pc["n_genes"][is.na(metadata_pc["n_genes"])] <- 0


## Merge metadata ####
metadata_merged <- rbind(metadata_ion, metadata_pc)

# create mixed variable age/celltype
metadata_merged <- metadata_merged %>% mutate(celltype_age = interaction(Age_dis, celltype, sep="_"))

colors_celltype <- c("Young_ION"="darkorange","Old_ION"="darkred", "Young_PC"="skyblue", "Old_PC"="darkblue")


## Scatter plot Nfeatures / Ncounts ####
ggplot(metadata_merged, aes(x = n_genes, y = log10(n_counts), color = celltype_age)) +
  geom_point(size = 3, alpha=.7) +
  scale_color_manual(values = colors_celltype) +
  labs(x = "Number of detected genes", y = "Total counts (log10)", title = "") +
  geom_hline(yintercept = 4, color="red", linetype =2)+
  geom_vline(xintercept = 10000, color="red", linetype=2)+
  theme_classic2(base_size = 35)


## Filter ####
metadata_merged <- metadata_merged %>% filter(n_counts > 10000, n_genes < 10000)

# ION counts
ionraw_filt <- ionraw[,colnames(ionraw)%in%metadata_merged$Library]

# PC counts 
pcraw_filt <- pcraw[,colnames(pcraw)%in%metadata_merged$Library]




# 3) Log-Normalize ####
rownames(metadata_merged) <- metadata_merged$Library
metadata_ion_merge <- metadata_merged %>% filter(celltype == "ION")
metadata_pc_merge <- metadata_merged %>% filter(celltype == "PC")


## Create Seurat object IONs and PCs####
ionseurat <- CreateSeuratObject(counts = ionraw_filt, min.cells = 0, min.features = 0, meta.data = as.data.frame(metadata_ion_merge))
ionseurat <- SetIdent(ionseurat, value=ionseurat@meta.data$Age_dis)

pcseurat <- CreateSeuratObject(counts = pcraw_filt, min.cells = 0, min.features = 0, meta.data = as.data.frame(metadata_pc_merge))
pcseurat <- SetIdent(pcseurat, value=pcseurat@meta.data$Age_dis)


## Log normalization and save ####
ionseurat <- NormalizeData(ionseurat, normalization.method = "LogNormalize", scale.factor = 1e6)
ionnorm <- GetAssayData(ionseurat, layer='data') %>% as.data.frame(.)
saveRDS(ionnorm,"ion_lognorm.RDS")

pcseurat <- NormalizeData(pcseurat, normalization.method = "LogNormalize", scale.factor = 1e6)
pcnorm <- GetAssayData(pcseurat, layer='data') %>% as.data.frame(.)
saveRDS(pcnorm,"pc_lognorm.RDS")



# 4) Marker genes expression ####
## Define marker genes ####
ionmarkers <- c("Calb1", "Pou4f1", "Foxp2", "Foxp1", "Slc17a6", "Crh", "Gap43", "Robo3", "Dcc", "Reln",
                "Grm5", "Olig3", "Igsf9", "Cadps2", "Pax6", "Rph3a", "C1ql1", "Nrcam", "Lgi2", "Sema4f")
pcmarkers <- c("Calb1", "Slc1a2", "Gabbr2", "Hspb1", "Ppp1r17", "Itpr1", "Prkcg", "Pcp2", "B3gat1", "Nptn", "Grm1", "Slc1a6",
               "Aldoc", "Sphk1", "Ebf2", "Pcp4", "Gabbr1", "Plcb3")
markers_pc_ion <- unique(c(ionmarkers, pcmarkers))


## Subset counts ####
ionmarkers_cts <- ionnorm[rownames(ionnorm) %in% markers_pc_ion,]
ionmarkers_cts <- as.matrix(ionmarkers_cts)
mode(ionmarkers_cts) <- "numeric"

pcmarkers_cts <- pcnorm[rownames(pcnorm) %in% markers_pc_ion,]
pcmarkers_cts <- as.matrix(pcmarkers_cts)
mode(pcmarkers_cts) <- "numeric"



## Dotplot ####
### Define dotplot function ####
make_dotplot_df_2 <- function(mat, metadata, celltype, age_col = "Age_dis", library_col = "Library") {
  
  tmp <- as.matrix(mat)
  mode(tmp) <- "numeric"
  metadata <- as.data.frame(metadata)
  
  meta_tmp <- metadata %>% select(Library = all_of(library_col), Age_dis = all_of(age_col)) %>%
    filter(Library %in% colnames(tmp)) %>% arrange(match(Library, colnames(tmp)))
  
  if (!all(meta_tmp$Library == colnames(tmp))) {
    stop("Metadata and matrix columns are not in the same order.")}
  
  tmp %>% as.data.frame() %>% tibble::rownames_to_column("gene") %>%
    tidyr::pivot_longer(cols = -gene, names_to = "Library", values_to = "expr") %>%
    left_join(meta_tmp, by = "Library") %>% group_by(gene, Age_dis) %>%
    summarise(pct_detected = mean(expr > 0, na.rm = TRUE) * 100, avg_expr = mean(expr, na.rm = TRUE), .groups = "drop") %>%
   mutate(Celltype = celltype)
}


### Create a df per cell type ####
dot_ion <- make_dotplot_df_2(mat = ionmarkers_cts, metadata = metadata_ion_merge, celltype = "ION",
                                  age_col = "Age_dis", library_col = "Library")

dot_pc <- make_dotplot_df_2(mat = pcmarkers_cts, metadata = metadata_pc_merge, celltype = "PC",
                                 age_col = "Age_dis", library_col = "Library")

# Bind
dot_ion_pc <- bind_rows(dot_pc, dot_ion)


### Compute z-score ####
dot_ion_pc <- dot_ion_pc %>% group_by(gene) %>% mutate(avg_expr_z = as.numeric(scale(avg_expr))) %>% ungroup()


### Order genes by average expression ####
gene_order_df <- dot_ion_pc %>% group_by(gene, Celltype) %>%
  summarise(mean_expr_celltype = mean(avg_expr, na.rm = TRUE), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = Celltype, values_from = mean_expr_celltype, values_fill = 0) %>%
  mutate(PC_vs_ION_score = PC - ION) %>% arrange(desc(PC_vs_ION_score))

gene_order <- gene_order_df$gene

dot_ion_pc <- dot_ion_pc %>% mutate(gene = factor(gene, levels = gene_order), Age_dis = factor(Age_dis, levels = c("Young", "Old")),
                                    Celltype = factor(Celltype, levels = c("PC", "ION")))


### Plot ####
ggplot(dot_ion_pc, aes(x = gene, y = Age_dis)) +
  geom_point(aes(size = pct_detected, color = avg_expr_z)) +
  facet_grid(Celltype ~ ., scales = "free_y", space = "free_y") +
  scale_size_continuous(range = c(0, 10), limits = c(0, 100)) +
  scale_color_gradientn(colors = c("grey25", "grey70","white","#3BC1A8","#061E29"),
                        values = scales::rescale(c(min(dot_ion_pc$avg_expr_z, na.rm = TRUE), -1,0,0.8,1.6, max(dot_ion_pc$avg_expr_z, na.rm = TRUE))),na.value = "grey90")+
  labs(x = NULL, y = NULL, size = "% cells",color = "Mean expr.\nZ-score") +
  theme_classic2(base_size = 27) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), strip.background = element_blank(), strip.text.y = element_text(angle = 0))



## Calb1 expression ####
### Subset ION counts ####
calb1_counts_ion <- ionnorm[rownames(ionnorm) %in% "Calb1",]
calb1_counts_ion <- calb1_counts_ion[,metadata_ion_merge$Library]

# make long format
calb1_counts_ion_long <- calb1_counts_ion %>%as.data.frame() %>%rownames_to_column("gene") %>%
  pivot_longer(cols = -gene, names_to = "Library", values_to = "expr")

# merge with metadata
calb1_counts_ion_long <- calb1_counts_ion_long %>% left_join(metadata_ion_merge, by = "Library")
calb1_counts_ion_long$Age_dis <- factor(calb1_counts_ion_long$Age_dis, levels = c("Young", "Old"))


### PC counts ####
# subset normalized counts
calb1_counts_pc <- pcnorm[rownames(pcnorm) %in% "Calb1",]
calb1_counts_pc <- calb1_counts_pc[,metadata_pc_merge$Library]

# make long format
calb1_counts_pc_long <- calb1_counts_pc %>%as.data.frame() %>%rownames_to_column("gene") %>%
  pivot_longer(cols = -gene, names_to = "Library", values_to = "expr")

# merge with metadata
calb1_counts_pc_long <- calb1_counts_pc_long %>% left_join(metadata_pc_merge, by = "Library")
calb1_counts_pc_long$Age_dis <- factor(calb1_counts_pc_long$Age_dis, levels = c("Young", "Old"))


### Merge ####
calb1_counts_both <- bind_rows(calb1_counts_pc_long, calb1_counts_ion_long)

# add combined variable
calb1_counts_both <- calb1_counts_both %>% left_join(metadata_merged %>% select(Library, celltype_age), by = "Library") %>%
  mutate(celltype_age = factor(celltype_age, levels = c("Young_ION", "Old_ION", "Young_PC", "Old_PC")))


### Plot ####
ggplot(as.data.frame(calb1_counts_both), aes(x = celltype_age, y = expr, fill = celltype_age)) +
  geom_half_violin(trim = FALSE, alpha = 0.5, color = "transparent") +
  geom_half_boxplot(width = 0.2, alpha = 0.6, outlier.shape = NA, side="r") +
  geom_jitter(width = 0.08, size = 3, alpha = 0.7) +
  scale_fill_manual(values = c("Young_PC" = "#4f97c7", "Old_PC" = "darkblue","Young_ION" = "darkorange", "Old_ION" = "darkred")) +
  coord_flip()+
  labs(x = NULL, y = "Norm. expression", title = "Calb1") +
  theme_classic2(base_size = 45)+
  theme(legend.position = "none", plot.title = element_text(face = "italic", hjust = .5))



# 5) Differential expression IONs ####
## LRT ####
library(DESeq2)
dds_ion <- DESeqDataSetFromMatrix(as.matrix(ionraw_filt), colData = metadata_ion_merge, design = ~Age_dis)
dds_ion <- DESeq(dds_ion, test = "LRT", useT = TRUE, minmu = 1e-6, minReplicatesForReplace = Inf, reduced = ~1)
res_lrt_ion <- results(dds_ion)
sig_res_lrt_ion <- res_lrt_ion %>% data.frame() %>% rownames_to_column(var="gene") %>% as_tibble() %>% filter(padj < 0.05)


## Extract up and downregulated for EnrichR ####
genes_up_age_ion <- sig_res_lrt_ion %>% filter(log2FoldChange > 0) %>% pull(gene)
genes_down_age_ion <- sig_res_lrt_ion %>% filter(log2FoldChange < 0) %>% pull(gene)

#remove pseudogenes
genes_up_age_ion <- genes_up_age_ion[!grepl("^Gm", genes_up_age_ion)]
genes_down_age_ion <- genes_down_age_ion[!grepl("^Gm", genes_down_age_ion)]


# The vectors genes_up_age_ion and genes_down_age_ion were fed to EnrichR to perform 
# pathway analysis online (https://maayanlab.cloud/Enrichr/)



## EnrichR results with significantly enriched pathways - upregulated genes IONs ####
### Biological process (UP) ####
#import enrichr results 
ionageup_enrichr_BP <- read.delim("ION_age_UP_GO_Biological_Process_2026_table.txt")
ionageup_enrichr_BP <- ionageup_enrichr_BP[ionageup_enrichr_BP$Adjusted.P.value<0.1,]
ionageup_enrichr_BP <- ionageup_enrichr_BP %>% add_column(Group=1)
ionageup_enrichr_BP$GenesList <- strsplit(ionageup_enrichr_BP$Genes, split = ";\\s*")

# add columns with log-adj. p-value and genes lists
ionageup_enrichr_BP <- ionageup_enrichr_BP %>% mutate(nGenes = sapply(GenesList, length), neglog10_p = -log10(Adjusted.P.value))

ionageup_enrichr_BP$Term <- str_wrap(ionageup_enrichr_BP$Term, width = 25)

# order by significance
ionageup_enrichr_BP <- ionageup_enrichr_BP %>% mutate(Term = reorder(Term, neglog10_p))

#plot
ggplot(ionageup_enrichr_BP, aes(x = neglog10_p, y = Term, size = nGenes)) +
  geom_point(show.legend = T, alpha=.8) +
  labs(y="", x = "-log(adj. p)", size = "Number of Genes", title = "Upregulated in Young IONs (BP)") +
  theme_classic2(base_size = 20) +
  theme(axis.text.y = element_text(size = 22), axis.text.x = element_text(size = 18))



## EnrichR results with significantly enriched pathways - downregulated genes IONs ####
### Biological process (DOWN) ####
#import enrichr results 
ionagedown_enrichr_BP <- read.delim("ION_age_down_GO_Biological_Process_2026_table.txt")
ionagedown_enrichr_BP <- ionagedown_enrichr_BP[ionagedown_enrichr_BP$Adjusted.P.value<0.1,]
ionagedown_enrichr_BP <- ionagedown_enrichr_BP %>% add_column(Grodown=1)
ionagedown_enrichr_BP$GenesList <- strsplit(ionagedown_enrichr_BP$Genes, split = ";\\s*")

# add columns with log-adj. p-value and genes lists
ionagedown_enrichr_BP <- ionagedown_enrichr_BP %>% mutate(nGenes = sapply(GenesList, length), neglog10_p = -log10(Adjusted.P.value))

ionagedown_enrichr_BP$Term <- str_wrap(ionagedown_enrichr_BP$Term, width = 25)

# order by significance
ionagedown_enrichr_BP <- ionagedown_enrichr_BP %>% mutate(Term = reorder(Term, neglog10_p))

#plot
ggplot(ionagedown_enrichr_BP, aes(x = neglog10_p, y = Term, size = nGenes)) +
  geom_point(show.legend = T, alpha=.8) +
  labs(y="", x = "-log(adj. p)", size = "Number of Genes", title = "downregulated in Young IONs (BP)") +
  theme_classic2(base_size = 20) +
  theme(axis.text.y = element_text(size = 22), axis.text.x = element_text(size = 18))


### Cellular component (DOWN) ####
#import enrichr results 
ionagedown_enrichr_CC <- read.delim("ION_age_down_GO_Cellular_Component_2026_table.txt")
ionagedown_enrichr_CC <- ionagedown_enrichr_CC[ionagedown_enrichr_CC$Adjusted.P.value<0.1,]
ionagedown_enrichr_CC <- ionagedown_enrichr_CC %>% add_column(Grodown=1)
ionagedown_enrichr_CC$GenesList <- strsplit(ionagedown_enrichr_CC$Genes, split = ";\\s*")

# add columns with log-adj. p-value and genes lists
ionagedown_enrichr_CC <- ionagedown_enrichr_CC %>% mutate(nGenes = sapply(GenesList, length), neglog10_p = -log10(Adjusted.P.value))

ionagedown_enrichr_CC$Term <- str_wrap(ionagedown_enrichr_CC$Term, width = 25)

# order by significance
ionagedown_enrichr_CC <- ionagedown_enrichr_CC %>% mutate(Term = reorder(Term, neglog10_p))

#plot
ggplot(ionagedown_enrichr_CC, aes(x = neglog10_p, y = Term, size = nGenes)) +
  geom_point(show.legend = T, alpha=.8) +
  labs(y="", x = "-log(adj. p)", size = "Number of Genes", title = "downregulated in Young IONs (CC)") +
  theme_classic2(base_size = 20) +
  theme(axis.text.y = element_text(size = 22), axis.text.x = element_text(size = 18))



## Extract genes from pathways ####
# Define function to extract genes from enriched pathways
prepare_pathway_gene_counts <- function(enrichr_df, norm_counts, metadata) {
  genes <- enrichr_df$GenesList %>% unlist() %>% str_to_title(tolower()) %>% unique()
  
  counts_sub <- norm_counts[rownames(norm_counts) %in% genes, , drop = FALSE]
  counts_sub <- counts_sub[, colnames(counts_sub) %in% metadata$Library, drop = FALSE]
  
  counts_long <- counts_sub %>% as.data.frame() %>% rownames_to_column("gene") %>%
    pivot_longer(cols = -gene, names_to = "Library", values_to = "expr") %>%
    left_join(metadata, by = "Library") %>%
    mutate(Age_dis = factor(Age_dis, levels = c("Young", "Old")))
  
  list(genes = genes, counts_long = counts_long)
}


# Define function to plot violins for genes
plot_pathway_violins <- function(counts_long, suffix, outdir = ".", colors = c("Young" = "darkorange", "Old" = "darkred")) {
  genes_to_plot <- unique(counts_long$gene)
  
  for (g in genes_to_plot) {
    df_gene <- counts_long %>%filter(gene == g)
    
    p <- ggplot(df_gene, aes(x = Age_dis, y = expr, fill = Age_dis)) +
      geom_half_violin(trim = FALSE, alpha = 0.5, color = "transparent") +
      geom_half_boxplot(width = 0.2, alpha = 0.6, outlier.shape = NA, side = "r") +
      geom_jitter(width = 0.08, size = 3, alpha = 0.7) +
      scale_fill_manual(values = colors) +
      labs(x = NULL, y = "Norm. expression", title = g) +
      theme_classic2(base_size = 45) +
      theme(legend.position = "none", plot.title = element_text(hjust = 0.5, face = "italic"))
    
    g_safe <- str_replace_all(g, "[^A-Za-z0-9_.-]", "_")
    
    ggsave(filename = file.path(outdir, paste0("Vln_", g_safe, "_", suffix, ".pdf")), plot = p, device = cairo_pdf, width = 8, height = 8)
  }
}


# Merge functions
process_pathway_violin_set <- function(enrichr_df, norm_counts, metadata, sig_res, suffix, outdir = ".",
                                       colors = c("Young" = "darkorange", "Old" = "darkred")) {
  prepared <- prepare_pathway_gene_counts(enrichr_df = enrichr_df, norm_counts = norm_counts, metadata = metadata)
  
  plot_pathway_violins(counts_long = prepared$counts_long, suffix = suffix, outdir = outdir, colors=colors)
  
  pvals <- sig_res %>% filter(gene %in% prepared$genes) %>% mutate(pathway_set = suffix)
  
  list(genes = prepared$genes, counts_long = prepared$counts_long, pvals = pvals
  )
}


### Genes from GO Biological Process (upregulated) ####
dir.create("Vln_genes_pathways_Age_ION", showWarnings = FALSE)

up_BP <- process_pathway_violin_set(enrichr_df = ionageup_enrichr_BP, norm_counts = ionnorm, 
                                    metadata = metadata_ion_merge, sig_res = sig_res_lrt_ion, 
                                    suffix = "upreg_age_BP",outdir = "Vln_genes_pathways_Age_ION")

### Genes from GO Biological Process (downregulated) ####
down_BP <- process_pathway_violin_set(enrichr_df = ionagedown_enrichr_BP, norm_counts = ionnorm, 
                                      metadata = metadata_ion_merge, sig_res = sig_res_lrt_ion,
                                      suffix = "downreg_age_BP", outdir = "Vln_genes_pathways_Age_ION")

### Genes from GO Cellular component (downregulated) ####
down_CC <- process_pathway_violin_set(enrichr_df = ionagedown_enrichr_CC, norm_counts = ionnorm,
                                      metadata = metadata_ion_merge, sig_res = sig_res_lrt_ion,
                                      suffix = "downreg_age_CC", outdir = "Vln_genes_pathways_Age_ION")


### P-values from LRT ####
pvalue_genes_pathways <- bind_rows(up_BP$pvals, down_BP$pvals, down_CC$pvals)




## Other GOIs from literature ####
goi_ion <- c("Gria3", "Atxn1", "Camk2b", "Foxp1", "Grin1", "Alkbh5", "Fto", "Nova2", "Adam22",
             "Grik4", "Bdnf", "Camk4", "Cacng2", "Cspg5", "Prkcg", "Camk2a", "Grm1", "Plcb4", "Cacna1a",
             "Adgrb3", "Arc", "Grn", "Sema3a", "Sema7a", "Htr3a", "Igf1", "Nrxn1", "Nrxn2", "Ntrk2", 
             "Plxna4", "Myo5a", "Camk2a", "Camk4", "Ptprd", "Fmr1", "C1ql1")

# subset normalized counts
goi_ion_counts <- ionnorm[rownames(ionnorm) %in% goi_ion,]
goi_ion_counts <- goi_ion_counts[,metadata_ion_merge$Library]

# make long format
goi_ion_counts_long <- goi_ion_counts %>%as.data.frame() %>% rownames_to_column("gene") %>%
  pivot_longer(cols = -gene, names_to = "Library", values_to = "expr")

# merge with metadata
goi_ion_counts_long <- goi_ion_counts_long %>% left_join(metadata_ion_merge, by = "Library")
goi_ion_counts_long$Age_dis <- factor(goi_ion_counts_long$Age_dis, levels = c("Young", "Old"))

# plot
genes_to_plot <- unique(goi_ion_counts_long$gene)
dir.create("Vln GOIs ION", showWarnings = FALSE)

for (g in genes_to_plot){
  df_gene <- goi_ion_counts_long %>% filter(gene == g)
  
  p <- ggplot(df_gene, aes(x = Age_dis, y = expr, fill = Age_dis)) +
    geom_half_violin(trim = FALSE, alpha = 0.5, color = "transparent") +
    geom_half_boxplot(width = 0.2, alpha = 0.6, outlier.shape = NA, side="r") +
    geom_jitter(width = 0.08, size = 3, alpha = 0.7) +
    scale_fill_manual(values = c("Young" = "darkorange", "Old" = "darkred")) +
    labs(x = NULL, y = "Norm. expression", title = g) +
    theme_classic2(base_size = 45)+
    theme(legend.position = "none", plot.title = element_text(hjust = .5, face = "italic"))
  
  ggsave(paste0("./Vln GOIs ION/Vln_GOI_ION_", g, "_.pdf"), plot = p, device = cairo_pdf, width = 8, height = 8)
}

pval_cf_ion <- sig_res_lrt_ion[sig_res_lrt_ion$gene %in% goi_ion_counts_long$gene,]




## Volcano plot DE IONs ####
#prepare data for volcano plot
volcano_df_ion <- res_lrt_ion %>% as.data.frame() %>% rownames_to_column("gene") %>% filter(!is.na(log2FoldChange), !is.na(padj)) %>%
  mutate(padj_plot = ifelse(padj == 0, min(padj[padj > 0], na.rm = TRUE), padj),
         neg_log10_padj = -log10(padj_plot), regulation = case_when(padj < 0.05 & log2FoldChange > 0.5  ~ "Up", 
                                                                    padj < 0.05 & log2FoldChange < -0.5 ~ "Down", TRUE ~ "NS"))
# select genes to highlight from pathways analysis
genes_to_label_ion <- c("Cntn2", "Agrn", "Sptbn4", "Nrcam", "Nlgn1", "Lrp4", "Prnp", "Cdkl5", "Pak1", "Grik2", "Atcay", "Ago2", "Prrt1", "Dagla", "Cspg5")

# plot
ggplot(volcano_df_ion, aes(x = log2FoldChange, y = neg_log10_padj)) +
  ggrastr::rasterise(geom_point(aes(color = regulation), alpha = 0.6, size = 2), dpi=300)+
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", color = "grey40") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40") +
  geom_text_repel(data = volcano_df_ion %>% filter(gene %in% genes_to_label_ion), aes(label = gene),
                  size = 6, box.padding = 0.8, point.padding = 0.6, force = 5, force_pull = 0.5, max.overlaps = Inf,
                  min.segment.length = 0, segment.color = "black", segment.size = .6) +
  scale_color_manual(values = c("Up" = "darkorange", "Down" = "darkred", "NS" = "grey70")) +
  labs(x = "log2 fold change", y = "-log10 adjusted p-value", color = NULL) +
  theme_classic2(base_size = 30)




# 6) Differential expression PCs ####
## LRT ####
dds_pc <- DESeqDataSetFromMatrix(as.matrix(pcraw_filt), colData = metadata_pc_merge, design = ~Age_dis)
dds_pc <- DESeq(dds_pc, test = "LRT", useT = TRUE, minmu = 1e-6, minReplicatesForReplace = Inf, reduced = ~1)
res_lrt_pc <- results(dds_pc)
sig_res_lrt_pc <- res_lrt_pc %>% data.frame() %>% rownames_to_column(var="gene") %>% as_tibble() %>% filter(padj < 0.05)


## Extract up and downregulated for EnrichR ####
genes_up_age_pc <- sig_res_lrt_pc %>% filter(log2FoldChange > 0) %>% pull(gene)
genes_down_age_pc <- sig_res_lrt_pc %>% filter(log2FoldChange < 0) %>% pull(gene)

#remove pseudogenes
genes_up_age_pc <- genes_up_age_pc[!grepl("^Gm", genes_up_age_pc)]
genes_down_age_pc <- genes_down_age_pc[!grepl("^Gm", genes_down_age_pc)]



## EnrichR results with significantly enriched pathways - upregulated genes PCs ####
### Biological process (UP) ####
#import enrichr results 
pcageup_enrichr_BP <- read.delim("PC_age_UP_GO_Biological_Process_2026_table.txt")
pcageup_enrichr_BP <- pcageup_enrichr_BP[pcageup_enrichr_BP$Adjusted.P.value<0.1,]
pcageup_enrichr_BP <- pcageup_enrichr_BP %>% add_column(Group=1)
pcageup_enrichr_BP$GenesList <- strsplit(pcageup_enrichr_BP$Genes, split = ";\\s*")

# add columns with log-adj. p-value and genes lists
pcageup_enrichr_BP <- pcageup_enrichr_BP %>% mutate(nGenes = sapply(GenesList, length), neglog10_p = -log10(Adjusted.P.value))

pcageup_enrichr_BP$Term <- str_wrap(pcageup_enrichr_BP$Term, width = 25)

# order by significance
pcageup_enrichr_BP <- pcageup_enrichr_BP %>% mutate(Term = reorder(Term, neglog10_p))

#plot
ggplot(pcageup_enrichr_BP, aes(x = neglog10_p, y = Term, size = nGenes)) +
  geom_point(show.legend = T, alpha=.8) +
  labs(y="", x = "-log(adj. p)", size = "Number of Genes", title = "Upregulated in Young PCs (BP)") +
  theme_classic2(base_size = 20) +
  theme(axis.text.y = element_text(size = 22), axis.text.x = element_text(size = 18))

ggsave("dotplot_DE_up_PC_age_BP.pdf", plot = ggplot2::last_plot(), device = cairo_pdf, width = 13, height =10)



### Cellular component (UP) ####
#import enrichr results 
pcageup_enrichr_CC <- read.delim("PC_age_UP_GO_Cellular_Component_2026_table.txt")
pcageup_enrichr_CC <- pcageup_enrichr_CC[pcageup_enrichr_CC$Adjusted.P.value<0.1,]
pcageup_enrichr_CC <- pcageup_enrichr_CC %>% add_column(Group=1)
pcageup_enrichr_CC$GenesList <- strsplit(pcageup_enrichr_CC$Genes, split = ";\\s*")

# add columns with log-adj. p-value and genes lists
pcageup_enrichr_CC <- pcageup_enrichr_CC %>% mutate(nGenes = sapply(GenesList, length), neglog10_p = -log10(Adjusted.P.value))

pcageup_enrichr_CC$Term <- str_wrap(pcageup_enrichr_CC$Term, width = 25)

# order by significance
pcageup_enrichr_CC <- pcageup_enrichr_CC %>% mutate(Term = reorder(Term, neglog10_p))

#plot
ggplot(pcageup_enrichr_CC, aes(x = neglog10_p, y = Term, size = nGenes)) +
  geom_point(show.legend = T, alpha=.8) +
  labs(y="", x = "-log(adj. p)", size = "Number of Genes", title = "Upregulated in Young PCs (CC)") +
  theme_classic2(base_size = 20) +
  theme(axis.text.y = element_text(size = 22), axis.text.x = element_text(size = 18))

ggsave("dotplot_DE_up_PC_age_CC.pdf", plot = ggplot2::last_plot(), device = cairo_pdf, width = 13, height =10)



## Extract and plot genes from enriched pathways ####
pc_age_colors <- c("Young" = "#4f97c7", "Old" = "darkblue")
dir.create("Vln_genes_pathways_Age_PC", showWarnings = FALSE)

### Genes from GO Biological Process (upregulated) ####
pc_up_BP <- process_pathway_violin_set(enrichr_df = pcageup_enrichr_BP, norm_counts = pcnorm,
                                       metadata = metadata_pc_merge, sig_res = sig_res_lrt_pc,
                                       suffix = "PC_upreg_age_BP", outdir = "Vln_genes_pathways_Age_PC",
                                       colors = pc_age_colors)

### Genes from GO Cellular Component (upregulated) ####
pc_up_CC <- process_pathway_violin_set(enrichr_df = pcageup_enrichr_CC, norm_counts = pcnorm,
                                       metadata = metadata_pc_merge, sig_res = sig_res_lrt_pc,
                                       suffix = "PC_upreg_age_CC", outdir = "Vln_genes_pathways_Age_PC",
                                       colors = pc_age_colors)

### P-values from LRT ####
pvalue_genes_pathways_pc <- bind_rows(pc_up_BP$pvals, pc_up_CC$pvals)



## Other GOIs from literature ####
goi_pc <- c("Gria3", "Atxn1", "Camk2b","Camk2a", "Foxp1", "Grin1", "Alkbh5", "Fto", "Nova2", "Adam22",
             "Grik4", "Bdnf", "Camk4", "Cacng2", "Cspg5", "Prkcg", "Camk2a", "Grm1", "Plcb4", "Cacna1a",
             "Adgrb3", "Arc", "Grn", "Sema3a", "Sema7a", "Grid2", "Htr3a", "Igf1", "Nlgn1", "Nlgn2", "Nlgn3",
             "Ntrk2", "Plxna4", "Myo5a", "Rora", "Ptprd", "Fmr1", "Auts2")

# subset normalized counts
goi_pc_counts <- pcnorm[rownames(pcnorm) %in% goi_pc,]
goi_pc_counts <- goi_pc_counts[,metadata_pc_merge$Library]

# make long format
goi_pc_counts_long <- goi_pc_counts %>% as.data.frame() %>% rownames_to_column("gene") %>%
  pivot_longer(cols = -gene, names_to = "Library", values_to = "expr")

# merge with metadata
goi_pc_counts_long <- goi_pc_counts_long %>% left_join(metadata_pc_merge, by = "Library")
goi_pc_counts_long$Age_dis <- factor(goi_pc_counts_long$Age_dis, levels = c("Young", "Old"))

# plot
dir.create("Vln GOIs PC", showWarnings = FALSE)
genes_to_plot <- unique(goi_pc_counts_long$gene)

for (g in genes_to_plot){
  df_gene <- goi_pc_counts_long %>% filter(gene == g)
  
  p <- ggplot(df_gene, aes(x = Age_dis, y = expr, fill = Age_dis)) +
    geom_half_violin(trim = FALSE, alpha = 0.5, color = "transparent") +
    geom_half_boxplot(width = 0.2, alpha = 0.6, outlier.shape = NA, side="r") +
    geom_jitter(width = 0.08, size = 3, alpha = 0.7) +
    scale_fill_manual(values = pc_age_colors) +
    labs(x = NULL, y = "Norm. expression", title = g) +
    theme_classic2(base_size = 45)+
    theme(legend.position = "none", plot.title = element_text(hjust = .5, face = "italic"))
  
  ggsave(paste0("./Vln GOIs PC/Vln_GOI_PC_", g, "_.pdf"), plot = p, device = cairo_pdf, width = 8, height = 8)
}

pval_goi_pc <- sig_res_lrt_pc[sig_res_lrt_pc$gene %in% goi_pc_counts_long$gene,]



## Volcano plot PCs ####
# Change padj = 0 by the minimum padj to allow graphic representation
volcano_df_pc <- res_lrt_pc %>% as.data.frame() %>% rownames_to_column("gene") %>% filter(!is.na(log2FoldChange), !is.na(padj)) %>%
  mutate(padj_plot = ifelse(padj == 0, min(padj[padj > 0], na.rm = TRUE), padj), neg_log10_padj = -log10(padj_plot), 
  regulation = case_when(padj < 0.05 & log2FoldChange > 0.5  ~ "Up", padj < 0.05 & log2FoldChange < -0.5 ~ "Down", TRUE ~ "NS"))


# select genes to highlight from pathways analysis
genes_to_label_pc <- c("Grik4", "Hdac4", "Foxc1", "Cacna1h", "Epha7", "Celsr3", "Cadm1", "Ahi1", "Ctnnb1",
                       "Prkn","Sirt1")

ggplot(volcano_df_pc, aes(x = log2FoldChange, y = neg_log10_padj)) +
  ggrastr::rasterise(geom_point(aes(color = regulation), alpha = 0.6, size = 2), dpi=300)+
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", color = "grey40") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40") +
  geom_text_repel(data = volcano_df_pc %>% filter(gene %in% genes_to_label_pc), aes(label = gene),
                  size = 6, box.padding = 0.8, point.padding = 0.6, force = 5, force_pull = 0.5, max.overlaps = Inf,
                  min.segment.length = 0, segment.color = "black", segment.size = .6) +
  scale_color_manual(values = c("Up" = "#4f97c7", "Down" = "darkblue", "NS" = "grey70")) +
  labs(x = "log2 fold change", y = "-log10 adjusted p-value", color = NULL) +
  theme_classic2(base_size = 30)
