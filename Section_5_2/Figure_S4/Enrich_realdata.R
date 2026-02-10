rm(list = ls())

library(reshape2)
library(knockoff)

### Enrichment analysis
## Install (first time only)
# BiocManager::install(c("clusterProfiler","org.Hs.eg.db","ReactomePA","enrichplot"))
# BiocManager::install("fgsea", force=TRUE)
# BiocManager::install("pathview")  # optional KEGG visualization

library(tidyverse)
library(data.table)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ReactomePA)
library(enrichplot)
library(msigdbr)
library(fgsea)

### Calculate q-value for LaCIS-DA

# load the workspace, with feature importance statistics
load("HuVascAD_cellastrocyte_methodmultidecomp_testUseMAST_downNULL_gene23537_PC43_10llam1_seed2017_maxIterImp100.rda")
source("UtilityFunctions.R")
q = 0.1

### Selections made by LaCIS-DA
W_imp = FinalResult[[1]]$W_imp_b# Always bonferroni corrected
feature_names = FinalResult[[1]]$feature.names.new
selected_method = multiknk(W_imp,q)$rej
method_cut = multiknk(W_imp,q)$cut
q_selected_method = qmulti(W_imp)

set.seed(1)
deg_df <- tibble(
  gene_symbol = FinalResult[[1]]$feature.names.new,
  qval_ko = q_selected_method,
) 

## Mark significant DEGs (q = 0.1)
sig_cutoff <- 0.1
deg_sig_ko <- deg_df %>% filter(qval_ko <= sig_cutoff)

## Background universe (all tested genes)
gene_universe_symbols <- deg_df$gene_symbol


############################################################
## 1) ID conversion: SYMBOL -> ENTREZ
## clusterProfiler expects ENTREZ for many functions
############################################################
id_map_all <- bitr(gene_universe_symbols,
                   fromType = "SYMBOL",
                   toType   = "ENTREZID",
                   OrgDb    = org.Hs.eg.db) %>% distinct(SYMBOL, .keep_all = TRUE)

id_map_sig <- bitr(deg_sig_ko$gene_symbol,
                   fromType = "SYMBOL",
                   toType   = "ENTREZID",
                   OrgDb    = org.Hs.eg.db) %>% distinct(SYMBOL, .keep_all = TRUE)

entrez_universe <- id_map_all$ENTREZID
entrez_sig      <- id_map_sig$ENTREZID

############################################################
## 2) ORA with GO / KEGG / Reactome
############################################################


## 2a) GO Biological Process ORA
ego_bp <- enrichGO(gene          = entrez_sig,
                   universe      = entrez_universe,
                   OrgDb         = org.Hs.eg.db,
                   ont           = "BP",
                   pAdjustMethod = "BH",
                   qvalueCutoff  = 0.25,  # FDR
                   readable      = TRUE)   # back-convert to SYMBOL

## Quick looks
print(head(as.data.frame(ego_bp)))
print(tail(as.data.frame(ego_bp)))
nrow(as.data.frame(ego_bp))

library(dplyr)
library(stringr)

## AD-related keyword filter 
ad_keywords <- c(
  "alzheimer", "amyloid", "tau",
  "astrocyte", "glial", "gliogenesis",
  "complement", "immune", "interferon", "cytokine", "inflamm",
  "lipid", "cholesterol", "lipoprotein", "sterol", "apoe",
  "lysosom", "autophag", "endosome", "proteasome", "ubiquitin",
  "oxidative", "reactive oxygen", "ros", "mitochond", "respiratory",
  "synap", "glutamate", "neurotransmitter",
  "vascular", "blood brain", "angiogen", "endothelial"
)

## Pull GO BP ORA results
ego_bp_df <- as.data.frame(ego_bp)

## Filter to AD-associated terms by keyword match on Description
ego_bp_ad <- ego_bp_df %>%
  filter(grepl(paste(ad_keywords, collapse="|"), Description, ignore.case = TRUE))

## Add a faceting category (theme) for AD biology
ego_bp_ad <- ego_bp_ad %>%
  mutate(
    Category = case_when(
      str_detect(tolower(Description), "complement|immune|interferon|cytokine|inflamm|antigen|mhc") ~ "Immune / Complement",
      str_detect(tolower(Description), "lipid|cholesterol|lipoprotein|sterol|apoe|fatty acid")     ~ "Lipid / Cholesterol",
      str_detect(tolower(Description), "lysosom|autophag|endosome|proteasome|ubiquitin")          ~ "Lysosome / Autophagy",
      str_detect(tolower(Description), "oxidative|reactive oxygen|ros|mitochond|respiratory")     ~ "Oxidative / Mitochondria",
      str_detect(tolower(Description), "synap|glutamate|neurotransmitter")                        ~ "Synaptic / Glutamate",
      str_detect(tolower(Description), "vascular|blood brain|angiogen|endothelial")               ~ "Vascular / BBB",
      str_detect(tolower(Description), "astrocyte|glial|gliogenesis|reactive|astrogliosis")       ~ "Astrocyte Reactivity",
      TRUE                                                                                        ~ "Other"
    )
  )

#install.packages("tidytext")  # if needed
library(ggplot2)
library(tidytext)

top_per_facet <- 6  # adjust

plot_df <- ego_bp_ad %>%
  arrange(GeneRatio) %>%
  group_by(Category) %>%
  slice_head(n = top_per_facet) %>%
  ungroup()

### Directly generate a dot plot
plot_df <- ego_bp_ad %>%
  filter(Category != "Other") %>%
  group_by(Category) %>%
  arrange(p.adjust) %>%
  slice_head(n = 8) %>%
  ungroup() %>%
  mutate(
    GeneRatio_num = as.numeric(sub("/.*", "", GeneRatio)) /
      as.numeric(sub(".*/", "", GeneRatio))
  )

# width = 14, height = 11
jpeg(filename = paste0("ko_enrichment_PC43_GO_BP_ORA.jpg"), width = 14, height = 11, units = "in",
     res = 800)
ggplot(plot_df,
       aes(x = GeneRatio_num,
           y = reorder_within(
             str_wrap(Description, width = 18), 
             GeneRatio_num, 
             Category))) +
  geom_point(aes(size = Count, color = -log10(p.adjust))) +
  facet_wrap(~ Category, 
             scales = "free_y", 
             labeller = labeller(
               Category = function(x) str_wrap(x, width = 16)
               )
             ) +
  scale_y_reordered() +
  scale_color_viridis_c(name = "-log10(adj P)") +
  labs(
    x = "Gene Ratio",
    y = NULL,
    title = "AD-associated astrocyte GO BP"
  ) +
  theme_bw() +
  theme(
    strip.text = element_text(size = 16),
    axis.text.y = element_text(size =16),
    axis.text.x = element_text(size = 12),
    axis.title.x = element_text(size = 16),
    plot.title = element_text(size = 18),
    legend.title = element_text(size = 16),
    legend.text  = element_text(size = 14)
  )
dev.off()
