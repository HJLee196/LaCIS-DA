rm(list=ls())

############################################################
## 0. Packages
############################################################

# BiocManager::install("AnnotationDbi")
# BiocManager::install("org.Hs.eg.db")
# BiocManager::install("dorothea")
# BiocManager::install("igraph")
# BiocManager::install("ggraph")

# Use org.Hs.eg.db version 3.21.0 for reproducibility.
# .libPaths(c("~/R/old_bioc_lib", .libPaths()))
library(org.Hs.eg.db)
# Verify the installed package version.
packageVersion("org.Hs.eg.db")
# Check the path from which the package is loaded, if needed.
# find.package("org.Hs.eg.db")

library(dplyr)
library(tibble)
library(AnnotationDbi)
# library(org.Hs.eg.db)
library(dorothea)
library(igraph)
library(ggraph)
library(ggplot2)

source("UtilityFunctions.R")

############################################################
## 1. Differential analysis table
############################################################

######## Step 1 ########
# Find the selected DEGs
########################
library(knockoff)
test =  "mast" #"LR" #
# test = "lrt"
# ll = c("ll00","ll01","ll03","ll05")
ll = "ll01"
# method = c("eBH","multi")
method = "multi"
# method = "eBH"
# knockoff = c("decomp","LR")
knockoff = "decomp" #"LR"

if(method == "eBH"){
  knk_method = paste0("e",knockoff)
} else {
  knk_method = paste0(method,knockoff) 
}

#lfc = c(0,0.05,0.1,0.25)*100 # log2 fold change
# lfc = c(0.07,0.14,0.36,0.72) # log2 fold change to log fold change
# lfc = 0.14*100
lfc = 0.25*100
# lfc = 0.25*100

# G_imp = c("G_imp")
# sgn = c("03")
# hippo_only = c(FALSE)

# PC = c(115)
# cellid = c("astrocyte", "oligo","pericyte")
# PC = c(30,50)
# cellid = c("microglia", "neuron") 
# PC = c(43) # c(60) c(115)
cellid = "astrocyte"
# cellid = "oligo"
# PC = 50
#cellid = "microglia"
PC = 43
# cellid = "neuron"

# q = c(0.05,0.1) # Target FDR
q = 0.1

if(knk_method == "eLR"){
  load(file = paste0("HuVascAD_downNULL_geneNULL_scsoftImpute_err_PC",PC,
                     "_expresc0_",
                     ll,"_",
                     "G_imp","_",
                     "multiLR_",
                     test,
                     "hippo_only",FALSE,"_",
                     cellid,"_.rda"))
} else if(knk_method == "multiLR"){
  load(file = paste0("HuVascAD_downNULL_geneNULL_scsoftImpute_err_PC",PC,
                     "_expresc0_",
                     ll,"_",
                     "G_imp","_",
                     "multiLR_",
                     test,
                     "hippo_only",FALSE,"_",
                     cellid,"_.rda"))
} else if(knk_method == "multidecomp"){
  load(file = paste0("HuVascAD_downNULL_geneNULL_scsoftImpute_err_PC",PC,
                     "_expresc0_",
                     ll,"_",
                     "G_imp","_",
                     "multidecomp_",
                     test,
                     "hippo_only",FALSE,"_",
                     cellid,"_.rda"))
} else if(knk_method == "edecomp"){
  if(cellid == "astrocyte" | cellid == "oligo"){
    load(file = paste0("HuVascAD_downNULL_geneNULL_scsoftImpute_err_PC",PC,
                       "_expresc0_",
                       ll,"_",
                       "G_imp","_",
                       "multidecomp_",
                       test,
                       "hippo_only",FALSE,"_",
                       cellid,"_.rda"))
  } else if(cellid == "microglia" | cellid == "neuron"){
  
  } else if(cellid == "pericyte"){

  } else {
    print("Error: edecomp - CellID issue?")
  }
} else {
  print("Error: knockoff method - no match?")
}

W_imp = FinalResult[[1]]$W_imp_b# Always bonferroni corrected
feature_names = FinalResult[[1]]$feature.names.new

m = 10 # since we will fix M at 10.

E <- matrix(NA, m, dim(W_imp)[2])
ind_sel = rep(NA,m)
ind_sel_which = matrix(0,nrow = m, ncol = dim(W_imp)[2])
for(w_it in 1:m){
  E[w_it,] <- W_imp[1,] - W_imp[(w_it+1),]
  tau <- knockoff.threshold(E[w_it,],fdr = q/2)
  ind_sel[w_it] = sum(E[w_it,] >= tau)
  ind_sel_which[w_it,which(E[w_it,] >= tau)] = 1
  
  E[w_it,] <- (E[w_it,] >= tau) / (1 + sum(E[w_it,] <= -tau))
}

ind_sel_which_sum = colSums(ind_sel_which)
E <- (dim(W_imp)[2]) * colMeans(E)

# Average knockoff values
K_imp <- matrix(NA, m, dim(W_imp)[2])
for(w_it in 1:m){
  K_imp[w_it,] <- W_imp[1,] - W_imp[(w_it+1),]
}

if(method == "eBH"){
  selected_method = ebh(E,q)$rej
  method_cut = ebh(E,q)$cut
  q_selected_method = qebh(W_imp)
  
} else if(method == "multi"){
  selected_method = multiknk(W_imp,q)$rej
  method_cut = multiknk(W_imp,q)$cut
  q_selected_method = qmulti(W_imp)
  
} else {
  print("Error: method - selection or qvalue error?")
}

deg_df <- tibble(
  gene = FinalResult[[1]]$feature.names.new,
  qval_ko = q_selected_method,
  padj = q_selected_method
)

## Load BH / MAST results
load("HuVascAD_downNULL_geneNULL_BH_PC43_hippo_onlyFALSE_astrocyte_.rda")

deg_df <- deg_df |>
  mutate(
    #row 1 is logFC =0
    #row 7 is logFC =0.1
    #row 10 is logFC = 0.25
    padj = FinalResult[[1]]$q_value_mast[10, ],#FinalResult[[1]]$q_value_mast[1, ],
    logFC = FinalResult[[1]]$selection_result_mast$avg_log2FC,
    selected_method = qval_ko < 0.1,
    selected_bh = padj < 0.1,
    selected_bh_logFC_0 = FinalResult[[1]]$q_value_mast[10, ] < 0.1
  )

#deg_df$selected_bh_logFC_0 = FinalResult[[1]]$q_value_mast[1, ] < 0.1
#deg_df$padj <- q_selected_method

############################################################
## 2. Extract genes in GO:0007169
## cell surface receptor protein tyrosine kinase signaling pathway
############################################################

go_id <- "GO:0007169"
#go_id <- "GO:0050808" #synapse organization
#go_id <- "GO:0001525"  # angiogenesis
#go_id <- "GO:0001568" #blood vessel development 
#go_id <- "GO:0030301" #cholesterol transport 
#go_id <- "GO:0034375" #lipoprotein particle remodeling 
#go_id <- "GO:0006869" #lipid transport 
#go_id <- "GO:0042632" #cholesterol homeostasis 

### For pathway with APOE

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

#BiocManager::install("KEGGREST")

# library(KEGGREST)
# ppar_pathway <- KEGGREST::keggGet("hsa03320")[[1]]
# 
# ppar_genes <- tibble(
#   gene = sub(
#     ";.*",
#     "",
#     ppar_pathway$GENE[seq(2, length(ppar_pathway$GENE), by = 2)]
#   )
# ) |>
#   distinct(gene)
# 
# pathway_gene_symbols <- unique(ppar_genes$gene)
# 
# length(pathway_gene_symbols)
# head(pathway_gene_symbols)

### For GO pathways
pathway_genes <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = go_id,
  keytype = "GOALL",
  columns = c("SYMBOL", "ENTREZID")
) |>
  as_tibble() |>
  dplyr::filter(!is.na(SYMBOL)) |>
  dplyr::distinct(SYMBOL, .keep_all = TRUE) |>
  dplyr::rename(gene = SYMBOL)

pathway_gene_symbols <- pathway_genes$gene

# Check For reproducibility
"ARID5B" %in% pathway_gene_symbols
"EBF1" %in% pathway_gene_symbols
"KLF6" %in% pathway_gene_symbols
"BCL6" %in% pathway_gene_symbols

# EBF1
"FYN" %in% pathway_gene_symbols
# KLF6
"IGF1R" %in% pathway_gene_symbols
# BCL6
"RHOQ" %in% pathway_gene_symbols


############################################################
## 3. Load DoRothEA TF-target database
############################################################

data(dorothea_hs)

# dorothea_hs |>
#   dplyr::filter(target == "APOE") |>
#   dplyr::count(confidence)

tf_target_all <- dorothea_hs |>
  as_tibble() |>
  dplyr::filter(confidence %in% c("A", "B", "C")) |> # How confident the data are about the relationships
  dplyr::select(tf, target, mor, confidence) |> # mor indicates a direction of an effect of each gene. 1: positive, -1: negative, 
  dplyr::distinct()

"ARID5B" %in% tf_target_all$target

tf_target_all[tf_target_all$target == "ARID5B",]

############################################################
## 4. Keep TF-target pairs where either TF or target is in pathway
############################################################

tf_target_pathway <- tf_target_all |>
  dplyr::filter(
    tf %in% pathway_gene_symbols |
      target %in% pathway_gene_symbols
  ) |>
  dplyr::mutate(
    tf_in_pathway = tf %in% pathway_gene_symbols,
    target_in_pathway = target %in% pathway_gene_symbols
  )

############################################################
## 5. Add differential analysis results to TFs and targets
############################################################

tf_target_annotated <- tf_target_pathway |>
  dplyr::left_join(
    deg_df |>
      dplyr::select(
        tf = gene,
        tf_logFC = logFC,
        tf_padj = padj,
        tf_qval_ko = qval_ko,
        tf_selected_method = selected_method,
        tf_selected_bh = selected_bh,
        tf_selected_bh_logFC_0 = selected_bh_logFC_0
      ),
    by = "tf"
  ) |>
  dplyr::left_join(
    deg_df |>
      dplyr::select(
        target = gene,
        target_logFC = logFC,
        target_padj = padj,
        target_qval_ko = qval_ko,
        target_selected_method = selected_method,
        target_selected_bh = selected_bh,
        target_selected_bh_logFC_0 = selected_bh_logFC_0
      ),
    by = "target"
  )

############################################################
## 6. Optional: keep only pairs where TF or target has DA result
############################################################

tf_target_annotated <- tf_target_annotated |>
  dplyr::filter(
    (tf_selected_method | tf_selected_bh_logFC_0 ) &
    (target_selected_method | target_selected_bh_logFC_0 ) &
    (!is.na(tf_logFC) | !is.na(target_logFC))
  )

############################################################
## 7. Build node annotation
############################################################

all_tfs <- unique(tf_target_annotated$tf)
all_targets <- unique(tf_target_annotated$target)
all_nodes <- unique(c(all_tfs, all_targets))

nodes <- tibble(name = all_nodes) |>
  dplyr::left_join(
    deg_df |>
      dplyr::select(
        name = gene,
        logFC,
        padj,
        qval_ko,
        selected_method,
        selected_bh
      ),
    by = "name"
  ) |>
  dplyr::mutate(
    in_pathway = name %in% pathway_gene_symbols,
    is_TF = name %in% all_tfs,
    is_target = name %in% all_targets,
    node_type = dplyr::case_when(
      is_TF & is_target ~ "TF and target",
      is_TF ~ "TF",
      is_target ~ "Target gene",
      TRUE ~ "Other"
    ),
    selection_status = dplyr::case_when(
      # name %in% c("VEGFA","ARID5B","NTRK3") ~ "Highlight",
      # TRUE ~"Not highlight"
      selected_method & selected_bh ~ "Both",
      selected_method & !selected_bh ~ "Method only",
      !selected_method & selected_bh ~ "BH only",
      TRUE ~ "Not selected"
    ),
    highlighted_node = dplyr::case_when(
      name %in% c("ARID5B", "VEGFA", "NTRK3", "EBF1", "EPAS1", "SOX2") ~ "highlighted",
      TRUE ~ "other"
    ),
    neg_log10_padj = -log10(padj)
  )

############################################################
## 8. Build edge table
############################################################

edges <- tf_target_annotated |>
  dplyr::transmute(
    from = tf,
    to = target,
    confidence = confidence,
    regulation = dplyr::case_when(
      mor == 1 ~ "activation",
      mor == -1 ~ "repression",
      TRUE ~ "unknown"
    ),
    highlighted_edge = dplyr::case_when(
      ((target == "VEGFA") & (tf == "EPAS1"))|
        ((target == "NTRK3") & (tf == "SOX2")|
           ((target == "ARID5B") & (tf == "EBF1"))) ~ "highlighted",
      TRUE ~ "other"
    ),
    tf_in_pathway,
    target_in_pathway
  )

############################################################
## 9. Construct graph
############################################################

g <- igraph::graph_from_data_frame(
  d = edges,
  vertices = nodes,
  directed = TRUE
)

############################################################
## 10. Force-directed network
############################################################

set.seed(1)

p1 <- ggraph(g, layout = "fr") +
  geom_edge_link(
    aes(linetype = regulation, color = highlighted_edge),
    arrow = arrow(length = unit(3, "mm"), type = "closed"),
    end_cap = circle(3, "mm"),
    alpha = 0.6
  ) +
  geom_node_point(
    aes(
      size = neg_log10_padj,
      color = selection_status,
      shape = node_type
    ),
    alpha = 0.9
  ) +
  geom_node_text(
    data = function(x) dplyr::filter(
      x,
      highlighted_node != "highlighted"
    ),
    aes(label = name),
    repel = TRUE,
    size = 3
  ) +
  scale_edge_linetype_manual(
    values = c(
      activation = "solid",
      repression = "dashed",
      unknown = "dotted"
    )
  ) +
  scale_color_manual(
    values = c(
      "Method only" = "#1B9E77",
      "BH only" = "#D95F02",
      "Both" = "#7570B3", #,"#3399FF"
      "Not selected" = "grey80"
    )
  ) + scale_edge_color_manual(
    values = c(
      highlighted = "red",
      other = "black"
    )
  ) + 
  # geom_node_point(
  # data = function(x) dplyr::filter(x, highlighted_node == "highlighted"),
  # aes(
  #   size = neg_log10_padj
  # ),
  # fill = NA,
  # color = "red",
  # stroke = 1.3,
  # show.legend = FALSE) + 
  geom_node_text(
    data = function(x) dplyr::filter(x, highlighted_node == "highlighted"),
    aes(label = name),
    repel = TRUE,
    size = 3.5,
    fontface = "bold",
    color = "red"
  ) +
  theme_void()

p1


# For illustration, not labled by the selection result

# p1 <- ggraph(g, layout = "fr") +
#   geom_edge_link(
#     aes(linetype = regulation),
#     arrow = arrow(length = unit(2, "mm"), type = "closed"),
#     end_cap = circle(3, "mm"),
#     alpha = 0.5
#   ) +
#   geom_node_point(
#     aes(
#       size = neg_log10_padj,
#       color = selection_status,
#       shape = node_type
#     ),
#     alpha = 0.6
#   ) +
#   geom_node_text(
#     aes(label = name),
#     repel = TRUE,
#     size = 3
#   ) +
#   scale_edge_linetype_manual(
#     values = c(
#       activation = "solid",
#       repression = "dashed",
#       unknown = "dotted"
#     )
#   )+
#   scale_color_manual(
#         values = c(
#           "Highlight" = "#0099CC",
#           "Not highlight"= "#006666"#"#0099CC"
#         ),
#         guide = "none"
#   ) +
#   theme_void()
# 
# p1

set.seed(1)

p1_subplot <- ggraph(g, layout = "fr") +
  geom_edge_link(
    aes(linetype = regulation, color = highlighted_edge),
    arrow = arrow(length = unit(2, "mm"), type = "closed"),
    end_cap = circle(3, "mm"),
    alpha = 0.5
  ) +
  geom_node_point(
    aes(
      size = neg_log10_padj,
      color = highlighted_node,
      shape = node_type
    ),
    alpha = 0.6
  ) +
  geom_node_text(
    data = function(x) dplyr::filter(
      x,
      highlighted_node != "highlighted"
    ),
    aes(label = name),
    repel = TRUE,
    size = 3
  ) +
  scale_edge_linetype_manual(
    values = c(
      activation = "solid",
      repression = "dashed",
      unknown = "dotted"
    )
  ) +
  # scale_color_manual(
  #   values = c(
  #     "Method only" = "#1B9E77",
  #     "BH only" = "#D95F02",
  #     "Both" = "#7570B3", #,"#3399FF"
  #     "Not selected" = "grey80"
  #   )
  # ) +
  scale_color_manual(
    values = c(
      "highlighted" = "#0099CC",
      "other" = "#006666"
    ),
    guide = "none"
  ) +
  scale_edge_color_manual(
    values = c(
      highlighted = "red",
      other = "black"
    ),
    guide = "none"
  ) + 
  # geom_node_point(
  # data = function(x) dplyr::filter(x, highlighted_node == "highlighted"),
  # aes(
  #   size = neg_log10_padj
  # ),
  # fill = NA,
  # color = "red",
  # stroke = 1.3,
  # show.legend = FALSE) + 
  geom_node_text(
    data = function(x) dplyr::filter(x, highlighted_node == "highlighted"),
    aes(label = name),
    repel = TRUE,
    size = 3.5,
    fontface = "bold",
    color = "red"
  ) +
  theme_void()

p1_subplot


############################################################
## 11. Hierarchical TF -> target network
############################################################

# p2 <- ggraph(g, layout = "sugiyama") +
#   geom_edge_link(
#     aes(linetype = regulation),
#     arrow = arrow(length = unit(3, "mm"), type = "closed"),
#     end_cap = circle(3, "mm"),
#     alpha = 0.6
#   ) +
#   geom_node_point(
#     aes(
#       size = neg_log10_padj,
#       color = logFC,
#       shape = node_type
#     )
#   ) +
#   geom_node_text(
#     aes(label = name),
#     repel = TRUE,
#     size = 3
#   ) +
#   scale_edge_linetype_manual(
#     values = c(
#       activation = "solid",
#       repression = "dashed",
#       unknown = "dotted"
#     )
#   ) +
#   theme_void()
# 
# p2

############################################################
## 12. Export results
############################################################

# write.csv(
#   pathway_genes,
#   "GO0007169_pathway_genes.csv",
#   row.names = FALSE
# )
# 
# write.csv(
#   tf_target_pathway,
#   "GO0007169_TF_target_pairs.csv",
#   row.names = FALSE
# )
# 
# write.csv(
#   tf_target_annotated,
#   "GO0007169_TF_target_pairs_with_DA.csv",
#   row.names = FALSE
# )
# 
# write.csv(
#   nodes,
#   "GO0007169_network_nodes_with_DA.csv",
#   row.names = FALSE
# )
# 
# write.csv(
#   edges,
#   "GO0007169_network_edges.csv",
#   row.names = FALSE
# )

############################################################
## 13. Save figures
############################################################
# setwd("/Users/linxiliu/Dropbox/Research/Knockoff_scRNAseq/RealData/direct_causal/")
# ggsave(
#   "GO0007169_TF_target_force_network_logFC_0_25_subplot.pdf",
#   p1,
#   width = 6,
#   height = 5
# )
# 
# ggsave(
#   "GO0007169_TF_target_force_network_decomp_mast_logFC_0_10_June25_2026.pdf",
#   p1,
#   width = 20,
#   height = 15
# )
# 
# ggsave(
#   "GO0050808_TF_target_force_network_logFC_0_25.pdf",
#   p1,
#   width = 20,
#   height = 15
# )
# 
# ggsave(
#   "GO0001525_TF_target_force_network_logFC_0_25.pdf",
#   p1,
#   width = 20,
#   height = 15
# )
# 
# ggsave(
#   "GO0007169_TF_target_hierarchical_network.pdf",
#   p2,
#   width = 10,
#   height = 8
# )

ggsave(
  "GO0007169_TF_target_force_network_decomp_mast_logFC_0_025_subplot_FigureS5_June29_2026.pdf",
  p1_subplot,
  width = 6.6,
  height = 5.5
)

############################################################
## 14 Count selected TFs and targets pathway by pathway
############################################################

go_terms <- tibble::tribble(
  ~go_id,      ~pathway,
  "GO:0007169", "Cell surface receptor protein tyrosine kinase signaling pathway",
  "GO:0050808", "Synapse organization",
  "GO:0001525", "Angiogenesis",
  "GO:0001568", "Blood vessel development",
  "GO:0030301", "Cholesterol transport",
  "GO:0034375", "Lipoprotein particle remodeling",
  "GO:0006869", "Lipid transport",
  "GO:0042632", "Cholesterol homeostasis"
)

summarize_one_go <- function(go_id, pathway_name) {
  
  pathway_genes <- AnnotationDbi::select(
    org.Hs.eg.db,
    keys = go_id,
    keytype = "GOALL",
    columns = c("SYMBOL", "ENTREZID")
  ) |>
    as_tibble() |>
    dplyr::filter(!is.na(SYMBOL)) |>
    dplyr::distinct(SYMBOL, .keep_all = TRUE) |>
    dplyr::rename(gene = SYMBOL)
  
  pathway_gene_symbols <- pathway_genes$gene
  
  tf_target_pathway <- tf_target_all |>
    dplyr::filter(
      tf %in% pathway_gene_symbols |
        target %in% pathway_gene_symbols
    )
  
  tf_target_annotated <- tf_target_pathway |>
    dplyr::left_join(
      deg_df |>
        dplyr::select(
          tf = gene,
          tf_selected_method = selected_method,
          tf_selected_bh = selected_bh
        ),
      by = "tf"
    ) |>
    dplyr::left_join(
      deg_df |>
        dplyr::select(
          target = gene,
          target_selected_method = selected_method,
          target_selected_bh = selected_bh
        ),
      by = "target"
    )
  
  tibble(
    GO_ID = go_id,
    Pathway = pathway_name,
    TF_method = tf_target_annotated |>
      dplyr::filter(tf_selected_method == TRUE) |>
      dplyr::pull(tf) |>
      unique() |>
      length(),
    TF_BH = tf_target_annotated |>
      dplyr::filter(tf_selected_bh == TRUE) |>
      dplyr::pull(tf) |>
      unique() |>
      length(),
    Target_method = tf_target_annotated |>
      dplyr::filter(target_selected_method == TRUE) |>
      dplyr::pull(target) |>
      unique() |>
      length(),
    Target_BH = tf_target_annotated |>
      dplyr::filter(target_selected_bh == TRUE) |>
      dplyr::pull(target) |>
      unique() |>
      length()
  )
}

pathway_summary <- purrr::map2_dfr(
  go_terms$go_id,
  go_terms$pathway,
  summarize_one_go
)

pathway_summary
