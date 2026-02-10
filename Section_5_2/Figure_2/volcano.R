rm(list = ls())

library(dplyr)
library(reshape2)
library(knockoff)
library(ggplot2)
library(Seurat)
library(ggrepel)
library(patchwork)

source("UtilityFunctions.R")

load("HuVascAD_cellmicroglia_testUseBH_downNULL_gene23537_PC11_10llam1_seed17_maxIterImp100_conditional_multidecomp_LR.rda")
# ko_selected <- c("SLC2A5",     "TNFRSF1B",   "SGIP1",      "MIR181A1HG", "RNF149",    
#                  "MAP4K4",     "TMEM163",    "ARHGAP15",   "SLC11A1",    "DOCK10",    
#                  "SLC4A7",     "ADAMTS9",    "FOXP1",      "AC092957.1", "P2RY12",    
#                  "MBNL1",      "TBC1D14",    "RBPJ",       "ELF2",       "ACSL1",     
#                  "LINC02211",  "F13A1",      "EPB41L2",    "TAB2",       "RPS6KA2",   
#                  "SDK1",       "HDAC9",      "NAMPT",      "AC016831.7", "CALD1",     
#                  "CTSB",       "BNIP3L",     "DENND3",     "DOCK8",      "FTH1",     
#                  "NEAT1",      "MALAT1",     "SFMBT2",     "PLXDC2",     "MAP3K8",    
#                  "AC074327.1", "ETV6",       "IRAK3",      "CHST11",     "FLT3",      
#                  "FLT1",       "PCDH9",      "HIF1A-AS2",  "GLDN",       "CHSY1",     
#                  "FAM157C",    "PRKCA",      "GRB2",       "RAB31",      "SOCS6",     
#                  "BSG",        "FP671120.1", "PCBP3"    )

### Conditional on knockoff-selected genes

deg <- FinalResult[[1]]$selection_result_mast_cond

deg[FinalResult[[1]]$ko_selected_ind, ]

deg[FinalResult[[1]]$ko_selected_ind, "p_val"]<-1

deg$p_val_adj <- p.adjust(deg$p_val, method = "BH")

#deg <- result.rand.mast.reordered
deg$gene <- rownames(deg)

# Avoid infinite -log10(0)
deg$p_val_adj[deg$p_val_adj == 0] <- .Machine$double.xmin

# Add significance categories
logfc_cutoff <- 0.0
padj_cutoff  <- 0.1 #1e-5#

deg <- deg %>%
  mutate(
    neg_log10_padj = -log10(p_val_adj),
    sig = case_when(
      avg_log2FC >= logfc_cutoff & p_val_adj < padj_cutoff ~ "Up",
      avg_log2FC <= -logfc_cutoff & p_val_adj < padj_cutoff ~ "Down",
      TRUE ~ "Not selected"
    )
  )

top_genes <- deg %>%
  arrange(p_val_adj) %>%
  head(20)

p_B<- ggplot(deg, aes(x = avg_log2FC, y = neg_log10_padj, color = sig)) +
  scale_color_manual(values = c("Down" = "royalblue3", "Not selected" = "grey70", "Up" = "violetred2"))+
  geom_point(alpha = 0.7, size = 1.2) +
  geom_text_repel(
    data = top_genes,
    aes(label = gene),
    size = 2.5,
    max.overlaps = Inf
  ) +
  geom_vline(xintercept = c(-logfc_cutoff, logfc_cutoff), linetype = "dashed", col="grey70") +
  geom_hline(yintercept = -log10(padj_cutoff), linetype = "dashed",col="grey70") +
  labs(
    x = "log2 fold change",
    y = "-log10(q-value)",
    color = ""
  ) +
  theme_classic()+
  theme(
    # axis titles
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12)
  )+
  coord_cartesian(ylim = c(0, 40))#+

p_B

### Conditional on random-selected genes

deg <- FinalResult[[1]]$selection_result_mast_rand

deg[FinalResult[[1]]$random_ind, ]

deg[FinalResult[[1]]$random_ind, "p_val"]<-1

deg$p_val_adj <- p.adjust(deg$p_val, method = "BH")

#deg <- result.rand.mast.reordered
deg$gene <- rownames(deg)

# Avoid infinite -log10(0)
deg$p_val_adj[deg$p_val_adj == 0] <- .Machine$double.xmin

# Add significance categories
logfc_cutoff <- 0.0
padj_cutoff  <- 0.1 #1e-5#

deg <- deg %>%
  mutate(
    neg_log10_padj = -log10(p_val_adj),
    sig = case_when(
      avg_log2FC >= logfc_cutoff & p_val_adj < padj_cutoff ~ "Up",
      avg_log2FC <= -logfc_cutoff & p_val_adj < padj_cutoff ~ "Down",
      TRUE ~ "Not selected"
    )
  )

top_genes <- deg %>%
  arrange(p_val_adj) %>%
  head(20)

p_C<- ggplot(deg, aes(x = avg_log2FC, y = neg_log10_padj, color = sig)) +
  scale_color_manual(values = c("Down" = "royalblue3", "Not selected" = "grey70", "Up" = "violetred2"))+
  geom_point(alpha = 0.7, size = 1.2) +
  geom_text_repel(
    data = top_genes,
    aes(label = gene),
    size = 2.5,
    max.overlaps = Inf
  ) +
  geom_vline(xintercept = c(-logfc_cutoff, logfc_cutoff), linetype = "dashed", col="grey70") +
  geom_hline(yintercept = -log10(padj_cutoff), linetype = "dashed",col="grey70") +
  labs(
    x = "log2 fold change",
    y = "-log10(q-value)",
    color = ""
  ) +
  theme_classic()+
  theme(
    # axis titles
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12)
  )+
  coord_cartesian(ylim = c(0, 40))#+
p_C

remove(FinalResult)

load("HuVascAD_cellmicroglia_testUseBH_downNULL_gene23537_PC11_10llam1_seed17_maxIterImp100.rda")

deg <- FinalResult[[1]]$selection_result_lr

deg$p_val_adj <- p.adjust(deg$p_val, method = "BH")

deg$gene <- rownames(deg)

# Avoid infinite -log10(0)
deg$p_val_adj[deg$p_val_adj == 0] <- .Machine$double.xmin

# Add significance categories
logfc_cutoff <- 0.0
padj_cutoff  <- 0.1 #1e-5#

deg <- deg %>%
  mutate(
    neg_log10_padj = -log10(p_val_adj),
    sig = case_when(
      avg_log2FC >= logfc_cutoff & p_val_adj < padj_cutoff ~ "Up",
      avg_log2FC <= -logfc_cutoff & p_val_adj < padj_cutoff ~ "Down",
      TRUE ~ "Not selected"
    )
  )

top_genes <- deg %>%
  arrange(p_val_adj) %>%
  head(20)

p_A<- ggplot(deg, aes(x = avg_log2FC, y = neg_log10_padj, color = sig)) +
  scale_color_manual(values = c("Down" = "royalblue3", "Not selected" = "grey70", "Up" = "violetred2"))+
  geom_point(alpha = 0.7, size = 1.2) +
  geom_text_repel(
    data = top_genes,
    aes(label = gene),
    size = 2.5,
    max.overlaps = Inf
  ) +
  geom_vline(xintercept = c(-logfc_cutoff, logfc_cutoff), linetype = "dashed", col="grey70") +
  geom_hline(yintercept = -log10(padj_cutoff), linetype = "dashed",col="grey70") +
  labs(
    x = "log2 fold change",
    y = "-log10(q-value)",
    color = ""
  ) +
  theme_classic()+
  coord_cartesian(ylim = c(0, 40))

p_A

(p_A | p_B | p_C) +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "A") &
  theme(legend.position = "bottom",
        legend.title = element_text(size = 12),
        legend.text  = element_text(size = 12),)

ggsave(
  "conditional_microglia.pdf",
  (p_A | p_B | p_C) +
    plot_layout(guides = "collect") +
    plot_annotation(tag_levels = "A") &
    theme(legend.position = "bottom",
          legend.title = element_text(size = 12),
          legend.text  = element_text(size = 12),),
  width = 10,
  height = 4.5
)



