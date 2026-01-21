rm(list = ls())

library(dplyr)
library(reshape2)
library(knockoff)
library(ggplot2)
library(ggrepel)

source("UtilityFunctions.R")

# test = c("mast","wrt","lrt")
test = "mast"
# ll = c("ll00","ll01","ll03","ll05")
ll = "ll01"

#lfc = c(0,0.05,0.1,0.25)*100 # log2 fold change
lfc = 0.25*100


PC = c(43) 

cellid = "astrocyte"

q = 0.1

# gene_name_compare = "Grubman651" Compare gene set
gene_name_compare = "Top45AD"

#### List of gene names to compare ####
top45AD = c("APOE","BIN1","MS4A6A","PICALM","CR1","CLU","TREM2","ABCA7","NYAP1","PTK2B","PLCG2",
            "SPI1","SORL1","HLA-DRB1","CD2AP","SLC24A4","RIN3","ADAMTS1","ADAMTS4","CASS4","ADAM10",
            "FERMT2","HAVCR2","SCIMP","CLNK","ECHDC3","TNIP1","ABCA1","CNTNAP2","USP6NL","INPP5D",
            "CD33","ACE","IQCK","WWOX","ABI3","HESX1","FHL2","APH1B","HS3ST1","CHRNE","CCDC6","AGRN",
            "KAT8","IL34")


# load the BH procedure results
load(file = "HuVascAD_cellastrocyte_testUseBH_downNULL_gene23537_PC43_10llam1_seed17_maxIterImp100.rda") 

### Need to change if use a different test
row_index <- which(lfc/100 == rownames(FinalResult[[1]]$q_value_mast))
p_val_match = FinalResult[[1]]$selection_result_mast$p_val # this one should be correct.
q_BH_lfc = FinalResult[[1]]$q_value_mast[row_index, ]

selected_match <- which(q_BH_lfc <=q)

### dataframe for the plot 
# The gene set to compare with
gene_name_list <- top45AD #Grubman651#
gg_dat = data.frame(lpval = (p_val_match),
                    qval_BH = q_BH_lfc
                    #eval = E,
                    #ksmean = apply(K_imp, 2, mean)
)


# All cellID, all knockoffs: 
# format: (list(W_imp, W_imp_b, feature.names.new))

### Load multidecomp results

# method = c("eBH","multi")
method = "multi"
# knockoff = c("decomp","LR")
knockoff = "decomp"#"LR"

if(method == "eBH"){
  knk_method = paste0("e",knockoff)
} else {
  knk_method = paste0(method,knockoff) 
}

load("HuVascAD_cellastrocyte_methodmultidecomp_testUseMAST_downNULL_gene23537_PC43_10llam1_seed2017_maxIterImp100.rda")

# if(knk_method == "eLR"){
#   load(file = paste0("HuVascAD_downNULL_geneNULL_scsoftImpute_err_PC",PC,
#                      "_expresc0_",
#                      ll,"_",
#                      "G_imp","_",
#                      "multiLR_",
#                      test,
#                      "hippo_only",FALSE,"_",
#                      cellid,"_.rda"))
# } else if(knk_method == "multiLR"){
#   load(file = paste0("HuVascAD_downNULL_geneNULL_scsoftImpute_err_PC",PC,
#                      "_expresc0_",
#                      ll,"_",
#                      "G_imp","_",
#                      "multiLR_",
#                      test,
#                      "hippo_only",FALSE,"_",
#                      cellid,"_.rda"))
# } else if(knk_method == "multidecomp"){
#   load(file = paste0("HuVascAD_downNULL_geneNULL_scsoftImpute_err_PC",PC,
#                      "_expresc0_",
#                      ll,"_",
#                      "G_imp","_",
#                      "multidecomp_",
#                      test,
#                      "hippo_only",FALSE,"_",
#                      cellid,"_.rda"))
# } else if(knk_method == "edecomp"){
#   load(file = paste0("HuVascAD_downNULL_geneNULL_scsoftImpute_err_PC",PC,
#                      "_expresc0_",
#                      ll,"_",
#                      "G_imp","_",
#                      "multidecomp_",
#                      test,
#                      "hippo_only",FALSE,"_",
#                      cellid,"_.rda"))
# } else {
#   print("Error: knockoff method - no match?")
# }

W_imp = FinalResult[[1]]$W_imp_b # extract feature importance statistics
feature_names = FinalResult[[1]]$feature.names.new # extract gene names

m = 10 # since we will fix M at 10.

#e-values
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

# Calculate feature importance statistics for conditional inference
K_imp <- matrix(NA, m, dim(W_imp)[2])
for(w_it in 1:m){
  K_imp[w_it,] <- W_imp[1,] - W_imp[(w_it+1),]
}

### Add e-values to the plot data
gg_dat$eval =E

# Select variables, and calculate q-values, saved in q_selected_method
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

# Add q-values to the plot data
varname = paste0("qval_", knk_method)
gg_dat[[varname]] = q_selected_method


#gg_dat$bh_sel[selected_match] = 1
#gg_dat$ebh_sel[selected_method] = 1

# sum(selected_match %in% selected_method)
# sum(selected_method %in% selected_match)

only_method = selected_method[-which(selected_method %in% selected_match)]
only_bh = selected_match[-which(selected_match %in% selected_method)]
both_method_bh = selected_match[which(selected_match %in% selected_method)]

# This show out of the gene list, which genes are only identified by the proposed method
gene_name_list[ gene_name_list %in% FinalResult[[1]]$feature.names.new[only_method]]
# Only BH
gene_name_list[ gene_name_list %in% FinalResult[[1]]$feature.names.new[only_bh]]
# Both 
gene_name_list[ gene_name_list %in% FinalResult[[1]]$feature.names.new[both_method_bh]]

# color for the points
varname = paste0("color_code_", knk_method)
gg_dat[[varname]] = 0
gg_dat[[varname]][only_method] = 1
gg_dat[[varname]][only_bh] = 2
gg_dat[[varname]][both_method_bh] = 3

gg_dat[[varname]] = factor(gg_dat[[varname]], levels=c(0,1,2,3))
#gg_dat$sel_times = factor(gg_dat$sel_times)

# Change to edecomp 
method <- "eBH"
if(method == "eBH"){
  knk_method = paste0("e",knockoff)
} else {
  knk_method = paste0(method,knockoff) 
}

# Select variables, and calculate q-values, saved in q_selected_method
if(method == "eBH"){
  selected_method = ebh(E,q)$rej
  method_cut = ebh(E,q)$cut
  system.time({
  q_selected_method = qebh(W_imp)})
} else if(method == "multi"){
  selected_method = multiknk(W_imp,q)$rej
  method_cut = multiknk(W_imp,q)$cut
  system.time({
    q_selected_method = qmulti(W_imp)
  })
  
} else {
  print("Error: method - selection or qvalue error?")
}

varname = paste0("qval_", knk_method)
gg_dat[[varname]] = q_selected_method


#gg_dat$bh_sel[selected_match] = 1
#gg_dat$ebh_sel[selected_method] = 1

# sum(selected_match %in% selected_method)
# sum(selected_method %in% selected_match)

only_method = selected_method[-which(selected_method %in% selected_match)]
only_bh = selected_match[-which(selected_match %in% selected_method)]
both_method_bh = selected_match[which(selected_match %in% selected_method)]

# This show out of the gene list, which genes are only identified by the proposed method
gene_name_list[ gene_name_list %in% FinalResult[[1]]$feature.names.new[only_method]]
# Only BH
gene_name_list[ gene_name_list %in% FinalResult[[1]]$feature.names.new[only_bh]]
# Both 
gene_name_list[ gene_name_list %in% FinalResult[[1]]$feature.names.new[both_method_bh]]

# color for the points
varname = paste0("color_code_", knk_method)
gg_dat[[varname]] = 0
gg_dat[[varname]][only_method] = 1
gg_dat[[varname]][only_bh] = 2
gg_dat[[varname]][both_method_bh] = 3

gg_dat[[varname]] = factor(gg_dat[[varname]], levels=c(0,1,2,3))


genecompare_index = match(gene_name_list,feature_names)
print(paste0("Number of the genes that don't exist in our data: ",
             sum(is.na(genecompare_index)))) # some of the identified genes don't even exist
genecompare_index_NA = is.na(genecompare_index)
genecompare_index = genecompare_index[!genecompare_index_NA] # and I need to exclude them

gg_dat$gene_name_list = ""
gg_dat$gene_name_list[genecompare_index] = gene_name_list[!genecompare_index_NA]

#### Plot ####

library(ggrepel)
library(patchwork)

gg_dat$qval_BH = q_BH_lfc
### for gene that has been filtered out by lfc, the p-values are set to 1
gg_dat$qval_BH[gg_dat$qval_BH==2] <- 1
### extremely small q-values are lowered truncated
gg_dat$qval_BH[gg_dat$qval_BH < 1e-52] <- 1e-52

# Plot for multiple knockoffs
y_var   <- paste0("qval_", "multi", knockoff)
col_var <- paste0("color_code_", "multi", knockoff)

label_dat <- gg_dat[
  gg_dat$gene_name_list != "" &
    !is.na(gg_dat$gene_name_list) &
    gg_dat[[col_var]] != "0",
]

p_A <- ggplot(
  gg_dat,
  aes(
    x = qval_BH,
    y = .data[[y_var]],
    color = factor(.data[[col_var]])
  )) +
    xlab("log-transformed q-value (BH)") +
    ylab(paste0("q-value (LaCIS, multi-", knockoff, ")")) +   
    scale_color_brewer(name="Selected in",
                       breaks=c("0","1", "2", "3"),
                       labels=c("None","only LaCIS","only BH", "BH and LaCIS"),
                       palette="Set2") +
  geom_point(size = 1, alpha = 0.8) +
  scale_x_continuous(trans = "log10") +
  geom_text_repel(
    data = label_dat,
    aes(label = gene_name_list),
    size = 3,
    fontface = "bold",
    box.padding = 1.5,
    point.padding = 0.2,
    show.legend = FALSE
  ) +
  theme_bw()

p_A

# Plot for the eBH procedure
y_var   <- paste0("qval_", "e", knockoff)
col_var <- paste0("color_code_", "e", knockoff)

label_dat <- gg_dat[
  gg_dat$gene_name_list != "" &
    !is.na(gg_dat$gene_name_list) &
    gg_dat[[col_var]] != "0",
]

p_B <- ggplot(
  gg_dat,
  aes(
    x = qval_BH,
    y = .data[[y_var]],
    color = factor(.data[[col_var]])
  )) +
  xlab("log-transformed q-value (BH)") +
  ylab(paste0("q-value (LaCIS, e-", knockoff, ")")) +   
  scale_color_brewer(name="Selected in",
                     breaks=c("0","1", "2", "3"),
                     labels=c("None","only LaCIS","only BH", "BH and LaCIS"),
                     palette="Set2") +
  geom_point(size = 1, alpha = 0.8) +
  scale_x_continuous(trans = "log10") +
  geom_text_repel(
    data = label_dat,
    aes(label = gene_name_list),
    size = 3,
    fontface = "bold",
    box.padding = 1.5,
    point.padding = 0.2,
    show.legend = FALSE
  ) +
  theme_bw()

p_B

(p_A | p_B) +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "A") &
  theme(legend.position = "bottom",
        legend.title = element_text(size = 12),
        legend.text  = element_text(size = 12),
  )

ggsave(
  "BH_vs_LaCIS_panels.pdf",
  (p_A | p_B) +
    plot_layout(guides = "collect") +
    plot_annotation(tag_levels = "A") &
    theme(legend.position = "bottom",
          legend.title = element_text(size = 12),
          legend.text  = element_text(size = 12),),
  width = 10,
  height = 6
)


