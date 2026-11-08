rm(list = ls())

library(tidyverse)

setwd("./")

# Read R script

# mother_file_name = "35_mother_code_decomp_confounding.R"
mother_file_name = "38_mother_code_BH_confounding.R"

Rcodetext <- readLines(mother_file_name)

testUse_set = c("LR", "MAST") #"LCD", "wilcox_limma", "LR", "MAST" decomp c("MAST", "LR") BH c("LR")
llam_set = c(0, 0.1, 0.2, 0.3) %>% as.character() #c(0, 0.1, 0.2, 0.3)
sign_strength_set = "2" #seq(from = 1, to = 3, length.out = 5) %>% as.character()

batch_option_set = c("with", "without", "permute")
CDR_impute_option_set = c("with", "without") # c("with", "without")
CDR_p_compute_option_set = c("with", "without")

nSim = 20
down = NULL
PC = 13
seed_num = 2017
# file_path = "path/to/HuVascAD_seuratcluster3.rds"
file_path = "./HuVascAD_seuratcluster3.rds"

max_iter_imp = 100

where_seed_num = grep("^seed_num[[:blank:]]*=[[:blank:]]*",Rcodetext)
where_down = grep("^down[[:blank:]]*=[[:blank:]]*",Rcodetext)
where_gene.index = grep("^gene.index[[:blank:]]*=[[:blank:]]*",Rcodetext)
where_sign_strength = grep("^sign_strength[[:blank:]]*=[[:blank:]]*",Rcodetext)
where_n_signal = grep("^n_signal[[:blank:]]*=[[:blank:]]*",Rcodetext)
where_PC = grep("^PC[[:blank:]]*=[[:blank:]]*",Rcodetext)
where_llam = grep("^llam[[:blank:]]*=[[:blank:]]*",Rcodetext)
where_nSim = grep("^nSim[[:blank:]]*=[[:blank:]]*",Rcodetext)
where_testUse = grep("^testUse[[:blank:]]*=[[:blank:]]*",Rcodetext)
where_max_iter_imp = grep("^max_iter_imp[[:blank:]]*=[[:blank:]]*",Rcodetext)

where_batch_option = grep("^batch_option[[:blank:]]*=[[:blank:]]*",Rcodetext)
where_CDR_impute_option = grep("^CDR_impute_option[[:blank:]]*=[[:blank:]]*",Rcodetext)
where_CDR_p_compute_option = grep("^CDR_p_compute_option[[:blank:]]*=[[:blank:]]*",Rcodetext)

# batch_option # "with" "without" "permute"
# CDR_impute_option # "with" "without" 
# CDR_p_compute_option # "with" "without" 

# where_libpath = grep("\\.libPaths",Rcodetext)
where_readRDS = 
  grep("^HuVascAD = readRDS", Rcodetext)

# Change Single

Rcodetext[where_nSim] =
  Rcodetext[where_nSim] %>%
  str_replace(pattern = "^nSim[[:blank:]]*=[[:blank:]]*[[:digit:]]*", 
              replacement = paste0("nSim = ", nSim))

# Rcodetext[where_llam] =
#   Rcodetext[where_llam] %>%
#   str_replace(pattern = "^llam[[:blank:]]*=[[:blank:]]*[[:digit:]\\.]*", replacement = "llam = 0.1")

Rcodetext[where_down] =
  Rcodetext[where_down] %>%
  str_replace(pattern = "(^down[[:blank:]]*=[[:blank:]]*)([[:alnum:]]*)",
              replacement = paste0("\\1", ifelse(is.null(down), yes = "NULL", no = down)))

# Rcodetext[where_libpath] =
#   Rcodetext[where_libpath] %>%
#   str_replace("^#[[:blank:]]*", replacement = "")

Rcodetext[where_PC] =
  Rcodetext[where_PC] %>%
  str_replace(pattern = "^(PC[[:blank:]]*=[[:blank:]]*)([[:digit:]\\.]*)", 
              replacement = str_c("\\1", PC))

Rcodetext[where_seed_num] =
  Rcodetext[where_seed_num] %>%
  str_replace(pattern = "^seed_num[[:blank:]]*=[[:blank:]]*[[:digit:]]*", 
              replacement = paste0("seed_num = ", seed_num))

Rcodetext[where_readRDS] = paste0("HuVascAD = readRDS(file = '", 
                                  file_path, 
                                  "')")

Rcodetext[where_max_iter_imp] = 
  Rcodetext[where_max_iter_imp] %>%
  str_replace(pattern = "(^max_iter_imp[[:blank:]]*=[[:blank:]]*)([[:digit:]]*)", 
              replacement = str_c("\\1", max_iter_imp))

# Read sh script
mother_file_name_sh =
  str_replace(mother_file_name, pattern = "\\.R", replacement = "\\.sh")
shcodetext <- readLines(mother_file_name_sh)

where_Rscript = grep("^Rscript",shcodetext)
where_jobName = grep("job-name",shcodetext)
where_mem = grep("mem", shcodetext)
where_time = grep("time", shcodetext)

# Change Single
shcodetext[where_mem] = 
  shcodetext[where_mem] %>%
  str_replace(pattern = "mem=[[:digit:]]*g",
              replacement = "mem=40g")

shcodetext[where_time] =
  shcodetext[where_time] %>%
  str_replace(pattern = "[[:digit:]]{2}:[[:digit:]]{2}:[[:digit:]]{2}",
              replacement = "23:00:00")

# Modify the parameters and save it.

if (str_detect(mother_file_name, pattern = "_BH"))
{
  testUse_set = setdiff(x = testUse_set, y = "LCD") # LCD cannot be applied to BH
}

# batch_option # "with" "without" "permute"
# CDR_impute_option # "with" "without" 
# CDR_p_compute_option # "with" "without" 

if (str_detect(mother_file_name, pattern = "_BH"))
{
  CDR_impute_option_set = setdiff(x = CDR_impute_option_set, y = "without") # We dont need this
  llam_set = 0 # We dont need this, neither.
}

for (i_testUse in 1:length(testUse_set))
{
  for (i_llam in 1:length(llam_set))
  {
    for (i_sign_strength in 1:length(sign_strength_set))
    {
      for (i_batch_option in 1:length(batch_option_set))
      {
        for (i_CDR_impute_option in 1:length(CDR_impute_option_set))
        {
          for (i_CDR_p_compute_option in 1:length(CDR_p_compute_option_set))
          {
            # Initialize
            Rcodetext_child = Rcodetext
            shcodetext_child = shcodetext
            
            # Change Parameters for the R script
            
            # testUse and Create Dir
            testUse_temp = testUse_set[i_testUse]
            Rcodetext_child[where_testUse] = str_c("testUse = \"", testUse_temp, "\"")
            
            folder_name = 
              mother_file_name %>%
              str_to_lower() %>%
              str_replace(pattern = ".r$", replacement = "")
            
            folder_name = 
              paste(folder_name, testUse_temp, sep = "_")
            
            if (!dir.exists(folder_name)){
              dir.create(folder_name)
            }
            
            # llam
            llam_temp = llam_set[i_llam]
            Rcodetext_child[where_llam] =
              Rcodetext_child[where_llam] %>% 
              str_replace(pattern = "(=[[:blank:]]*)([[:digit:]\\.]*)", replacement = str_c("\\1", llam_temp))
            
            # sign_strength
            sign_strength_temp = sign_strength_set[i_sign_strength]
            Rcodetext_child[where_sign_strength] = 
              Rcodetext_child[where_sign_strength] %>%
              str_replace(pattern = "(=[[:blank:]]*)([[:digit:]\\.]*)", replacement = str_c("\\1", sign_strength_temp))
            
            # batch_option
            batch_option_temp = batch_option_set[i_batch_option]
            Rcodetext_child[where_batch_option] =
              Rcodetext_child[where_batch_option] %>%
              str_replace(pattern = "(=[[:blank:]]*)([[:alpha:]\"]*)", 
                          replacement = str_c("\\1", '\\\"', batch_option_temp, '\\\"'))
            
            # CDR_impute_option
            CDR_impute_option_temp = CDR_impute_option_set[i_CDR_impute_option]
            Rcodetext_child[where_CDR_impute_option] =
              Rcodetext_child[where_CDR_impute_option] %>%
              str_replace(pattern = "(=[[:blank:]]*)([[:alpha:]\"]*)", 
                          replacement = str_c("\\1", '\\\"', CDR_impute_option_temp, '\\\"'))
            
            # CDR_p_compute_option
            CDR_p_compute_option_temp = CDR_p_compute_option_set[i_CDR_p_compute_option]
            Rcodetext_child[where_CDR_p_compute_option] =
              Rcodetext_child[where_CDR_p_compute_option] %>%
              str_replace(pattern = "(=[[:blank:]]*)([[:alpha:]\"]*)", 
                          replacement = str_c("\\1", '\\\"', CDR_p_compute_option_temp, '\\\"'))
            
            if ((CDR_impute_option_temp == "with") & 
                (CDR_p_compute_option_temp == "without") & 
                (!str_detect(mother_file_name, pattern = "_BH"))){
              next
            }
            
            if ( !(testUse_temp %in% c("negbinom", "poisson", "MAST", "LR")) )
            {
              where_latent.vars = grep(pattern = "latent.vars", x = Rcodetext)
              for (i_latent in 1:length(where_latent.vars))
              {
                where_latent.vars_temp = where_latent.vars[i_latent]
                Rcodetext_child[where_latent.vars_temp] =
                  Rcodetext_child[where_latent.vars_temp] %>% 
                  str_replace(pattern = "latent\\.vars[[:blank:]]*=[^\\)]*\\)",
                              replacement = "latent\\.vars = NULL\\)")
              }
            }
            
            # Save the R file.
            file_name_R =
              str_replace(mother_file_name, pattern = "mother", replacement = "child") %>%
              str_replace(pattern = "\\.R", replacement = str_c("_down", ifelse(is.null(down), yes = "NULL", no = down), "\\.R")) %>%
              str_replace(pattern = "\\.R", replacement = str_c("_testUse", testUse_temp, "\\.R")) %>%
              str_replace(pattern = "\\.R", replacement = str_c("_10llam", 10*as.numeric(llam_temp), "\\.R")) %>%
              str_replace(pattern = "\\.R", replacement = str_c("_10sgnStrength", 10*as.numeric(sign_strength_temp), "\\.R")) %>%
              str_replace(pattern = "\\.R", replacement = str_c("_batch", batch_option_temp, "\\.R")) %>%
              str_replace(pattern = "\\.R", replacement = str_c("_impCDR", CDR_impute_option_temp, "\\.R")) %>%
              str_replace(pattern = "\\.R", replacement = str_c("_compCDR", CDR_p_compute_option_temp, "\\.R")) 
            
            folder_file_name_R = paste(folder_name, file_name_R, sep = "/")
            
            write(Rcodetext_child,
                  file = folder_file_name_R)
            
            # Change the sh file
            log_name_R = str_replace(file_name_R, pattern = "\\.R", replacement = "\\.log")
            
            Rscript_line = 
              str_c("Rscript ", file_name_R, 
                    " > ", log_name_R, 
                    " 2>&1")
            
            jobName_line = 
              str_replace(shcodetext[where_jobName], 
                          pattern = "job-name=[[:graph:]]*", 
                          replacement = str_c("job-name=", file_name_R)) %>%
              str_replace(pattern = "\\.R", replacement = "")
            
            shcodetext_child[where_Rscript] = Rscript_line
            shcodetext_child[where_jobName] = jobName_line
            
            file_name_sh = 
              str_replace(file_name_R, pattern = "\\.R", replacement = "\\.sh")
            
            folder_file_name_sh = paste(folder_name, file_name_sh, sep = "/")
            
            write(shcodetext_child, 
                  file = folder_file_name_sh) 
          }
        }
      }
    }
  }
}


