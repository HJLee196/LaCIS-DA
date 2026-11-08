rm(list = ls())

library(tidyverse)

setwd("./")

# Read R script
# mother_file_name = "01_mother_code_asdp.R"
# mother_file_name = "02_mother_code_LR.R"
# mother_file_name = "03_04_mother_code_multiLR.R"
# mother_file_name = "05_mother_code_decomp.R"
# mother_file_name = "06_07_mother_code_multidecomp.R"
# mother_file_name = "08_mother_code_BH.R"
mother_file_name = "09_mother_code_naive.R"

Rcodetext <- readLines(mother_file_name)

testUse_set = c("LCD") #"LCD", "wilcox_limma", "LR", "MAST"
llam_set = c(0, 0.1, 0.2, 0.3) %>% as.character() # 0 0.1 0.2 0.3
sign_strength_set = c(2) %>% as.character() #seq(from = 1, to = 3, length.out = 5) %>% as.character()

nSim = 20
down = NULL
PC = 30
seed_num = 2017
# file_path = "path/to/combined_superior_parietal_lobe.rds"
file_path = "./combined_superior_parietal_lobe.rds"

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

# where_libpath = grep("\\.libPaths",Rcodetext)
where_readRDS = 
  grep("^superior_parietal_lobe = readRDS", Rcodetext)

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

Rcodetext[where_readRDS] = paste0("superior_parietal_lobe = readRDS(file = '", 
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


if (str_detect(mother_file_name, pattern = "_BH"))
{
  testUse_set = setdiff(x = testUse_set, y = "LCD") # LCD cannot be applied to BH
}

if (str_detect(mother_file_name, pattern = "_BH") | str_detect(mother_file_name, pattern = "_naive"))
{
  llam_set = 0 # We dont need this, neither.
}

for (i_testUse in 1:length(testUse_set))
{
  for (i_llam in 1:length(llam_set))
  {
    for (i_sign_strength in 1:length(sign_strength_set))
    {
      # Initialize
      Rcodetext_child = Rcodetext
      shcodetext_child = shcodetext
      
      folder_name = 
        mother_file_name %>%
        str_to_lower() %>%
        str_replace(pattern = ".r$", replacement = "")
      
      if (!dir.exists(folder_name)){
        dir.create(folder_name)
      }
      
      # Change Parameters for the R script
      
      # testUse
      testUse_temp = testUse_set[i_testUse]
      Rcodetext_child[where_testUse] = str_c("testUse = \"", testUse_temp, "\"")
      
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
        str_replace(pattern = "\\.R", replacement = str_c("_testUse", testUse_temp, "\\.R")) %>%
        str_replace(pattern = "\\.R", replacement = str_c("_10llam", 10*as.numeric(llam_temp), "\\.R")) %>%
        str_replace(pattern = "\\.R", replacement = str_c("_10sgnStrength", 10*as.numeric(sign_strength_temp), "\\.R"))
      
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


