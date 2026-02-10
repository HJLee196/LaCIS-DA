rm(list = ls())

library(tidyverse)

setwd("./")

# Read R script 43_44
# mother_file_name = "43_44_mother_code_multiLR_realdata-AllPart12.R"

# Read R script 46_47
mother_file_name = "46_47_mother_code_multidecomp_realdata-AllPart12.R"

Rcodetext <- readLines(mother_file_name)

cell_name_set = c("astrocyte", "oligo", "pericyte", "microglia", "neuron")
PC_set = c("43", "64", "76", "11", "13")

testUse_set = c("LCD") # c("LCD", "wilcox_limma", "LR", "MAST")
llam_set = 0.1 %>% as.character() #c(0, 0.1, 0.2, 0.3)

nSim = 1
down = NULL
seed_num = 2017
# file_path = "path/to/your/rds_folder"
file_path = "./"

max_iter_imp = 100

part1 = FALSE
part2 = TRUE

if (part1) {
  part1_name = "1"
} else {
  part1_name = ""
} 

if (part2) {
  part2_name = "2"
} else {
  part2_name = ""
} 

part_name = paste0("Part", part1_name, part2_name)

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

where_part1 = grep("^part1[[:blank:]]*=[[:blank:]]*",Rcodetext)
where_part2 = grep("^part2[[:blank:]]*=[[:blank:]]*",Rcodetext)

# batch_option # "with" "without" "permute"
# CDR_impute_option # "with" "without" 
# CDR_p_compute_option # "with" "without" 

# where_libpath = grep("\\.libPaths",Rcodetext)
where_readRDS = 
  grep("^HuVascAD = readRDS", Rcodetext)

where_file_path = 
  grep("^file_path", Rcodetext)
where_cell_name =
  grep("^cell_name", Rcodetext)

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

Rcodetext[where_file_path] = 
  paste0("file_path = '",
         file_path,
         "'")

Rcodetext[where_seed_num] =
  Rcodetext[where_seed_num] %>%
  str_replace(pattern = "^seed_num[[:blank:]]*=[[:blank:]]*[[:digit:]]*", 
              replacement = paste0("seed_num = ", seed_num))

Rcodetext[where_max_iter_imp] = 
  Rcodetext[where_max_iter_imp] %>%
  str_replace(pattern = "(^max_iter_imp[[:blank:]]*=[[:blank:]]*)([[:digit:]]*)", 
              replacement = str_c("\\1", max_iter_imp))

Rcodetext[where_part1] = 
  Rcodetext[where_part1] %>% 
  str_replace(pattern = "(^part1[[:blank:]]*=[[:blank:]]*)([[:alnum:]]*)",
              replacement = str_c("\\1", part1))

Rcodetext[where_part2] = 
  Rcodetext[where_part2] %>% 
  str_replace(pattern = "(^part2[[:blank:]]*=[[:blank:]]*)([[:alnum:]]*)",
              replacement = str_c("\\1", part2))

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
              replacement = "mem=256g")

shcodetext[where_time] =
  shcodetext[where_time] %>%
  str_replace(pattern = "[[:digit:]]*:[[:digit:]]*:[[:digit:]]*",
              replacement = "23:00:00")

# Modify the parameters and save it.

 # c("LCD", "wilcox_limma", "LR", "MAST") #"LCD", "wilcox_limma", "LR", "MAST"
if (str_detect(mother_file_name, pattern = "_BH"))
{
  testUse_set = setdiff(x = testUse_set, y = "LCD") # LCD cannot be applied to BH
}

# signal_strength_set = "2" #seq(from = 1, to = 3, length.out = 5) %>% as.character()

batch_option_set = "with"
CDR_impute_option_set = "with" # c("with", "without")
CDR_p_compute_option_set = "with"

for (i_testUse in 1:length(testUse_set))
{
  for (i_llam in 1:length(llam_set))
  {
    for (i_batch_option in 1:length(batch_option_set))
    {
      for (i_CDR_impute_option in 1:length(CDR_impute_option_set))
      {
        for (i_CDR_p_compute_option in 1:length(CDR_p_compute_option_set))
        {
          for (i_PC_cell_name in 1:length(cell_name_set))
          {
            # Initialize
            Rcodetext_child = Rcodetext
            shcodetext_child = shcodetext
            
            folder_name = 
              mother_file_name %>%
              str_to_lower() %>%
              str_replace(pattern = ".r$", replacement = "")
            
            folder_name = str_replace(string = folder_name, pattern = "allpart12", replacement = part_name)
            
            if (!dir.exists(folder_name)){
              dir.create(folder_name)
            }
            
            # Change Parameters for the R script
            
            # testUse
            testUse_temp = testUse_set[i_testUse]
            Rcodetext_child[where_testUse] = str_c("testUse = \"", testUse_temp, "\"")
            
            # if (str_extract(mother_file_name, pattern = "Part[[:digit:]]*") == "Part12"){
            #   print("AllPart12")
            # } else if (i_testUse == 1 & str_extract(mother_file_name, pattern = "-Part[[:digit:]]") == "-Part1"){
            #   testUse_temp = "Part1"
            # } else if (i_testUse != 1 & str_extract(mother_file_name, pattern = "-Part[[:digit:]]") == "-Part1"){
            #   next
            # }
            
            # llam
            llam_temp = llam_set[i_llam]
            Rcodetext_child[where_llam] =
              Rcodetext_child[where_llam] %>% 
              str_replace(pattern = "(=[[:blank:]]*)([[:digit:]\\.]*)", replacement = str_c("\\1", llam_temp))
            
            # signal_strength
            # signal_strength_temp = signal_strength_set[i_signal_strength]
            # Rcodetext_child[where_sign_strength] = 
            #   Rcodetext_child[where_sign_strength] %>%
            #   str_replace(pattern = "(=[[:blank:]]*)([[:digit:]\\.]*)", replacement = str_c("\\1", signal_strength_temp))
            
            # batch_option
            batch_option_temp = batch_option_set[i_batch_option]
            Rcodetext_child[where_batch_option] =
              Rcodetext_child[where_batch_option] %>%
              str_replace(pattern = "(=[[:blank:]]*)([[:alpha:]\"]*)", replacement = str_c("\\1", '\\\"', batch_option_temp, '\\\"'))
            
            # CDR_impute_option
            CDR_impute_option_temp = CDR_impute_option_set[i_CDR_impute_option]
            Rcodetext_child[where_CDR_impute_option] =
              Rcodetext_child[where_CDR_impute_option] %>%
              str_replace(pattern = "(=[[:blank:]]*)([[:alpha:]\"]*)", replacement = str_c("\\1", '\\\"', CDR_impute_option_temp, '\\\"'))
            
            # CDR_p_compute_option
            CDR_p_compute_option_temp = CDR_p_compute_option_set[i_CDR_p_compute_option]
            Rcodetext_child[where_CDR_p_compute_option] =
              Rcodetext_child[where_CDR_p_compute_option] %>%
              str_replace(pattern = "(=[[:blank:]]*)([[:alpha:]\"]*)", replacement = str_c("\\1", '\\\"', CDR_p_compute_option_temp, '\\\"'))
            
            # PC and Cell_name
            cell_name_temp = cell_name_set[i_PC_cell_name]
            PC_temp = PC_set[i_PC_cell_name]
            Rcodetext_child[where_PC] = 
              Rcodetext_child[where_PC] %>%
              str_replace(pattern = "^(PC[[:blank:]]*=[[:blank:]]*)([[:digit:]\\.]*)",
                          replacement = str_c("\\1", PC_temp))
            
            Rcodetext_child[where_cell_name] =
              Rcodetext[where_cell_name] %>%
              str_replace(pattern = "^(cell_name[[:blank:]]*=[[:blank:]]*)([[:alpha:][:punct:]]*)", 
                          replacement = str_c("\\1", '\\\"', cell_name_temp, '\\\"'))
            
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
                  str_replace(pattern = "latent\\.vars[[:blank:]]*=[[:blank:]]*covariate_names",
                              replacement = "latent\\.vars = NULL")
              }
            }
            
            # Save the R file.
            file_name_R =
              str_replace(mother_file_name, pattern = "mother", replacement = "child") %>%
              str_replace(pattern = "\\.R", replacement = str_c("_cell", cell_name_temp, "\\.R")) %>%
              str_replace(pattern = "\\.R", replacement = str_c("_down", ifelse(is.null(down), yes = "NULL", no = down), "\\.R")) %>%
              str_replace(pattern = "\\.R", replacement = str_c("_testUse", testUse_temp, "\\.R")) %>%
              str_replace(pattern = "\\.R", replacement = str_c("_10llam", 10*as.numeric(llam_temp), "\\.R")) %>%
              str_replace(pattern = "\\.R", replacement = str_c("_PC", PC_temp, "\\.R")) 
            # str_replace(pattern = "\\.R", replacement = str_c("_10sgnStrength", 10*as.numeric(signal_strength_temp), "\\.R")) %>%
            # str_replace(pattern = "\\.R", replacement = str_c("_batch", batch_option_temp, "\\.R")) %>%
            # str_replace(pattern = "\\.R", replacement = str_c("_impCDR", CDR_impute_option_temp, "\\.R")) %>%
            # str_replace(pattern = "\\.R", replacement = str_c("_compCDR", CDR_p_compute_option_temp, "\\.R")) 
            
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


