library(DrugTargetMR)
library(stringr)

setwd("F:/DrugTargetMR_Examples/")

####### Meta-analysis of multiple cohorts.
files <- c("./data_prepare/ieu-a-7.txt",
           "./data_prepare/finngen_R8_I9_CHD.txt")

postgwas_meta(
  metal_exe_path = "e:/Windows-metal/metal.exe",
  file_path =files,
  out_path = "./results/work_flow3/metal/"
)



####### Cross-trait meta-analysis.
?postgwas_mtag
files <- c("./data_prepare/IBD_ebi-a-GCST004131.txt", "./data_prepare/PBC_34033851-GCST90061440-EFO_1001486.h.txt" )
postgwas_mtag(gwas_files = files,
              maf = 0.01,
              ancestry="EUR", 
              out_path = "./results/work_flow3/mtag"
)

?postgwas_cpassoc
postgwas_cpassoc( files_path = files,
                  traits = c("IBD","PBC"),
                  sample_sizes = c(59957,24510),
                  method = 1,
                  out_path = './results/work_flow4/cpassoc/',
                  out_prefix = "cpassoc")
