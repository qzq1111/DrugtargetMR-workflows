setwd("F:/DrugTargetMR_Examples/")


# Load the DrugTargetMR package
library(DrugTargetMR)


# WORKFLOW 2: Effector-based MR analysis


################### example 1
# e.g., IL6 → CRP (C-reactive protein) → AIDs 


# bfile is the linkage disequilibrium reference file corresponding to the study population, 
# here it refers to the European population
bfile <- "f:/data_ref/1kg.v3/EUR"

#############  STEP1: IL6 instrumental variables was selected based on gwas summary-level data of CRP
IL6 <- Ensembl_GRCh37[Ensembl_GRCh37$gene_name =="IL6",]

# GCST90029070 was downloaded from GWAS Catalog, and prepared with prepare_others()
extract_target_ivs(gwas_path = "./data_prepare/GCST90029070.txt",
                   pval = 5e-08,
                   target_info = IL6,
                   up_num = 100 * 1000,
                   down_num = 100 * 1000,
                   clump = T,
                   clump_r2 = 0.001,
                   clump_kb = 10000,
                   bfile = bfile,
                   out_path = "./results/work_flow2/")


# exp <- extract_target_ivs(gwas_path = "data_clean/GCST90029070.txt",
#                           rsid = "rs2228145",
#                           out_path = "./result_example2")



#############  STEP2: GWAS
prepare_finngen(file_path = "./data/finngen_R8_T2D_WIDE.gz",
                out_path = "./data_prepare/",
                generate_mr = T,
                generate_smr = T)



#############  STEP3: MR analysis based on clumped ivs
files_ivs <- dir( "./results/work_flow2/",pattern = ".csv$",full.names = T)
mr_local(local_exp_path = files_ivs,
         local_out_path = "./data_prepare/finngen_R8_T2D_WIDE.txt",
         
         clump = T,
         clump_kb = 10000,
         clump_r2 = 0.001,
         bfile = bfile,
         exp_p = 5e-08,
         out_p = 1e-05,
         out_path = "./results/work_flow2/mr/")




################### example 2
# e.g.,targets → LDL-C→ CAD; targets → TG → CAD 

#############  STEP1: targets instrumental variables was selected based on gwas summary-level data of LDL-c or TG
ldl <- Ensembl_GRCh37[Ensembl_GRCh37$gene_name %in% c("LDLR","HMGCR","NPC1L1","PCSK9","APOB","ABCG5","ABCG8",
                                                      "LPL","PPARA","ANGPTL3","APOC3"),]
extract_target_ivs(gwas_path = "data_prepare/ieu-a-300.txt",
                   pval = 5e-08,
                   target_info = ldl,
                   clump = T,
                   clump_r2 = 0.2,
                   clump_kb = 250,
                   bfile = bfile,
                   out_path = "./results/work_flow2/lipid_ldl_c")


tg <- Ensembl_GRCh37[Ensembl_GRCh37$gene_name %in% c("LPL","PPARA","ANGPTL3","APOC3"),]
extract_target_ivs(gwas_path = "data_prepare/ieu-a-302.txt",
                   pval = 5e-08,
                   target_info = tg,
                   clump = T,
                   clump_r2 = 0.2,
                   clump_kb = 250,
                   bfile = bfile,
                   out_path = "./results/work_flow2/lipid_TG")

#############  STEP2: mr analysis
mr_local(local_exp_path = c(dir("./results/work_flow2/lipid_ldl_c/",full.names = T),
                            dir("./results/work_flow2//lipid_TG/",full.names = T)),
         local_out_path = c("./data_prepare/GCST90132315.h.txt","./data_prepare/ebi-a-GCST009722.txt"),
         out_path = "./results/work_flow2/lipid_mr",
         clump = F,
         bfile = bfile,
         exp_p = 5e-08,
         out_p = 1e-5,
         index_snp = F,
         find_proxy = T)

