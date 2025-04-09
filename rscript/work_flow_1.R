setwd("F:/DrugTargetMR_Examples/")

# Load the DrugTargetMR package
library(DrugTargetMR)

# WORKFLOW 1: Multi-omics QTLs based MR and colocalization analysis 

#############  STEP1: Multi-omics QTLs
# e.g..  Preparation of protein cis-pQTLs from DECODE (Ferkingstad, E. et al. Large-scale integration of the plasma proteome with genetics and disease.)
# dir1 is the address where all protein_MRinput.csv is stored

# and the data contains the following columns
#            SNP effect_allele other_allele     eaf   beta       se      pval samplesize chr       pos             Phenotype   gene gene_chr gene_start gene_end cis_trans_hg38      chr_pos
# 1   rs76185902             A            G 0.35813 0.1055 0.008331 9.448e-37      34954  14 106014331 10000_28_CRYBB2_CRBB2 CRYBB2       22   25212564 25231870          trans 14:106014331
# 2 rs1312340709             G            A 0.35941 0.1056 0.008343 1.018e-36      34954  14 106027837 10000_28_CRYBB2_CRBB2 CRYBB2       22   25212564 25231870          trans 14:106027837
# 3   rs35132310             T            A 0.35843 0.1055 0.008336 1.031e-36      34954  14 106009931 10000_28_CRYBB2_CRBB2 CRYBB2       22   25212564 25231870          trans 14:106009931
# 4   rs34118325             G            A 0.35838 0.1055 0.008336 1.045e-36      34954  14 106009601 10000_28_CRYBB2_CRBB2 CRYBB2       22   25212564 25231870          trans 14:106009601
# 5   rs12892841             T            C 0.35834 0.1055 0.008337 1.049e-36      34954  14 106009594 10000_28_CRYBB2_CRBB2 CRYBB2       22   25212564 25231870          trans 14:106009594
# 6   rs12891681             C            T 0.35837 0.1054 0.008332 1.108e-36      34954  14 106009651 10000_28_CRYBB2_CRBB2 CRYBB2       22   25212564 25231870          trans 14:106009651


dir1 <- "d:/MR/data/deCODE_2021/deCODE-2021-MRinput-clean/"
all_files <- dir(path = dir1,full.names = T)

# bfile is the linkage disequilibrium reference file corresponding to the study population, 
# here it refers to the European population
bfile <- "f:/data_ref/1kg.v3/EUR"


# filter and clump
mr_clump(local_exp_path = all_files,

         exp_p = 5e-8,
         maf_filter =T,
         maf = 0.01,
         Fvalue_filter = T,
         Fvalue = 10,
         cis = T ,
         cis_trans_col = "cis_trans_hg38",
         pos_col = "pos",
         clump = T,
         clump_kb = 10000,
         clump_r2 = 0.001,
         bfile = bfile,
         out_path = "./data/deCODE_2021/",
         out_prefix = "decode2021")


#############  STEP2: GWAS
prepare_finngen(file_path = "./data/finngen_R8_T2D_WIDE.gz",
                out_path = "./data_prepare/",
                generate_mr = T,
                generate_smr = T)

#############  STEP3: MR analysis based on clumped ivs
mr_local_clumped(local_exp_path = "./data/deCODE_2021/decode2021-5e-08-clumped-10000kb-0.001r2-20240307.csv",
                 local_out_path = "./data_prepare/finngen_R8_T2D_WIDE.txt",
                 index_snp =  F,
                 out_p = 1e-5,
                 find_proxy = T ,
                 ld_r2 = 0.8,
                 bfile = bfile,
                 out_path = "./results/work_flow1/mr/")


#############  STEP4: coloc analysis
## Select proteins that are significant in MR analysis for colocalization analysis
mr_res <- read.csv("./results/work_flow1/mr/table_s2_mr_results.csv",row.names = 1)
mr_res <- mr_res[mr_res$method %in% c("Inverse variance weighted","Wald ratio"),]
Bonferroni_p <- 0.05/length(unique(mr_res$exposure))
mr_res_sig <- mr_res[mr_res$pval_mr < Bonferroni_p,]


# The function of coloc_batch_decode is to use the  protein GWAS summary RAWdata from DECODE
# to perform colocalization analysis.
?coloc_batch_decode
gwas_path2 <- dir("./decode_rawdata/",full.names = T)
coloc_batch_decode(gwas_path1 = "./data_prepare/finngen_R8_T2D_WIDE.txt",
                   gwas_annotation = gwas_annotation_decode,
                   gwas_path2 = gwas_path2,
                   prepare_method = 2,
                   type1 = "quant",
                   SS1 = 35559,
                   NC1 = 0,
                   type2 = "cc",
                   SS2 = 0,
                   NC2 = 0,
                   bfile = bfile,
                   coloc_plot = T,
                   coloc_plot_pph4 = 0.6,
                   out_path = "./results/work_flow1/coloc/")
