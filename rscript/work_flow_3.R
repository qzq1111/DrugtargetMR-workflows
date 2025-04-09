library(DrugTargetMR)
library(stringr)



# Load the DrugTargetMR package
library(DrugTargetMR)

# WORKFLOW 3: Multi-omics-wide association studies
setwd("F:/DrugTargetMR_Examples/")


############ UTMOST ############
############ step1: Performing single tissue utmost analysis ############
# Note 1: gwas_file, provide the corresponding GWAS file in standard format. Use prepare_others() for data preprocessing.
# Note 2: utmost_single_tissue performs single-tissue utmost analysis: this step will organize the result data by tissue and summarize the merged results and gene information.
# Note 3: For function usage, please execute: ?utmost_single_tissue

prepare_finngen(file_path = "./data/finngen_R8_T2D_WIDE.gz",
                out_path = "./data_prepare/",
                generate_mr = T,
                generate_smr = T)

utmost <- utmost_single_tissue(gwas_file = "./data_prepare/finngen_R8_T2D_WIDE.txt",
                               cov_folder  = "F:/data_utmost/GTEX_v8_20230609/covariance_GTEx8_normalized_pruned/",
                               weight_db_folder  = "F:/data_utmost/GTEX_v8_20230609/database_normalized_pruned/",
                               workers = 6,
                               out_path = "./results/work_flow3/utmost/" ,
                               out_prefix =  "utmost")


############ Step 2: Perform cross-organizational utmost analysis ############
# Note 1: utmost_joint_GBJ performs cross-tissue utmost analysis. The input data used in this step is the output file from the single-tissue utmost analysis.


# ?utmost_joint_GBJ
files <- dir(path = "./results/work_flow3/utmost/",pattern = "^utmost_summary",recursive = T,full.names =T)
utmost_joint_GBJ(gene_info_file =files[i],
                 input_file = files[i+1],
                 cov_rds = "F:/data_utmost/GTEX_v8_20230609/covariance_joint_GTEx8_normalized_pruned.rds",
                 weight_db_rds =  "F:/data_utmost/GTEX_v8_20230609/database_normalized_pruned.rds",
                 out_path = "./results/work_flow3/utmost/",
                 out_prefix ="finngen_R8_T2D_WIDE",
                 workers = 6)



############ FUSION ############
## Note 1: Convert the GWAS standard file into a .sumstats file.
## Note 2: The file_sumstat returned by the fusion_prepare_sumstat() function is the path of the output file.
file_sumstat <- fusion_prepare_sumstat(file_path = "./data_prepare/finngen_R8_T2D_WIDE.txt",
                                       out_path ="./data_prepare/" )

?fusion_assoc
fusion_assoc(sumstat_path = "./data_prepare/finngen_R8_T2D_WIDE.sumstats",
             weights = "f:/data_fusion/Mancuso-GTEx多组织基因权重/GTExv8.ALL2/GTExv8.ALL.Whole_Blood/GTExv8.ALL.Whole_Blood.pos",
             weights_dir = "f:/data_fusion/Mancuso-GTEx多组织基因权重/GTExv8.ALL2/GTExv8.ALL.Whole_Blood/",
             ref_ld_chr_prefix = "f:/data_fusion/LDREF_hg19/1000G.EUR.",
             ref_ld_chr_num = 1:22,
             out_path = "./results/work_flow3/fusion/",
             workers = 4)





############ SMR ############
# data clean
smr_prepare_ma(file_path ="data_prepare/finngen_R8_T2D_WIDE.ma",out_path = "./data_prepare/")

# smr_qtl2gwas analysis for Single tissues:
smr_qtl2gwas(smr_exe_path = "e:/smr-1.3.1-win-x86_64/smr-1.3.1-win.exe",
             bfile = bfile,
             gwas_path = "data_prepare/finngen_R8_T2D_WIDE.ma",
             qtls_path = "f:/data_smr/GTEx_V8_cis_eqtl_summary_lite/Whole_Blood.lite",
             
             smr_multi_snp = T,
             smr_peqtl = 5e-08,
             thread_num = 10,
             smr_cis_wind = 1000,  
             
             out_path = "./results/work_flow3/smr/",
             out_prefix = "gtex_blood")



############ GCTA-fastbat ############
postgwas_gcta_fastBAT(
  gcta_exe = "e:/gcta-1.94.1-Win-x86_64/exe/gcta-win-1.94.1.exe",
  file_path = "./data_prepare/finngen_R8_T2D_WIDE.ma",
  bfile = bfile,
  glist = "y:/data_smr/glist-hg19",
  out_path = "./results/work_flow3/GCTA_fastbat/",
  out_prefix = "finngen_R8_T2D_WIDE"
)


############ MAGMA ############
## Coordinate transformation: Convert coordinates from the hg38 to the hg19
transform_hg38ToHg19(file_path = "./data_prepare/finngen_R8_T2D_WIDE.txt",
                     out_path = "./data_prepare/")

## Sample size information: 33,043 cases and 284,971 controls.
append_col_customs(file_path = "./data_prepare/finngen_R8_T2D_WIDE_hg19.txt",
                   add_col_names = "samplesize",
                   add_col_value = 33043+284971,
                   out_path = "./data_prepare/")

## Perform MAGMA analysis.
postgwas_magma(magma_exe = "f:/data_magma/magma.exe",
               gwas_path ="./data_prepare/finngen_R10_G6_MIGRAINE_hg19.txt",
               bfile = "f:/data_magma/ld/g1000_eur/g1000_eur",
               geneloc_file = "f:/data_magma/MAGMA_gene_boundary_files/ENSGv110.coding.genes.txt", # Gene annotation is based on the hg19 coordinates
               geneset_anno = T,
               geneset_file = "f:/data_magma/Gene_set_files/MSigDB_20231Hs_MAGMA.txt",
               genecovar_file = "f:/data_magma/Gene_expression_files/gtex_v8_ts_avg_log2TPM_54tissues.txt",
               out_path = "./results/work_flow3/magma/",
               out_prefix = "magma")



