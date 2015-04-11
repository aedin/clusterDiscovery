setwd("~/CloudChen/TCGA-Assembler/")
source("Module_A.r")

# Traverse directories if needed
TraverseAllDirectories(entryPoint = "https://tcga-data.nci.nih.gov/tcgafiles/ftp_auth/distro_ftpusers/anonymous/tumor/blca/",
                       fileLabel = "Traverse_BLCA_Directories_Sept_26_2014")

# 
# =================================================================================
# =                                                                               =
# =                               data download by assembler                      =
# =                                                                               =
# =================================================================================
DownloadClinicalData(traverseResultFile = "./Traverse_BLCA_Directories_Sept_26_2014.rda", 
                     saveFolderName = "~/CloudChen/Projects/MOGSA/Dat/BLCA/Clinical/", 
                     cancerType = "BLCA", 
                     clinicalDataType = c("patient", "drug", "follow_up", "radiation"), 
                     outputFileName = "")


DownloadRNASeqData(traverseResultFile = "./Traverse_BLCA_Directories_Sept_26_2014.rda", 
                   saveFolderName = "~/CloudChen/Projects/MOGSA/Dat/BLCA/RNASeqV2/", 
                   cancerType = "BLCA", 
                   dataType = "rsem.genes.normalized_results",
                   assayPlatform = "RNASeqV2", 
                   # tissueType = , 
                   # inputPatientIDs = , 
                   outputFileName = "")

DownloadCNAData(traverseResultFile = "./Traverse_BLCA_Directories_Sept_26_2014.rda", 
                saveFolderName = "~/CloudChen/Projects/MOGSA/Dat/BLCA/CNV/", 
                cancerType = "BLCA", 
                # assayPlatform = , 
                # tissueType = , 
                # inputPatientIDs = , 
                outputFileName = "")

# =================================================================================
# =                                                                               =
# =                               data prepressing by assembler                   =
# =                                                                               =
# =================================================================================
source("Module_B.r")
setwd("~/CloudChen/TCGA-Assembler/")

cnv <- ProcessCNAData(inputFilePath = "~/CloudChen/Projects/MOGSA/Dat/BLCA/CNV/BLCA__broad.mit.edu__genome_wide_snp_6__nocnv_hg19__2014.txt", 
                      outputFileName = "cnv_snp_6_nosnp_hg19.txt", 
                      outputFileFolder = "~/CloudChen/Projects/MOGSA/Dat/BLCA/Assembly_process/", 
                      refGenomeFile = "./SupportingFiles/Hg19GenePosition.txt")

rnaseqv2 <- ProcessRNASeqData(inputFilePath = "~/CloudChen/Projects/MOGSA/Dat/BLCA/RNASeqV2/BLCA__unc.edu__illuminahiseq_rnaseqv2__rsem.genes.normalized_results__2014.txt", 
                              outputFileName = "rnaseqv2.txt", 
                              outputFileFolder = "~/CloudChen/Projects/MOGSA/Dat/BLCA/Assembly_process/", 
                              dataType = "GeneExp", 
                              verType = "RNASeqV2")



