# ------------------------------------------------------------------------------
# Run a cis meQTL association analysis with MatrixEQTL.
#
# The script prepares the genotype, CpG methylation and covariate matrices,
# aligns the sample order across all of them, and then runs MatrixEQTL.
# ------------------------------------------------------------------------------

suppressMessages(library(dplyr))
suppressMessages(library(data.table))
suppressMessages(library(stringr))
suppressMessages(library(MatrixEQTL))
suppressMessages(library(argparser))

##### Parse command-line arguments #####
parser <- arg_parser("Run the cis meQTL association analysis with MatrixEQTL.")
parser <- add_argument(parser, "--genotype", help = "Genotype matrix file.")
parser <- add_argument(parser, "--genotype_cvrt", help = "Genotype PCA covariates file (eigenvec).")
parser <- add_argument(parser, "--genotype_cvrt_K", help = "Number of genotype covariates.", default = 5)
parser <- add_argument(parser, "--expression", help = "Methylation beta matrix file.")
parser <- add_argument(parser, "--exp_cvrt", help = "Methylation PCA covariates file.")
parser <- add_argument(parser, "--exp_cvrt_K", help = "Number of expression covariates.", default = 5)
parser <- add_argument(parser, "--pvOutputThreshold_cis", help = "Cis p-value threshold.", default = 1e-2)
parser <- add_argument(parser, "--cisDist", help = "Cis distance in bp.", default = 1e6)
parser <- add_argument(parser, "--noFDRsaveMemory", help = "Do not store the FDR table to save memory.", default = FALSE)
parser <- add_argument(parser, "--clinic", help = "Clinical data file.")
parser <- add_argument(parser, "--prefix", help = "Prefix for the output files.")
parser <- add_argument(parser, "--output_dir", short = "-o", help = "Output directory.", default = ".")

argv <- parse_args(parser)

dir.create(argv$output_dir, showWarnings = FALSE, recursive = TRUE)

genotype_covariate_count <- as.numeric(argv$genotype_cvrt_K)
expression_covariate_count <- as.numeric(argv$exp_cvrt_K)

##### Load input data #####
message("Load genotype data.")
genotype_snp_loc <- as.data.frame(fread(
    argv$genotype, header = TRUE, sep = "\t", data.table = FALSE, nThread = 100
))

message("Load genotype PCA covariates.")
genotype_cvrt <- as.data.frame(fread(
    cmd = str_c("cat ", argv$genotype_cvrt, " | sed '1d'"), header = FALSE, sep = " "
))

message("Load CpG methylation data (WGBS).")
methylation <- as.data.frame(fread(
    argv$expression, header = TRUE, sep = "\t", data.table = FALSE, nThread = 100
))

message("Load methylation PCA covariates.")
expression_cvrt <- read.table(
    argv$exp_cvrt, header = TRUE, sep = "\t", row.names = 1,
    check.names = FALSE, stringsAsFactors = FALSE
)

message("Load clinical data.")
clinical <- read.table(
    argv$clinic, header = TRUE, sep = "\t", check.names = FALSE,
    stringsAsFactors = FALSE
)

##### Prepare the genotype matrix #####
message("Prepare the genotype matrix.")
genotype <- genotype_snp_loc %>%
    select(-Chr, -Pos)
rownames(genotype) <- genotype$RSID
genotype$RSID <- NULL
snp_location <- genotype_snp_loc %>%
    select(RSID, Chr, Pos)

##### Prepare the methylation matrix and CpG locations #####
message("Prepare the methylation matrix and CpG locations.")
sample_list <- colnames(methylation)[3:length(colnames(methylation))]
colnames(methylation) <- c("Chr", "Pos", sample_list)
methylation$ID <- str_c(methylation$Chr, ":", methylation$Pos)

cpg_location <- methylation %>%
    select(ID, Chr, Pos)
cpg_location$start <- as.numeric(cpg_location$Pos) - 1
cpg_location$end <- as.numeric(cpg_location$Pos)
cpg_location <- cpg_location[, c("ID", "Chr", "start", "end")]

methylation <- methylation[, c("ID", sample_list)]
rownames(methylation) <- methylation$ID
methylation$ID <- NULL

##### Prepare the covariate matrix #####
message("Prepare the covariate matrix.")
clinical$gender_factor <- case_when(
    clinical$Gender == "Female" ~ 0,
    clinical$Gender == "Male" ~ 1
)

stage_factor <- factor(clinical$Stage)
stage_dummy <- as.data.frame(model.matrix(~ stage_factor - 1))
stage_dummy <- stage_dummy[, c("stage_factorI", "stage_factorII", "stage_factorIII"),
                           drop = FALSE]
colnames(stage_dummy) <- c("Stage_I", "Stage_II", "Stage_III")
clinical <- cbind(clinical, stage_dummy)

# Keep Purity (and optionally CNVcount) as additional covariates.
clinical <- clinical %>%
    select(Sample, Normal, Age, gender_factor, Stage_I, Stage_II, Stage_III) %>%
    rename(Gender = gender_factor)

clinical[clinical == "UNKNOW"] <- NA

# Optional: impute missing covariate values if needed.
# clinical$Purity[is.na(clinical$Purity)] <- mean(clinical$Purity, na.rm = TRUE)
# clinical$CNVcount[is.na(clinical$CNVcount)] <- median(clinical$CNVcount, na.rm = TRUE)

colnames(genotype_cvrt) <- c(
    "Normal",
    paste("geno_cvrt", seq_len(length(colnames(genotype_cvrt)) - 2), sep = "_"),
    "class"
)
genotype_cvrt$class <- NULL
covariates <- merge(clinical, genotype_cvrt[, 1:(genotype_covariate_count + 1)],
                    by = "Normal", all.x = TRUE)

expression_cvrt <- as.data.frame(t(expression_cvrt))
expression_cvrt <- cbind(rownames(expression_cvrt), expression_cvrt)
colnames(expression_cvrt) <- replace(colnames(expression_cvrt), 1, "Sample")
covariates <- merge(covariates, expression_cvrt[, 1:(expression_covariate_count + 1)],
                    by = "Sample", all.x = TRUE)

covariates <- covariates %>%
    select(-Normal)
rownames(covariates) <- covariates$Sample
covariates <- covariates %>%
    select(-Sample)
covariates <- as.data.frame(t(covariates))

##### Align the sample order across all matrices #####
sample_order <- sample_list[sample_list %in% colnames(genotype)]
methylation <- methylation[, sample_order]
genotype <- genotype[, sample_order]
covariates <- covariates[, sample_order]
covariates <- cbind(rownames(covariates), covariates)
colnames(covariates) <- replace(colnames(covariates), 1, "covariates")

message("genotype shape:   ", nrow(genotype), " x ", ncol(genotype))
message("expression shape: ", nrow(methylation), " x ", ncol(methylation))
message("covariates shape: ", nrow(covariates), " x ", ncol(covariates))
message("Sample order matches between genotype and methylation: ",
        identical(colnames(genotype), colnames(methylation)))
message("Sample order matches between methylation and covariates: ",
        identical(colnames(methylation), colnames(covariates)[2:length(colnames(covariates))]))

##### Write the covariate matrix that MatrixEQTL reads back from disk #####
covariate_file <- file.path(argv$output_dir, paste0(argv$prefix, ".covariates.merge.txt"))
if (!file.exists(covariate_file)) {
    write.table(covariates, covariate_file,
                quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
}

##### Build the MatrixEQTL SlicedData objects #####
message("Build the MatrixEQTL input objects.")
genotype_sliced <- SlicedData$new(as.matrix(genotype))
genotype_sliced$ResliceCombined(sliceSize = 50000)

methylation_sliced <- SlicedData$new(as.matrix(methylation))
methylation_sliced$ResliceCombined(sliceSize = 5000)

covariate_sliced <- SlicedData$new()
covariate_sliced$fileDelimiter <- "\t"
covariate_sliced$fileSliceSize <- 5000
covariate_sliced$ResliceCombined()
covariate_sliced$fileSkipRows <- 1
covariate_sliced$fileSkipColumns <- 1
covariate_sliced$LoadFile(covariate_file)

##### Run MatrixEQTL #####
message("Run MatrixEQTL.")
output_file_cis <- paste0(argv$output_dir, "/", argv$prefix, "_cis.txt")

meqtl_result <- Matrix_eQTL_main(
    snps = genotype_sliced,
    gene = methylation_sliced,
    cvrt = covariate_sliced,
    pvOutputThreshold = argv$pvOutputThreshold_tra,
    useModel = modelLINEAR,
    errorCovariance = numeric(),
    verbose = TRUE,
    output_file_name.cis = output_file_cis,
    pvOutputThreshold.cis = argv$pvOutputThreshold_cis,
    snpspos = snp_location,
    genepos = cpg_location,
    cisDist = argv$cisDist,
    pvalue.hist = "qqplot",
    min.pv.by.genesnp = TRUE,
    noFDRsaveMemory = argv$noFDRsaveMemory
)

saveRDS(meqtl_result, paste0(argv$output_dir, "/", argv$prefix, ".rds"))
message("Finished: results written to ", argv$output_dir)
