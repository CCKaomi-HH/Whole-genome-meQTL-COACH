# ------------------------------------------------------------------------------
# Run a cis eQTL association analysis with MatrixEQTL.
#
# The script prepares the genotype, gene-expression and covariate matrices,
# aligns the sample order across all of them, and then runs MatrixEQTL.
# ------------------------------------------------------------------------------

suppressMessages(library(dplyr))
suppressMessages(library(data.table))
suppressMessages(library(stringr))
suppressMessages(library(MatrixEQTL))
suppressMessages(library(argparser))

##### Parse command-line arguments #####
parser <- arg_parser("Run the cis eQTL association analysis with MatrixEQTL.")
parser <- add_argument(parser, "--genotype", help = "Genotype matrix file.")
parser <- add_argument(parser, "--genotype_cvrt", help = "Genotype PCA covariates file (eigenvec).")
parser <- add_argument(parser, "--genotype_cvrt_K", help = "Number of genotype covariates.", default = 5)
parser <- add_argument(parser, "--expression", help = "Gene expression matrix file.")
parser <- add_argument(parser, "--gene_loc", help = "Gene location file.")
parser <- add_argument(parser, "--exp_cvrt", help = "Expression PCA covariates file.")
parser <- add_argument(parser, "--exp_cvrt_K", help = "Number of expression covariates.", default = 5)
parser <- add_argument(parser, "--pvOutputThreshold_cis", help = "Cis p-value threshold.", default = 1e-2)
parser <- add_argument(parser, "--cisDist", help = "Cis distance in bp.", default = 1e6)
parser <- add_argument(parser, "--noFDRsaveMemory", help = "Do not store the FDR table to save memory.", default = FALSE)
parser <- add_argument(parser, "--clinic", help = "Clinical data file.")
parser <- add_argument(parser, "--prefix", help = "Prefix for the output files.")
parser <- add_argument(parser, "--output_dir", short = "-o", help = "Output directory.", default = ".")

argv <- parser$parse_args()

dir.create(argv$output_dir, showWarnings = FALSE, recursive = TRUE)

genotype_covariate_count <- argv$genotype_cvrt_K
expression_covariate_count <- argv$exp_cvrt_K

##### Load input data #####
message("Load genotype data.")
genotype_snp_loc <- as.data.frame(fread(
    argv$genotype, header = TRUE, sep = "\t", data.table = FALSE, nThread = 50
))

message("Load genotype PCA covariates.")
genotype_cvrt <- as.data.frame(fread(
    cmd = str_c("cat ", argv$genotype_cvrt, " | sed '1d'"),
    header = FALSE, sep = " ", data.table = FALSE, nThread = 50
))

message("Load gene expression data (RNA).")
expression <- as.data.frame(fread(
    argv$expression, header = TRUE, sep = "\t", data.table = FALSE, nThread = 50
))

message("Load gene locations.")
gene_location <- as.data.frame(fread(
    argv$gene_loc, header = TRUE, sep = "\t", data.table = FALSE, nThread = 50
))

message("Load expression PCA covariates.")
expression_cvrt <- read.table(
    argv$exp_cvrt, header = FALSE, sep = "\t", row.names = 1,
    stringsAsFactors = FALSE
)

message("Load clinical data.")
clinical <- read.table(
    argv$clinic, header = TRUE, sep = "\t", stringsAsFactors = FALSE
)

##### Prepare the genotype matrix #####
message("Prepare the genotype matrix.")
genotype <- genotype_snp_loc %>%
    select(-Chr, -Pos)
snp_location <- genotype_snp_loc %>%
    select(RSID, Chr, Pos)

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

# additional covariates
clinical <- clinical %>%
    select(Sample, Normal, Age, gender_factor, Stage_I, Stage_II, Stage_III) %>%
    rename(Gender = gender_factor)

clinical[clinical == "UNKNOW"] <- NA

colnames(genotype_cvrt) <- c(
    "Normal",
    paste("geno_cvrt", seq_len(length(colnames(genotype_cvrt)) - 2), sep = "_"),
    "class"
)
genotype_cvrt$class <- NULL
covariates <- merge(clinical, genotype_cvrt[, 1:(genotype_covariate_count + 1)],
                    by = "Normal", all.x = TRUE)

expression_cvrt <- as.data.frame(t(expression_cvrt))
colnames(expression_cvrt) <- c(
    "Sample",
    paste("exp_cvrt", seq_len(length(colnames(expression_cvrt)) - 1), sep = "_")
)
covariates <- merge(covariates, expression_cvrt[, 1:(expression_covariate_count + 1)],
                    by = "Sample", all.y = TRUE)

covariates <- covariates %>%
    select(-Normal)
rownames(covariates) <- covariates$Sample
covariates <- covariates %>%
    select(-Sample)
covariates <- as.data.frame(t(covariates))

##### Align the sample order across all matrices #####
sample_order <- colnames(expression)[2:length(colnames(expression))]
sample_order <- sample_order[sample_order %in% colnames(genotype)]

expression <- expression[, c(colnames(expression)[1], sample_order)]
colnames(expression) <- c("gene_id", sample_order)
rownames(expression) <- expression$gene_id
expression$gene_id <- NULL

genotype <- genotype[, c("RSID", sample_order)]
rownames(genotype) <- genotype$RSID
genotype$RSID <- NULL

covariates <- covariates[, sample_order]
covariates <- cbind(rownames(covariates), covariates)
colnames(covariates) <- replace(colnames(covariates), 1, "covariates")

message("genotype shape:   ", nrow(genotype), " x ", ncol(genotype))
message("expression shape: ", nrow(expression), " x ", ncol(expression))
message("covariates shape: ", nrow(covariates), " x ", ncol(covariates))
message("Sample order matches between genotype and expression: ",
        identical(colnames(genotype), colnames(expression)))

##### Write the covariate matrix that MatrixEQTL reads back from disk #####
covariate_file <- file.path(argv$output_dir, paste0(argv$prefix, ".covariates.merge.txt"))
write.table(covariates, covariate_file,
            quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

##### Build the MatrixEQTL SlicedData objects #####
message("Build the MatrixEQTL input objects.")
genotype_sliced <- SlicedData$new(as.matrix(genotype))
genotype_sliced$ResliceCombined(sliceSize = 10000)

expression_sliced <- SlicedData$new(as.matrix(expression))
expression_sliced$ResliceCombined(sliceSize = 5000)

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

eqtl_result <- Matrix_eQTL_main(
    snps = genotype_sliced,
    gene = expression_sliced,
    cvrt = covariate_sliced,
    pvOutputThreshold = argv$pvOutputThreshold_tra,
    useModel = modelLINEAR,
    errorCovariance = numeric(),
    verbose = TRUE,
    output_file_name.cis = output_file_cis,
    pvOutputThreshold.cis = argv$pvOutputThreshold_cis,
    snpspos = snp_location,
    genepos = gene_location,
    cisDist = argv$cisDist,
    pvalue.hist = "qqplot",
    min.pv.by.genesnp = TRUE,
    noFDRsaveMemory = argv$noFDRsaveMemory
)

saveRDS(eqtl_result, paste0(argv$output_dir, "/", argv$prefix, ".rds"))
message("Finished: results written to ", argv$output_dir)
