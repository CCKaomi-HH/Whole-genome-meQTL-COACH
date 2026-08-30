#' -----------------------------------------------------------------------------
#' eQTL and meQTL clumping using Plink for per gene and per CpG.
#' -----------------------------------------------------------------------------
suppressMessages(library(qvalue))
suppressMessages(library(data.table))
suppressMessages(library(stringr))
suppressMessages(library(foreach))
suppressMessages(library(parallel))
suppressMessages(library(reshape2))
suppressMessages(library(argparser, quietly = TRUE))
parser <- arg_parser("QTL clumping using Plink for per gene and per CpG")
parser <- add_argument(parser, "--eQTL", help = "eQTL cis file")
parser <- add_argument(parser, "--meQTL", help = "meQTL cis file")
parser <- add_argument(parser, "--fdr", help = "fdr value")
parser <- add_argument(parser, "--bfile", help = "plink file")
parser <- add_argument(parser, "--cpu", help = "the number of cpu")
parser <- add_argument(parser, "--prefix", help = "the prefixion of file name")
parser <- add_argument(parser, "--output_dir", short = "-o", help = "output directory", default = ".")
argv <- parse_args(parser)
# ------------------------------------------------------------------------------
##### Load data
# ------------------------------------------------------------------------------
if (!dir.exists(file.path(argv$output_dir, "/tmp"))) {
    dir.create(file.path(argv$output_dir, "/tmp"))
}
message(str_c("Importing eQTL data....."))
cis_1 <- as.data.frame(fread(argv$cis1, sep = "\t", header = TRUE, data.table = FALSE, nThread = 50))
cis_1 <- cis_1[cis_1$FDR < 0.05, ]
cis_1 <- cis_1[, c("SNP", "gene", "FDR")]
colnames(cis_1) <- c("SNP", "gene", "FDR.qtl_1")

message(str_c("Importing meQTL data....."))
cis_2 <- as.data.frame(fread(argv$cis2, sep = "\t", header = TRUE, data.table = FALSE, nThread = 50))
colnames(cis_2) <- c("SNP", "gene", "beta", "t-stat", "p-value", "chr")
cis_2 <- cis_2[cis_2$`p-value` < 2e-12, ]
cis_2 <- cis_2[, c("SNP", "gene", "p-value")]
colnames(cis_2) <- c("SNP", "CpG", "P.qtl_2")

xQTLs <- merge(cis_1, cis_2, by = "SNP", all.x = TRUE)
xQTLs <- xQTLs[!is.na(xQTLs$CpG), ]
print(dim(xQTLs))
xQTLs <- xQTLs[xQTLs$FDR.qtl_1 < 0.05 & xQTLs$P.qtl_2 < 2e-12, ]
xQTLs$ID <- str_c(xQTLs$gene, "=", xQTLs$CpG)
print(dim(xQTLs))
overlap.gene_cpgs <- unique(xQTLs$ID)
message(str_c("length overlap.gene_cpgs: ", length(overlap.gene_cpgs)))

pvalue <- as.numeric(argv$fdr)
bfile <- argv$bfile
output_dir <- argv$output_dir
prefix <- argv$prefix
data <- NULL
process_gene <- function(i) {
    gene_cpg <- overlap.gene_cpgs[i]
    stats <- xQTLs[xQTLs$ID == gene_cpg, ]
    stats1 <- stats[, c("SNP", "gene", "FDR.qtl_1")]
    stats1$`p-value.qtl_1` <- stats1$FDR.qtl_1
    stats1$FDR2.qtl_1 <- stats1$FDR.qtl_1
    colnames(stats1) <- c("SNP", "gene", "P", "p-value.qtl_1", "FDR.eQTL")

    stats2 <- stats[, c("SNP", "CpG", "P.qtl_2")]
    stats2$`p-value.qtl_2` <- stats2$P.qtl_2
    stats2$FDR.qtl_2 <- stats2$P.qtl_2
    colnames(stats2) <- c("SNP", "CpG", "P", "p-value.qtl_2", "P.meQTL")
    write.table(stats1, file = paste0(output_dir, "/tmp/", gene_cpg, "_", prefix, ".qtl_1.assoc"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
    write.table(stats2, file = paste0(output_dir, "/tmp/", gene_cpg, "_", prefix, ".qtl_2.assoc"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

    tryCatch(
        {
            system(paste0(
                "cd ", output_dir, "/tmp", "\n",
                "/data/work/COACH_QTL_analysis/Tools/plink --bfile ", bfile,
                " --clump ", gene_cpg, "_", prefix, ".qtl_1.assoc", ",", gene_cpg, "_", prefix, ".qtl_2.assoc",
                " --clump-p1 ", pvalue,
                " --clump-p2 1",
                " --clump-r2 0.5",
                " --clump-kb 250",
                " --allow-extra-chr",
                " --out ", gene_cpg, "_", prefix
            ))
            clump <- read.table(file = paste0(output_dir, "/tmp/", gene_cpg, "_", prefix, ".clumped"), stringsAsFactors = FALSE, h = TRUE)
            clump$Gene <- strsplit(gene_cpg, "=")[[1]][1]
            clump$CpG <- strsplit(gene_cpg, "=")[[1]][2]
            clump <- merge(clump, stats1[, c("SNP", "FDR.eQTL")], by = "SNP", all.x = TRUE)
            clump <- merge(clump, stats2[, c("SNP", "P.meQTL")], by = "SNP", all.x = TRUE)
            clump <- clump[, c("Gene", "CpG", "SNP", "BP", "P", "FDR.eQTL", "P.meQTL", "TOTAL", "NSIG", "S05", "S01", "S001", "S0001", "SP2")]
            colnames(clump) <- c("Gene", "CpG", "Lead_SNP", "BP", "P", "FDR.eQTL", "P.meQTL", "TOTAL", "NSIG", "S05", "S01", "S001", "S0001", "SP2")
            return(clump)
        },
        error = function(e) {
            error_data <- data.frame(
                "Gene" = strsplit(gene_cpg, "=")[[1]][1],
                "Gene" = strsplit(gene_cpg, "=")[[1]][2],
                matrix(NA,
                    nrow = 1,
                    ncol = 5 + 7
                ),
                stringsAsFactors = FALSE
            )
            colnames(error_data) <- c("Gene", "CpG", "Lead_SNP", "BP", "P", "FDR.eQTL", "P.meQTL", "TOTAL", "NSIG", "S05", "S01", "S001", "S0001", "SP2")
            return(error_data)
        }
    )
}
# Use parallel processing.
cl <- makeCluster(as.numeric(argv$cpu))
clusterExport(cl, c("xQTLs", "overlap.gene_cpgs", "pvalue", "bfile", "output_dir", "prefix"))
clusterEvalQ(cl, {
    library(stringr)
    library(dplyr)
})
res <- parLapply(cl, 1:length(overlap.gene_cpgs), process_gene)
stopCluster(cl)

res <- do.call(rbind, res)
data <- as.data.frame(res)
data <- as.data.frame(data)
write.table(data, file = file.path(argv$output_dir, str_c(argv$prefix, ".txt")), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
