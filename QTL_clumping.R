#' -----------------------------------------------------------------------------
#' QTL clumping using Plink for per gene/CpG.
#' -----------------------------------------------------------------------------
suppressMessages(library(ggplot2))
suppressMessages(library(RColorBrewer))
suppressMessages(library(qvalue))
suppressMessages(library(data.table))
suppressMessages(library(stringr))
suppressMessages(library(foreach))
suppressMessages(library(parallel))
suppressMessages(library(reshape2))
suppressMessages(library(argparser, quietly = TRUE))
parser <- arg_parser("QTL clumping using Plink for per gene/CpG")
parser <- add_argument(parser, "--cis", help = "MatrixEQTL cis file")
parser <- add_argument(parser, "--bfile", help = "plink file")
parser <- add_argument(parser, "--fdr", help = "fdr value")
parser <- add_argument(parser, "--cpu", help = "the number of cpu")
parser <- add_argument(parser, "--prefix", help = "the prefixion of file name")
parser <- add_argument(parser, "--output_dir", short = "-o", help = "output directory", default = ".")
argv <- parse_args(parser)
# ------------------------------------------------------------------------------
##### Load data
# ------------------------------------------------------------------------------
pvalue <- as.numeric(argv$fdr)
message("Importing data.....")
if (!grepl("meQTL", argv$prefix)) {
    cis <- as.data.frame(fread(argv$cis, sep = "\t", header = TRUE, data.table = FALSE, nThread = 50))
} else {
    message("input meQTL cis file....")
    cis <- as.data.frame(fread(argv$cis, sep = "\t", header = FALSE, data.table = FALSE, nThread = 70))
    colnames(cis) <- c("SNP", "gene", "beta", "t-stat", "p-value", "chr")
    cis$chr <- NULL
    cis$FDR <- cis$`p-value`
    print(dim(cis))
}
message("Importing data accomplish.....")
if (!c("FDR") %in% colnames(cis) & !grepl("meQTL", argv$prefix)) {
    cis$FDR <- p.adjust(cis$`p-value`, "BH")
}
cis <- cis[cis$FDR < pvalue, ]
print(head(cis, 4))
if (!dir.exists(file.path(argv$output_dir, "/tmp"))) {
    dir.create(file.path(argv$output_dir, "/tmp"))
}
data <- NULL
process_gene <- function(i) {
    gene <- gene_list[i]
    message(paste0("gene ", gene, " (", which(gene_list == gene), "/", length(gene_list), "):"))
    stats <- cis[cis$gene == gene, ]
    stats1 <- stats[, c("SNP", "FDR")]
    colnames(stats1) <- c("SNP", "P")
    if (grepl("[|]", gene)) {
        gene_name <- str_split(gene, "[|]")[[1]][2]
    } else if (grepl(":", gene)) {
        gene_name <- gsub(":", "-", gene)
    } else {
        gene_name <- gene
    }
    tryCatch(
        {
            write.table(stats1, file = paste0(output_dir, "/tmp/", gene_name, "_", prefix, ".assoc"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
            system(paste0(
                "cd ", output_dir, "/tmp", "\n",
                "/data/work/COACH_QTL_analysis/Tools/plink --bfile ", bfile,
                " --clump ", gene_name, "_", prefix, ".assoc",
                " --clump-p1 ", pvalue,
                " --clump-p2 1",
                " --clump-r2 0.5",
                " --clump-kb 250",
                " --allow-extra-chr",
                # " --clump-best --clump-annotate pval.eqtl,FDR.eqtl",
                # " --clump-verbose",
                " --out ", gene_name, "_", prefix
            ))
            clump <- read.table(file = paste0(output_dir, "/tmp/", gene_name, "_", prefix, ".clumped"), stringsAsFactors = FALSE, h = TRUE)
            clump$SP2 <- unlist(sapply(clump$SP2, function(x) {
                x <- gsub("\\(1\\)", "", x)
            }))
            clump$Gene <- gene
            clump <- clump[, c("Gene", "SNP", "BP", "P", "TOTAL", "NSIG", "S05", "S01", "S001", "S0001", "SP2")]
            colnames(clump) <- c("Gene", "Lead_SNP", "BP", "P", "TOTAL", "NSIG", "S05", "S01", "S001", "S0001", "SP2")
            return(clump)
        },
        error = function(e) {
            return(data.frame("Gene" = gene, "Lead_SNP" = NA, "BP" = NA, "P" = NA, "TOTAL" = NA, "NSIG" = NA, "S05" = NA, "S01" = NA, "S001" = NA, "S0001" = NA, "SP2" = NA, check.names = FALSE))
        }
    )
}
# Use parallel processing.
cl <- makeCluster(as.numeric(argv$cpu))
gene_list <- unique(cis$gene)
bfile <- argv$bfile
output_dir <- argv$output_dir
prefix <- argv$prefix
clusterExport(cl, c("cis", "gene_list", "pvalue", "bfile", "output_dir", "prefix"))
clusterEvalQ(cl, {
    library(stringr)
})
res <- parLapply(cl, 1:length(gene_list), process_gene)
stopCluster(cl)

res <- do.call(rbind, res)
data <- as.data.frame(res)
data <- as.data.frame(data)
write.table(data, file = file.path(argv$output_dir, str_c(argv$prefix, "_clump.txt")), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
