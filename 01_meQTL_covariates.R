# ------------------------------------------------------------------------------
# Run PCAForQTL: estimate principal components of the methylation matrix, choose
# the number of PCs with the elbow method, and export the top PCs.
# ------------------------------------------------------------------------------
suppressWarnings(suppressMessages(library(PCAForQTL)))
suppressWarnings(suppressMessages(library(stringr)))
suppressWarnings(suppressMessages(library(matrixStats)))
suppressWarnings(suppressMessages(library(torch)))
suppressWarnings(suppressMessages(library(qrpca)))
suppressWarnings(suppressMessages(library(data.table)))
suppressWarnings(suppressMessages(library(argparser)))

parser <- arg_parser("Run PCAForQTL PCA estimation.")
parser <- add_argument(parser, "--expr.file", help = "Path of the input molecular phenotype matrix.")
parser <- add_argument(parser, "--prefix", help = "Prefix for the output file names.")
parser <- add_argument(parser, "--output_dir", short = "-o", help = "Output directory.", default = ".")
argv <- parse_args(parser)

runBE <- function(X,
                  prcomp_result = NULL,
                  n_perm = 20,
                  alpha = 0.05,
                  mc_cores = min(n_perm, parallel::detectCores() - 1),
                  verbose = FALSE) {
    if (alpha < 0 || alpha > 1) {
        stop("alpha must be between 0 and 1.")
    }

    n <- nrow(X)  # Number of observations.
    p <- ncol(X)  # Number of features.
    d <- min(n, p)  # Total number of PCs.

    if (verbose) {
        cat("Running PCA on permuted data...\n")
    }

    perm_pve <- lapply(seq_len(n_perm), function(b) {
        if (verbose) {
            cat("Permutation ", b, " out of ", n_perm, "...\n", sep = "")
        }
        # Permute the observations within each feature.
        x_permuted <- apply(X, 2, function(x) sample(x, size = n, replace = FALSE))
        prcomp_perm <- qrpca::qrpca(x_permuted)
        importance_perm <- summary(prcomp_perm)$importance
        return(importance_perm[2, ])
    })
    perm_stats <- matrix(data = unlist(perm_pve), nrow = d, byrow = FALSE)  # PC by permutation.

    if (is.null(prcomp_result)) {
        if (verbose) {
            cat("Running PCA on the unpermuted data...\n")
        }
        prcomp_result <- qrpca::qrpca(X)
    }

    importance <- summary(prcomp_result)$importance
    pve <- importance[2, ]

    # P-value for the j-th PC: proportion of permutations whose PVE >= observed PVE.
    p_values <- (rowSums(perm_stats >= pve) + 1) / (n_perm + 1)

    # Enforce a monotone increase of the p-values.
    for (j in 2:d) {
        if (p_values[j] < p_values[j - 1]) {
            p_values[j] <- p_values[j - 1]
        }
    }

    list(p_values = p_values,
         alpha = alpha,
         num_pcs_chosen = sum(p_values <= alpha))
}

if (!file.exists(argv$output_dir)) {
    dir.create(argv$output_dir, recursive = TRUE)
}

##### Load and prepare the expression/methylation matrix #####
cat("PCAForQTL: loading data ...\n")
expr <- as.data.frame(fread(argv$expr.file, header = TRUE, sep = "\t", data.table = FALSE, nThread = 50))
expr$chrom <- NULL
expr$pos <- NULL
expr$sd <- NULL
sample_id <- colnames(expr)

##### Run PCA and estimate the number of PCs #####
expr <- as.data.frame(t(expr))
prcomp_result <- qrpca(expr)
importance <- summary(prcomp_result)$importance
pve <- importance[2, ]
message("Sum of PVE: ", sum(pve))

# Elbow method for automatic selection of the number of PCs.
result_elbow <- runElbow(prcomp_result = prcomp_result)
message("Elbow method selected K: ", result_elbow)

k_elbow <- result_elbow

makeScreePlot(prcomp_result,
              labels = c("Elbow"),
              values = c(k_elbow),
              titleText = argv$prefix)
ggplot2::ggsave(file.path(argv$output_dir, str_c(argv$prefix, ".PCA.pdf")),
                width = 16, height = 11, unit = "cm")

##### Export the top 200 principal components #####
pc_scores <- prcomp_result
top_pcs <- as.data.frame(pc_scores[, 1:200, drop = FALSE])
pc_names <- colnames(top_pcs)
top_pcs$sample <- rownames(expr)
top_pcs <- top_pcs[, c("sample", pc_names)]
rownames(top_pcs) <- top_pcs$sample
top_pcs$sample <- NULL
top_pcs <- as.data.frame(t(top_pcs))
pc_names <- colnames(top_pcs)
top_pcs$ID <- rownames(top_pcs)
top_pcs <- top_pcs[, c("ID", pc_names)]

write.table(top_pcs,
            file = file.path(argv$output_dir, str_c(argv$prefix, ".PCA.matrix.txt")),
            quote = FALSE,
            sep = "\t",
            row.names = FALSE,
            col.names = TRUE)
