# ------------------------------------------------------------------------------
# Mediation analysis script.
# ------------------------------------------------------------------------------
suppressMessages(library(stringr))
suppressMessages(library(data.table))
suppressMessages(library(dplyr))
suppressMessages(library(parallel))
suppressMessages(library(ggplot2))
suppressMessages(library(purrr))
suppressMessages(library(tidyr))
suppressMessages(library(bmediatR))
suppressMessages(library(mediation))
suppressMessages(library(bda))
suppressMessages(library(fdrtool))

# ==============================================================================
# Section 1. mediation-package bootstrap + bda Sobel mediation
# ==============================================================================

##### Load input data #####
coloc_result_pp4 <- as.data.table(fread(
    "/data/work/COACH_QTL_analysis/mediation/data/eQTL_meQTL_clump_coloc.summary.txt",
    sep = "\t", header = TRUE, nThread = 50
))

expr_mat <- as.matrix(fread("/data/work/COACH_QTL_analysis/mediation/data/gene_expression.txt", header = TRUE))
colnames_expr <- colnames(expr_mat)

methy_mat <- as.matrix(fread("/data/work/COACH_QTL_analysis/mediation/data/CpG_methylation.txt", header = TRUE))
colnames_methy <- colnames(methy_mat)

geno_mat <- as.matrix(fread("/data/work/COACH_QTL_analysis/mediation/data/SNP_genotype.txt", header = TRUE))
colnames_geno <- colnames(geno_mat)

cov_mat <- as.matrix(fread("/data/work/COACH_QTL_analysis/mediation/data/covariates.txt", header = TRUE))
colnames_cov <- colnames(cov_mat)

stopifnot(nrow(expr_mat) == nrow(methy_mat),
          nrow(expr_mat) == nrow(geno_mat),
          nrow(expr_mat) == nrow(cov_mat))

message("Data loaded:")
message("  samples: ", nrow(expr_mat))
message("  genes: ", ncol(expr_mat))
message("  CpGs: ", ncol(methy_mat))
message("  SNPs: ", ncol(geno_mat))
message("  colocalization results: ", nrow(coloc_result_pp4))

##### Per-trio mediation function (bootstrap + Sobel) #####
run_mediation_analysis <- function(i) {
    snp_id <- coloc_result_pp4$leadSNP[i]
    gene_id <- coloc_result_pp4$gene[i]
    cpg_id <- coloc_result_pp4$CpG[i]

    expr <- expr_mat[, gene_id]
    methy <- methy_mat[, cpg_id]
    geno <- geno_mat[, snp_id]
    covs <- cov_mat

    datt <- data.frame(expr = expr, methy = methy, geno = geno)
    datt <- cbind(datt, covs)
    cov_formula <- paste(colnames_cov, collapse = " + ")

    result <- list(
        sme_prop = NA, sme_pval = NA, sme_ci_low = NA, sme_ci_high = NA,
        sem_prop = NA, sem_pval = NA, sem_ci_low = NA, sem_ci_high = NA,
        sobel_p_sme = NA, sobel_p_sem = NA
    )

    # SME model: SNP -> methylation -> expression (bootstrap).
    tryCatch({
        b_sme <- lm(as.formula(paste("methy ~ geno +", cov_formula)), data = datt)
        c_sme <- lm(as.formula(paste("expr ~ methy + geno +", cov_formula)), data = datt)
        med_sme <- mediate(b_sme, c_sme, sims = 1000, treat = "geno", mediator = "methy")
        summary_result <- summary(med_sme)
        result$sme_prop <- summary_result$n1
        result$sme_pval <- summary_result$n1.p
        result$sme_ci_low <- summary_result$n1.ci[1]
        result$sme_ci_high <- summary_result$n1.ci[2]
    }, error = function(e) {
        message(sprintf("Bootstrap SME failed [%d]: %s", i, e$message))
    })

    # SEM model: SNP -> expression -> methylation (bootstrap).
    tryCatch({
        b_sem <- lm(as.formula(paste("expr ~ geno +", cov_formula)), data = datt)
        c_sem <- lm(as.formula(paste("methy ~ expr + geno +", cov_formula)), data = datt)
        med_sem <- mediate(b_sem, c_sem, sims = 1000, treat = "geno", mediator = "expr")
        summary_result <- summary(med_sem)
        result$sem_prop <- summary_result$n1
        result$sem_pval <- summary_result$n1.p
        result$sem_ci_low <- summary_result$n1.ci[1]
        result$sem_ci_high <- summary_result$n1.ci[2]
    }, error = function(e) {
        message(sprintf("Bootstrap SEM failed [%d]: %s", i, e$message))
    })

    # Sobel test, SME model.
    tryCatch({
        m1_sme <- lm(as.formula(paste("methy ~ geno +", cov_formula)), data = datt)
        m2_sme <- lm(as.formula(paste("expr ~ methy + geno +", cov_formula)), data = datt)
        res_sme <- mediation.test(m1_sme$residuals + m1_sme$fitted.values,
                                  geno,
                                  m2_sme$residuals + m2_sme$fitted.values)
        p_sme <- res_sme[2, 1]
        result$sobel_p_sme <- min(max(p_sme, 0), 1)
    }, error = function(e) {
        message(sprintf("Sobel SME failed [%d]: %s", i, e$message))
    })

    # Sobel test, SEM model.
    tryCatch({
        m1_sem <- lm(as.formula(paste("expr ~ geno +", cov_formula)), data = datt)
        m2_sem <- lm(as.formula(paste("methy ~ expr + geno +", cov_formula)), data = datt)
        res_sem <- mediation.test(m1_sem$residuals + m1_sem$fitted.values,
                                  geno,
                                  m2_sem$residuals + m2_sem$fitted.values)
        p_sem <- res_sem[2, 1]
        result$sobel_p_sem <- min(max(p_sem, 0), 1)
    }, error = function(e) {
        message(sprintf("Sobel SEM failed [%d]: %s", i, e$message))
    })

    return(result)
}

##### Parallel computation (chunked) #####
message("Starting parallel computation...")
cluster <- makeCluster(20)
clusterExport(cluster, c("coloc_result_pp4", "expr_mat", "methy_mat", "geno_mat",
                         "cov_mat", "colnames_cov", "run_mediation_analysis"))
clusterEvalQ(cluster, {
    library(mediation)
    library(bda)
    library(data.table)
})

chunk_size <- 1000
n_tests <- nrow(coloc_result_pp4)
chunks <- split(seq_len(n_tests), ceiling(seq_along(seq_len(n_tests)) / chunk_size))

all_results <- list()
for (chunk in chunks) {
    chunk_results <- parLapply(cluster, chunk, run_mediation_analysis)
    all_results <- c(all_results, chunk_results)
    message("Processed ", max(chunk), " / ", n_tests, " tests")
}
stopCluster(cluster)
message("Parallel computation complete.")

##### Format results #####
results_list <- list(
    sme_prop = sapply(all_results, `[[`, "sme_prop"),
    sme_ci_low = sapply(all_results, `[[`, "sme_ci_low"),
    sme_ci_high = sapply(all_results, `[[`, "sme_ci_high"),
    sme_pval = sapply(all_results, `[[`, "sme_pval"),
    sem_prop = sapply(all_results, `[[`, "sem_prop"),
    sem_ci_low = sapply(all_results, `[[`, "sem_ci_low"),
    sem_ci_high = sapply(all_results, `[[`, "sem_ci_high"),
    sem_pval = sapply(all_results, `[[`, "sem_pval"),
    sobel_p_sme = sapply(all_results, `[[`, "sobel_p_sme"),
    sobel_p_sem = sapply(all_results, `[[`, "sobel_p_sem")
)
results_dt <- as.data.table(results_list)

##### FDR correction for the Sobel p-values #####
safe_fdrtool <- function(pvals) {
    if (!is.vector(pvals)) {
        warning("Input is not a vector; coercing.")
        pvals <- as.vector(pvals)
    }

    pvals_clean <- na.omit(pvals)
    pvals_clean <- pvals_clean[pvals_clean >= 0 & pvals_clean <= 1]

    if (length(pvals_clean) == 0) {
        warning("No valid p-values available for FDR correction.")
        return(rep(NA, length(pvals)))
    }

    if (length(pvals_clean) < 3) {
        warning(sprintf("Too few valid p-values (%d); skipping FDR correction.", length(pvals_clean)))
        return(rep(NA, length(pvals)))
    }

    tryCatch({
        fdr_result <- fdrtool(pvals_clean, statistic = "pvalue", plot = FALSE)
        qvals <- rep(NA, length(pvals))
        qvals[pvals %in% pvals_clean] <- fdr_result$qval
        return(qvals)
    }, error = function(e) {
        warning("FDR correction failed: ", e$message)
        return(rep(NA, length(pvals)))
    })
}

results_dt[, `:=`(
    sobel_q_sme = safe_fdrtool(sobel_p_sme),
    sobel_q_sem = safe_fdrtool(sobel_p_sem)
)]

##### Merge results back into the colocalization table #####
coloc_result_pp4[, `:=`(
    Prop_each_SME = results_dt$sme_prop,
    SME_ci_low = results_dt$sme_ci_low,
    SME_ci_high = results_dt$sme_ci_high,
    Prop_each_SME_pval = results_dt$sme_pval,
    Prop_each_SEM = results_dt$sem_prop,
    SEM_ci_low = results_dt$sem_ci_low,
    SEM_ci_high = results_dt$sem_ci_high,
    Prop_each_SEM_pval = results_dt$sem_pval,
    sobel_p_SME = results_dt$sobel_p_sme,
    sobel_q_SME = results_dt$sobel_q_sme,
    sobel_p_SEM = results_dt$sobel_p_sem,
    sobel_q_SEM = results_dt$sobel_q_sem
)]

##### Write results #####
output_file <- "/data/work/COACH_QTL_analysis/mediation/result/mediation_result.txt"
fwrite(coloc_result_pp4, output_file, sep = "\t", row.names = FALSE, quote = FALSE, nThread = 50)
message("\nAnalysis complete.")

# ==============================================================================
# Section 2. bmediatR mediation (SME: SNP -> methylation -> expression)
# ==============================================================================

##### Load input data #####
coloc_data <- as.data.table(fread(
    "/data/work/COACH_QTL_analysis/colocalization/eQTL_meQTL_coloc/result/eQTL_meQTL_clump_coloc.summary.txt",
    sep = "\t", header = TRUE, nThread = 5))

sample_names <- scan(
    "/data/work/COACH_QTL_analysis/mediation/data/head.txt",
    what = character(), quiet = TRUE)

check_dimensions <- function(mat, expected_rows) {
    if (nrow(mat) != expected_rows) {
        stop("Matrix has ", nrow(mat), " rows, but ", expected_rows, " samples are expected.")
    }
}

n_samples <- length(sample_names)
check_dimensions(geno_mat, n_samples)
check_dimensions(methy_mat, n_samples)
check_dimensions(expr_mat, n_samples)
check_dimensions(cov_mat, n_samples)

add_rownames <- function(mat, names) {
    rownames(mat) <- names
    mat
}

geno_mat <- add_rownames(geno_mat, sample_names)
methy_mat <- add_rownames(methy_mat, sample_names)
expr_mat <- add_rownames(expr_mat, sample_names)
cov_mat <- add_rownames(cov_mat, sample_names)

run_bmediatr <- function(i) {
    snp_id <- coloc_data$leadSNP[i]
    cpg_id <- coloc_data$CpG[i]
    gene_id <- coloc_data$gene[i]

    x_all <- geno_mat[, snp_id]
    valid_samples <- !is.na(x_all)
    valid_sample_names <- names(x_all)[valid_samples]

    if (sum(valid_samples) < 10) {
        return(list(index = i,
                    SNP = snp_id,
                    CpG = cpg_id,
                    Gene = gene_id,
                    Error = "Insufficient valid samples"))
    }

    x <- as.matrix(x_all[valid_samples])
    m <- as.matrix(methy_mat[valid_sample_names, cpg_id])
    y <- as.matrix(expr_mat[valid_sample_names, gene_id])
    z <- as.matrix(cov_mat[valid_sample_names, ])

    tryCatch(
        {
            result <- bmediatR(y = y, M = m, X = x, Z = z,
                               ln_prior_c = "complete",
                               options_X = list(sum_to_zero = FALSE, center = FALSE, scale = FALSE))

            list(index = i,
                 SNP = snp_id,
                 CpG = cpg_id,
                 Gene = gene_id,
                 N_samples = sum(valid_samples),
                 ln_post_c = result$ln_post_c,
                 ln_post_odds = result$ln_post_odds,
                 ln_prior_c = result$ln_prior_c,
                 ln_prior_odds = result$ln_prior_odds,
                 ln_ml = result$ln_ml,
                 Error = NULL)
        },
        error = function(e) {
            list(index = i,
                 SNP = snp_id,
                 CpG = cpg_id,
                 Gene = gene_id,
                 Error = e$message)
        }
    )
}

cluster <- makeCluster(20)
clusterExport(cluster, varlist = c("coloc_data",
                                   "geno_mat",
                                   "methy_mat",
                                   "expr_mat",
                                   "cov_mat",
                                   "run_bmediatr"))
clusterEvalQ(cluster, {
    library(bmediatR)
    library(data.table)
})

bmediatr_results <- parLapply(cluster, seq_len(nrow(coloc_data)), function(i) {
    tryCatch(
        run_bmediatr(i),
        error = function(e) {
            list(index = i,
                 SNP = coloc_data$leadSNP[i],
                 CpG = coloc_data$CpG[i],
                 Gene = coloc_data$gene[i],
                 Error = e$message)
        }
    )
})
stopCluster(cluster)

process_bmediatr_results <- function(results) {
    processed <- lapply(results, function(x) {
        result_df <- data.frame(
            Index = x$index,
            SNP = x$SNP,
            CpG = x$CpG,
            Gene = x$Gene,
            N_samples = ifelse(is.null(x$N_samples), NA, x$N_samples),
            Error = ifelse(is.null(x$Error), NA, x$Error),
            stringsAsFactors = FALSE
        )

        if (is.null(x$Error)) {
            post_prob <- exp(x$ln_post_c)

            result_df$complete_mediation_prob <- ifelse(is.na(post_prob[1, "1,1,0"]), 0, post_prob[1, "1,1,0"])
            result_df$partial_mediation_prob <- ifelse(is.na(post_prob[1, "1,1,1"]), 0, post_prob[1, "1,1,1"])
            result_df$colocalization_prob <- ifelse(is.na(post_prob[1, "1,0,1"]), 0, post_prob[1, "1,0,1"])
            result_df$partial_mediation_reactive_prob <- ifelse(is.na(post_prob[1, "1,*,1"]), 0, post_prob[1, "1,*,1"])
            result_df$complete_mediation_reactive_prob <- ifelse(is.na(post_prob[1, "0,*,1"]), 0, post_prob[1, "0,*,1"])
            result_df$other_non_mediation_prob <- max(0, 1 - sum(
                result_df$complete_mediation_prob,
                result_df$partial_mediation_prob,
                result_df$colocalization_prob,
                result_df$partial_mediation_reactive_prob,
                result_df$complete_mediation_reactive_prob,
                na.rm = TRUE
            ))
        }

        return(result_df)
    })

    rbindlist(processed, fill = TRUE)
}

bmediatr_final <- process_bmediatr_results(bmediatr_results)

success_count <- sum(is.na(bmediatr_final$Error))
fail_count <- sum(!is.na(bmediatr_final$Error))
message("Completed bmediatR analyses: ", success_count)
message("Failed bmediatR analyses: ", fail_count)

if (fail_count > 0) {
    message("\nIndices and error messages for failed bmediatR analyses:")
    print(bmediatr_final[!is.na(bmediatr_final$Error),
                         c("Index", "SNP", "CpG", "Gene", "Error")])
}

fwrite(bmediatr_final,
       "/data/work/COACH_QTL_analysis/mediation/result/bmediatR_coloc_SME.txt",
       sep = "\t", na = "NA")