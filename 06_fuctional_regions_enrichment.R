# ------------------------------------------------------------------------------
# Functional regions Fisher enrichment test for meQTL and eQTL sets.
# ------------------------------------------------------------------------------

# Settings previously exported by the Bash wrapper.
.libPaths(c("/data/work/COACH_QTL_analysis/Tools/xQTLbiolinks", .libPaths()))
Sys.setenv(TMPDIR = "/data/work/COACH_QTL_analysis/enrichment/result/tmp")

suppressWarnings(suppressMessages(library(data.table)))
suppressWarnings(suppressMessages(library(parallel)))

base_fisher_dir <- "/data/work/COACH_QTL_analysis/enrichment/data"

run_fisher_enrichment <- function(qtl_type, n_to_sample, output_dir,
                                  n_sets = 1000, n_cores = 10) {
    sig_snp_file <- file.path(
        base_fisher_dir,
        sprintf("COACH_250s_%s_cis_significant.prune.in.genomic.region.txt", qtl_type)
    )
    nonsig_snp_file <- file.path(
        base_fisher_dir,
        sprintf("COACH_250s_%s_cis_non_significant.prune.in.genomic.region.txt", qtl_type)
    )
    type_label <- sprintf("250s_%s_genetic_regulatory_set", qtl_type)

    ##### Load the full SNP annotation tables once #####
    sig_snps_full <- as.data.frame(fread(sig_snp_file, sep = "\t", header = TRUE,
                                         data.table = FALSE, nThread = 50))
    nonsig_snps_full <- as.data.frame(fread(nonsig_snp_file, sep = "\t", header = TRUE,
                                            data.table = FALSE, nThread = 50))

    annotation_elements <- colnames(sig_snps_full)[-1]

    ##### Fisher enrichment test for one random set #####
    process_set <- function(set_index) {
        sig_set <- sig_snps_full[sample(seq_len(nrow(sig_snps_full)),
                                        min(n_to_sample, nrow(sig_snps_full))), , drop = FALSE]
        nonsig_set <- nonsig_snps_full[sample(seq_len(nrow(nonsig_snps_full)),
                                              min(n_to_sample, nrow(nonsig_snps_full))), , drop = FALSE]

        result_list <- list()
        for (element in annotation_elements) {
            sig_in <- sum(sig_set[[element]] == 1)
            sig_out <- nrow(sig_set) - sig_in
            nonsig_in <- sum(nonsig_set[[element]] == 1)
            nonsig_out <- nrow(nonsig_set) - nonsig_in

            fisher_result <- fisher.test(matrix(c(sig_in, sig_out, nonsig_in, nonsig_out), nrow = 2),
                                         workspace = 1e8)
            result_list[[element]] <- c(type = paste0(type_label, set_index),
                                        element = element,
                                        pvalue = fisher_result$p.value,
                                        or = unname(fisher_result$estimate),
                                        sig_in = sig_in,
                                        sig_out = sig_out,
                                        nonsig_in = nonsig_in,
                                        nonsig_out = nonsig_out)
        }

        result_df <- do.call(rbind, result_list)
        colnames(result_df) <- c("Type", "Element", "Pvalue", "OR",
                                 "Sig_in", "Sig_out", "nonSig_in", "nonSig_out")

        set_dir <- file.path(output_dir, paste0("set_", set_index))
        dir.create(set_dir, showWarnings = FALSE, recursive = TRUE)
        write.table(result_df,
                    file = file.path(set_dir, paste0(type_label, set_index, ".fisher.txt")),
                    sep = "\t", quote = FALSE, row.names = FALSE)
    }

    ##### Run the random sets in parallel #####
    dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
    message("Running ", qtl_type, " Fisher enrichment (", n_to_sample,
            " SNPs per set, ", n_sets, " sets)...")
    invisible(mclapply(seq_len(n_sets), process_set, mc.cores = n_cores, mc.set.seed = TRUE))
    message("Finished ", qtl_type, " Fisher enrichment: ", n_sets,
            " random sets written to ", output_dir, "\n")
}

##### meQTL, then eQTL #####
run_fisher_enrichment(qtl_type = "meQTL",
                      n_to_sample = 100000,
                      output_dir = file.path(base_fisher_dir, "result_meQTL"))

run_fisher_enrichment(qtl_type = "eQTL",
                      n_to_sample = 10000,
                      output_dir = file.path(base_fisher_dir, "result_eQTL"))
