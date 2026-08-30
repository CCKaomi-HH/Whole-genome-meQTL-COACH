#' -----------------------------------------------------------------------------
#' Colocalization analysis
#' -----------------------------------------------------------------------------
suppressMessages(library(tibble))
suppressMessages(library(dplyr))
suppressMessages(library(qvalue))
suppressMessages(library(data.table))
suppressMessages(library(stringr))
suppressMessages(library(coloc))
suppressMessages(library(foreach))
suppressMessages(library(parallel))
suppressMessages(library(ggplot2))
suppressMessages(library(tidyr))
suppressMessages(library(reshape2))
suppressMessages(library(argparser, quietly = TRUE))
parser <- arg_parser("Colocalization analysis for per gene/CpG")
parser <- add_argument(parser, "--GWAS", help = "Path of GWAS summary file")
parser <- add_argument(parser, "--GWAS_name", help = "Name of GWAS cohort")
parser <- add_argument(parser, "--case_ratio", help = "The proportion of cases among all the samples in GWAS")
parser <- add_argument(parser, "--dataset", help = "A type of MatrixEQTL cis file (eQTL/meQTL)")
parser <- add_argument(parser, "--dataset_clump", help = "Clumping file path correspond to dataset")
parser <- add_argument(parser, "--dataset_exp", help = "Expression file path correspond to dataset")
parser <- add_argument(parser, "--dataset_name", help = "Name of QTL (eQTL/meQTL)")
parser <- add_argument(parser, "--dataset_N", help = "The sample size of dataset")
parser <- add_argument(parser, "--cpu", help = "The number of cpu")
parser <- add_argument(parser, "--prefix", help = "Prefix for output file names")
parser <- add_argument(parser, "--output_dir", short = "-o", help = "output directory", default = ".")
argv <- parse_args(parser)

# ------------------------------------------------------------------------------
############    GWAS summary file
# ------------------------------------------------------------------------------
if(argv$GWAS_name=="GCST90129505"){
  GWAS = as.data.frame(fread(argv$GWAS, header=T, sep="\t",data.table=FALSE,nThread=50))
  GWAS = GWAS[,c("SNP","b","se","p")]
  colnames(GWAS) = c("SNP","beta","se","P")
  GWAS$N = 254791
  GWAS = GWAS[,c("SNP","P","beta","se","N")]
  colnames(GWAS) = c("SNP",str_c(c("P","beta","se","N"),".GWAS"))
  print("dim GCST90129505: ")
  print(dim(GWAS))
}

if(argv$GWAS_name=="TPMI_BBJ_meta"){
  GWAS = as.data.frame(fread(argv$GWAS, header=T, sep="\t",data.table=FALSE,nThread=50))
  GWAS = GWAS[,c("MarkerName","Effect","StdErr","P-value")]
  colnames(GWAS) = c("SNP","beta","se","P")
  GWAS$N = 478192
  GWAS = GWAS[,c("SNP","P","beta","se","N")]
  colnames(GWAS) = c("SNP",str_c(c("P","beta","se","N"),".GWAS"))
  print("dim TPMI_BBJ_meta: ")
  print(dim(GWAS))
}

# ------------------------------------------------------------------------------
############    QTL summary file
# ------------------------------------------------------------------------------
if (grepl("meQTL", argv$dataset_name)) {
  # dataset: cis-meQTL file (Note: For simplicity, 'FDR' denotes P-values in meQTL analysis.)
  dataset <- as.data.frame(fread(argv$dataset, sep = "\t", header = TRUE, data.table = FALSE, nThread = 20))
  colnames(dataset) <- c("SNP", "gene", "beta", "t-stat", "FDR")
  tmp <- dataset[, c("SNP", "gene", "FDR")]
  tmp <- merge(tmp, GWAS[, c("SNP", "P.GWAS")], by = "SNP", all.x = TRUE)
  tmp <- tmp[!is.na(tmp$`P.GWAS`), ]
  tmp <- tmp[tmp$P.GWAS < 1e-5, ]
  print("The dim of tmp in meQTL:")
  print(dim(tmp))
  # dataset_clump
  dataset_clump <- as.data.frame(fread(cmd = str_c("cat ", argv$dataset_clump, "|cut -f 1,2,4,11"), sep = "\t", header = TRUE, data.table = FALSE, nThread = 50))
  dataset_clump <- dataset_clump[dataset_clump$P < 2e-12, ]
  dataset_clump$SP2 <- ifelse(dataset_clump$SP2 != "NONE", str_c(dataset_clump$Lead_SNP, ",", dataset_clump$SP2), dataset_clump$Lead_SNP)
  dataset_clump_longer <- dataset_clump %>% separate_longer_delim(SP2, ",")
  colnames(dataset_clump_longer) <- c("Gene", str_c("Lead_SNP.", argv$dataset_name), "leadSNP.FDR", "SNP")
  colnames(dataset_clump) <- c("Gene", "Lead_SNP", "leadSNP.FDR", "SNP")
  dataset_clump <- dataset_clump[dataset_clump$Gene %in% unique(tmp$gene), ]
  rownames(dataset_clump) <- paste0(dataset_clump$Gene, "_", dataset_clump$Lead_SNP)
  dataset_exp <- NULL
  tmp <- NULL
} else {
  # dataset: cis-eQTL file
  dataset <- as.data.frame(fread(argv$dataset, sep = "\t", header = TRUE, data.table = FALSE, nThread = 50))
  tmp <- dataset[, c("SNP", "gene", "FDR")]
  tmp <- merge(tmp, GWAS[, c("SNP", "P.GWAS")], by = "SNP", all.x = TRUE)
  tmp <- tmp[!is.na(tmp$`P.GWAS`), ]
  tmp <- tmp[tmp$P.GWAS < 1e-5, ]
  print("The dim of tmp in eQTL:")
  print(dim(tmp))
  # dataset_clump
  dataset_clump <- as.data.frame(fread(cmd = str_c("cat ", argv$dataset_clump, "|cut -f 1,2,4,11"), sep = "\t", header = TRUE, data.table = FALSE, nThread = 50))
  dataset_clump$SP2 <- ifelse(dataset_clump$SP2 != "NONE", str_c(dataset_clump$Lead_SNP, ",", dataset_clump$SP2), dataset_clump$Lead_SNP)
  dataset_clump_longer <- dataset_clump %>% separate_longer_delim(SP2, ",")
  colnames(dataset_clump_longer) <- c("Gene", str_c("Lead_SNP.", argv$dataset_name), "leadSNP.FDR", "SNP")
  colnames(dataset_clump) <- c("Gene", "Lead_SNP", "leadSNP.FDR", "SNP")
  dataset_clump <- dataset_clump[dataset_clump$Gene %in% unique(tmp$gene), ]
  rownames(dataset_clump) <- paste0(dataset_clump$Gene, "_", dataset_clump$Lead_SNP)
  if (grepl("meQTL", argv$dataset_name)) {
    dataset_exp <- NULL
  } else {
    dataset_exp <- as.data.frame(fread(argv$dataset_exp, sep = "\t", header = TRUE, data.table = FALSE, nThread = 50))
    colnames(dataset_exp) <- replace(colnames(dataset_exp), 1, "gene_id")
    dataset_exp <- dataset_exp %>% column_to_rownames("gene_id")
  }
  tmp <- NULL
}

print("dim dataset: ")
print(dim(dataset))

process_clump <- function(i) {
  coloc.abf.res <- list()
  target <- target_list[i]
  coloc.summary <- data.frame(leadSNP = dataset_clump[rownames(dataset_clump)==target,"Lead_SNP"],
                              gene = dataset_clump[rownames(dataset_clump)==target,"Gene"],
                              leadSNP.FDR =  dataset_clump[rownames(dataset_clump)==target,"leadSNP.FDR"],
                              matrix(NA,
                                     nrow = 1,
                                     ncol = 6),
                              stringsAsFactors = F,
                              row.names = target)
  colnames(coloc.summary) <- c("leadSNP", "gene","leadSNP.FDR","N.target",
                               "H0", str_c("H1.",dataset_name), str_c("H2.",GWAS_name),"H3", "H4")
  gene <- coloc.summary[target, "gene"]
  lead <- coloc.summary[target, "leadSNP"]
  dataset_tmp <- dataset[dataset$gene == gene,]
  GWAS_tmp = GWAS[GWAS$SNP %in% dataset_tmp$SNP,]
  
  QTL <- merge(dataset_tmp, GWAS_tmp, by = "SNP",all.x=T)
  QTL <- QTL[!is.na(QTL$`se.GWAS`), ]
  if(grepl("meQTL",dataset_name)){
    qtl.set <- list(N = as.numeric(dataset_N),
                    beta = QTL$beta,
                    varbeta = (QTL$beta / QTL$`t-stat`)^2,
                    type = "quant",
                    sdY = 1,
                    snp = QTL$SNP)
    GWAS.set <- list(N = as.numeric(QTL$N.GWAS),
                     beta = as.double(QTL$`beta.GWAS`),
                     varbeta =  as.double(QTL$`se.GWAS`)^2,
                     pvalues = QTL$`P.GWAS`,    
                     type = "cc",
                     s = case_ratio,
                     snp = QTL$SNP)
  }else{
    qtl.set <- list(N = as.numeric(dataset_N),
                    beta = QTL$beta,
                    varbeta = (QTL$beta / QTL$`t-stat`)^2,
                    type = "quant",
                    sdY = sd(dataset_exp[gene, ]),
                    snp = QTL$SNP)
    GWAS.set <- list(N = as.numeric(QTL$N.GWAS),
                     beta = as.double(QTL$`beta.GWAS`),
                     varbeta =  as.double(QTL$`se.GWAS`)^2,
                     pvalues = QTL$`P.GWAS`,          
                     type = "cc",
                     s = case_ratio,  
                     snp = QTL$SNP)
  }
  coloc.abf.res[[target]] <- coloc.abf(qtl.set, GWAS.set)
  print("accomplish coloc.abf......")  
  coloc.summary[target, c("N.target","H0", str_c("H1.",dataset_name), str_c("H2.",GWAS_name),"H3", "H4")] <- coloc.abf.res[[target]]$summary
  return(coloc.summary)
}

target_list <- rownames(dataset_clump)
print("length of target_list: ")
print(length(target_list))

# Use parallel processing.
cl <- makeCluster(as.numeric(argv$cpu))
dataset_name <- argv$dataset_name
GWAS_name <- argv$GWAS_name
dataset_N <- argv$dataset_N
case_ratio <- as.numeric(argv$case_ratio)
clusterExport(cl, c("dataset", "dataset_exp", "GWAS", "dataset_clump", "target_list", "case_ratio", "dataset_N", "dataset_name", "GWAS_name"))
clusterEvalQ(cl, {
  library(stringr)
  library(dplyr)
  library(tidyr)
  library(data.table)
  library(coloc)
})
res <- parLapply(cl, 1:length(target_list), process_clump)
stopCluster(cl)

res <- do.call(rbind, res)
data <- as.data.frame(res)
data <- as.data.frame(data)

print("Writing summary file.....")
write.table(data, file = file.path(argv$output_dir, str_c(argv$prefix, "_clump_coloc.summary.txt")), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)