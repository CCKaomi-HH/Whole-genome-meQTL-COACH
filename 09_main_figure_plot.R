#' -----------------------------------------------------------------------------
#' Main Figure - plot Script
#' -----------------------------------------------------------------------------

#######################################################
# Load libraries
#######################################################
suppressMessages(library(data.table))
suppressMessages(library(ggplot2))
suppressMessages(library(dplyr))
suppressMessages(library(tibble))
suppressMessages(library(reshape2))
suppressMessages(library(ggsci))
suppressMessages(library(stringr))
suppressMessages(library(paletteer))
suppressMessages(library(latex2exp))
suppressMessages(library(forcats))
suppressMessages(library(ggthemes))
suppressMessages(library(ggrepel))
suppressMessages(library(lemon))
suppressMessages(library(ggpubr))
suppressMessages(library(ggprism))
suppressMessages(library(aplot))
suppressMessages(library(patchwork))
suppressMessages(library(purrr))
suppressMessages(library(scales))
suppressMessages(library(tidyverse))
suppressMessages(library(ggforce))
suppressMessages(library(Gviz))
suppressMessages(library(tidyr))
suppressMessages(library(GenomicRanges))
suppressMessages(library(rtracklayer))
suppressMessages(library(bedtoolsr))

#######################################################
##### Fig.2B
#######################################################
setwd("/data/work/COACH_QTL_analysis/functional_annotation/result")
genomic_fisher <- as.data.frame(fread("mean_meQTL_SNP_genetic_regulatory.fisher.txt", sep = "\t", header = TRUE,data.table = FALSE, nThread = 50))

y_max <- max(genomic_fisher$OR_up_diff, na.rm = TRUE) + 0.02
y_min <- min(genomic_fisher$OR_down_diff, na.rm = TRUE) - 0.02
y_breaks <- seq(floor(y_min / 0.5) * 0.5, ceiling(y_max / 0.5) * 0.5, by = 0.5)

fold_enrichment_plot <- ggplot(genomic_fisher, aes(x = Element, y = OR_diff)) +
    geom_bar(
        stat = "identity",
        fill = ifelse(genomic_fisher$OR > 1, "#F48FB1FF", "#b7cbe9"),
        width = 0.7,
        position = position_dodge(width = 0.8)
    ) +
    geom_errorbar(
        aes(ymin = OR_down_diff, ymax = OR_up_diff),
        width = 0.15,
        color = "black",
        position = position_dodge(width = 0.8)
    ) +
    geom_hline(yintercept = 0, linetype = "solid", color = "grey") +
    scale_y_continuous(
        name = "Fold enrichment",
        breaks = y_breaks,
        labels = function(x) round(1 + x, 1),
        limits = c(y_min, y_max),
        expand = expansion(mult = c(0, 0))
    ) +
    labs(x = NULL, title = NULL) +
    theme_minimal() +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, colour = "black"),
        axis.text.y = element_text(colour = "black"),
        panel.border = element_rect(colour = "black", fill = NA, size = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks.y = element_line(color = "black"),
        axis.ticks.length = unit(0.15, "cm"),
        axis.ticks.x = element_line(color = "black")
    )

ggsave("/data/work/COACH_QTL_analysis/functional_annotation/plot/Fold_enrichment_of_meQTL_SNP_in_genetic_elements.pdf", plot = fold_enrichment_plot, width = 6, height = 4, dpi = 300, device = "pdf")


#######################################################
##### Fig.2C
#######################################################
setwd("/data/work/COACH_QTL_analysis/meQTL_clumping/data")
meQTL_cis_clumping = as.data.frame(fread("COACH_meQTL_CpG_clumping_clump.txt",sep = "\t",header = T,data.table=FALSE,nThread=10))
#SNP
SNP.path = "/data/work/COACH_QTL_analysis/genotype/250s.snp.vqsr_pass.plink.vcf.alt.txt"
plink.vcf.alt = as.data.frame(fread(SNP.path,sep = "\t",header = T,data.table=FALSE,nThread=10))
colnames(plink.vcf.alt) = c("CHR","POS","SNP","REF","ALT")
plink.vcf.alt$Start = plink.vcf.alt$POS
plink.vcf.alt$CHR = str_c("chr",plink.vcf.alt$CHR)
#chromatin states
chromatin.path = "/data/work/COACH_QTL_analysis/dataset/E075_15_coreMarks_hg38lift_mnemonics.bed.gz"
chromHMM_E075 = as.data.frame(fread(chromatin.path,sep = "\t",header = F,data.table=FALSE,nThread=50))
colnames(chromHMM_E075) = c("CHR","Start","End","states")
tmp = bedtoolsr::bt.intersect(a = plink.vcf.alt[,c("CHR","Start","POS","SNP")], b = chromHMM_E075,wo = T)
#Distribution of chromatin states for meQTLs, grouped by absolute effect sizes
data = tmp[tmp$V4 %in% meQTL_cis_clumping$Lead_SNP,c("V4","V8")]
colnames(data) = c("SNP","status")
data = merge(data,meQTL_cis_significant,by="SNP",all.x=T)

data$Class <- case_when(
  abs(data$beta) >= 0 & abs(data$beta) <= 0.6 ~ "0-0.6",
  abs(data$beta) > 0.6 & abs(data$beta) <= 0.7 ~ "0.6-0.7",
  abs(data$beta) > 0.7 & abs(data$beta) <= 0.8 ~ "0.7-0.8",
  abs(data$beta) > 0.8 & abs(data$beta) <= 0.9 ~ "0.8-0.9",
  abs(data$beta) > 0.9 & abs(data$beta) <= 1 ~ "0.9-1",
  abs(data$beta) > 1 & abs(data$beta) <= 1.1 ~ "1-1.1",
  abs(data$beta) > 1.1 & abs(data$beta) <= 1.2 ~ "1.1-1.2",
  abs(data$beta) > 1.2 ~ ">1.2")
data$status2 <- case_when(
  data$status == "1_TssA" ~ "TssA",
  data$status == "2_TssAFlnk" ~ "TssAFlnk",
  data$status == "3_TxFlnk" ~ "TxFlnk",
  data$status == "4_Tx" ~ "Tx",
  data$status == "5_TxWk" ~ "TxWk",
  data$status == "6_EnhG" ~ "EnhG",
  data$status == "7_Enh" ~ "Enh",
  data$status == "8_ZNF/Rpts" ~ "ZNF/Rpts",
  data$status == "9_Het" ~ "Het",
  data$status == "10_TssBiv" ~ "TssBiv",
  data$status == "11_BivFlnk" ~ "BivFlnk",
  data$status == "12_EnhBiv" ~ "EnhBiv",
  data$status == "13_ReprPC" ~ "ReprPC",
  data$status == "14_ReprPCWk" ~ "ReprPCWk",
  data$status == "15_Quies" ~ "Quies"
)
data$Class = factor(data$Class, levels=c("0-0.6","0.6-0.7","0.7-0.8","0.8-0.9","0.9-1","1-1.1","1.1-1.2",">1.2"))
data$status2 = factor(data$status2, levels=c("Quies","Tx","TxWk","Enh","EnhG","EnhBiv","TssAFlnk","TssA",
                                             "BivFlnk","ReprPCWk","ReprPC","Het","ZNF/Rpts",
                                             "TssBiv","TxFlnk"))
table(data$Class)
Ratio <- data %>% group_by(Class,status2) %>%
  count() %>%
  group_by(Class) %>%
  mutate(Freq = n/sum(n)*100)
Ratio = Ratio[!Ratio$status2 %in% c("TssBiv","TxFlnk"),]

color =c("Quies"="#b0b1b1","Tx"="#95a9d6","TxWk"="#b2e2f7","Enh"="#b0d077","EnhG"="#91a76d","EnhBiv"="#cfc177",
         "TssAFlnk"="#f5a172","TssA"="#f2727d","BivFlnk"="#f7c19e","ReprPCWk"="#98b0c0","ReprPC"="#6d8fa6",
         "Het"="#a57cb4","ZNF/Rpts"="#c28ca5")
ggplot(Ratio, aes(x = Class, y = Freq, fill = status2))+
  geom_col()+
  #geom_text(aes(label = paste(round(Freq, 1),"%")),
  #          position = position_stack(vjust = 0.5))+
  theme_classic()+
  labs(y = "Proportion",x="meQTL beta value range",fill = "")+
  scale_y_continuous(expand=c(0.01,0)) +
  scale_fill_manual(values=color)+
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black", fill = NA),
        axis.ticks = element_line(size=0.6,colour = "black"),
        axis.text.x = element_text(size = 13,angle = 90,vjust=0.5,hjust = 0.5,colour = "black"),
        axis.title.x = element_text(size=16),
        axis.text.y = element_text(size=13,colour = "black"),
        axis.title.y = element_text(size=16),
        legend.key.size = unit(1, 'cm'), #change legend key size
        legend.key.height = unit(0.8, 'cm'), #change legend key height
        legend.key.width = unit(0.8, 'cm'), #change legend key width
        legend.title = element_text(size=10), #change legend title font size
        legend.text = element_text(size=10),
        legend.position = "right")
ggsave("/data/work/COACH_QTL_analysis/meQTL_clumping/plot/Distribution_chromatin_states_meQTLs_grouped_absolute_effect_sizes_hg38.pdf", height = 6, width = 8,dpi = 800)


#######################################################
##### Fig.3B
#######################################################
meqtl_file <- "/data/work/COACH_QTL_analysis/meQTL_mapping/result/COACH_250s_meQTL_cis_2e12_significant.txt"
sig_meqtl <- fread(meqtl_file, sep = "\t", header = TRUE, nThread = 50)

chrom_state_list <- readLines("/data/work/COACH_QTL_analysis/dataset/chromHMM_states_hg38/file_list.txt")
state_tables <- list()

for (chrom_state_file in chrom_state_list) {
    chrom_state <- fread(chrom_state_file, col.names = c("chr", "start", "end", "state"))

    sig_meqtl_bed <- sig_meqtl[, .(chr, start = pos, end = pos)]

    matched <- bedtoolsr::bt.intersect(a = sig_meqtl_bed, b = chrom_state, wa = TRUE, wb = TRUE)
    colnames(matched) <- c("chr", "start", "end", "chr_state", "start_state", "end_state", "state")

    file_id <- sub("_.*", "", basename(chrom_state_file))
    state_counts <- matched %>%
        group_by(state) %>%
        summarise(count = n(), .groups = "drop") %>%
        mutate(proportion = count / nrow(sig_meqtl))

    colnames(state_counts)[2:3] <- paste0(file_id, "_", colnames(state_counts)[2:3])
    state_tables[[file_id]] <- state_counts
}

state_counts_wide <- reduce(state_tables, full_join, by = "state")
state_counts_wide <- state_counts_wide %>%
    mutate(state_num = as.numeric(sub("_.*", "", state))) %>%
    arrange(state_num)
state_counts_wide <- state_counts_wide %>%
    filter(!is.na(state))

# Write the combined chromHMM state count table.
print(state_counts_wide)
fwrite(state_counts_wide,
       "/data/work/COACH_QTL_analysis/dataset/chromHMM_states_hg38/15_chromHMM_states.txt",
       sep = "\t", row.names = FALSE, quote = FALSE, nThread = 50)

state_counts_wide$state2 <- case_when(
    state_counts_wide$state == "1_TssA" ~ "TssA",
    state_counts_wide$state == "2_TssAFlnk" ~ "TssAFlnk",
    state_counts_wide$state == "3_TxFlnk" ~ "TxFlnk",
    state_counts_wide$state == "4_Tx" ~ "Tx",
    state_counts_wide$state == "5_TxWk" ~ "TxWk",
    state_counts_wide$state == "6_EnhG" ~ "EnhG",
    state_counts_wide$state == "7_Enh" ~ "Enh",
    state_counts_wide$state == "8_ZNF/Rpts" ~ "ZNF/Rpts",
    state_counts_wide$state == "9_Het" ~ "Het",
    state_counts_wide$state == "10_TssBiv" ~ "TssBiv",
    state_counts_wide$state == "11_BivFlnk" ~ "BivFlnk",
    state_counts_wide$state == "12_EnhBiv" ~ "EnhBiv",
    state_counts_wide$state == "13_ReprPC" ~ "ReprPC",
    state_counts_wide$state == "14_ReprPCWk" ~ "ReprPCWk",
    state_counts_wide$state == "15_Quies" ~ "Quies"
)
state_counts_wide[is.na(state_counts_wide)] <- 0

state_counts_wide <- state_counts_wide %>%
    mutate(state = factor(state, levels = unique(state[order(as.numeric(sub("_.*", "", state)))])))

state_counts_long <- state_counts_wide %>%
    pivot_longer(cols = -c(state, state_num, state2),
                 names_to = c("file", ".value"),
                 names_sep = "_")

state_counts_long <- state_counts_long %>%
    mutate(point_color = ifelse(file %in% c("E075", "E076", "E101", "E102", "E103", "E106"),
                                "Colorectum", "Other"))

state_counts_long$state2 <- factor(
    state_counts_long$state2,
    levels = unique(state_counts_long$state2[order(state_counts_long$state_num)])
)

state_proportion_plot <- ggplot(state_counts_long, aes(x = state2, y = proportion)) +
    geom_boxplot(fill = "lightblue", color = "black", alpha = 0.3, outlier.shape = NA) +
    geom_jitter(aes(color = point_color), width = 0.2, height = 0, size = 1.5, alpha = 0.7) +
    scale_color_manual(values = c("Colorectum" = "red", "Other" = "lightgray")) +
    labs(x = "Region State", y = "Proportion of CpGs", title = NULL) +
    theme_minimal(base_size = 14) +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1, size = 12, color = "black"),
        axis.text.y = element_text(size = 12, color = "black"),
        axis.title.x = element_text(size = 14, face = "bold", margin = margin(t = 10)),
        axis.title.y = element_text(size = 14, face = "bold", margin = margin(r = 10)),
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        legend.position = c(0.086, 0.917),
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        legend.background = element_rect(fill = "white", color = "black")
    ) +
    scale_y_continuous(expand = c(0, 0),
                       limits = c(0, max(state_counts_long$proportion) * 1.1))

ggsave("/data/work/COACH_QTL_analysis/meQTL_mapping/plot/CpG_Proportion_by_Region_State_Scatter.pdf",
       plot = state_proportion_plot, width = 9, height = 6, device = "pdf")


#######################################################
##### Fig.3D
#######################################################
cpg_sd_file <- "/data/work/COACH_QTL_analysis/molecular_phenotype/methylation/result/250s_T_methy_CG_10x_filtered.matirx.RN.SD.txt"
cpg_sd <- as.data.frame(fread(cpg_sd_file, sep = "\t", header = TRUE,
                              data.table = FALSE, nThread = 50))
meqtl_bed <- cpg_sd[, c("chr", "start", "end", "CpG")]

cpg_value_file <- "/data/work/COACH_QTL_analysis/molecular_phenotype/methylation/result/250s_T_methy_CG_10x.matirx.filter.txt"
cpg_value <- as.data.frame(fread(cpg_value_file, sep = "\t", header = TRUE,
                                 data.table = FALSE, nThread = 50))

add_mean_expression <- function(expression_data) {
    expression_data$mean_expression <- rowMeans(
        expression_data[, 3:(ncol(expression_data) - 1)], na.rm = TRUE
    )
    return(expression_data)
}

cpg_expression <- add_mean_expression(cpg_value)
cpg_expression <- cpg_expression[, c("chrom", "pos", "CpG", "mean_expression")]

cgi <- as.data.frame(fread("/data/work/COACH_QTL_analysis/functional_annotation/data/CpG_Island.bed",
                           sep = "\t", header = FALSE, data.table = FALSE, nThread = 50))
cgi_shore <- as.data.frame(fread("/data/work/COACH_QTL_analysis/functional_annotation/data/CGI_Shore.bed",
                                 sep = "\t", header = FALSE, data.table = FALSE, nThread = 50))
cgi_shelf <- as.data.frame(fread("/data/work/COACH_QTL_analysis/functional_annotation/data/CGI_Shelf.bed",
                                 sep = "\t", header = FALSE, data.table = FALSE, nThread = 50))
open_sea <- as.data.frame(fread("/data/work/COACH_QTL_analysis/functional_annotation/data/Open_Sea.bed",
                                sep = "\t", header = FALSE, data.table = FALSE, nThread = 50))

find_overlap_cpgs <- function(query_bed, subject_bed, bed_name) {
    overlap_result <- bedtoolsr::bt.intersect(a = query_bed, b = subject_bed, wa = TRUE, wb = TRUE)
    message("Number of CpGs overlapping with ", bed_name, ": ", nrow(overlap_result))
    message("First 10 overlapping CpGs with ", bed_name, ":")
    print(head(overlap_result, 10))
    return(overlap_result)
}

cgi_overlaps <- find_overlap_cpgs(meqtl_bed, cgi, "CGI")
cgi_shore_overlaps <- find_overlap_cpgs(meqtl_bed, cgi_shore, "CGI_Shore")
cgi_shelf_overlaps <- find_overlap_cpgs(meqtl_bed, cgi_shelf, "CGI_Shelf")
open_sea_overlaps <- find_overlap_cpgs(meqtl_bed, open_sea, "Open_Sea")

cgi_expression <- cpg_expression %>% filter(CpG %in% cgi_overlaps$V4)
cgi_shore_expression <- cpg_expression %>% filter(CpG %in% cgi_shore_overlaps$V4)
cgi_shelf_expression <- cpg_expression %>% filter(CpG %in% cgi_shelf_overlaps$V4)
open_sea_expression <- cpg_expression %>% filter(CpG %in% open_sea_overlaps$V4)

plot_expression_distribution <- function(expression_data, title, file_name) {
    bins_per_major <- 15
    major_breaks_pct <- seq(0, 100, by = 20)
    total_bins <- length(major_breaks_pct) * bins_per_major
    breaks <- seq(0, 1, length.out = total_bins + 1)

    expression_data$expression_bin <- cut(expression_data$mean_expression,
                                          breaks = breaks,
                                          include.lowest = TRUE)

    plot_data <- as.data.frame(table(expression_data$expression_bin))
    colnames(plot_data) <- c("Methylation", "Frequency")
    plot_data$Frequency_scaled <- plot_data$Frequency / 1000

    n_bins <- nrow(plot_data)
    n_major_breaks <- ceiling(n_bins / bins_per_major)
    major_breaks_pos <- round(seq(1, n_bins, length.out = n_major_breaks))
    major_breaks_pct <- round(seq(0, 100, length.out = n_major_breaks))

    p <- ggplot(plot_data, aes(x = as.numeric(Methylation), y = Frequency_scaled)) +
        geom_bar(stat = "identity",
                 width = 1,
                 color = "black",
                 fill = "skyblue",
                 size = 0.3) +
        labs(x = "Methylation %",
             y = expression("Frequency (10"^3*")"),
             title = title) +
        scale_y_continuous(expand = expansion(mult = c(0.02, 0.05)),
                           limits = c(0, max(plot_data$Frequency_scaled) * 1.05)) +
        scale_x_continuous(
            breaks = major_breaks_pos,
            labels = major_breaks_pct,
            limits = c(0.5, n_bins + 0.5),
            expand = expansion(mult = c(0.01, 0.02))
        ) +
        theme_classic(base_rect_size = 0) +
        theme(
            panel.background = element_blank(),
            axis.ticks = element_line(size = 0.6, colour = "black"),
            axis.text.x = element_text(size = 16, colour = "black"),
            axis.title.x = element_text(size = 18),
            axis.text.y = element_text(size = 16, colour = "black"),
            axis.title.y = element_text(size = 18),
            plot.margin = unit(c(1, 1, 0.5, 0.5), "cm"),
            plot.title = element_text(size = 20, hjust = 0.5)
        )
    ggsave(file_name, p, width = 7, height = 4, dpi = 300)
}

plot_expression_distribution(cgi_expression, "CpG Islands", "CGI_Expression_Distribution.pdf")
plot_expression_distribution(cgi_shore_expression, "CpG Shores", "CGI_Shore_Expression_Distribution.pdf")
plot_expression_distribution(cgi_shelf_expression, "CpG Shelves", "CGI_Shelf_Expression_Distribution.pdf")
plot_expression_distribution(open_sea_expression, "Open Sea", "Open_Sea_Expression_Distribution.pdf")


#######################################################
##### Fig.4A
#######################################################
setwd("/data/work/COACH_QTL_analysis/colocalization/eQTL_meQTL_coloc/data")
data = as.data.frame(fread("Standardized_Fig.4A_data_input.txt", header=T, sep="\t",data.table=FALSE))
colnames(data) = c("gene","gene_list","eQTL_gene_count","meQTL_CpG_count","ALL")
data$eQTL_gene_count = as.numeric(data$eQTL_gene_count)
data$meQTL_CpG_count = as.numeric(data$meQTL_CpG_count)

data$class <- case_when(
  data$meQTL_CpG_count == 1 ~ "1",
  data$meQTL_CpG_count == 2 ~ "2",
  data$meQTL_CpG_count == 3 ~ "3",
  data$meQTL_CpG_count == 4 ~ "4",
  data$meQTL_CpG_count == 5 ~ "5",
  data$meQTL_CpG_count == 6 ~ "6",
  data$meQTL_CpG_count >= 7 & data$meQTL_CpG_count <= 8 ~ "7-8",
  data$meQTL_CpG_count >= 9 & data$meQTL_CpG_count <= 10 ~ "9-10",
  data$meQTL_CpG_count >= 11 & data$meQTL_CpG_count <= 12 ~ "11-12",
  data$meQTL_CpG_count >= 13 & data$meQTL_CpG_count <= 14 ~ "13-14",
  data$meQTL_CpG_count >= 15 & data$meQTL_CpG_count <= 16 ~ "15-16",
  data$meQTL_CpG_count >= 17 & data$meQTL_CpG_count <= 18 ~ "17-18",
  data$meQTL_CpG_count >= 19 & data$meQTL_CpG_count <= 20 ~ "19-20",
  data$meQTL_CpG_count >= 21 & data$meQTL_CpG_count <= 25 ~ "21-25",
  data$meQTL_CpG_count >= 26 & data$meQTL_CpG_count <= 30 ~ "26-30",
  data$meQTL_CpG_count >= 31  ~ ">30")

table(data$class)

inpurt = data %>% group_by(class) %>% summarise(Count = length(class))
colnames(inpurt) = c("class","count")
inpurt$class = factor(inpurt$class, levels=c("1","2","3","4","5","6","7-8","9-10","11-12","13-14","15-16","17-18","19-20","21-25","26-30",">30"))
ggplot(inpurt, aes(y=count, x=class)) +
  geom_bar(stat="identity",width = 0.8,color="#9E80BC",fill="#9E80BC")+
  geom_text(aes(x = class, y = count-12, label = count),
            size = 5,
            fontface = "bold",
            color = 'white')+
  labs(x="Number of meCpG", y="Number of colocalized eGene-meCpG groups") +
  theme_classic(base_rect_size = 0)+
  scale_y_continuous(expand = c(0.01, 0),limit=c(0,650)) +
  scale_x_discrete(expand = expansion(mult = c(0.035, 0.02))) +
  theme(panel.background = element_blank(),
        axis.ticks = element_line(size=0.6,colour = "black"),
        axis.text.x = element_text(size = 10,colour = "black"),
        axis.title.x = element_text(size=16),
        axis.text.y = element_text(size=10,colour = "black"),
        axis.title.y = element_text(size=16),
        plot.margin=unit(c(1,1,0.5,0.5),units=,"cm"),
        legend.direction = "horizontal",legend.position = "top")+
  guides(fill=guide_legend(title = ""))
ggsave("/data/work/COACH_QTL_analysis/colocalization/eQTL_meQTL_coloc/plot/Number_meQTL_eQTL-colocalization_histogram.pdf", height = 6, width = 10,dpi = 800)


#######################################################
##### Fig.4B
#######################################################
beta_08_file <- "/data/work/COACH_QTL_analysis/colocalization/eQTL_meQTL_coloc/result/gene_CpG_beta_ppi_08_new.txt"
beta_08 <- as.data.frame(fread(beta_08_file, sep = "\t", header = TRUE,
                               data.table = FALSE, nThread = 50))

beta_1_file <- "/data/work/COACH_QTL_analysis/colocalization/eQTL_meQTL_coloc/result/gene_CpG_beta_ppi_1_new.txt"
beta_1 <- as.data.frame(fread(beta_1_file, sep = "\t", header = TRUE,
                              data.table = FALSE, nThread = 50))

prepare_direction_data <- function(df, file_type) {
    df %>%
        mutate(
            original_direction = case_when(
                beta_meQTL * beta_eQTL > 0 ~ "same",
                beta_meQTL * beta_eQTL < 0 ~ "different",
                TRUE ~ NA_character_
            ),
            eqtm_direction = case_when(
                beta_eQTM > 0 ~ "same",
                beta_eQTM < 0 ~ "different",
                TRUE ~ NA_character_
            )
        ) %>%
        pivot_longer(
            cols = c(original_direction, eqtm_direction),
            names_to = "analysis_type",
            values_to = "direction"
        ) %>%
        mutate(
            category = case_when(
                analysis_type == "original_direction" & Promoter_region ~ paste0(file_type, "_Promoter"),
                analysis_type == "original_direction" & !Promoter_region ~ paste0(file_type, "_NonPromoter"),
                analysis_type == "eqtm_direction" & Promoter_region ~ paste0(file_type, "_eQTM_Promoter"),
                analysis_type == "eqtm_direction" & !Promoter_region ~ paste0(file_type, "_eQTM_NonPromoter")
            )
        ) %>%
        filter(!is.na(direction))
}

beta_1_data <- prepare_direction_data(beta_1, "All")
beta_08_data <- prepare_direction_data(beta_08, "0.8")
combined_data <- bind_rows(beta_1_data, beta_08_data)

rename_categories <- function(category) {
    case_when(
        category == "All_Promoter" ~ "All \n(Promoter)",
        category == "0.8_Promoter" ~ "0.8 \n(Promoter)",
        category == "All_NonPromoter" ~ "All \n(Non-promoter)",
        category == "0.8_NonPromoter" ~ "0.8 \n(Non-promoter)",
        category == "All_eQTM_Promoter" ~ "All \n(Promoter)",
        category == "0.8_eQTM_Promoter" ~ "0.8 \n(Promoter)",
        category == "All_eQTM_NonPromoter" ~ "All \n(Non-promoter)",
        category == "0.8_eQTM_NonPromoter" ~ "0.8 \n(Non-promoter)",
        TRUE ~ category
    )
}

direction_summary <- combined_data %>%
    mutate(
        display_category = rename_categories(category),
        analysis_group = ifelse(
            grepl("eQTM", category),
            "gene expression and DNAm levels",
            "eQTL and meQTL"
        )
    ) %>%
    group_by(display_category, analysis_group, direction) %>%
    summarise(count = n(), .groups = "drop") %>%
    group_by(display_category, analysis_group) %>%
    mutate(
        total = sum(count),
        proportion = count / total
    ) %>%
    ungroup() %>%
    mutate(
        display_category = factor(display_category, levels = c(
            "All \n(Promoter)", "0.8 \n(Promoter)",
            "All \n(Non-promoter)", "0.8 \n(Non-promoter)"
        )),
        direction = factor(direction, levels = c("same", "different"))
    )

direction_plot <- direction_summary %>%
    filter(analysis_group == "eQTL and meQTL") %>%
    ggplot(aes(x = display_category, y = proportion, fill = direction)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7,
             color = "black", linewidth = 0.3) +
    scale_fill_manual(
        values = c("same" = "#4E79A7", "different" = "#E15759"),
        labels = c("Same direction", "Different direction")
    ) +
    scale_y_continuous(
        labels = function(x) x * 100,
        breaks = seq(0, 1, by = 0.2),
        limits = c(0, max(direction_summary$proportion) * 1.1)
    ) +
    labs(
        x = "PP4 in coloc",
        y = "Percentage (%)",
        fill = NULL,
        title = "eQTL and meQTL"
    ) +
    theme_bw(base_size = 12) +
    theme(
        axis.text.x = element_text(angle = 0, size = 16, color = "black"),
        axis.text.y = element_text(size = 16, color = "black"),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 0.8),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        legend.position = c(0.95, 0.95),
        legend.direction = "horizontal",
        legend.justification = c("right", "top"),
        legend.box.just = "right",
        legend.margin = margin(6, 6, 6, 6),
        legend.background = element_rect(color = NA, fill = "white"),
        legend.text = element_text(size = 16)
    ) +
    geom_text(
        aes(label = round(proportion * 100)),
        position = position_dodge(width = 0.8),
        vjust = -0.5,
        size = 7
    )
ggsave("/data/work/COACH_QTL_analysis/colocalization/eQTL_meQTL_coloc/plot/Histogram_percentage_of_directions_eGene–meCpG_colocalization.pdf", plot = direction_plot, width = 12, height = 6, dpi = 300)


#######################################################
##### Fig.5A
#######################################################
promoter_file <- "/data/work/COACH_QTL_analysis/dataset/SCREEN/GRCh38-PLS.bed"
promoter <- as.data.frame(fread(promoter_file, sep = "\t", header = FALSE,
                                data.table = FALSE, nThread = 50))

proximal_enhancer_file <- "/data/work/COACH_QTL_analysis/dataset/SCREEN/GRCh38-cCREs.pELS.bed"
proximal_enhancer <- as.data.frame(fread(proximal_enhancer_file, sep = "\t", header = FALSE,
                                         data.table = FALSE, nThread = 50))

distal_enhancer_file <- "/data/work/COACH_QTL_analysis/dataset/SCREEN/GRCh38-cCREs.dELS.bed"
distal_enhancer <- as.data.frame(fread(distal_enhancer_file, sep = "\t", header = FALSE,
                                       data.table = FALSE, nThread = 50))

insulator_file <- "/data/work/COACH_QTL_analysis/dataset/SCREEN/GRCh38-cCREs.CTCF-only.bed"
insulator <- as.data.frame(fread(insulator_file, sep = "\t", header = FALSE,
                                 data.table = FALSE, nThread = 50))

gene_body_file <- "/data/work/COACH_QTL_analysis/dataset/Gencode_v32/gencode.v32.gene.body.bed"
gene_body <- as.data.frame(fread(gene_body_file, sep = "\t", header = FALSE,
                                 data.table = FALSE, nThread = 50))

bim_file <- "/data/work/COACH_QTL_analysis/genotype/250s.snp.vqsr_pass.plink.bim"
bim <- as.data.frame(fread(bim_file, sep = "\t", header = FALSE,
                           data.table = FALSE, nThread = 50))

coloc_file <- "/data/work/COACH_QTL_analysis/colocalization/eQTL_meQTL_coloc/result/eQTL_meQTL_clump_coloc.summary.txt"
coloc <- as.data.frame(fread(coloc_file, sep = "\t", header = TRUE,
                             data.table = FALSE, nThread = 50))

h4_threshold <- 0.8
coloc <- coloc[coloc$H4 > h4_threshold, ]
coloc_cpgs <- unique(coloc$CpG)

target_genes <- data.frame(
    gene = c("ABHD12B", "NIN", "PYGL"),
    color = c("#595959", "#595959", "#595959")
)

info_file <- "/data/work/COACH_QTL_analysis/meQTL_clumping/data/PYGL_loci_CpGs_summary.txt"
info <- as.data.frame(fread(info_file, sep = "\t", header = TRUE,data.table = FALSE, nThread = 50))

methylation_file <- "/data/work/COACH_QTL_analysis/molecular_phenotype/methylation/split_chr/250s_methy_CG_10x.matirx.filter.RN_chr14.txt"
methylation_rn <- as.data.frame(fread(methylation_file, sep = "\t", header = TRUE,
                                      data.table = FALSE, nThread = 50))
rownames(methylation_rn) <- methylation_rn$CpG
methylation_rn <- methylation_rn[, -1]

methylation_rn_t <- as.data.frame(t(methylation_rn))
rownames(methylation_rn_t) <- NULL
correlation_rn <- as.data.frame(cor(methylation_rn_t))

snp_correlation_data <- data.frame(
    CpG = rownames(correlation_rn),
    Lead_CpG_correlation = correlation_rn$chr14_50899858
)

snp_correlation_data <- merge(info, snp_correlation_data, by.x = "TargetID", by.y = "CpG", all.x = TRUE)
colnames(snp_correlation_data) <- c("SNP", "CHR", "MAPINFO", "Pval", "Lead_CpG_correlation")
snp_correlation_data$Lead_CpG_correlation <- abs(snp_correlation_data$Lead_CpG_correlation)
print(head(snp_correlation_data, 3))

suppressWarnings(suppressMessages(library(EnsDb.Hsapiens.v86)))
locus_result <- locus(snp_correlation_data,
                      gene = target_genes$gene,
                      pos = "MAPINFO",
                      p = "Pval",
                      LD = "Lead_CpG_correlation",
                      flank = 2.5e5,
                      ens_db = "EnsDb.Hsapiens.v86")

exon_df <- as.data.frame(locus_result$EX)
exon_df$gene_name <- mapIds(EnsDb.Hsapiens.v86,
                            keys = exon_df$gene_id,
                            column = "GENENAME",
                            keytype = "GENEID",
                            multiVals = "first")
exon <- exon_df[exon_df$gene_name %in% target_genes$gene, ]
exon <- exon[, c("seqnames", "start", "end", "strand", "gene_name")]
colnames(exon) <- c("chr", "start", "end", "strand", "gene")
exon$chr <- paste0("chr", exon$chr)

target_snp <- "rs6572708"
target_snp_pos <- 50899857

target_genes_info <- gene_body %>%
    dplyr::filter(gene %in% target_genes$gene) %>%
    left_join(target_genes, by = "gene")
target_genes_info <- target_genes_info %>%
    mutate(gene_id = match(gene, target_genes$gene))
print(target_genes_info)

target_data <- coloc %>% dplyr::filter(gene %in% target_genes$gene)

clump_file <- "/data/work/COACH_QTL_analysis/eQTL_meQTL_clumping/result/250s_eQTL_meQTL_clumping.txt"
clump <- as.data.frame(fread(clump_file, sep = "\t", header = TRUE,
                             data.table = FALSE, nThread = 50))

snp_rows <- clump %>%
    dplyr::filter(if_any(everything(), ~ str_detect(., target_snp)))
snp_cpgs <- unique(snp_rows$CpG)
message("CpGs linked to ", target_snp, ": ", length(snp_cpgs))

overlap_cpgs <- intersect(coloc_cpgs, snp_cpgs)
message("CpGs shared by coloc and the clump file: ", length(overlap_cpgs))

target_data <- target_data %>% dplyr::filter(CpG %in% overlap_cpgs)

target_snp_info <- bim %>%
    dplyr::filter(SNP == target_snp) %>%
    mutate(chr = paste0("chr", chr))

target_cpg_info <- target_data %>%
    mutate(
        cpg_chr = paste0("chr", sub(":.*", "", CpG)),
        cpg_pos = as.numeric(sub(".*:", "", CpG))
    ) %>%
    dplyr::filter(cpg_chr == unique(target_genes_info$chr))

region_chr <- unique(target_genes_info$chr)
region_start <- min(c(target_genes_info$start, target_snp_info$pos, target_cpg_info$cpg_pos)) - 20000
region_end <- max(c(target_genes_info$end, target_snp_info$pos, target_cpg_info$cpg_pos)) + 20000

region_start_line <- min(c(target_genes_info$start, target_snp_info$pos, target_cpg_info$cpg_pos)) - 10000
region_end_line <- max(c(target_genes_info$end, target_snp_info$pos, target_cpg_info$cpg_pos)) + 10000

promoter_in_region <- promoter %>%
    dplyr::filter(chr == region_chr, end >= region_start_line, start <= region_end_line) %>%
    mutate(type = "promoter")

proximal_enhancer_in_region <- proximal_enhancer %>%
    dplyr::filter(chr == region_chr, end >= region_start_line, start <= region_end_line) %>%
    mutate(type = "proximal_enhancer")

distal_enhancer_in_region <- distal_enhancer %>%
    dplyr::filter(chr == region_chr, end >= region_start_line, start <= region_end_line) %>%
    mutate(type = "distal_enhancer")

insulator_in_region <- insulator %>%
    dplyr::filter(chr == region_chr, end >= region_start_line, start <= region_end_line) %>%
    mutate(type = "insulator")

create_gene_line <- function(start, end, ymin, ymax, strand) {
    ycenter <- (ymin + ymax) / 2
    data.frame(
        x = c(start, end),
        y = c(ycenter, ycenter),
        type = "gene_line"
    )
}

create_exons <- function(exon_df, gene_info, ymin, ymax) {
    body_width <- (ymax - ymin) * 0.5
    ycenter <- (ymin + ymax) / 2

    gene_exons <- exon_df %>% dplyr::filter(gene == gene_info$gene)

    if (nrow(gene_exons) > 0) {
        gene_exons %>%
            mutate(
                ymin = ycenter - body_width,
                ymax = ycenter + body_width,
                type = "exon"
            ) %>%
            dplyr::select(xmin = start, xmax = end, ymin, ymax, type)
    } else {
        NULL
    }
}

gene_shapes <- target_genes_info %>%
    group_by(gene) %>%
    do({
        gene_row <- .

        line_data <- create_gene_line(
            gene_row$start,
            gene_row$end,
            0.04 + 0.015 * (gene_row$gene_id - 1),
            0.03 + 0.015 * (gene_row$gene_id - 1)
        )
        line_data$gene <- gene_row$gene
        line_data$gene_id <- gene_row$gene_id
        line_data$color <- gene_row$color

        exon_data <- create_exons(
            exon_df = exon,
            gene_info = gene_row,
            ymin = 0.04 + 0.015 * (gene_row$gene_id - 1),
            ymax = 0.03 + 0.015 * (gene_row$gene_id - 1)
        )

        if (!is.null(exon_data)) {
            exon_data$gene <- gene_row$gene
            exon_data$gene_id <- gene_row$gene_id
            exon_data$color <- gene_row$color
            bind_rows(line_data, exon_data)
        } else {
            line_data
        }
    })

all_curves <- map_df(seq_len(nrow(target_cpg_info)), function(i) {
    cpg_pos <- target_cpg_info$cpg_pos[i]
    data.frame(
        x = c(
            min(target_snp_pos, cpg_pos),
            (target_snp_pos + cpg_pos) / 2,
            max(target_snp_pos, cpg_pos)
        ),
        y = c(0.162, 0.222, 0.162),
        group = i
    )
})

all_gene_curves <- map_df(seq_len(nrow(target_genes_info)), function(i) {
    gene_pos <- if (target_genes_info$strand[i] == "+") {
        target_genes_info$start[i]
    } else {
        target_genes_info$end[i]
    }
    data.frame(
        x = c(
            min(target_snp_pos, gene_pos),
            (target_snp_pos + gene_pos) / 2,
            max(target_snp_pos, gene_pos)
        ),
        y = c(0.138, 0.078, 0.138),
        group = paste0("gene_", i)
    )
})

genomic_region_plot <- ggplot() +
    # Gene bodies.
    lapply(unique(gene_shapes %>% dplyr::filter(type == "gene_line") %>% pull(gene)), function(g) {
        geom_line(
            data = gene_shapes %>% dplyr::filter(gene == g, type == "gene_line"),
            aes(x = x, y = y),
            color = "#595959",
            linewidth = 0.5
        )
    }) +
    # Exons.
    lapply(unique(gene_shapes %>% dplyr::filter(type == "exon") %>% pull(gene)), function(g) {
        geom_rect(
            data = gene_shapes %>% dplyr::filter(gene == g, type == "exon"),
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = "#595959",
            color = "#595959",
            linewidth = 0.5
        )
    }) +
    # Gene labels and SNP-to-gene curves / dashed connectors.
    lapply(seq_len(nrow(target_genes_info)), function(i) {
        gene_info <- target_genes_info[i, ]
        dash_pos <- if (gene_info$strand == "+") {
            gene_info$start
        } else {
            gene_info$end
        }
        list(
            geom_text(
                data = gene_info,
                aes(x = (start + end) / 2,
                    y = 0 + 0.015 * (gene_id - 1),
                    label = gene),
                size = 3
            ),
            geom_bezier(
                data = all_gene_curves,
                aes(x = x, y = y, group = group),
                color = "blue",
                linewidth = 0.3,
                alpha = 0.8
            ),
            geom_segment(
                data = gene_info,
                aes(x = dash_pos, xend = dash_pos,
                    y = 0, yend = 0.14),
                linetype = "dashed",
                color = "red",
                linewidth = 0.5,
                alpha = 0.6
            )
        )
    }) +
    geom_rect(
        data = promoter_in_region,
        aes(xmin = start, xmax = end, ymin = 0.14, ymax = 0.16, fill = "Promoter"),
        alpha = 1
    ) +
    geom_rect(
        data = proximal_enhancer_in_region,
        aes(xmin = start, xmax = end, ymin = 0.14, ymax = 0.16, fill = "Proximal Enhancer"),
        alpha = 1
    ) +
    geom_rect(
        data = distal_enhancer_in_region,
        aes(xmin = start, xmax = end, ymin = 0.14, ymax = 0.16, fill = "Distal Enhancer"),
        alpha = 1
    ) +
    geom_rect(
        data = insulator_in_region,
        aes(xmin = start, xmax = end, ymin = 0.14, ymax = 0.16, fill = "Insulator"),
        alpha = 1
    ) +
    geom_point(
        data = target_snp_info,
        aes(x = target_snp_pos, y = 0.15),
        shape = 23,
        fill = "white",
        size = 2,
        stroke = 0.5
    ) +
    geom_text(
        aes(x = target_snp_pos, y = 0.23, label = target_snp),
        size = 3, vjust = 0, color = "black"
    ) +
    geom_point(
        data = target_cpg_info,
        aes(x = cpg_pos, y = 0.162),
        shape = 24,
        fill = "blue",
        size = 0.5,
        stroke = 0.5
    ) +
    geom_bezier(
        data = all_curves,
        aes(x = x, y = y, group = group),
        color = "#e31a1c",
        linewidth = 0.3,
        alpha = 0.8
    ) +
    scale_fill_manual(
        name = "Regulatory\nElements",
        values = c("Promoter" = "#f81212", "Proximal Enhancer" = "#e49c3d",
                   "Distal Enhancer" = "#1c9e74", "Insulator" = "#00b3fb"),
        guide = guide_legend(
            override.aes = list(alpha = 1, size = 2),
            title.position = "top",
            title.theme = element_text(
                size = 8,
                face = "bold",
                margin = margin(b = 3)
            ),
            label.theme = element_text(size = 7),
            direction = "vertical",
            ncol = 1,
            keywidth = unit(0.3, "cm"),
            keyheight = unit(0.3, "cm"),
            label.position = "right",
            label.margin = margin(l = 2, r = 2)
        )
    ) +
    # Manual SNP / CpG legend markers.
    geom_point(
        aes(x = region_start + 122000, y = 0.375),
        shape = 23, fill = "white", size = 2
    ) +
    annotate(
        "text", x = region_start + 117000, y = 0.375,
        label = "SNP", size = 3
    ) +
    geom_point(
        aes(x = region_start + 122000, y = 0.345),
        shape = 24, fill = "blue", size = 0.5
    ) +
    annotate(
        "text", x = region_start + 117000, y = 0.345,
        label = "CpG", size = 3
    ) +
    scale_x_continuous(
        name = paste0(region_chr, " (kb)"),
        limits = c(region_start, region_end),
        breaks = seq(
            round(region_start / 10000) * 10000,
            region_end,
            by = 25000
        ),
        labels = scales::comma_format(scale = 1 / 1000)
    ) +
    scale_y_continuous(breaks = NULL, limits = c(0, 0.5)) +
    labs(title = paste("Genomic Region around", target_snp)) +
    theme_minimal() +
    theme(
        axis.title.y = element_blank(),
        panel.grid = element_blank(),
        plot.title = element_text(size = 10),
        legend.position = c(0.05, 0.85),
        legend.justification = c(0, 1),
        legend.box = "horizontal",
        legend.background = element_rect(
            fill = "white",
            color = "gray50",
            size = 0.2
        ),
        legend.margin = margin(3, 5, 5, 5),
        axis.line.x = element_line(color = "black", linewidth = 0.5),
        axis.ticks.x = element_line(color = "black", linewidth = 0.5),
        axis.ticks.length.x = unit(2, "pt"),
        axis.title.x = element_text(margin = margin(t = 10))
    )

ggsave(paste0(target_snp, "_plot_with_legend.pdf"),plot = genomic_region_plot, width = 10, height = 4, dpi = 300, device = "pdf")


#######################################################
##### Fig.6C
#######################################################
setwd("/data/work/COACH_QTL_analysis/mediation/result")
bmediatR_coloc_SME = read.table(file.path(PATH,"bmediatR_coloc_SME_08.txt"),sep = "\t",header = T,check.names = F)
bmediatR_coloc_SME$ID = str_c(bmediatR_coloc_SME$SNP,"_",bmediatR_coloc_SME$CpG,"_",bmediatR_coloc_SME$Gene)
bmediatR_coloc_SME$class <- case_when(
  bmediatR_coloc_SME$complete_mediation_prob > 0.5 ~ "Complete mediation",
  bmediatR_coloc_SME$partial_mediation_prob > 0.5 ~ "Partial mediation",
  bmediatR_coloc_SME$colocalization_prob > 0.5 ~ "Co-local",
  TRUE ~ "Other")
table(bmediatR_coloc_SME$class)
Ratio <- bmediatR_coloc_SME %>% 
  group_by(class) %>% 
  count()
Ratio = Ratio[Ratio$class %in% c("Complete mediation","Partial mediation","Co-local","Other"),]
Ratio$Freq = round(Ratio$n/sum(Ratio$n)*100,1)
Ratio$class = factor(Ratio$class, levels = rev(c("Complete mediation","Partial mediation","Co-local","Other")))
Ratio$label = str_c(Ratio$Freq,"% (",Ratio$n,"/",sum(Ratio$n),")")
Ratio = Ratio[Ratio$class!="Other",]

ggplot(data=Ratio,aes(x=n,y=class))+
  geom_col(aes(fill=class))+
  scale_x_continuous(expand=expansion(mult=c(0,0)))+
  theme_classic()+
  geom_text(aes(label = label),nudge_x = -2400, size=5) +
  scale_fill_manual(values = c("Complete mediation"="#fc8eb7",
                               "Partial mediation"="#F1CAA8FF",
                               "Co-local"="#98a2c6",
                               "Other"="#d8dddc"))+
  #scale_y_discrete(labels=c("MRT67307","GW6741","Vehicle","Vehicle"))+
  labs(y=NULL,x="Triples")+
  theme(legend.position = "none",
        axis.title = element_text(size = 15,colour = "black"),
        axis.text = element_text(size = 12,colour = "black"),
        plot.margin = unit(c(0.1,2,0.1,0.1),'cm'),
        ggh4x.axis.ticks.length.minor= rel(0.9),
        axis.ticks.length.x = unit(0.4,'lines'))
ggsave("/data/work/COACH_QTL_analysis/mediation/plot/Proportion_different_models_each_SME_triplets.pdf", height = 6, width = 10)

#######################################################
##### Fig.6D
#######################################################
library(ggtern)
setwd("/data/work/COACH_QTL_analysis/mediation/data")
data = as.data.frame(fread("Standardized_Fig.6D_data_input.txt", header=T, sep="\t",data.table=FALSE))
ggtern()+
  geom_mask()+
  geom_point(data=data[is.na(data$label),],
             aes(x=complete_mediation_prob,
                 y=partial_mediation_prob,
                 z=colocalization_prob,
                 color=class,
                 size=log10_P.meQTL),stroke = 0.1,alpha=0.4)+
  geom_point(data=data[!is.na(data$label),],
             aes(x=complete_mediation_prob,
                 y=partial_mediation_prob,
                 z=colocalization_prob,
                 color=class,
                 size=log10_P.meQTL),stroke = 0,alpha=0.8)+
  geom_text(data=data[!is.na(data$label),],
            aes(x=complete_mediation_prob,
                y=partial_mediation_prob,
                z=colocalization_prob,label = label),color="black", size=3) + 
  theme_classic() +
  theme_showarrows()+
  labs(x="Complete mediation",y="Partial mediation",z="Co-local","color"="Model")+
  scale_color_manual(values = c("Complete mediation"="#fc8eb7",
                                "Partial mediation"="#F1CAA8FF",
                                "Co-local"="#98a2c6",
                                "Other"="#d8dddc",
                                "show"="black"))+
  guides(size = guide_legend(override.aes = list(size =c(1,2,3))),
         color=guide_legend(override.aes = list(size=3)),
         fill = "none"# Remove background tag
  )+
  scale_size("-log10(P-value)",range=c(0,3),breaks=c(25,50,75,100),labels=c("25","50","75","100"))+
  theme(panel.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 0.6),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        
        plot.margin = unit(c(1, 1, 1, 2), "cm"),
        plot.title = element_text(size = 28, hjust = 0.5),
        legend.position = "right",
        strip.background = element_rect(colour = "black", fill = "white", size = 0.3),
        strip.text = element_text(size = 14))
ggsave("/data/work/COACH_QTL_analysis/mediation/plot/Posterior_probabilities_different_models_each_SME_triplets.pdf", height = 8, width = 10)
