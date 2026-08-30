#' -----------------------------------------------------------------------------
#' Supplementary Figure - plot Script
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

#######################################################
##### Additional file 1: Fig. S1D
#######################################################
setwd("/data/work/COACH_QTL_analysis/covariates/data")
data = as.data.frame(fread("Standardized_Fig.S1B_data_input.txt", header=T, sep="\t",data.table=FALSE,nThread=5))
data$omics = factor(data$omics,levels = c("RNA","Methylation"))
data$X = factor(data$X,levels = c(str_c("RNA PC",1:30),str_c("Methylation PC",1:7)))
data$feature = factor(data$feature,levels = clinical_feacture)
p <- ggplot() +
  geom_tile(
    data = data,
    aes(feature,X, fill = R2), 
    colour = "white", size = 0.5) +
  scale_fill_gradient2(low = "#003366",high = "firebrick3",midpoint=0,na.value = "white",guide = "colourbar",name = expression(R^2))+
  labs(x = "",y = "") +
  scale_x_discrete(position = "bottom") +
  ggtitle("Correlation of omics PEER factors with known covariates")+
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    plot.title = element_text(size = 15, face = "bold", hjust = 0.5),
    axis.text.x = element_text(color="black",size=12,vjust=0.5,hjust=1,angle = 90),
    axis.text.y = element_text(color="black",size=11,margin=margin(0.5,0.4,0,0.5,'cm')),
    legend.position = "right", 
    axis.ticks.y = element_blank(),
    legend.direction="vertical",
    legend.title = element_text(hjust=0.5),
    legend.box.just = "left",
    plot.margin = margin( t = 0.2,
                          r = 0.9,
                          b = 0.2,
                          l = 0.2,
                          unit = "cm")
  )
omics.class <- data$X %>% as.data.frame() %>%
  mutate(group=data$omics) %>%
  mutate(p="") %>%
  ggplot(aes(.,p,fill=group))+
  geom_tile() + 
  scale_x_discrete(position="right") +
  scale_fill_manual(values = c("RNA"="#4A90E2","Methylation"="#FB8D62"))+
  theme_minimal()+xlab(NULL) + ylab(NULL) +
  theme(panel.grid = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x =element_blank(),
        axis.ticks.x = element_blank(),
        #legend.direction = "horizontal",
        legend.position = c(50,0.8)
  )+
  labs(fill = "omics")+
  coord_flip()

bottom <- ggplotGrob(omics.class)
p+annotation_custom(bottom,xmin=-0.9,xmax=0.7,ymin=-0.4,ymax=38.4)
ggsave("/data/work/COACH_QTL_analysis/covariates/plot/COACH_Correlation_Fig.S1B.plot.pdf", p, width=8, height=8)


#######################################################
##### Additional file 1: Fig. S1H
#######################################################
overlap_file <- "/data/work/COACH_QTL_analysis/replication/result/overlap_meQTL_TCGA_all.txt"
overlap_data <- as.data.frame(fread(overlap_file, sep = "\t", header = TRUE,
                                    data.table = FALSE, nThread = 50))

overlap_data$beta_scaled <- rescale(overlap_data$beta, to = c(-1, 1))
overlap_data$tcga_beta_scaled <- rescale(overlap_data$TCGA_beta, to = c(-1, 1))

direction_plot <- ggplot(overlap_data, aes(x = beta_scaled, y = tcga_beta_scaled)) +
  rasterise(
    geom_pointdensity(size = 0.01),
    dpi = 300
  ) +
  scale_color_viridis(
    name = NULL,
    labels = NULL,
    guide = guide_colorbar(
      ticks = TRUE,
      ticks.colour = "white",
      frame.colour = "white"
    )
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.5) +
  geom_smooth(
    method = "lm", se = FALSE,
    color = "red2", linewidth = 1.2, linetype = "dashed"
  ) +
  theme(
    panel.background = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, size = 0.6),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks = element_line(size = 0.6, colour = "black"),
    axis.text = element_blank(),
    axis.title = element_blank(),
    plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
    legend.position = "right",
    legend.text = element_blank(),
    legend.title = element_blank()
  )

ggsave("/data/work/COACH_QTL_analysis/replication/plot/Direction_of_FDR_005_meQTL_between_COCC_and_TCGA_plot_rasterized_new.png",
       plot = direction_plot,
       width = 7,
       height = 6,
       dpi = 300,
       bg = "white")


#######################################################
##### Additional file 1: Fig. S2A
#######################################################
y2_plot <- function(p1, p2) {
  # Obtain the underlying grobs of the two plots.
  g1 <- ggplotGrob(p1)
  g2 <- ggplotGrob(p2)
  pos <- c(subset(g1$layout, name == "panel", select = t:r))
  
  # Add the panel body of g2 onto the panel of g1.
  g <- gtable_add_grob(
    g1,
    g2$grobs[[which(g2$layout$name == "panel")]],
    pos$t, pos$l, pos$b, pos$l
  )
  
  # Move the y-axis (and label) of g2 to the right side of the merged plot.
  index <- which(g2$layout$name == "ylab-l")
  ylab <- g2$grobs[[index]]
  
  # Swap edge widths so the label sits on the right.
  widths <- ylab$widths
  ylab$widths[1] <- widths[3]
  ylab$widths[3] <- widths[1]
  ylab$vp[[1]]$layout$widths[1] <- widths[3]
  ylab$vp[[1]]$layout$widths[3] <- widths[1]
  
  # Flip the text alignment for the right-hand label.
  ylab$children[[1]]$hjust <- 1 - ylab$children[[1]]$hjust
  ylab$children[[1]]$vjust <- 1 - ylab$children[[1]]$vjust
  ylab$children[[1]]$x <- unit(1, "npc") - ylab$children[[1]]$x
  
  g <- gtable_add_cols(g, g2$widths[g2$layout[index, ]$l], pos$r)
  g <- gtable_add_grob(g, ylab, pos$t, pos$r + 1, pos$b, pos$r + 1,
                       clip = "off", name = "ylab-r")
  
  # Move the y-axis tick marks and tick labels to the right.
  index <- which(g2$layout$name == "axis-l")
  yaxis <- g2$grobs[[index]]
  ticks <- yaxis$children[[2]]
  
  # Swap the tick and tick-label positions.
  ticks$widths <- rev(ticks$widths)
  ticks$grobs <- rev(ticks$grobs)
  
  # Shift the tick marks.
  ticks$grobs[[1]]$x <- ticks$grobs[[1]]$x - unit(1, "npc") + unit(3, "pt")
  
  # Mirror the tick-label width and alignment for the right side.
  ticks$grobs[[2]]$widths[1] <- widths[3]
  ticks$grobs[[2]]$widths[3] <- widths[1]
  ticks$grobs[[2]]$vp[[1]]$layout$widths[1] <- widths[3]
  ticks$grobs[[2]]$vp[[1]]$layout$widths[3] <- widths[1]
  ticks$grobs[[2]]$children[[1]]$hjust <- 1 - ticks$grobs[[2]]$children[[1]]$hjust
  ticks$grobs[[2]]$children[[1]]$vjust <- 1 - ticks$grobs[[2]]$children[[1]]$vjust
  ticks$grobs[[2]]$children[[1]]$x <- unit(1, "npc") - ticks$grobs[[2]]$children[[1]]$x
  yaxis$children[[2]] <- ticks
  
  g <- gtable_add_cols(g, g2$widths[g2$layout[index, ]$l], pos$r)
  g <- gtable_add_grob(g, yaxis, pos$t, pos$r + 1, pos$b, pos$r + 1,
                       clip = "off", name = "axis-r")
  
  return(g)
}

##### Load SNP coordinates #####
bim_file <- "/data/work/COACH_QTL_analysis/genotype/250s.snp.vqsr_pass.plink.bim"
bim <- as.data.frame(fread(
  cmd = str_c("cat ", bim_file, " | cut -f 1,2,4,5,6"),
  sep = "\t", header = FALSE, data.table = FALSE, nThread = 30
))
colnames(bim) <- c("Chr", "SNP", "Pos", "Ref", "Alt")
bim$Chr <- str_c("chr", bim$Chr)

##### Load gene coordinates (TSS) #####
gene_loc_file <- "/data/work/COACH_QTL_analysis/molecular_phenotype/RNA/result/COACH_250sample_TPM.expression.gene.loc.txt.gz"
eqtl_gene_loc <- as.data.frame(fread(
  gene_loc_file, sep = "\t", header = TRUE, data.table = FALSE, nThread = 50
))
eqtl_gene_loc <- eqtl_gene_loc[, c("gene_id", "end")]
colnames(eqtl_gene_loc) <- c("gene", "gene_TSS")

##### Load eQTL results and compute the distance to TSS #####
eqtl_file <- "/data/work/COACH_QTL_analysis/eQTL_mapping/result/COACH_250_eQTL_Tumor_cis.txt"
eqtl <- as.data.frame(fread(
  eqtl_file, sep = "\t", header = TRUE, data.table = FALSE, nThread = 50
))

fdr_threshold <- 0.05
eqtl$FDR <- p.adjust(eqtl$`p-value`, "BH")
eqtl <- eqtl[eqtl$FDR < fdr_threshold, ]

eqtl <- merge(bim, eqtl, by = "SNP", all.y = TRUE)
eqtl <- merge(eqtl, eqtl_gene_loc, by = "gene", all.x = TRUE)

eqtl$distance <- eqtl$Pos - eqtl$gene_TSS
eqtl$log10_fdr <- -log10(eqtl$FDR)
eqtl$distance <- eqtl$distance / 1000

##### Build the two panels and overlaid them with the secondary y-axis #####
message("Drawing eQTL distance distribution...")

eqtl_point_plot <- ggplot(data = eqtl, aes(x = distance, y = log10_fdr)) +
  # Rasterise the point layer to keep the large PDF small.
  rasterise(
    geom_point(size = 0.5, alpha = 0.7, color = "#9181BB"),
    dpi = 300
  ) +
  coord_cartesian(xlim = c(-1000, 1000)) +
  scale_y_continuous(n.breaks = 5) +
  scale_x_continuous(n.breaks = 5) +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(linewidth = 0.8, fill = NA),
    panel.background = element_blank(),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    axis.text = element_text(size = 10)
  ) +
  labs(y = "-Log10 p-value of eQTLs", x = "eQTL:distance to TSS (kb)")

eqtl_density_plot <- ggplot(data = eqtl, aes(x = distance)) +
  geom_density(aes(x = distance, y = after_stat(density)), color = "red") +
  coord_cartesian(xlim = c(-1000, 1000)) +
  scale_y_continuous(n.breaks = 5) +
  scale_x_continuous(n.breaks = 5) +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(linewidth = 0.8, fill = NA),
    panel.background = element_blank(),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    axis.text = element_text(size = 10)
  ) +
  labs(y = "Density", x = "eQTL:distance to TSS (kb)")

# The two panels share the x-axis; the first keeps its y-axis on the left and
# the second adds a y-axis on the right.
merged_plot <- y2_plot(eqtl_point_plot, eqtl_density_plot)

ggsave("/data/work/COACH_QTL_analysis/eQTL_mapping/plot/Distance_to_TSS_eQTLs_rasterized.pdf",plot = merged_plot, width = 6, height = 6, dpi = 300, device = "pdf")


#######################################################
##### Additional file 1: Fig. S2C
#######################################################
setwd("/data/work/COACH_QTL_analysis/functional_annotation/result")
genomic_fisher <- as.data.frame(fread("mean_eQTL_SNP_genetic_regulatory.fisher.txt", sep = "\t", header = TRUE,data.table = FALSE, nThread = 5))

element_order <- c("Missense", "5'UTR", "Splice region", "Synonymous",
                   "3'UTR", "Non coding exon", "Intron", "Intergenic")
genomic_fisher$Element <- factor(genomic_fisher$Element, levels = element_order)
genomic_fisher$OR_diff <- genomic_fisher$OR - 1
genomic_fisher$OR_down_diff <- genomic_fisher$OR_down - 1
genomic_fisher$OR_up_diff <- genomic_fisher$OR_up - 1

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
print(fold_enrichment_plot)
ggsave("/data/work/COACH_QTL_analysis/functional_annotation/plot/Fold_enrichment_of_eQTL_SNP_in_genetic_elements.pdf",
       plot = fold_enrichment_plot, width = 6, height = 4, dpi = 300, device = "pdf")


#######################################################
##### Additional file 1: Fig. S2I
#######################################################
setwd("/data/work/COACH_QTL_analysis/GWAS/data")
GWAS = fread(file.path(Path,"Standardized_QTL_GWAS_pvalue_different_thresholds.txt"),header = T,data.table=FALSE)
GWAS = GWAS %>% pivot_longer(cols = -c("QTL_name","threshold"),names_to = "cohort", values_to = "number")
GWAS$QTL_name = factor(GWAS$QTL_name,levels = c("eQTL","meQTL"))

GWAS$cohort = case_when(
  GWAS$cohort == "UKBB_GERA" ~ "UK Biobank and GERA",
  GWAS$cohort == "BBJ" ~ "Biobank Japan",
  GWAS$cohort == "finngen" ~ "FinnGen",
  GWAS$cohort == "PLCO_East" ~ "PLCO")
GWAS$ID = str_c(GWAS$QTL_name,":",GWAS$cohort)
GWAS$ID = factor(GWAS$ID,levels = unique(GWAS$ID))
GWAS$cohort = factor(GWAS$cohort,levels = c("UK Biobank and GERA","Biobank Japan","FinnGen","PLCO"))

GWAS$log10_P_value = -log10(GWAS$threshold)
GWAS$threshold = case_when(
  GWAS$threshold == 5e-02 ~ "5e-02",
  GWAS$threshold == 5e-03 ~ "5e-03",
  GWAS$threshold == 5e-04 ~ "5e-04",
  GWAS$threshold == 5e-05 ~ "5e-05",
  GWAS$threshold == 5e-06 ~ "5e-06",
  GWAS$threshold == 5e-07 ~ "5e-07",
  GWAS$threshold == 5e-08 ~ "5e-08")
GWAS$threshold = factor(GWAS$threshold,levels = c("5e-02","5e-03","5e-04","5e-05","5e-06","5e-07","5e-08"))
GWAS[GWAS$number==0,"number"] = NA
colnames(GWAS)

panel1 <- ggplot(GWAS,aes(threshold,ID)) +
  scale_color_gradientn(colours =  colorRampPalette(c("#F9D8B1","#881912"))(8))+
  theme_bw()+
  geom_point(aes(size=number,color=log10_P_value))+
  geom_text(aes(label = number),nudge_y = 0.2) +
  xlab("The statistical significance threshold (P-value)") + ylab(NULL)+
  labs(colour = expression(Log[10](P)))+
  theme(panel.border = element_rect(fill=NA,color="black", linewidth=1, linetype="solid"),
        #panel.grid = element_blank(),
        panel.grid.major = element_line(colour = "gray"),
        axis.text.x =element_text(colour = 'black',size =12,angle =0,hjust =0.5,vjust = 0.5),
        axis.title.x = element_text(size=16),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        legend.text = element_text(colour = 'black', size = 12),
        legend.title = element_text(colour = 'black', size = 12)
  )+
  scale_size("xQTLs count",range=c(2,15),breaks=rev(c(5,10,15)),labels=c("5000","10000","15000"))+
  guides(size = guide_legend(reverse=F,override.aes = list(size =c(6,7,8))),
         fill = "none"
  )

sidebar_a <- GWAS[,c("ID","cohort")]%>%
  ggplot(aes("Cohort",ID,fill=cohort))+
  geom_tile() + 
  scale_y_discrete(position="left") +
  theme_minimal()+xlab(NULL) + ylab(NULL) +
  theme(panel.grid = element_blank(),
        panel.border =element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        legend.text = element_text(colour = 'black', size = 12),
        legend.title = element_text(colour = 'black', size = 12),
        axis.text.x =element_text(colour = 'black',size =12,angle =90,hjust =0,vjust = 0.5)
  )+
  labs(x=NULL,fill = "Cohort",y=NULL)+
  scale_fill_manual(values=c("UK Biobank and GERA"="#449DB3FF","Biobank Japan"="#60BFAEFF","FinnGen"="#CFE2B4","PLCO"="#FBD6DD"))

sidebar_b <- GWAS[,c("ID","QTL_name")]%>%
  ggplot(aes("xQTL",ID,fill=QTL_name))+
  geom_tile() + 
  scale_y_discrete(position="left") +
  theme_minimal()+xlab(NULL) + ylab(NULL)+
  theme(panel.grid = element_blank(),
        panel.border =element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        legend.text = element_text(colour = 'black', size = 12),
        legend.title = element_text(colour = 'black', size = 12),
        axis.text.x =element_text(colour = 'black',size =12,angle =90,hjust =0,vjust = 0.5)
  )+
  labs(x=NULL,fill = "",y=NULL)+
  scale_fill_manual(values=c("eQTL"="#bebada","meQTL"="#80b1d3"))

P <- panel %>% insert_left(sidebar_a, width = .05) %>% insert_left(sidebar_b,width=.05)
setwd(Path)
ggsave("/data/work/COACH_QTL_analysis/GWAS/plot/GWAS_eQTL_meQTL_overlap_different_thresholds.pdf", P, height = 8, width = 8)



#######################################################
##### Additional file 1: Fig. S3A
#######################################################
#Density plot shows the distance from meQTL to the corresponding meCpG
setwd("/data/work/COACH_QTL_analysis/meQTL_clumping/result")
meQTL_cis_clumping = as.data.frame(fread("COACH_meQTL_CpG_clumping_clump.txt",sep = "\t",header = T,data.table=FALSE,nThread=10))
meQTL_cis_clumping$CpG_BP = as.numeric(unlist(lapply(strsplit(meQTL_cis_clumping$CpG,":"),function(x) x[[2]]))) 
meQTL_cis_clumping$distance = abs(meQTL_cis_clumping$CpG_BP-meQTL_cis_clumping$BP)/1000
data = meQTL_cis_clumping %>% group_by(CpG) %>% summarise(SNP_Count = length(unique(Lead_SNP)))
data = merge(meQTL_cis_clumping, data, by="CpG", all.x=T)

data$Class <- case_when(
  data$SNP_Count == 1 ~ "1",
  data$SNP_Count == 2 ~ "2",
  data$SNP_Count == 3 ~ "3",
  data$SNP_Count == 4 ~ "4",
  data$SNP_Count >= 5 & data$SNP_Count < 10 ~ ">=5",
  data$SNP_Count >= 10 ~ ">=10")
data$Class = factor(data$Class, levels=c("1","2","3","4",">=5",">=10"))
data_median = data %>% 
  group_by(Class) %>% 
  summarise(median_distance = median(distance))
table(data$Class)
color = c("1"="#3F88C5","2"="#F3AC66","3"="#019875FF","4"="#9041ba",">=5"="#E84A5FFF",">=10"="#fec10b")
ggplot(data[data$Class==c("1","2","3","4",">=5",">=10"),],aes(x=distance,fill=Class)) +
  geom_density(color="#e9ecef", alpha=0.5) +
  geom_vline(xintercept = data_median$median_distance,
             color=color,linetype="dashed") +
  theme_classic(base_rect_size = 0)+
  labs(x = "Distance to CpG site (kb)",y = "Density",fill="Number of lead meQTLs") +
  ggtitle("") +
  scale_y_continuous(expand=c(0.01,0),limits = c(0,0.2)) +
  scale_x_continuous(expand=c(0.01,0),limits = c(0,60)) +
  scale_fill_manual(values=color)+
  theme(panel.background = element_blank(),
        #panel.grid.major = element_blank(),
        #panel.grid.minor = element_blank(),
        #panel.border = element_rect(color = "black", fill = NA, size = 1),
        panel.grid = element_blank(),
        panel.spacing.x = unit(0,"cm"),
        panel.border = element_blank(),
        panel.spacing = unit(0,"lines"),
        axis.ticks = element_line(size=0.6,colour = "black"),
        axis.text.x = element_text(size = 10,colour = "black"),
        axis.title.x = element_text(size=16),
        axis.text.y = element_text(size=10,colour = "black"),
        axis.title.y = element_text(size=16),
        #legend.direction = "horizontal",
        legend.position = "right")
ggsave("/data/work/COACH_QTL_analysis/meQTL_clumping/plot/distance_meQTL_corresponding_meCpG_density.pdf", height = 6, width = 10)



#######################################################
##### Additional file 1: Fig. S4B
#######################################################
setwd("/data/work/COACH_QTL_analysis/colocalization/eQTL_meQTL_coloc/result")
eQTL_meQTL_clump_coloc = read.table("eQTL_meQTL_clump_coloc.summary.txt",sep = "\t",header = T,check.names = F)
eQTL_meQTL_clump_coloc$T_label <- case_when(
  eQTL_meQTL_clump_coloc$Assigned.Tier == "T1" ~ str_c("T1: ",nrow(eQTL_meQTL_clump_coloc[!is.na(eQTL_meQTL_clump_coloc$Assigned.Tier) & eQTL_meQTL_clump_coloc$Assigned.Tier == "T1",])),
  eQTL_meQTL_clump_coloc$Assigned.Tier == "T2" ~ str_c("T2: ",nrow(eQTL_meQTL_clump_coloc[!is.na(eQTL_meQTL_clump_coloc$Assigned.Tier) & eQTL_meQTL_clump_coloc$Assigned.Tier == "T2",])),
  eQTL_meQTL_clump_coloc$Assigned.Tier == "T3" ~ str_c("T3: ",nrow(eQTL_meQTL_clump_coloc[!is.na(eQTL_meQTL_clump_coloc$Assigned.Tier) & eQTL_meQTL_clump_coloc$Assigned.Tier == "T3",])),
  eQTL_meQTL_clump_coloc$Assigned.Tier == "T4" ~ str_c("T4: ",nrow(eQTL_meQTL_clump_coloc[!is.na(eQTL_meQTL_clump_coloc$Assigned.Tier) & eQTL_meQTL_clump_coloc$Assigned.Tier == "T4",])),
  eQTL_meQTL_clump_coloc$Assigned.Tier == "T5" ~ str_c("T5: ",nrow(eQTL_meQTL_clump_coloc[!is.na(eQTL_meQTL_clump_coloc$Assigned.Tier) & eQTL_meQTL_clump_coloc$Assigned.Tier == "T5",])),
  TRUE ~ NA)
ggplot() +
  geom_point(data = eQTL_meQTL_clump_coloc[is.na(eQTL_meQTL_clump_coloc$Assigned.Tier),],
             aes(x = Signed_log10_pvalue_meQTL, y = Signed_log10_adjust_pvalue_geneEffectScore),
             color = "grey",
             size = 1,
             show.legend = NA)+ 
  geom_point(data = eQTL_meQTL_clump_coloc[!is.na(eQTL_meQTL_clump_coloc$Assigned.Tier),],
             aes(x = Signed_log10_pvalue_meQTL, 
                 y = Signed_log10_adjust_pvalue_geneEffectScore,
                 color = T_label),
             size = 2.2,
             show.legend = NA)+ 
  scale_color_manual(values=c("T1: 8"="#8f2638","T2: 11"="#fd77ab","T3: 14"="#ff7323","T4: 38"="#2082a2","T5: 18"="#3a5fb1"))+
  geom_hline(yintercept = c(log10(0.01)),size = 0.5,color = "red",lty = "dashed") + #Horizontal threshold line
  geom_vline(xintercept = c(-log10(2e-10)),size = 0.5,color = "red",lty = "dashed") +
  geom_vline(xintercept = c(log10(2e-10)),size = 0.5,color = "red",lty = "dashed") +
  scale_x_continuous(limits = c(-100,100))+
  scale_y_continuous(limits = c(-50,0))+
  geom_label_repel(data = eQTL_meQTL_clump_coloc[!is.na(eQTL_meQTL_clump_coloc$label_gene),],
                   aes(x = Signed_log10_pvalue_meQTL,
                       y = Signed_log10_adjust_pvalue_geneEffectScore,
                       label = label_gene,
                       fill = Assigned.Tier),
                   color = "white",
                   parse=T,#color=`Up/Down-Regulation`
                   box.padding = 0.5, label.padding = 0.25, 
                   point.padding = 1e-06, label.r = 0.15, label.size = 0.1, min.segment.length = 0,
                   arrow=NULL, force=20, force_pull=1, nudge_x = 0, nudge_y = 0, seed = 42, xlim = c(NA, NA),
                   segment.size = 0.6, segment.color = "black", show_guide = FALSE) +
  scale_fill_manual(values=c("T1"="#8f2638","T2"="#fd77ab","T3"="#ff7323","T4"="#2082a2","T5"="#3a5fb1"))+
  labs(y = expression(Signed -log[10]("adjusted p-value CRISPR gene effect")),
       x = expression(Signed -log[10]("p-value of meQTL in both meQTL-eQTL coloc")),
       color = str_c("Total: ",length(unique(eQTL_meQTL_clump_coloc$gene))),
       fill = "")+
  ggtitle("Prioritization of targetable mGenes in both meQTL-eQTL coloc")+
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 14),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16),
        plot.title = element_text(size = 16, face = "plain", hjust = 0.5),)
ggsave("/data/work/COACH_QTL_analysis/colocalization/eQTL_meQTL_coloc/plot/Prioritization_targetable_mGenes_both_meQTL-eQTL_coloc.pdf", height = 8, width = 7)


#######################################################
##### Additional file 1: Fig. S4F
#######################################################
suppressMessages(library(ggunchained))
suppressMessages(library(gghalves))
setwd("/data/work/COACH_QTL_analysis/molecular_phenotype/methylation/data")
data = as.data.frame(fread("250s_T_methy_CG_10x.matirx_Standardized_PYGL.txt",header = T,sep = "\t",data.table=FALSE,nThread=10))
ggplot() +
  # Draw a split violin plot
  ###Left half: Group 1
  geom_half_violin(data = data[data$class == "Normal",],
                   aes(x = CpG,y = value),colour = "white",fill = "#b75555",
                   side = "l",nudge = 0.01) +
  ###Right half: Group 2
  geom_half_violin(data = data[data$class == "Tumor",],
                   aes(x = CpG,y = value),colour = "white",fill = "#4576af",
                   side = "r",nudge = 0.01) + 
  ###Add the mean points
  geom_point(data = data, aes(x = CpG,y = value, fill = class),
             stat = 'summary', fun = mean,size = 1,
             show.legend = T,
             position = position_dodge(width = 0.4))+
  ###Add errorbars
  stat_summary(data = data, aes(x = CpG, y = value, fill = class),
               fun.min = function(x){quantile(x)[2]},
               fun.max = function(x){quantile(x)[4]},
               geom = 'errorbar', color = 'black',
               width = 0.1,size = 0.5,
               position = position_dodge(width = 0.4)) + 
  stat_compare_means(data = data,aes(x = CpG,y=value,group = class),
                     method="wilcox.test",
                     show.legend = NA,
                     label="p.signif",
                     label.y=2.5) +
  scale_fill_manual(values = c("Tumor"="#4576af","Normal"="#b75555")) +
  scale_y_continuous(limits = c(-2.5,2.7))+
  labs(x=NULL,y=NULL) +
  guides(fill=guide_legend(position = "top")) +
  theme_test() +
  theme(strip.placement = "outside",
        strip.background = element_blank(),
        strip.text=element_text(color="black"),
        axis.text.x=element_text(color="black",size=10,angle = 90,vjust=0.5,hjust=0.5),
        axis.text.y=element_text(color="black",size=10),
        legend.title = element_blank(),
        legend.direction="horizontal",
        legend.key.height = unit(0.4,"cm"),
        legend.key.width = unit(0.4,"cm"),
        legend.background = element_blank(),
        legend.position.inside = c(0.5,0.95))
ggsave("/data/work/COACH_QTL_analysis/molecular_phenotype/methylation/plot/COACH_250_PYGL_41CpGs_methy_diff_T-vs-N_wilcox.pdf", height = 4, width = 12)


#######################################################
##### Additional file 1: Fig. S5A
#######################################################
setwd("/data/work/COACH_QTL_analysis/mediation/data")
mediation_result = as.data.frame(fread("mediation_result_sobel_08.txt",sep = "\t",header = T,data.table=FALSE,nThread=10))
ggplot(mediation_result, aes(x=proportion, y=pvalue,group=Type)) +
  geom_point_rast(aes(color=Type),size=0.7)+
  scale_color_manual(values=c("SME"='#c55e79',"SEM"='#4173b2'))+
  labs(x="Mediation proportion", y = "-log10(Sobel pvalues)",color="")+
  #theme_bw(base_size=22)+
  #theme_classic(base_size = 14)+
  geom_vline(aes(xintercept = 0), size = 0.4,color = "black", linetype = "dashed")+
  geom_hline(aes(yintercept = 0), size = 0.4,color = "black", linetype = "dashed")+
  #
  # geom_hline(aes(yintercept = log10(0.05)),size = 0.5,color = "red",linetype = "dashed")+ 
  # geom_hline(aes(yintercept = -log10(0.05)),size = 0.5,color = "red",linetype = "dashed")+ 
  theme(panel.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 0.8),
        axis.ticks = element_line(size=0.6,colour = "black"),
        axis.title = element_text(size=16),
        axis.text = element_text(size=14,colour = "black"),
        plot.margin=unit(c(0.5,0.5,0.5,0.5),units=,"cm"),
        legend.key.size = unit(0.5, "cm"),
        legend.background = element_blank(),
        legend.key = element_blank(),
        legend.text = element_text(colour = 'black', size = 12),
        legend.direction = "horizontal",legend.position = c(0.2,0.94))
ggsave("/data/work/COACH_QTL_analysis/mediation/plot/Mediation_analysis_for_SME_and_SEM_Sobel_tests.pdf", height = 6, width = 6,dpi = 800)


#######################################################
##### Additional file 1: Fig. S5B
#######################################################
setwd("/data/work/COACH_QTL_analysis/mediation/data")
mediation_result = as.data.frame(fread("mediation_result_sobel_08.txt",sep = "\t",header = T,data.table=FALSE,nThread=10))
ggplot(mediation_result, aes(x=proportion,fill=Type)) +
  geom_histogram(bins = 40,alpha=0.5)+
  geom_vline(aes(xintercept=median(mediation_result[mediation_result$Type=="SME","proportion"],na.rm = TRUE)),size=0.7,linetype="dashed",color="#0072B1")+ 
  geom_vline(aes(xintercept=median(mediation_result[mediation_result$Type=="SEM","proportion"],na.rm = TRUE)),size=0.7,linetype="dashed",color="#FB8D62")+ 
  #geom_vline(aes(xintercept=0),color="black")+
  scale_fill_manual(values=c("SME"='#0072B1',"SEM"='#FB8D62'))+
  labs(fill="",x=paste0("Mediation Proportion"), y = "Count")+
  #theme_bw(base_size = 22)+
  xlim(0,1.5)+
  #geom_density(alpha=0.6)+
  theme(panel.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 0.8),
        axis.ticks = element_line(size=0.6,colour = "black"),
        axis.title = element_text(size=16),
        axis.text = element_text(size=14,colour = "black"),
        plot.margin=unit(c(0.5,0.5,0.5,0.5),units=,"cm"),
        legend.key.size = unit(0.5, "cm"),
        legend.background = element_blank(),
        legend.key = element_blank(),
        legend.text = element_text(colour = 'black', size = 14),
        legend.direction = "horizontal",legend.position = c(0.7,0.94))
ggsave("/data/work/COACH_QTL_analysis/mediation/plot/Histograms_mediation_proportion_significant_SEM_and_SME.pdf", height = 6, width = 6,dpi = 800)


#######################################################
##### Additional file 1: Fig. S5E
#######################################################
setwd("/data/work/COACH_QTL_analysis/mediation/result")
counts <- read.csv("SME_Tumor_Normal_triplets_comparison_counts.csv", stringsAsFactors = FALSE)

draw_triplet_barplot <- function() {
  bar_data <- as.matrix(counts[, c("Tumor", "Normal", "Overlap")])
  rownames(bar_data) <- counts$Category
  colnames(bar_data) <- c("Tumor", "Normal", "Overlap")
  
  group_colors <- c(
    "Tumor"   = "#E64B35",
    "Normal"  = "#4DBBD5",
    "Overlap" = "#7E57C2"
  )
  
  par(mar = c(5, 5, 4, 2))
  bar_positions <- barplot(
    t(bar_data),
    beside    = TRUE,
    names.arg = rownames(bar_data),
    col       = group_colors,
    border    = NA,
    ylim      = c(0, max(bar_data) * 1.35),
    ylab      = "Number of SNP-CpG-Gene Triplets",
    xlab      = "SME Mediation Category",
    main      = "Tumor vs Normal SME Triplet Comparison",
    cex.main  = 1.4,
    cex.lab   = 1.2,
    cex.axis  = 1.1,
    cex.names = 1.2
  )
  
  # Add value labels on top of each bar.
  for (i in seq_len(nrow(bar_data))) {
    for (j in seq_len(ncol(bar_data))) {
      val <- bar_data[i, j]
      x_pos <- bar_positions[j, i]
      text(x_pos, val + max(bar_data) * 0.03, labels = val,
           cex = 1.1, font = 2, col = "black")
    }
  }
  
  legend("topright",
         legend = c("Tumor", "Normal", "Overlap"),
         fill   = group_colors,
         border = NA,
         bty    = "n",
         cex    = 1.1,
         pt.cex = 2)
  
  # Add a subtle grid, then redraw the bars on top of it.
  grid(nx = NA, ny = NULL, col = "grey90", lty = "dashed")
  
  bar_positions <- barplot(
    t(bar_data),
    beside    = TRUE,
    names.arg = rownames(bar_data),
    col       = group_colors,
    border    = NA,
    ylim      = c(0, max(bar_data) * 1.35),
    ylab      = "Number of SNP-CpG-Gene Triplets",
    xlab      = "SME Mediation Category",
    main      = "Tumor vs Normal SME Triplet Comparison",
    cex.main  = 1.4,
    cex.lab   = 1.2,
    cex.axis  = 1.1,
    cex.names = 1.2,
    add       = TRUE
  )
  
  # Redraw the value labels on top of the redrawn bars.
  for (i in seq_len(nrow(bar_data))) {
    for (j in seq_len(ncol(bar_data))) {
      val <- bar_data[i, j]
      x_pos <- bar_positions[j, i]
      text(x_pos, val + max(bar_data) * 0.03, labels = val,
           cex = 1.1, font = 2, col = "black")
    }
  }
  
  box(lwd = 2, col = "black")
}

pdf("/data/work/COACH_QTL_analysis/mediation/plot/Tumor-vs-Normal_SME_triplets_comparison_figure.pdf", width = 10, height = 7, useDingbats = FALSE)
draw_triplet_barplot()
dev.off()


#######################################################
##### Additional file 1: Fig. S6G
#######################################################
suppressMessages(library("ggVennDiagram"))
#GME
setwd("/data/work/COACH_QTL_analysis/colocalization/xQTL_GWAS_moloc/result")
xQTL_GWAS_moloc_clump_moloc = read.table("xQTL_GWAS_clump_moloc.summary.txt",sep = "\t",header = T,check.names = F)

#EM
setwd("/data/work/COACH_QTL_analysis/colocalization/eQTL_meQTL_coloc/result")
eQTL_meQTL_clump_coloc = read.table("eQTL_meQTL_clump_coloc.summary.txt",sep = "\t",header = T,check.names = F)

#GE
setwd("/data/work/COACH_QTL_analysis/colocalization/eQTL_GWAS_coloc/result")
eQTL_GCST90129505_clump_coloc = read.table("eQTL_GCST90129505_clump_coloc.summary.txt",sep = "\t",header = T,check.names = F)

#GM
setwd("/data/work/COACH_QTL_analysis/colocalization/meQTL_GWAS_coloc/result")
meQTL_GCST90129505_clump_coloc = read.table("meQTL_GCST90129505_clump_coloc.summary.txt",sep = "\t",header = T,check.names = F)

#CpG
ggVennDiagram(list(GME = sample(xQTL_GWAS_moloc_clump_moloc$CpG),
                   EM = sample(eQTL_meQTL_clump_coloc$CpG),
                   GM = sample(meQTL_GCST90129505_clump_coloc$CpG)),
              label = "count",
              label_alpha = 0,
              label_size = 6,
              category.names = c("GME", "EM", "GM"),
              edge_lty = 1,
              edge_size = 2,
              alpha =0.7,
              #fill = c("#586fc0", "#e07777", "#8ace8a"),
              set_size = 7,
              set_color = c("#586fc0", "#e07777", "#8ace8a")) +
  scale_fill_gradient(low = "white", high =  "white") +
  ggtitle("CpGs")+
  theme_void()+
  theme(legend.position = "none",
        plot.title = element_text(size = 20, hjust = 0.5))
ggsave("/data/work/COACH_QTL_analysis/colocalization/plot/GEM_EM_GM_colocalization_CpGs_overlap.pdf", height = 6, width = 8,dpi = 800)


#genes
ggVennDiagram(list(GME = sample(xQTL_GWAS_moloc_clump_moloc$gene),
                   EM = sample(eQTL_meQTL_clump_coloc$gene),
                   GE = sample(eQTL_GCST90129505_clump_coloc$gene)),
              label = "count",
              label_alpha = 0,
              label_size = 6,
              category.names = c("GME", "EM", "GE"),
              edge_lty = 1,
              edge_size = 2,
              alpha =0.7,
              #fill = c("#586fc0", "#e07777", "#8ace8a"),
              set_size = 7,
              set_color = c("#586fc0", "#e07777", "#8ace8a")) +
  scale_fill_gradient(low = "white", high =  "white") +
  theme_void()+
  ggtitle("Genes")+
  theme(legend.position = "none",
        plot.title = element_text(size = 20, hjust = 0.5))
ggsave("/data/work/COACH_QTL_analysis/colocalization/plot/GEM_EM_GE_colocalization_genes_overlap.pdf", height = 6, width = 8,dpi = 800)





