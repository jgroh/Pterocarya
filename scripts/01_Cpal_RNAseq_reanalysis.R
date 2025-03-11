library(data.table)
library(tximport)
library(forcats)
library("DESeq2")
library(ggplot2)
library(EnhancedVolcano)
library(tidyverse)
library(apeglm)
setwd("~/workspace/Pterocarya/Cpal_expression/")
tx2gene_pr <- fread("Chen_etal_2019/tx2gene_paralogs.txt", header = F, col.names = c("tx", "gene"))

# -------- Qu et al. 2023 -----
# https://doi.org/10.1016/j.gpb.2023.02.001.
# "Male and female floral buds in various PA and PG individuals 
# were collected at five different stages including 
# (1) S0, physiological differentiation period; 
# (2) S1, dormancy period; 
# (3) S2, germination period; 
# (4) S3, inflorescence elongation period; and 
# (5) S4, maturation period. 
# Three biological replicates of the extracted RNA from flowers 
# (in PA and PG individuals) were sequenced on Illumina NovaSeq platform"

Qu23.meta <- fread("Qu_etal_2023/RNAseq_meta.txt", header = T, col.names = c("Run", "sample"))
Qu23.meta[, c("stage", "type", "fl", "rep") := tstrsplit(sample, split = "")]
Qu23.meta[fl == 'F', fl2 := 'Pistillate (F)']
Qu23.meta[fl == 'M', fl2 := 'Staminate (M)']
Qu23.meta[type == 'A', type2 := 'Protogynous (F 1st)']
Qu23.meta[type == 'B', type2 := 'Protandrous (M 1st)']



# transcript quantifications, annotation lacking both alleles 
# (originally labelled as paralogs)
Qu23.quants_files_pr <- list.files("Qu_etal_2023/quants_paralogs/", recursive = T, pattern = ".*quant\\.sf$", full.names = T)
Qu23.all_names_pr <- gsub("_quant/quant.sf", "", gsub("quants_paralogs/", "", Qu23.quants_files_pr))
names(Qu23.quants_files_pr) <- Qu23.all_names_pr


Qu23.txi_pr <- tximport(Qu23.quants_files_pr, type = "salmon", tx2gene = tx2gene_pr, countsFromAbundance = 'no')
head(Qu23.txi_pr$counts)
head(Qu23.txi_pr$abundance)

Qu23.counts_pr <- data.table(Qu23.txi_pr$counts)
Qu23.counts_pr[, gene := rownames(Qu23.txi_pr$counts)]
Qu23.counts_pr <- melt(Qu23.counts_pr, id.vars = 'gene', variable.name = 'Run', value.name = 'count')

Qu23.counts_pr[, Run := gsub("Qu_etal_2023//", "", Run)]
Qu23.counts_pr <- merge(Qu23.counts_pr, Qu23.meta)


Qu23.counts.FAFL_pr <- Qu23.counts_pr[grepl("GFAFL", gene)]
Qu23.counts.FAFL_pr[, gene := gsub("2PA_", "", gene)]
Qu23.counts.FAFL_pr[, gene := gsub("FL", "FL.", gene)]


# zero counts from type A (protogynous)
Qu23.counts.FAFL_pr[type2 == 'Protogynous (F 1st)' & gene == 'GFAFL.2', count := NA]
#Qu23.counts.FAFL_pr[gene == 'GFAFL1' & stage == 4 & type2 == 'protandrous' & fl == 'M']

Qu23.counts.FAFL_pr[, stage := as.numeric(stage)]
Qu23_means <- Qu23.counts.FAFL_pr[, .(count = mean(count)), by = .(type2, fl2, gene, stage)]

# ----- plot using raw transcript counts from salmon -----
ggplot(Qu23.counts.FAFL_pr) + 
  geom_point(aes(x = stage, y = count, group = interaction(type2, gene), color = gene), size = 1.5, shape = 21, stroke = 1)  + 
  geom_line(data = Qu23_means, aes(x = stage, y = count, group = interaction(type2, gene, fl2), color = gene)) +
  geom_segment(data = Qu23_means, aes(x = stage-0.1, xend = as.numeric(stage)+0.1, y = count, yend = count, group = interaction(type2, gene, stage, fl2), color = gene)) +
  
  scale_color_manual(values = c("turquoise", "darkblue")) +
  facet_grid(fl2~type2, scales = 'free_y') + 
  #facet_grid(fl2~type2, scales = 'free_y') + 
  
  labs(color = '', x = 'Developmental stage',y= 'Estimated transcript count') +
  theme_bw() + 
  theme(aspect.ratio = 1,
        plot.margin = margin(t = 10, r = 10, b = 10, l = 10, unit = 'pt'),
        strip.background = element_rect(fill = NA), 
        legend.text = element_text(face = 'italic'))


# ----- plot using library-size corrected transcript counts from DESeq2 -------
# From DESeqDataSetFromTximport help page:
# Rows of colData correspond to columns of countData. Here, this is unecessary, but included for completenes.

all_order_indices <- match(gsub("Qu_etal_2023//","", names(as.data.table(Qu23.txi_pr$counts))), Qu23.meta$Run)
Qu23.meta <-  Qu23.meta[all_order_indices]                  



ddsTxi.all <- DESeqDataSetFromTximport(Qu23.txi_pr,
                                       colData = Qu23.meta,
                                       design = ~ type + fl + stage)

dds.all <- DESeq(ddsTxi.all)


# counts normalized by library size
G1_counts <- plotCounts(dds.all, gene='2PA_GFAFL1', intgroup=c("type", "fl", 'stage'), pc = 0, returnData = T)
G2_counts <- plotCounts(dds.all, gene='2PA_GFAFL2', intgroup=c("type", "fl", 'stage'), pc = 0, returnData = T)

# reformat to get metadata back
G1_counts <- cbind(Run = rownames(G1_counts),data.table(G1_counts))
G2_counts <- cbind(Run = rownames(G2_counts),data.table(G2_counts))

G1_counts[, Run := gsub("Qu_etal_2023//", "", Run)]
G2_counts[, Run := gsub("Qu_etal_2023//", "", Run)]

G1_counts <- merge(Qu23.meta, G1_counts, by = c("Run", "stage", "type", "fl"))
G2_counts <- merge(Qu23.meta, G2_counts, by = c("Run", "stage", "type", "fl"))

G1_counts[, allele := 'G1']
G2_counts[, allele := 'G2']

GFAFL_counts_final <- rbind(G1_counts, G2_counts)
GFAFL_counts_final[, stage:= as.numeric(as.character(stage))]
GFAFL_counts_final[, log.count := log(count + 0.5)]

GFAFL_mean_counts_raw <- GFAFL_counts_final[, .(count = mean(count)), by = .(stage, fl2, type2, allele)]
GFAFL_mean_counts_log <- GFAFL_counts_final[, .(log.count = mean(log.count+0.5)), by = .(stage, fl2, type2, allele)]
GFAFL_mean_counts_raw[, stage := as.numeric(as.character(stage))]
GFAFL_mean_counts_log[, stage := as.numeric(as.character(stage))]

GFAFL_counts_final[type == 'A' & allele == 'G2', count := NA]
GFAFL_mean_counts_raw[type2 == 'Protogynous (F 1st)' & allele == 'G2', count := NA]
GFAFL_mean_counts_log[type2 == 'Protogynous (F 1st)' & allele == 'G2', log.count := NA]
# plot on raw scale
ggplot(GFAFL_counts_final) + 
  geom_point(aes(x = stage, y = count, group = interaction(type2, allele), color = allele), size = 1.5, shape = 21, stroke = 1)  + 
  geom_line(data = GFAFL_mean_counts_raw, aes(x = stage, y = count, group = interaction(type2, allele, fl2), color = allele)) +
  geom_segment(data = GFAFL_mean_counts_raw, aes(x = stage-0.12, xend = stage+0.12, y = count, yend = count, group = interaction(type2, allele, stage, fl2), color = allele)) +
  
  scale_color_manual(values = c("turquoise", "darkblue")) +
  facet_grid(fl2~type2, scales = 'free_y') + 
  #facet_grid(fl2~type2, scales = 'free_y') + 
  
  labs(color = '', x = 'Developmental stage',y= 'Corrected transcript count') +
  theme_bw() + 
  theme(aspect.ratio = 1,
        plot.margin = margin(t = 10, r = 10, b = 10, l = 10, unit = 'pt'),
        strip.background = element_rect(fill = NA), 
        legend.text = element_text(face = 'italic'))


# plot on log scale
ggplot(GFAFL_counts_final) + 
  geom_point(aes(x = stage, y = log.count, group = interaction(type2, allele), color = allele), size = 1.5, shape = 21, stroke = 1)  + 
  geom_line(data = GFAFL_mean_counts_log, aes(x = stage, y = log.count, group = interaction(type2, allele, fl2), color = allele)) +
  geom_segment(data = GFAFL_mean_counts_log, aes(x = stage-0.1, xend = stage+0.1, y = log.count, yend = log.count, group = interaction(type2, allele, stage, fl2), color = allele)) +
  
  scale_color_manual(values = c("turquoise", "darkblue")) +
  facet_grid(fl2~type2, scales = 'free_y') + 
  #facet_grid(fl2~type2, scales = 'free_y') + 
  
  labs(color = '', x = 'Developmental stage',y= 'Estimated transcript count') +
  theme_bw() + 
  theme(aspect.ratio = 1,
        plot.margin = margin(t = 10, r = 10, b = 10, l = 10, unit = 'pt'),
        strip.background = element_rect(fill = NA), 
        legend.text = element_text(face = 'italic'))






# ------- Orthologs of other genes -------

# TPPD1 ortholog is CpaM1st31681
# counts normalized by library size
TPPD1.counts <- plotCounts(dds.all, gene='CpaM1st31681', intgroup=c("type", "fl", 'stage'), pc = 0, returnData = T)

# reformat to get metadata back
TPPD1.counts <- cbind(Run = rownames(TPPD1.counts),data.table(TPPD1.counts))
TPPD1.counts[, Run := gsub("Qu_etal_2023//", "", Run)]
TPPD1.counts <- merge(Qu23.meta, TPPD1.counts, by = c("Run", "stage", "type", "fl"))
TPPD1.counts[, stage:= as.numeric(as.character(stage))]
TPPD1.counts.mean <- TPPD1.counts[, .(count = mean(count)), by = .(stage, fl2, type2)]

TPPD1_plt <- ggplot(TPPD1.counts) + 
  geom_point(aes(x = stage, y = count, group = type2, color = type2), size = 1.5, shape = 21, stroke = 1)  + 
  geom_line(data = TPPD1.counts.mean, aes(x = stage, y = count, group = interaction(type2, fl2), color = type2)) +
  geom_segment(data = TPPD1.counts.mean, aes(x = stage-0.12, xend = stage+0.12, y = count, yend = count, group = interaction(type2, stage, fl2), color = type2)) +
  scale_color_manual(values = c("orange", 'darkmagenta')) +
  facet_wrap(~fl2, scales = 'free_y') + 
  labs(title = 'TPPD1 (CpaM1st31681)', color = '', x = 'Developmental stage', 
       y= 'Corrected transcript count') +
  theme_bw() + 
  theme(aspect.ratio = 1,
        plot.margin = margin(t = 10, r = 10, b = 10, l = 10, unit = 'pt'),
        strip.background = element_rect(fill = NA),
        plot.title = element_text(face = 'italic', size = 10))


# NDR1 CpaM1st31682
NDR1.counts <- plotCounts(dds.all, gene='CpaM1st31682', intgroup=c("type", "fl", 'stage'), pc = 0, returnData = T)

# reformat to get metadata back
NDR1.counts <- cbind(Run = rownames(NDR1.counts),data.table(NDR1.counts))
NDR1.counts[, Run := gsub("Qu_etal_2023//", "", Run)]
NDR1.counts <- merge(Qu23.meta, NDR1.counts, by = c("Run", "stage", "type", "fl"))
NDR1.counts[, stage:= as.numeric(as.character(stage))]
NDR1.counts.mean <- NDR1.counts[, .(count = mean(count)), by = .(stage, fl2, type2)]

NDR1_plt <- ggplot(NDR1.counts) + 
  geom_point(aes(x = stage, y = count, group = type2, color = type2), size = 1.5, shape = 21, stroke = 1)  + 
  geom_line(data = NDR1.counts.mean, aes(x = stage, y = count, group = interaction(type2, fl2), color = type2)) +
  geom_segment(data = NDR1.counts.mean, aes(x = stage-0.12, xend = stage+0.12, y = count, yend = count, group = interaction(type2, stage, fl2), color = type2)) +
  scale_color_manual(values = c("orange", 'darkmagenta')) +
  facet_wrap(~fl2, scales = 'free_y') + 
  labs(title = 'NDR1/HIN1-like (CpaM1st31682)', color = '', x = 'Developmental stage', 
       y= 'Corrected transcript count') +
  theme_bw() + 
  theme(aspect.ratio = 1,
        plot.margin = margin(t = 10, r = 10, b = 10, l = 10, unit = 'pt'),
        strip.background = element_rect(fill = NA),
        plot.title = element_text(face = 'italic', size = 10))


library(cowplot)
plot_grid(TPPD1_plt, NDR1_plt, ncol = 1, align = 'v')


# ------ Carya G-locus gene orthologs -----
# EMS1 CpaM1st35237
EMS1.counts <- plotCounts(dds.all, gene='CpaM1st35237', intgroup=c("type", "fl", 'stage'), pc = 0, returnData = T)

# reformat to get metadata back
EMS1.counts <- cbind(Run = rownames(EMS1.counts),data.table(EMS1.counts))
EMS1.counts[, Run := gsub("Qu_etal_2023//", "", Run)]
EMS1.counts <- merge(Qu23.meta, EMS1.counts, by = c("Run", "stage", "type", "fl"))
EMS1.counts[, stage:= as.numeric(as.character(stage))]
EMS1.counts.mean <- EMS1.counts[, .(count = mean(count)), by = .(stage, fl2, type2)]

EMS1_plt <- ggplot(EMS1.counts) + 
  geom_point(aes(x = stage, y = count, group = type2, color = type2), size = 1.5, shape = 21, stroke = 1)  + 
  geom_line(data = EMS1.counts.mean, aes(x = stage, y = count, group = interaction(type2, fl2), color = type2)) +
  geom_segment(data = EMS1.counts.mean, aes(x = stage-0.12, xend = stage+0.12, y = count, yend = count, group = interaction(type2, stage, fl2), color = type2)) +
  scale_color_manual(values = c("orange", 'darkmagenta')) +
  facet_wrap(~fl2, scales = 'free_y') + 
  labs(title = 'EMS1 (CpaM1st35237)', color = '', x = 'Developmental stage', 
       y= 'Corrected transcript count') +
  theme_bw() + 
  theme(aspect.ratio = 1,
        plot.margin = margin(t = 10, r = 10, b = 10, l = 10, unit = 'pt'),
        strip.background = element_rect(fill = NA),
        legend.position = 'none',
        plot.title = element_text(face = 'italic', size = 10))

EMS1_plt



# SLK2 CpaM1st35250
SLK2.counts <- plotCounts(dds.all, gene='CpaM1st35250', intgroup=c("type", "fl", 'stage'), pc = 0, returnData = T)

# reformat to get metadata back
SLK2.counts <- cbind(Run = rownames(SLK2.counts),data.table(SLK2.counts))
SLK2.counts[, Run := gsub("Qu_etal_2023//", "", Run)]
SLK2.counts <- merge(Qu23.meta, SLK2.counts, by = c("Run", "stage", "type", "fl"))
SLK2.counts[, stage:= as.numeric(as.character(stage))]
SLK2.counts.mean <- SLK2.counts[, .(count = mean(count)), by = .(stage, fl2, type2)]

SLK2_plt <- ggplot(SLK2.counts) + 
  geom_point(aes(x = stage, y = count, group = type2, color = type2), size = 1.5, shape = 21, stroke = 1)  + 
  geom_line(data = SLK2.counts.mean, aes(x = stage, y = count, group = interaction(type2, fl2), color = type2)) +
  geom_segment(data = SLK2.counts.mean, aes(x = stage-0.12, xend = stage+0.12, y = count, yend = count, group = interaction(type2, stage, fl2), color = type2)) +
  scale_color_manual(values = c("orange", 'darkmagenta')) +
  facet_wrap(~fl2, scales = 'free_y') + 
  labs(title = 'SLK2 (CpaM1st35250)', color = '', x = 'Developmental stage', 
       y= 'Corrected transcript count') +
  theme_bw() + 
  theme(aspect.ratio = 1,
        plot.margin = margin(t = 10, r = 10, b = 10, l = 10, unit = 'pt'),
        strip.background = element_rect(fill = NA),
        legend.position = 'none',
        plot.title = element_text(face = 'italic', size = 10))

SLK2_plt


# ---- Rhomboid-like CpaM1st35241
rhom.counts <- plotCounts(dds.all, gene='CpaM1st35241', intgroup=c("type", "fl", 'stage'), pc = 0, returnData = T)
# reformat to get metadata back
rhom.counts <- cbind(Run = rownames(rhom.counts),data.table(rhom.counts))
rhom.counts[, Run := gsub("Qu_etal_2023//", "", Run)]
rhom.counts <- merge(Qu23.meta, rhom.counts, by = c("Run", "stage", "type", "fl"))
rhom.counts[, stage:= as.numeric(as.character(stage))]
rhom.counts.mean <- rhom.counts[, .(count = mean(count)), by = .(stage, fl2, type2)]

rhom_plt <- ggplot(rhom.counts) + 
  geom_point(aes(x = stage, y = count, group = type2, color = type2), size = 1.5, shape = 21, stroke = 1)  + 
  geom_line(data = rhom.counts.mean, aes(x = stage, y = count, group = interaction(type2, fl2), color = type2)) +
  geom_segment(data = rhom.counts.mean, aes(x = stage-0.12, xend = stage+0.12, y = count, yend = count, group = interaction(type2, stage, fl2), color = type2)) +
  scale_color_manual(values = c("orange", 'darkmagenta')) +
  facet_wrap(~fl2, scales = 'free_y') + 
  labs(title = 'RHOMBOID-like (CpaM1st35241)', color = '', x = 'Developmental stage', 
       y= 'Corrected transcript count') +
  theme_bw() + 
  theme(aspect.ratio = 1,
        plot.margin = margin(t = 10, r = 10, b = 10, l = 10, unit = 'pt'),
        strip.background = element_rect(fill = NA),
        legend.position = 'none',
        plot.title = element_text(face = 'italic', size = 10))
rhom_plt

# ----- IQ Domain CpaM1st35244
iq.counts <- plotCounts(dds.all, gene='CpaM1st35244', intgroup=c("type", "fl", 'stage'), pc = 0, returnData = T)
# reformat to get metadata back
iq.counts <- cbind(Run = rownames(iq.counts),data.table(iq.counts))
iq.counts[, Run := gsub("Qu_etal_2023//", "", Run)]
iq.counts <- merge(Qu23.meta, iq.counts, by = c("Run", "stage", "type", "fl"))
iq.counts[, stage:= as.numeric(as.character(stage))]
iq.counts.mean <- iq.counts[, .(count = mean(count)), by = .(stage, fl2, type2)]

iq_plt <- ggplot(iq.counts) + 
  geom_point(aes(x = stage, y = count, group = type2, color = type2), size = 1.5, shape = 21, stroke = 1)  + 
  geom_line(data = iq.counts.mean, aes(x = stage, y = count, group = interaction(type2, fl2), color = type2)) +
  geom_segment(data = iq.counts.mean, aes(x = stage-0.12, xend = stage+0.12, y = count, yend = count, group = interaction(type2, stage, fl2), color = type2)) +
  scale_color_manual(values = c("orange", 'darkmagenta')) +
  facet_wrap(~fl2, scales = 'free_y') + 
  labs(title = 'IQ-DOMAIN 23-like (CpaM1st35244)', color = '', x = 'Developmental stage', 
       y= 'Corrected transcript count') +
  theme_bw() + 
  theme(aspect.ratio = 1,
        plot.margin = margin(t = 10, r = 10, b = 10, l = 10, unit = 'pt'),
        strip.background = element_rect(fill = NA),
        legend.position = 'none',
        plot.title = element_text(face = 'italic', size = 10))
iq_plt

# ------ generic for test  -------
# CpaM1st35251 ccb2
# CpaM1st35255 detox shows nice pattern
# CpaM1st35254 BAG
# CpaM1st35252 dihydro
# CpaM1st35249 At4g00755
# CpaM1st35248 Pp2-A12
# CpaM1st35247 RNA pol subunit 6A
# CpaM1st35246 PMT16
# CpaM1st35243 non-spec lipd transfer 14
# CpaM1st35242 LOC122306943 uncharacterized
# CpaM1st35240 Root primordium defective
# CpaM1st35239 XBAT33
# CpaM1st35238 CEN
# CpaM1st35235 LOC122306891 unchar
# CpaM1st35234 phos inhibitor 2


X.counts <- plotCounts(dds.all, gene='CpaM1st35255', intgroup=c("type", "fl", 'stage'), pc = 0, returnData = T)
# reformat to get metadata back
X.counts <- cbind(Run = rownames(X.counts),data.table(X.counts))
X.counts[, Run := gsub("Qu_etal_2023//", "", Run)]
X.counts <- merge(Qu23.meta, X.counts, by = c("Run", "stage", "type", "fl"))
X.counts[, stage:= as.numeric(as.character(stage))]
X.counts.mean <- X.counts[, .(count = mean(count)), by = .(stage, fl2, type2)]

X_plt <- ggplot(X.counts) + 
  geom_point(aes(x = stage, y = count, group = type2, color = type2), size = 1.5, shape = 21, stroke = 1)  + 
  geom_line(data = X.counts.mean, aes(x = stage, y = count, group = interaction(type2, fl2), color = type2)) +
  geom_segment(data = X.counts.mean, aes(x = stage-0.12, xend = stage+0.12, y = count, yend = count, group = interaction(type2, stage, fl2), color = type2)) +
  scale_color_manual(values = c("orange", 'darkmagenta')) +
  facet_wrap(~fl2, scales = 'free_y') + 
  labs(title = 'DETOXIFICATION 49-like (CpaM1st35255)', color = '', x = 'Developmental stage', 
       y= 'Corrected transcript count') +
  theme_bw() + 
  theme(aspect.ratio = 1,
        plot.margin = margin(t = 10, r = 10, b = 10, l = 10, unit = 'pt'),
        strip.background = element_rect(fill = NA),
        legend.position = 'none',
        plot.title = element_text(face = 'italic', size = 10))
X_plt




plot_grid(EMS1_plt, X_plt, ncol = 1)









# # FIL1 ortholog CpaM1st23605
# Qu23.counts.FIL1 <- Qu23.counts_pr[gene == 'CpaM1st23605']
# 
# ggplot(Qu23.counts.FIL1) + 
#   geom_line(aes(x = stage, y = count, group = interaction(type2, rep), color = type2))  + 
#   geom_point(aes(x = stage, y = count, group = type2, color = type2))  + 
#   scale_color_manual(values = c("orange", 'darkmagenta')) +
#   facet_wrap(~fl2, scales = 'free_y') + 
#   labs(title = 'FIL1', x = 'Developmental stage', 
#        color = '',
#        y= 'Estimated transcript count') +
#   theme_bw() + 
#   theme(aspect.ratio = 1,
#         plot.margin = margin(t = 10, r = 10, b = 10, l = 10, unit = 'pt'),
#         strip.background = element_rect(fill = NA),
#         plot.title = element_text(face = 'italic'))





# -------- Chen et al 2019 -----
Chen.meta <- fread("Chen_etal_2019/sra.txt", header = F, col.names = c("Run", "sample"))

# transcript quantifications, annotation lacking both paralogs
Chen.quants_files <- list.files("Chen_etal_2019/quants/", recursive = T, pattern = ".*quant\\.sf$", full.names = T)
Chen.all_names <- gsub("_quant/quant.sf", "", gsub("quants//", "", Chen.quants_files))
names(Chen.quants_files) <- Chen.all_names

tx2gene <- fread("Chen_etal_2019/tx2gene.txt", header = F, col.names = c("tx", "gene"))

# The tximport package has a single function for importing transcript-level estimates. 
# The type argument is used to specify what software was used for estimation. 
# A simple list with matrices, "abundance", "counts", and "length", is returned, 
# where the transcript level information is summarized to the gene-level.

Chen.txi <- tximport(Chen.quants_files, type = "salmon", tx2gene = tx2gene)
head(Chen.txi$counts)

Chen.counts <- data.table(Chen.txi$counts)
Chen.counts[, gene := rownames(Chen.txi$counts)]
Chen.counts <- melt(Chen.counts, id.vars = 'gene', variable.name = 'Run', value.name = 'count')
Chen.counts[, Run := gsub("Chen_etal_2019/", "", Run)]
Chen.counts <- merge(Chen.counts, Chen.meta)

Chen.counts[, c("type", "fl_rep") := tstrsplit(sample, split = "-")]
Chen.counts[, c("fl", "rep") := tstrsplit(fl_rep, split = "")]

Chen.counts.FAFL <- Chen.counts[ gene == 'CpaM1st14242']
Chen.counts.FAFL

# --- with paralogs --- 

# transcript quantifications, annotation lacking both paralogs
Chen.quants_files_pr <- list.files("Chen_etal_2019/quants_paralogs/", recursive = T, pattern = ".*quant\\.sf$", full.names = T)
Chen.all_names_pr <- gsub("_quant/quant.sf", "", gsub("quants_paralogs/", "", Chen.quants_files_pr))
names(Chen.quants_files_pr) <- Chen.all_names_pr

tx2gene_pr <- fread("Chen_etal_2019/tx2gene_paralogs.txt", header = F, col.names = c("tx", "gene"))

Chen.txi_pr <- tximport(Chen.quants_files_pr, type = "salmon", tx2gene = tx2gene_pr, countsFromAbundance = 'no')
head(Chen.txi_pr$counts)

Chen.counts_pr <- data.table(Chen.txi_pr$counts)
Chen.counts_pr[, gene := rownames(Chen.txi_pr$counts)]
Chen.counts_pr <- melt(Chen.counts_pr, id.vars = 'gene', variable.name = 'Run', value.name = 'Count')

Chen.counts_pr[, Run := gsub("Chen_etal_2019//", "", Run)]
Chen.counts_pr <- merge(Chen.counts_pr, Chen.meta)
Chen.counts_pr[, c("type", "fl_rep") := tstrsplit(sample, split = "-")]
Chen.counts_pr[, c("fl", "rep") := tstrsplit(fl_rep, split = "")]

Chen.counts.FAFL_pr <- Chen.counts_pr[grepl("GFAFL", gene)]
Chen.counts.FAFL_pr




# expression of GFAFL2 is zero in these samples.
Chen.counts.FAFL_pr[type == 'PG' & gene == '2PA_GFAFL2', Count := NA]

Chen.counts.FAFL_pr[type == 'PA', type2 := 'Protandrous']
Chen.counts.FAFL_pr[type == 'PG', type2 := 'Protogynous']
Chen.counts.FAFL_pr[,gene := gsub("2PA_", "", gene)]


ggplot(Chen.counts.FAFL_pr) + 
  geom_bar(aes(x = rep, y = Count, fill = gene), stat = 'identity') +
  scale_fill_manual(values = c("blue", 'red')) +
  facet_grid(fl~type2, scales = 'free_y') + 
  labs(x = 'Replicate', fill = 'allele') +
  theme_bw() + 
  theme(aspect.ratio = 1)
  



# -------- Zhang et al. 2024 -----
Zhang.meta <- fread("Zhang_etal_2024/meta.txt", header = F, col.names = c("Run", "sample"))

# transcript quantifications
Zhang.quants_files <- list.files("Zhang_etal_2024/quants/", recursive = T, pattern = ".*quant\\.sf$", full.names = T)
Zhang.all_names <- gsub("_quant/quant.sf", "", gsub("quants//", "", Zhang.quants_files))
names(Zhang.quants_files) <- Zhang.all_names

Zhang.txi <- tximport(Zhang.quants_files, type = "salmon", tx2gene = tx2gene, countsFromAbundance = 'scaledTPM')
Zhang.counts <- data.table(Zhang.txi$counts)
Zhang.counts[, gene := rownames(Zhang.txi$counts)]
Zhang.counts <- melt(Zhang.counts, id.vars = 'gene', variable.name = 'Run', value.name = 'TPM')
Zhang.counts[, Run := gsub("Zhang_etal_2024/", "", Run)]
Zhang.counts <- merge(Zhang.counts, Zhang.meta)
Zhang.counts.FAFL <- Zhang.counts[ gene == 'CpaM1st14242']
Zhang.counts.FAFL

Zhang.counts.FAFL[, c("tissue", "replicate") := tstrsplit(sample, ".replicate")]

ggplot(Zhang.counts.FAFL, aes(x = replicate, y = TPM)) + 
  geom_bar(stat = 'identity') + 
  facet_wrap(~tissue, nrow = 1) + 
  labs(x = "Replicate", y = 'Transcripts per million') +
  theme_classic() + 
  theme(aspect.ratio = 1)




# -------- Qu et al. 2021 -----
Qu21.meta <- fread("Qu_etal_2021/meta.txt", header = F, col.names = c("Run", "sample"))

# transcript quantifications
Qu21.quants_files <- list.files("Qu_etal_2021/quants/", recursive = T, pattern = ".*quant\\.sf$", full.names = T)
Qu21.all_names <- gsub("_quant/quant.sf", "", gsub("quants//", "", Qu21.quants_files))
names(Qu21.quants_files) <- Qu21.all_names

Qu21.txi <- tximport(Qu21.quants_files, type = "salmon", tx2gene = tx2gene, countsFromAbundance = 'no')
Qu21.counts <- data.table(Qu21.txi$counts)
Qu21.counts[, gene := rownames(Qu21.txi$counts)]
Qu21.counts <- melt(Qu21.counts, id.vars = 'gene', variable.name = 'Run', value.name = 'count')
Qu21.counts[, Run := gsub("Qu_etal_2021/", "", Run)]
Qu21.counts <- merge(Qu21.counts, Qu21.meta)

Qu21.counts[, c("stage", "type", "fl", "rep") := tstrsplit(sample, split = "")]

Qu21.counts.FAFL <- Qu21.counts[ gene == 'CpaM1st14242']
Qu21.counts.FAFL



# -------- match A-B to PA-PG -----
pcdat.long <- rbind(Qu23.counts_pr[stage == 2, .(Run, sample, type = type, fl, gene, log.count=log10(count+0.5))], 
      Chen.counts_pr[, .(Run, sample, type = type, fl, gene, log.count = log10(Count+0.5))])


pcdat.long[, type := paste0(type, "_", fl )]

pcdat<-  dcast(pcdat.long, Run + sample + type ~ gene, value.var = 'log.count')

pca <- prcomp(pcdat[, -(1:4)])
pc12 <- cbind(pcdat[, 1:4], pca$x[, 1:2])

library(pals)
ggplot(pc12, aes(x = PC1, y = PC2, color = type)) +
  geom_point(size = 3) + 
  scale_color_manual(values = tol(8)) + 
  labs(color = 'metadata') +
  theme_classic() + 
  theme(aspect.ratio = 1)


# Chen 2019 vs Qu 2021
pcdat2.long <- rbind(Qu21.counts[stage == 2, .(Run, sample, type, fl, gene, log.count=log10(count+0.5))], 
                    Chen.counts[, .(Run, sample, type, fl, gene, log.count = log10(count+0.5))])

pcdat2 <-  dcast(pcdat2.long, Run + sample + type ~ gene, value.var = 'log.count')

pca2 <- prcomp(pcdat2[, -(1:4)])
pc12.2 <- cbind(pcdat2[, 1:4], pca2$x[, 1:2])

pc12.2[PC1 > 0]
ggplot(pc12.2, aes(x = PC1, y = PC2, color = type)) +
  geom_point(size = 3) + 
  scale_color_manual(values = tol(4)) + 
  labs(color = 'metadata') +
  theme_classic() + 
  theme(aspect.ratio = 1)



