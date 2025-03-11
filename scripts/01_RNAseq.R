library(data.table)
library(tximport)
library(forcats)
library("DESeq2")
library(ggplot2)
library(tidyverse)
library(ggpubfigs)
library(lubridate) # for manipulating POSIX
library(ggpattern)

setwd("~/workspace/Pterocarya/RNAseq/")

# ---- read phenotypes -----
meta <- fread("~/workspace/Pterocarya/Pste_WGS_phenotypes.txt", col.names = c("sample", "pheno"))
meta[pheno==0, phenotype := 'male-first']
meta[pheno==1, phenotype := 'female-first']

meta_RNAseq <- fread("meta.txt")
meta_RNAseq[Type == 'female-first' & ID != "Pste_WS_17.01_mbuds", geno := 'Gg']
meta_RNAseq[ID == "Pste_WS_17.01_mbuds", geno := 'GG']
meta_RNAseq[ID == "DV_145", geno := 'Gg']
meta_RNAseq[Type == 'male-first', geno := 'gg']


# ----- read transcript quantification files -----
# this is from a competitive mapping approach, where both alleles of SDMYB were included

# get names of quant files, all samples
quants_files <- list.files("quants", recursive = T, pattern = ".*quant\\.sf$", full.names = T)
all_names <- gsub("_quant/quant.sf", "", gsub("quants/", "", quants_files))
names(quants_files) <- all_names

quants_files

# files that are specifically just male bud tissue
mbuds_quant_files <- quants_files[names(quants_files) %in% meta_RNAseq[Tissue == 'male buds pre-dormancy', ID]]

# files that are specifically just male bud tissue
f_quant_files <- quants_files[names(quants_files) %in% meta_RNAseq[grepl("female", Tissue), ID]]

# files that are specifically just mature male pre-anthesis
mm_quant_files <- quants_files[names(quants_files) %in% meta_RNAseq[Tissue == 'mature male pre-anthesis', ID]]

# gene to species mapping
tx2gene <- fread("tx2gene_competitive.txt", header = F, col.names = c("tx", "gene"))

# tximport
txi.all <- tximport(quants_files, type = "salmon", tx2gene = tx2gene, countsFromAbundance = 'scaledTPM')
head(txi.all$counts)

txi.mbuds <- tximport(mbuds_quant_files, type = "salmon", tx2gene = tx2gene, countsFromAbundance = 'scaledTPM')
txi.f <- tximport(f_quant_files, type = "salmon", tx2gene = tx2gene, countsFromAbundance = 'scaledTPM')
txi.mm <- tximport(mm_quant_files, type = "salmon", tx2gene = tx2gene, countsFromAbundance = 'scaledTPM')

# gene counts
counts <- data.table(txi.all$counts)
counts[, gene := rownames(txi.all$counts)]

# ------ PCA, all samples -----
# *** critical to make rowns of metadata match columns of count matrix
# From DESeqDataSetFromTximport help page:
# Rows of colData correspond to columns of countData
all_order_indices <- match(names(as.data.table(txi.all$counts)), meta_RNAseq$ID)
meta_RNAseq.all <- meta_RNAseq[all_order_indices]

# 
ddsTxi.all <- DESeqDataSetFromTximport(txi.all,
                                       colData = meta_RNAseq.all,
                                       design = ~ Tissue + Type)


dds.all <- DESeq(ddsTxi.all)
vsd.all <- vst(dds.all)
#

#
plotPCA(vsd.all, intgroup=c("Type", "Tissue"), returnData = F)

pca_data_all <- setDT(plotPCA(vsd.all, intgroup=c("Type", "Tissue"), returnData = T))

pca_data_all[Type == 'female-first', geno := 'Gg']
pca_data_all[name == 'Pste_WS_17.01_mbuds', geno := 'GG']
pca_data_all[Type == 'male-first', geno := 'gg']
pca_data_all[name == 'DV_145', geno := 'Gg']

ggplot(pca_data_all, aes(x = PC1, y = PC2)) + 
  geom_point(aes(shape = Tissue, color = geno), 
             size= 1.5, stroke = 1) + 
  labs(x = 'PC1: 75% variance', y = 'PC2: 16% variance', color = 'Genotype') +
  scale_shape_manual(values = c(0,1,2,5),
                     labels = c("immature female", 
                                "male buds pre-dormancy", 
                                "mature female", 
                                "mature male pre-anthesis")) +
  scale_color_manual(
    values = c('gg' = '#DDAA33', 'Gg' = '#004488', 
               'GG' = '#DC267F'),
    guide =
      guide_legend(label.theme = element_text(angle = 0, face = "italic"))) +
              #'PSTE_W1_5' = 'turquoise1', 'PSTE_DV_134' = '#DC267F')
  theme_classic() + 
  theme(aspect.ratio = 1) 

pca_data_all


# see what naive pca of raw counts looks like

# all.raw.cnts <- dcast(melt(counts, id.vars = 'gene', value.name = 'counts', variable.name = 'sample'),
#                          sample ~ gene, value.var = 'counts')
# 
# all.pca.raw.cnts <- cbind(all.raw.cnts[,1], prcomp(all.raw.cnts[,-1])$x)
# 
# all.pca.raw.cnts <- merge(all.pca.raw.cnts, meta_RNAseq[, .(sample = ID, Tissue, geno, Type)])
# 
# # naive pca of counts shows broadly similar pattern, albeit somewhat weaker clustering by tissue
# # the DeSeq PCA seems more sensible
# ggplot(all.pca.raw.cnts, aes(x = PC1, y = PC2, color = geno, shape = Tissue))+ geom_point()
# 


# ----- PCA, just male buds samples -----
mbuds_order_indices <- match(names(as.data.table(txi.mbuds$counts)), meta_RNAseq$ID)
meta_RNAseq.mbuds <- meta_RNAseq[mbuds_order_indices]

ddsTxi.mbuds <- DESeqDataSetFromTximport(txi.mbuds,
                                       colData = meta_RNAseq.mbuds,
                                       design = ~ Type)

dds.mbuds <- DESeq(ddsTxi.mbuds)
vsd.mbuds <- vst(dds.mbuds)
plotPCA(vsd.mbuds, intgroup=c("Type"))

pca_data_mbuds <- setDT(plotPCA(vsd.mbuds, intgroup=c("Type", "Tissue"), returnData = T))


pca_data_mbuds[Type == 'female-first', geno := 'Gg']
pca_data_mbuds[name == 'Pste_WS_17.01_mbuds', geno := 'GG']
pca_data_mbuds[Type == 'male-first', geno := 'gg']
pca_data_mbuds[name == 'DV_145', geno := 'Gg']

mbuds <- ggplot(pca_data_mbuds, aes(x = PC1, y = PC2)) + 
  geom_point(aes(color = geno), size= 2, stroke = 1, shape =1) + 
  labs(x = 'PC1: 26% variance', 
       y = 'PC2: 15% variance', color = 'Genotype') +
  scale_color_manual(
    values = c('gg' = '#DDAA33', 'Gg' = '#004488', 
               'GG' = '#DC267F')) +
    
  #'PSTE_W1_5' = 'turquoise1', 'PSTE_DV_134' = '#DC267F')
  theme_classic() + 
  theme(aspect.ratio = 1,
        plot.title = element_text(size = 8),
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 8),
        legend.position = 'none') 
mbuds



# ------ PCA using ascertained genes that are SCOs to DEGs in Jhin OR Jreg -----
# Juglans_DEG_Orthologs <- fread("~/workspace/Pterocarya/RNAseq/Pste_candidate_pollen_genes.txt", col.names = 'gene')
# 
# mbuds.counts <- data.table(txi.mbuds$counts)
# mbuds.counts[, gene := rownames(txi.mbuds$counts)]
# 
# mbuds.counts.Jug.DEG <- mbuds.counts[gene %in% Juglans_DEG_Orthologs$gene]
# select.cols <- names(mbuds.counts.Jug.DEG)[names(mbuds.counts.Jug.DEG) != 'gene']
# 
# mbuds.counts.Jug.DEG <- cbind(mbuds.counts.Jug.DEG[, gene],
#       mbuds.counts.Jug.DEG[, lapply(.SD, function(x){log(x+0.5)}), 
#                            .SDcols = select.cols])
# setnames(mbuds.counts.Jug.DEG, 'V1', 'gene')
# 
# mbuds.counts.Jug.DEG <- dcast( melt(mbuds.counts.Jug.DEG, 
#                                    id.vars = 'gene', value.name = 'log.counts', 
#                                    variable.name = 'sample'),
#       sample ~ gene, value.var = 'log.counts')
# 
# 
# 
# 
# mbuds.Jug.DEG.PC <- cbind(mbuds.counts.Jug.DEG[, 1],prcomp(mbuds.counts.Jug.DEG[,-1])$x)
# mbuds.Jug.DEG.PC
# 
# mbuds.Jug.DEG.PC[sample %in% meta_RNAseq[geno == 'gg', ID], geno := 'gg']
# mbuds.Jug.DEG.PC[sample %in% meta_RNAseq[geno == 'Gg', ID], geno := 'Gg']
# mbuds.Jug.DEG.PC[sample %in% meta_RNAseq[geno == 'GG', ID], geno := 'GG']
# 
# ggplot(mbuds.Jug.DEG.PC, aes(x = PC1, y = PC2, color = geno))+ geom_point()
# 
# mbuds.pca.sub <- cbind(mbuds.sub[,1], prcomp(mbuds.sub[,-1])$x, center = T, scale. = T)
# mbuds.pca.sub[sample %in% meta_RNAseq[geno == 'gg', ID], geno := 'gg']
# mbuds.pca.sub[sample %in% meta_RNAseq[geno == 'Gg', ID], geno := 'Gg']
# mbuds.pca.sub[sample %in% meta_RNAseq[geno == 'GG', ID], geno := 'GG']
# 
# ggplot(mbuds.pca.sub, aes(x = PC1, y = PC2, color = geno))+ geom_point()




# ----- PCA, just female samples -----
f_order_indices <- match(names(as.data.table(txi.f$counts)), meta_RNAseq$ID)
meta_RNAseq.f <- meta_RNAseq[f_order_indices]

ddsTxi.f <- DESeqDataSetFromTximport(txi.f,
                                         colData = meta_RNAseq.f,
                                         design = ~ Type)

dds.f <- DESeq(ddsTxi.f)
vsd.f <- vst(dds.f, blind = FALSE)
plotPCA(vsd.f, intgroup=c("Type"))

pca_data_f <- setDT(plotPCA(vsd.f, intgroup=c("Type", "Tissue"), returnData = T))

pca_data_f[Type == 'female-first', geno := 'Gg']
pca_data_f[name == 'Pste_WS_17.01_mbuds', geno := 'GG']
pca_data_f[Type == 'male-first', geno := 'gg']
pca_data_f[name == 'DV_145', geno := 'Gg']

fm <- ggplot(pca_data_f, aes(x = PC1, y = PC2)) + 
  geom_point(aes(shape = Tissue, color = geno), size= 2, stroke = 1) + 
  scale_shape_manual(values = c(0,2),
                     labels = c("immature female", 
                                "mature female")) +
  labs(x = 'PC1: 50% variance', 
       y = 'PC2: 28% variance', color = 'Genotype') +
  scale_color_manual(
    values = c('gg' = '#DDAA33', 'Gg' = '#004488', 
               'GG' = '#DC267F'),
    guide = guide_legend(label.theme = element_text(angle = 0, face = "italic"))) +
  #'PSTE_W1_5' = 'turquoise1', 'PSTE_DV_134' = '#DC267F')
  theme_classic() + 
  theme(aspect.ratio = 1,
        plot.title = element_text(size = 8),
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 8),
        legend.position = 'none') 

fm


# ----- PCA, just mature male -----
mm_order_indices <- match(names(as.data.table(txi.mm$counts)), meta_RNAseq$ID)
meta_RNAseq.mm <- meta_RNAseq[mm_order_indices]

ddsTxi.mm <- DESeqDataSetFromTximport(txi.mm,
                                     colData = meta_RNAseq.mm,
                                     design = ~ Type)

dds.mm <- DESeq(ddsTxi.mm)
vsd.mm <- vst(dds.mm)
plotPCA(vsd.mm, intgroup=c("Type"))


pca_data_mm <- setDT(plotPCA(vsd.mm, intgroup=c("Type", "Tissue"), returnData = T))

pca_data_mm[Type == 'female-first', geno := 'Gg']
pca_data_mm[name == 'Pste_WS_17.01_mbuds', geno := 'GG']
pca_data_mm[Type == 'male-first', geno := 'gg']
pca_data_mm[name == 'DV_145', geno := 'Gg']

mm <- ggplot(pca_data_mm, aes(x = PC1, y = PC2)) + 
  geom_point(aes(color = geno), size= 2, shape = 5, stroke = 1) + 
  labs(x = 'PC1: 74% variance', 
       y = 'PC2: 16% variance', color = 'Genotype') +
  scale_color_manual(
    values = c('gg' = '#DDAA33', 'Gg' = '#004488', 
               'GG' = '#DC267F')) + 
  #'PSTE_W1_5' = 'turquoise1', 'PSTE_DV_134' = '#DC267F')
  theme_classic() + 
  theme(aspect.ratio = 1,
        plot.title = element_text(size = 8),
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 8),
        legend.position = 'none') 


mm

plot_grid(mbuds, mm, fm, ncol = 3)



# ----- select top DEGs ----
res.mbuds <- results(dds.mbuds)
topn <- 10
rownames(head(res.mbuds[order(res.mbuds$padj),], topn))

res.mbuds[grepl("g21834", rownames(res.mbuds)),]








# ---- FAFL gene counts -----
# extract focal gene and reformat

# from DESeq estimated counts
FAFL_hap1_cnts.all <- plotCounts(dds.all, gene='P.stenopera_UCD1470_hap1.g21834.FAFL', intgroup=c("Type", 'Tissue'), returnData = T, pc = 0)
FAFL_hap1_cnts.all$sample.ID <- rownames(FAFL_hap1_cnts.all)
setDT(FAFL_hap1_cnts.all)
FAFL_hap1_cnts.all[, allele := 'g']

FAFL_hap2_cnts.all <- plotCounts(dds.all, gene='P.stenoptera_UCD1470_hap2.g21575.FAFL', intgroup=c("Type", "Tissue"), returnData = T, pc = 0) 
FAFL_hap2_cnts.all$sample.ID <- rownames(FAFL_hap2_cnts.all)
setDT(FAFL_hap2_cnts.all)
FAFL_hap2_cnts.all[, allele := 'G']

FAFL_counts.all <- rbind(FAFL_hap1_cnts.all, FAFL_hap2_cnts.all)

# using counts from scaled TPM
# dom <- melt(counts[grepl("hap2", gene)], id.vars = 'gene', value.name = 'count', variable.name = 'sample.ID')
# rec <- melt(counts[grepl("hap1", gene)], id.vars = 'gene', value.name = 'count', variable.name = 'sample.ID')
# dom[, gene := NULL][, allele := 'G'][]
# rec[, gene := NULL][, allele := 'g'][]
# 
# candidate <- rbind(dom, rec)
# candidate <- merge(candidate, meta_RNAseq[, .(sample.ID = ID, Type, Tissue)], by = 'sample.ID')
# candidate


# show all
ggplot(FAFL_counts.all, aes(x = sample.ID)) +
  facet_grid(Tissue~Type, scales = 'free_x') + 
  geom_bar(stat = 'identity', aes(y = count, fill = allele)) + 
  labs(x = '') + 
  theme_classic() + 
  theme(axis.text.x = element_blank())

# ---- plot only tissues where FAFL is expressed ---

ggplot(FAFL_counts.all[Tissue %in% c("male buds pre-dormancy", "mature male pre-anthesis") &
                   Type != 'ambiguous'], aes(x = sample.ID)) +
  facet_grid(Tissue~Type, scales = 'free_x') + 
  geom_bar(stat = 'identity', aes(y = count, fill = allele)) + 
  labs(x = '') + 
  theme_classic() + 
  theme(axis.text.x = element_blank())

#candidate$Type <- as.factor(candidate$Type, levels = c("male-first", "female-first", "ambiguous"))
#candidate$sample.ID <- fct_reorder2(candidate$sample.ID, candidate$Type, candidate$count)

# reorder IDs
FAFL_counts.all$sample.ID <- factor(FAFL_counts.all$sample.ID, 
                                 levels = c(
                                   # first by tissue, then type then order
                                   # PA
                                   "Pste_WS_10.01_mbuds", "Pste_WS_10.02_mbuds","PSTE_WS_2.10","Pste_WS_10.05_mbuds",
                                   "Pste_WS_02.06_mbuds",'PSTE_WS_2.08',
                                   #PG , increasing order
                                   "PSTE_UCD1", "PSTE_WS_11.01","Pste_WS_01.01_mbuds",
                                   "Pste_WS_17.01_mbuds","Pste_WS_01.05_mbuds",
                              
                                   # mature male
                                   'DV_136', "DV_145", "DV_149_M", 'DV_138_M' ))

FAFL_counts.all[Tissue == "mature male pre-anthesis", tissuelabel := 'mature male fls']
FAFL_counts.all[Tissue == "male buds pre-dormancy", tissuelabel := 'pre-dormant male fls']
FAFL_counts.all$tissuelabel <- factor(FAFL_counts.all$tissuelabel, levels = c("pre-dormant male fls", "mature male fls"))

FAFL_counts.all[Type == 'male-first', geno := 'gg']
FAFL_counts.all[Type == 'female-first', geno := 'Gg']
FAFL_counts.all[grepl('17', sample.ID), geno := 'GG']
FAFL_counts.all[grepl('DV_145', sample.ID), geno := 'Gg']


FAFL_counts.all
FAFL_counts.ml <- FAFL_counts.all[Tissue %in% c("male buds pre-dormancy", "mature male pre-anthesis")]
FAFL_counts.mbud <- FAFL_counts.ml[Tissue %in% c("male buds pre-dormancy")]
ggplot(FAFL_counts.ml, 
       aes(x = sample.ID, fill = allele)) + 
  # bar chart
  geom_bar(stat = 'identity', position = 'stack', aes(y = count)) + 
  facet_wrap(~tissuelabel, scales = 'free', ncol = 1) + 
  
  # customize labels, colors 
  labs(x = '', y = 'Normalized transcript count') +
  scale_fill_manual(values = c(#'g' = '#DDCC10', 
                                 'g' = 'lightsalmon', 
                                 'G' = 'darkmagenta')) +
  scale_x_discrete(labels = c("Pste_WS_10.01_mbuds" = "gg", 
                              "Pste_WS_10.02_mbuds" = "gg",
                              "PSTE_WS_2.10" = "gg",
                              "Pste_WS_10.05_mbuds" = "gg",
                              "Pste_WS_02.06_mbuds" = "gg",
                              'PSTE_WS_2.08' = "gg",
                              "PSTE_UCD1" = "Gg",  
                              "PSTE_WS_11.01" = "Gg",
                              "Pste_WS_01.01_mbuds" = "Gg",
                              "Pste_WS_17.01_mbuds" = "GG",
                              "Pste_WS_01.05_mbuds" = 'Gg',
                              'DV_136' = "gg", "DV_145" = "Gg", 
                              "DV_149_M" = "Gg", 'DV_138_M' ="Gg"
                               )) +
  # customize theme
  labs(fill = 'Expressed\nallele') +
  theme_classic() + 
  theme(aspect.ratio = 1,
        legend.text = element_text(face = 'italic'),
        axis.text.x = element_text(angle = 90, vjust = 0.5, face = 'italic'),
        legend.title = element_text(size = 9))
     





# ------ P values -------
# 1. allele-specific expression in buds of heterozygotes



# normalized counts from DESeq2
x <- FAFL_counts.mbud[Type == 'female-first' & 
                  !grepl("17.01", sample.ID) & 
                  allele == 'g',
                count, by = sample.ID]
x

y <- FAFL_counts.mbud[Type == 'female-first' & 
                       Tissue == 'male buds pre-dormancy' &
                       !grepl("17.01", sample.ID) & 
                       allele == 'G',
                     count, by = sample.ID]

x
y
# round down to be conservative
prop.test(x = 0, n = 9)
prop.test(x = 0, n = 16.66)

# counts from salmon (to check that the size correction counts from DESeq2 are not biasing result)
counts[grepl("FAF", gene), 'Pste_WS_01.01_mbuds']
counts[grepl("FAF", gene), 'Pste_WS_01.05_mbuds']
# for this sample, the depth correction is considerably lowering the transcript count, 
# so the result is actually much more significant. 



# mature fls
x1 <- FAFL_counts.all[Type %in% c('female-first', 'ambiguous') &
                       Tissue == 'mature male pre-anthesis' &
                       !grepl("17.01", sample.ID) & 
                       allele == 'g',
                     count]
y1 <- FAFL_counts.all[Type %in% c('female-first', 'ambiguous') & 
                       Tissue == 'mature male pre-anthesis' &
                       !grepl("17.01", sample.ID) & 
                       allele == 'G',
                     count]

prop.test(x = x1[1], n = x1[1] + y1[1])
prop.test(x = x1[2], n = x1[2] + y1[2])
prop.test(x = x1[3], n = x1[3] + y1[3])


# 2. test difference in g expression between homozygotes and heterozygotes in buds
s.gg <- FAFL_counts.mbud[geno == 'gg' & allele == 'g', count/2]
s.Gg <- FAFL_counts.mbud[geno == 'Gg' & allele == 'g', count]

library(MASS)  # For glm.nb()
library(stats)

# Example data: counts and group labels
set.seed(123)
group <- factor(c(rep("gg", length(s.gg)), rep('Gg', length(s.Gg))))
counts <- c(s.gg, s.Gg)

# Fit a Negative Binomial model
nb.model <- glm.nb(counts ~ group)

# Perform Likelihood Ratio Test (LRT)
nb.reduced <- glm.nb(counts ~ 1)  # Reduced model without group effect
anova(nb.model, nb.reduced, test = "Chisq")



# fit poisson model
p.model <- glm(counts ~ group, family = poisson())

# Perform Likelihood Ratio Test (LRT)
p.reduced <-glm(counts ~ 1, family = poisson())
anova(p.model, p.reduced, test = "Chisq")




# ---- TPPD1, EMS1, FIL ortholog gene counts -----
# extract focal gene and reformat
#g11858

# from DESeq estimated counts
TPPD1_cnts.all <- plotCounts(dds.all, gene='Pste.g11858', intgroup=c("Type", 'Tissue'), returnData = T, pc = 0)
TPPD1_cnts.all$sample.ID <- rownames(TPPD1_cnts.all)
setDT(TPPD1_cnts.all)

ggplot(TPPD1_cnts.all, aes(x = sample.ID, y = count, fill = Type)) + 
  geom_bar(stat = 'identity')+ 
  facet_wrap(~Tissue)
res.mbuds[rownames(res.mbuds) == 'Pste.g11858',]

# check correlation between FAFL expression and orthologs (answer: none)
FAFL_counts_combined <- FAFL_counts.all[Tissue == 'male buds pre-dormancy', .(FAFL_count = sum(count)), by = .(Type, Tissue, sample.ID)]
TPPD1_mbuds <- TPPD1_cnts.all[Tissue == 'male buds pre-dormancy', .(TPPD1_count = count, sample.ID)]
ggplot(merge(FAFL_counts_combined[Type == 'female-first'], TPPD1_mbuds, by = 'sample.ID'), 
       aes(x = FAFL_count, y = TPPD1_count)) + geom_point()

EMS1_cnts.all <- plotCounts(dds.all, gene='Pste.g13028', intgroup=c("Type", 'Tissue'), returnData = T, pc = 0)
EMS1_cnts.all$sample.ID <- rownames(EMS1_cnts.all)
setDT(EMS1_cnts.all)
ggplot(EMS1_cnts.all, aes(x = sample.ID, y = count, fill = Type)) + 
  geom_bar(stat = 'identity')+ 
  facet_wrap(~Tissue, scales = 'free_y')
res.mbuds[rownames(res.mbuds) == 'Pste.g13028',]

# EMS1_mbuds <- EMS1_cnts.all[Tissue == 'male buds pre-dormancy', .(EMS1_count = count, sample.ID)]
# ggplot(merge(FAFL_counts_combined, EMS1_mbuds, by = 'sample.ID'), 
#        aes(x = FAFL_count, y = EMS1_count)) + geom_point()


FIL1_cnts.all <- plotCounts(dds.all, gene='Pste.g13037', intgroup=c("Type", 'Tissue'), returnData = T, pc = 0)
FIL1_cnts.all$sample.ID <- rownames(FIL1_cnts.all)
setDT(FIL1_cnts.all)
ggplot(FIL1_cnts.all, aes(x = sample.ID, y = count, fill = Type)) + 
  geom_bar(stat = 'identity')+ 
  facet_wrap(~Tissue, scales = 'free_y')
res.mbuds[rownames(res.mbuds) == 'Pste.g13037',]


# ------- Plot other genes in region ------
# from DESeq estimated counts
g21833 <- plotCounts(dds.all, gene='Pste.g21833', intgroup=c("Type", 'Tissue'), returnData = T, pc = 0)
g21833$sample.ID <- rownames(g21833)
setDT(g21833)

ggplot(g21833[Type != 'ambiguous'], aes(x = sample.ID)) +
  facet_grid(Tissue~Type, scales = 'free_x') + 
  geom_bar(stat = 'identity', aes(y = count)) + 
  labs(x = '') + 
  theme_classic() + 
  theme(axis.text.x = element_blank())


g21835 <- plotCounts(dds.all, gene='Pste.g21835', intgroup=c("Type", 'Tissue'), returnData = T, pc = 0)
g21835$sample.ID <- rownames(g21835)
setDT(g21835)
ggplot(g21835[Type != 'ambiguous'], aes(x = sample.ID)) +
  facet_grid(Tissue~Type, scales = 'free_x') + 
  geom_bar(stat = 'identity', aes(y = count)) + 
  labs(x = '') + 
  theme_classic() + 
  theme(axis.text.x = element_blank())

g21836 <- plotCounts(dds.all, gene='Pste.g21836', intgroup=c("Type", 'Tissue'), returnData = T, pc = 0)
g21836$sample.ID <- rownames(g21836)
setDT(g21836)
ggplot(g21836[Type != 'ambiguous'], aes(x = sample.ID)) +
  facet_grid(Tissue~Type, scales = 'free_x') + 
  geom_bar(stat = 'identity', aes(y = count)) + 
  labs(x = '') + 
  theme_classic() + 
  theme(axis.text.x = element_blank())







