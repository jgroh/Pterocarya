library(data.table)
library(tidyverse)
library(pals)
library(ggpubfigs)

# read in phenotypes
meta <- fread("~/workspace/Pterocarya/Pste_WGS_phenotypes.txt", col.names = c("ID", "Type"))
meta[, Type := as.character(Type)]
meta[Type == 0, Type := 'protandrous']
meta[Type == 1, Type := 'protogynous']
meta <- meta[Type != '-9']

# read in pca data
pca <- fread("~/workspace/Pterocarya/calls/PSTE_UCD1470_HAP2/P.ste_pruned.eigenvec")
eigenval <- scan("~/workspace/Pterocarya/calls/PSTE_UCD1470_HAP2/P.ste_pruned.eigenval")

pve <- eigenval/sum(eigenval)*100

# reformat pca data
pca <- pca[, -c(1)]
setnames(pca, c("ID", paste0("PC", 1:(ncol(pca)-1))))

pca <- merge(pca, meta, by = 'ID')


pca[grepl("PSTE_SA", ID), Locality := 'Sacramento']
pca[grepl("PSTE_DV", ID), Locality := 'Davis']
pca[grepl("PSTE_UCD", ID), Locality := 'Davis']
pca[grepl("PSTE_W", ID), Locality := 'NCGR Wolfskill']
pca[grepl("PSTE_SO", ID), Locality := 'Sonoma Bot. Garden']

#pca[grepl("SRR", ID), Locality := 'SRA']

pca[grepl("PSTE_SA", ID), Locality_abbrv := '1']
pca[grepl("PSTE_DV", ID), Locality_abbrv := '2']
pca[grepl("PSTE_UCD", ID), Locality_abbrv := '2']
pca[grepl("PSTE_W", ID), Locality_abbrv := '3']
pca[grepl("PSTE_SO", ID), Locality_abbrv := '4']
#pca[grepl("SRR", ID), Locality := 'SRA']




# add variable for putative families
pca[, unique(ID)]

pca[grepl("PSTE_SO_1998", ID), Family := 'f1']
pca[grepl("PSTE_SO_2003", ID), Family := 'f2']

pca[ID %in% c("PSTE_W1_1", "PSTE_W1_5", "PSTE_WS_1_6"), Family := 'f3']
pca[ID %in% c("PSTE_WS_2_8", "PSTE_WS_2_10", "PSTE_WS_2_6"), Family := 'f4']
pca[ID %in% c("PSTE_WS_10_5", "PSTE_W10_1", "PSTE_W10_5"), Family := 'f5']
pca[ID %in% c("PSTE_W11_1", "PSTE_W11_2"), Family := 'f6']
pca[is.na(Family), Family := 'unknown']





ggplot(pca[Type != 'unknown'], aes(x = PC1, y = PC2)) + 
  geom_jitter(aes(color = Type, shape = Family), size = 1.5, stroke=0.8, width = 0.005) + 
  
  scale_shape_manual(values = c('unknown' = 16, 'f1' = 0, 
                                'f2' = 2, 'f3' = 1, 'f4' = 5,
                                'f5' = 6, 'f6' = 4)) +
  
  scale_color_manual(values = c('darkmagenta','lightsalmon'),
                     labels = c("protogynous", 'protandrous')) +
  #geom_text(aes(label = Locality_abbrv), nudge_y = 0.03) +
  
  #scale_color_manual(values = friendly_pal(name = 'zesty_four', n = 4)) +
  
  # customize labels
  xlab(paste0("PC1 (", signif(pve[1], 3), "%)")) + 
  ylab(paste0("PC2 (", signif(pve[2], 3), "%)")) + 
  theme_classic() + 
  theme(aspect.ratio = 1,
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        axis.title.x = element_text(margin=margin(t=10)),
        axis.title.y = element_text(margin=margin(r=10)),
        legend.key.size = unit(0.3, 'cm'),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 8)
  ) 


