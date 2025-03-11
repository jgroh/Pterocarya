library(ggplot2)
library(data.table)
library(seqinr)
library(ape)

##########################################
# Analysis of non-singleton biallelic snps 
# phased across the PsFAFL transcript.
##########################################


# ----- SNP positions ----- 

leg <- fread("~/workspace/Pterocarya/calls/PSTE_UCD1470_HAP2/PsFAFL.impute.legend")

# ID column in this table is both redundant and unecessary
leg[, ID := NULL]

# ----- haps ----- Each column is a biallelic SNP position. Columns are haplotypes.
haps0 <- fread("~/workspace/Pterocarya/calls/PSTE_UCD1470_HAP2/PsFAFL.impute.hap")
colnames(haps0)
nc <- ncol(haps0)
setnames(haps0, paste0("V", 1:nc), paste0(rep(seq(1:(nc/2)), each = 2), c("_1", "_2")))

d0 <- melt(cbind(leg, haps0), id.vars = c("pos", "allele0", "allele1"), value.name = 'allele')
d0[, c("indiv.id", "hap.id") := tstrsplit(variable, split = "_", fixed = T)]
d0[, variable := NULL]

# ----- individual IDs -----
indiv <- fread("~/workspace/Pterocarya/calls/PSTE_UCD1470_HAP2/PsFAFL.impute.hap.indv", header = F, col.names = c("ID"))
indiv[, indiv.id := as.character(seq(1, .N))]

d1 <- merge(d0, indiv, by = 'indiv.id', all = T)
d1[, indiv.id := NULL]


# ---- merge with genotype (assigned by phenotype or coverage) -----

geno <- fread("~/workspace/Pterocarya/P.stenoptera_MK_genotypes.txt", header = F, col.names = c("ID", 'genotype'))
haps <- merge(geno, d1, by = 'ID', all = T)

haps

# ----- investigate distribution of nonreference alleles among haplotypes
a <- haps[hap.id == 1, .(nonref_alleles_hap1 = sum(allele)), by = .(ID,genotype)]
b <- haps[hap.id == 2, .(nonref_alleles_hap2 = sum(allele)), by = .(ID,genotype)]
nonref_alleles <- merge(a,b)
nonref_allhaps <- melt(nonref_alleles, id.vars = c("ID", 'genotype'), value.name = "Non ref alleles")

ggplot(nonref_allhaps, aes(x = `Non ref alleles`)) + geom_histogram() + 
  theme_classic() + 
  theme(aspect.ratio = 1)

# specifically look at heterozygotes (expect bimodal)
nonref_alleles[genotype == 'Gg', hap_more := max(nonref_alleles_hap1, nonref_alleles_hap2), by = .(ID,genotype)]
nonref_alleles[genotype == 'Gg', hap_less := min(nonref_alleles_hap1, nonref_alleles_hap2), by = .(ID,genotype)]

nonref <- melt(nonref_alleles[, .(ID, hap_more, hap_less, genotype)], id.vars = c('ID', 'genotype'), value.name = 'nonref', variable.name = 'haplotype')
ggplot(nonref[genotype == 'Gg'], aes(x = nonref)) + geom_histogram()
ggplot(nonref[genotype == 'Gg'], aes(x = haplotype, y = nonref, group = ID)) + geom_line()

# ------ Assign tentative haplotype identities -----
# first we assign the G haplotypes of heterozygotes (reference is G haplotype, so hap with more alleles is g)
tmp <- haps[genotype == 'Gg', .(nonref = sum(allele)), by = .(ID, hap.id)]
tmp <- tmp[, .SD[which.min(nonref)], by = ID]
tmp[, nonref := NULL]
tmp[, haplotype := 'G']
haps <- merge(haps, tmp, all.x = T, by = c("ID", "hap.id"))

# now assign g haplotypes to homozygotes
haps[genotype == 'gg', haplotype := 'g']
haps[genotype == 'GG', haplotype := 'G']

# remaining unassigned ones are g haplotypes of heterozygotes
haps[is.na(haplotype), haplotype := 'g']



# ----- plot haplotypes -----
haps[, SNPrank := seq(1, .N), by = .(ID, hap.id)]

setkey(haps, haplotype, pos)

haps[, hapRank := seq(1, .N), by = .(pos)]

# look at all haplotypes. Phasing looks correct here.
ggplot(haps[], 
       aes(x = SNPrank, y = hapRank)) + 
  facet_wrap(~haplotype, scales = 'free') +

  geom_tile(aes(fill = as.factor(allele)), width = 1, height = 1) + 
  scale_fill_manual(values = c('cornsilk', 'chocolate4')) + 
  theme_classic() + 
  labs(x = '', y = '', fill = 'Allele') + 
  theme(axis.text = element_blank(),
        legend.position = 'none')


# look at heterozygotes by individual - 
# would be possible to spot phasing errors this way, but none seen.
ggplot(haps[genotype == 'Gg'], 
       aes(x = SNPrank, y = hap.id)) + 
  facet_wrap(~ID) + 
  geom_tile(aes(fill = as.factor(allele)), width = 1, height = 1) + 
  scale_fill_manual(values = c('cornsilk', 'chocolate4')) + 
  theme_classic() + 
  labs(x = '', y = '', fill = 'Allele') + 
  theme(axis.text = element_blank(),
        legend.position = 'none', aspect.ratio = .5)  
  



# ----- MK test -----

# sequence
PsFAFL_seq <- read.fasta("~/workspace/Pterocarya/01_Gloc/genes/FAF/CDS/Pste_dom.fasta")
PsFAFL_seq <- as.character(PsFAFL_seq[[1]])
length(PsFAFL_seq) #948 (including stop codon), 315 AA

# coding sequence coordinates
PsFAFL_cds_coords <- fread("~/workspace/Pterocarya/01_Gloc/genes/FAF/CDS/PsFAFL_CDS_coords_hap2.gff", header = F, sep = '\t')
PsFAFL_cds_coords <- PsFAFL_cds_coords[, c(4,5)]
setnames(PsFAFL_cds_coords, c("start", "end"))
setkey(PsFAFL_cds_coords, "start")

# combine positions and CDS sequence
# *** If adapting this code, double check result. positions may need to be reversed depending on strand
PsFAFL_cds <- data.table(PsFAFL_cds_coords[, seq(start,end), by = start][, .(pos = V1)], refseq=PsFAFL_seq)
PsFAFL_cds[, seq_rank := seq(1, .N)]

PsFAFL_cds[seq_rank %% 3 == 1, codon_position := 1]
PsFAFL_cds[seq_rank %% 3 == 2, codon_position := 2]
PsFAFL_cds[seq_rank %% 3 == 0, codon_position := 3]
# verify positions against annotation

# remove stop codon
PsFAFL_cds <- PsFAFL_cds[1:(.N-3)]


# get allele counts by haplotype at each site
PsFAFL_cnts <- haps[, .(count = sum(allele)), by = .(pos, haplotype, allele0, allele1)]
PsFAFL_cnts <- dcast(PsFAFL_cnts, pos + allele0 + allele1 ~ haplotype, value.var = 'count')


# how many of each haplotype?
haps[pos == min(pos), table(haplotype)]
# 56 g, 
# 22 G

# UTR
PsFAFL_UTR <- PsFAFL_cnts[pos < 4366881 | pos > 4367828]

# coding seq final data set
PsFAFL_final <- merge(PsFAFL_cnts[pos >= 4366881 & pos < 4367828],PsFAFL_cds, by = 'pos')
setkey(PsFAFL_final, seq_rank)

# function to classify synonymous vs nonsynonymous
syn_vs_nonsyn <- function(CDS_dt, rnk, nuc, codon_position){
  CDS <- copy(CDS_dt)
  if(codon_position == 1){
    refaa <- translate(CDS[seq_rank %in% rnk:(rnk+2), refseq])
    qryaa <- translate(c(nuc, CDS[seq_rank %in% (rnk+1):(rnk+2), refseq]))
  }
  if(codon_position == 2){
    refaa <- translate(CDS[seq_rank %in% (rnk-1):(rnk+1), refseq])
    qryaa <- translate( c( CDS[seq_rank %in% (rnk-1), refseq], nuc, CDS[seq_rank %in% (rnk+1), refseq]) )
  }
  if(codon_position == 3){
    refaa <- translate(CDS[seq_rank %in% (rnk-2):(rnk), refseq])
    qryaa <- translate( c( CDS[seq_rank %in% (rnk-2):(rnk-1), refseq], nuc ) )
  }
  if(refaa == qryaa){
    return('S')
  } else if(refaa != qryaa){
    return('N')
  }
}

syn_vs_nonsyn_v <- Vectorize(syn_vs_nonsyn)
vals <- vector()

for(i in 1:nrow(PsFAFL_final)){
  vals[i] <- syn_vs_nonsyn(CDS_dt = PsFAFL_cds, rnk = PsFAFL_final$seq_rank[i], nuc = PsFAFL_final$allele1[i], codon_position = PsFAFL_final$codon_position[i])
}

PsFAFL_final[, type := vals]
PsFAFL_final[,1:9]
# make sure you double check final numbers when assigning these classes
PsFAFL_final[(G >= 21 & g <= 1) | (G <= 1 & g >= 55), fixed_poly := 'fixed']
PsFAFL_final[g > 1 & g < 55 & G > 1 & G < 21, fixed_poly := 'polymorphic_both']

PsFAFL_final[g > 1  & g < 55 & (G >= 21 | G <= 1), fixed_poly := 'polymorphic_g']
PsFAFL_final[G < 21 & G > 1 & (g <= 1 | g >= 55), fixed_poly := 'polymorphic_G']
PsFAFL_final[, fixed_poly_MK := fixed_poly]
PsFAFL_final[fixed_poly == 'polymorphic_both', fixed_poly_MK := 'polymorphic_g']

PsFAFL_tbl <- table(PsFAFL_final$fixed_poly,PsFAFL_final$type)
PsFAFL_MK_tbl <- table(PsFAFL_final$fixed_poly_MK, PsFAFL_final$type)
PsFAFL_MK_tbl

fisher.test(PsFAFL_MK_tbl[1:2,])

PsFAFL_final[fixed_poly == 'fixed' & type == 'N']
PsFAFL_final[pos == 4367112]


fxd <-  ggplot(PsFAFL_final[fixed_poly == 'fixed'], aes(x = type)) + 
  geom_bar(stat = 'count')  + 
  theme_classic() + 
  #scale_y_continuous(limits = c(0, 20), breaks = 1:10) + 
  labs(title = 'Fixed', x = '', y = 'Count')  + 
  theme(aspect.ratio = 1)

polymorphic_G <- ggplot(PsFAFL_final[fixed_poly == 'polymorphic_G'], aes(x = type)) + 
  geom_bar(stat = 'count')  + 
  theme_classic() + 
  #scale_y_continuous(limits = c(0, 20), breaks = 1:10) + 
  labs(title = 'Polymorphic in G', x = '', y = '') + 
  theme(aspect.ratio = 1)

polymorphic_g <- ggplot(PsFAFL_final[fixed_poly == 'polymorphic_g'], aes(x = type)) + 
  geom_bar(stat = 'count')  + 
  theme_classic() + 
  labs(title = 'Polymorphic in g', x = '', y = '') + 
  #scale_y_continuous(limits = c(0, 20), breaks = 1:10) + 
  theme(aspect.ratio = 1)

library(gridExtra)
grid.arrange(fxd,polymorphic_G, polymorphic_g, ncol = 3)



# ----- FIL1 -----
# sequence
FIL1_seq <- read.fasta("~/workspace/heterodichogamy/HC_locus_structure/Pawnee_genes/stamen-specific-protein-FIL1-like/Pawnee_FIL1_CDS.fasta")
FIL1_seq <- as.character(FIL1_seq[[1]])
length(FIL1_seq) #300

# coding sequence coordinates
FIL1_cds_coords <- fread("~/workspace/heterodichogamy/HC_locus_structure/Pawnee_genes/stamen-specific-protein-FIL1-like/FIL1_cds.gff")
FIL1_cds_coords <- FIL1_cds_coords[, c(4,5)]
setnames(FIL1_cds_coords, c("start", "end"))
setkey(FIL1_cds_coords, "start")

# combine positions and CDS sequence
FIL1_cds <- data.table(FIL1_cds_coords[, seq(start,end), by = start][, .(pos = rev(V1))], refseq=FIL1_seq)
FIL1_cds[, seq_rank := seq(1, .N)]

FIL1_cds[seq_rank %% 3 == 1, codon_position := 1]
FIL1_cds[seq_rank %% 3 == 2, codon_position := 2]
FIL1_cds[seq_rank %% 3 == 0, codon_position := 3]
# verify positions against annotation

# are phasing errors evident for this gene? No. Assigned haplotypes should be correct
haps[pos > 6534375 & pos < 6534778]
ggplot(haps[genotype == 'Hh'], 
       aes(x = SNPrank, y = hap.id)) + 
  facet_wrap(~run) + 
  geom_tile(aes(fill = as.factor(allele)), width = 1, height = 1) + 
  scale_fill_manual(values = c('cornsilk', 'chocolate4')) + 
  theme_classic() + 
  labs(x = '', y = '', fill = 'Allele') + 
  theme(axis.text = element_blank(),
        legend.position = 'none', aspect.ratio = .1)  + 
  geom_vline(xintercept = 285, color = 'turquoise2') + 
  geom_vline(xintercept = 293, color = 'turquoise2') 


# get allele counts by haplotype at each site
FIL1_cnts <- haps[pos >= 6534375 & pos <= 6534778, .(count = sum(allele)), by = .(pos, haplotype, allele0, allele1)]
FIL1_cnts <- dcast(FIL1_cnts, pos + allele0 + allele1 ~ haplotype, value.var = 'count')

# **reverse complement these to be on the correct strand**
FIL1_cnts[, c('allele0', 'allele1') := lapply(.SD, comp), .SDcols = c('allele0', 'allele1')]

# how many of each haplotype?
haps[pos == min(pos), table(haplotype)]
# 47, 19

FIL1_final <- merge(FIL1_cnts, FIL1_cds)
setkey(FIL1_final, seq_rank)

# calculate 
vals <- vector()
for(i in 1:nrow(FIL1_final)){
  vals[i] <- syn_vs_nonsyn(CDS_dt = FIL1_cds, rnk = FIL1_final$seq_rank[i], nuc = FIL1_final$allele1[i], codon_position = FIL1_final$codon_position[i])
}

FIL1_final[, type := vals]

# make sure you double check final numbers when assigning these classes
FIL1_final[(H >= 18 & h <= 1) | (H <= 1 & h >= 46), fixed_poly := 'fixed']
FIL1_final[h > 1 & h < 46 & H > 1 & H < 18, fixed_poly := 'polymorphic_both']

FIL1_final[h > 1  & h < 46 & (H >= 18 | H <= 1), fixed_poly := 'polymorphic_h']
FIL1_final[H < 18 & H > 1 & (h <= 1 | h >= 46), fixed_poly := 'polymorphic_H']
FIL1_final[, fixed_poly_MK := fixed_poly]
FIL1_final[fixed_poly == 'polymorphic_both', fixed_poly_MK := 'polymorphic_h']

FIL1_tbl <- table(FIL1_final$fixed_poly, FIL1_final$type)
FIL1_MK_tbl <- table(FIL1_final$fixed_poly_MK, FIL1_final$type)
fisher.test(FIL1_MK_tbl)





# ----- EMS1 -----
# sequence
EMS1_seq <- read.fasta("~/workspace/heterodichogamy/HC_locus_structure/Pawnee_genes/EMS1-like/Pawnee_EMS1_cds.fasta")
EMS1_seq <- as.character(EMS1_seq[[1]])
length(EMS1_seq) #3885

# coding sequence coordinates
EMS1_cds_coords <- fread("~/workspace/heterodichogamy/HC_locus_structure/Pawnee_genes/EMS1-like/EMS1_cds.gff")
EMS1_cds_coords <- EMS1_cds_coords[, c(4,5)]
setnames(EMS1_cds_coords, c("start", "end"))
setkey(EMS1_cds_coords, "start")

# combine positions and CDS sequence
EMS1_cds <- data.table(EMS1_cds_coords[, seq(start,end), by = start][, .(pos = rev(V1))], refseq=EMS1_seq)
EMS1_cds[, seq_rank := seq(1, .N)]

EMS1_cds[seq_rank %% 3 == 1, codon_position := 1]
EMS1_cds[seq_rank %% 3 == 2, codon_position := 2]
EMS1_cds[seq_rank %% 3 == 0, codon_position := 3]
# verify positions against annotation

# are phasing errors evident for this gene? Yes, leave these individuals out to be safe
haps[pos > 6471551 & pos < 6475435]
ggplot(haps[genotype == 'Hh'], 
       aes(x = SNPrank, y = hap.id)) + 
  facet_wrap(~run) + 
  geom_tile(aes(fill = as.factor(allele)), width = 1, height = 1) + 
  scale_fill_manual(values = c('cornsilk', 'chocolate4')) + 
  theme_classic() + 
  labs(x = '', y = '', fill = 'Allele') + 
  theme(axis.text = element_blank(),
        legend.position = 'none', aspect.ratio = .1)  + 
  geom_vline(xintercept = 10, color = 'turquoise2') + 
  geom_vline(xintercept = 74, color = 'turquoise2') 

# get rid of individuals with ambiguous phasing
haps_EMS1 <- haps[!run %in% c("SRR15911530", "SRR15911531", "SRR15911546", "SRR15911528")]

# get allele counts by haplotype at each site
EMS1_cnts <- haps_EMS1[pos >= 6471551 & pos <= 6475435, .(count = sum(allele)), by = .(pos, haplotype, allele0, allele1)]
EMS1_cnts <- dcast(EMS1_cnts, pos + allele0 + allele1 ~ haplotype, value.var = 'count')

# **reverse complement these to be on the correct strand**
EMS1_cnts[, c('allele0', 'allele1') := lapply(.SD, comp), .SDcols = c('allele0', 'allele1')]

# how many of each haplotype?
haps_EMS1[pos == min(pos), table(haplotype)]
# 43, 15

EMS1_final <- merge(EMS1_cnts, EMS1_cds)
setkey(EMS1_final, seq_rank)

# calculate 
vals <- vector()
for(i in 1:nrow(EMS1_final)){
  vals[i] <- syn_vs_nonsyn(CDS_dt = EMS1_cds, rnk = EMS1_final$seq_rank[i], nuc = EMS1_final$allele1[i], codon_position = EMS1_final$codon_position[i])
}

EMS1_final[, type := vals]

# make sure you double check final numbers when assigning these classes
EMS1_final[(H >= 14 & h <= 1) | (H <= 1 & h >= 42), fixed_poly := 'fixed']
EMS1_final[h > 1 & h < 42 & H > 1 & H < 14, fixed_poly := 'polymorphic_both']

EMS1_final[h > 1  & h < 42 & (H >= 14 | H <= 1), fixed_poly := 'polymorphic_h']
EMS1_final[H < 14 & H > 1 & (h <= 1 | h >= 42), fixed_poly := 'polymorphic_H']
EMS1_final[, fixed_poly_MK := fixed_poly]
EMS1_final[fixed_poly == 'polymorphic_both', fixed_poly_MK := 'polymorphic_h']

EMS1_tbl <- table(EMS1_final$fixed_poly, EMS1_final$type)
EMS1_MK_tbl <- table(EMS1_final$fixed_poly_MK, EMS1_final$type)
fisher.test(EMS1_MK_tbl)

fxd <-  ggplot(EMS1_final[fixed_poly == 'fixed'], aes(x = type)) + 
  geom_bar(stat = 'count')  + 
  theme_classic() + 
  #scale_y_continuous(limits = c(0, 20), breaks = 1:10) + 
  labs(title = 'Fixed', x = '', y = 'Count')  + 
  theme(aspect.ratio = 1)

polymorphic_H <- ggplot(EMS1_final[fixed_poly == 'polymorphic_H'], aes(x = type)) + 
  geom_bar(stat = 'count')  + 
  theme_classic() + 
  #scale_y_continuous(limits = c(0, 20), breaks = 1:10) + 
  labs(title = 'Polymorphic in H', x = '', y = '') + 
  theme(aspect.ratio = 1)

polymorphic_h <- ggplot(EMS1_final[fixed_poly == 'polymorphic_h'], aes(x = type)) + 
  geom_bar(stat = 'count')  + 
  theme_classic() + 
  labs(title = 'Polymorphic in h', x = '', y = '') + 
  #scale_y_continuous(limits = c(0, 20), breaks = 1:10) + 
  theme(aspect.ratio = 1)

grid.arrange(fxd,polymorphic_H, polymorphic_h, ncol = 3)





# ----- RHOMBOID -----
# sequence
RHOM_seq <- read.fasta("~/workspace/heterodichogamy/HC_locus_structure/Pawnee_genes/RHOMBOID-like-protein2/Pawnee_RHOMBOID_cds.fasta")
RHOM_seq <- as.character(RHOM_seq[[1]])
length(RHOM_seq) #990

# coding sequence coordinates
RHOM_cds_coords <- fread("~/workspace/heterodichogamy/HC_locus_structure/Pawnee_genes/RHOMBOID-like-protein2/RHOMBOID_cds.gff")
RHOM_cds_coords <- RHOM_cds_coords[, c(4,5)]
setnames(RHOM_cds_coords, c("start", "end"))
setkey(RHOM_cds_coords, "start")

# combine positions and CDS sequence
RHOM_cds <- data.table(RHOM_cds_coords[, seq(start,end), by = start][, .(pos = V1)], refseq=RHOM_seq)
RHOM_cds[, seq_rank := seq(1, .N)]

RHOM_cds[seq_rank %% 3 == 1, codon_position := 1]
RHOM_cds[seq_rank %% 3 == 2, codon_position := 2]
RHOM_cds[seq_rank %% 3 == 0, codon_position := 3]
# verify positions against annotation

# are phasing errors evident for this gene? No. Assigned haplotypes should be correct here
haps[pos >= 6514141 & pos <= 6517293]
ggplot(haps[genotype == 'Hh'], 
       aes(x = SNPrank, y = hap.id)) + 
  facet_wrap(~run) + 
  geom_tile(aes(fill = as.factor(allele)), width = 1, height = 1) + 
  scale_fill_manual(values = c('cornsilk', 'chocolate4')) + 
  theme_classic() + 
  labs(x = '', y = '', fill = 'Allele') + 
  theme(axis.text = element_blank(),
        legend.position = 'none', aspect.ratio = .1)  + 
  geom_vline(xintercept = 176, color = 'turquoise2') + 
  geom_vline(xintercept = 195, color = 'turquoise2') 

# get allele counts by haplotype at each site
RHOM_cnts <- haps[pos >= 6514141 & pos <= 6517293, .(count = sum(allele)), by = .(pos, haplotype, allele0, allele1)]
RHOM_cnts <- dcast(RHOM_cnts, pos + allele0 + allele1 ~ haplotype, value.var = 'count')


# how many of each haplotype?
haps[pos == min(pos), table(haplotype)]
# 47, 19

RHOM_final <- merge(RHOM_cnts, RHOM_cds, by = 'pos')
setkey(RHOM_final, seq_rank)

# calculate 
vals <- vector()
for(i in 1:nrow(RHOM_final)){
  vals[i] <- syn_vs_nonsyn(CDS_dt = RHOM_cds, rnk = RHOM_final$seq_rank[i], nuc = RHOM_final$allele1[i], codon_position = RHOM_final$codon_position[i])
}

RHOM_final[, type := vals]

# make sure you double check final numbers when assigning these classes
RHOM_final[(H >= 18 & h <= 1) | (H <= 1 & h >= 46), fixed_poly := 'fixed']
RHOM_final[h > 1 & h < 46 & H > 1 & H < 18, fixed_poly := 'polymorphic_both']

RHOM_final[h > 1  & h < 46 & (H >= 18 | H <= 1), fixed_poly := 'polymorphic_h']
RHOM_final[H < 18 & H > 1 & (h <= 1 | h >= 46), fixed_poly := 'polymorphic_H']
RHOM_final[, fixed_poly_MK := fixed_poly]
RHOM_final[fixed_poly == 'polymorphic_both', fixed_poly_MK := 'polymorphic_h']

RHOM_tbl <- table(RHOM_final$fixed_poly, RHOM_final$type)
RHOM_MK_tbl <- table(RHOM_final$fixed_poly_MK, RHOM_final$type)
fisher.test(RHOM_MK_tbl)



# ----- CEN1 -----

# sequence
CEN1_seq <- read.fasta("~/workspace/heterodichogamy/HC_locus_structure/Pawnee_genes/CEN-like-protein1/Pawnee_CEN1_cds.fasta")
CEN1_seq <- as.character(CEN1_seq[[1]])
length(CEN1_seq) #510

# coding sequence coordinates
CEN1_cds_coords <- fread("~/workspace/heterodichogamy/HC_locus_structure/Pawnee_genes/CEN-like-protein1/CEN1_cds.gff")
CEN1_cds_coords <- CEN1_cds_coords[, c(4,5)]
setnames(CEN1_cds_coords, c("start", "end"))
setkey(CEN1_cds_coords, "start")
CEN1_cds_coords[start == 6497654, start := 6497666]

# combine positions and CDS sequence
CEN1_cds <- data.table(CEN1_cds_coords[, seq(start,end), by = start][, .(pos = V1)], refseq=CEN1_seq)
CEN1_cds[, seq_rank := seq(1, .N)]

CEN1_cds[seq_rank %% 3 == 1, codon_position := 1]
CEN1_cds[seq_rank %% 3 == 2, codon_position := 2]
CEN1_cds[seq_rank %% 3 == 0, codon_position := 3]
# verify positions against annotation

# are phasing errors evident for this gene? No. Assigned haplotypes should be correct here
haps[pos >= 6497666 & pos <= 6498545]
ggplot(haps[genotype == 'Hh'], 
       aes(x = SNPrank, y = hap.id)) + 
  facet_wrap(~run) + 
  geom_tile(aes(fill = as.factor(allele)), width = 1, height = 1) + 
  scale_fill_manual(values = c('cornsilk', 'chocolate4')) + 
  theme_classic() + 
  labs(x = '', y = '', fill = 'Allele') + 
  theme(axis.text = element_blank(),
        legend.position = 'none', aspect.ratio = .1)  + 
  geom_vline(xintercept = 75, color = 'turquoise2') + 
  geom_vline(xintercept = 90, color = 'turquoise2') 

# get rid of individuals with ambiguous phasing
haps_CEN1 <- haps[!run %in% c("SRR15911530", "SRR15911531", "SRR15911546", "SRR15911528")]


# get allele counts by haplotype at each site
CEN1_cnts <- haps_CEN1[pos >= 6497666 & pos <= 6498545, .(count = sum(allele)), by = .(pos, haplotype, allele0, allele1)]
CEN1_cnts <- dcast(CEN1_cnts, pos + allele0 + allele1 ~ haplotype, value.var = 'count')


# how many of each haplotype?
haps_CEN1[pos == min(pos), table(haplotype)]
# 43, 15

CEN1_final <- merge(CEN1_cnts, CEN1_cds, by = 'pos')
setkey(CEN1_final, seq_rank)

# calculate 
vals <- vector()
for(i in 1:nrow(CEN1_final)){
  vals[i] <- syn_vs_nonsyn(CDS_dt = CEN1_cds, rnk = CEN1_final$seq_rank[i], nuc = CEN1_final$allele1[i], codon_position = CEN1_final$codon_position[i])
}

CEN1_final[, type := vals]

# make sure you double check final numbers when assigning these classes
CEN1_final[(H >= 14 & h <= 1) | (H <= 1 & h >= 42), fixed_poly := 'fixed']
CEN1_final[h > 1 & h < 42 & H > 1 & H < 14, fixed_poly := 'polymorphic_both']

CEN1_final[h > 1  & h < 42 & (H >= 14 | H <= 1), fixed_poly := 'polymorphic_h']
CEN1_final[H < 14 & H > 1 & (h <= 1 | h >= 42), fixed_poly := 'polymorphic_H']
CEN1_final[, fixed_poly_MK := fixed_poly]
CEN1_final[fixed_poly == 'polymorphic_both', fixed_poly_MK := 'polymorphic_h']

CEN1_tbl <- table(CEN1_final$fixed_poly, CEN1_final$type)
CEN1_MK_tbl <- table(CEN1_final$fixed_poly_MK, CEN1_final$type)
fisher.test(CEN1_MK_tbl)



# ----- XBAT33 -----

# sequence
XBAT33_seq <- read.fasta("~/workspace/heterodichogamy/HC_locus_structure/Pawnee_genes/XBAT33-like/Pawnee_XBAT33_cds.fasta")
XBAT33_seq <- as.character(XBAT33_seq[[1]])
length(XBAT33_seq) #1539

# coding sequence coordinates
XBAT33_cds_coords <- fread("~/workspace/heterodichogamy/HC_locus_structure/Pawnee_genes/XBAT33-like/Pawnee_XBAT33_cds.gff")
XBAT33_cds_coords <- XBAT33_cds_coords[, c(4,5)]
setnames(XBAT33_cds_coords, c("start", "end"))
setkey(XBAT33_cds_coords, "start")

# combine positions and CDS sequence
XBAT33_cds <- data.table(XBAT33_cds_coords[, seq(start,end), by = start][, .(pos = rev(V1))], refseq=XBAT33_seq)
XBAT33_cds[, seq_rank := seq(1, .N)]

XBAT33_cds[seq_rank %% 3 == 1, codon_position := 1]
XBAT33_cds[seq_rank %% 3 == 2, codon_position := 2]
XBAT33_cds[seq_rank %% 3 == 0, codon_position := 3]
# verify positions against annotation

# are phasing errors evident for this gene? No. Assigned haplotypes should be correct here
haps[pos >= 6501075 & pos <= 6505794]
ggplot(haps[genotype == 'Hh'], 
       aes(x = SNPrank, y = hap.id)) + 
  facet_wrap(~run) + 
  geom_tile(aes(fill = as.factor(allele)), width = 1, height = 1) + 
  scale_fill_manual(values = c('cornsilk', 'chocolate4')) + 
  theme_classic() + 
  labs(x = '', y = '', fill = 'Allele') + 
  theme(axis.text = element_blank(),
        legend.position = 'none', aspect.ratio = .1)  + 
  geom_vline(xintercept = 91, color = 'turquoise2') + 
  geom_vline(xintercept = 133, color = 'turquoise2') 

# get rid of individuals with ambiguous phasing
haps_XBAT33 <- haps[!run %in% c("SRR15911530", "SRR15911531", "SRR15911546", "SRR15911528")]


# get allele counts by haplotype at each site
XBAT33_cnts <- haps_XBAT33[pos >= 6501075 & pos <= 6505794, .(count = sum(allele)), by = .(pos, haplotype, allele0, allele1)]
XBAT33_cnts <- dcast(XBAT33_cnts, pos + allele0 + allele1 ~ haplotype, value.var = 'count')

# **reverse complement these to be on the correct strand**
XBAT33_cnts[, c('allele0', 'allele1') := lapply(.SD, comp), .SDcols = c('allele0', 'allele1')]

# how many of each haplotype?
haps_XBAT33[pos == min(pos), table(haplotype)]
# 43, 15

XBAT33_final <- merge(XBAT33_cnts, XBAT33_cds, by = 'pos')
setkey(XBAT33_final, seq_rank)

# calculate 
vals <- vector()
for(i in 1:nrow(XBAT33_final)){
  vals[i] <- syn_vs_nonsyn(CDS_dt = XBAT33_cds, rnk = XBAT33_final$seq_rank[i], nuc = XBAT33_final$allele1[i], codon_position = XBAT33_final$codon_position[i])
}

XBAT33_final[, type := vals]

# make sure you double check final numbers when assigning these classes
XBAT33_final[(H >= 14 & h <= 1) | (H <= 1 & h >= 42), fixed_poly := 'fixed']
XBAT33_final[h > 1 & h < 42 & H > 1 & H < 14, fixed_poly := 'polymorphic_both']

XBAT33_final[h > 1  & h < 42 & (H >= 14 | H <= 1), fixed_poly := 'polymorphic_h']
XBAT33_final[H < 14 & H > 1 & (h <= 1 | h >= 42), fixed_poly := 'polymorphic_H']
XBAT33_final[, fixed_poly_MK := fixed_poly]
XBAT33_final[fixed_poly == 'polymorphic_both', fixed_poly_MK := 'polymorphic_h']

XBAT33_tbl <- table(XBAT33_final$fixed_poly, XBAT33_final$type)
XBAT33_MK_tbl <- table(XBAT33_final$fixed_poly_MK, XBAT33_final$type)
fisher.test(XBAT33_MK_tbl)


# ----- BAG1 -----

# sequence
BAG1_seq <- read.fasta("~/workspace/heterodichogamy/HC_locus_structure/Pawnee_genes/BAG-family-molecular-chaperone-regulator-1-like/Pawnee_BAG1_cds.fasta")
BAG1_seq <- as.character(BAG1_seq[[1]])
length(BAG1_seq) #1050

# coding sequence coordinates
BAG1_cds_coords <- fread("~/workspace/heterodichogamy/HC_locus_structure/Pawnee_genes/BAG-family-molecular-chaperone-regulator-1-like/Pawnee_BAG1_cds.gff")
BAG1_cds_coords <- BAG1_cds_coords[, c(4,5)]
setnames(BAG1_cds_coords, c("start", "end"))
setkey(BAG1_cds_coords, "start")

# combine positions and CDS sequence
BAG1_cds <- data.table(BAG1_cds_coords[, seq(start,end), by = start][, .(pos = V1)], refseq=BAG1_seq)
BAG1_cds[, seq_rank := seq(1, .N)]

BAG1_cds[seq_rank %% 3 == 1, codon_position := 1]
BAG1_cds[seq_rank %% 3 == 2, codon_position := 2]
BAG1_cds[seq_rank %% 3 == 0, codon_position := 3]
# verify positions against annotation

# are phasing errors evident for this gene? No. Assigned haplotypes should be correct here
haps[pos >= 6631726 & pos <= 6633223]
ggplot(haps[genotype == 'Hh'], 
       aes(x = SNPrank, y = hap.id)) + 
  facet_wrap(~run) + 
  geom_tile(aes(fill = as.factor(allele)), width = 1, height = 1) + 
  scale_fill_manual(values = c('cornsilk', 'chocolate4')) + 
  theme_classic() + 
  labs(x = '', y = '', fill = 'Allele') + 
  theme(axis.text = element_blank(),
        legend.position = 'none', aspect.ratio = .1)  + 
  geom_vline(xintercept = 563, color = 'turquoise2') + 
  geom_vline(xintercept = 581, color = 'turquoise2') 



# get allele counts by haplotype at each site
BAG1_cnts <- haps[pos >= 6631726 & pos <= 6633223, .(count = sum(allele)), by = .(pos, haplotype, allele0, allele1)]
BAG1_cnts <- dcast(BAG1_cnts, pos + allele0 + allele1 ~ haplotype, value.var = 'count')

# how many of each haplotype?
haps[pos == min(pos), table(haplotype)]
# 47, 19

BAG1_final <- merge(BAG1_cnts, BAG1_cds, by = 'pos')
setkey(BAG1_final, seq_rank)

# calculate 
vals <- vector()
for(i in 1:nrow(BAG1_final)){
  vals[i] <- syn_vs_nonsyn(CDS_dt = BAG1_cds, rnk = BAG1_final$seq_rank[i], nuc = BAG1_final$allele1[i], codon_position = BAG1_final$codon_position[i])
}

BAG1_final[, type := vals]

# make sure you double check final numbers when assigning these classes
BAG1_final[(H >= 18 & h <= 1) | (H <= 1 & h >= 46), fixed_poly := 'fixed']
BAG1_final[h > 1 & h < 46 & H > 1 & H < 18, fixed_poly := 'polymorphic_both']

BAG1_final[h > 1  & h < 46 & (H >= 18 | H <= 1), fixed_poly := 'polymorphic_h']
BAG1_final[H < 18 & H > 1 & (h <= 1 | h >= 46), fixed_poly := 'polymorphic_H']
BAG1_final[, fixed_poly_MK := fixed_poly]
BAG1_final[fixed_poly == 'polymorphic_both', fixed_poly_MK := 'polymorphic_h']

BAG1_tbl <- table(BAG1_final$fixed_poly, BAG1_final$type)
BAG1_MK_tbl <- table(BAG1_final$fixed_poly_MK, BAG1_final$type)
fisher.test(BAG1_MK_tbl)


# ----- ROOT PRIMORDIUM DEFECTIVE -----

# sequence
RPD_seq <- read.fasta("~/workspace/heterodichogamy/HC_locus_structure/Pawnee_genes/ROOT-PRIMORDIUM-DEFECTIVE-1/Pawnee_ROOT_PRIMORDIUM_DEFECTIVE_cds.fasta")
RPD_seq <- as.character(RPD_seq[[1]])
length(RPD_seq) #1401

# coding sequence coordinates
RPD_cds_coords <- fread("~/workspace/heterodichogamy/HC_locus_structure/Pawnee_genes/ROOT-PRIMORDIUM-DEFECTIVE-1/Pawnee_RPD_cds.gff")
RPD_cds_coords <- RPD_cds_coords[, c(4,5)]
setnames(RPD_cds_coords, c("start", "end"))
setkey(RPD_cds_coords, "start")

# combine positions and CDS sequence
RPD_cds <- data.table(RPD_cds_coords[, seq(start,end), by = start][, .(pos = rev(V1))], refseq=RPD_seq)
RPD_cds[, seq_rank := seq(1, .N)]

RPD_cds[seq_rank %% 3 == 1, codon_position := 1]
RPD_cds[seq_rank %% 3 == 2, codon_position := 2]
RPD_cds[seq_rank %% 3 == 0, codon_position := 3]
# verify positions against annotation

# are phasing errors evident for this gene? No. Assigned haplotypes should be correct here
haps[pos >= 6508200 & pos <= 6509600]
ggplot(haps[genotype == 'Hh'], 
       aes(x = SNPrank, y = hap.id)) + 
  facet_wrap(~run) + 
  geom_tile(aes(fill = as.factor(allele)), width = 1, height = 1) + 
  scale_fill_manual(values = c('cornsilk', 'chocolate4')) + 
  theme_classic() + 
  labs(x = '', y = '', fill = 'Allele') + 
  theme(axis.text = element_blank(),
        legend.position = 'none', aspect.ratio = .1)  + 
  geom_vline(xintercept = 134, color = 'turquoise2') + 
  geom_vline(xintercept = 175, color = 'turquoise2') 

# get rid of individuals with ambiguous phasing
haps_RPD <- haps[!run %in% c("SRR15911546")]


# get allele counts by haplotype at each site
RPD_cnts <- haps_RPD[pos >= 6508200 & pos <= 6509600, .(count = sum(allele)), by = .(pos, haplotype, allele0, allele1)]
RPD_cnts <- dcast(RPD_cnts, pos + allele0 + allele1 ~ haplotype, value.var = 'count')


# **reverse complement these to be on the correct strand**
RPD_cnts[, c('allele0', 'allele1') := lapply(.SD, comp), .SDcols = c('allele0', 'allele1')]


# how many of each haplotype?
haps_RPD[pos == min(pos), table(haplotype)]
# 46, 18

RPD_final <- merge(RPD_cnts, RPD_cds, by = 'pos')
setkey(RPD_final, seq_rank)

# calculate 
vals <- vector()
for(i in 1:nrow(RPD_final)){
  vals[i] <- syn_vs_nonsyn(CDS_dt = RPD_cds, rnk = RPD_final$seq_rank[i], nuc = RPD_final$allele1[i], codon_position = RPD_final$codon_position[i])
}

RPD_final[, type := vals]

# make sure you double check final numbers when assigning these classes
RPD_final[(H >= 17 & h <= 1) | (H <= 1 & h >= 45), fixed_poly := 'fixed']
RPD_final[h > 1 & h < 45 & H > 1 & H < 17, fixed_poly := 'polymorphic_both']

RPD_final[h > 1  & h < 45 & (H >= 17 | H <= 1), fixed_poly := 'polymorphic_h']
RPD_final[H < 17 & H > 1 & (h <= 1 | h >= 45), fixed_poly := 'polymorphic_H']
RPD_final[, fixed_poly_MK := fixed_poly]
RPD_final[fixed_poly == 'polymorphic_both', fixed_poly_MK := 'polymorphic_h']

RPD_tbl <- table(RPD_final$fixed_poly, RPD_final$type)
RPD_MK_tbl <- table(RPD_final$fixed_poly_MK, RPD_final$type)
fisher.test(RPD_MK_tbl)


# ----- PMT16 -----

# sequence
PMT_seq <- read.fasta("~/workspace/heterodichogamy/HC_locus_structure/Pawnee_genes/PMT16/Pawnee_PMT16_cds.fasta")
PMT_seq <- as.character(PMT_seq[[1]])
length(PMT_seq) #1752

# coding sequence coordinates
PMT_cds_coords <- fread("~/workspace/heterodichogamy/HC_locus_structure/Pawnee_genes/PMT16/Pawnee_PMT16_cds.gff")
PMT_cds_coords <- PMT_cds_coords[, c(4,5)]
setnames(PMT_cds_coords, c("start", "end"))
setkey(PMT_cds_coords, "start")

# combine positions and CDS sequence
PMT_cds <- data.table(PMT_cds_coords[, seq(start,end), by = start][, .(pos = rev(V1))], refseq=PMT_seq)
PMT_cds[, seq_rank := seq(1, .N)]

PMT_cds[seq_rank %% 3 == 1, codon_position := 1]
PMT_cds[seq_rank %% 3 == 2, codon_position := 2]
PMT_cds[seq_rank %% 3 == 0, codon_position := 3]
# verify positions against annotation

# are phasing errors evident for this gene? No. Assigned haplotypes should be correct here
haps[pos >= 6541912 & pos <= 6545860]
ggplot(haps[genotype == 'Hh'], 
       aes(x = SNPrank, y = hap.id)) + 
  facet_wrap(~run) + 
  geom_tile(aes(fill = as.factor(allele)), width = 1, height = 1) + 
  scale_fill_manual(values = c('cornsilk', 'chocolate4')) + 
  theme_classic() + 
  labs(x = '', y = '', fill = 'Allele') + 
  theme(axis.text = element_blank(),
        legend.position = 'none', aspect.ratio = .1)  + 
  geom_vline(xintercept = 294, color = 'turquoise2') + 
  geom_vline(xintercept = 366, color = 'turquoise2') 


# get allele counts by haplotype at each site
PMT_cnts <- haps[pos >= 6541912 & pos <= 6545860, .(count = sum(allele)), by = .(pos, haplotype, allele0, allele1)]
PMT_cnts <- dcast(PMT_cnts, pos + allele0 + allele1 ~ haplotype, value.var = 'count')

# **reverse complement these to be on the correct strand**
PMT_cnts[, c('allele0', 'allele1') := lapply(.SD, comp), .SDcols = c('allele0', 'allele1')]


# how many of each haplotype?
haps[pos == min(pos), table(haplotype)]
# 47, 19

PMT_final <- merge(PMT_cnts, PMT_cds, by = 'pos')
setkey(PMT_final, seq_rank)

# calculate 
vals <- vector()
for(i in 1:nrow(PMT_final)){
  vals[i] <- syn_vs_nonsyn(CDS_dt = PMT_cds, rnk = PMT_final$seq_rank[i], nuc = PMT_final$allele1[i], codon_position = PMT_final$codon_position[i])
}

PMT_final[, type := vals]

# make sure you double check final numbers when assigning these classes
PMT_final[(H >= 18 & h <= 1) | (H <= 1 & h >= 46), fixed_poly := 'fixed']
PMT_final[h > 1 & h < 46 & H > 1 & H < 18, fixed_poly := 'polymorphic_both']

PMT_final[h > 1  & h < 46 & (H >= 18 | H <= 1), fixed_poly := 'polymorphic_h']
PMT_final[H < 18 & H > 1 & (h <= 1 | h >= 46), fixed_poly := 'polymorphic_H']
PMT_final[, fixed_poly_MK := fixed_poly]
PMT_final[fixed_poly == 'polymorphic_both', fixed_poly_MK := 'polymorphic_h']

PMT_tbl <- table(PMT_final$fixed_poly, PMT_final$type)
PMT_MK_tbl <- table(PMT_final$fixed_poly_MK, PMT_final$type)
fisher.test(PMT_MK_tbl)





# ----- F-box PP2-A12 -----

# sequence
PP2_seq <- read.fasta("~/workspace/heterodichogamy/HC_locus_structure/Pawnee_genes/F-box-protein-PP2-A12-like/Pawnee_PP2-A12.fasta")
PP2_seq <- as.character(PP2_seq[[1]])
length(PP2_seq) #897

# coding sequence coordinates
PP2_cds_coords <- fread("~/workspace/heterodichogamy/HC_locus_structure/Pawnee_genes/F-box-protein-PP2-A12-like/Pawnee_PP2-A12_cds.gff")
PP2_cds_coords <- PP2_cds_coords[, c(4,5)]
setnames(PP2_cds_coords, c("start", "end"))
setkey(PP2_cds_coords, "start")

# combine positions and CDS sequence
PP2_cds <- data.table(PP2_cds_coords[, seq(start,end), by = start][, .(pos = V1)], refseq=PP2_seq)
PP2_cds[, seq_rank := seq(1, .N)]

PP2_cds[seq_rank %% 3 == 1, codon_position := 1]
PP2_cds[seq_rank %% 3 == 2, codon_position := 2]
PP2_cds[seq_rank %% 3 == 0, codon_position := 3]
# verify positions against annotation

# are phasing errors evident for this gene? No. Assigned haplotypes should be correct here
haps[pos >= 6571273 & pos <= 6573700]
ggplot(haps[genotype == 'Hh'], 
       aes(x = SNPrank, y = hap.id)) + 
  facet_wrap(~run) + 
  geom_tile(aes(fill = as.factor(allele)), width = 1, height = 1) + 
  scale_fill_manual(values = c('cornsilk', 'chocolate4')) + 
  theme_classic() + 
  labs(x = '', y = '', fill = 'Allele') + 
  theme(axis.text = element_blank(),
        legend.position = 'none', aspect.ratio = .1)  + 
  geom_vline(xintercept = 375, color = 'turquoise2') + 
  geom_vline(xintercept = 389, color = 'turquoise2') 


# get allele counts by haplotype at each site
PP2_cnts <- haps[pos >=6571273 & pos <= 6573700, .(count = sum(allele)), by = .(pos, haplotype, allele0, allele1)]
PP2_cnts <- dcast(PP2_cnts, pos + allele0 + allele1 ~ haplotype, value.var = 'count')


# how many of each haplotype?
haps[pos == min(pos), table(haplotype)]
# 47, 19

PP2_final <- merge(PP2_cnts, PP2_cds, by = 'pos')
setkey(PP2_final, seq_rank)

# calculate 
vals <- vector()
for(i in 1:nrow(PP2_final)){
  vals[i] <- syn_vs_nonsyn(CDS_dt = PP2_cds, rnk = PP2_final$seq_rank[i], nuc = PP2_final$allele1[i], codon_position = PP2_final$codon_position[i])
}

PP2_final[, type := vals]

# make sure you double check final numbers when assigning these classes
PP2_final[(H >= 18 & h <= 1) | (H <= 1 & h >= 46), fixed_poly := 'fixed']
PP2_final[h > 1 & h < 46 & H > 1 & H < 18, fixed_poly := 'polymorphic_both']

PP2_final[h > 1  & h < 46 & (H >= 18 | H <= 1), fixed_poly := 'polymorphic_h']
PP2_final[H < 18 & H > 1 & (h <= 1 | h >= 46), fixed_poly := 'polymorphic_H']
PP2_final[, fixed_poly_MK := fixed_poly]
PP2_final[fixed_poly == 'polymorphic_both', fixed_poly_MK := 'polymorphic_h']

PP2_tbl <- table(PP2_final$fixed_poly, PP2_final$type)
PP2_MK_tbl <- table(PP2_final$fixed_poly_MK, PP2_final$type)
fisher.test(PP2_MK_tbl)



# ----- F-box At4g00755 -----

# sequence
at4g_seq <- read.fasta("~/workspace/heterodichogamy/HC_locus_structure/Pawnee_genes/F-box-protein-At4g00755-like/Pawnee_At4g_cds.fasta")
at4g_seq <- as.character(at4g_seq[[1]])
length(at4g_seq) #1080

# coding sequence coordinates
at4g_cds_coords <- fread("~/workspace/heterodichogamy/HC_locus_structure/Pawnee_genes/F-box-protein-At4g00755-like/Pawnee_At4g_cds.gff")
at4g_cds_coords <- at4g_cds_coords[, c(4,5)]
setnames(at4g_cds_coords, c("start", "end"))
setkey(at4g_cds_coords, "start")

# combine positions and CDS sequence
at4g_cds <- data.table(at4g_cds_coords[, seq(start,end), by = start][, .(pos = V1)], refseq=at4g_seq)
at4g_cds[, seq_rank := seq(1, .N)]

at4g_cds[seq_rank %% 3 == 1, codon_position := 1]
at4g_cds[seq_rank %% 3 == 2, codon_position := 2]
at4g_cds[seq_rank %% 3 == 0, codon_position := 3]
# verify positions against annotation

# are phasing errors evident for this gene? No. Assigned haplotypes should be correct here
haps[pos >= 6583289 & pos <= 6591527]
ggplot(haps[genotype == 'Hh'], 
       aes(x = SNPrank, y = hap.id)) + 
  facet_wrap(~run) + 
  geom_tile(aes(fill = as.factor(allele)), width = 1, height = 1) + 
  scale_fill_manual(values = c('cornsilk', 'chocolate4')) + 
  theme_classic() + 
  labs(x = '', y = '', fill = 'Allele') + 
  theme(axis.text = element_blank(),
        legend.position = 'none', aspect.ratio = .1)  + 
  geom_vline(xintercept = 390, color = 'turquoise2') + 
  geom_vline(xintercept = 434, color = 'turquoise2') 


# get allele counts by haplotype at each site
at4g_cnts <- haps[pos >=6583289 & pos <= 6591527, .(count = sum(allele)), by = .(pos, haplotype, allele0, allele1)]
at4g_cnts <- dcast(at4g_cnts, pos + allele0 + allele1 ~ haplotype, value.var = 'count')


# how many of each haplotype?
haps[pos == min(pos), table(haplotype)]
# 47, 19

at4g_final <- merge(at4g_cnts, at4g_cds, by = 'pos')
setkey(at4g_final, seq_rank)

# calculate 
vals <- vector()
for(i in 1:nrow(at4g_final)){
  vals[i] <- syn_vs_nonsyn(CDS_dt = at4g_cds, rnk = at4g_final$seq_rank[i], nuc = at4g_final$allele1[i], codon_position = at4g_final$codon_position[i])
}

at4g_final[, type := vals]

# make sure you double check final numbers when assigning these classes
at4g_final[(H >= 18 & h <= 1) | (H <= 1 & h >= 46), fixed_poly := 'fixed']
at4g_final[h > 1 & h < 46 & H > 1 & H < 18, fixed_poly := 'polymorphic_both']

at4g_final[h > 1  & h < 46 & (H >= 18 | H <= 1), fixed_poly := 'polymorphic_h']
at4g_final[H < 18 & H > 1 & (h <= 1 | h >= 46), fixed_poly := 'polymorphic_H']
at4g_final[, fixed_poly_MK := fixed_poly]
at4g_final[fixed_poly == 'polymorphic_both', fixed_poly_MK := 'polymorphic_h']

at4g_tbl <- table(at4g_final$fixed_poly, at4g_final$type)
at4g_MK_tbl <- table(at4g_final$fixed_poly_MK, at4g_final$type)
fisher.test(at4g_MK_tbl)


# ----- IQ-Domain -----

# sequence
IQ_seq <- read.fasta("~/workspace/heterodichogamy/HC_locus_structure/Pawnee_genes/IQ-DOMAIN-23-like/Pawnee_IQ_cds.fasta")
IQ_seq <- as.character(IQ_seq[[1]])
length(IQ_seq) #1371

# coding sequence coordinates
IQ_cds_coords <- fread("~/workspace/heterodichogamy/HC_locus_structure/Pawnee_genes/IQ-DOMAIN-23-like/Pawnee_IQ_cds.gff")
IQ_cds_coords <- IQ_cds_coords[, c(4,5)]
setnames(IQ_cds_coords, c("start", "end"))
setkey(IQ_cds_coords, "start")

# combine positions and CDS sequence
IQ_cds <- data.table(IQ_cds_coords[, seq(start,end), by = start][, .(pos = rev(V1))], refseq=IQ_seq)
IQ_cds[, seq_rank := seq(1, .N)]

IQ_cds[seq_rank %% 3 == 1, codon_position := 1]
IQ_cds[seq_rank %% 3 == 2, codon_position := 2]
IQ_cds[seq_rank %% 3 == 0, codon_position := 3]
# verify positions against annotation

# are phasing errors evident for this gene? No. Assigned haplotypes should be correct here
haps[pos >=6527616 & pos <= 6531831]
ggplot(haps[genotype == 'Hh'], 
       aes(x = SNPrank, y = hap.id)) + 
  facet_wrap(~run) + 
  geom_tile(aes(fill = as.factor(allele)), width = 1, height = 1) + 
  scale_fill_manual(values = c('cornsilk', 'chocolate4')) + 
  theme_classic() + 
  labs(x = '', y = '', fill = 'Allele') + 
  theme(axis.text = element_blank(),
        legend.position = 'none', aspect.ratio = .1)  + 
  geom_vline(xintercept = 224, color = 'turquoise2') + 
  geom_vline(xintercept = 284, color = 'turquoise2') 


# get allele counts by haplotype at each site
IQ_cnts <- haps[pos >=6527616 & pos <=  6531831, .(count = sum(allele)), by = .(pos, haplotype, allele0, allele1)]
IQ_cnts <- dcast(IQ_cnts, pos + allele0 + allele1 ~ haplotype, value.var = 'count')

# **reverse complement these to be on the correct strand**
IQ_cnts[, c('allele0', 'allele1') := lapply(.SD, comp), .SDcols = c('allele0', 'allele1')]


# how many of each haplotype?
haps[pos == min(pos), table(haplotype)]
# 47, 19

IQ_final <- merge(IQ_cnts, IQ_cds, by = 'pos')
setkey(IQ_final, seq_rank)

# calculate 
vals <- vector()
for(i in 1:nrow(IQ_final)){
  vals[i] <- syn_vs_nonsyn(CDS_dt = IQ_cds, rnk = IQ_final$seq_rank[i], nuc = IQ_final$allele1[i], codon_position = IQ_final$codon_position[i])
}

IQ_final[, type := vals]

# make sure you double check final numbers when assigning these classes
IQ_final[(H >= 18 & h <= 1) | (H <= 1 & h >= 46), fixed_poly := 'fixed']
IQ_final[h > 1 & h < 46 & H > 1 & H < 18, fixed_poly := 'polymorphic_both']

IQ_final[h > 1  & h < 46 & (H >= 18 | H <= 1), fixed_poly := 'polymorphic_h']
IQ_final[H < 18 & H > 1 & (h <= 1 | h >= 46), fixed_poly := 'polymorphic_H']
IQ_final[, fixed_poly_MK := fixed_poly]
IQ_final[fixed_poly == 'polymorphic_both', fixed_poly_MK := 'polymorphic_h']

IQ_tbl <- table(IQ_final$fixed_poly, IQ_final$type)
IQ_MK_tbl <- table(IQ_final$fixed_poly_MK, IQ_final$type)
fisher.test(IQ_MK_tbl)



# ----- DETOX49 -----

# sequence
DX_seq <- read.fasta("~/workspace/heterodichogamy/HC_locus_structure/Pawnee_genes/protein-DETOXIFICATION-49-like/Pawnee_DETOX_cds.fasta")
DX_seq <- as.character(DX_seq[[1]])
length(DX_seq) #1629

# coding sequence coordinates
DX_cds_coords <- fread("~/workspace/heterodichogamy/HC_locus_structure/Pawnee_genes/protein-DETOXIFICATION-49-like/Pawnee_DETOX_cds.gff")
DX_cds_coords <- DX_cds_coords[, c(4,5)]
setnames(DX_cds_coords, c("start", "end"))
setkey(DX_cds_coords, "start")

# combine positions and CDS sequence
DX_cds <- data.table(DX_cds_coords[, seq(start,end), by = start][, .(pos = rev(V1))], refseq=DX_seq)
DX_cds[, seq_rank := seq(1, .N)]

DX_cds[seq_rank %% 3 == 1, codon_position := 1]
DX_cds[seq_rank %% 3 == 2, codon_position := 2]
DX_cds[seq_rank %% 3 == 0, codon_position := 3]
# verify positions against annotation

# are phasing errors evident for this gene? No. Assigned haplotypes should be correct here
haps[pos >= 6656922 & pos <= 6658550]


# get allele counts by haplotype at each site
DX_cnts <- haps[pos >= 6656922 & pos <= 6658550, .(count = sum(allele)), by = .(pos, haplotype, allele0, allele1)]
DX_cnts <- dcast(DX_cnts, pos + allele0 + allele1 ~ haplotype, value.var = 'count')

# **reverse complement these to be on the correct strand**
DX_cnts[, c('allele0', 'allele1') := lapply(.SD, comp), .SDcols = c('allele0', 'allele1')]


# how many of each haplotype?
haps[pos == min(pos), table(haplotype)]
# 47, 19

DX_final <- merge(DX_cnts, DX_cds, by = 'pos')
setkey(DX_final, seq_rank)

# calculate 
vals <- vector()
for(i in 1:nrow(DX_final)){
  vals[i] <- syn_vs_nonsyn(CDS_dt = DX_cds, rnk = DX_final$seq_rank[i], nuc = DX_final$allele1[i], codon_position = DX_final$codon_position[i])
}

DX_final[, type := vals]

# make sure you double check final numbers when assigning these classes
DX_final[(H >= 18 & h <= 1) | (H <= 1 & h >= 46), fixed_poly := 'fixed']
DX_final[h > 1 & h < 46 & H > 1 & H < 18, fixed_poly := 'polymorphic_both']

DX_final[h > 1  & h < 46 & (H >= 18 | H <= 1), fixed_poly := 'polymorphic_h']
DX_final[H < 18 & H > 1 & (h <= 1 | h >= 46), fixed_poly := 'polymorphic_H']
DX_final[, fixed_poly_MK := fixed_poly]
DX_final[fixed_poly == 'polymorphic_both', fixed_poly_MK := 'polymorphic_h']

DX_tbl <- table(DX_final$fixed_poly, DX_final$type)
DX_MK_tbl <- table(DX_final$fixed_poly_MK, DX_final$type)
fisher.test(DX_MK_tbl)




# ----- CRR1 -----

# sequence
CR_seq <- read.fasta("~/workspace/heterodichogamy/HC_locus_structure/Pawnee_genes/dihydrodipicolinate-reductase-like-protein-CRR1-chloroplastic/Pawnee_CRR1_cds.fasta")
CR_seq <- as.character(CR_seq[[1]])
length(CR_seq) #903

# coding sequence coordinates
CR_cds_coords <- fread("~/workspace/heterodichogamy/HC_locus_structure/Pawnee_genes/dihydrodipicolinate-reductase-like-protein-CRR1-chloroplastic/Pawnee_CRR1_cds.gff")
CR_cds_coords <- CR_cds_coords[, c(4,5)]
setnames(CR_cds_coords, c("start", "end"))
setkey(CR_cds_coords, "start")

# combine positions and CDS sequence
CR_cds <- data.table(CR_cds_coords[, seq(start,end), by = start][, .(pos = rev(V1))], refseq=CR_seq)
CR_cds[, seq_rank := seq(1, .N)]

CR_cds[seq_rank %% 3 == 1, codon_position := 1]
CR_cds[seq_rank %% 3 == 2, codon_position := 2]
CR_cds[seq_rank %% 3 == 0, codon_position := 3]
# verify positions against annotation

# are phasing errors evident for this gene? No. Assigned haplotypes should be correct here
haps[pos >= 6622893 & pos <= 6625554]


# get allele counts by haplotype at each site
CR_cnts <- haps[pos >= 6622893 & pos <= 6625554, .(count = sum(allele)), by = .(pos, haplotype, allele0, allele1)]
CR_cnts <- dcast(CR_cnts, pos + allele0 + allele1 ~ haplotype, value.var = 'count')

# **reverse complement these to be on the correct strand**
CR_cnts[, c('allele0', 'allele1') := lapply(.SD, comp), .SDcols = c('allele0', 'allele1')]


# how many of each haplotype?
haps[pos == min(pos), table(haplotype)]
# 47, 19

CR_final <- merge(CR_cnts, CR_cds, by = 'pos')
setkey(CR_final, seq_rank)

# calculate 
vals <- vector()
for(i in 1:nrow(CR_final)){
  vals[i] <- syn_vs_nonsyn(CDS_dt = CR_cds, rnk = CR_final$seq_rank[i], nuc = CR_final$allele1[i], codon_position = CR_final$codon_position[i])
}

CR_final[, type := vals]

# make sure you double check final numbers when assigning these classes
CR_final[(H >= 18 & h <= 1) | (H <= 1 & h >= 46), fixed_poly := 'fixed']
CR_final[h > 1 & h < 46 & H > 1 & H < 18, fixed_poly := 'polymorphic_both']

CR_final[h > 1  & h < 46 & (H >= 18 | H <= 1), fixed_poly := 'polymorphic_h']
CR_final[H < 18 & H > 1 & (h <= 1 | h >= 46), fixed_poly := 'polymorphic_H']
CR_final[, fixed_poly_MK := fixed_poly]
CR_final[fixed_poly == 'polymorphic_both', fixed_poly_MK := 'polymorphic_h']

CR_tbl <- table(CR_final$fixed_poly, CR_final$type)
CR_MK_tbl <- table(CR_final$fixed_poly_MK, CR_final$type)
fisher.test(CR_MK_tbl)




# ----- CCB2 -----

# sequence
CCB2_seq <- read.fasta("~/workspace/heterodichogamy/HC_locus_structure/Pawnee_genes/protein-COFACTOR-ASSEMBLY-OF-COMPLEX-C-SUBUNIT-B-CCB2-chloroplastic-like/Pawnee_CCB2.fasta")
CCB2_seq <- as.character(CCB2_seq[[1]])
length(CCB2_seq) #867

# coding sequence coordinates
CCB2_cds_coords <- fread("~/workspace/heterodichogamy/HC_locus_structure/Pawnee_genes/protein-COFACTOR-ASSEMBLY-OF-COMPLEX-C-SUBUNIT-B-CCB2-chloroplastic-like/Pawnee_CCB2_cds.gff")
CCB2_cds_coords <- CCB2_cds_coords[, c(4,5)]
setnames(CCB2_cds_coords, c("start", "end"))
setkey(CCB2_cds_coords, "start")

# combine positions and CDS sequence
CCB2_cds <- data.table(CCB2_cds_coords[, seq(start,end), by = start][, .(pos = V1)], refseq=CCB2_seq)
CCB2_cds[, seq_rank := seq(1, .N)]

CCB2_cds[seq_rank %% 3 == 1, codon_position := 1]
CCB2_cds[seq_rank %% 3 == 2, codon_position := 2]
CCB2_cds[seq_rank %% 3 == 0, codon_position := 3]
# verify positions against annotation

# are phasing errors evident for this gene? No. Assigned haplotypes should be correct here
haps[pos >= 6618426 & pos <= 6622821]


# get allele counts by haplotype at each site
CCB2_cnts <- haps[pos >= 6618426 & pos <= 6622821, .(count = sum(allele)), by = .(pos, haplotype, allele0, allele1)]
CCB2_cnts <- dcast(CCB2_cnts, pos + allele0 + allele1 ~ haplotype, value.var = 'count')

# how many of each haplotype?
haps[pos == min(pos), table(haplotype)]
# 47, 19

CCB2_final <- merge(CCB2_cnts, CCB2_cds, by = 'pos')
setkey(CCB2_final, seq_rank)

# calculate 
vals <- vector()
for(i in 1:nrow(CCB2_final)){
  vals[i] <- syn_vs_nonsyn(CDS_dt = CCB2_cds, rnk = CCB2_final$seq_rank[i], nuc = CCB2_final$allele1[i], codon_position = CCB2_final$codon_position[i])
}

CCB2_final[, type := vals]

# make sure you double check final numbers when assigning these classes
CCB2_final[(H >= 18 & h <= 1) | (H <= 1 & h >= 46), fixed_poly := 'fixed']
CCB2_final[h > 1 & h < 46 & H > 1 & H < 18, fixed_poly := 'polymorphic_both']

CCB2_final[h > 1  & h < 46 & (H >= 18 | H <= 1), fixed_poly := 'polymorphic_h']
CCB2_final[H < 18 & H > 1 & (h <= 1 | h >= 46), fixed_poly := 'polymorphic_H']
CCB2_final[, fixed_poly_MK := fixed_poly]
CCB2_final[fixed_poly == 'polymorphic_both', fixed_poly_MK := 'polymorphic_h']

CCB2_tbl <- table(CCB2_final$fixed_poly, CCB2_final$type)
CCB2_MK_tbl <- table(CCB2_final$fixed_poly_MK, CCB2_final$type)
fisher.test(CCB2_MK_tbl)


# ----- PPI2 -----

# sequence
PPI2_seq <- read.fasta("~/workspace/heterodichogamy/HC_locus_structure/Pawnee_genes/Protein_phosphatase_inhibitor_2-like/Pawnee_PPI2_cds.fasta")
PPI2_seq <- as.character(PPI2_seq[[1]])
length(PPI2_seq) #564

# coding sequence coordinates
PPI2_cds_coords <- fread("~/workspace/heterodichogamy/HC_locus_structure/Pawnee_genes/Protein_phosphatase_inhibitor_2-like/Pawnee_PPI2_cds.gff")
PPI2_cds_coords <- PPI2_cds_coords[, c(4,5)]
setnames(PPI2_cds_coords, c("start", "end"))
setkey(PPI2_cds_coords, "start")

# combine positions and CDS sequence
PPI2_cds <- data.table(PPI2_cds_coords[, seq(start,end), by = start][, .(pos = V1)], refseq=PPI2_seq)
PPI2_cds[, seq_rank := seq(1, .N)]

PPI2_cds[seq_rank %% 3 == 1, codon_position := 1]
PPI2_cds[seq_rank %% 3 == 2, codon_position := 2]
PPI2_cds[seq_rank %% 3 == 0, codon_position := 3]
# verify positions against annotation

# are phasing errors evident for this gene? Yes
haps[pos >= 6462468 & pos <= 6464730]


# get rid of individuals with ambiguous phasing
haps_PPI2 <- haps[!run %in% c("SRR15911530", "SRR15911531", "SRR15911546", "SRR15911528")]


# get allele counts by haplotype at each site
PPI2_cnts <- haps_PPI2[pos >= 6462468 & pos <= 6464730, .(count = sum(allele)), by = .(pos, haplotype, allele0, allele1)]
PPI2_cnts <- dcast(PPI2_cnts, pos + allele0 + allele1 ~ haplotype, value.var = 'count')

# how many of each haplotype?
haps_PPI2[pos == min(pos), table(haplotype)]
# 43, 15

PPI2_final <- merge(PPI2_cnts, PPI2_cds, by = 'pos')
setkey(PPI2_final, seq_rank)

# calculate 
vals <- vector()
for(i in 1:nrow(PPI2_final)){
  vals[i] <- syn_vs_nonsyn(CDS_dt = PPI2_cds, rnk = PPI2_final$seq_rank[i], nuc = PPI2_final$allele1[i], codon_position = PPI2_final$codon_position[i])
}

PPI2_final[, type := vals]

# make sure you double check final numbers when assigning these classes
PPI2_final[(H >= 14 & h <= 1) | (H <= 1 & h >= 42), fixed_poly := 'fixed']
PPI2_final[h > 1 & h < 42 & H > 1 & H < 14, fixed_poly := 'polymorphic_both']

PPI2_final[h > 1  & h < 42 & (H >= 14 | H <= 1), fixed_poly := 'polymorphic_h']
PPI2_final[H < 14 & H > 1 & (h <= 1 | h >= 42), fixed_poly := 'polymorphic_H']
PPI2_final[, fixed_poly_MK := fixed_poly]
PPI2_final[fixed_poly == 'polymorphic_both', fixed_poly_MK := 'polymorphic_h']

PPI2_tbl <- table(PPI2_final$fixed_poly, PPI2_final$type)
PPI2_MK_tbl <- table(PPI2_final$fixed_poly_MK, PPI2_final$type)
fisher.test(PPI2_MK_tbl)






# ----- nsLTP14 -----

# sequence
nsLTP14_seq <- read.fasta("~/workspace/heterodichogamy/HC_locus_structure/Pawnee_genes/non-specific-lipid-transfer-protein14/Pawnee_nsLTP14.fasta")
nsLTP14_seq <- as.character(nsLTP14_seq[[1]])
length(nsLTP14_seq) #378

# coding sequence coordinates
nsLTP14_cds_coords <- fread("~/workspace/heterodichogamy/HC_locus_structure/Pawnee_genes/non-specific-lipid-transfer-protein14/Pawnee_nsLTP14_cds.gff")
nsLTP14_cds_coords <-nsLTP14_cds_coords[, c(4,5)]
setnames(nsLTP14_cds_coords, c("start", "end"))
setkey(nsLTP14_cds_coords, "start")

# combine positions and CDS sequence
nsLTP14_cds <- data.table(nsLTP14_cds_coords[, seq(start,end), by = start][, .(pos = V1)], refseq=nsLTP14_seq)
nsLTP14_cds[, seq_rank := seq(1, .N)]

nsLTP14_cds[seq_rank %% 3 == 1, codon_position := 1]
nsLTP14_cds[seq_rank %% 3 == 2, codon_position := 2]
nsLTP14_cds[seq_rank %% 3 == 0, codon_position := 3]
# verify positions against annotation

# are phasing errors evident for this gene? No
haps[pos >= 6524724 & pos <= 6525368]

# are phasing errors evident for this gene? No. 
ggplot(haps[genotype == 'Hh'], 
       aes(x = SNPrank, y = hap.id)) + 
  facet_wrap(~run) + 
  geom_tile(aes(fill = as.factor(allele)), width = 1, height = 1) + 
  scale_fill_manual(values = c('cornsilk', 'chocolate4')) + 
  theme_classic() + 
  labs(x = '', y = '', fill = 'Allele') + 
  theme(axis.text = element_blank(),
        legend.position = 'none', aspect.ratio = .1)  + 
  geom_vline(xintercept = 215, color = 'turquoise2') + 
  geom_vline(xintercept = 223, color = 'turquoise2') 



# get allele counts by haplotype at each site
nsLTP14_cnts <- haps[pos >= 6524724 & pos <= 6525368, .(count = sum(allele)), by = .(pos, haplotype, allele0, allele1)]
nsLTP14_cnts <- dcast(nsLTP14_cnts, pos + allele0 + allele1 ~ haplotype, value.var = 'count')

# how many of each haplotype?
haps[pos == min(pos), table(haplotype)]
# 47, 19

nsLTP14_final <- merge(nsLTP14_cnts, nsLTP14_cds, by = 'pos')
setkey(nsLTP14_final, seq_rank)

# calculate 
vals <- vector()
for(i in 1:nrow(nsLTP14_final)){
  vals[i] <- syn_vs_nonsyn(CDS_dt = nsLTP14_cds, rnk = nsLTP14_final$seq_rank[i], nuc = nsLTP14_final$allele1[i], codon_position = nsLTP14_final$codon_position[i])
}

nsLTP14_final[, type := vals]

# make sure you double check final numbers when assigning these classes
nsLTP14_final[(H >= 18 & h <= 1) | (H <= 1 & h >= 46), fixed_poly := 'fixed']
nsLTP14_final[h > 1 & h < 46 & H > 1 & H < 18, fixed_poly := 'polymorphic_both']

nsLTP14_final[h > 1  & h < 46 & (H >= 18 | H <= 1), fixed_poly := 'polymorphic_h']
nsLTP14_final[H < 18 & H > 1 & (h <= 1 | h >= 46), fixed_poly := 'polymorphic_H']
nsLTP14_final[, fixed_poly_MK := fixed_poly]
nsLTP14_final[fixed_poly == 'polymorphic_both', fixed_poly_MK := 'polymorphic_h']

nsLTP14_tbl <- table(nsLTP14_final$fixed_poly, nsLTP14_final$type)
nsLTP14_MK_tbl <- table(nsLTP14_final$fixed_poly_MK, nsLTP14_final$type)
fisher.test(nsLTP14_MK_tbl)





# ----- DNA-directed RNA polymerases II, IV and V subunit 6A-like -----

# sequence
RNApol_seq <- read.fasta("~/workspace/heterodichogamy/HC_locus_structure/Pawnee_genes/DNA-directed-RNA-polymerase-II-IV-V-subunit-6A-like/Pawnee_RNApol_cds.fasta")
RNApol_seq <- as.character(RNApol_seq[[1]])
length(RNApol_seq) #429

# coding sequence coordinates
RNApol_cds_coords <- fread("~/workspace/heterodichogamy/HC_locus_structure/Pawnee_genes/DNA-directed-RNA-polymerase-II-IV-V-subunit-6A-like/Pawnee_RNApol_cds.gff")
RNApol_cds_coords <- RNApol_cds_coords[, c(4,5)]
setnames(RNApol_cds_coords, c("start", "end"))
setkey(RNApol_cds_coords, "start")

# combine positions and CDS sequence
RNApol_cds <- data.table(RNApol_cds_coords[, seq(start,end), by = start][, .(pos = V1)], refseq=RNApol_seq)
RNApol_cds[, seq_rank := seq(1, .N)]

RNApol_cds[seq_rank %% 3 == 1, codon_position := 1]
RNApol_cds[seq_rank %% 3 == 2, codon_position := 2]
RNApol_cds[seq_rank %% 3 == 0, codon_position := 3]
# verify positions against annotation

# are phasing errors evident for this gene? No
haps[pos >= 6564067 & pos <= 6567952]


# get allele counts by haplotype at each site
RNApol_cnts <- haps[pos >= 6564067 & pos <= 6567952, .(count = sum(allele)), by = .(pos, haplotype, allele0, allele1)]
RNApol_cnts <- dcast(RNApol_cnts, pos + allele0 + allele1 ~ haplotype, value.var = 'count')

# how many of each haplotype?
haps[pos == min(pos), table(haplotype)]
# 47, 19

RNApol_final <- merge(RNApol_cnts, RNApol_cds, by = 'pos')
setkey(RNApol_final, seq_rank)

# calculate 
vals <- vector()
for(i in 1:nrow(RNApol_final)){
  vals[i] <- syn_vs_nonsyn(CDS_dt = RNApol_cds, rnk = RNApol_final$seq_rank[i], nuc = RNApol_final$allele1[i], codon_position = RNApol_final$codon_position[i])
}

RNApol_final[, type := vals]

# make sure you double check final numbers when assigning these classes
RNApol_final[(H >= 18 & h <= 1) | (H <= 1 & h >= 46), fixed_poly := 'fixed']
RNApol_final[h > 1 & h < 46 & H > 1 & H < 18, fixed_poly := 'polymorphic_both']

RNApol_final[h > 1  & h < 46 & (H >= 18 | H <= 1), fixed_poly := 'polymorphic_h']
RNApol_final[H < 18 & H > 1 & (h <= 1 | h >= 46), fixed_poly := 'polymorphic_H']
RNApol_final[, fixed_poly_MK := fixed_poly]
RNApol_final[fixed_poly == 'polymorphic_both', fixed_poly_MK := 'polymorphic_h']

RNApol_tbl <- table(RNApol_final$fixed_poly, RNApol_final$type)
RNApol_MK_tbl <- table(RNApol_final$fixed_poly_MK, RNApol_final$type)
fisher.test(RNApol_MK_tbl)







# ----- LOC LOC122306891 -----

# sequence
loc1_seq <- read.fasta("~/workspace/heterodichogamy/HC_locus_structure/Pawnee_genes/LOC122306891_uncharacterized/Pawnee_uncharacterized1.fasta")
loc1_seq <- as.character(loc1_seq[[1]])
length(loc1_seq) #522

# coding sequence coordinates
loc1_cds_coords <- fread("~/workspace/heterodichogamy/HC_locus_structure/Pawnee_genes/LOC122306891_uncharacterized/Pawnee_uncharacterized1_cds.gff")
loc1_cds_coords <- loc1_cds_coords[, c(4,5)]
setnames(loc1_cds_coords, c("start", "end"))
setkey(loc1_cds_coords, "start")

# combine positions and CDS sequence
loc1_cds <- data.table(loc1_cds_coords[, seq(start,end), by = start][, .(pos = V1)], refseq=loc1_seq)
loc1_cds[, seq_rank := seq(1, .N)]

loc1_cds[seq_rank %% 3 == 1, codon_position := 1]
loc1_cds[seq_rank %% 3 == 2, codon_position := 2]
loc1_cds[seq_rank %% 3 == 0, codon_position := 3]
# verify positions against annotation

# are phasing errors evident for this gene? Yes 
haps[pos >= 6468176 & pos <= 6469520]

ggplot(haps[genotype == 'Hh'], 
       aes(x = SNPrank, y = hap.id)) + 
  facet_wrap(~run) + 
  geom_tile(aes(fill = as.factor(allele)), width = 1, height = 1) + 
  scale_fill_manual(values = c('cornsilk', 'chocolate4')) + 
  theme_classic() + 
  labs(x = '', y = '', fill = 'Allele') + 
  theme(axis.text = element_blank(),
        legend.position = 'none', aspect.ratio = .1)  + 
  geom_vline(xintercept = 7, color = 'turquoise2') + 
  geom_vline(xintercept = 9, color = 'turquoise2') 

# get rid of individuals with ambiguous phasing
haps_loc1 <- haps[!run %in% c("SRR15911530", "SRR15911531", "SRR15911546", "SRR15911528")]


# get allele counts by haplotype at each site
loc1_cnts <- haps[pos >= 6468176 & pos <= 6469520, .(count = sum(allele)), by = .(pos, haplotype, allele0, allele1)]
loc1_cnts <- dcast(loc1_cnts, pos + allele0 + allele1 ~ haplotype, value.var = 'count')

# how many of each haplotype?
haps_loc1[pos == min(pos), table(haplotype)]
# 43, 15

loc1_final <- merge(loc1_cnts, loc1_cds, by = 'pos')
setkey(loc1_final, seq_rank)

# calculate 
vals <- vector()
for(i in 1:nrow(loc1_final)){
  vals[i] <- syn_vs_nonsyn(CDS_dt = loc1_cds, rnk = loc1_final$seq_rank[i], nuc = loc1_final$allele1[i], codon_position = loc1_final$codon_position[i])
}

loc1_final[, type := vals]

# make sure you double check final numbers when assigning these classes
loc1_final[(H >= 14 & h <= 1) | (H <= 1 & h >= 42), fixed_poly := 'fixed']
loc1_final[h > 1 & h < 42 & H > 1 & H < 14, fixed_poly := 'polymorphic_both']

loc1_final[h > 1  & h < 42 & (H >= 14 | H <= 1), fixed_poly := 'polymorphic_h']
loc1_final[H < 14 & H > 1 & (h <= 1 | h >= 42), fixed_poly := 'polymorphic_H']
loc1_final[, fixed_poly_MK := fixed_poly]
loc1_final[fixed_poly == 'polymorphic_both', fixed_poly_MK := 'polymorphic_h']

loc1_tbl <- table(loc1_final$fixed_poly, loc1_final$type)
loc1_MK_tbl <- table(loc1_final$fixed_poly_MK, loc1_final$type)
fisher.test(loc1_MK_tbl)







# ----- LOC122306943 -----

# sequence
loc2_seq <- read.fasta("~/workspace/heterodichogamy/HC_locus_structure/Pawnee_genes/LOC122306943/Pawnee_loc2_cds.fasta")
loc2_seq <- as.character(loc2_seq[[1]])
length(loc2_seq) #771

# coding sequence coordinates
# note that we left off 2nd exon as reading frame gets ambiguous
loc2_cds_coords <- fread("~/workspace/heterodichogamy/HC_locus_structure/Pawnee_genes/LOC122306943/Pawnee_loc2_cds.gff")
loc2_cds_coords <- loc2_cds_coords[, c(4,5)]
setnames(loc2_cds_coords, c("start", "end"))
setkey(loc2_cds_coords, "start")

# combine positions and CDS sequence
loc2_cds <- data.table(loc2_cds_coords[1, seq(start,end), by = start][, .(pos = V1)], refseq=loc2_seq)
loc2_cds[, seq_rank := seq(1, .N)]

loc2_cds[seq_rank %% 3 == 1, codon_position := 1]
loc2_cds[seq_rank %% 3 == 2, codon_position := 2]
loc2_cds[seq_rank %% 3 == 0, codon_position := 3]
# verify positions against annotation

# are phasing errors evident for this gene? No
haps[pos >= 6522904 & pos <= 6523674]

ggplot(haps[genotype == 'Hh'], 
       aes(x = SNPrank, y = hap.id)) + 
  facet_wrap(~run) + 
  geom_tile(aes(fill = as.factor(allele)), width = 1, height = 1) + 
  scale_fill_manual(values = c('cornsilk', 'chocolate4')) + 
  theme_classic() + 
  labs(x = '', y = '', fill = 'Allele') + 
  theme(axis.text = element_blank(),
        legend.position = 'none', aspect.ratio = .1)  + 
  geom_vline(xintercept = 196, color = 'turquoise2') + 
  geom_vline(xintercept = 214, color = 'turquoise2') 


# get allele counts by haplotype at each site
loc2_cnts <- haps[pos >= 6522904 & pos <= 6523674, .(count = sum(allele)), by = .(pos, haplotype, allele0, allele1)]
loc2_cnts <- dcast(loc2_cnts, pos + allele0 + allele1 ~ haplotype, value.var = 'count')

# how many of each haplotype?
haps[pos == min(pos), table(haplotype)]
# 47,19

loc2_final <- merge(loc2_cnts, loc2_cds, by = 'pos')
setkey(loc2_final, seq_rank)

# calculate 
vals <- vector()
for(i in 1:nrow(loc2_final)){
  vals[i] <- syn_vs_nonsyn(CDS_dt = loc2_cds, rnk = loc2_final$seq_rank[i], nuc = loc2_final$allele1[i], codon_position = loc2_final$codon_position[i])
}

loc2_final[, type := vals]

# make sure you double check final numbers when assigning these classes
loc2_final[(H >= 18 & h <= 1) | (H <= 1 & h >= 46), fixed_poly := 'fixed']
loc2_final[h > 1 & h < 46 & H > 1 & H < 18, fixed_poly := 'polymorphic_both']

loc2_final[h > 1  & h < 46 & (H >= 18 | H <= 1), fixed_poly := 'polymorphic_h']
loc2_final[H < 18 & H > 1 & (h <= 1 | h >= 46), fixed_poly := 'polymorphic_H']
loc2_final[, fixed_poly_MK := fixed_poly]
loc2_final[fixed_poly == 'polymorphic_both', fixed_poly_MK := 'polymorphic_h']

loc2_tbl <- table(loc2_final$fixed_poly, loc2_final$type)
loc2_MK_tbl <- table(loc2_final$fixed_poly_MK, loc2_final$type)
fisher.test(loc2_MK_tbl)








