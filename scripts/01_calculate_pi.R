library(data.table)

# calculate pi within FAF coding sequence
# for P. stenoptera, this includes the mapping population and Geng et al. data
# data for P. macroptera and P. hupehensis were phased separately.

# ----- SNP positions ----- 
leg.ste <- fread("~/workspace/Pterocarya/calls/PSTE_UCD1470_HAP2/PsFAFL.impute.legend")
leg.hup <- fread("~/workspace/Pterocarya/Geng_etal_2024_data/P.hupehensis.impute.legend")
leg.mac <- fread("~/workspace/Pterocarya/Geng_etal_2024_data/P.macroptera.impute.legend")

# ID column in this table is both redundant and unecessary
leg.ste[, ID := NULL]
leg.hup[, ID := NULL]
leg.mac[, ID := NULL]

# ----- haps ----- Each column is a biallelic SNP position. Columns are haplotypes.
haps0.ste <- fread("~/workspace/Pterocarya/calls/PSTE_UCD1470_HAP2/PsFAFL.impute.hap")
nc.ste <- ncol(haps0.ste)
setnames(haps0.ste, paste0("V", 1:nc.ste), paste0(rep(seq(1:(nc.ste/2)), each = 2), c("_1", "_2")))

haps0.hup <- fread("~/workspace/Pterocarya/Geng_etal_2024_data/P.hupehensis.impute.hap")
nc.hup <- ncol(haps0.hup)
setnames(haps0.hup, paste0("V", 1:nc.hup), paste0(rep(seq(1:(nc.hup/2)), each = 2), c("_1", "_2")))

haps0.mac <- fread("~/workspace/Pterocarya/Geng_etal_2024_data/P.macroptera.impute.hap")
nc.mac <- ncol(haps0.mac)
setnames(haps0.mac, paste0("V", 1:nc.mac), paste0(rep(seq(1:(nc.mac/2)), each = 2), c("_1", "_2")))


d0.ste <- melt(cbind(leg.ste, haps0.ste), id.vars = c("pos", "allele0", "allele1"), value.name = 'allele')
d0.ste[, c("indiv.id", "hap.id") := tstrsplit(variable, split = "_", fixed = T)]
d0.ste[, variable := NULL]

d0.hup <- melt(cbind(leg.hup, haps0.hup), id.vars = c("pos", "allele0", "allele1"), value.name = 'allele')
d0.hup[, c("indiv.id", "hap.id") := tstrsplit(variable, split = "_", fixed = T)]
d0.hup[, variable := NULL]

d0.mac <- melt(cbind(leg.mac, haps0.mac), id.vars = c("pos", "allele0", "allele1"), value.name = 'allele')
d0.mac[, c("indiv.id", "hap.id") := tstrsplit(variable, split = "_", fixed = T)]
d0.mac[, variable := NULL]


# ----- individual IDs -----
indiv.ste <- fread("~/workspace/Pterocarya/calls/PSTE_UCD1470_HAP2/PsFAFL.impute.hap.indv", header = F, col.names = c("ID"))
indiv.ste[, indiv.id := as.character(seq(1, .N))]

indiv.hup <- fread("~/workspace/Pterocarya/Geng_etal_2024_data/P.hupehensis.impute.hap.indv", header = F, col.names = c("ID"))
indiv.hup[, indiv.id := as.character(seq(1, .N))]

indiv.mac <- fread("~/workspace/Pterocarya/Geng_etal_2024_data/P.macroptera.impute.hap.indv", header = F, col.names = c("ID"))
indiv.mac[, indiv.id := as.character(seq(1, .N))]

d1.ste <- merge(d0.ste, indiv.ste, by = 'indiv.id', all = T)
d1.ste[, indiv.id := NULL]

d1.hup <- merge(d0.hup, indiv.hup, by = 'indiv.id', all = T)
d1.hup[, indiv.id := NULL]

d1.mac <- merge(d0.mac, indiv.mac, by = 'indiv.id', all = T)
d1.mac[, indiv.id := NULL]



# ---- merge with genotype (for stenoptera, assigned by phenotype or coverage) -----
geno.ste <- fread("~/workspace/Pterocarya/P.stenoptera_MK_genotypes.txt", header = F, col.names = c("ID", 'genotype'))
haps.ste <- merge(geno.ste, d1.ste, by = 'ID', all = T)

# ---- subset to FAF (phased across broader region, which resolved phase switching in heterozygotes)----
haps.hup <- d1.hup
haps.hup <- haps.hup[pos >= 4366881 & pos <= 4367828]

haps.mac <- d1.mac
haps.mac <- haps.mac[pos >= 4366881 & pos <= 4367828]

# ----- assign ranks -----
haps.ste[, hapRank := seq(1, .N), by = .(pos)]
haps.hup[, hapRank := seq(1, .N), by = .(pos)]
haps.mac[, hapRank := seq(1, .N), by = .(pos)]


# ------ Assign haplotype identities -----
# first we assign the G haplotypes of heterozygotes 
g.hap.ids.ste <- haps.ste[pos == 4367555 & allele == 1, hapRank]
G.hap.ids.ste <- haps.ste[pos == 4367555 & allele == 0, hapRank]
haps.ste[hapRank %in% g.hap.ids.ste, haplotype := 'g']
haps.ste[hapRank %in% G.hap.ids.ste, haplotype := 'G']

g.hap.ids.hup <- haps.hup[pos == 4367555 & allele == 1, hapRank]
G.hap.ids.hup <- haps.hup[pos == 4367555 & allele == 0, hapRank]
haps.hup[hapRank %in% g.hap.ids.hup, haplotype := 'g']
haps.hup[hapRank %in% G.hap.ids.hup, haplotype := 'G']

g.hap.ids.mac <- haps.mac[pos == 4367555 & allele == 1, hapRank]
G.hap.ids.mac <- haps.mac[pos == 4367555 & allele == 0, hapRank]
haps.mac[hapRank %in% g.hap.ids.mac, haplotype := 'g']
haps.mac[hapRank %in% G.hap.ids.mac, haplotype := 'G']

# confirm that the assignment by the SNP above agrees with the other CDS SNP
# haps.ste[pos == 4367603 & haplotype == 'g']
# haps.ste[pos == 4367603 & haplotype == 'G']
# haps.hup[pos == 4367603 & haplotype == 'g']
# haps.hup[pos == 4367603 & haplotype == 'G']
# haps.mac[pos == 4367603 & haplotype == 'g']
# haps.mac[pos == 4367603 & haplotype == 'G']


# ---- inspect phasing results -----
# assign genotypes

hup.hets <- haps.hup[haplotype == 'G', unique(ID)]
mac.hets <- haps.mac[haplotype == 'G', unique(ID)]

haps.hup[ID %in% hup.hets, genotype := 'Gg']
haps.hup[!ID %in% hup.hets, genotype := 'gg']
haps.mac[ID %in% mac.hets, genotype := 'Gg']
haps.mac[!ID %in% mac.hets, genotype := 'gg']


# look at all haplotypes, check phasing result in hupehensis
haps.hup[, SNPrank := seq(1, .N), by = .(ID, hap.id)]
haps.mac[, SNPrank := seq(1, .N), by = .(ID, hap.id)]


# distribution of non-reference alleles. Not informative for macroptera, given little shared polymorphisms
a.hup <- haps.hup[hap.id == 1, .(nonref_alleles_hap1 = sum(allele)), by = .(ID,genotype)]
b.hup <- haps.hup[hap.id == 2, .(nonref_alleles_hap2 = sum(allele)), by = .(ID,genotype)]
nonref_alleles.hup <- merge(a.hup,b.hup)
nonref_allhaps.hup <- melt(nonref_alleles.hup, id.vars = c("ID", 'genotype'), value.name = "Non ref alleles")

ggplot(nonref_allhaps.hup, aes(x = `Non ref alleles`)) + geom_histogram() + 
  theme_classic() + 
  theme(aspect.ratio = 1)



# specifically look at heterozygotes in hupehensis (expect bimodal)
nonref_alleles.hup[genotype == 'Gg', hap_more := max(nonref_alleles_hap1, nonref_alleles_hap2), by = .(ID,genotype)]
nonref_alleles.hup[genotype == 'Gg', hap_less := min(nonref_alleles_hap1, nonref_alleles_hap2), by = .(ID,genotype)]

nonref.hup <- melt(nonref_alleles.hup[, .(ID, hap_more, hap_less, genotype)], id.vars = c('ID', 'genotype'), value.name = 'nonref', variable.name = 'haplotype')
ggplot(nonref.hup[genotype == 'Gg'], aes(x = nonref)) + geom_histogram()
ggplot(nonref.hup[genotype == 'Gg'], aes(x = haplotype, y = nonref, group = ID)) + geom_line()


#  visualize haplotypes 

# hupehensis
ggplot(haps.hup, 
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
ggplot(haps.hup[ID %in% hup.hets ], #& ID == 'SRR25637677'
       aes(x = SNPrank, y = hap.id)) + 
  facet_wrap(~ID) + 
  geom_tile(aes(fill = as.factor(allele)), width = 1, height = 1) + 
  scale_fill_manual(values = c('cornsilk', 'chocolate4')) + 
  theme_classic() + 
  labs(x = '', y = '', fill = 'Allele') + 
  theme(axis.text = element_blank(),
        legend.position = 'none', aspect.ratio = .5)  

# macroptera
ggplot(haps.mac, 
       aes(x = SNPrank, y = hapRank)) + 
  facet_wrap(~haplotype, scales = 'free') +
  
  geom_tile(aes(fill = as.factor(allele)), width = 1, height = 1) + 
  scale_fill_manual(values = c('cornsilk', 'chocolate4')) + 
  theme_classic() + 
  labs(x = '', y = '', fill = 'Allele') + 
  theme(axis.text = element_blank(),
        legend.position = 'none')

ggplot(haps.mac[ID %in% mac.hets ], #& ID == 'SRR25637677'
       aes(x = SNPrank, y = hap.id)) + 
  facet_wrap(~ID) + 
  geom_tile(aes(fill = as.factor(allele)), width = 1, height = 1) + 
  scale_fill_manual(values = c('cornsilk', 'chocolate4')) + 
  theme_classic() + 
  labs(x = '', y = '', fill = 'Allele') + 
  theme(axis.text = element_blank(),
        legend.position = 'none', aspect.ratio = .5)  



# ------ diversity ---------

haps.FAF.ste <- haps.ste[pos >= 4366881 & pos <= 4367828]
haps.FAF.hup <- haps.hup[pos >= 4366881 & pos <= 4367828]
haps.FAF.mac <- haps.mac[pos >= 4366881 & pos <= 4367828]

g.pairs.ste <- t(combn(g.hap.ids.ste, 2))
G.pairs.ste <- t(combn(G.hap.ids.ste, 2))

g.pairs.hup <- t(combn(g.hap.ids.hup, 2))
G.pairs.hup <- t(combn(G.hap.ids.hup, 2))

g.pairs.mac <- t(combn(g.hap.ids.mac, 2))
G.pairs.mac <- t(combn(G.hap.ids.mac, 2))


# ----- pi in g haps -----

# stenoptera
g.pi.vals.ste <- vector()
for(k in 1:nrow(g.pairs.ste)){
  i <- g.pairs.ste[k, 1]
  j <- g.pairs.ste[k, 2]
  ij <- merge(haps.FAF.ste[hapRank == i, .(pos, i.allele = allele)], haps.FAF.ste[hapRank == j, .(pos, j.allele = allele)])
  g.pi.vals.ste[k] <- ij[, length(which(i.allele != j.allele))/(4367828-4366881+1)]
}
pste.g <- mean(g.pi.vals.ste) #0.003340594

# hupehensis
g.pi.vals.hup <- vector()
for(k in 1:nrow(g.pairs.hup)){
  i <- g.pairs.hup[k, 1]
  j <- g.pairs.hup[k, 2]
  ij <- merge(haps.FAF.hup[hapRank == i, .(pos, i.allele = allele)], haps.FAF.hup[hapRank == j, .(pos, j.allele = allele)])
  g.pi.vals.hup[k] <- ij[, length(which(i.allele != j.allele))/(4367828-4366881+1)]
}
phup.g <- mean(g.pi.vals.hup) # 0.004718067

# macroptera
g.pi.vals.mac <- vector()
for(k in 1:nrow(g.pairs.mac)){
  i <- g.pairs.mac[k, 1]
  j <- g.pairs.mac[k, 2]
  ij <- merge(haps.FAF.mac[hapRank == i, .(pos, i.allele = allele)], haps.FAF.mac[hapRank == j, .(pos, j.allele = allele)])
  g.pi.vals.mac[k] <- ij[, length(which(i.allele != j.allele))/(4367828-4366881+1)]
}
pmac.g <- mean(g.pi.vals.mac) # 0.002329949


# ----- pi in G haps -----

# stenoptera
G.pi.vals.ste <- vector()
for(k in 1:nrow(G.pairs.ste)){
  i <- G.pairs.ste[k, 1]
  j <- G.pairs.ste[k, 2]
  ij <- merge(haps.FAF.ste[hapRank == i, .(pos, i.allele = allele)], haps.FAF.ste[hapRank == j, .(pos, j.allele = allele)])
  G.pi.vals.ste[k] <- ij[, length(which(i.allele != j.allele))/((4367828-4366881+1))]
}
pste.G <- mean(G.pi.vals.ste) #0.001657625

# hupehensis
G.pi.vals.hup <- vector()
for(k in 1:nrow(G.pairs.hup)){
  i <- G.pairs.hup[k, 1]
  j <- G.pairs.hup[k, 2]
  ij <- merge(haps.FAF.hup[hapRank == i, .(pos, i.allele = allele)], haps.FAF.hup[hapRank == j, .(pos, j.allele = allele)])
  G.pi.vals.hup[k] <- ij[, length(which(i.allele != j.allele))/((4367828-4366881+1))]
}
phup.G <- mean(G.pi.vals.hup) # 0.002310629

# macroptera
G.pi.vals.mac <- vector()
for(k in 1:nrow(G.pairs.mac)){
  i <- G.pairs.mac[k, 1]
  j <- G.pairs.mac[k, 2]
  ij <- merge(haps.FAF.mac[hapRank == i, .(pos, i.allele = allele)], haps.FAF.mac[hapRank == j, .(pos, j.allele = allele)])
  G.pi.vals.mac[k] <- ij[, length(which(i.allele != j.allele))/(4367828-4366881+1)]
}
pmac.G <- mean(G.pi.vals.mac) # 0.0003516174


N <- sapply(
  list(g.hap.ids.ste, G.hap.ids.ste, 
     g.hap.ids.hup, G.hap.ids.hup,
     g.hap.ids.mac, G.hap.ids.mac),
  length)


final <- data.table(species = rep(c("P. stenoptera", "P. hupehensis", "P. macroptera"), each = 2),
           haplotype = rep(c("g", "G"), 3),
          N = paste0("N=",N),
           pi = c(pste.g, pste.G, phup.g, phup.G, pmac.g, pmac.G))

final[, species := factor(species, levels = c("P. stenoptera", "P. hupehensis", "P. macroptera"))]
final[, textpos := c(rep(0.00025, 2), rep(0.0003, 2), rep(0.00015, 2))]

final

ggplot(final, aes(x = haplotype, y = pi)) + 
  geom_bar(stat = 'identity', color = 'black', fill = 'gray') + 
  #geom_text(aes(label = N, y = textpos), size = 2.5) +
  scale_y_continuous(position = 'right') +
  labs(y = expression(paste("Pairwise diversity (", pi, ")")), x = '') +
  facet_wrap(~species, ncol = 1, scales = 'free_y') + 
  theme_classic() + 
  theme(aspect.ratio = 0.7, 
        strip.text = element_text(face = 'italic'),
        axis.text.x = element_text(face = 'italic', size = 12), 
        axis.title.y = element_text(size = 10, vjust=-3))
