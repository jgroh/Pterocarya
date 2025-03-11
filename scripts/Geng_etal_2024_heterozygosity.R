library(data.table)
library(ggplot2)

sra <- fread("~/workspace/Pterocarya/Geng_etal_2024_data/Geng_etal_2024_SRA_samples.txt", header = F, col.names = c("ID", "run", "species"))
#het <- fread("~/workspace/Pterocarya/Geng_etal_2024_data/Geng_etal_het_fltvcf.txt")

# Pste - Pmac TSP ascertained snps
het1 <- fread("~/workspace/Pterocarya/Geng_etal_2024_data/Geng_etal_het_fltvcf_Pste_Pmac_sites.txt")
het1[, SNPset := 'Pste_Pmac']

# Pste G-loc SNPs
het2 <- fread("~/workspace/Pterocarya/Geng_etal_2024_data/Geng_etal_het_Pste_Gloc_SNPs.txt")
het2[, SNPset := 'Pste']

het <- rbind(het1, het2)

setnames(het, c("V1", "V2", "V3", "V4"), c("chr", "pos", "G", "g"))

het <- melt(het, id.vars = c("chr", "pos", "G", 'g', 'SNPset'), value.name = 'sample_gt')
het[, variable := NULL]
het[, c("run", "geno") := tstrsplit(sample_gt, split = '=')]
het[, sample_gt := NULL]

het[geno== '0/0', het := 0]
het[geno == '0/1' | geno == '1/0', het := 1]
het[geno == '1/1', het := 0]
het[geno == './.', het := NA]

het <- merge(het, sra, by = 'run')


# ------- combine with frax, rho data ------
het.fra.rho <- fread("~/workspace/Pterocarya/calls/PSTE_UCD1470_HAP2/nonPste_samples_TSP_SNPs_het.txt")
het.fra.rho[, SNPset := 'Pste_Pmac']
setnames(het.fra.rho, c("V1", "V2", "V3", "V4"), c("chr", "pos", "G", "g"))

het.fra.rho <- melt(het.fra.rho, id.vars = c("chr", "pos", "G", 'g', 'SNPset'), value.name = 'sample_gt')
het.fra.rho[, variable := NULL]
het.fra.rho[, c("run", "geno") := tstrsplit(sample_gt, split = '=')]
het.fra.rho[, sample_gt := NULL]

het.fra.rho[geno== '0/0', het := 0]
het.fra.rho[geno == '0/1' | geno == '1/0', het := 1]
het.fra.rho[geno == '1/1', het := 0]
het.fra.rho[geno == './.', het := NA]
het.fra.rho[, ID := NA]
het.fra.rho[grepl("RHO", run), species := 'rhoifolia']
het.fra.rho[grepl("FRA", run), species := 'fraxinifolia']
het.fra.rho[, ID := run]

het <- rbind(het, het.fra.rho[, c("run", 'chr', 'pos', 'G', 'g', 'SNPset', 'geno', 'het', 'ID', 'species')])



# subest to 2 coding SNPs
het_2cds_SNPs <- het[SNPset == 'Pste_Pmac' & pos %in% c(4367555,4367603)]

het[, nsites := sum(!is.na(het)), by = .(ID, SNPset)]
het_2cds_SNPs[, nsites := sum(!is.na(het)), by = ID]

# code individuals as heterozygote based on 2 CDS SNPs
het_2cds_SNPs[, type := ifelse(sum(het) == 2, 'het', 'hom_gg'), by = ID]

# inspecting genotypes, this individual appears GG homozygote
het_2cds_SNPs[run == 'SRR25637527', type := 'hom_GG']

# write out genotypes 
Geng_etal_genotypes <- het_2cds_SNPs[pos == 4367603 & species == 'stenoptera', .(run, type)]
Geng_etal_genotypes[type == 'het', geno := 'Gg']
Geng_etal_genotypes[type == 'hom_GG', geno := 'GG']
Geng_etal_genotypes[type == 'hom_gg', geno := 'gg']
Geng_etal_genotypes[,type := NULL]
# fwrite(Geng_etal_genotypes, file = '~/workspace/Pterocarya/Geng_etal_2024_data/Geng_etal_stenoptera_Gloc_genotypes.txt',
#        quote = F, col.names = F, row.names = F, sep = '\t')

het[run %in% het_2cds_SNPs[type == 'het', run], type := 'het']
het[run %in% het_2cds_SNPs[type == 'hom_GG', run], type := 'hom']
het[run %in% het_2cds_SNPs[type == 'hom_gg', run], type := 'hom']

het[run %in% het_2cds_SNPs[type == 'het', run], Gloc_geno := '01']
het[run %in% het_2cds_SNPs[type == 'hom_GG', run], Gloc_geno := '00']
het[run %in% het_2cds_SNPs[type == 'hom_gg', run], Gloc_geno := '11']

# check hets per species
het[species == 'stenoptera', 
    .(het = sum(het, na.rm = T), 
      N = sum(!is.na(het))), by = .(run,type, SNPset)]

het[species == 'macroptera', 
    .(het = sum(het, na.rm = T), 
      N = sum(!is.na(het))), by = .(run,type, SNPset)]

het[species == 'hupehensis', 
    .(het = sum(het, na.rm = T), 
      N = sum(!is.na(het))), by = .(run,type, SNPset)]


# --- check boundaries of heterozygosity ----

het[species == 'stenoptera' & Gloc_geno == '01' & pos ==4363775 ]
het[species == 'hupehensis' & Gloc_geno == '01' & pos ==4363775 ]
het[species == 'macroptera' & Gloc_geno == '01' & pos ==4363775 ]


# define heterozygosity
het_p <- het[, .(X = length(which(geno == '0/1')), N = length(which(geno!='./.'))), by = .(ID, species, SNPset)]
het_p[, p := X/N, by = SNPset]
het_p[SNPset == 'Pste_Pmac']

het_p <- merge(het_p[SNPset== 'Pste_Pmac'], 
               unique(het[, .(ID, species, Gloc_geno)]), by = c('ID', 'species'))


# ----- make plot ------

het_p[, species := paste0("P. ", species)]
#het_p$species = factor(het_p$species, levels=c("P. stenoptera","P. hupehensis",'P. fraxinifolia', 'P. macroptera', 'P. rhoifolia'))
het_p$species = factor(het_p$species, levels=c("P. fraxinifolia","P. stenoptera","P. hupehensis", 'P. macroptera', 'P. rhoifolia'))

het_p[, y := as.numeric(as.factor(ID)), by = species]

# Agresti-Coull 95% CI
het_p[, pprime := (X+2)/(N+4)]
het_p[, CI.lower := pprime - 1.96*sqrt((pprime*(1-pprime))/(N+4))]
het_p[, CI.upper := pprime + 1.96*sqrt((pprime*(1-pprime))/(N+4))]

ggplot(het_p[SNPset == 'Pste_Pmac' & species != 'P. rhoifolia'], aes(x = p, y = y, color = Gloc_geno)) +

  facet_wrap(~species, scales = 'free_y', ncol = 1) +
  scale_color_manual(
    values = c('11' = '#DDAA33', '01' = '#004488', 
               '00' = '#DC267F'),
    labels = c('11' = 'gg', '01' = 'Gg','00' = 'GG')) +
  geom_vline(xintercept = 0, linetype = 2, color = 'gray') + 
  geom_vline(xintercept = 1, linetype = 2, color = 'gray') +
    
  #geom_point(position=position_dodge(width=10), size = 2)  + 
  geom_point() +
  geom_errorbar(aes(xmin = CI.lower,
                    xmax = CI.upper),
                width = 0, 
                position=position_dodge(width=1)) +

  scale_x_continuous(breaks = c(0, .2, 0.4,0.6, 0.8, 1)) +
  scale_y_continuous(expand = c(0,1)) +
  labs(x = 'Heterozygosity at\nshared G-locus SNPs',
       y = '', color = 'Inferred\nG-locus\ngenotype', 
       title = '') +
  theme_classic() + 
  theme(
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text = element_text(size = 8),
        axis.title.y = element_text(vjust = 6),
        axis.title.x = element_text(size = 9, vjust = -1),
        plot.margin= margin(t=0,r=4,b=1,l=4,unit = 'cm'),
        
        plot.title = element_text(face = 'italic', size = 14),
        strip.text = element_text(face = 'italic'),
        legend.text = element_text(face = 'italic'),
        legend.title = element_text(size = 8)) +
  ggforce::facet_col(vars(species),
                     space="free",
                     scale="free_y")





