library(data.table)
library(ggplot2)

# ---- read putative genotypes based on coverage -----
geno <- fread("~/workspace/Pterocarya/Cyclocarya_paliurus/wgs_samples_geno.tsv")


#  ------ allelic depth, Illumina data -----
# *** this must include ALL sites, even non-variant, for correction to work 

# short read diploids
adsr <- fread("~/workspace/Pterocarya/Cyclocarya_paliurus/results/allelic_depth_Cpal_hap1_GFAFL1_diploids.txt", header = F, sep = ' ')
setnames(adsr, paste0("V", 1:3), c("pos", "ref", "alt"))
adsr <- melt(adsr, id.vars = c("pos", "ref", "alt"))

adsr[, c("sample", "AD") := tstrsplit(value, split = "=")]
adsr[, value := NULL]
adsr[, variable := NULL]

# merge with putative genotypes 
adsr <- merge(adsr, geno[, .(sample = run, geno)], by = 'sample', all = T)
adsr[sample == 'SRR23378891', type := 'PA']
adsr[sample == 'SRR23378890', type := 'PG']

# short read tetraploids
adsr.tetra <- fread("~/workspace/Pterocarya/Cyclocarya_paliurus/results/allelic_depth_Cpal_hap1_GFAFL1_tetraploids.txt", header = F, sep = ' ')
setnames(adsr.tetra, paste0("V", 1:3), c("pos", "ref", "alt"))
adsr.tetra <- melt(adsr.tetra, id.vars = c("pos", "ref", "alt"))

adsr.tetra[, c("sample", "AD") := tstrsplit(value, split = "=")]
adsr.tetra[, value := NULL]
adsr.tetra[, variable := NULL]
adsr.tetra[, geno := NA]
adsr.tetra[sample == 'SRR23378879', type := 'PA']

adsr <- rbind(adsr[, .(sample, pos, ref, alt, AD, geno, type)],
              adsr.tetra[, .(sample, pos, ref, alt, AD, geno, type)])

adsr <- merge(adsr, geno[, .(sample = run, ploidy, study)])


# ----- pacbio allelic depth -----
# This only worked for the HiFi sample. The Qu2023 data is P6-C4 chemistry and error rate is too high,
# so variant calling failed. So will use short read samples from assembly inds for this purpose. 

adpb <- fread("~/workspace/Pterocarya/Cyclocarya_paliurus/results/allelic_depth_Cpal_hap1_GFAFL1_pacbio.txt", header = F, sep = ' ')
setnames(adpb, paste0("V", 1:3), c("pos", "ref", "alt"))
adpb <- melt(adpb, id.vars = c("pos", "ref", "alt"))

adpb[, c("sample", "AD") := tstrsplit(value, split = "=")]
adpb[, value := NULL]
adpb[, variable := NULL]
adpb[, sample := gsub("/group/gmcoopgrp3/jgroh/Pterocarya2/Cyclocarya_paliurus/alignment_files/Cpal_hap1/","",sample)]
adpb[, sample := gsub(".bam","",sample)]
adpb[sample == 'CpalSBG', geno := 'G1G1']
adpb[sample == 'CpalSBG', type := 'PG']
adpb[sample == 'CpalSBG', ploidy := 'diploid']

adpb[sample == 'CRR309098', type := 'PA']
adpb[sample == 'CRR309098', geno := 'G1G2']
adpb[sample == 'CRR309098', ploidy := 'diploid']

adpb[sample == 'CRR309099', type := 'PG']
adpb[sample == 'CRR309099', geno := 'G1G1']
adpb[sample == 'CRR309099', ploidy := 'diploid']

adpb[sample == 'CRR309100', type := 'PA']
adpb[sample == 'CRR309100', geno := 'G1G1G1G2']
adpb[sample == 'CRR309100', ploidy := 'tetraploid']

adpb[grepl("CRR", sample), study := 'Qu2023']
adpb[!grepl("CRR", sample), study := 'Groh2025']

ad <- rbind(adsr[, .(sample, pos, ref, alt, AD,type, geno, ploidy, study)], 
            adpb[, .(sample, pos, ref, alt, AD,type, geno, ploidy, study)])


ad[, c("ref_cnt", "alt_cnt") := tstrsplit(AD, split = ",")]
ad[, ref_cnt := as.numeric(ref_cnt)]
ad[, alt_cnt := as.numeric(alt_cnt)]

# for invariant sites alt_cnt is zero according to above formatting procedure, set alt_cnt to zero here
ad[is.na(alt_cnt), alt_cnt := 0]
ad[, AD := NULL]
ad <- ad[!is.na(ref_cnt) & !is.na(alt_cnt)]
ad[, pos := gsub("CM056974.1:", "", pos)]

ad

# ------ Apply binomial correction -------

# assign unique identifier to alleles, (rather than positions) as there may be some multi-allelic sites.
ad[, allele.id := seq_len(.N), by = sample]

ad[, minor_allele_count := min(ref_cnt, alt_cnt), by = .(sample, allele.id)]

ad[ sample == 'CRR309100' & pos == 26902560, unique(ref_cnt)]

ad[, site_depth_pt1 := unique(ref_cnt), by = .(sample, pos)]
ad[, site_depth_pt2 := sum(alt_cnt), by = .(sample, pos)]
ad[, site_depth := site_depth_pt1+site_depth_pt2]

# ascertain sites. criterion: 0.05 freq among GFAFL copies within individual.
ad[, maf := minor_allele_count/(site_depth)]


# are there some multiallelic sites within individuals? 
# yes, for many individuals (63%), but at a tiny fraction sites within individuals (~1%, calculation below)
# nonetheless in the calculation here these are included as separate segregating sites
# to conform within infinite-sites model
ad[maf >= 0.05, .(Nseg = .N, npos =length(unique(pos))), by = sample]

# check that site depth was calculated correctly, should be sum across alleles at multi-allelic sites 
# without double counting ref allele
ad[maf >= 0.05, .SD[.N > 1], by = .(pos,sample)]


# calculate observed number of segregating sites
S_obs <- ad[maf >= 0.05, .(segregating_sites = .N), by = .(sample)]


#  write function to compute depth correction factor from binomial sampling

SingleSiteCorrectionFactor <- function(site_depth, threshold_maf=0.05, freq){
  
  prob.seg.given.depth <- vector()
  sum_index <- ceiling(threshold_maf*site_depth):floor((1-threshold_maf)*site_depth)
  
  # if read depth is zero or one, there is zero probability of ascertaining conditioning on site segregating
  if(site_depth %in% c(0,1)){
    prob.seg.given.depth <- 0
  } else{
    
    # otherwise loop over possible minor allel counts that satisfy criterion for an allele to count as a segregating site
    for(i in 1:length(sum_index)){
      prob.seg.given.depth[i] <- dbinom(x = sum_index[i], size = site_depth, prob = freq)
    }
  }
  return(sum(prob.seg.given.depth))
}

# check values are correct for test cases
# SingleSiteCorrectionFactor(site_depth = 0, freq = 0.5)
# SingleSiteCorrectionFactor(site_depth = 1, freq = 0.5)
# SingleSiteCorrectionFactor(site_depth = 2, freq = 0.5)
# SingleSiteCorrectionFactor(site_depth = 3, freq = 0.5)
# SingleSiteCorrectionFactor(site_depth = 4, freq = 0.5)

# test out in data.table using sample with lowest average depth
# ad[sample == 'SRR23378852', SingleSiteCorrectionFactor(site_depth, freq = 0.5), by = allele.id]

# compute depth correction factor

ad[ploidy == 'diploid', DepthCorrectionFactor_tetra0.25 := SingleSiteCorrectionFactor(site_depth = site_depth, freq = 0.5), by = .(sample, site_depth, allele.id)]
ad[ploidy == 'tetraploid', DepthCorrectionFactor_tetra0.25 := SingleSiteCorrectionFactor(site_depth = site_depth, freq = 0.25), by = .(sample, site_depth, allele.id)]

ad[ploidy == 'diploid', DepthCorrectionFactor_all0.5 := SingleSiteCorrectionFactor(site_depth = site_depth, freq = 0.5), by = .(sample, site_depth, allele.id)]
ad[ploidy == 'tetraploid', DepthCorrectionFactor_all0.5 := SingleSiteCorrectionFactor(site_depth = site_depth, freq = 0.5), by = .(sample, site_depth, allele.id)]


DepthCorrectionFactors <- ad[, .(DepthCorrectionFactor_tetra0.25 = mean(DepthCorrectionFactor_tetra0.25), 
                                  DepthCorrectionFactor_all0.5 = mean(DepthCorrectionFactor_all0.5)), 
                              by = .(sample, type, geno, ploidy, study)]

seg_final <- merge(S_obs[, .(S_obs = segregating_sites, sample)], DepthCorrectionFactors)

seg_final[, S_hat := S_obs/DepthCorrectionFactor_tetra0.25]

# check correspondence between correction factors in tetraploid
ggplot(DepthCorrectionFactors[],
       aes(x = DepthCorrectionFactor_tetra0.25, 
           y = DepthCorrectionFactor_all0.5, color = ploidy)) + 
  geom_point() + theme_classic() +
  geom_abline()


# ----- read coverage over GFAFL copies -----
G1 <- rbindlist(lapply(list.files("~/workspace/Pterocarya/Cyclocarya_paliurus/results/coverage_Cyclo2PA_GFAFL1/",
                                           pattern = '*.txt.gz', full.names = T),
                                function(x){
                                  z <- fread(x, select = 2:3, col.names = c("position", "coverage"))
                                  z[, sample := gsub(".txt.gz", "", basename(x))]
                                  return(z)
                                }))
G1_cvg <- G1[, .(G1_cvg = mean(coverage)), by = sample]

G2 <- rbindlist(lapply(list.files("~/workspace/Pterocarya/Cyclocarya_paliurus/results/coverage_Cyclo2PA_GFAFL2/",
                                  pattern = '*.txt.gz', full.names = T),
                       function(x){
                         z <- fread(x, select = 2:3, col.names = c("position", "coverage"))
                         z[, sample := gsub(".txt.gz", "", basename(x))]
                         return(z)
                       }))
G2_cvg <- G2[, .(G2_cvg = mean(coverage)), by = sample]

GFAFL_cvg <- merge(G1_cvg, G2_cvg)

# for normalization
cvg_2PA_nrm <- rbindlist(lapply(list.files("~/workspace/Pterocarya/Cyclocarya_paliurus/results/coverage_Cyclo2PA_GFAFL1/norm/",
                                           pattern = '*_norm.txt.gz', full.names = T),
                                function(x){
                                  z <- fread(x, col.names = c("avg_cvg"))
                                  z[, sample := gsub("_norm.txt.gz", "", basename(x))]
                                  return(z)
                                }))

cvg_2PA_nrm


GFAFL_cvg <- merge(GFAFL_cvg, cvg_2PA_nrm, by = 'sample')
GFAFL_cvg[, G1_cvg_nrm := G1_cvg/avg_cvg]
GFAFL_cvg[, G2_cvg_nrm := G2_cvg/avg_cvg]
GFAFL_cvg[, G2_prp := G2_cvg_nrm/(G1_cvg_nrm + G2_cvg_nrm), by = sample]


fin <- merge(seg_final, GFAFL_cvg, all = T)

fin


# final data clean up
# the samples CRR309098 and SRR23378891 should be the same sample, one long read (older pacbio) and the latter short read.
# the pacbio is not suitable for calling segregating sites. On the other hand the SRR23378891 sample
# shows highly variable coverage patterns across the whole genome, suggesting there could be 
# pcr duplication or some other artefact of library prep. indeed, it appears as a bit of an outlier in the read depth ratio here
# compared to all other diploid samples that have non-zero coverage over GFAFL2.
# Consolidate these two data sources for this individual, 
# as a compromise take the segregating sites from the short read but the depth from the long read

fin[sample %in% c('CRR309098','SRR23378891')]
fin[sample == 'CRR309098', S_hat := fin[sample == 'SRR23378891', S_hat]]
fin[sample == 'CRR309098', type := 'PA']
fin[sample == 'CRR309098', ploidy := 'diploid']
fin[sample == 'CRR309098', study := 'Qu2023']


fin <- fin[sample != 'SRR23378891']


# double check that known individuals have phenotype information
fin[is.na(type), type := 'unknown']
fin[type != 'unknown']

fin <- fin[!is.na(ploidy) & !is.na(type)]


ggplot(fin[study %in% c("Qu2023", 'Groh2025')], 
       aes(x = S_hat, y = G2_prp, size = avg_cvg, color = ploidy, shape = type)) + 
  geom_point() + 
  labs(x = 'Segregating sites', y = 'Proportion of reads mapping to G2') + 
  theme_classic() + 
  theme(aspect.ratio = 1)



# ---- last step: correction for copy number (Watterson's theta) ----
# denominator of Watterson's theta for tetraploids is 1.833 (1 + 1/2 + 1/3)

fin[, S_hat2 := S_obs/DepthCorrectionFactor_all0.5]

fin[ploidy == 'diploid', ThetaW := S_hat]
fin[ploidy == 'tetraploid', ThetaW := S_hat/1.833]

fin[ploidy == 'diploid', ThetaW.2 := S_hat2]
fin[ploidy == 'tetraploid', ThetaW.2 := S_hat2/1.833]

fin[type != 'unknown']
fin[, avg_cvg := as.numeric(avg_cvg)]

ggplot(fin[ ],#study %in% c("Qu2023", 'Groh2025')], 
       aes(x = S_hat, y = G2_prp, color = ploidy,  shape = type, size = type, alpha = log10(avg_cvg))) + 
  geom_jitter(width = 0, height = 0,stroke=1.1) + 
  
  scale_shape_manual(values = c(24,25,21)) +
  #scale_shape_manual(values = c(18,17,16)) +
  
  scale_size_manual(values = c(4,4,2), guide = 'none') +
  scale_color_manual(values = c("cyan4", "chocolate3")) +
  
  scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75), limits = c(-0.02, 0.76)) +
  
  labs(x = 'Corrected no. segregating sites', 
       y = 'Relative depth: G2/(G1+G2)', 
       color = 'Ploidy', shape = 'Type', size = NULL,
       alpha = expression(Log[10] ~ '(depth)')) + 
  #guides(fill = guide_legend(override.aes = list(shape = 1))) +
  theme_bw() + 
  theme(aspect.ratio = 1,
        plot.margin = margin(t=30,r=30,b=30,l=30, unit = 'pt')) 

# which sample
#fin[S_hat > 250 & ploidy == 'tetraploid' & G2_prp < 0.125]
#fin[sample == 'SRR23378874']

fin[G2_prp < 0.375 & G2_prp > 0.125 & ploidy == 'tetraploid', geno := 'G1G2']
fin[G2_prp < 0.125 & ploidy == 'tetraploid', geno := 'G1G1']



# ------- Ratio tests ----- 
fin[grepl("SRR", sample) & ploidy == 'tetraploid' & type == 'unknown' & study == 'Qu2023', table(geno)]
fin[grepl("SRR", sample) & ploidy == 'diploid' & type == 'unknown' & study == 'Qu2023', table(geno)]

prop.test(13, 13+22)
prop.test(8, 10)


# fwrite(merge(geno[, .(run, ID, ploidy, phenotype, study)],
#              fin[, .(run = sample, geno)]),
#        file = "~/workspace/Pterocarya/Cyclocarya_paliurus/wgs_samples_geno.tsv",
#        quote = F, col.names = T, row.names = F, sep = "\t")

# fwrite(x = merge(sra, genotypes, all = T), file = '~/workspace/Pterocarya/Cyclocarya_paliurus/wgs_samples_geno.tsv',
#        quote = F, col.names = T, row.names = F, sep = "\t")
# 
# ggplot(fin[ploidy == 'diploid' & geno == 'G1G2'], 
#        aes(x = S_obs, y = S_hat, size = avg_cvg)) + geom_point() + 
#   geom_abline()
# 
# ggplot(fin[ploidy == 'tetraploid' & geno == 'G1G2'], 
#        aes(x = S_obs, y = S_hat, size = avg_cvg)) + geom_point() + 
#   geom_abline()

