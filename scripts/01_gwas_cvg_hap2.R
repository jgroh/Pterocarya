library(data.table)
library(ggplot2)

# ---- read phenotypes -----
meta <- fread("~/workspace/Pterocarya/Pste_WGS_phenotypes.txt", col.names = c("sample", "pheno"))
meta[pheno==0, phenotype := 'male-first']
meta[pheno==1, phenotype := 'female-first']

meta_RNAseq <- fread("~/workspace/Pterocarya/RNAseq/meta.txt")
#meta_RNAseq[pheno==0, phenotype := 'male-first']
#meta_[pheno==1, phenotype := 'female-first']

# ----- read p vals ------
p2 <- fread("~/workspace/Pterocarya/output/PSTE_UCD1470_HAP2/Pste.assoc.txt")
p2[, N:= seq_len(.N)]
# 
nchr <- p2[, length(unique(chr))]
# 
# # ---- read gff -----
#gff2 <- fread("~/workspace/Pterocarya/denovo_assembly/scaffolded/P.stenoptera_UCD1470/hap2/P.stenoptera_UCD1470_hap2_liftoff.gff", fill = T, skip = 3)
gff2 <- fread("~/workspace/Pterocarya/denovo_assembly/scaffolded/P.stenoptera_UCD1470/hap2/P.stenoptera_UCD1470_hap2.gff3", fill = T)

gff2 <- gff2[, .(V1, V3, V4, V5, V7, V9)]
setnames(gff2, c('chr', 'feature', 'start', 'end', 'strand', 'info'))
# 
# # ----read TE annotation ----
cnames <- unlist(strsplit("Chr,Source,Type,Start,End,Score,Strand,Phase,Attributes", split = ','))
TElines <- readLines("~/workspace/Pterocarya/denovo_assembly/scaffolded/P.stenoptera_UCD1470/hap2/P.stenoptera_UCD1470_hap2.fasta.mod.EDTA.TEanno.gff3")
TElines_fltd <- TElines[!grepl( "#", TElines)]
TEs <- fread(text = paste(TElines_fltd, collapse = "\n"), col.names = cnames)

# ----- read LD -----
#ld <- fread("~/workspace/Pterocarya/out.geno.ld", col.names = c("chr", 'pos1', 'pos2', 'n_indv', 'r2'))

# ----- read TSP SNPs ----
#TSP_SNPs <- fread("~/workspace/Pterocarya/01_Gloc/LHS_of_SV/TSP_SNPs_stenoptera_macroptera.txt", col.names ='pos')
# filtered SNP table from minimap2 alignments on LHS of SV
TSP_tbl1 <- fread("~/workspace/Pterocarya/01_Gloc/LHS_wide/SNP_table_fltd.txt")

# ---- read coverage ------
cvg2 <- rbindlist(lapply(list.files("~/workspace/Pterocarya/results/coverage_UCD1470_HAP2//",
                                            pattern = '*_Gloc.txt.gz', full.names = T),
                                 function(x){
                                   z <- fread(x, select = 2:3, col.names = c("position", "coverage"))
                                   z[, sample := gsub("_Gloc.txt.gz", "", basename(x))]
                                   return(z)
                                 }))

# ----- genome-wide average -----
cvg_nrm2 <- rbindlist(lapply(list.files("~/workspace/Pterocarya/results/coverage_IMNU/",
                                          pattern = '*norm.txt', full.names = T),
                               function(x){
                                 z <- fread(x, col.names = c("avg_cvg"))
                                 z[, sample := gsub("_norm.txt", "", basename(x))]
                                 return(z)
                               }))

# combine
cvg2 <- merge(cvg2, cvg_nrm2, all.x = T)

# combine with phenotypes
cvg2 <- merge(cvg2, meta, all.x = T)
cvg2[, pheno := as.factor(pheno)]



# ----- read IsoSeq coverage -----
# IsoSeq_cvg <- rbindlist(lapply(list.files("~/workspace/Pterocarya/results/coverage_IsoSeq_IMNU/",
#                                    pattern = '*_Gloc.txt.gz', full.names = T),
#                         function(x){
#                           z <- fread(x, select = 2:3, col.names = c("position", "coverage"))
#                           z[, sample := gsub("_Gloc.txt.gz", "", basename(x))]
#                           return(z)
#                         }))
# 
# IsoSeq_cvg[sample == 'SRR25617144', tissue := 'catkin']
# IsoSeq_cvg[sample == 'SRR25617145', tissue := 'female_fl']
# IsoSeq_cvg[sample == 'SRR25617146', tissue := 'bark2']
# IsoSeq_cvg[sample == 'SRR25617147', tissue := 'bark1']
# IsoSeq_cvg[sample == 'SRR25617149', tissue := 'fruit']
# IsoSeq_cvg[sample == 'SRR25617150', tissue := 'leaf']
# IsoSeq_cvg[sample == 'SRR26994970', tissue := 'mixed']



# ----- read RNA-seq coverage -----
RNAseq_cvg <- rbindlist(lapply(list.files("~/workspace/Pterocarya/results/coverage_hap2_RNAseq/", recursive = T,
                                          pattern = '*_Gloc_RNAseq.txt.gz', full.names = T),
                               function(x){
                                 z <- fread(x, select = 2:3, col.names = c("position", "coverage"))
                                 z[, sample := gsub("_Gloc_RNAseq.txt.gz", "", basename(x))]
                                 return(z)
                               }))


RNAseq_cvg <- merge(meta_RNAseq[, .(sample = ID, Type, Tissue, Location, Date, WGS)], RNAseq_cvg, by = 'sample')




# ----- read divergence (pending. currently plotted by running plot Dxy script in parallel and sharing environment -----
# Dxy <- rbindlist(lapply(list.files("~/workspace/Pterocarya/anchorwave-alignments/QRY_vs_Pste_hap2_alignments/",
#                                    recursive = T, full.names = T, pattern = "Dxy_Chr11.txt.gz"),
#                         fread))
# # set window size (Dxy already calculated for this window size)
# window_size <- 1e3




# ----- plot gwas -----
#whole genome
nchr
p2[-log10(p_lrt) == max(-log10(p_lrt), na.rm=T)]

p2[!is.na(p_lrt), min(p_lrt)]
p2_chrs <- p2[!grepl("h2t", chr)]
p2_chrs[, chr := as.numeric(chr)]
setkey(p2_chrs, chr)
p2_chrs[, N:= seq_len(.N)]

# takes a while to render, uncomment to run
# ggplot(p2_chrs[], 
#        aes(x = N, y = -log10(p_lrt))) +
#   geom_point(size = 1, aes(color = as.factor(chr))) +
#   scale_color_manual(values = rep(c('black', 'gray'), 8)) +
#   theme_classic() +
#   labs(y = expression(-log[10](P)),
#        x = '') +
#   theme(
#     aspect.ratio = .3,
#     plot.margin = margin(30,30,30,30, "pt"),
#     axis.ticks.length.x = unit(0, 'cm'),
#     axis.text.x = element_blank(),
#     text = element_text(size = 12),
#     axis.text = element_text(size = 12),
#     #axis.text.x = element_blank(),
#     axis.title.x = element_text(size = 10, vjust = -2),
#     axis.title.y = element_text(size = 12, vjust = 2),
#     legend.position = 'none'
#   )


# ----- chromosome -----
p2[p_lrt == min(p_lrt, na.rm=T), chr]

# ggplot(p2[chr == "1"], aes(x = ps, y = -log10(p_lrt))) + 
#   geom_point(size = 1, color = 'black') + 
#   theme_classic() + 
#   labs(y = expression(-log[10](P)), 
#        x = '') +
#   theme(
#     aspect.ratio = .3,
#     plot.margin = margin(30,30,30,30, "pt"),
#     axis.ticks.length.x = unit(0, 'cm'),
#     text = element_text(size = 12),
#     axis.text = element_text(size = 12),
#     axis.title.x = element_text(size = 10, vjust = -2),
#     axis.title.y = element_text(size = 12, vjust = 2),
#     legend.position = 'none'
#   )


# what's the peak on chr1
p2[p_lrt != min(p_lrt, na.rm=T)][p_lrt == min(p_lrt, na.rm=T)]

# ---- focal peak ----
p2[chr == "1"][p_lrt == min(p_lrt)]

mn <- 4266881
mx <- 4534617

p2[p_lrt == min(p_lrt, na.rm = T)]

utr <- fread("~/workspace/Pterocarya/01_Gloc/genes/FAF/UTR/UTR_coords.txt")

gwas_plt_focal <- ggplot(p2[chr == "11" & ps > mn & ps < mx]) + 
  geom_point(size = 0.4, aes(x = ps, y = -log10(p_lrt))) + 
  theme_classic() + 
  labs(y = expression(-log[10](P)), x = ''
       #x = 'Position (Mb)'
       ) +
  scale_x_continuous(breaks = seq(4.3e6, 4.5e6, length.out = 5), labels = sprintf("%.2f", seq(4.3, 4.5, length.out = 5))) +
  scale_y_continuous(expand = c(0.01,0), limits = c(-2.6,20)) +
  theme(
    aspect.ratio = .2,
    plot.margin = margin(5,5,5,5, "pt"),
    text = element_text(size = 8),
    axis.text = element_text(size = 8),
    plot.title = element_text(face = 'italic'),
    axis.title.x = element_text(size = 10, vjust = -2),
    axis.title.y = element_text(size = 10),
    legend.position = c(0.9, 1.05),
    legend.text = element_text(size = 10),
    legend.key.size = unit(.5, 'cm')
  )  +
  #geom_vline(aes(xintercept = 4364658)) +
  #geom_vline(aes(xintercept = 4428444)) +
  geom_rect(data = gff2[chr == 'Chr11' & feature == 'exon' & start > mn & end < mx], aes(xmin = start, xmax = end, ymin=-2.5, ymax = -.5), fill = 'gray40')   +
  geom_segment(data = gff2[chr == 'Chr11' & feature == 'gene' & start > mn & end < mx], aes(x = start, xend = end, y=-1.5, yend = -1.5), color = 'gray40')  +
  
  annotate('rect', xmin = 4364001,xmax= 4364074, ymin=-2.5, ymax = -.5, fill = 'firebrick1') +
  geom_segment(data = gff2[grepl("g21575", info)], aes(x = 4364001, xend = 4364074, y=-1.5, yend = -1.5), color = 'red')  
  
  #geom_vline(aes(xintercept = 4483150), linetype = 2, color = 'red')+
  #geom_vline(aes(xintercept = 4561150), linetype = 2, color = 'red')  
gwas_plt_focal

# 



# ----- plot IsoSeq from Stevens et al. 2018 -----
# 
# IsoSeq_cvg[, cvg_nrm := coverage/mean(coverage), by = sample]
# IsoSeq_cvg_plt <- ggplot(IsoSeq_cvg[position > mn & position < mx 
#                                    # & tissue != 'mixed'
#                                    # & tissue %in% c('catkin','female_fl')
#                                     ]) +
#   geom_line(linewidth = 0.5, alpha = 0.9,aes(x = position, y = cvg_nrm, group = sample, color = tissue) )  +
# 
#   theme_classic() + 
#   labs(x = '', y =  'Read depth', color = '') +
#   theme(
#     aspect.ratio = .3,
#     plot.margin = margin(30,30,30,30, "pt"),
#     text = element_text(size = 12),
#     axis.text = element_text(size = 10),
#     axis.title.x = element_text(size = 10, vjust = -2),
#     axis.title.y = element_text(size = 12, vjust = 2),
#     legend.position = c(0.6, 0.8)
#   ) +
#   geom_rect(data = gff[chr == 'Chr11' & feature == 'exon' & start > mn & end < mx], aes(xmin = start, xmax = end, ymin=-3, ymax = -1))   +
#   geom_rect(data = gff[chr == 'Chr11' & feature == 'gene' & start > mn & end < mx], aes(xmin = start, xmax = end, ymin=-2.1, ymax = -1.9)) #+
# 

# IsoSeq_cvg_plt
# 
# library(cowplot)
# plot_grid(gwas_plt_IMNU, IsoSeq_cvg_plt, ncol = 1)

# ----- RNAseq plot -----

Gloc_LHS <- 4364658
Gloc_RHS <- 4428444

mn2.1 <- 4364658 - 20e3
mx2.1 <- 4428444 + 20e3

RNAseq_cvg[coverage != 0, cvg_nrm := coverage/mean(coverage), by = sample]
RNAseq_cvg[coverage != 0, sqrtcvg := sqrt(coverage), by = sample]

# these are from GWAS, but note separatedly ascertained fixed sites afterward (marginally longer region)
assoc_sites <- data.table(pos = c(4364658,4366662,4366877,4366886,4366912,4366930,
                 4366933,4366945,4367051,4428444,4428472))
#mac_het <- fread("~/Desktop/P.mac_het_sites.txt")
#fra_het <- fread("~/Desktop/P.fra_het_sites.txt")

ggplot(RNAseq_cvg[position > mn2.1 & position < mx2.1 & Tissue == 'mature male pre-anthesis' #& Type == 'female-first'
                  # & tissue %in% c('catkin','female_fl')
]) +
  facet_grid(Tissue~Type) +
  #geom_point(aes(x = position, y = coverage)) +
  geom_ribbon(aes(ymin = 0, ymax = sqrtcvg, x = position), fill = 'purple3') + 
  
  theme_classic() + 
  labs(x = '', y =  expression(sqrt('Read depth')), color = '') +
  #scale_y_continuous(limits = c(-3.5, 10)) +
  theme(
    aspect.ratio = .3,
    plot.margin = margin(30,30,30,30, "pt"),
    text = element_text(size = 12),
    axis.text = element_text(size = 10),
    axis.title.x = element_text(size = 10, vjust = -2),
    axis.title.y = element_text(size = 12, vjust = 2),
    #legend.position = c(0.6, 0.8)
  ) +
  #geom_rug(data = assoc_sites, aes(x = pos, y = NULL), length = unit(.5, 'cm'))   +

  geom_vline(aes(xintercept = Gloc_LHS), linetype = 2) +
  geom_vline(aes(xintercept = Gloc_RHS), linetype = 2) +
  geom_rect(data = gff2[chr == 'Chr11' & feature == 'exon' & start > mn2.1 & end < mx2.1], aes(xmin = start, xmax = end, ymin=-9, ymax = -2))   +
  geom_segment(data = gff2[chr == 'Chr11' & feature == 'gene' & start > mn2.1 & end < mx2.1], aes(x = start, xend = end, y=-5, yend = -5)) 



# ----- examine candidate gene -----
# view with leftmost assoc site in center
#mn2.2 <- 4366881 - 8000
#mx2.2 <- 4367828 + 3000

# view with candidate gene in center
mn2.2 <- 4366881 - 4000
mx2.2 <- 4367828 + 1000

ggplot(RNAseq_cvg[position > mn2.2 & position < mx2.2 & 
                  #  Tissue == 'male buds pre-dormancy' #& Type == 'male-first'
                  Tissue == 'mature male pre-anthesis'
                  # & tissue %in% c('catkin','female_fl')
]) +
  facet_grid(~sample) +
  #geom_point(aes(x = position, y = coverage)) +
  geom_ribbon(aes(ymin = 0, ymax = coverage, x = position), fill = 'purple3') + 
  
  theme_classic() + 
  labs(x = '', y =  'Read depth', color = '') +
  #scale_y_continuous(limits = c(-10, 70)) +
  theme(
    aspect.ratio = .3,
    plot.margin = margin(30,30,30,30, "pt"),
    text = element_text(size = 12),
    axis.text = element_text(size = 10),
    axis.title.x = element_text(size = 10, vjust = -2),
    axis.title.y = element_text(size = 12, vjust = 2),
    #legend.position = c(0.6, 0.8)
  ) +
  #geom_rug(data = fra_het[pos > mn2.2 & pos < mx2.2], aes(x = pos, y = NULL), color = 'green', length = unit(.6, 'cm'))   +
  
  #geom_rug(data = fra_het[pos > mn2.2 & pos < mx2.2], aes(x = pos, y = NULL), color = 'orange', length = unit(.5, 'cm'))   +
  
  geom_vline(aes(xintercept = Gloc_LHS), linetype = 2) +
 # geom_vline(aes(xintercept = Gloc_RHS), linetype = 2) +
  geom_rect(data = gff2[chr == 'Chr11' & feature == 'exon' & start > mn2.2 & end < mx2.2], aes(xmin = start, xmax = end, ymin=-9, ymax = -2))   +
  geom_segment(data = gff2[chr == 'Chr11' & feature == 'gene' & start > mn2.2 & end < mx2.2], aes(x = start, xend = end, y=-5, yend = -5))  + 
  geom_rug(data = assoc_sites[ pos < mx2.2 & pos > mn2.2], aes(x = pos, y = NULL), length = unit(.5, 'cm'), color = 'green')   
  

gff2[chr == 'Chr11' & feature == 'gene' & start > mn2.2 & end < mx2.2]








# ----- calculate coverage in windows -----


# for figure, use 1000
# for below windows to exclude based on CNV, use 100
winsize <- 1000
cvg2[, window := cut(position, breaks = seq(mn, mx+winsize, by = winsize), labels = seq(mn, mx, by = winsize), include.lowest =T), by = sample]
cvg_win2 <- cvg2[, .(coverage = mean(coverage)), by = .(sample, window, avg_cvg)]
cvg_win2[, window := as.numeric(as.character((window)))]

# normalize
cvg_win2[, nrm_cvg := coverage/avg_cvg]

# combine with phenotypes
cvg_win2 <- merge(cvg_win2, meta, all.x = T)
cvg_win2[, pheno := as.factor(pheno)]

# normalize and combine with phenotypes for non-windowed raw coverage
cvg2[, nrm_cvg := coverage/avg_cvg]



# ----- plot windowed coverage focal peak -----

# which male-first looks like it might have a dominant allele
#cvg_win2[window == 4395881 & phenotype == 'male-first']
# DV_134

#which female first has zero cvg over 1/2 of indel
cvg_win2[window == 4395881 & phenotype == 'female-first']


cvg_plt2 <- ggplot(cvg_win2[window > mn & window < mx  
                 # P.ste sample w/o phenotype
                 #& sample != 'PSTE_W19_1'
               # 2 fraxinifolia samples
               # & !sample %in% c("PFRA_W14_2", "PFRA_W18_2")
               # 2 rhoifolia samples
               #& !sample %in% c("PRHO_SO1_S61", "PRHO_SO6_S64")
               # macroptera sample
               & !sample %in% c("PMAC_SO4_S63")
               & phenotype %in% c('male-first', 'female-first')
               #& phenotype %in% c('female-first')
               
               # exclude sample with putative GG genotype, added separately
               #& sample %in% c("PSTE_W17_1")
               #& sample %in% c("PSTE_W1_5") | phenotype %in% c('male-first')
               #& (sample == 'PSTE_DV_134' | phenotype == 'female-first')
               
               
              # PacBio
              & !sample=='SRR26994977'
               ]) +
  #geom_rect(data = TEs[Chr == 'Chr11' & Start > mn & End < mx],
  #          aes(xmin = Start, xmax = End, ymin = 0, ymax = 2), fill = 'lightgray', alpha = 0.5) +
  scale_x_continuous(breaks = seq(4.3e6, 4.5e6, length.out = 5), labels = sprintf("%.2f", seq(4.3, 4.5, length.out = 5))) +
  #scale_y_continuous(limits = c(0, 3)) +
  geom_line(linewidth = 0.5, alpha = 0.5,aes(x = window, y =nrm_cvg, group = sample, color = phenotype))  +
  
  # specific samples
  #geom_line(data = cvg_win2[window > mn & window < mx & sample == 'PSTE_W17_1'],
  #                         linewidth = 1, alpha = 0.9, aes(x = window, y = nrm_cvg, group = sample, color = sample))  +
  #geom_line(data = cvg_win2[window > mn & window < mx & sample == 'PSTE_W1_5'],
  #          linewidth = 1, alpha = 0.9, aes(x = window, y = nrm_cvg, group = sample, color = sample))  +
  #geom_line(data = cvg_win2[window > mn & window < mx & sample == 'PSTE_DV_134'],
  #          linewidth = 1, alpha = 0.9,  aes(x = window, y = nrm_cvg, group = sample, color = sample))  +
  
  
 # '#DC267F','turquoise2','green2'
  
  #scale_color_manual(values = c('#004488','#DDAA33')) +
  scale_color_manual(values = c('darkmagenta','lightsalmon'),
                     labels = c("protogynous", 'protandrous')) +
  
  # Custom colors for phenotype and specific samples
  #scale_color_manual(
  #  values = c('male-first' = '#DDAA33', 'female-first' = '#004488', 
  #             'PSTE_W17_1' = 'green2', 'PSTE_W1_5' = 'turquoise1', 'PSTE_DV_134' = '#DC267F'),
  #  labels = c('male-first' = 'male-first', 'female-first' = 'female-first', 
  #             'PSTE_W17_1' = '17.1 (female-first)', 'PSTE_W1_5' = '1.5 (female-first)', 'PSTE_DV_134' = '134 (male-first)')
  #) +
  
  scale_y_continuous(limits = c(0, 2)) +
  #geom_vline(aes(xintercept = 4364658), linetype = 2) +
  #geom_vline(aes(xintercept = 4428444), linetype = 2) +
  theme_classic() + 
  labs(x = '', y =  'Normalized read depth', color = '') +
  theme(
    aspect.ratio = .2,
    plot.margin = margin(5,5,5,5, "pt"),
    text = element_text(size = 8),
    axis.text = element_text(size = 8),
    plot.title = element_text(face = 'italic'),
    axis.title.x = element_text(size = 10, vjust = -2),
    axis.title.y = element_text(size = 10),
    legend.position = c(0.9, 1.1),
    legend.text = element_text(size = 8),
    legend.key.size = unit(.4, 'cm')
  )+ 
  guides(color = guide_legend(override.aes = list(linewidth = 1, alpha = 1))) #+
   #geom_rect(data = TEs[Chr == 'Chr11' & Start > mn & End < mx ],
          #  aes(ymin = -Inf, ymax = Inf, xmin = Start, xmax = End), fill = 'lightgray') +
  #geom_line(linewidth = 0.5, alpha = 0.8,aes(x = window, y = nrm_cvg, group = sample, color = phenotype) )  
  
  
  
cvg_plt2
library(cowplot)
plot_grid(gwas_plt_focal, cvg_plt2, dxy_plt, ncol = 1, align = 'v')


# ----- coverage anomalies -----
# boundaries for plotting coverage anomalies showing
#mn <- 4364658 - 10000
#mx <-  4428444 + 10000

# ------ Fig 4 ---------
# main view with candidate gene in center
mn2.4 <- 4358000
mx2.4 <- 4370000 #+ 5000

# read depth from long reads
lrdepth <- rbindlist(lapply(list.files("~/workspace/Pterocarya/01_Gloc/LHS_wide/depth/", pattern = "*.txt.gz",
           full.names = T, recursive = T), function(x){
             z <- fread(x)
             z[, id := gsub(".txt.gz", "", basename(x))]
             setnames(z, c("Chr", "pos", "depth", "id"))
           }))

# read UTR coords
utr_coords <- fread("~/workspace/Pterocarya/01_Gloc/genes/FAF/UTR/UTR_coords.txt",
                    col.names = c("assembly", "gene", "utr", "chr", "start", "end"), fill = T)
# format TSP table
TSP_tbl1
# write out TSP sites, for measuring short read Dxy
# fwrite(TSP_tbl1[pos > 4364000 & 
#            Pste_hap2 == Pmac_hap1 &
#            Pste_hap1 == Pmac_hap2 &
#            Pste_hap1 != Pste_hap2, .(pos)],
#        file = '~/workspace/Pterocarya/TSP_test_sites.txt',
#        row.names = F, col.names = F)

# data subsets for plotting
cvg_candidate_region <- cvg2[position > mn2.4 & position < mx2.4 & phenotype %in% c("male-first", "female-first")]

# normalize within individuals
cvg_candidate_region[, cvg_norm := coverage/avg_cvg]

# average across individuals
cvg_candidate_region[, cvg_nrm := mean(cvg_norm), by = .(position, phenotype)]

RNA_seq_candidate_region <- RNAseq_cvg[position > mn2.4 & position < mx2.4 & Tissue == 'mature male pre-anthesis']

RNA_seq_candidate_region[position < 4.363e6, gene_region := 1]
RNA_seq_candidate_region[position > 4.363e6, gene_region := 2]

RNA_seq_candidate_region[, nrm := coverage/quantile(coverage, 0.99), by = .(Type, gene_region)]

# make plot
lroffset <- 0.2

ggplot() + 

  # customize axes 
  #scale_y_continuous(expand = c(0,0), breaks = c(-1,0,1), labels = c(1,0,1)) +
  scale_y_continuous(expand = c(0,0), breaks = NULL) +
  
  scale_x_continuous(expand = c(0,0), breaks = seq(4.358e6, 4.37e6, by = 2e3),
                     labels = seq(4358, 4370, by = 2)) +
  
  # labels
  labs(x = expression(italic("P. stenoptera ") * 'Chr 11 ' * italic("G") * ' haplotype (kb)'), 
       #y = 'Scaled\nRNAseq depth', title = 'P.stenoptera'
       y = '') +
  
  # gene expression
  #geom_ribbon(data = RNA_seq_candidate_region[Type %in% c("male-first", 'female-first')], 
  #            aes(ymin = 0, ymax = nrm/2, x = position, group = Type, fill = Type), alpha = 0.85) +
  #geom_ribbon(data = RNA_seq_candidate_region[Type %in% c('female-first')], 
  #            aes(ymin = 0, ymax = nrm/2, x = position, fill = 'darkmagenta'), alpha = 0.85) +
  
  # pick colors
  #scale_color_manual(
  #scale_color_manual(
  #  values = c('male-first' = 'darkgoldenrod1', 'female-first' = 'darkmagenta')) +
  #scale_fill_manual(
  #  values = c('male-first' = 'darkgoldenrod1', 'female-first' = 'darkmagenta')) +
  
  
  # add long read depths. 4 colors for site classification. 
  
  # Pste_hap2
  geom_ribbon(data = lrdepth[id == 'Pste_hap2' & pos >= mn2.4 & pos <= mx2.4], 
              aes(x = pos, ymin = -2*lroffset, ymax = depth*lroffset-2*lroffset), fill = 'gray') +
 
   geom_segment(data = TSP_tbl1[Pste_hap2 != Pste_hap1 & pos >= mn2.4 & pos <= mx2.4], 
              aes(x = pos, xend= pos, y = -2*lroffset, yend = 1*lroffset-2*lroffset), color = 'darkmagenta') +
  # matches both (and not macroptera)
  geom_segment(data = TSP_tbl1[(Pste_hap2 == Pste_hap1) &
                                 (Pste_hap2 != Pmac_hap2) & 
                                 (Pste_hap2 != Pmac_hap1) & 
                                 pos >= mn2.4 & pos <= mx2.4], 
               aes(x = pos, xend= pos, y = -2*lroffset, yend = 1*lroffset-2*lroffset), color = 'yellow') +
  
  
  # Pste_hap1
   geom_ribbon(data = lrdepth[id == 'Pste_hap1' & pos >= mn2.4 & pos <= mx2.4], 
              aes(x = pos, ymin = -4*lroffset, ymax = depth*lroffset-4*lroffset), fill = 'gray') +
   
   geom_segment(data = TSP_tbl1[Pste_hap1 != Pste_hap2 & pos >= mn2.4 & pos <= mx2.4], 
               aes(x = pos, xend= pos, y = -4*lroffset, yend = 1*lroffset-4*lroffset), color = 'lightsalmon') +
  # matches both (and not macroptera)
  geom_segment(data = TSP_tbl1[(Pste_hap1 == Pste_hap2) &
                                 (Pste_hap1 != Pmac_hap2) & 
                                 (Pste_hap1 != Pmac_hap1) & 
                                 pos >= mn2.4 & pos <= mx2.4], 
               aes(x = pos, xend= pos, y = -4*lroffset, yend = 1*lroffset-4*lroffset), color = 'yellow') +
  
  # Pmac_hap1
  geom_ribbon(data = lrdepth[id == 'Pmac_hap1' & pos >= mn2.4 & pos <= mx2.4], 
              aes(x = pos, ymin = -6*lroffset, ymax = depth*lroffset-6*lroffset), fill = 'gray') +
  # matches Pste G but not Pste g
  geom_segment(data = TSP_tbl1[(Pmac_hap1 == Pste_hap2) & 
                                 (Pmac_hap1 != Pste_hap1) & 
                                 pos >= mn2.4 & pos <= mx2.4], 
               aes(x = pos, xend= pos, y = -6*lroffset, yend = 1*lroffset-6*lroffset), color = 'darkmagenta') +
  # matches neither 
  geom_segment(data = TSP_tbl1[(Pmac_hap2 == Pmac_hap1) &
                                 (Pmac_hap1 != Pste_hap2) & 
                                 (Pmac_hap1 != Pste_hap1) & 
                                 pos >= mn2.4 & pos <= mx2.4], 
               aes(x = pos, xend= pos, y = -6*lroffset, yend = 1*lroffset-6*lroffset), color = 'cyan') +
  
  # Pmac_hap2
  geom_ribbon(data = lrdepth[id == 'Pmac_hap2' & pos >= mn2.4 & pos <= mx2.4], 
              aes(x = pos, ymin = -8*lroffset, ymax = depth*lroffset-8*lroffset), fill = 'gray') +
  # matches Pste g but not Pste G
  geom_segment(data = TSP_tbl1[(Pmac_hap2 == Pste_hap1) & 
                                 (Pmac_hap2 != Pste_hap2) & 
                                 pos >= mn2.4 & pos <= mx2.4], 
               aes(x = pos, xend= pos, y = -8*lroffset, yend = 1*lroffset-8*lroffset), color = 'lightsalmon') +
  # matches neither 
  geom_segment(data = TSP_tbl1[(Pmac_hap2 == Pmac_hap1) &
                                 (Pmac_hap2 != Pste_hap2) &
                                 (Pmac_hap2 != Pste_hap1) & 
                                 pos >= mn2.4 & pos <= mx2.4], 
               aes(x = pos, xend= pos, y = -8*lroffset, yend = 1*lroffset-8*lroffset), color = 'cyan') +
  
  # first position of leftmost fixed difference in P.stenoptera
  #geom_vline(aes(xintercept = 4363775), linetype = 2) + 
  
  # G-0 start vertical line
  geom_vline(aes(xintercept = 4369829), linetype = 2) + 
  # G-0 end
  #geom_vline(aes(xintercept = 4371311))  +
  
  
  

  # annotated genes
  geom_rect(data = gff2[chr == 'Chr11' & feature == 'exon' & start > mn2.4 & end < mx2.4], aes(xmin = start, xmax = end, ymin=0.1, ymax = 0.6), fill = 'maroon')   +
  
  geom_rect(data = utr_coords[assembly == 'P.stenoptera_UCD1470_hap2'], aes(xmin = start, xmax = end, ymin=0.1, ymax = 0.6), fill = 'maroon', alpha = 0.2)   +
  geom_segment(data = utr_coords[assembly == 'P.stenoptera_UCD1470_hap2' & gene == 'g21574'], 
               aes(x = min(start), xend = max(end), y=0.35, yend = 0.35), arrow = arrow(length = unit(0.05, "npc")))  + 
  geom_segment(data = utr_coords[assembly == 'P.stenoptera_UCD1470_hap2' & gene == 'g21575'], 
               aes(x = min(start), xend = max(end), y=0.35, yend = 0.35), arrow = arrow(length = unit(0.05, 'npc')))  + 
  
  
  # add diagnostic SNPs
  #geom_rug(data = assoc_sites[ pos < mx2.4 & pos > mn2.4], aes(x = pos, y = NULL), length = unit(.2, 'cm'), color = 'blue') +

  # add TSP SNPs
  #geom_rug(data = TSP_SNPs, aes(x = pos, y = NULL), length = unit(.2, 'cm'), color = 'red') +
  
  # add TEs
  #geom_segment(data = TEs[Chr == 'Chr11' & Start >= mn2.4 & End <= mx2.4], 
  #            aes(x = Start, xend = End, y = 0.5, yend = 0.5)) +
  # customize theme
  #geom_hline(aes(yintercept = 0)) +
  theme_classic() + 
  theme(aspect.ratio = 0.2, 
        plot.margin = margin(1,1,1,1, 'line'),
        plot.title = element_text(face = 'italic', size = 12),
        legend.position = 'none',
        axis.line.y = element_blank()) #+
  #annotate("segment", x = mn2.4, xend = mn2.4, y = 0, yend = 1, color = "black", size = 0.8)   # Y-axis line + 


# the two coding SNPs that are TSP
#4367555,4367603





# ----- plot anomalous samples -----
# boundaries for plotting coverage anomalies showing
mn <- 4364658 - 10000
mx <-  4428444 + 10000

cvg_anom2 <- ggplot(cvg_win2[window > mn & window < mx  
                            # P.ste sample w/o phenotype
                            #& sample != 'PSTE_W19_1'
                            # 2 fraxinifolia samples
                            # & !sample %in% c("PFRA_W14_2", "PFRA_W18_2")
                            # 2 rhoifolia samples
                            #& !sample %in% c("PRHO_SO1_S61", "PRHO_SO6_S64")
                            # macroptera sample
                            & !sample %in% c("PMAC_SO4_S63")
                            & phenotype %in% c('male-first', 'female-first')
                            #& phenotype %in% c('female-first')
                            
                            # exclude sample with putative GG genotype, added separately
                            #& sample %in% c("PSTE_W17_1")
                            #& sample %in% c("PSTE_W1_5") | phenotype %in% c('male-first')
                            #& (sample == 'PSTE_DV_134' | phenotype == 'female-first')
                            
                            
                            # PacBio
                            & !sample=='SRR26994977'
]) +
  #geom_rect(data = TEs[Chr == 'Chr11' & Start > mn & End < mx],
  #          aes(xmin = Start, xmax = End, ymin = 0, ymax = 2), fill = 'lightgray', alpha = 0.5) +
  #scale_x_continuous(breaks = seq(4.3e6, 4.5e6, length.out = 5), labels = sprintf("%.2f", seq(4.3, 4.5, length.out = 5))) +
  scale_x_continuous(breaks = seq(4.36e6, 4.44e6, length.out = 5), labels = sprintf("%.2f", seq(4.36, 4.44, length.out = 5))) +
  scale_y_continuous(breaks = c(0,1,2)) +
  
  #scale_y_continuous(limits = c(0, 3)) +
  geom_line(linewidth = 0.5, alpha = 0.5,aes(x = window, y =nrm_cvg, group = sample, color = phenotype))  +
  
  # specific samples
  geom_line(data = cvg_win2[window > mn & window < mx & sample == 'PSTE_W17_1'],
            linewidth = 0.7, alpha = 0.9, aes(x = window, y = nrm_cvg, group = sample, color = sample))  +
  geom_line(data = cvg_win2[window > mn & window < mx & sample == 'PSTE_W1_5'],
            linewidth = 0.7, alpha = 0.9, aes(x = window, y = nrm_cvg, group = sample, color = sample))  +
  #geom_line(data = cvg_win2[window > mn & window < mx & sample == 'PSTE_DV_134'],
  #          linewidth = 1, alpha = 0.9,  aes(x = window, y = nrm_cvg, group = sample, color = sample))  +
  
  
  # '#DC267F','turquoise2','green2'
  
  #scale_color_manual(values = c('#004488','#DDAA33')) +
  
  # Custom colors for phenotype and specific samples
  scale_color_manual(
    values = c('male-first' = '#DDAA33', 'female-first' = '#004488', 
               'PSTE_W17_1' = 'green2', 'PSTE_W1_5' = 'turquoise1', 'PSTE_DV_134' = '#DC267F'),
    labels = c('male-first' = 'PA', 'female-first' = 'PG', 
               'PSTE_W17_1' = 'DPTE_17.1 (PG)', 'PSTE_W1_5' = 'DPTE_1.5 (PG)', 'PSTE_DV_134' = '134 (male-first)')
  ) +
  
  #scale_y_continuous(limits = c(0, 10)) +
  geom_vline(aes(xintercept = 4363775), linetype = 2) +
  geom_vline(aes(xintercept = 4429295), linetype = 2) +
  theme_classic() + 
  labs(x = 'Chr 11 hap2 (Mb)', y =  'Normalized read depth', color = '') +
  theme(
    aspect.ratio = .6,
    plot.margin = margin(30,30,30,30, "pt"),
    text = element_text(size = 12),
    axis.text = element_text(size = 8),
    plot.title = element_text(face = 'italic'),
    axis.title.x = element_text(size = 9, vjust = -2),
    axis.title.y = element_text(size = 9, vjust = 2),
    legend.position = c(.35, 0.9),
    legend.text = element_text(size = 8),
    legend.key.size = unit(.2, 'cm')
  )+ 
  guides(color = guide_legend(override.aes = list(linewidth = 1, alpha = 1))) #+
#geom_rect(data = TEs[Chr == 'Chr11' & Start > mn & End < mx ],
#  aes(ymin = -Inf, ymax = Inf, xmin = Start, xmax = End), fill = 'lightgray') +
#geom_line(linewidth = 0.5, alpha = 0.8,aes(x = window, y = nrm_cvg, group = sample, color = phenotype) )  

cvg_anom2





# ----- plot average by genotype, exclude CNV regions -----
avg_cvg_by_geno <- cvg_win[grepl("PSTE", sample) & 
                             phenotype %in% c("male-first", "female-first") & 
                             sample != "PSTE_DV_134", 
                           .(nrm_cvg = mean(nrm_cvg)),
                           by = .(window, phenotype)]

# without exclude
ggplot(avg_cvg_by_geno, aes(x = window, y = nrm_cvg, color = phenotype)) + 
  geom_line() + 
  scale_y_continuous(limits = c(0,2))

# plot masking regions with CNV`
cvg_diff <- dcast(avg_cvg_by_geno, window~phenotype)[, .(diff = abs(`female-first` - `male-first`)), by = window]

cnv_windows <- cvg_diff[diff > 0.2, window]

low_depth_windows <- dcast(avg_cvg_by_geno, window~phenotype)[`female-first` < 0.4 & `male-first` < 0.4, window]

exclude_windows <- c(cnv_windows, low_depth_windows)
# exclude anything in CNV region
exclude_windows <- c(exclude_windows, seq(4500000, 4550000, by = 100))

# excluding
# fixed sites '~/workspace/Pterocarya/Pste_fixed_sites_IMNU.txt'

# the het sites within CNV regions were exluded for checking heterozygosity in other species
UCD2_het_sites <- fread("~/workspace/Pterocarya/calls/IMNU/UCD2_heterozygous_SNPs.bed", col.names = c("chr", "start", "pos"))
UCD2_het_sites_CNV_excluded <- fread("~/workspace/Pterocarya/calls/IMNU/UCD2_heterozygous_SNPs_CNV_excluded.txt", col.names = c("chr", "pos"))

# confirm indexing was done correctly
length(UCD2_het_sites$pos)
length(UCD2_het_sites_CNV_excluded$pos)
length(intersect(UCD2_het_sites$pos,UCD2_het_sites_CNV_excluded$pos))

ggplot(avg_cvg_by_geno[!window %in% exclude_windows & window > 4400000 & window < 4600000]) + 
  geom_point(aes(x = window, y = nrm_cvg, color = phenotype)) + 
  #scale_y_continuous(limits = c(0,2)) + 
  theme_classic() + 
  theme(aspect.ratio = .2) + 
  geom_point(data = data.table(window = as.numeric(Pste_fixed_sites), nrm_cvg = 0.5),
             aes(x = window, y=nrm_cvg)) + 
  geom_point(data = data.table(window = as.numeric(GG_het_sites), nrm_cvg = 0.5),
             aes(x = window, y=nrm_cvg), color = 'blue') + 
  geom_point(data = data.table(window = as.numeric(DV_134_sites), nrm_cvg = 0.5),
             aes(x = window, y=nrm_cvg), color = 'gray') + 
  geom_point(data = UCD2_het_sites[pos > 4400000 & pos < 4600000], aes(x = pos, y = 0.2), color = 'red') +
  geom_point(data = UCD2_het_sites_CNV_excluded[pos > 4400000 & pos < 4600000], aes(x = pos, y = 0.2), color = 'green2') +
  
  geom_rect(data = gff[chr == 'Chr11' & feature == 'exon' & start > 4400000 & end < 4600000], aes(xmin = start, xmax = end, ymin=-0.4, ymax = 0))   +
  geom_segment(data = gff[chr == 'Chr11' & feature == 'gene' & start > 4400000 & end < 4600000], aes(x = start, xend = end, y=-.2, yend = -.2))
  
Pste_fixed_sites

# write out exclude windows
exclude_windows_out <- data.table(chr = 'Chr11', start = exclude_windows, end = exclude_windows + winsize)
setkey(exclude_windows_out, start)
exclude_windows_out

fwrite(exclude_windows_out, file = '~/workspace/Pterocarya/CNV_regions_to_exclude_IMNU.bed',
       quote = F, col.names = F, row.names = F, sep = "\t")

# ------ plot coverage other peaks -----
chr12_cvg <- ggplot(cvg_win_chr12[window > 1432500 & window < 1437500]) +
  geom_line(linewidth = 0.5, alpha = 0.9,aes(x = window, 
                                             y = nrm_cvg,
                                             # y = coverage,
                                             group = sample, color = phenotype) )  +
  geom_line(data = cvg_win_chr12[window > 1432500 & window < 1437500 & sample == 'PSTE_W17_1'],
            linewidth = 0.5, alpha = 0.9, color = 'green', 
            aes(x = window, 
                y = nrm_cvg,
                #y = coverage, 
                group = sample, color = phenotype))  +
  scale_color_manual( values = c('darkmagenta','lightsalmon')) +
  
  #scale_y_continuous(limits = c(0, 10)) +
  theme_classic() + 
  labs(x = '', y =  'Normalized\nread depth', color = '') +
  theme(
    aspect.ratio = .3,
    plot.margin = margin(30,30,30,30, "pt"),
    text = element_text(size = 12),
    axis.text = element_text(size = 12),
    plot.title = element_text(face = 'italic'),
    axis.title.x = element_text(size = 12, vjust = -2),
    axis.title.y = element_text(size = 12, vjust = 2),
    legend.position = c(0.9, 1.05),
    legend.text = element_text(size = 10),
    legend.key.size = unit(.5, 'cm')
  )

library(cowplot)
plot_grid(chr12_pvals, chr12_cvg, ncol = 1)

chr14_cvg <- ggplot(cvg_win_chr14) +
  geom_line(linewidth = 0.5, alpha = 0.9,aes(x = window, 
                                             y = nrm_cvg,
                                             # y = coverage,
                                             group = sample, color = phenotype) )  +
  geom_line(data = cvg_win_chr14[sample == 'PSTE_W17_1'],
            linewidth = 0.5, alpha = 0.9, color = 'green', 
            aes(x = window, 
                y = nrm_cvg,
                #y = coverage, 
                group = sample, color = phenotype))  +
  scale_color_manual( values = c('darkmagenta','lightsalmon')) +
  
  #scale_y_continuous(limits = c(0, 10)) +
  theme_classic() + 
  labs(x = '', y =  'Normalized\nread depth', color = '') +
  theme(
    aspect.ratio = .3,
    plot.margin = margin(30,30,30,30, "pt"),
    text = element_text(size = 12),
    axis.text = element_text(size = 12),
    plot.title = element_text(face = 'italic'),
    axis.title.x = element_text(size = 12, vjust = -2),
    axis.title.y = element_text(size = 12, vjust = 2),
    legend.position = c(0.9, 1.05),
    legend.text = element_text(size = 10),
    legend.key.size = unit(.5, 'cm')
  )

plot_grid(chr14_pvals, chr14_cvg, ncol = 1)

chr16_cvg <- ggplot(cvg_win_chr16) +
  geom_line(linewidth = 0.5, alpha = 0.9,aes(x = window, 
                                             y = nrm_cvg,
                                             # y = coverage,
                                             group = sample, color = phenotype) )  +
  geom_line(data = cvg_win_chr16[sample == 'PSTE_W17_1'],
            linewidth = 0.5, alpha = 0.9, color = 'green', 
            aes(x = window, 
                y = nrm_cvg,
                #y = coverage, 
                group = sample, color = phenotype))  +
  scale_color_manual( values = c('darkmagenta','lightsalmon')) +
  
  #scale_y_continuous(limits = c(0, 10)) +
  theme_classic() + 
  labs(x = '', y =  'Normalized\nread depth', color = '') +
  theme(
    aspect.ratio = .3,
    plot.margin = margin(30,30,30,30, "pt"),
    text = element_text(size = 12),
    axis.text = element_text(size = 12),
    plot.title = element_text(face = 'italic'),
    axis.title.x = element_text(size = 12, vjust = -2),
    axis.title.y = element_text(size = 12, vjust = 2),
    legend.position = c(0.9, 1.05),
    legend.text = element_text(size = 10),
    legend.key.size = unit(.5, 'cm')
  )

plot_grid(chr16_pvals, chr16_cvg, ncol = 1)



# ---- examine coverage in more detail in most obvious indel ---- 

# whole thing
mn1 <- 4480000 ; mx1 <- 4515000

# part 1
mn1 <- 4485000 + 10000 ; mx1 <- 4490000 + 15000


ggplot(cvg[
  position > mn1 & position < mx1
  # P.ste sample w/o phenotype
  & sample != 'PSTE_W19_1'
  # 2 fraxinifolia samples
  & !sample %in% c("PSTE_W14_2", "PSTE_W18_2")
  # 2 rhoifolia samples
  & !sample %in% c("PRHO_SO1_S61", "PRHO_SO6_S64")
  & phenotype %in% c('male-first', 'female-first')
  # exclude sample with homozygous GG genotype, added separately
  & !sample %in% c("PSTE_W17_1")
]) +
  geom_line(linewidth = 0.5, alpha = 0.9,aes(x = position, y = nrm_cvg, group = sample, color = phenotype) )  +
  geom_line(data = cvg[position > mn1 & position < mx1 & sample == 'PSTE_W17_1'],
            linewidth = 0.5, alpha = 0.9, color = 'green', aes(x = position, y = nrm_cvg, group = sample, color = phenotype))  +
  geom_line(data = cvg[position > mn1 & position < mx1 & sample == 'PSTE_DV_134'],
            linewidth = 0.5, alpha = 0.9, color = 'black', aes(x = position, y = nrm_cvg, group = sample, color = phenotype))  +
  scale_color_manual( values = c('darkmagenta','lightsalmon')) +
  
  #scale_y_continuous(limits = c(0, 2.5)) +
  theme_classic() + 
  labs(x = '', y =  'Normalized\nread depth', color = '') +
  theme(
    aspect.ratio = .3,
    plot.margin = margin(30,30,30,30, "pt"),
    text = element_text(size = 12),
    axis.text = element_text(size = 12),
    plot.title = element_text(face = 'italic'),
    axis.title.x = element_text(size = 12, vjust = -2),
    axis.title.y = element_text(size = 12, vjust = 2),
    legend.position = c(0.9, 1.05),
    legend.text = element_text(size = 10),
    legend.key.size = unit(.5, 'cm')
  ) 



# ----- examine coverage over Pst11G003420 -----
# whole gene
mn1 <- 4483161; mx1 <- 4484089

# indel 1
#mn1 <- 4483403 -10; mx1 <- 4483411 + 10

# indel 2
# mn1 <- 4483680; mx1 <- 4483700


ggplot(cvg[
  position > mn1 & position < mx1
               # P.ste sample w/o phenotype
               & sample != 'PSTE_W19_1'
               # 2 fraxinifolia samples
               & !sample %in% c("PSTE_W14_2", "PSTE_W18_2")
               # 2 rhoifolia samples
               & !sample %in% c("PRHO_SO1_S61", "PRHO_SO6_S64")
               & phenotype %in% c('male-first', 'female-first')
               # exclude sample with homozygous GG genotype, added separately
               & !sample %in% c("PSTE_W17_1")
]) +
  geom_line(linewidth = 0.5, alpha = 0.9,aes(x = position, y = nrm_cvg, group = sample, color = phenotype) )  +
  geom_line(data = cvg[position > mn1 & position < mx1 & sample == 'PSTE_W17_1'],
            linewidth = 0.5, alpha = 0.9, color = 'green', aes(x = position, y = nrm_cvg, group = sample, color = phenotype))  +
  geom_line(data = cvg[position > mn1 & position < mx1 & sample == 'PSTE_DV_134'],
            linewidth = 0.5, alpha = 0.9, color = 'black', aes(x = position, y = nrm_cvg, group = sample, color = phenotype))  +
  scale_color_manual( values = c('darkmagenta','lightsalmon')) +

  scale_y_continuous(limits = c(0, 2.5)) +
  theme_classic() + 
  labs(x = '', y =  'Normalized\nread depth', color = '') +
  theme(
    aspect.ratio = .3,
    plot.margin = margin(30,30,30,30, "pt"),
    text = element_text(size = 12),
    axis.text = element_text(size = 12),
    plot.title = element_text(face = 'italic'),
    axis.title.x = element_text(size = 12, vjust = -2),
    axis.title.y = element_text(size = 12, vjust = 2),
    legend.position = c(0.9, 1.05),
    legend.text = element_text(size = 10),
    legend.key.size = unit(.5, 'cm')
  ) + 
  geom_vline(aes(xintercept = 4483161), linetype = 2, color = 'red')

# # ---- combine GWAS and coverage plots ----
# library(cowplot)
# plot_grid(gwas_plt_BNU, cvg_plt_BNU, ncol = 1)
# 
# 
# 
# # ----- plot LD -----
# library(viridis)
# ggplot(ld[!is.na(r2) & 
#             pos1 > 1750000 & pos2 > 1750000 &
#             pos1 < 1900000 & pos2 < 1900000], aes(x = pos2, y = pos1, fill = r2)) + 
#   geom_tile(width = 400, height = 400) +
#   scale_fill_viridis(option = 'D') + 
#   labs(x = "Position (Mb)", y = "Position (Mb)", fill = expression(r^2)) +
#   theme_classic() + 
#   theme(aspect.ratio = 1, 
#         #legend.position = c(0.1, 0.7), 
#         axis.text = element_text(size = 10)) + 
#   annotate("rect", xmin = 1785000, xmax = 1865000,
#            ymin = 1785000, ymax = 1865000,
#            alpha = .9,fill = 'NA', color = 'red') 



