library(data.table) 
library(ggplot2)
library(cowplot)

# read metadata, samples from Qu et al. 
sra <- fread("~/workspace/Pterocarya/Cyclocarya_paliurus/wgs_samples.tsv")

# PG
# SRR23378890
# PA
# SRR23378891	 

# annotation
gff <- fread("grep -v '^#' ~/workspace/Pterocarya/Cyclocarya_genome_assemblies/GWHBKKX00000000.gff", fill = T)
setnames(gff, c("V1", "V3", "V4", "V5", "V7"), c("Chr", "Type", "Start", "End", "Strand"))

gff_curated <- fread("~/workspace/Pterocarya/Cyclocarya_genome_assemblies/2PA_focal_region_curated_ncbi_chrnames.gff")
setnames(gff_curated, c("V1", "V3", "V4", "V5", "V7"), c("Chr", "Type", "Start", "End", "Strand"))


# ---- read coverage, 2PA assembly ------
cvg_2PA <- rbindlist(lapply(list.files("~/workspace/Pterocarya/Cyclocarya_paliurus/results/coverage_Cyclo2PA",
                                   pattern = '*Pterocarya_locus.txt.gz', full.names = T),
                        function(x){
                          z <- fread(x, select = 2:3, col.names = c("position", "coverage"))
                          z[, sample := gsub("_Pterocarya_locus.txt.gz", "", basename(x))]
                          return(z)
                        }))

# normalization, 2PA assembly
cvg_2PA_nrm <- rbindlist(lapply(list.files("~/workspace/Pterocarya/Cyclocarya_paliurus/results/coverage_Cyclo2PA",
                                       pattern = '*norm.txt', full.names = T),
                            function(x){
                              z <- fread(x, col.names = c("avg_cvg"))
                              z[, sample := gsub("_norm.txt", "", basename(x))]
                              return(z)
                            }))

# combine
cvg_2PA <- merge(cvg_2PA, cvg_2PA_nrm, by = 'sample')


# normalize
cvg_2PA[, nrm_cvg := coverage/avg_cvg]






# ----- pacbio coverage -----
cvg_2PA_pb <- rbindlist(lapply(list.files(
  #"~/workspace/Pterocarya/Cyclocarya_paliurus/results/coverage_Cyclo2PA_pacbio//",
  #pattern = '*Pterocarya_locus.txt.gz', full.names = T),
  "~/workspace/Pterocarya/Cpal_expression/Qu_etal_2023/pacbio/results/coverage_2PA/",
  pattern = '*.txt.gz', full.names = T),
                            function(x){
                              z <- fread(x, select = 2:3, col.names = c("position", "coverage"))
                              z[, sample := gsub(".txt.gz", "", basename(x))]
                              return(z)
                            }))

cvg_2PA_pb_nrm <- rbindlist(lapply(list.files(
  #"~/workspace/Pterocarya/Cyclocarya_paliurus/results/coverage_Cyclo2PA_pacbio/",
  "~/workspace/Pterocarya/Cpal_expression/Qu_etal_2023/pacbio/results/coverage_2PA/",
                                           pattern = '*norm.txt', full.names = T),
                                function(x){
                                  z <- fread(x, col.names = c("avg_cvg"))
                                  z[, sample := gsub("_norm.txt", "", basename(x))]
                                  return(z)
                                }))

cvg_2PA_pb <- merge(cvg_2PA_pb, cvg_2PA_pb_nrm, by = "sample")
cvg_2PA_pb[, cvg_nrm := coverage/avg_cvg]




# ------ read blast alignment -----
blast <- fread("~/workspace/Pterocarya/Cyclocarya_genome_assemblies/dotplot/2PA_vs_Cpal_hap2.csv")


# ----- calculate coverage in windows -----
st_2PA <- min(cvg_2PA$position)
en_2PA <- max(cvg_2PA$position)

winsize <- 1000
cvg_2PA[, window := cut(position, breaks = seq(st_2PA, en_2PA+winsize, by = winsize), labels = seq(st_2PA, en_2PA, by = winsize), include.lowest =T), by = sample]
cvg_win_2PA <- cvg_2PA[, .(coverage = mean(coverage)), by = .(sample, window, avg_cvg)]
cvg_win_2PA[, window := as.numeric(as.character((window)))]

# normalize
cvg_win_2PA[, nrm_cvg := coverage/avg_cvg]


# ----- calculate pacbio coverage in windows -----

# for broader view
mn_pb <- cvg_2PA_pb[, min(position)]
mx_pb <- cvg_2PA_pb[, max(position)]


winsize <- 1000
cvg_2PA_pb[, window := cut(position, breaks = seq(mn_pb, mx_pb+winsize, by = winsize), labels = seq(mn_pb, mx_pb, by = winsize), include.lowest =T), by = sample]
cvg_win_2PA_pb <- cvg_2PA_pb[, .(coverage = mean(coverage)), by = .(sample, window, avg_cvg)]
cvg_win_2PA_pb[, window := as.numeric(as.character((window)))]

# normalize
cvg_win_2PA_pb[, nrm_cvg := coverage/avg_cvg]

# ------ plot pacbio coverage -----
mn_pb_plt <- 25.05e6
mx_pb_plt <- 25.160e6

ggplot(cvg_win_2PA_pb[window > mn_pb_plt & window < mx_pb_plt]) + 
  geom_line(aes(x = window , y = nrm_cvg, group = sample, color = sample)) + 
  theme_classic() + 
  #scale_x_reverse(breaks = rev(seq(25e6, 25.2e6, by = 5e4)),
  #               labels = sprintf('%.2f', rev(seq(25, 25.2, by = 0.05))),
  #               expand = c(0,0))+
  scale_y_continuous(limits = c(-0.3, 5)) +
  labs(x = 'Chr 13 (Mb)', y =  'Normalized\nread depth', 
       color = 'CNV genotype', title = '') +
  
  theme(
    aspect.ratio = .3,
    plot.margin = margin(l = 20, r = 20),
    text = element_text(size = 12),
    axis.text = element_text(size = 10),
    plot.title = element_text(face = 'italic'),
    axis.title.x = element_text(size = 10, vjust = -2),
    axis.title.y = element_text(size = 10, vjust = 2),
    legend.position = c(0.9, 1.05),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 9),
    legend.key.size = unit(.2, 'cm')
  ) + 
  
  geom_segment(data = gff_curated[Start > mn_pb_plt & End < mx_pb_plt], 
               aes(x = Start, xend = End, y = -0.1, yend = -0.1)) +
  
  # forward strand, transcribed regions
  geom_segment(data = gff[Chr == 'GWHBKKX00000005' & Start >= mn_pb_plt & End <= mx_pb_plt & Type == 'gene' & Strand == '+'],
               aes(x = Start, xend = End, y = -0.1, yend =-0.1), 
               #arrow = arrow(length = unit(0.1, "cm")), linewidth = 1
  ) + 
  # forward strand, exons
  geom_rect(data = gff[Chr == 'GWHBKKX00000005' & Start >= mn_pb_plt & End <= mx_pb_plt & Type == 'CDS' & Strand == '+'],
            aes(xmin = End, xmax = Start, ymin = -0.21, ymax = -0.01)) +
  
  # reverse strand, transcribed regions
  geom_segment(data = gff[Chr == 'GWHBKKX00000005' & Start >= mn_pb_plt & End <= mx_pb_plt & Type == 'gene' & Strand == '-'],
               aes(x = End, xend = Start, y = -0.1, yend =-0.1), 
               #arrow = arrow(length = unit(0.1, "cm")), linewidth = 1
  ) + 
  # forward strand, exons
  geom_rect(data = gff[Chr == 'GWHBKKX00000005' & Start >= mn_pb_plt & End <= mx_pb_plt & Type == 'CDS' & Strand == '-'],
            aes(xmin = Start, xmax = End, ymin = -0.21, ymax = -0.01)) + 
  
  # Gloc-FAFlike, annotated copy 
  geom_rect(data = gff[grepl("CpaM1st14242", V9)],
            aes(xmin = Start, xmax = End, ymin = -0.21, ymax = -0.01), fill = 'red') +
  
  # Gloc-FAF-like duplicate (only in this haplotype, but not annotated)
  annotate("rect", xmin = 25140259, xmax = 25141220,
           ymin = -0.21, ymax = -0.01, fill = 'red') 





# ----- assign putative genotypes of short read data  -----

geno11 <- cvg_win_2PA[window == 25.128e6
                      & grepl("SRR", sample)
                      & nrm_cvg < 0.1, sample]
geno12 <- cvg_win_2PA[window == 25.128e6
                      & grepl("SRR", sample)
                      & nrm_cvg > 0.1, sample]
geno22 <- cvg_win_2PA[window == 25.128e6
                      & nrm_cvg > 1, sample]
geno12 <- geno12[!geno12 %in% geno22]


genotypes <- data.table(run = c(geno11, geno12, geno22), 
           geno = c(rep("G1G1", length(geno11)),
                    rep("G1G2", length(geno12)),
                    rep("G2G2", length(geno22))
           ))

# fwrite(x = merge(sra, genotypes, all = T), file = '~/workspace/Pterocarya/Cyclocarya_paliurus/wgs_samples_geno.tsv',
#        quote = F, col.names = T, row.names = F, sep = "\t")

cvg_win_2PA[sample %in% geno11, putative_geno := '00 (Illumina)']
cvg_win_2PA[sample %in% geno12, putative_geno := '01 (Illumina)']
cvg_win_2PA[sample %in% geno22, putative_geno := '11 (Illumina)']
cvg_win_2PA[sample == 'SRR23378890', putative_geno := 'known PG']
cvg_win_2PA[sample == 'SRR23378891', putative_geno := 'known PA']
cvg_win_2PA[sample == 'CpalSBG', putative_geno := 'known PG (CCS)']

cvg_win_2PA_pb[sample == 'CRR309098', putative_geno := 'known PA (P6-C4)']
cvg_win_2PA_pb[sample == 'CRR309099', putative_geno := 'known PG (P6-C4)']
cvg_win_2PA_pb[sample == 'CpalSBG', putative_geno := 'known PG (CCS)']

cvg_2PA[sample %in% geno11, putative_geno := '00']
cvg_2PA[sample %in% geno12, putative_geno := '01']
cvg_2PA[sample %in% geno22, putative_geno := '11']
cvg_2PA[sample == 'CpalSBG', putative_geno := 'known PG']





# This sample has dramatically lower coverage than expected
# given , as the coverage is dramatically lower than lister
#cvg_2PA[sample == '2PA_pacbio', putative_geno := 'known PA']

# ------ Plot coverage -----
mn <- 25000000 
mx <- 25200000 

# try alternative normalization
cvg_win_2PA[window > mn & window < mx, nrm2 := .SD[coverage != 0, mean(coverage)], by = sample]
cvg_win_2PA_pb[window > mn & window < mx, nrm2 := .SD[coverage != 0, mean(coverage)], by = sample]

cvg_win_2PA[, nrm2_cvg := coverage/nrm2, by = sample]
cvg_win_2PA_pb[, nrm2_cvg := coverage/nrm2, by = sample]


cvg_win_2PA_sub <- cvg_win_2PA[ sample %in% c(sra[study == 'Qu2023', run], "CpalSBG")]

plt_data <- rbind(cvg_win_2PA_pb[, .(sample, window, nrm_cvg = nrm2_cvg, putative_geno)],
      cvg_win_2PA_sub[!sample %in% c('SRR23378890',  'SRR23378891'), .(sample, window, nrm_cvg = nrm2_cvg, putative_geno)])


library(pals)

# CpaM1st14242 is the annotated ortholog of Gloc-FAF-like. The other copy is not annotated. Coordinates determined by blast.

cvg_plt <- ggplot(plt_data[window > mn & window < mx]) +
  geom_line(aes(x = window, y = nrm_cvg, group = sample, color = putative_geno)) + 
  geom_line(data = plt_data[window > mn & window < mx & putative_geno == '01 (Illumina)'],
            aes(x = window, y = nrm_cvg, group = sample), 
               color = "#DDCC77", linewidth = 0.5, alpha = 0.9)+
    geom_line(data = plt_data[window > mn & window < mx & putative_geno == '00 (Illumina)'],
              aes(x = window, y = nrm_cvg, group = sample), 
              color = "#332288",linewidth = 0.5, alpha = 0.9)+
    geom_line(data = plt_data[window > mn & window < mx & putative_geno == 'known PA (P6-C4)'],
              aes(x = window, y = nrm_cvg, group = sample), 
              color = "maroon", linewidth = 0.5, alpha = 0.9)+
    geom_line(data = plt_data[window > mn & window < mx & putative_geno == 'known PG (P6-C4)'],
              aes(x = window, y = nrm_cvg, group = sample), 
                  color = "green", linewidth = 0.5, alpha = 0.9)+
    geom_line(data = plt_data[window > mn & window < mx & putative_geno == 'known PG (CCS)'],
              aes(x = window, y = nrm_cvg, group = sample), 
                  color = "cyan", linewidth = 0.5, alpha = 0.9) +
    
  scale_x_reverse(breaks = rev(seq(25e6, 25.2e6, by = 5e4)), 
                  labels = sprintf('%.2f', rev(seq(25, 25.2, by = 0.05))),
                  expand = c(0,0)) +
  scale_y_continuous(limits = c(-0.5, 3), expand = c(0,0))+
  scale_color_manual(values = c("01 (Illumina)" = "#DDCC77",
                                "00 (Illumina)" = "#332288",
                                "known PA (P6-C4)" = "maroon",
                                "known PG (P6-C4)" = 'green',
                                "known PG (CCS)" = 'cyan'),
                     labels = c(
                       "01 (Illumina)" = "01 (Illumina)",
                         "00 (Illumina)" = "00 (Illumina)",
                         "known PA (P6-C4)" = "known PA (P6-C4)",
                         "known PG (P6-C4)" = "known PG (P6-C4)",
                         "known PG (CCS)" = "known PG (CCS)")
                     ) +


  theme_classic() + 
    labs(x = 'Chr 13 (Mb)', y =  'Normalized read depth', 
       color = '', title = '') +

  theme(
    aspect.ratio = .3,
    plot.margin = margin(l = 20, r = 20),
    text = element_text(size = 12),
    axis.text = element_text(size = 10),
    plot.title = element_text(face = 'italic'),
    axis.title.x = element_text(size = 10, vjust = -2),
    axis.title.y = element_text(size = 10, vjust = 2),
    legend.position = c(0.9, 1.05),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 9),
    legend.key.size = unit(.2, 'cm')
  ) + 
  # transcribed regions, non-dup region 
  geom_segment(data = gff[Chr == 'GWHBKKX00000005' & ((Start >= mn & End <= 25.05e6) | (Start >=25.155e6 & End <= mx) ) & Type == 'gene'],
               aes(x = Start, xend = End, y = -0.2, yend =-0.2), 
               ) + 
  # exons
  geom_rect(data = gff[Chr == 'GWHBKKX00000005' & ((Start >= mn & End <= 25.05e6) | (Start >=25.15e6 & End <= mx) ) & Type == 'CDS'],
               aes(xmin = End, xmax = Start, ymin = -0.35, ymax = -0.05)) +
  
  # transcribed regions, dup region 1
  geom_segment(data = gff[Chr == 'GWHBKKX00000005' & (Start >= 25.05e6 & End <= 25.1e6) & Type == 'gene'],
               aes(x = Start, xend = End, y = -0.2, yend =-0.2), color = 'turquoise') + 
  geom_rect(data = gff[Chr == 'GWHBKKX00000005' & (Start >= 25.05e6 & End <= 25.1e6) & Type == 'CDS'],
            aes(xmin = End, xmax = Start, ymin = -0.35, ymax = -0.05), fill = 'turquoise') +
  
  # transcribed regions, dup region 2
  geom_segment(data = gff[Chr == 'GWHBKKX00000005' & (Start >= 25.1e6 & End <= 25.16e6) & Type == 'gene'],
               aes(x = Start, xend = End, y = -0.2, yend =-0.2), color = 'blue') + 
  geom_rect(data = gff[Chr == 'GWHBKKX00000005' & (Start >= 25.1e6 & End <= 25.16e6) & Type == 'CDS'],
            aes(xmin = End, xmax = Start, ymin = -0.35, ymax = -0.05), fill = 'blue') +

  
  # Gloc-FAF-like duplicate (only in this haplotype, but not annotated)
  annotate("rect", xmin = 25140259, xmax = 25141220,
           ymin = -0.3, ymax = -0.05, fill = 'blue')  +
  # contig break
  geom_vline(aes(xintercept = 25.080e6), linetype = 2) + 
  annotate('text', x = 25140500, y = -.45, label = '*')

cvg_plt




# ----- plot dotplot --------
dotplot <-  ggplot(blast) + 
  geom_segment(aes(x = Query_start, xend = Query_end, 
                   y = Subject_start, yend = Subject_end)) + 
  scale_x_reverse(breaks = rev(seq(25e6, 25.2e6, by = 5e4)), 
                  labels = sprintf('%.2f', rev(seq(25, 25.2, by = 0.05))),
                  expand = c(0,0))+
  scale_y_reverse(breaks = rev(seq(26.15e6, 26.25e6, by = 5e4)), 
                  labels = sprintf('%.2f', rev(seq(26.15, 26.25, by = 0.05))),
                  expand = c(0.1,0))+
  labs(x = '2PA Chr 13 (Mb)', y = '2PG SBG.hap2 Chr 13 (Mb)') +

  theme_classic() + 
  theme(aspect.ratio = 0.3, 
        plot.margin = margin(l = 20, r = 20, t = -50)) + 
  annotate('rect', ymin = 26187313,ymax = 26188281, xmin = -Inf, xmax = Inf, alpha = 0.1) + 
  # contig break
  geom_vline(aes(xintercept = 25.080e6), linetype = 2)


#dotplot
# ---- combine -----

plot_grid(cvg_plt, dotplot, ncol = 1, align = 'v')




# -------- fine scale coverage, CpGFAFL1 ------
ggplot(cvg_2PA[position > 25055078 - 100 & position < 25056040 + 100
                                  #& sample != 'SRR23378891'
]) +
  geom_line(aes(x = position, y = nrm_cvg, group = sample, color = putative_geno, linewidth = putative_geno), alpha = 0.9)  +
  scale_x_reverse() +
  # scale_x_reverse(breaks = rev(seq(25e6, 25.2e6, by = 5e4)), 
  #                 labels = sprintf('%.2f', rev(seq(25, 25.2, by = 0.05))),
  #                 expand = c(0,0))+
  scale_y_continuous(limits = c(-0.21, 2))+
  
  scale_color_manual( values = c(tol()[c(1,7)], 'turquoise1', 'maroon', 'darkgreen') ) +
  scale_linewidth_manual(values = c(0.7,0.5,0.4, 0.4, 0.5)) +
  guides(linewidth = 'none') +
  #scale_x_continuous(breaks = rev(seq(25e6, 25.2e6, by = 5e4)), 
  #                   labels = sprintf('%.2f', rev(seq(25, 25.2, by = 0.05)))
  #                   ) +
  # geom_line(data = cvg_win_2PA[window > 25028000 & window < 25148000 & sample == 'SRR23378891'], color = 'green') +
  
  #geom_line(data = cvg_win_2PA[window >mn & window < mx & 
  #                               sample == 'SRR23378890'],
  #                             aes(x = window, y = nrm_cvg, group = sample, color = putative_geno), 
  #          color = 'turquoise1', linewidth = 0.5) +
  theme_classic() + 
  labs(x = 'Chr 13 (Mb)', y =  'Normalized\nread depth', 
       color = 'CNV genotype', title = '') +
  
  theme(
    aspect.ratio = .3,
    plot.margin = margin(l = 20, r = 20),
    text = element_text(size = 12),
    axis.text = element_text(size = 10),
    plot.title = element_text(face = 'italic'),
    axis.title.x = element_text(size = 10, vjust = -2),
    axis.title.y = element_text(size = 10, vjust = 2),
    legend.position = c(0.9, 1.05),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 9),
    legend.key.size = unit(.2, 'cm')
  ) 
  

# -------- fine scale coverage, CpGFAFL2 ------
#annotate("rect", xmin = 25140259, xmax = 25141220,
 
ggplot(cvg_2PA[position > 25140259 & position < 25141220
               & sample != 'SRR23378891'
]) +
  geom_line(aes(x = position, y =nrm_cvg, group = sample, color = putative_geno, linewidth = putative_geno), alpha = 0.9)  +
  scale_x_reverse() +
  # scale_x_reverse(breaks = rev(seq(25e6, 25.2e6, by = 5e4)), 
  #                 labels = sprintf('%.2f', rev(seq(25, 25.2, by = 0.05))),
  #                 expand = c(0,0))+
  #scale_y_continuous(limits = c(-0.21, 2))+
  
  scale_color_manual( values = c(tol()[c(1,7)], 'turquoise1', 'maroon', 'darkgreen') ) +
  scale_linewidth_manual(values = c(0.7,0.5,0.4, 0.4, 0.5)) +
  guides(linewidth = 'none') +
  #scale_x_continuous(breaks = rev(seq(25e6, 25.2e6, by = 5e4)), 
  #                   labels = sprintf('%.2f', rev(seq(25, 25.2, by = 0.05)))
  #                   ) +
  geom_line(data = cvg_2PA[position > 25140259 & position < 25141220 & sample == 'SRR23378891'], 
            aes(x = position , y = nrm_cvg), color = 'green') +
  
  #geom_line(data = cvg_win_2PA[window >mn & window < mx & 
  #                               sample == 'SRR23378890'],
  #                             aes(x = window, y = nrm_cvg, group = sample, color = putative_geno), 
  #          color = 'turquoise1', linewidth = 0.5) +
  theme_classic() + 
  labs(x = 'Chr 13 (Mb)', y =  'Normalized\nread depth', 
       color = 'CNV genotype', title = '') +
  
  theme(
    aspect.ratio = .3,
    plot.margin = margin(l = 20, r = 20),
    text = element_text(size = 12),
    axis.text = element_text(size = 10),
    plot.title = element_text(face = 'italic'),
    axis.title.x = element_text(size = 10, vjust = -2),
    axis.title.y = element_text(size = 10, vjust = 2),
    legend.position = c(0.9, 1.05),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 9),
    legend.key.size = unit(.2, 'cm')
  )  



  
