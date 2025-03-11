library(data.table) 
library(ggplot2)
library(cowplot)

# ----- metadata -----
sra <- fread("~/workspace/Pterocarya/Cyclocarya_paliurus/wgs_samples_geno.tsv")
diploids <- sra[ploidy == 'diploid', unique(run)]
tetraploids <- sra[ploidy == 'tetraploid', unique(run)]

# ---- read coverage, 4PA assembly ------
cvg_4PA_1bp <- rbindlist(lapply(list.files("~/workspace/Pterocarya/Cyclocarya_paliurus/results/coverage_tetraPA/",
                                       pattern = '*.txt.gz', full.names = T),
                            function(x){
                              z <- fread(x, select = 2:3, col.names = c("position", "coverage"))
                              z[, sample := gsub(".txt.gz", "", basename(x))]
                              return(z)
                            }))

cvg_4PA_1bp[coverage >= 1, avg := median(coverage), by = sample]
cvg_4PA_1bp[, nrm_cvg := coverage/avg]


# ------ coverage in windows ------
#mn <- cvg_4PA[, min(position)]
#mx <- cvg_4PA[, max(position)]
mn <- cvg_4PA_1bp[, min(position)]
mx <- cvg_4PA_1bp[, max(position)]

mn_faf <- 32070679-100e3
mx_faf <- 32071533+100e3

# 100bp
# cvg_4PA_1bp[, window_100bp := cut(position, breaks = seq(mn, mx+1e2, by = 1e2), labels = seq(mn, mx, by = 1e2), include.lowest =T), by = sample]
# cvg_win_4PA_100bp <- cvg_4PA_1bp[, .(coverage = mean(coverage)), by = .(sample, window_100bp)]
# cvg_win_4PA_100bp[, window := as.numeric(as.character((window_100bp)))]
# 
# cvg_win_4PA_100bp[, win_avg_cvg := .SD[coverage != 0, median(coverage)], by = sample]
# cvg_win_4PA_100bp[, nrm_cvg := coverage/win_avg_cvg]

# 1kb
# cvg_4PA_1bp[, window_1kb := cut(position, breaks = seq(mn, mx+1e3, by = 1e3), labels = seq(mn, mx, by = 1e3), include.lowest =T), by = sample]
# cvg_win_4PA_1kb <- cvg_4PA_1bp[, .(coverage = mean(coverage)), by = .(sample, window_1kb)]
# cvg_win_4PA_1kb[, window := as.numeric(as.character((window_1kb)))]
# 
# cvg_win_4PA_1kb[, win_avg_cvg := .SD[coverage != 0, median(coverage)], by = sample]
# cvg_win_4PA_1kb[, nrm_cvg := coverage/win_avg_cvg]

# 5kb
cvg_4PA_1bp[, window_5kb := cut(position, breaks = seq(mn, mx+5e3, by = 5e3), labels = seq(mn, mx, by = 5e3), include.lowest =T), by = sample]

cvg_win_4PA_5kb <- cvg_4PA_1bp[, .(coverage = mean(coverage)), by = .(sample, window_5kb)]
cvg_win_4PA_5kb[, window := as.numeric(as.character((window_5kb)))]
cvg_win_4PA_5kb[, win_avg_cvg := .SD[coverage != 0, mean(coverage)], by = sample]
cvg_win_4PA_5kb[, nrm_cvg := coverage/win_avg_cvg]

# 10kb
# cvg_4PA_1bp[, window_10kb := cut(position, breaks = seq(mn, mx+10e3, by = 10e3), labels = seq(mn, mx, by = 10e3), include.lowest =T), by = sample]
# cvg_win_4PA_10kb <- cvg_4PA_1bp[, .(coverage = mean(coverage)), by = .(sample, window_10kb)]
# cvg_win_4PA_10kb[, window := as.numeric(as.character((window_10kb)))]
# cvg_win_4PA_10kb[, win_avg_cvg := .SD[coverage != 0, median(coverage)], by = sample]
# cvg_win_4PA_10kb[, nrm_cvg := coverage/win_avg_cvg]


# cvg_win_4PA_100bp <- merge(sra[, .(sample=run, geno, ploidy)], cvg_win_4PA_100bp, all = T, by = 'sample')
# cvg_win_4PA_1kb <- merge(sra[, .(sample=run, geno, ploidy)], cvg_win_4PA_1kb, all = T, by = 'sample')
cvg_win_4PA_5kb <- merge(sra[, .(sample=run, geno, ploidy)], cvg_win_4PA_5kb, all = T, by = 'sample')
# cvg_win_4PA_10kb <- merge(sra[, .(sample=run, geno, ploidy)], cvg_win_4PA_10kb, all = T, by = 'sample')

# ------ read blast -----
# ------ read blast alignment -----
blast <- fread("~/workspace/Pterocarya/01_Gloc/Cpal_haplotypes/4PA.13C_vs_Cpal.hap1.csv")

# ----- plots -----
pacbio <- cvg_4PA_1bp[!grepl("SRR", sample), unique(sample)]

# plot just pacbio samples
ggplot(cvg_win_4PA_5kb[sample %in% pacbio], aes(x = window, y = nrm_cvg)) + 
  geom_line(aes(group = sample, color = sample)) + 
  theme_classic() + 
  theme(aspect.ratio = 0.3) +
  geom_vline(aes(xintercept = 32070679)) 

# # pacbio just over FAFL
# ggplot(cvg_win_4PA_100bp[sample %in% pacbio & window >= mn_faf & window <= mx_faf], aes(x = window, y = nrm_cvg)) + 
#   geom_line(aes(group = sample, color = sample)) + 
#   theme_classic() + 
#   theme(aspect.ratio = 0.3) +
#   geom_vline(aes(xintercept = 32070679)) 


# plot just tetraploids
cvg_win_4PA_5kb[sample == 'CRR309100', geno := 'tetraploid PA']
cvg_win_4PA_5kb[sample == 'CRR309100', ploidy := 'tetraploid']

cvg_win_4PA_5kb[sample == 'CRR309098', geno := 'diploid PA']
cvg_win_4PA_5kb[sample == 'CRR309099', geno := 'diploid PG']
cvg_win_4PA_5kb[sample %in% c('CRR309099','CRR309098'), ploidy := 'diploid']


#cvg_win_4PA_10kb[sample == 'CRR309100', geno := '4PA']
#cvg_win_4PA_10kb[sample == 'CRR309098', geno := '2PA']
#cvg_win_4PA_10kb[sample == 'CRR309099', geno := '2PA']


# ----- main coverage plot 13D -----
cvg_win_4PA_5kb[sample == 'SRR23378879', geno := 'tetraploid PA']
cvg_win_4PA_5kb[sample == 'SRR23378891', geno := 'diploid PA']
cvg_win_4PA_5kb[sample == 'CpalSBG', geno := 'diploid PG']

cvg_win_4PA_5kb[, unique(geno)]
cvg_win_4PA_5kb[is.na(geno), ]

cvg_win_4PA_5kb[, geno := factor(geno, levels = c('G1G1', 'G1G2', 'G2G2', 'tetraploid PA', 'diploid PG', 'diploid PA'))]
# renormalize for plot window

cvg_plt_data <- cvg_win_4PA_5kb[window >= 31900000 -10e3 & window <= 32130000 + 10e3]
cvg_plt_data[, nrm_cvg2 := nrm_cvg/mean(nrm_cvg)]

mn2 <- cvg_plt_data[, min(window)]
mx2 <- cvg_plt_data[, max(window)]
cvg_plt_data[, geno := as.factor(geno)]

sra[study %in% c("Qu2023", "Groh"), run]
cvg_plt <- ggplot(cvg_plt_data[! sample %in% c('SRR23378879','SRR23378890','SRR23378891') & 
                                 !sample %in% sra[!study %in% c("Qu2023", "Groh"), run] & 
                                 geno %in% c("diploid PA", "diploid PG")], 
                  aes(x = window, y = nrm_cvg2, color = geno)) +
  
  geom_line(aes(group = sample)) +
  geom_line(data = cvg_plt_data[geno == 'G1G2' ], aes(group = sample, color = 'G1G2'), alpha = 0.6, color = "#DDCC77") + 
  geom_line(data = cvg_plt_data[geno == 'G1G1' ], aes(group = sample, color = 'G1G1'), alpha = 0.6, color = "darkmagenta") + 
  geom_line(data = cvg_plt_data[sample == 'CRR309100' ], aes(group = sample, color = 'CRR309100'), alpha = 1, color = "red", linewidth = 1) + 
  geom_line(data = cvg_plt_data[sample == 'CRR309098' ], aes(group = sample, color = 'CRR309098'), alpha = 1, color = "lightsalmon", linewidth = 1) + 
  
  geom_line(data = cvg_plt_data[sample %in% c("CpalSBG", 'CRR309099') ], aes(group = sample), alpha = 1, color = "cyan", linewidth = 1) + 
  

  scale_x_reverse(
    breaks = rev(seq(31.9e6, 32.1e6, by = 5e4)), 
                  labels = sprintf('%.2f', rev(seq(31.9, 32.1, by = .05))),
                  expand = c(0,0)) +
  scale_color_manual(
  values = c(
    "G1G2" = "#DDCC77",
    "G1G1" = "darkmagenta",
    "G2G2" = NULL,
    "tetraploid PA" = "red",
    "diploid PA" = "lightsalmon",
    "diploid PG" = "turquoise"
  ),
  labels = c(
    "G1G2" = "Tetraploid G1/G1/G1/G2",
    "G1G1" = "Tetraploid G1/G1/G1/G1",
    "G2G2" = '',
    "tetraploid PA" = "Tetraploid PA (G1/G1/G1/G2)",
    "diploid PA" = "Diploid PA (G1/G2)",
    "diploid PG" = "Diploid PG (G1/G1)"
  )) + 

  
  scale_y_continuous( expand = c(0,0), limits =c(-0.7,10)) +
  labs(color = '', x = '', y = 'Normalized depth') +
  theme_classic() + 
  theme(aspect.ratio = 0.3,
        legend.position = 'none'
        #legend.text = element_text(),
        #legend.position = c(0.9,0.9)
  ) +
  # 5' UTR exon 1
  annotate('rect', xmin = 32073938, xmax = 32074033, ymin = -0.5, ymax = -0.1) + 
  # boundaries of exon 2
  annotate('rect', xmin = 32070489, xmax = 32071916, ymin = -0.5, ymax = -0.1) +
  # CDS
  annotate('rect', xmin = 32070679, xmax = 32071648, ymin = -0.5, ymax = -0.1) + 
  annotate('segment', x = 32070489, xend = 32074033, y = -0.3, yend = -0.3)

  # 32073938-32074033 is the 1st exon of 5' UTR
  #32071916 is the start of the 2nd exon (5' UTR)
  # coordinates of CDS 32070679-32071648
  #approximate end of 3' UTR (from blast) 32070489


cvg_plt
cvg_plt_data[, unique(geno)]

# ----- plot dotplot --------
blast[, min(Subject_end)]
dotplot <-  ggplot(blast[Query_end < 32.1e6 & Subject_end < 26.925e6]) + 
  geom_segment(aes(x = Query_start, xend = Query_end, 
                   y = Subject_start, yend = Subject_end), linewidth = 0.6) + 
  scale_x_reverse(
    breaks = rev(seq(31.9e6, 32.1e6, by = 5e4)), 
    labels = sprintf('%.2f', rev(seq(31.9, 32.1, by = .05))),
    expand = c(0,0)) +
  
  scale_y_reverse(expand=c(0,0), 
                  breaks = c(26.9e6, 26.88e6, 26.86e6),
                  labels = sprintf('%.2f', c(26.9, 26.88, 26.86)))+
                  #, limits = c(26853226,26.91e6))+
  labs(x = 'Chr13 G2 haplotype (Mb)', y = 'Chr13 G1 haplotype (Mb)') +
  theme_classic() + 
  theme(aspect.ratio = 1.1, 
        plot.margin = margin(l = 30, r = 30, t = 30, b= 30, unit = 'pt'),
        axis.text = element_text(size = 10)
        ) + 
annotate('rect', ymin = 26898850, ymax = 26902615, 
           xmin = 32070489, xmax = 32074033, alpha = 0.5, fill = 'salmon')

dotplot
plot_grid(cvg_plt, dotplot, ncol = 1, align = 'v')
dotplot

# ------ dotplot finer scale, extract coordinates ------
ggplot(blast[Query_start > 31.95e6 & Query_end < 32.085e6 & Subject_start > 26.875e6 & Subject_end < 26.91e6]) + 
  geom_segment(aes(x = Query_start, xend = Query_end, 
                   y = Subject_start, yend = Subject_end)) + 
  scale_x_reverse(expand = c(0,0)) +
  scale_y_reverse() +
  #scale_y_reverse(expand = c(0,0), limits = c(26.875e6,26.905e6)) + 
  
  annotate('rect', ymin = 26898850, ymax = 26902615, 
           xmin = 32070489, xmax = 32074033, alpha = 0.5, fill = 'salmon') + 
  
  theme_bw()
  # boundaries of repeat region
  #geom_vline(aes(xintercept =32066202)) + 
  #geom_vline(aes(xintercept = 31968996))
  
  # G1-0 
  #geom_hline(aes(yintercept =   26896960)) +
  #geom_hline(aes(yintercept =   26890608)) +
  
  # G2-0
  #geom_vline(aes(xintercept = 32.0686e6))+
  #geom_vline(aes(xintercept = 32066202)) + 
  
  # G2-1
  #geom_vline(aes(xintercept = 32056412)) +
  #geom_vline(aes(xintercept = 32033270))
  
  # G2-2
  #geom_vline(aes(xintercept = 32012706)) +
  #geom_vline(aes(xintercept = 32031869))
  
  # G-3
  #geom_vline(aes(xintercept = 32011304))  +
 # geom_vline(aes(xintercept = 31968200))  
  




# get coordinates of entire repeat region
blast[Query_start > 32.065e6, min(Query_start)]
blast[Query_start < 31.97e6, max(Query_end)]

# get repeat coordinates of G1-0
# blast[Subject_start > 26.895e6 & Subject_end < 26.9e6 & Query_start < 32.012e6, max(Query_end)]
# #32011304
# # blast[Query_end == 32011304]
# 
# blast[Subject_start > 26.89e6 & Subject_end < 26.9e6 & Query_start < 32.013e6, max(Query_start)]
# blast[Query_start == 32012706]

# get repeat coordinates of G2-0
#blast[Subject_start > 26.893e6 & Subject_end < 26.895e6 & Query_start > 32.065e6, min(Query_start)]

# get repeat coordinates of G2-1
#blast[Subject_start > 26.895e6 & Subject_end < 26.898e6 & Query_start > 32.055e6, min(Query_end)]
#blast[Subject_start > 26.89e6 & Subject_start < 26.893e6 & Query_start > 32.033e6, min(Query_start)]

# get repeat coordinates of G2-2
#blast[Subject_start > 26.895e6 & Subject_start < 26.900e6 & Query_start > 32.025e6, min(Query_end)]


# get repeat coordinates of G2-3
#blast[Subject_start > 26.895e6 & Subject_start < 26.900e6 & Query_start > 32.025e6, min(Query_end)]







