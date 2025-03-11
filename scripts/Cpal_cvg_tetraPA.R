library(data.table)
library(ggplot)


# ---- read coverage, 2PA assembly ------
cvg_4PA <- rbindlist(lapply(list.files("~/workspace/Pterocarya/Cyclocarya_paliurus/results/coverage_tetraPA/",
                                       pattern = '*.txt.gz', full.names = T),
                            function(x){
                              z <- fread(x, select = 2:3, col.names = c("position", "coverage"))
                              z[, sample := gsub(".txt.gz", "", basename(x))]
                              return(z)
                            }))

# normalization, 2PA assembly
cvg_4PA_nrm <- rbindlist(lapply(list.files("~/workspace/Pterocarya/Cyclocarya_paliurus/results/coverage_Cyclo2PA/",
                                           pattern = '*norm.txt', full.names = T),
                                function(x){
                                  z <- fread(x, col.names = c("avg_cvg"))
                                  z[, sample := gsub("_norm.txt", "", basename(x))]
                                  return(z)
                                }))

# combine
cvg_4PA <- merge(cvg_4PA, cvg_4PA_nrm, by = 'sample')



# ----- calculate coverage in windows -----
st_4PA <- min(cvg_4PA$position)
en_4PA <- max(cvg_4PA$position)

winsize <- 10000
cvg_4PA[, window := cut(position, breaks = seq(st_4PA, en_4PA+winsize, by = winsize), labels = seq(st_4PA, en_4PA, by = winsize), include.lowest =T), by = sample]
cvg_win_4PA <- cvg_4PA[, .(coverage = mean(coverage)), by = .(sample, window, avg_cvg)]
cvg_win_4PA[, window := as.numeric(as.character((window)))]

# normalize
cvg_win_4PA[, nrm_cvg := coverage/avg_cvg]


ggplot(cvg_win_4PA) +
  geom_line(aes(x = window, y = nrm_cvg, group = sample)) +
  # scale_x_reverse(breaks = rev(seq(25e6, 25.2e6, by = 5e4)), 
  #                 labels = sprintf('%.2f', rev(seq(25, 25.2, by = 0.05))),
  #                 expand = c(0,0))+
   scale_y_continuous(limits = c(-0.21, 2))
  
  # scale_color_manual(values = c(tol()[c(1,7)], 'turquoise1')) +
  # scale_linewidth_manual(values = c(0.7,0.5,0.4)) +
  # guides(linewidth = 'none') +
  #scale_x_continuous(breaks = rev(seq(25e6, 25.2e6, by = 5e4)), 
  #                   labels = sprintf('%.2f', rev(seq(25, 25.2, by = 0.05)))
  #                   ) +
  # geom_line(data = cvg_win_2PA[window > 25028000 & window < 25148000 & sample == 'SRR23378891'], color = 'green') +
  
  #geom_line(data = cvg_win_2PA[window >mn & window < mx & 
  #                               sample == 'SRR23378890'],
  #                             aes(x = window, y = nrm_cvg, group = sample, color = putative_geno), 
  #          color = 'turquoise1', linewidth = 0.5) +
  theme_classic() 
  












