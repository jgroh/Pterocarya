library(data.table)
library(ggplot2)
library(viridis)
library(cowplot)

# ----- read blast read alignments ------

# stenoptera
P.ste <- fread("~/workspace/Pterocarya/01_Gloc/10kb_border/P.stenoptera/P.ste_megablast.csv")
P.ste[, Dxy := (100-Similarity)/100]
#P.ste[, Query_start := Query_start + 4354657]
#P.ste[, Query_end := Query_end + 4354657]

#P.ste[, Subject_start := Subject_start + 3869234]
#P.ste[, Subject_end := Subject_end + 3869234]

# macroptera
P.mac <- fread("~/workspace/Pterocarya/01_Gloc/10kb_border/P.macroptera/P.macroptera_megablast.csv")
P.mac[, Dxy := (100-Similarity)/100]


# ----- read gene annotations -----
# P.ste_G_gff <- fread("~/workspace/Pterocarya/denovo_assembly/scaffolded/P.stenoptera_UCD1470/hap2/P.stenoptera_UCD1470_hap2.gff3", fill = T)
# P.ste_G_gff <- P.ste_G_gff[, .(V1, V3, V4, V5, V7, V9)]
# setnames(P.ste_G_gff, c('chr', 'feature', 'start', 'end', 'strand', 'info'))
# 
# P.ste_g_gff <- fread("~/workspace/Pterocarya/denovo_assembly/scaffolded/P.stenoptera_UCD1470/hap1/P.stenoptera_UCD1470_hap1.gff3", fill = T)
# P.ste_g_gff <- P.ste_g_gff[, .(V1, V3, V4, V5, V7, V9)]
# setnames(P.ste_g_gff, c('chr', 'feature', 'start', 'end', 'strand', 'info'))
# 
# # create table for gene annotations
# P.ste_genes_anno <- cbind(
#   P.ste_G_gff[chr == 'Chr11' & 
#                 start >= P.ste[, min(Query_start)] & 
#                 end <= P.ste[, max(Query_end)] &
#                 feature == 'gene', .(xstart = start, xend = end)],
#   P.ste_g_gff[chr == 'Chr11' & 
#                 start >= P.ste[, min(Subject_start)] & 
#                 end <= P.ste[, max(Subject_end)] &
#                 feature == 'gene', .(ystart = start, yend = end)]
# )

# macroptera annotations (use independent braker annotations)
# P.mac_G_gff <- fread("~/workspace/Pterocarya/denovo_assembly/scaffolded/P.macroptera_SBG1992.306B/hap1/P.macroptera_SBG1992.306B_hap1_braker.gff", fill = T)
# P.mac_G_gff <- P.mac_G_gff[, .(V1, V3, V4, V5, V7, V9)]
# setnames(P.mac_G_gff, c('chr', 'feature', 'start', 'end', 'strand', 'info'))
# 
# P.mac_g_gff <- fread("~/workspace/Pterocarya/denovo_assembly/scaffolded/P.macroptera_SBG1992.306B/hap2/P.macroptera_SBG1992.306B_hap2_braker.gff", fill = T)
# P.mac_g_gff <- P.mac_g_gff[, .(V1, V3, V4, V5, V7, V9)]
# setnames(P.mac_g_gff, c('chr', 'feature', 'start', 'end', 'strand', 'info'))
# 
# # create table for gene annotations
# P.mac_genes_anno <- cbind(
#   P.mac_G_gff[chr == 'Chr11' & 
#                 start >= P.mac[, min(Query_start)] & 
#                 end <= P.mac[, max(Query_end)] &
#                 feature == 'gene', .(xstart = start, xend = end)],
#   P.mac_g_gff[chr == 'Chr11' & 
#                 start >= P.mac[, min(Subject_start)] & 
#                 end <= P.mac[, max(Subject_end)] &
#                 feature == 'gene', .(ystart = start, yend = end)]
# )


# coords in Pmac_hap1
# LHS gene: g35298.1, Chr11:4496681-4499939
# FAF-like1: g35299.1, Chr11:4504089-4505027

P.mac_genes_anno <- data.table(xstart = 4504089, xend = 4505027, 
                               ystart = 3882194, yend = 3883144)

# ----- read TE annotations -----
# TE proportion
# Pmac_dom_Chr11_TEprp <- fread("~/workspace/Pterocarya/denovo_assembly/scaffolded/P.macroptera_SBG1992.306B/hap1/Chr11_TEprp_1kb.txt")
# Pmac_rec_Chr11_TEprp <- fread("~/workspace/Pterocarya/denovo_assembly/scaffolded/P.macroptera_SBG1992.306B/hap2/Chr11_TEprp_1kb.txt")
# 
# cnames <- unlist(strsplit("Chr,Source,Type,Start,End,Score,Strand,Phase,Attributes", split = ','))
# P.ste_G_TElines <- readLines("~/workspace/Pterocarya/denovo_assembly/scaffolded/P.stenoptera_UCD1470/hap2/P.stenoptera_UCD1470_hap2.fasta.mod.EDTA.TEanno.gff3")
# P.ste_G_TElines_fltd <- P.ste_G_TElines[!grepl( "#", P.ste_G_TElines)]
# P.ste_G_TEs <- fread(text = paste(P.ste_G_TElines_fltd, collapse = "\n"), col.names = cnames)
# 
# P.ste_g_TElines <- readLines("~/workspace/Pterocarya/denovo_assembly/scaffolded/P.stenoptera_UCD1470/hap1/P.stenoptera_UCD1470_hap1.fasta.mod.EDTA.TEanno.gff3")
# P.ste_g_TElines_fltd <- P.ste_g_TElines[!grepl( "#", P.ste_g_TElines)]
# P.ste_g_TEs <- fread(text = paste(P.ste_g_TElines_fltd, collapse = "\n"), col.names = cnames)
# 
# P.mac_G_TElines <- readLines("~/workspace/Pterocarya/denovo_assembly/scaffolded/P.macroptera_SBG1992.306B/hap1/P.macroptera_SBG1992.306B_hap1.fasta.mod.EDTA.TEanno.gff3")
# P.mac_G_TElines_fltd <- P.mac_G_TElines[!grepl( "#", P.mac_G_TElines)]
# P.mac_G_TEs <- fread(text = paste(P.mac_G_TElines_fltd, collapse = "\n"), col.names = cnames)
# 
# P.mac_g_TElines <- readLines("~/workspace/Pterocarya/denovo_assembly/scaffolded/P.macroptera_SBG1992.306B/hap2/P.macroptera_SBG1992.306B_hap2.fasta.mod.EDTA.TEanno.gff3")
# P.mac_g_TElines_fltd <- P.mac_g_TElines[!grepl( "#", P.mac_g_TElines)]
# P.mac_g_TEs <- fread(text = paste(P.mac_g_TElines_fltd, collapse = "\n"), col.names = cnames)





# ------ Plot stenoptera -----

pste_mn_dom <- P.ste[, min(Query_start)]
pste_mx_dom <- P.ste[, max(Query_start)]
pste_mn_rec <- P.ste[, min(Subject_start)]
pste_mx_rec <- P.ste[, max(Subject_start)]

# SNPs in Pste based on ascertained fixed sites
pste_Gloc_mn_dom <- 4363775 #4364658 based on GWAS  most significant p val
pste_Gloc_mx_dom <- 4429295 # 4428472 based on GWAS most significant p val
pste_Gloc_mn_rec <- 3878403 # 3879235 based on GWAS most significant p val
pste_Gloc_mx_rec <- 3955586 # 3955295 based on GWAS most significant p val

P.ste_plt <- ggplot(P.ste) + 

  # annotate region of association
  annotate("rect", xmin = pste_Gloc_mn_dom,
           xmax = pste_Gloc_mx_dom, 
           ymin=pste_Gloc_mn_rec, 
           ymax = pste_Gloc_mx_rec, 
           #alpha = .3,
           #fill = '#bfbaa1') +
           fill = 'gray80', alpha = 0.3) +

  # annotate TEs at sides 
  # geom_rect(data = P.ste_G_TEs[Chr == 'Chr11' & Start > P.ste[, min(Query_start)] & End < P.ste[, max(Query_start)] ],
  #           aes(xmin = Start, xmax = End, ymin = P.ste[, min(Subject_start)-3000],
  #               ymax = P.ste[, min(Subject_start)]),
  #               #ymax = Inf),
  #           fill = 'darkblue', alpha = 0.5) +
  # geom_rect(data = P.ste_g_TEs[Chr == 'Chr11' & Start > P.ste[, min(Subject_start)] & End < P.ste[, max(Subject_start)]],
  #           aes(ymin = Start, ymax = End, xmin = P.ste[, min(Query_start)-2000], 
  #               xmax = P.ste[, min(Query_start)]), 
  #           #xmax = Inf),
  #           fill = 'darkblue', alpha = 0.5) +
  
  # annotated genes at sides
  #geom_rect(data = P.ste_G_gff[chr == 'Chr11' & feature == 'exon' & start > pste_mn_dom & end < pste_mx_dom], aes(xmin = start, xmax = end, ymin=pste_mn_rec-1e3, ymax = pste_mn_rec+1e3, fill = 'maroon'))   +
  

 

  # add alignment segments
  geom_segment(aes(x = Query_start, xend = Query_end, 
                   y = Subject_start, yend = Subject_end), #color = Dxy),
               linewidth = 0.5) +
  # annotate genes predicted from BRAKER3
  #geom_rect(data = P.ste_genes_anno, aes(xmin=xstart, xmax=xend,
  #         ymin=ystart,ymax=yend), fill = 'green') +
  # geom_rect(data = P.ste_genes_anno, aes(xmin=xstart, xmax=xend,
  #                                        ymin=ystart,ymax=yend), fill = 'turquoise') +
  #geom_rect(data = P.ste_genes_anno, aes(xmin=-Inf, xmax=Inf,
  #                                       ymin=ystart,ymax=yend), fill = 'green') 
  #scale_color_viridis(option = 'A', 
  #                    limits = c( min(0, P.ste[, min(Dxy)], P.mac[, min(Dxy)]),
  #                                max(P.ste[, max(Dxy)], P.mac[, max(Dxy)])))  + 

  
  # annotate FAFL using manually defined UTR coords
  annotate("rect", xmin = 4364001,
           xmax = 4368049, 
           ymin=3878625, 
           ymax = 3882202,
           fill = 'salmon', alpha = 0.6) +
# customize axes
  scale_x_continuous(breaks = seq(4.36e6, 4.44e6, length.out = 5), labels = sprintf("%.2f", seq(4.36, 4.44, length.out = 5) ), expand = c(0,0)) + 
  scale_y_continuous(breaks = seq(3.88e6, 3.96e6, length.out = 5), labels = sprintf("%.2f", seq(3.88, 3.96, length.out = 5) ), expand = c(0,0)) +
  labs(x = expression(italic("G") * " haplotype (Mb)"), y = expression(italic("g") * " haplotype (Mb)"), title = 'P. stenoptera', color = expression(D[XY])) +
  
  # customize theme
  theme_classic() + 
  theme(aspect.ratio = 1.14,
        #legend.position = c(0.95, 0.4),
        legend.position = 'none',
        
        legend.margin = margin(t=0, r=0,b=0,l=0, unit = 'pt'),
        #legend.position = c(1, 0.5),
        plot.margin = unit(c(0.5,0.5,0.5,0.5), 'lines'),
        axis.text = element_text(size=8),
        axis.title = element_text(size = 10),
        #axis.title.y = element_text(size=14, vjust = 2),
        #axis.title.x = element_text(size = 14, vjust = -2),
        #text = element_text(family = "", face = 'bold'), 
        plot.title = element_text(face = 'italic', size = 10),
        #plot.margin = unit(c(.1, 0, 0, 0), "cm")
        #panel.background = element_rect(fill = NA),
        #panel.ontop = TRUE
  ) 
  
  # check positions of misc. features

  # end of Grepeat-0 + 
  #geom_vline(aes(xintercept = 4371311))+
  #geom_hline(aes(yintercept = 3885177))

  # last repeat
  # annotate("rect", xmin = 4398002,
  #          xmax = 4399194, 
  #          ymin=3883305, 
  #          ymax = 3884499 , 
  #          #alpha = .3,
  #          #fill = '#bfbaa1') +
  #          fill = 'green') 

P.ste_plt
  
# what is large TE?
# P.ste_G_TEs[Chr == 'Chr11' & Start < 4.42e6 & End > 4.42e6]
# P.ste_g_TEs[Chr == 'Chr11' & Start < 3.94e6 & End > 3.94e6]
# P.ste_g_TEs[Chr == 'Chr11' & Start > 3.92e6 & End < 3.95e6 & End - Start > 1000]



# ----- Plot macroptera -------

Pmac_G_TE_prp_plt <- ggplot(Pmac_dom_Chr11_TEprp[V2 > P.mac[,min(Query_start)] & V3 < P.mac[, max(Query_end)]], 
       aes(x = (V2+V3)/2)) +
  scale_x_continuous(breaks = seq(4.5e6, 4.58e6, length.out = 5), labels = sprintf("%.2f", seq(4.5, 4.58, length.out = 5) ), expand = c(0,0)) + 
  scale_y_continuous(breaks = c(0,1), expand = c(0,0)) + 
  
  geom_ribbon(aes(ymin = 0, ymax = V4)) + 
  labs(x = '') +
  theme_classic() + 
  theme(aspect.ratio = 0.1)

Pmac_g_TE_prp_plt <-  ggplot(Pmac_rec_Chr11_TEprp[V2 > P.mac[,min(Subject_start)] & V3 < P.mac[, max(Subject_end)]], 
       aes(x = (V2+V3)/2)) +
  scale_x_continuous(breaks = seq(3.88e6, 3.96e6, length.out = 5),
                     labels = sprintf("%.2f", seq(3.88, 3.96, length.out = 5) ), expand = c(0,0)) +
  scale_y_continuous(breaks = c(0,1), expand = c(0,0)) + 
  
  geom_ribbon(aes(ymin = 0, ymax = V4)) + 
  labs(x = '') +
  theme_classic() + 
  theme(aspect.ratio = 0.1)
  
plot_grid(Pmac_G_TE_prp_plt, Pmac_g_TE_prp_plt, ncol = 1)


P.mac_plt <- ggplot(P.mac) + 
  
  # customize axes
  scale_x_continuous(breaks = seq(4.5e6, 4.58e6, length.out = 5), labels = sprintf("%.2f", seq(4.5, 4.58, length.out = 5) ), expand = c(0,0)) + 
  scale_y_continuous(breaks = seq(3.88e6, 3.96e6, length.out = 5),
                     labels = sprintf("%.2f", seq(3.88, 3.96, length.out = 5) ), expand = c(0,0)) +
  labs(x = expression(italic("G") * " haplotype (Mb)"), y = expression(italic("g") * " haplotype (Mb)"), title = 'P. macroptera', color = expression(D[XY])) +
  
  # annotate TEs at sides 
  geom_rect(data = P.mac_G_TEs[Chr == 'Chr11' & Start > P.mac[, min(Query_start)] & End < P.mac[, max(Query_start)] & Type == 'LTR_retrotransposon'],
            aes(xmin = Start, xmax = End, ymin = P.mac[, min(Subject_start)-3000],
                ymax = P.mac[, min(Subject_start)]),
            #ymax = Inf),
            fill = 'darkblue', alpha = 0.5) +
  # geom_rect(data = P.ste_g_TEs[Chr == 'Chr11' & Start > P.ste[, min(Subject_start)] & End < P.ste[, max(Subject_start)]],
  #           aes(ymin = Start, ymax = End, xmin = P.ste[, min(Query_start)-2000], 
  #               xmax = P.ste[, min(Query_start)]), 
  #           #xmax = Inf),
  #           fill = 'darkblue', alpha = 0.5) +
  
  
  
  
  # customize theme
  theme_classic() + 
  theme(aspect.ratio = 1.14,
        legend.position = c(0.9, 0.4),
        legend.margin = margin(t=0, r=0,b=0,l=, unit = 'pt'),
        #legend.position = c(1, 0.5),
        plot.margin = unit(c(0.5,0.5,0.5,0.5), 'lines'),
        axis.text = element_text(size=8),
        axis.title = element_text(size = 10),
        #axis.title.y = element_text(size=14, vjust = 2),
        #axis.title.x = element_text(size = 14, vjust = -2),
        #text = element_text(family = "", face = 'bold'), 
        plot.title = element_text(face = 'italic', size = 10),
        #plot.margin = unit(c(.1, 0, 0, 0), "cm")
        #panel.background = element_rect(fill = NA),
        #panel.ontop = TRUE
  ) + 
  #geom_vline(aes(xintercept = 4553373 )) +
  
  # annotate FAFL using manually defined UTR coords
  # note that these were not determined from RNAseq, but from identifying the 5' UTR splice motif and adding same # of bp
  # upstream, and similar for 3' UTR
  
  annotate("rect", xmin = 4501100,
           xmax = 4505227, 
           ymin=3879871,
           ymax = 3883364,
           fill = 'salmon', alpha = 0.6) +

  # add alignment segments
  geom_segment(aes(x = Query_start, xend = Query_end, 
                   y = Subject_start, yend = Subject_end), #, color = Dxy),
               linewidth = 0.5)  
  # # annotate TEs at sides 
  # geom_rect(data = P.mac_G_TEs[Chr == 'Chr11' & Start > P.mac[, min(Query_start)] & End < P.mac[, max(Query_start)]],
  #           aes(xmin = Start, xmax = End, ymin = P.mac[, min(Subject_start)-3000],
  #               ymax = P.mac[, min(Subject_start)]),
  #               #ymax = Inf),
  #           fill = 'darkblue', alpha = 0.5) +
  # geom_rect(data = P.mac_g_TEs[Chr == 'Chr11' & Start > P.mac[, min(Subject_start)] & End < P.mac[, max(Subject_start)]],
  #           aes(ymin = Start, ymax = End, xmin = P.mac[, min(Query_start)-2000],
  #               xmax = P.mac[, min(Query_start)]),
  #           #xmax = Inf),
  #           fill = 'darkblue', alpha = 0.5) 



  #geom_vline(aes(xintercept = 4501100)) + geom_vline(aes(xintercept = 4571100)) + 
  #geom_hline(aes(yintercept = 3879871)) + geom_hline(aes(yintercept = 3943500))
  #scale_color_viridis(option = 'A', 
  #                    limits = c( min(0, P.ste[, min(Dxy)], P.mac[, min(Dxy)]),
   #                               max(P.ste[, max(Dxy)], P.mac[, max(Dxy)])))  

P.mac_plt


# ----- multipanel of whole region -----

P.ste_P.mac_plt <- plot_grid(P.ste_plt, P.mac_plt, ncol = 1)
P.ste_P.mac_plt
#ggsave(P.ste_P.mac_plt, file = '~/workspace/Pterocarya/01_Figures/dotplot.pdf')





# ------- P. stenoptera closeup ------


P.ste.repeat.coords <- data.table(
  x = c(4369829,4371311,4372970,4375123,4376799,4378980, 
        4380087,4382384,4384925,4385815,4388015,4389646,
        4392388,4395199, 4398005),
  y = rep(3.8830e6, 15), 
  label = 0:14)



P.ste_Grepeat_aln_blastn <- fread("~/workspace/Pterocarya/01_Gloc/Grepeats/G-repeats-dotplot/Pste_Grepeat_blastn-short.csv")

P.ste_Grepeat_aln_blastn[, range(abs(Query_start-Query_end) +1)]



# P.ste_repeats_closeup <- ggplot(P.ste[Query_start > 4.36e6 & 
#                Query_start < 4.405e6 & 
#                Subject_start > 3.877e6 &
#                Subject_end < 3.89e6]) + 

P.ste_repeats_closeup <- ggplot(P.ste_Grepeat_aln_blastn) +
  # annotate TEs at sides 
  # geom_rect(data = P.ste_G_TEs[Chr == 'Chr11' & Start > P.ste[, min(Query_start)] & End < P.ste[, max(Query_start)] & Type == 'helitron'],
  #           aes(xmin = Start, xmax = End, ymin = P.ste_Grepeat_aln_blastn[, min(Subject_start)-1000],
  #               ymax = P.ste_Grepeat_aln_blastn[, min(Subject_start)]),
  #           #ymax = Inf),
  #           fill = 'darkblue', alpha = 0.5) +
  
  # annotate FAFL using manually defined UTR coords
 

  
  # just intron
  annotate("rect", xmin = 4364075,
           xmax = 4366579, 
           ymin=3878720	, 
           ymax = 3880667,
           fill = 'gray', alpha = 0.6) +
  
  # just 5' UTRs
  # 1
  annotate("rect", xmin = 4364001,
           xmax = 4364074,
           ymin=3878625, 
           ymax = 3878719,
           fill = 'turquoise', alpha = 1) +
  # 2
  annotate("rect", xmin = 4366580,
           xmax = 4366880,
           ymin=3880668,
           ymax = 3880924,
           fill = 'turquoise', alpha = 1) +
  

  # just 3' UTR
  annotate("rect", xmin = 4367829,
           xmax = 4368049,
           ymin=3881879,
           ymax = 3882202,
           fill = 'turquoise', alpha = 1) +
  
  
  # just CDS
  annotate("rect", xmin = 4366881,
           xmax = 4367828, 
           ymin=3880925, 
           ymax = 3881878,
           fill = 'maroon', alpha = 0.6) +
  
  
  
  # add alignment segments
  geom_segment(aes(x = Query_start, xend = Query_end, 
                   y = Subject_start, yend = Subject_end),
               linewidth = 0.5) + 
  # scale_color_viridis(option = 'A', 
  #                     limits = c( min(0, P.ste[, min(Dxy)], P.mac[, min(Dxy)]),
  #                                 max(P.ste[, max(Dxy)], P.mac[, max(Dxy)])))   +
  # 
  # re add small 5' UTR
  annotate("rect", xmin = 4364001,
           xmax = 4364074,
           ymin=3878625, 
           ymax = 3878719,
           fill = 'turquoise', alpha = 1) +
  # customize axes
  scale_x_continuous(expand = c(0,0), 
                     breaks = seq(4.365e6, 4.4e6, by = 5e3),
                     labels = sprintf("%.3f", seq(4.365e6, 4.4e6, by = 5e3)/1e6)) +
  scale_y_continuous(expand = c(0,0), 
                     breaks = c(3.878e6, 3.88e6, 3.882e6,3.884e6), labels = sprintf("%.3f", c(3.878e6,3.88e6, 3.882e6,3.884e6)/1e6)) +
  labs(x = expression(italic("G") * " haplotype (Mb)"), y = expression(italic("g") * " haplotype (Mb)"), title = 'P. stenoptera', color = expression(D[XY])) +
  
  # customize theme
  theme_classic() + 
  theme(aspect.ratio = .3,
        legend.position = 'none',
        
        #legend.position = c(0.96, 0.26),
        legend.key.size = unit(8, 'pt'),
        #legend.position = c(1, 0.5),
        plot.margin = unit(c(1,1,1,1), 'lines'),
        axis.text = element_text(size=8),
        axis.title = element_text(size=10),
        #text = element_text(family = "", face = 'bold'), 
        plot.title = element_text(face = 'italic', size = 10),
        #plot.margin = unit(c(.1, 0, 0, 0), "cm")
  )   + 
  
  # ----- check repeat coords -----
  # approximate end of UTR
  #geom_vline(aes(xintercept = 4.368e6)) + 
  
  # G-0 start
  #geom_vline(aes(xintercept = 4369829)) + 
  # G-0 end
  #geom_vline(aes(xintercept = 4371311))  + 
  
  # G-1 start
  #geom_vline(aes(xintercept = 4371311)) + 
  # G-1 end
  #geom_vline(aes(xintercept = 4372970))  + 
  
  # G-2 start
  #geom_vline(aes(xintercept = 4372970))  + 
  # G-2 end
  #geom_vline(aes(xintercept = 4375123)) + 

  # G-3 start
  #geom_vline(aes(xintercept = 4375123)) + 
  # G-3 end
  #geom_vline(aes(xintercept = 4376799))  +
  
  # G-4 start
  #geom_vline(aes(xintercept = 4376799)) + 
  # G-4 end
  #geom_vline(aes(xintercept = 4378980))  +
  
  # G-5 start
  #geom_vline(aes(xintercept = 4378980)) + 
  # G-5 end
  #geom_vline(aes(xintercept = 4380087))  +
  
  # G-6 start
  #geom_vline(aes(xintercept = 4380087)) + 
  # G-6 end
  #geom_vline(aes(xintercept = 4382379))  +
  
  # G-7 start
  #geom_vline(aes(xintercept = 4382384)) + 
  # G-7 end
  #geom_vline(aes(xintercept = 4384925))  +
  
  # G-8 start
  #geom_vline(aes(xintercept = 4384925)) + 
  # G-8 end
  #geom_vline(aes(xintercept = 4385579))  +
  
  # G-9 start
  #geom_vline(aes(xintercept = 4385815)) + 
  # G-9 end
  #geom_vline(aes(xintercept = 4387402))  +
  
  # G-10 start
  #geom_vline(aes(xintercept = 4388015)) + 
  # G-10 end
  #geom_vline(aes(xintercept = 4389646))  +
  
  # G-11 start
  #geom_vline(aes(xintercept = 4389646)) + 
  # G-11 end
  #geom_vline(aes(xintercept = 4392388))  +
  
  # G-12 start
  #geom_vline(aes(xintercept = 4392388)) + 
  # G-12 end
  #geom_vline(aes(xintercept = 4395199))  +
  
  # G-13 start
  #geom_vline(aes(xintercept = 4395199)) + 
  # G-13 end
  #geom_vline(aes(xintercept = 4398005))  +
  
  # G-14 start
  #geom_vline(aes(xintercept = 4398005)) + 
  # G-14 end
  #geom_vline(aes(xintercept = 4400096))  +

  # g-0 start
  #geom_hline(aes(yintercept = 3883305)) + 
  # g-0 end
  #geom_hline(aes(yintercept = 3885177)) + 
  
  geom_text(data = P.ste.repeat.coords, 
            aes(x = x, y = y, label = label), size = 3)  +


  # locations of motif that matches between G repeats and g intron
  
  # location within intron
  geom_segment(aes(y = 3878625+1485+18, x = 4.3625e6+1, xend = 4.3625e6, yend = 3878625+1485+18),
               arrow = arrow(length = unit(0.2, "cm")), color = 'red') +
  
  #locations within repeats
  geom_segment(aes(x = 4369829+3350, xend = 4369829+3350, y = 3.877e6+1, yend = 3.877e6),
               arrow = arrow(length = unit(0.2, "cm"), type = 'closed'), color = 'red') +
  
  geom_segment(aes(x = 4369829+7198, xend = 4369829+7198, y = 3.877e6+1, yend = 3.877e6),
               arrow = arrow(length = unit(0.2, "cm"), type = 'closed'), color = 'red') +
  
  geom_segment(aes(x = 4369829+20032, xend = 4369829+20032, y = 3.877e6+1, yend = 3.877e6),
               arrow = arrow(length = unit(0.2, "cm"), type = 'closed'), color = 'red')  + 
  geom_hline(aes(yintercept = 3884726)) + 
  geom_vline(aes(xintercept = 4372538))
  


P.ste_repeats_closeup

# find inverted repeats
P.ste_Grepeat_aln_blastn[Subject_end > Subject_start & 
                           Subject_end > 3.884e6 & 
                           Subject_end < 3.885e6]

P.ste_Grepeat_aln_blastn[abs(Subject_start - 3884726) < 100 & Subject_start > Subject_end]




# ------even closer, find inverted repeats, stenoptera -----
ggplot(P.ste_Grepeat_aln_blastn[Subject_start > 3.883e6 & Query_start > 4370000]) + 
  # add alignment segments
  geom_segment(aes(x = Query_start, xend = Query_end, 
                   y = Subject_start, yend = Subject_end),
               linewidth = 0.5) + 
  scale_x_continuous(expand = c(0,0), 
                     breaks = seq(4.370e6, 4.4e6, by = 5e3),
                     labels = sprintf("%.3f", seq(4.370e6, 4.4e6, by = 5e3)/1e6))  + 
  scale_y_continuous(expand = c(0,0), 
                     breaks = c(3884000, 3885000),
                     labels = c(3.884, 3.885))  + 
  theme_classic() + 
  theme(aspect.ratio = 0.4,
        axis.text.y = element_text(angle = 90, hjust = 0.5),
        plot.margin = margin(t = 20, r = 20, b = 20, l = 20, unit = 'pt')) + 
  # G-0 start
  # geom_vline(aes(xintercept = 4369829)) + 
  # # G-0 end
  # geom_vline(aes(xintercept = 4371311))  
  # 3' UTR end
  #geom_vline(aes(xintercept = 4368049), linetype = 2, color = 'blue')  + 
  
  geom_segment(aes(y = 3884835, yend = 3884835, x = 4370000, xend = 4370000),
               arrow = arrow(length = unit(0.2, "cm")), color = 'red') 
  
  
  # # just 3' UTR
  # annotate("rect", xmin = 4367829,
  #          xmax = 4368049,
  #          ymin=3881879,
  #          ymax = 3882202,
  #          fill = 'turquoise', alpha = 1) 
  

# isolate a single repeat
ggplot(P.ste_Grepeat_aln_blastn[Subject_start > 3.883e6 & Query_start > 4.38e6 & Query_end <4.385e6 & Query_end < 4.3825e6]) + 
  # add alignment segments
  geom_segment(aes(x = Query_start, xend = Query_end, 
                   y = Subject_start, yend = Subject_end),
               linewidth = 0.5)  + 
  # G-14 IR
  # geom_vline(aes(xintercept = 4399591), color = 'red') + 
  # geom_vline(aes(xintercept = 4399961), color = 'red')
  geom_vline(aes(xintercept = 4381668), color = 'red') +
  geom_vline(aes(xintercept =4382025), color = 'red')


# G-14
P.ste_Grepeat_aln_blastn[Subject_start > 3.883e6 & Query_start > 4.398e6 & Subject_start > Subject_end]
  
# G-7
P.ste_Grepeat_aln_blastn[Subject_start > 3.883e6 & Query_start > 4.38e6 & Query_end <4.385e6 & Subject_start > Subject_end]
  

# for getting coordinates
# P.ste[Query_start > 4.2e6 & 
#         Query_start < 4.405e6 & 
#         Subject_start > 3.878e6 &
#         Subject_end < 3.89e6][7]






# ------ P. mac closeup -----
P.mac.repeat.coords <- data.table(
  x = c(4507000, 4509356, 4511744, 4514072, 4516636, 4518595, 4520067, 4522389, 
        4524698, 4527249, 4529552,4531569, 4533817, 4536135, 4537463, 4545603), #, partial repeat off to the side4553373),
  y = rep(3.8840e6, 16), 
  label = 0:15)


P.mac_Grepeat_blastn <- fread("~/workspace/Pterocarya/01_Gloc/Grepeats/G-repeats-dotplot/Pmac_Grepeat_blastn-short.csv")
P.mac_Grepeat_blastn[, range(abs(Query_start-Query_end)+1)]

# macroptera
# P.mac_repeats_closeup <- ggplot(P.mac[Query_start > 4.5e6 & 
#                                         Query_end < 4.55e6 & # 4.56e6 shows another partial repeat
#                                         Subject_start > 3.879e6 &
#                                         Subject_end < 3.888e6]) +


P.mac_repeats_closeup <- ggplot(P.mac_Grepeat_blastn) +
  
  # just intron
  annotate("rect", xmin = 4501198,
           xmax = 4503803,
           ymin=3879971,
           ymax = 3881928,
           fill = 'gray', alpha = 0.6) +
  
  # just 5' UTRs
  # 1
  annotate("rect", xmin = 4501100,
           xmax = 4501197,
           ymin=3879871, 
           ymax = 3879970,
           fill = 'turquoise', alpha = 1) +
  
  # 2
  annotate("rect", xmin = 4503804,
           xmax = 4504088,
           ymin=3881929,
           ymax = 3882193,
           fill = 'turquoise', alpha = 1) +
  
  # just CDS
  annotate("rect", xmin = 4504089,
           xmax = 4505027,
           ymin=3882194,
           ymax = 3883144,
           fill = 'maroon', alpha = 0.6) +
  
  # just 3' UTR
  annotate("rect", xmin = 4505028,
           xmax = 4505227,
           ymin=3883145,
           ymax = 3883364,
           fill = 'turquoise', alpha = 1) +
  
  # add alignment segments
  geom_segment(aes(x = Query_start, xend = Query_end, 
                   y = Subject_start, yend = Subject_end),
               linewidth = 0.5) + 
  #scale_color_viridis(option = 'A', 
  #                    limits = c( min(0, P.ste[, min(Dxy)], P.mac[, min(Dxy)]),
  #                                max(P.ste[, max(Dxy)], P.mac[, max(Dxy)])))  +
  
  # re add small 5' UTR exon
  annotate("rect", xmin = 4501100,
           xmax = 4501197,
           ymin=3879871, 
           ymax = 3879970,
           fill = 'turquoise', alpha = 1) +
  #scale_x_continuous(breaks = seq(min(Query_start))
  # customize axes
  scale_x_continuous(expand = c(0,0), 
                     breaks = seq(4.5e6, 4.555e6, by = 5e3), 
                     labels = sprintf("%.3f", seq(4.5e6, 4.555e6, by = 5e3)/1e6)) +
  scale_y_continuous(expand = c(0,0),
    breaks = seq(3.881e6, 3.887e6, by = 2e3), 
                     labels = sprintf("%.3f", seq(3.881e6, 3.887e6, by = 2e3)/1e6)) +
   labs(x = expression(italic("G") * " haplotype (Mb)"), y = expression(italic("g") * " haplotype (Mb)"), title = 'P. macroptera', color = expression(D[XY])) +
  
  # customize theme
  theme_classic() + 
  theme(aspect.ratio = .3,
        legend.position = 'none',
        #legend.position = c(1, 0.5),
        plot.margin = unit(c(1,1,1,1), 'lines'),
        axis.text = element_text(size=8),
        axis.title = element_text(size=10),
        #axis.title.y = element_text(size=14, vjust = 2),
        #axis.title.x = element_text(size = 14, vjust = -2),
        #text = element_text(family = "", face = 'bold'), 
        plot.title = element_text(face = 'italic', size = 10),
        #plot.margin = unit(c(.1, 0, 0, 0), "cm")
  )   +
  geom_text(data = P.mac.repeat.coords, 
            aes(x = x, y = y, label = label), size = 2.5) +
  
  # add alignment segments
  geom_segment(aes(x = Query_start, xend = Query_end, 
                   y = Subject_start, yend = Subject_end), #, color = Dxy),
               linewidth = 0.5)  +
  
  # # annotate TEs at sides 
  # geom_rect(data = P.mac_G_TEs[Type == 'L1_LINE_retrotransposon' & Chr == 'Chr11' & Start > P.mac_Grepeat_blastn[, min(Query_start)] & End < P.mac_Grepeat_blastn[, max(Query_start)]],
  #           aes(xmin = Start, xmax = End, ymin = P.mac_Grepeat_blastn[, min(Subject_start)],
  #               ymax = P.mac_Grepeat_blastn[, min(Subject_start)]-1000),
  #           #ymax = Inf),
  #           fill = 'darkblue', alpha = 0.5) +
  # geom_rect(data = P.mac_g_TEs[Type == 'L1_LINE_retrotransposon' & Chr == 'Chr11' & Start > P.mac_Grepeat_blastn[, min(Subject_start)] & End < P.mac_Grepeat_blastn[, max(Subject_start)]],
  #           aes(ymin = Start, ymax = End, xmin = P.mac_Grepeat_blastn[, min(Query_start)],
  #               xmax = P.mac_Grepeat_blastn[, min(Query_start)]-1000),
  #           #xmax = Inf),
  #           fill = 'darkblue', alpha = 0.5) +

  # ----- check repeat coords -----
  # g-0 start
  #geom_hline(aes(yintercept =  3884549)) +
  # g-0 end
  #geom_hline(aes(yintercept = 3886645))  

  # G-0 start
  #geom_vline(aes(xintercept = 4507000)) + 
  # G-0 end
  #geom_vline(aes(xintercept = 4509352))  
  
  # G-1 start
  #geom_vline(aes(xintercept = 4509355)) + 
  # G-1 end
  #geom_vline(aes(xintercept = 4511649))  
  
  # G-2 start
  #geom_vline(aes(xintercept = 4511744))  + 
  # G-2 end
  #geom_vline(aes(xintercept = 4514071)) 

  # G-3 start
  #geom_vline(aes(xintercept = 4514071)) + 
  # G-3 end
  #geom_vline(aes(xintercept = 4516635))  
  
  # G-4 start
  #geom_vline(aes(xintercept = 4516635)) + 
  # G-4 end
  #geom_vline(aes(xintercept = 4518469))  
  
  # G-5 start
  #geom_vline(aes(xintercept = 4518595)) + 
  # G-5 end
  #geom_vline(aes(xintercept = 4519793))  
  
  # G-6 start
  #geom_vline(aes(xintercept = 4520067)) + 
  # G-6 end
  #geom_vline(aes(xintercept = 4522389))  
  
  # G-7 start
  #geom_vline(aes(xintercept = 4522386)) + 
  # G-7 end
  #geom_vline(aes(xintercept = 4524697))  
  
  # G-8 start
  #geom_vline(aes(xintercept = 4524694)) + 
  # G-8 end
  #geom_vline(aes(xintercept = 4527248))  
  
  # G-9 start
  #geom_vline(aes(xintercept = 4527248)) + 
  # G-9 end
  #geom_vline(aes(xintercept = 4529551))  
  
  # G-10 start
  #geom_vline(aes(xintercept = 4529551)) + 
  # G-10 end
  #geom_vline(aes(xintercept = 4531568))  
  
  # G-11 start
  #geom_vline(aes(xintercept = 4531568)) + 
  # G-11 end
  #geom_vline(aes(xintercept = 4533816))  
  
  # G-12 start
  #geom_vline(aes(xintercept = 4533816)) + 
  # G-12 end
  #geom_vline(aes(xintercept = 4536134))  
  
  # G-13 start
  #geom_vline(aes(xintercept = 4536134)) + 
  # G-13 end
  #geom_vline(aes(xintercept = 4537462))  
  
  # G-14 start
  #geom_vline(aes(xintercept = 4537462)) + 
  # G-14 end
  #geom_vline(aes(xintercept = 4539285))  
  
  # G-15 start
  #geom_vline(aes(xintercept = 4545603)) + 
  # G-15 end
  #geom_vline(aes(xintercept = 4547853))  
  
  # G-16 start
 # geom_vline(aes(xintercept = 4553373)) + 
  # G-16 end
  #geom_vline(aes(xintercept = 4554323))  

  # add locations of matching motifs between g intron and G repeats


  # location within intron
  geom_segment(aes(y = 3879871+1519, x = 4.5e6+1, xend = 4.5e6, yend = 3879871+1519),
              arrow = arrow(length = unit(0.2, "cm")), color = 'red') +
  geom_segment(aes(x = 4501100 + 2082 , xend = 4501100 + 2082 , y = 3.879e6+1, yend = 3.879e6),
               arrow = arrow(length = unit(0.2, "cm")), color = 'blue') +
  
  #locations within repeats
  geom_segment(aes(x = 4507000+5180, xend = 4507000+5180, y = 3.879e6+1, yend = 3.879e6),
             arrow = arrow(length = unit(0.2, "cm"), type = 'closed'), color = 'red')  + 
  
  geom_segment(aes(x = 4507000+13500, xend = 4507000+13500, y = 3.879e6+1, yend = 3.879e6),
               arrow = arrow(length = unit(0.2, "cm"), type = 'closed'), color = 'red') +
  geom_segment(aes(x = 4507000+483, xend = 4507000+483, y = 3.879e6+1, yend = 3.879e6),
               arrow = arrow(length = unit(0.2, "cm"), type = 'closed'), color = 'red')  +
  geom_segment(aes(x = 4507000+25017, xend = 4507000+25017, y = 3.879e6+1, yend = 3.879e6),
               arrow = arrow(length = unit(0.2, "cm"), type = 'closed'), color = 'red')   +
  geom_segment(aes(x = 4507000+27253, xend = 4507000+27253, y = 3.879e6+1, yend = 3.879e6),
               arrow = arrow(length = unit(0.2, "cm"), type = 'closed'), color = 'red')  +
  
  geom_segment(aes(x = 4507000+39117, xend = 4507000+39117, y = 3.879e6+1, yend = 3.879e6),
               arrow = arrow(length = unit(0.2, "cm"), type = 'closed'), color = 'blue') +
  geom_segment(aes(x = 4507000+15666, xend = 4507000+15666, y = 3.879e6+1, yend = 3.879e6),
               arrow = arrow(length = unit(0.2, "cm"), type = 'closed'), color = 'blue') 



P.mac_repeats_closeup

library(cowplot)
plot_grid(P.ste_repeats_closeup, P.mac_repeats_closeup, ncol = 1)


# ------ even closer in macroptera, look for inverted repeats -----
ggplot(P.mac_Grepeat_blastn[Subject_start > 3.8851e6 & Subject_end < 3.887e6 & Query_start > 4.506e6]) +   # add alignment segments
  geom_segment(aes(x = Query_start, xend = Query_end, 
                   y = Subject_start, yend = Subject_end), #, color = Dxy),
               linewidth = 0.5)  + 
  theme_classic() + 
  theme(aspect.ratio = .4,
        axis.text.y = element_text(angle = 90, hjust = 0.5),
        plot.margin = margin(t = 20, r = 20, b = 20, l = 20, unit = 'pt')) + 
  # add end of 3' utr
  #geom_vline(aes(xintercept = 4505227), linetype = 2, color = 'blue') + 
  
  scale_x_continuous(expand = c(0.0,0.0), 
                     breaks = seq(4.5e6, 4.555e6, by = 5e3), 
                     labels = sprintf("%.3f", seq(4.5e6, 4.555e6, by = 5e3)/1e6))  + 
    scale_y_continuous(expand = c(0.0,0.0), 
                        breaks = c(3.8856e6, 3.8864e6), 
                       labels = sprintf("%.4f", c(3.8856, 3.8864))) + 
  labs(x = '', y = '') +

  geom_segment(aes(y = 3886178, yend = 3886178, x =4507754, xend = 4507754),
             arrow = arrow(length = unit(0.2, "cm")), color = 'red') 
#3886178
  # first inverted repeat (G-0)
  #geom_vline(aes(xintercept = 4508662), linetype = 2) + 
  #@geom_vline(aes(xintercept = 4508936), linetype = 2)
  

# isolate a single inverted repeat
P.mac_Grepeat_blastn[Subject_start > 3.8855e6 & 
                       Subject_end < 3.887e6 & 
                       Query_start < 4.510e6 & Subject_start > Subject_end]


# -------- Pste G-14 self-alignment ------
cr <- fread("~/workspace/Pterocarya/01_Gloc/Grepeats/P.stenoptera/Pste_G14_self_alignment.csv")
cr[, Dxy := (100-Similarity)/100]

repeat.self.aln <- ggplot(cr) + 
  labs(x = "Repeated motif (kb)", 
       y = "Repeated motif (kb)", 
       title = 'P. stenoptera') +
  
  #scale_x_continuous(breaks = seq(4.398e6, 4.4e6, by = 500), labels = sprintf("%.1f", seq(4398, 4400, by = 0.5))) +
  #scale_y_continuous(breaks = seq(4.398e6, 4.4e6, by = 500), labels = sprintf("%.1f", seq(4398, 4400, by = 0.5))) +
  
  theme_classic() + 
  theme(aspect.ratio = 1,
        plot.title = element_text(face = 'italic')) + 
  annotate("rect", xmin = 178,
           xmax = 361, 
           ymin=178, 
           ymax = 361, 
           fill = 'gray80') +
  annotate("rect", xmin = 385,
           xmax = 555, 
           ymin=385, 
           ymax = 555 , 
           fill = 'gray80') +
  geom_segment(aes(x = Query_start, 
                   xend = Query_end, 
                   y = Subject_start, 
                   yend = Subject_end, color = Dxy), linewidth = 0.8) +
  scale_color_viridis(option = 'A', 
                      limits = c( min(0, P.ste[, min(Dxy)], P.mac[, min(Dxy)]),
                                  max(P.ste[, max(Dxy)], P.mac[, max(Dxy)])))  





# ----- final plots ----
Grepeats_plt <- plot_grid(P.ste_repeats_closeup, P.mac_repeats_closeup, ncol = 1, align = 'v')
Grepeats_plt

P.ste_repeats_closeup
P.mac_repeats_closeup
repeat.self.aln

#ggsave(Grepeats_plt, file = "~/workspace/Pterocarya/01_Figures/Main/dotplot_repeats.pdf")








# --------- close up of UTR intron -----
# 
# intron_view_mn <- 4364020 - 1000
# intron_view_mx <- 4367007
# 
# P.ste_utr_intron <- ggplot(P.ste[Query_start > intron_view_mn & 
#                                         Query_start < intron_view_mx #& 
#                                         #Subject_start > 3.878e6 &
#                                         #Subject_end < 3.89e6
#                                       ]) + 
#   
#   # annotated genes
#   #geom_rect(data = P.ste_genes_anno[2], aes(xmin=xstart, xmax=xend,
#     #                                        ymin=ystart,ymax=yend), fill = 'green') +
#   # add alignment segments
#   geom_segment(aes(x = Query_start, xend = Query_end, 
#                    y = Subject_start, yend = Subject_end, color = Dxy),
#                linewidth = 0.8) + 
#   scale_color_viridis(option = 'A', 
#                       limits = c( min(0, P.ste[, min(Dxy)], P.mac[, min(Dxy)]),
#                                   max(P.ste[, max(Dxy)], P.mac[, max(Dxy)])))   +
#   
#   # customize axes
#  # scale_x_continuous(expand = c(0,1000), 
#  #                    breaks = seq(4.365e6, 4.4e6, by = 5e3),
#  #                    labels = sprintf("%.3f", seq(4.365e6, 4.4e6, by = 5e3)/1e6)) +
#  # scale_y_continuous(breaks = c(3.88e6, 3.882e6,3.884e6), labels = sprintf("%.3f", c(3.88e6, 3.882e6,3.884e6)/1e6)) +
#  # labs(x = expression(italic("G") * " haplotype (Mb)"), y = expression(italic("g") * " haplotype (Mb)"), title = 'P. stenoptera', color = expression(D[XY])) +
#   
#   # customize theme
#  # theme_bw() + 
#   theme(aspect.ratio = .2,
#         legend.position = c(0.96, 0.35),
#         legend.key.size = unit(8, 'pt'),
#         #legend.position = c(1, 0.5),
#         #plot.margin = unit(c(1,1,1,1), 'lines'),
#         axis.text = element_text(size=10),
#         #axis.title.y = element_text(size=14, vjust = 2),
#         #axis.title.x = element_text(size = 14, vjust = -2),
#         #text = element_text(family = "", face = 'bold'), 
#         plot.title = element_text(face = 'italic', size = 10),
#         #plot.margin = unit(c(.1, 0, 0, 0), "cm")
#   )   + 
#   geom_vline(aes(xintercept = 4365148)) + 
#   geom_vline(aes(xintercept = 4365700)) + 
#   # annotate TEs
#   geom_rect(data = P.ste_G_TEs[Chr == 'Chr11' & Start > intron_view_mn & End < intron_view_mx],
#             aes(xmin = Start, xmax = End, ymin = -Inf, 
#                 ymax = Inf), 
#             #ymax = Inf),
#             fill = 'darkblue', alpha = 0.5) 
# 
# 
# P.ste_utr_intron
# 
# 
# 
