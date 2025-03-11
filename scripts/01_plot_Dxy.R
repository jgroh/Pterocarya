library(data.table)
library(ggplot2)

# set window size (Dxy already calculated for this window size)
window_size <- 1e3

# ----- read divergence -----
Dxy <- rbindlist(lapply(list.files("~/workspace/Pterocarya/anchorwave-alignments/QRY_vs_Pste_hap2_alignments/",
           recursive = T, full.names = T, pattern = ".txt.gz"),
       fread))

# set viewing interval

mn <- 4266881
mx <- 4534617

# expected Dxy for completely random alignment
#Dxy[, Dxy_misalign := 1-(frq.A^2 + frq.T^2 + frq.C^2 + frq.G^2), by = window]


# -----read TE annotation ----
cnames <- unlist(strsplit("Chr,Source,Type,Start,End,Score,Strand,Phase,Attributes", split = ','))
P.ste_G_TElines <- readLines("~/workspace/Pterocarya/denovo_assembly/scaffolded/P.stenoptera_UCD1470/hap2/P.stenoptera_UCD1470_hap2.fasta.mod.EDTA.TEanno.gff3")
P.ste_G_TElines_fltd <- P.ste_G_TElines[!grepl( "#", P.ste_G_TElines)]
P.ste_G_TEs <- fread(text = paste(P.ste_G_TElines_fltd, collapse = "\n"), col.names = cnames)

# TEs annotated by repeatModeler + repeatMasker
RM <- fread("~/workspace/Pterocarya/denovo_assembly/scaffolded/P.stenoptera_UCD1470/hap2/P.stenoptera_UCD1470_hap2.fasta.out", fill = T, skip = 3)
setnames(RM, c("V5", "V6", "V7"), c('Chr', 'Start', 'End'))


# ----- read vcf of long read alignments from Pste and macroptera -----
# ---- LHS TSP -----
map1 <- fread("~/workspace/Pterocarya/01_Gloc/LHS_wide/toPsteHap2.map")
ped1 <- fread("~/workspace/Pterocarya/01_Gloc/LHS_wide/toPsteHap2.ped")
pos1 <- map1[, .(pos = V4)]

nuc1 <- ped1[,-c(2:6)]
setnames(nuc1, "V1", "ID")

# each column is duplicated immediately after because it thinks they're diploid
# how many SNPs?
mxn1 <- max( as.numeric( gsub("V", "", names(nuc1)[-1] ) ) )
SNP_nx2.1 <- mxn1 - 7 + 1
nSNP1 <- SNP_nx2.1/2

selectcols1 <- c("ID", paste0("V", seq(7, mxn1, by = 2)))
nuc1 <- nuc1[, ..selectcols1]

setnames(nuc1, paste0("V", seq(7, mxn1, by = 2)), as.character(1:(SNP_nx2.1/2)))

# to long, clean up for plotting
nuc1 <- melt(nuc1, id.vars = "ID", variable.name = "rank", value.name = "base")
nuc1[, pos := pos1[, pos], by = ID]
nuc1[,ID := gsub(".bam", "", basename(ID))]

nuc1[ID == 'P.ste_hap1', ID := 'g Pste']
nuc1[ID == 'P.ste_hap2', ID := 'G Pste']
nuc1[ID == 'P.mac_hap1', ID := 'G Pmac']
nuc1[ID == 'P.mac_hap2', ID := 'g Pmac']

nuc1[, base := factor(base, levels = c("0", "A", "T", 'C', "G"))]

# spread by ind
TSP_tbl1 <- dcast(nuc1, pos~ID, value.var = 'base')
TSP_tbl1

# confirm only biallelic sites
#nuc[, length(unique(base)), by = pos][, unique(V1)] # I think filtered in plink conversion
TSP_sites1 <- TSP_tbl1[ (Pmac_hap1 != Pmac_hap2) & 
                          (Pste_hap1 != Pste_hap2) &
                          (Pmac_hap1 == Pste_hap2) &
                          (Pmac_hap2 == Pste_hap1), .(pos)]
TSP_sites1

# Examine some of the TSP SNPs manually in IGV.
# 4362822 is of particular interest as it's the TSP SNP that seems to leapfrom species tree SNPs near 3UTR of L30
# but seems likely error as called at edge of A-repeat, where repeat is extended only in P. stenoptera g.
# (variant calling didn't call indels bc called on single read haplotypes. in the future could consider padding
# to simulate multiple reads to force indel call and filter based on distance.
TSP_sites1 <- TSP_sites1[pos != 4362822]

# other SNPs which are adjacent to repeat-indels (within 3 bp)
rmTSP <- c(4293290,4328072,4351494, 4354866)
TSP_sites1 <- TSP_sites1[!pos %in% rmTSP]

SP_sites1 <- TSP_tbl1[(Pste_hap1 == Pste_hap2) &
                        (Pmac_hap1 == Pmac_hap2) &
                        (Pste_hap1 != Pmac_hap1), .(pos)]

# likewise, manually check in IGV the rightmost of these which appears within the TSP region
# 4364812 is also a SNP called at the edge of an indel and appears to be error, remove. 
# other SNPs not manually checked bc this is the null expectation outside TSP region 
# (there are probably some errors within but the signal is overwhelming from plotting)
SP_sites1 <- SP_sites1[pos != 4364812]

# fwrite(TSP_tbl1[pos %in% TSP_sites1[,pos]], file = "~/workspace/Pterocarya/01_Gloc/LHS_wide/TSP_SNPs_fltd.txt",
#        row.names = F, quote = F, col.names = T, sep = "\t")

# fwrite(TSP_tbl1[pos %in% c(TSP_sites1[,pos], SP_sites1[, pos])], file = "~/workspace/Pterocarya/01_Gloc/LHS_wide/SNP_table_fltd.txt",
#        row.names = F, quote = F, col.names = T, sep = "\t")




# ----- RHS TSP -----
map2 <- fread("~/workspace/Pterocarya/01_Gloc/RHS_wide/toPsteHap2.map")
ped2 <- fread("~/workspace/Pterocarya/01_Gloc/RHS_wide/toPsteHap2.ped")
pos2 <- map2[, .(pos = V4)]

nuc2 <- ped2[,-c(2:6)]
setnames(nuc2, "V1", "ID")

# each column is duplicated immediately after because it thinks they're diploid
# how many SNPs?
mxn2 <- max( as.numeric( gsub("V", "", names(nuc2)[-1] ) ) )
SNP_nx2.2 <- mxn2 - 7 + 1
nSNP2 <- SNP_nx2.2/2

selectcols2 <- c("ID", paste0("V", seq(7, mxn2, by = 2)))
nuc2 <- nuc2[, ..selectcols2]

setnames(nuc2, paste0("V", seq(7, mxn2, by = 2)), as.character(1:(SNP_nx2.2/2)))

# to long, clean up for plotting
nuc2 <- melt(nuc2, id.vars = "ID", variable.name = "rank", value.name = "base")
nuc2[, pos := pos2[, pos], by = ID]
nuc2[,ID := gsub(".bam", "", basename(ID))]

nuc2[ID == 'P.ste_hap1', ID := 'g Pste']
nuc2[ID == 'P.ste_hap2', ID := 'G Pste']
nuc2[ID == 'P.mac_hap1', ID := 'G Pmac']
nuc2[ID == 'P.mac_hap2', ID := 'g Pmac']

nuc2[, base := factor(base, levels = c("0", "A", "T", 'C', "G"))]

# spread by ind
TSP_tbl2 <- dcast(nuc2, pos~ID, value.var = 'base')
TSP_tbl2

# confirm only biallelic sites
#nuc[, length(unique(base)), by = pos][, unique(V1)] # I think filtered in plink conversion
TSP_sites2 <- TSP_tbl2[ (Pmac_hap1 != Pmac_hap2) & 
                          (Pste_hap1 != Pste_hap2) &
                          (Pmac_hap2 == Pste_hap1) &
                          (Pmac_hap1 == Pste_hap2), .(pos)]
TSP_sites2
# These are not obvious calling errors but are in contexts that seem like they might be repeat mutations
# for example 4464431 is a SNP right next to a A-repeat, 4464340 is right next to a T-repeat and immediately
# adjacent to 2 other SNPs, 4446426 also occurs immediately adjacent to another SNP

SP_sites2 <- TSP_tbl2[(Pste_hap1 == Pste_hap2) &
                        (Pmac_hap1 == Pmac_hap2) &
                        (Pste_hap1 != Pmac_hap1), .(pos)]
                        





# ----- Calculate TE content per window -----
Dxy[, window_start := window]
Dxy[, window_end := window_start + window_size - 1]

TEs_tst <- P.ste_G_TEs[, .(TE_chr = Chr, TE_start = Start, TE_end = End)]

#TEs_tst <- RM[, .(TE_chr = Chr, TE_start = Start, TE_end = End)]

Dxy_tst <- Dxy[!is.na(window_start) & qryGnom == 'Pste_hap1', .(alignmentID, Dxy_chr = refChr, 
                                                                Dxy_num, Dxy_denom, Dxy, 
                                                                Dxy_window_start = window_start, 
                                                                Dxy_window_end = window_end)]

# find all fragments that overlap grid windows
setkey(TEs_tst, TE_chr, TE_start, TE_end)

TE_overlaps <- foverlaps(Dxy_tst, TEs_tst, by.x=c("Dxy_chr","Dxy_window_start","Dxy_window_end"), by.y=c("TE_chr","TE_start","TE_end"))

# For each TE, calculate proportion of interval it covers
TE_overlaps[, prp := (pmin(Dxy_window_end, TE_end) - pmax(Dxy_window_start, TE_start))/window_size, by=.(Dxy_chr,Dxy_window_start,Dxy_window_end)]
TE_overlaps[is.na(prp), prp := 0]

# inspect window with multiple TEs
#TE_overlaps[, .N, by = Dxy_window_start]

# sum across proportion of TEs in window, if greater than 1 (ie. nested TEs), take 1
TE_content_windows <- TE_overlaps[, .(tot.prp = min(1,sum(prp))), by = .(Dxy_chr, Dxy_window_start, Dxy_window_end)]

# merge back with Dxy data
Dxy <- merge(Dxy, 
             TE_content_windows[, .(refChr = Dxy_chr, window_start = Dxy_window_start, window_end = Dxy_window_end, TE.prp = tot.prp)],
             by = c("refChr", "window_start", "window_end"))

# ----- calculate mean & CI -----

Dxy[, avgDxy := mean(Dxy, na.rm = T), by = qryGnom]
Dxy[, qnt99 := quantile(Dxy, 0.99, na.rm = T), by = qryGnom]
Dxy[, qnt95 := quantile(Dxy, 0.95, na.rm = T), by = qryGnom]

#ggplot(Dxy, aes(x = TE.prp, y = Dxy_denom)) + geom_point()

# ----- plot -----
# set filters. (1) exclude windows with fewer than 500 bp
# and that are >90% primarily TE
dnm <- 500
Dxy_flt <- Dxy[TE.prp < 0.5 & Dxy_denom >= dnm]
Dxy_flt[, avgDxy := mean(Dxy, na.rm = T)]
Dxy_flt[, qnt99 := quantile(Dxy, 0.99, na.rm = T)]
Dxy_flt[, qnt95 := quantile(Dxy, 0.95, na.rm = T)]

# finally, on LHS trim to alignments before the Grepeats
# this value is the END of the G-0 (determined from BLAST alignment)
# after which alignment is likely to be garbage a priori
# note Dxy windows are labelled by start position

mx_wn <- 4371311

Dxy_flt <- Dxy_flt[window <= mx_wn | window > 4.4e6]
# set filter for min # of aligned bp in reference
dxy_plt <- ggplot(Dxy_flt[qryGnom %in% c('Pste_hap1') & 
             refChr == 'Chr11' &
             window > mn  & 
             window < mx & 
             Dxy_denom >= dnm]) + 
  #geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = 0, ymax = qnt99), fill = 'lightgray') +
  #annotate('rect', xmin = -Inf, xmax = Inf, ymin = 0, ymax = Dxy_flt[qryGnom == 'Pste_hap1', qnt99][1], fill = 'gray80')+
  geom_hline(aes(yintercept = avgDxy), linetype = 2, color = 'gray')   + 
  geom_hline(aes(yintercept = qnt99), linetype = 3, color = 'gray')   + 
  
  
  geom_rug(data = SP_sites1, aes(x = pos, color = "Sp"), length = unit(.2, 'cm')) + 
  geom_rug(data = SP_sites2, aes(x = pos, color = 'Sp'), length = unit(.2, 'cm')) +
  
  geom_rug(data = TSP_sites1, aes(x = pos, color = 'TSP'), length = unit(.2, 'cm')) + 
  geom_rug(data = TSP_sites2, aes(x = pos, color = 'TSP'), length = unit(.2, 'cm')) +
  
  geom_line(data = Dxy_flt[qryGnom %in% c('Pste_hap1') & 
                         refChr == 'Chr11' &
                         window < 4.4e6  & 
                         window > mn],
            aes(x = window, y = Dxy, group = qryGnom)) + 
  geom_line(data = Dxy_flt[qryGnom %in% c('Pste_hap1') & 
                         refChr == 'Chr11' &
                         window > 4.4e6  & 
                         window < mx],
            aes(x = window, y = Dxy, group = qryGnom)) +
  
  
  #geom_ribbon(data = helitron[start > mn & end < mx], aes(x = mid, ymin = 0, ymax = prp)) + 
  
  
  theme_classic() + 
  labs(y = expression(D[XY]), x = 'Position (Mb)') +
  scale_color_manual(values = c("Sp" = "gold", "TSP" = "blue"), name = NULL) + 
  scale_x_continuous(expand = c(0,0), breaks = seq(4.3e6, 4.5e6, length.out = 5), 
                     labels = sprintf("%.2f", seq(4.3, 4.5, length.out = 5))) +
  scale_y_continuous(expand= c(0,0)) + # , breaks = c(0, 0.05, .1, .15)) +
  theme(
    aspect.ratio = .2,
    plot.margin = margin(5,5,5,5, "pt"),
    text = element_text(size = 10),
    axis.text = element_text(size = 8),
    plot.title = element_text(face = 'italic'),
    axis.title.x = element_text(size = 10, vjust = 1),
    axis.title.y = element_text(size = 10),
    legend.position = c(0.2, 0.6),
    legend.text = element_text(size = 8),
    legend.margin = margin(0,0,0,0),
    legend.key.size = unit(.4, 'cm')
  )  
dxy_plt



# ---------- Helitron density -------
helitron <- fread("~/workspace/Pterocarya/denovo_assembly/scaffolded/P.stenoptera_UCD1470/hap2/Pste_dom_Chr11_helitron_density_5kb.txt")
setnames(helitron, c("V2", "V3", "V7"), c("start", "end", "prp"))
helitron[, mid := start + 500]


ggplot()











      