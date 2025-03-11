library(data.table)
library(ggplot2)
library(ape)
library('ggtree')
library(phytools)
library(pals)
library(ggnewscale)
library(dplyr)

setwd("~/workspace/Pterocarya/01_Gloc/Grepeats/Cyclo_Ptero")

# read tree
tree <- ape::read.tree("aligned-curated.fasta.treefile")
tipdata <- fread("tipdata.csv")

# format labels
tipdata[, label_expression := paste0(
  "italic('", species, "')",  # Italic species name
  " * '' * ",             # Another separator
  "' ", ind, "'",            # Individual
  " * ' | ' * ",             # Separator
  "italic('", allele, "')",
  " * '-' * ",             
  "italic('", no, "')"
)]
tipdata[species == 'J. regia', label_expression := "italic('J. regia')"]

# root 
tree <- phytools::reroot(tree, node.number = which(tree$tip.label == 'Jreg_0'), position = .01)

# quick look
ggtree(tree)  %<+% tipdata + geom_tiplab() + geom_text(aes(label=node), hjust=-0.3) 


# use Arial font
library(showtext)
showtext_auto()
library(sysfonts)
font_add("Arial", regular = "/System/Library/Fonts/Supplemental/Arial.ttf",
         italic = "/System/Library/Fonts/Supplemental/Arial Italic.ttf")

tree$node.label[tree$node.label == 'Root'] <- ""

p <- ggtree(tree) %<+% tipdata 

full.phylo <- p + xlim(0, 0.25) + 
  geom_highlight(mapping=aes(subset = node %in% c(60, 62, 65, 67)), 
                 type = 'gradient', 
                 fill = c('darkmagenta','lightsalmon', 'darkblue', 'turquoise'),
                 gradient.direction = 'rt', alpha = 0.7) + 
  geom_tree() +
  geom_tippoint(aes(shape = allele)) +
  scale_shape_manual(values = c(NA, NA, NA, NA, 20)) +
  geom_treescale(fontsize = 2, x = 0, offset = 0.2) + 
  geom_nodelab(size = 2, nudge_x = -.0035, nudge_y = 0.27)  + 
  geom_tiplab(aes(label = label_expression), 
              parse = T, size = 3, family = 'Arial', offset = 0.001)  + 
  theme(legend.position = 'none', aspect.ratio = 1.3) 

# geom_strip('Cpal_4PA.13C_G2-0', 'Cpal_4PA.13C_G2-2', barsize=1, color='darkblue', 
  #            label = "G2", offset.text=0.005) + 
  # geom_strip('Cpal_4PA.13B_G1-0', 'Cpal_4PA.13D_G1-0', barsize=1, color='turquoise', 
  #            label = "G1", offset.text=0.005) +  
  # geom_strip('Pmac_g-0', 'Pste_g-0', barsize=1, color='lightsalmon', 
  #            label = "g", offset.text=0.005)   +
  # geom_strip('Pmac_G-0', 'Pste_G-13', barsize=1, color='darkmagenta', 
  #            label = "g", offset.text=0.005, offset = 0.06)   
full.phylo

