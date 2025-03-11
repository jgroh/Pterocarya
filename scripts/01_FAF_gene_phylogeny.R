library(data.table)
library(ggplot2)
library(ape)
library('ggtree')
library(phytools)
library(pals)
library(ggnewscale)
library(dplyr)

setwd("~/workspace/Pterocarya/01_Gloc/genes/FAF/CDS/")


# read tree
tree <- ape::read.tree("FAF_domain/codon_align.fasta.treefile")

# read tip data
tipdata <- setDT(read.csv("FAF_domain/tipdata.csv", header = F,
                    col.names = c("tip", "species", "locus", "genename", 'shape')))
tipdata[grepl("stenoptera", species), focal := 'A']
tipdata[!grepl("stenoptera", species), focal := 'B']

# format labels
tipdata[, label_expression := paste0(
  "italic('", species, "')",  # Italic species name
  " * ' |  ' * ",             # Separator
  "' ", locus, "'",            # Plain locus
  " * '' * ",             # Another separator
  "italic('", genename, "')"
)]
# root using Marchantia as outgroup
tree <- phytools::reroot(tree, node.number = which(tree$tip.label == 'Marchantia_polymorpha|Mapoly0824s0001'),
                         position = 1)


tree$node.label[tree$node.label == 'Root'] <- ""
p <- ggtree(tree) %<+% tipdata
p + geom_text(aes(label=node), hjust=-0.3) + 
  geom_tiplab() #+ geom_nodelab()

# q <- rotate(p, 28) 
# q <- rotate(q, 44)
# q <- rotate(q, 48)
# q <- rotate(q, 33)
# q <- rotate(q, 37)
# q <- rotate(q, 42)

# q + xlim(0, 20) + 
#   geom_highlight(mapping=aes(subset = node %in% c(44, 29)), 
#                  type = 'gradient', 
#                  fill = c('orange','mediumvioletred'),
#                  gradient.direction = 'rt', alpha = 0.7) + 
#   geom_tree() +
#   geom_treescale(x = 0, width = 1) + 
#   geom_nodelab(size = 2.5, nudge_x = -0.3, nudge_y = 0.25)  + 
#   geom_tiplab(aes(label =label_expression), parse = T, size = 3) 
  

# use Arial font
library(showtext)
showtext_auto()
library(sysfonts)
#font_add("Arial", regular = "/System/Library/Fonts/Supplemental/Arial.ttf",
#         italic = "/System/Library/Fonts/Supplemental/Arial Italic.ttf")

p
p + xlim(0,30) + 
  geom_highlight(mapping=aes(subset = node %in% c(50, 35)), 
                 type = 'gradient', 
                 fill = c(#'#C1CA2D',
                    '#55AACC',
                          '#882255'),
                 gradient.direction = 'rt', alpha = 0.7) + 
  geom_tree() +
  geom_tippoint(aes(shape = shape, color = shape, size = shape)) +
  scale_shape_manual(values = c(18, NA, 20)) +
  scale_size_manual(values = c(2, NA, 1.5)) +
  
  
  geom_treescale(x = 0, width = 1,fontsize = 3, offset=0.1) + 
  geom_nodelab(size = 2, nudge_x = -0.35, nudge_y = 0.28)  + 
  geom_tiplab(aes(label =label_expression), 
              parse = T, size = 2.5, offset = 0.2) + 
              #family = 'Arial')  + 
  scale_color_manual(values = c("red", 'black', 'black')) +
  theme(legend.position = 'none')




