library(data.table)
library(ggplot)
library(pals)
library(forcats)

reg <- fread("~/workspace/Pterocarya/RNAseq/candidate_expression_in_Juglans_regia.csv")

reg$Sample <- forcats::fct_reorder(reg$Sample, reg$Expression)
ggplot(reg, aes(x = Sample, y = Expression)) + 
  geom_bar(stat = 'identity', fill = 'firebrick') + 
  #scale_fill_viridis_b(option = 'G') +
  #scale_fill_manual(values = rev(stepped(20))) + 
  labs(x = '', y = 'FPKM', title = 'LOC109012520') +
  theme_classic() +
  theme(aspect.ratio = 1,
        axis.text.x = element_text(size = 11, angle = 90, vjust = 0.5, hjust=1),
        axis.text.y = element_text(size = 11),
        legend.position= 'none')
