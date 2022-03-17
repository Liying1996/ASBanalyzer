library(ggplot2)
library(ggseqlogo)
library(ggpubr)
library(Unicode)

args <- commandArgs(trailingOnly=TRUE)

freq_file <- args[1]
ppm_file <- args[2]
output_file <- args[3]
biosample <- args[4]

disrupt_pos <- read.table(freq_file, sep = '\t', header=F)
colnames(disrupt_pos) <- c('type', 'pos', 'count', 'total', 'freq')

g1 <- ggplot() + geom_line(data = disrupt_pos, mapping = aes(x = pos, y = freq, color = type, linetype=type),size=1.4) +
  geom_point(data = disrupt_pos, mapping = aes(x = pos, y = freq, color = type, shape=type),size=3) + 
  geom_bar(data = disrupt_pos, mapping = aes(x = pos, y = freq, fill = type),position="dodge",stat="identity",alpha=0.7,color='gray') + 
  scale_color_manual(values = c('#eb8a3c', '#2b9464')) + scale_fill_manual(values = c('#eb8a3c', '#2b9464')) + 
  scale_x_continuous(breaks=seq(1,19,1))  + 
  labs(x = 'Position', y = 'Fraction', title = biosample) + 
  theme_classic() + theme(legend.position = c(1, 1),legend.justification = c(1,1)) + 
  theme(plot.title = element_text(hjust=0.5)) + 
  guides(color=guide_legend(title=NULL), fill=guide_legend(title=NULL), shape=guide_legend(title=NULL), linetype=guide_legend(title=NULL))

# MOTIF LOGO
ppm <- read.table(ppm_file, header=F)
ppm <- t(ppm) 
rownames(ppm) <- c("A", "C", "G", "T")
g2 <- ggseqlogo(ppm)
print(g2)

cairo_pdf(output_file, height=7, width=10)
ggarrange(g2, g1, 
          ncol = 1, nrow = 2, 
          widths = c(1, 3), heights = c(1,2))
dev.off()


