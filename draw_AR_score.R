library(ggpubr)
library(ggplot2)

args <- commandArgs(trailingOnly=TRUE)

AR_score_file <- args[1]
biosample <- args[2]
output_file <- args[3]

AR_score <- read.csv(AR_score_file, header=F, sep = '\t', stringsAsFactors = FALSE)
AR_score$diff <- AR_score$V14 - AR_score$V15


AR_score$V13[which(AR_score$V13 != 'significant')] <- 'nonAS'
AR_score$V13[which(AR_score$V13 == 'significant')] <- 'AS'

AR_score$V9 <- as.numeric(as.character(AR_score$V9))

pdf(output_file, width = 9, height=7)
ggplot(AR_score, aes(x=diff, y=V9, color=V13, shape=V13)) + geom_point() + 
  scale_color_manual(values = c('#eb8a3c', '#2b9464')) + 
  geom_smooth(method = 'lm') + 
  coord_cartesian(ylim = c(0,1))  + 
  stat_cor(method = 'pearson', aes(x=diff, y=V9, color=V13))  + 
  labs(x = 'Motif score change', y = 'Ref allele Ratio', title = biosample) +
  theme_classic() + theme(plot.title = element_text(hjust=0.5)) +
  guides(color=guide_legend(title=NULL), shape=guide_legend(title=NULL))

dev.off()
