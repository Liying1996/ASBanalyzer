library(ggpubr)

args <- commandArgs(trailingOnly=TRUE)

AR_score_file <- args[1]
biosample <- args[2]
output_file <- args[3]

AR_score <- read.csv(AR_score_file, header=F, sep = '\t')
AR_score$diff <- AR_score$V14 - AR_score$V18

AS_score <- AR_score[AR_score$V21=="AS",]
nonAS_score <- AR_score[AR_score$V21=="nonAS",]
sampled_nonAS <- nonAS_score[sample(nrow(nonAS_score), nrow(AS_score)),]
merged_AR_score <- rbind(AS_score, sampled_nonAS)
#merged_AR_score[merged_AR_score$V21 == 'nonAS',]$V21 = 'control'
merged_AR_score$V21 <- c(rep('AS', nrow(AS_score)), rep('control', nrow(AS_score)))

pdf(output_file, width = 9, height=7)
ggplot(merged_AR_score, aes(x=diff, y=V9, color=V21, shape=V21)) + geom_point() + 
  scale_color_manual(values = c('#eb8a3c', '#2b9464')) + 
  geom_smooth(method = 'lm') + 
  coord_cartesian(ylim = c(0,1))  + 
  stat_cor(method = 'pearson', aes(x=diff, y=V9, color=V21))  + 
  labs(x = 'Motif score change', y = 'Ref allele Ratio', title = biosample) +
  theme_classic() + theme(plot.title = element_text(hjust=0.5)) +
  guides(color=guide_legend(title=NULL), shape=guide_legend(title=NULL))

dev.off()
