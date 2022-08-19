library(ggpubr)
library(ggplot2)

args <- commandArgs(trailingOnly=TRUE)

AR_score_file <- args[1]
biosample <- args[2]
output_file <- args[3]


pwm_score <- read.csv(AR_score_file, header=F, sep = '\t')
pwm_score$diff <- pwm_score$V14 - pwm_score$V15

AS_score <- pwm_score[pwm_score$V13=="significant",]
nonAS_score <- pwm_score[pwm_score$V13!="significant",]
#sampled_nonAS <- nonAS_score[sample(nrow(nonAS_score), nrow(AS_score)),]
#merged_pwm_score <- rbind(AS_score, sampled_nonAS)
merged_AR_score <- rbind(AS_score, nonAS_score)
merged_pwm_score$V13 <- c(rep('AS', nrow(AS_score)), rep('control', nrow(nonAS_score)))

pdf(output_file, width = 9, height=7)
ggplot(merged_pwm_score, aes(x=diff, y=V9, color=V13, shape=V13)) + geom_point() + 
  scale_color_manual(values = c('#eb8a3c', '#2b9464')) + 
  geom_smooth(method = 'lm') + 
  coord_cartesian(ylim = c(0,1))  + 
  stat_cor(method = 'pearson', aes(x=diff, y=V9, color=V13))  + 
  labs(x = 'PWM score', y = 'Ref allele Ratio', title = biosample) +
  theme_classic() + theme(plot.title = element_text(hjust=0.5)) +
  guides(color=guide_legend(title=NULL), shape=guide_legend(title=NULL))

dev.off()
