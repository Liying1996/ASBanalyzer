library(ggsci)
library(ggplot2)

args <- commandArgs(trailingOnly=TRUE)

anno_file <- args[1]
biosample <- args[2]
output_file <- args[3]

anno <- read.csv(anno_file, header=F, sep = '\t')
anno$V8 <- gsub(pattern = "ANN=.\\|", replacement = '', anno$V8)
anno$V8 <- gsub(pattern = "_variant", replacement = '', anno$V8)

var <- c('intron', 'intergenic_region', 'downstream_gene', 'upstream_gene', '3_prime_UTR', '5_prime_UTR', 'missense')
anno$group <- anno$V8
anno[!anno$V8 %in% var, 'group'] <- 'Other'

pdf(output_file, height = 5, width = 10)


anno$group <- factor(anno$group, levels = c(var, 'Other'))
# num <- length(unique(anno$V8))
#num <- length(names)
num <- 8
ggplot(data=anno, aes(x = V3, fill=group)) + 
  theme_classic() + 
  geom_bar(stat="count",width=0.5,position='fill') + coord_flip() + 
  theme(legend.position = 'top') + scale_fill_manual(values=pal_jco('default')(num)) + 
  labs(x = '' ,y = 'Frequency', title=paste('Variant annotation (', biosample, ')', sep='')) + 
  theme(plot.title=element_text(hjust=0.5)) + 
  scale_x_discrete(breaks=unique(anno$V3), labels = c('AS', 'nonAS')) + 
  guides(fill=guide_legend(title=NULL))

dev.off()


