library(ggplot2)

args <- commandArgs(trailingOnly=TRUE)

chrom_file <- args[1]
AS_file <- args[2]
output <- args[3]

chrom_info <- read.table(chrom_file, header=F)
colnames(chrom_info) <- c('chr', 'len')
chrom_info$chr <- factor(chrom_info$chr, levels=paste('chr', c(seq(1,22), 'X'), sep=''))

snps <- read.csv(AS_file, header=T, sep = '\t')
snps$betabin.FDR <- as.character(snps$betabin.FDR)
snps$betabin.FDR[which(snps$betabin.FDR=='significant')] <- 'ASB'
snps$betabin.FDR[which(snps$betabin.FDR=='not significant')] <- 'non ASB'

pdf(output, width=12, height=7)

ggplot() + geom_bar(data = chrom_info, aes(x = chr, y = len/1000000), stat = "identity", fill='white', color='gray', width = 0.5) + 
  theme_classic() + geom_point(data = snps, aes(x=Chr, y=Pos/1000000, color=betabin.FDR, size=betabin.FDR, alpha=betabin.FDR), shape=95) +
  scale_color_manual(values=c('#791E94','#ffa200')) + scale_alpha_manual(values=c(1, 0.3)) + scale_size_manual(values=c(12,7)) +
  guides(alpha='none') +  guides(color='none') + guides(color=guide_legend(title=NULL)) + guides(size=guide_legend(title=NULL)) +
  guides(alpha=guide_legend(title=NULL)) + theme(legend.position = 'top') + 
  theme(legend.text = element_text(size=15, face="bold")) +
  theme(axis.text.x = element_text(size = 15, face = "bold", vjust = 0.5, hjust = 0.5, angle = 45)) +
  labs(x = '', y='len(Mb)')
dev.off()

