library(ggsci)
library(ggplot2)

args <- commandArgs(trailingOnly=TRUE)

anno_file <- args[1]
biosample <- args[2]
output_file <- args[3]

anno <- read.csv(anno_file, header=F, sep = '\t')
anno$V8 <- gsub(pattern = "ANN=.\\|", replacement = '', anno$V8)
anno$V8 <- gsub(pattern = "_variant", replacement = '', anno$V8)

names <- unique(anno$V8)
var <- rownames(table(anno$V8)[table(anno$V8) > 5])
anno <- anno[anno$V8 %in% var, ]

pdf(output_file, height = 5, width = 10)

# levels(anno$V8) <- c('intron', 'upstream_gene', 'missense', 'intergenic_region', 'downstream_gene', 'synonymous', '3_prime_UTR', '5_prime_UTR', 'splice_region&intron', 'sequence_feature', 'non_coding_transcript_exon', '5_prime_UTR_premature_start_codon_gain','test')

levels(anno$V8) <- names
# num <- length(unique(anno$V8))
num <- length(names)
ggplot(data=anno, aes(x = V3, fill=V8)) + 
  theme_classic() + 
  geom_bar(stat="count",width=0.5,position='fill') + coord_flip() + 
  theme(legend.position = 'top') + scale_fill_manual(values=pal_igv('default')(num)) + 
  labs(x = '' ,y = 'Frequency', title=paste('Variant annotation (', biosample, ')', sep='')) + 
  theme(plot.title=element_text(hjust=0.5)) + 
  scale_x_discrete(breaks=unique(anno$V3), labels = c('AS', 'nonAS')) + 
  guides(fill=guide_legend(title=NULL))

dev.off()


