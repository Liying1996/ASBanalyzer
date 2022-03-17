library(ggplot2)

args <- commandArgs(trailingOnly=TRUE)

ccre_file <- args[1]
biosample <- args[2]
output_file <- args[3]

ccre <- read.csv(ccre_file, header=F, sep = '\t')
ccre$V20 <- gsub(pattern = ',CTCF-bound', replacement = '', ccre$V20)

# ggplot(data=ccre, aes(x = V14, fill=V20)) +
#   geom_bar(stat="count",width=0.5,position='fill') +
#   geom_text(stat='count',aes(label=scales::percent(..count../sum(..count..)), group=V20) ,color="white", size=3.5,position=position_fill(0.5))

ccre$V20 <- factor(ccre$V20, levels = c('dELS', 'pELS', 'PLS', 'DNase-H3K4me3', 'CTCF-only', 'Unclassified'))

pdf(output_file, width=10, height=5)
ggplot(data=ccre, aes(x = V14, fill=V20)) +
  geom_bar(stat="count",width=0.5,position='fill') +
  guides(fill=guide_legend(title=NULL)) + coord_flip() +
  labs(x = '', y = 'Frequency', title=paste('cCREs distribution (', biosample, ')', sep='')) + 
  theme_classic() + 
  scale_fill_manual(values=c('#FFCD00', '#FFA700', '#FF0000', '#ffaaaa', '#00B0F0', '#8c8c8c')) + 
  theme(plot.title=element_text(hjust=0.5)) + 
  scale_x_discrete(breaks=unique(ccre$V14), labels = c('AS', 'nonAS')) 

 dev.off()
