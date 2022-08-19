# Draw the distribution of 10000 nonAS control simulations
# example : Rscript draw_motif.R ENCFF001HIA_inmotif_control_count.txt 10  ENCFF001HIA ENCFF001HIA_dist.pdf 

library(ggplot2)
library(Unicode)

args <- commandArgs(trailingOnly=TRUE)

controls <- read.table(args[1], header=F)
AS_num <- args[2]
name <- args[3]

colnames(controls) <- 'control_num'

AS_num <- as.numeric(AS_num)
zscore <- (AS_num - mean(controls$control_num)) / sd(controls$control_num)
p.val <- min(2*pnorm(q=zscore, lower.tail=FALSE),1)
p.val <- signif(p.val, 3)
write.table(p.val, 'tmp', quote = FALSE, col.names = FALSE,row.names = FALSE)
p.val <- paste(' = ', p.val, sep='')

#italic_p <- u_char_inspect(u_char_from_name("MATHEMATICAL ITALIC SMALL P"))["Char"]
pdf(args[4], height=5, width=8)
#cairo_pdf(args[4], height=5, width=8)
ggplot(data=controls, aes(x = control_num)) + geom_histogram(stat="count", fill='#888888') + labs(x='Number of control SNPs in motifs', y = 'count', title = name) +
  geom_vline(xintercept = AS_num, color='#FF5511') +
  geom_vline(xintercept = mean(controls$control_num), color='#0066FF') +
   theme_bw() + theme(plot.title = element_text(hjust=0.5)) +
  annotate("text", x = -Inf, y = Inf, label = substitute(paste(italic('p'), x), list(x=p.val)), hjust = -.2, vjust = 2)
dev.off()	
