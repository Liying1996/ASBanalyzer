# Convert ensembleID to gene symbol
# R --slave --no-restore --file=${basepath}/convertID.R --args ${output_dir}/html_summary/${name}_all_summary.txt ${output_dir}/html_summary/${name}_all_summary_ID.txt

library(stringr)
library(org.Hs.eg.db)

args <- commandArgs(trailingOnly=TRUE)

summary_file_name <- args[1]
output_file <- args[2]

summary_file <- read.csv(summary_file_name, header=FALSE, sep='\t')
k <- keys(org.Hs.eg.db,keytype = "ENSEMBL")

match_list <- select(org.Hs.eg.db,keys=k,columns = c("ENTREZID","SYMBOL"), keytype="ENSEMBL")

gene_symbols <- c()
for (gene in as.character(summary_file$V19)){
	if (gene != '-'){
		genes <- str_split(gene, ',')
		new_genes1 <- as.vector(unlist(genes[1]))
		new_genes2 <- sub("[.][0-9]*", "", new_genes1)

		res <- ""
		ID_list <- match_list[match(new_genes2, match_list[,"ENSEMBL"]), 3]

		for (i in 1:length(new_genes1)){
			if (is.na(ID_list[i])){
				res2 <- paste(new_genes1[i], '(-)', sep='')
			}
			else{
				res2 <- paste(new_genes1[i], '(', ID_list[i], ')', sep='')
			}
			if (res != ""){
				res <- paste(res, res2, sep=',')
			}else{
				res <- res2
			}
		}
	}
	else{res <- '-'}

	gene_symbols <- c(gene_symbols, res)
}

new_summary_file <- summary_file
new_summary_file$V19 <- gene_symbols
write.table(new_summary_file, file=output_file, sep='\t', row.names=FALSE, col.names=FALSE, quote=FALSE)

