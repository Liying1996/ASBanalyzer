### **ASBanalyzer manual**

---

#### Introduction

ASBanalyzer is a pipeline of analyzing allele-specific binding (ASB) sites.

The following directories and files are included with ASBanalyzer:

get_AS.sh: The main pipeline;

remove_bias_SE.sh/remove_bias_PE.sh: WASP to remove mapping bias;

betabinomial.R: Use Chen's method to obtain ASB sites;

draw_cCRE_pie.py: Draw pie charts of cCREs distribution;

draw_cCRE_hist.R: Draw histograms of cCREs distribution;

draw_freq.py: Plot the positions distribution of SNPs that disrupt the motif sequence；

draw_annotaion.R: Draw  histograms of snpEFF annotated regions；

draw_motif.R: Draw histograms of motif enrichment; 

draw_AR_score.R: Plot scatter plots of reference allele ratio and motifs scores of 2 alleles;

best_match.py: Select best Fimo outputs; 

motif_disrupt.py: Obtain SNPs' positions of motifs and calculate the frequency of pos of disruptions;

trans_alt.py: Convert the fasta with reference allele to fasta with alternative allele;

motif_score.sh: Calculate the motif scores of 2 alleles;

data: human and mouse cCREs from SCREEN database;

example_data:  Example data files that can be used to test this pipeline, including single-end and paired-end reads.

---

#### Usage and parameters

-s: Single-end reads, cannot specify option -s after specifying option -p;

-p: Paired-end reads, the 2 files are separated by a comma (','), cannot specify option -p after specifying option -s;

-w: The directory of WASP;

-i: The BWA index (Not provided in the example_data/);

-c: The ChromInfo file;

-h: The HDF5 files containing SNPs information converted from VCF, cannot specify option -h after specifying option -v;

-v: The VCF file containing SNPs, cannot specify option -v after specifying option -h;

-f: The peak file of ChIP-seq/DNase-seq data;

-j: The path of snpEff.jar;

-a: Species, include human and mouse. Human: hg38; mouse: mm10;

-o: Output folder;

-m: The motifs of a transcript factor (you can download from [here](https://meme-suite.org/meme/meme-software/Databases/motifs/motif_databases.12.21.tgz));

-g: Significant variant-gene associations downloaded from [GTEx](https://storage.googleapis.com/gtex_analysis_v8/single_tissue_qtl_data/GTEx_Analysis_v8_eQTL.tar); 

-r: Reference fasta file.


For example, single-end reads:

```shell
path=~/ASBanalyzer/
bash get_AS.sh \
	-s ${path}/example_data/single_end/ENCFF000OCP.fastq.gz \
	-w ~/WASP/ \
	-i ~/hg38_index/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta \
	-c  ${path}/example_data/chromeinfo/GRCh38_EBV.chrom.sizes.tsv \
	-h ${path}/example_data/h5 \
	-f ${path}/example_data/ENCFF801BDJ.bed \
	-j ${path}/data/snpEFF/snpEff/snpEff.jar \
	-a human \
	-o ${path}/example_data/output/single/
	-m ${path}/example_data/CTCF.meme \
	-g ~/GTEx/GTEx_Analysis_v8_eQTL/Whole_Blood.v8.signif_variant_gene_pairs.txt \
	-r ~/hg38_index/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta 
```

Paired-end reads:

```shell
path=~/ASBanalyzer/
bash get_AS.sh \
	-p ${path}/example_data/paired_end/ENCFF340SQP.fastq.gz,${path}/example_data/paired_end/ENCFF587OVW.fastq.gz \
	-w ~/WASP/ \
	-i ~/hg38_index/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta \
	-c  ${path}/example_data/chromeinfo/GRCh38_EBV.chrom.sizes.tsv \
	-h ${path}/example_data/h5 \
	-f ${path}/example_data/ENCFF801BDJ.bed \
	-j ${path}/data/snpEFF/snpEff/snpEff.jar \
	-a human \
	-o ${path}/example_data/output/paired/
	-m ${path}/example_data/CTCF.meme \
	-g ~/GTEx/GTEx_Analysis_v8_eQTL/Whole_Blood.v8.signif_variant_gene_pairs.txt \
	-r ~/hg38_index/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta 
```

---

#### Conda Environment
You can directly download the conda environment in: []

After downloaded, you'll get a file named 'ASB_env.tar.gz'

Then Use `tar -zxvf` to unzip;
Finally use the following code to activate the conda environment:
`source ./ASB_env/bin/activate`


#### Code and Data
You can direcly download from: []

Thus, you didn't need to download WASP and other datasets by yourself.


#### Dependencies

If you want to build an environment by yourself, the following are dependencies:

Python3 (Highly recommend Anaconda):

- numpy
- scipy
- pysam
- pytables
- matplotlib
- pandas
- argparse
- All [WASP](https://github.com/bmvdgeijn/WASP) needs

R:

- VGAM
- Unicode
- ggplot2
- ggpubr
- ggsci


Softwares:

- fastp
- Picard
- tabix
- bedtools intersect (aka. intersectBed)
- snpEFF
- fimo

Conda can install these softare and packages directly:

```shell
conda install picard, fastp,tabix, bedtools
```
