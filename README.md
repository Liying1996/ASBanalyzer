### **ASBanalyzer manual**

---

#### Introduction

ASB-analyzer is a pipeline for visualizing and analyzing allele-specific binding (ASB) SNPs.

The following directories and files are included with ASBfinder:

get_AS.sh: The main pipeline;

remove_bias_SE.sh/remove_bias_PE.sh: WASP to remove mapping bias;

betabinomial.R: Use Chen's method to obtain ASB SNPs;

cal_motif.sh: Use Fimo to scan motif sequences and calculate the motif enrichment of ASB and control (non-ASB) SNPs;

Chromosome.R: Draw the  distribution of ASB SNPs across 23 chromosomes；

convertID.R: Convert the Ensembl gene ID to gene symbol (for better display);

draw_cCRE.py: Draw pie charts of cCREs distribution SNPs enriched in;

draw_cCREs_hist.R: Draw histograms of cCREs distribution SNPs enriched in;

draw_annotation.R: Draw the genomic regions annotated by snpEFF;

draw_motif.R: Draw histograms of motif enrichment;

draw_PWM_score.R: Draw  scatter plots of  allele ratio and PWM score change of SNPs;

draw_freq.R: Draw the motif information content and the frequency of each postion of motif disrupted by SNPs;

best_match.py: Select best Fimo outputs;

motif_disrupt.py: Obtain SNPs' positions of motifs disrupted and calculate the frequency of pos of disruption;

motif_score.py: Get the PWM scores of ref and alt alleles;

trans_alt.py: Convert the fasta with reference allele to fasta with alternative allele;

data/ : human and mouse cCREs from SCREEN database;

example_data/ :  Example data files that can be used to test this pipeline, including single-end and paired-end reads.

---

#### Installation

You can download the whole software and related datasets via:

Then,

***Activate the conda enviroment***

```shell
cd ASBanalyzer/
source ASB_env/bin/activate
```

***Deactivate***

```
source ASB_env/bin/activate
```

---

#### Usage and parameters

***Required***:

-s: Single-end reads, cannot specify option -s after specifying option -p;

-p: Paired-end reads, the 2 files are separated by a comma (','), cannot specify option -p after specifying option -s;

-h: The HDF5 files containing SNPs information converted from VCF, cannot specify option -h after specifying option -v;

-v: The VCF file containing heterozygous SNPs, cannot specify option -v after specifying option -h;

-m: The motifs of a transcript factor (you can download from [here](https://meme-suite.org/meme/meme-software/Databases/motifs/motif_databases.12.21.tgz));

-o: Output folder;

***Optional***:

-i: The BWA index;

-c: The ChromInfo file;

-f: The peak file of ChIP-seq/DNase-seq data (default: called by macs2);

-g: Significant variant-gene associations downloaded from [GTEx](https://storage.googleapis.com/gtex_analysis_v8/single_tissue_qtl_data/GTEx_Analysis_v8_eQTL.tar);

-r: Reference fasta file.



***Example***


single-end reads (only required parameters):

```shell
path=~/test_ASB/ASB-analyzer/
sample=ENCFF001HIA
bash get_AS.sh \
  -s ${path}/example_data/single_end/${sample}.fastq.gz  \
	-h ${path}/example_data/h5 \
	-m ${path}/data/CTCF.meme \
	-o ${path}/example_data/output/test/
```

Paired-end reads:

```shell
bash get_AS.sh \
  -p ${path}/example_data/paired_end/ENCFF340SQP.fastq.gz,${path}/example_data/paired_end/ENCFF587OVW.fastq.gz  \
	-h ${path}/example_data/h5 \
	-m ${path}/data/CTCF.meme \
	-o ${path}/example_data/output/test/
```

---

#### Output Results

QC/: The QC results from fastp;

find_intersecting_snps/; map/; remap/; filter_remapped_reads/:  Intermediate files generated by WASP (for removing mapping bias);

dup/: Picard was implemented to remove the duplicates;

counts/: The read counts of ref and alt alleles were obtained by WASP-bam2h5.py;

inpeak/ : Only the SNPs in peaks were considered;

annotation/: Output and visualize the genomic regions and cCREs (candidate cis-regulatory elements) SNPs enriched in;

motif/: The results of motif enrichment, the associations of positions of motif disruption and motif information content, the correlations of allele ratio and PWM score change of ref and alt alleles.

screenshots: Genome browser screenshots of ASB SNPs (+-20bp);

html_summary/: An HTML summary of all results and plots above.
![]('example.png')





