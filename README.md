### **ASBanalyzer manual**

---

#### 1. Introduction

Allele-specific binding (ASB) events occur when transcription factors (TFs) bind more favorably to one of the two parental alleles at heterozygous single  nucleotide polymorphisms (SNPs).

ASB-analyzer is a pipeline for visualizing and analyzing allele-specific binding (ASB) SNPs.  It is a software platform that enables the users to quickly and efficiently input raw sequencing data to generate individual reports containing the cytogenetic map of ASB SNPs and their associated phenotypes.  This interactive tool thereby combines ASB SNP identification, biological annotation, motif analysis, GWAS  associations and report summary in one pipeline.

---

#### 2. Installation

The full-version software and related datasets can be downloaded via this [link](https://drive.google.com/file/d/1bNkO9bY2AIW-Lw80Y6S53Sgtm8Gh0OFX/view?usp=sharing).

Next, please proceed according to the following steps:


***File decompression*** 

```
 tar -jxvf ASBanalyzer.tar.bz2 -C ASBanalyzer_env/
```

***Activate the conda enviroment***

```shell
cd ASBanalyzer/
mkdir ASBanalyzer_env/
# tar zxvf ASBanalyzer_env.tar.gz -C ASBanalyzer_env/
source ASBanalyzer_env/bin/activate
```

***Build BWA index***

If you want to use your own reference genome and index, please skip this step.

```
bwa index GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta
```

***Build annovar database***

```
perl data/annovar/annotate_variation.pl --downdb --webfrom annovar --buildver hg38 avsnp150 data/annovar/humandb/
```
---

#### 3. Usage and parameters

***Required***:

-s: Single-end reads, cannot specify option -s after specifying option -p;

-p: Paired-end reads, the 2 files are separated by a comma (','), cannot specify option -p after specifying option -s;

-h: The HDF5 files containing SNPs information converted from VCF (users can also directly provide the VCF file using the `-v` option, and the built-in WASP in the pipeline will convert it to the HDF5 format), cannot specify option -h after specifying option -v;

-v: The VCF file containing heterozygous SNPs, cannot specify option -v after specifying option -h;

-m: The motif file of a transcript factor (users can download [here](https://meme-suite.org/meme/meme-software/Databases/motifs/motif_databases.12.21.tgz));

-o: Output folder;

***Optional***:

-i: The BWA index;

-c: The ChromInfo file (The workflow utilize `example_data/chromeinfo/GRCh38_EBV.chrom.sizes.tsv` as a default file. Alternatively, users have the option to directly access the data from the UCSC database.);

-f: The peak file of ChIP-seq/DNase-seq data (default: called by macs2 with default parameters);

-g: Significant variant-gene associations downloaded from [GTEx](https://storage.googleapis.com/gtex_analysis_v8/single_tissue_qtl_data/GTEx_Analysis_v8_eQTL.tar);

-r: Reference fasta file.



***Examples***


Single-end reads (required parameters):

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

#### 4. Detailed descriptions of each script

The following directories and scripts are included with ASBanalyzer:

get_AS.sh: The main pipeline. Users can directly obtain all analysis results using this main script. For the usage instructions, please refer to the "3. Usage and Parameters" section. ;

remove_bias_SE.sh/remove_bias_PE.sh: The scripts for mapping bias correction with WASP and duplicate removal with Picard. If users want to use them directly, follow the examples below:

```
# For paired-end data
sh ${basepath}/remove_bias_PE.sh \
-w ${wasp_path} \ # Path to WASP
-i ${bwa_index} \ # BWA index
-1 ${output_dir}/QC/${name}.fq1 \
-2 ${output_dir}/QC/${name}.fq2 \
-h ${h5_dir} \  # The HDF5 files containing SNPs information
-o ${output_dir}

# For single-end data
sh ${basepath}/remove_bias_SE.sh \
-w ${wasp_path} \
-i ${bwa_index} \
-f ${output_dir}/QC/${name}.fq \
-h ${h5_dir} -o ${output_dir}
```

betabinomial.R: The script to obtain ASB SNPs using Chen's Method.  If users want to use it directly, follow the examples below: 

```
R --slave --no-restore --file=${basepath}/betabinomial.R --args ${output_dir}/inpeak/${name}_counts.txt ${name} 0.1
```

cal_motif.sh: The script to scan motif sequences with Fimo and calculate the motif enrichment of ASB and control (non-ASB) SNPs. If users intend to directly employ the script, please refer to the provided examples:

```
sh ${basepath}/cal_motif.sh \
-i ${output_dir}/inpeak/${name}_counts_AS_0.1.txt \ # The ASB results obtained from betabinomial.R
-n ${name} \ # sample name
-f ${fasta} \ # Reference genome (.fa)
-m ${motif_file} \ # The .meme file downloaded from HOCOMOCO/JASPAR/CISBP
-o ${output_dir}
```

Chromosome.R: Draw the  distribution of ASB SNPs across 23 chromosomes with ggplot2. If users intend to directly use the script to draw the distribution, please refer to the following examples:

```
R --slave --no-restore --file=${basepath}/chromosome.R --args ${output_dir}/html_summary/${name}.chrom.info.txt \ # A file containing chromosomes and their respective lengths
${output_dir}/inpeak/${name}_counts_AS_0.1.txt \ # The ASB results obtained from betabinomial.R
${output_dir}/html_summary/${name}_chrom_dist.pdf # Output Fig
```

convertID.R: Perform the transformation of Ensembl gene IDs to gene symbols (for an improved display).

draw_cCRE.py: Generate pie charts of the distribution of cCREs (candidate cis-regulatory elements) SNPs enriched in.

draw_cCREs_hist.R: Generate histograms depicting the distribution of cCREs, specifically highlighting the enrichment of SNPs.

draw_annotation.R: Generate visual representations of genomic regions annotated by snpEFF tool.

draw_motif.R: Construct histograms illustrating the enrichment of motifs.

draw_AR_score.R: Create scatter plots depicting the correlation between the reference allele ratio and PWM score change of SNPs.

draw_freq.R: Generate graphical depictions showcasing the motif information content and the positional frequency of motif disruptions caused by SNPs.

best_match.py: Identify the best Fimo hits.

motif_disrupt.py: Get the positions of motifs disrupted by SNPs and calculate the frequency of disrupted position occurrences.

motif_score.py: Determine the positional weight matrix (PWM) scores for both reference and alternative alleles.

trans_alt.py: Convert the input FASTA file containing the reference alleles  to a FASTA file with alternative alleles.

data/: A directory containing human cCREs obtained from the SCREEN database.

example_data/: A directory containing sample names for testing the pipeline, encompassing both single-end and paired-end read data.

---

#### 5. Output Results

QC/: The quality control outcomes implemented by `fastp`.

find_intersecting_snps/; map/; remap/; filter_remapped_reads/:  Intermediate files generated by `WASP` (for removing mapping bias).

dup/: The results files after duplicate removal with Picard.

counts/: The read counts of ref and alt alleles performed by `WASP-bam2h5.py`. In this folder, the `${name}_counts.txt` file represents the subsequent read counts file for processing. 

inpeak/ : SNPs exclusively within ChIP-seq peaks. This folder contains three files: `${name}_counts.txt`, `${name}_counts_0.1.pdf`, and `${name}_counts_AS_0.1.txt`. Among them, `${name}_counts.txt` includes the read count information of SNPs in peaks; `${name}_counts_0.1.pdf` displays the  ref allele ratio distribution of ASB SNPs and non-ASB SNPs; `${name}_counts_AS_0.1.txt` encompasses the statistical information of all ASB SNPs and non-ASB SNPs. Examples of the last two files is provided below.

The distribution of ref allele ratio:

![${name}_counts_0.1.pdf](https://github.com/Liying1996/ASBanalyzer/raw/main/example_data/example_outputs/example_figs/example_counts_0.1.jpg)

The ASB SNPs:

![${name}_counts_0.1.txt](https://github.com/Liying1996/ASBanalyzer/raw/main/example_data/example_outputs/example_figs/example_counts_0.1_txt.jpg)

annotation/: Generate outputs and visualize the genomic regions as well as cCREs in which SNPs are enriched. This directory primarily includes five files: `${name}.ann.txt`,  `${name}.ann.distrubution.pdf`, `${name}.ccre.txt`, `${name}.cCREs.hist.pdf`, and `${name}.cCREs_pie.pdf`. These files consist of the annotation file obtained through snpEFF and the cCREs file in which the SNPs are enriched. Moreover, distribution charts are generated based on these two files.

Annotation:

![](https://github.com/Liying1996/ASBanalyzer/raw/main/example_data/example_outputs/example_figs/example.ann.distrubution.jpg)

Example of cCREs enrichment (Pie chart version):

![](https://github.com/Liying1996/ASBanalyzer/raw/main/example_data/example_outputs/example_figs/example_cCREs_pie.jpg)

Example of cCREs enrichment (Barplot chart version):

![](https://github.com/Liying1996/ASBanalyzer/raw/main/example_data/example_outputs/example_figs/example_cCREs_hist.jpg)

motif/: This directory contains the outcomes of motif enrichment, the associations involving positions of motif disruption and motif information content, as well as the correlations regarding the allele ratio and PWM score change of reference and alternative alleles. The results of motif analysis primarily includes the `${name}_AR_score.pdf`, `${name}_dist.pdf`, `${name}_disrupt_pos.pdf`, and `${name}_results_inmotif.txt` files. The former files provide visualizations, while the last offers statistical outcomes.. Additionally, this directory encompasses the GTEx_GWAS/ subfolder, within which results are presented for SNPs and their associations with corresponding GWAS phenotypes and GTEx eQTLs. Among them, the results are divided into SNPs disrupting motif recognition sequences and the outcomes for all SNPs associated with GWAS phenotypes and GTEx. 

Motif analysis results (${name}\_dist.pdf, ${name}\_AR_score.pdf,  ${name}_disrupt_pos.pdf): 

![motif analysis results](https://github.com/Liying1996/ASBanalyzer/raw/main/example_data/example_outputs/example_figs/example_motif_analysis.jpg)

Example GWAS Associations result: 

![example_gwas_associations](https://github.com/Liying1996/ASBanalyzer/raw/main/example_data/example_outputs/example_figs/example_gwas_associations.jpg)

Example GTEx genes result: 

![example_gtex_genepairs](https://github.com/Liying1996/ASBanalyzer/raw/main/example_data/example_outputs/example_figs/example_gtex_genepairs.jpg)

screenshots: Within this directory, an `${name}_screenshots_links.txt` file will be generated, containing screenshot links for each ASB SNP with a range of approximately ±20 base pairs corresponding to the Genome Browser and Variant Viewer. An example is illustrated in the figure below.

![example_screenshots_links](https://github.com/Liying1996/ASBanalyzer/raw/main/example_data/example_outputs/example_figs/example_screenshots_links.jpg)

html_summary/: This directory comprises a comprehensive HTML-based summary incorporating all the preceding results and graphical representations. Additionally, users can also exclusively review the `${name}_all_summary_ID.txt` file, which is a text-based summary excluding visualizations.

example_all_summary_ID:

![example_all_summary_ID.txt](https://github.com/Liying1996/ASBanalyzer/raw/main/example_data/example_outputs/example_figs/example_all_summary_ID.jpg)

The homepage of HTML-based summary:

![example_font.html](https://github.com/Liying1996/ASBanalyzer/raw/main/example_data/example_outputs/example_figs/example_font_html.jpg)

Upon clicking on a particular chromosome in the navigation bar of the page, the entirety of information regarding ASB SNPs on that chromosome will be presented on the right-hand side:

![](https://github.com/Liying1996/ASBanalyzer/raw/main/example_data/example_outputs/example_figs/example_font_chr.jpg)

The summary also includes the annotation and motif analysis results, an example (example_data/example_outputs/html_summary/ENCFF001HIA_font.html) is also presented for reference.  Users can download the [example output folder](https://github.com/Liying1996/ASBanalyzer/tree/main/example_data/example_outputs/html_summary) and simply open the ENCFF001HIA_font.html file to view the output results. (Please note that users must download all the files in this folder for the HTML to display correctly.)

Furthermore, the aforementioned example output files can all be found in the example_data/example_outputs/ directory.

