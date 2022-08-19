# Main function to analyze ASB SNPs


die () {
    echo "ERROR: $*. " >&2
    exit 1
}

opts=false
optp=false
opth=false
optv=false

# optional
basepath=$(cd `dirname $0`; pwd)
snpeff_jar=${basepath}/data/snpEFF/snpEff/snpEff.jar
wasp_path=${basepath}/WASP/
bwa_index=${basepath}/data/hg38_index/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta  
chrom_info=${basepath}/data/GRCh38_EBV.chrom.sizes.tsv 
file_peak=default # optioanl defalt:calling by macs2
gtex_file=${basepath}/data/GTEx_Analysis_v8_eQTL/Whole_Blood.v8.signif_variant_gene_pairs.txt 
fasta=${basepath}/data/hg38_index/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta

while getopts 's:p:h:v:m:o:i::c::f::g::r::' arg
do
	case $arg in

	s)
		$optp && die "Cannot specify option -s after specifying option -p"
	    opts=true
	    fastq=$OPTARG
	    name=`echo $fastq | awk -F '/' '{print $NF}' | sed "s/.\(fq\|fastq\)[0-9]*\(\.gz\)*//g"`
	    ;;
	p)
		$opts && die "Cannot specify option -p after specifying option -s"
		optp=true
		fastq1=`echo $OPTARG|awk -F ',' '{print $1}'`
		fastq2=`echo $OPTARG|awk -F ',' '{print $2}'`
		name1=`echo $fastq1 | awk -F '/' '{print $NF}' | sed "s/.\(fq\|fastq\)[0-9]*\(\.gz\)*//g"`
		name2=`echo $fastq2 | awk -F '/' '{print $NF}' | sed "s/.\(fq\|fastq\)[0-9]*\(\.gz\)*//g"`
		if [ "$name1" != "$name2" ];then
			name=${name1}_${name2}
		fi
		;;
     h)
		$optv && die "Cannot specify option -h after specifying option -v"
	    opth=true
		h5_dir=$OPTARG
		;;
	 v)
		$opth && die "Cannot specify option -v after specifying option -h"
		optv=true
		vcf_file=$OPTARG
		;;
	 m)
		motif_file=$OPTARG
		;;
     o)
		output_dir=$OPTARG
		;;
     i)
        bwa_index=$OPTARG
        ;;
     c)
		chrom_info=$OPTARG
		;;

	 f)
		file_peak=$OPTARG
		;;

	 g)
		gtex_file=$OPTARG
		;;
	 r)
		fasta=$OPTARG
		;;
	\?) die "Invalid option!"
	  	;;
	esac
done




# Make HDF5 files
if $optv;then
	echo "Convert VCF to HDF5 files ..."
	mkdir ${output_dir}/h5/

	bgzip -c ${vcf_file} > ${vcf_file}.gz

	tabix -p vcf ${vcf_file}.gz

    for chrom in `seq 1 22` X;do
        tabix -h  ${vcf_file}.gz chr${chrom} > ${output_dir}/h5/chr${chrom}.vcf
    done

    gzip ${output_dir}/h5/chr*.vcf
	${wasp_path}/snp2h5/snp2h5 --chrom ${chrom_info} \
	       --format vcf \
	       --haplotype ${output_dir}/h5/haplotypes.h5 \
	       --snp_index ${output_dir}/h5/snp_index.h5 \
	       --snp_tab   ${output_dir}/h5/snp_tab.h5 \
	       ${output_dir}/h5/chr*.vcf.gz
	h5_dir=${output_dir}/h5/
fi

# QC; Mapping and Remove mapping bias

echo "QC and Remove mapping bias now ... "
mkdir ${output_dir}/QC/

if $opts;then
	echo "input single-end reads: $fastq"
	fastp -i ${fastq} -o ${output_dir}/QC/${name}.fq -j ${output_dir}/QC/${name}_fastp.json -h ${output_dir}/QC/${name}_fastp.html
	sh ${basepath}/remove_bias_SE.sh -w ${wasp_path} -i ${bwa_index} -f ${output_dir}/QC/${name}.fq -h ${h5_dir} -o ${output_dir}

else
	echo "input paired-end reads: $fastq1 $fastq2"
	fastp -i ${fastq1} -o ${output_dir}/QC/${name}.fq1 -I ${fastq2} -O ${output_dir}/QC/${name}.fq2 -j ${output_dir}/QC/${name}_fastp.json -h ${output_dir}/QC/${name}_fastp.html
	sh ${basepath}/remove_bias_PE.sh -w ${wasp_path} -i ${bwa_index} -1 ${output_dir}/QC/${name}.fq1 -2 ${output_dir}/QC/${name}.fq2 -h ${h5_dir} -o ${output_dir}

fi

# Get SNPs' reads
echo "Try to count SNPs ... "
mkdir ${output_dir}/counts/
python ${wasp_path}/CHT/bam2h5.py \
	  --chrom $chrom_info \
	  --data_type uint16 \
     --snp_index ${h5_dir}/snp_index.h5 \
     --snp_tab ${h5_dir}/snp_tab.h5 \
     --ref_as_counts ${output_dir}/counts/ref_as_counts.${name}.h5 \
     --alt_as_counts ${output_dir}/counts/alt_as_counts.${name}.h5 \
     --other_as_counts ${output_dir}/counts/other_as_counts.${name}.h5 \
     --read_counts ${output_dir}/counts/read_counts.${name}.h5 \
     --txt_counts ${output_dir}/counts/${name}_counts.txt \
	   ${output_dir}/dup/${name}_qual_sort.bam


# Peak calling (Optional)

if [ ${file_peak} == "default" ];then
	echo "Peak Calling using MACS2..."
	mkdir ${output_dir}/peak/
	macs2 callpeak -t ${output_dir}/map/${name}.bam \
	               -f BAM \
	               -n ${name} \
	               -g hs \
	               --outdir ${output_dir}/peak/
	file_peak=${output_dir}/peak/${name}_peaks.narrowPeak

fi

# Get SNPs in peaks
echo "Obtain SNPs in peaks ..."
mkdir ${output_dir}/inpeak/
awk -F ' ' '{print $1"\t"$2-1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}' ${output_dir}/counts/${name}_counts.txt > ${output_dir}/counts/${name}_counts.bed # 0-base
sed -i 's/ /\t/g' ${output_dir}/counts/${name}_counts.bed
bedtools intersect -a ${output_dir}/counts/${name}_counts.bed -b ${file_peak} | sort | uniq > ${output_dir}/counts/${name}_tmp.txt
cut -f1,3,4,5,6,7 ${output_dir}/counts/${name}_tmp.txt > ${output_dir}/inpeak/${name}_counts.txt
rm ${output_dir}/counts/${name}_tmp.txt ${output_dir}/counts/${name}_counts.bed


# Use binomial model to obtain allele specific sites
echo 'Use betabinomial model to find AS sites ...'
# Rscript ${basepath}/betabinomial.R ${output_dir}/inpeak/${name}_counts.txt ${name} 0.1
R --slave --no-restore --file=${basepath}/betabinomial.R --args ${output_dir}/inpeak/${name}_counts.txt ${name} 0.1


# Annotation

ccre_file=${basepath}/data/cCREs/GRCh38-cCREs.bed
snpeff=GRCh38.p7.RefSeq


# cCREs
mkdir ${output_dir}/annotation
awk -F '\t' '{print $1"\t"$2-1"\t"$2}' ${output_dir}/inpeak/${name}_counts_AS_0.1.txt >  ${output_dir}/annotation/tmp1 # 0-based
sed -i '1s/1$/Pos2/' ${output_dir}/annotation/tmp1
cut -f3-13 ${output_dir}/inpeak/${name}_counts_AS_0.1.txt >  ${output_dir}/annotation/tmp2
paste ${output_dir}/annotation/tmp1 ${output_dir}/annotation/tmp2 > ${output_dir}/annotation/${name}_counts_AS_0.1.bed
sed -i '1d'  ${output_dir}/annotation/${name}_counts_AS_0.1.bed
bedtools intersect -a ${output_dir}/annotation/${name}_counts_AS_0.1.bed -b $ccre_file -wa -wb >  ${output_dir}/annotation/${name}_ccre.txt

bedtools intersect -a ${output_dir}/annotation/${name}_counts_AS_0.1.bed -b $ccre_file -v >  ${output_dir}/annotation/tmp3
awk -F '\t' '{print $0"\t-\t-\t-\t-\t-\tUnclassified"}'  ${output_dir}/annotation/tmp3 >> ${output_dir}/annotation/${name}_ccre.txt

rm ${output_dir}/annotation/tmp*

python ${basepath}/draw_cCREs_pie.py -file ${output_dir}/annotation/${name}_ccre.txt -name ${name}

R --slave --no-restore --file=${basepath}/draw_cCREs_hist.R --args ${output_dir}/annotation/${name}_ccre.txt ${name}  ${output_dir}/annotation/${name}.cCREs.hist.pdf

# snpEFF Annotation
cd ${output_dir}/annotation/
awk -F '\t' '{print $1"\t"$3"\t"$NF"\t"$4"\t"$5}' ${output_dir}/annotation/${name}_counts_AS_0.1.bed > ${output_dir}/annotation/${name}.txt
java -jar ${snpeff_jar} ${snpeff} ${output_dir}/annotation/${name}.txt >  ${output_dir}/annotation/${name}.ann.txt
rm ${output_dir}/annotation/${name}.txt

# Annotation plots

grep -v '##' ${output_dir}/annotation/${name}.ann.txt | cut -d '|' -f1-2 > ${output_dir}/annotation/${name}.ann.tmp.txt

R --slave --no-restore --file=${basepath}/draw_annotation.R --args ${output_dir}/annotation/${name}.ann.tmp.txt ${name} ${output_dir}/annotation/${name}.ann.distrubution.pdf
rm ${output_dir}/annotation/${name}.ann.tmp.txt

# Motif analysis

sh ${basepath}/cal_motif.sh -i ${output_dir}/inpeak/${name}_counts_AS_0.1.txt -n ${name} -f ${fasta} -m ${motif_file} -o ${output_dir}

if [ -f "${output_dir}/motif/temp_files/${name}_pos_freq.txt" ]; then
	rm ${output_dir}/motif/temp_files/${name}_pos_freq.txt
fi

nrow1=`grep -n 'letter-probability' ${motif_file} | cut -d ':' -f1`
nrow1=`echo $nrow1 | awk '{print $1+1}'`
nrow2=`grep -n "URL" ${motif_file} | cut -d ':' -f1`
nrow2=`echo $nrow2 | awk '{print $1-1}'`
sed -n "${nrow1}, ${nrow2}p" ${motif_file} | sed 's/^ //g' | sed 's/  /\t/g' > ${output_dir}/motif/temp_files/${name}.pfm.txt

len=`wc -l ${output_dir}/motif/temp_files/${name}.pfm.txt | cut -d ' ' -f1`

python ${basepath}/motif_disrupt.py -n ${name} -m ${output_dir}/motif/ -l ${len} -o ${output_dir}/motif/temp_files/


R --slave --no-restore --file=${basepath}/draw_freq.R --args ${output_dir}/motif/temp_files/${name}_pos_freq.txt ${output_dir}/motif/temp_files/${name}.pfm.txt ${output_dir}/motif/${name}_disrupt_pos.pdf ${name}



python ${basepath}/motif_score.py -n ${name} -f ${output_dir}/inpeak/${name}_counts_AS_0.1.txt -p ${output_dir}/motif/temp_files/${name}.pfm.txt -m ${output_dir}/motif/

R --slave --no-restore --file=${basepath}/draw_AR_score.R --args ${output_dir}/motif/${name}_score.txt ${name} ${output_dir}/motif/${name}_AR_score.pdf


# compare with GWAS and GTEx (only human hg38)

mkdir ${output_dir}/motif/GTEx_GWAS/
awk -F '[\t+]' 'FNR==NR {x["chr"$12"_"$13];next} ($1"_"$2 in x)' ${basepath}/data/GWAS/gwas_all_associations.txt ${output_dir}/inpeak/${name}_counts_AS_0.1.txt > ${output_dir}/motif/GTEx_GWAS/${name}_all_gwas_snps.txt
awk -F '[\t+]' 'FNR==NR {x[$1"_"$2];next} ("chr"$12"_"$13 in x)' ${output_dir}/inpeak/${name}_counts_AS_0.1.txt ${basepath}/data/GWAS/gwas_all_associations.txt  > ${output_dir}/motif/GTEx_GWAS/${name}_all_gwas_associations.txt
awk -F '[\t|_]' 'FNR==NR {x[$1"_"$2];next} ($1"_"$2 in x)' ${output_dir}/inpeak/${name}_counts_AS_0.1.txt ${gtex_file} > ${output_dir}/motif/GTEx_GWAS/${name}_all_genepairs.txt
awk -F '[\t|_]' 'FNR==NR {x[$1"_"$2];next} ($1"_"$2 in x)' ${gtex_file} ${output_dir}/inpeak/${name}_counts_AS_0.1.txt  > ${output_dir}/motif/GTEx_GWAS/${name}_all_eqtl_snps.txt


awk -F '[\t+]' 'FNR==NR {x["chr"$12"_"$13];next} ($1"_"$4 in x)' ${basepath}/data/GWAS/gwas_all_associations.txt ${output_dir}/motif/${name}_AS_inmotif.txt > ${output_dir}/motif/GTEx_GWAS/${name}_AS_inmotif_gwas_snps.txt
awk -F '[\t+]' 'FNR==NR {x[$1"_"$4];next} ("chr"$12"_"$13 in x)' ${output_dir}/motif/${name}_AS_inmotif.txt ${basepath}/data/GWAS/gwas_all_associations.txt > ${output_dir}/motif/GTEx_GWAS/${name}_AS_inmotif_gwas_disease.txt

awk -F '[\t|_]' 'FNR==NR {x[$1"_"$2];next} ($1"_"$4 in x)' ${gtex_file} ${output_dir}/motif/${name}_AS_inmotif.txt > ${output_dir}/motif/GTEx_GWAS/${name}_AS_inmotif_gtex_snps.txt
awk -F '[\t|_]' 'FNR==NR {x[$1"_"$4];next} ($1"_"$2 in x)' ${output_dir}/motif/${name}_AS_inmotif.txt ${gtex_file} > ${output_dir}/motif/GTEx_GWAS/${name}_AS_inmotif_gtex_genepairs.txt


# Screen shots
mkdir ${output_dir}/screenshots/
if [ -f "${output_dir}/genomebrowser_screenshot/${name}_screenshots_links.txt" ]; then
	rm ${output_dir}/genomebrowser_screenshot/${name}_screenshots_links.txt
fi

for line in `awk -F '\t' '$NF=="significant"' ${output_dir}/inpeak/${name}_counts_AS_0.1.txt | sed 's/\t/+/g'`;do
	chrom=`echo $line | cut -d '+' -f1`
	pos=`echo $line | cut -d '+' -f2` # 1-based
	pos1=`awk "BEGIN {print ${pos}-20}"`
	pos2=`awk "BEGIN {print ${pos}+20}"`
	# pdf_page='http://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr1%3A1477833%2D1477873&hgsid=1280864071_rnlawaLfUQQ1i4quzzQsWWEDlzt1&hgt.psOutput=on'
	pdf_page="http://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=${chrom}%3A${pos1}%2D${pos2}&hgsid=1280864071_rnlawaLfUQQ1i4quzzQsWWEDlzt1&hgt.psOutput=on"

	pdf_url=$(curl -s "$pdf_page" | grep "the current browser graphic in PDF" | grep -E -o "\".+\"" | tr -d "\"" | sed 's/../https:\/\/genome.ucsc.edu/')
	# curl -s -o ${output_dir}/screenshots/${name}_${chrom}_${pos}.pdf ${pdf_url} -k
	chr_num=`echo ${chrom} | sed -s 's/chr//g'`
	variantviewer_link="http://www.ncbi.nlm.nih.gov/variation/view/?assm=GCF_000001405.28&chr=${chr_num}&from=${pos1}&to=${pos2}&mk=${pos}|SNP"
	echo ${line}"+"${pdf_page}"+"${variantviewer_link} >> ${output_dir}/screenshots/${name}_screenshots_links.txt

done

sed -i 's/+/\t/g' ${output_dir}/screenshots/${name}_screenshots_links.txt


# html summary
echo 'Generating html summary...'

mkdir ${output_dir}/html_summary/

cp ${basepath}/html_sample/font.html ${output_dir}/html_summary/${name}_font.html
cp ${basepath}/html_sample/style.css ${output_dir}/html_summary/${name}_style.css
cp ${basepath}/html_sample/func.js ${output_dir}/html_summary/${name}_func.js
cp ${output_dir}/annotation/${name}*pdf ${output_dir}/html_summary/
cp ${output_dir}/motif/${name}*pdf ${output_dir}/html_summary/
cp ${output_dir}/screenshots/${name}_screenshots_links.txt ${output_dir}/html_summary/
head -n23 ${chrom_info} > ${output_dir}/html_summary/${name}.chrom.info.txt
R --slave --no-restore --file=${basepath}/chromosome.R --args ${output_dir}/html_summary/${name}.chrom.info.txt ${output_dir}/inpeak/${name}_counts_AS_0.1.txt ${output_dir}/html_summary/${name}_chrom_dist.pdf


sed -i "s/ENCFF001HIA/${name}/g" ${output_dir}/html_summary/${name}_font.html

sed -i "s/func\.js/${name}_func\.js/g"  ${output_dir}/html_summary/${name}_font.html

sed -i "s/style\.css/${name}_style\.css/g"  ${output_dir}/html_summary/${name}_font.html


# rsID
annovar=${basepath}/data/annovar/annotate_variation.pl
awk -F '\t' '{print $1"\t"$2"\t"$2"\t"$3"\t"$4}' ${output_dir}/html_summary/${name}_screenshots_links.txt > ${output_dir}/html_summary/${name}_as_snps.txt
${annovar} ${output_dir}/html_summary/${name}_as_snps.txt ${basepath}/data/annovar/humandb/ -filter -build hg38 -dbtype avsnp150

mv ${output_dir}/html_summary/${name}_as_snps.txt.hg38_avsnp150_dropped ${output_dir}/html_summary/${name}_rsID.txt


python ${basepath}/summary.py -n ${name} -dir ${output_dir}/ -o ${output_dir}/html_summary/
sed -i 's/ /_/g' ${output_dir}/html_summary/${name}_all_summary.txt
R --slave --no-restore --file=${basepath}/convertID.R --args ${output_dir}/html_summary/${name}_all_summary.txt ${output_dir}/html_summary/${name}_all_summary_ID.txt
echo $(sed "$ ! s/$/\\\n/" ${output_dir}/html_summary/${name}_all_summary.txt) | sed 's/ chr/chr/g' | sed '1s/^/var fileText="/g' | sed 's/$/"/g' | sed 's/_rnlawaLfUQQ1i4quzzQsWWEDlzt1&hgt\.psOutput=on//g' > ${output_dir}/html_summary/tmp
cat ${output_dir}/html_summary/tmp >> ${output_dir}/html_summary/${name}_func.js
rm ${output_dir}/html_summary/tmp

echo 'End!'

