die () {
    echo "ERROR: $*. " >&2
    exit 1
}

opts=false
optp=false
opth=false 
optv=false 
while getopts 's:p:w:i:c:h:v:f:j:a:o:m:g:r:' arg
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
     w)
        wasp_path=$OPTARG 
        ;;
     i)
        bwa_index=$OPTARG
        ;;
     c)
		chrom_info=$OPTARG
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
	 f)
		file_peak=$OPTARG
		;;
	 j)
		snpeff_jar=$OPTARG
		;;
	 a)
		annotation=$OPTARG 
		;;
     o)
		output_dir=$OPTARG
		;;
	 m)
		motif_file=$OPTARG
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

basepath=$(cd `dirname $0`; pwd)

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
Rscript ${basepath}/betabinomial.R ${output_dir}/inpeak/${name}_counts.txt ${name} 0.1


# Annotation

if [ "$annotation" == "human" ];then
	ccre_file=${basepath}/data/cCREs/GRCh38-ccREs.bed
	snpeff=GRCh38.p7.RefSeq
elif [ "$annotation" == "mouse" ];then
	ccre_file=${basepath}/data/cCREs/mm10-ccREs.bed
	snpeff=GRCm38.86
else
	echo 'End!'
	exit 1
fi

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

python ${basepath}/draw_cCRE.py -file ${output_dir}/annotation/${name}_ccre.txt -name ${name} 

# snpEFF Annotation
cd ${output_dir}/annotation/
awk -F '\t' '{print $1"\t"$3"\t"$NF"\t"$4"\t"$5}' ${output_dir}/annotation/${name}_counts_AS_0.1.bed > ${output_dir}/annotation/${name}.txt 
java -jar ${snpeff_jar} ${snpeff} ${output_dir}/annotation/${name}.txt >  ${output_dir}/annotation/${name}.ann.txt 
rm ${output_dir}/annotation/${name}.txt 

# Motif analysis
sh ${basepath}/cal_motif.sh -i ${output_dir}/inpeak/${name}_counts_AS_0.1.txt -n ${name} -f ${fasta} -m ${motif_file} -o ${output_dir}

# compare with GWAS and GTEx (only human hg38)
if [ "$annotation" == "human" ];then 
	mkdir ${output_dir}/motif/GTEx_GWAS/
	awk -F '[\t+]' 'FNR==NR {x["chr"$12"_"$13];next} ($1"_"$4 in x)' ${basepath}/data/GWAS/gwas_all_associations.txt ${output_dir}/motif/${name}_AS_inmotif.txt > ${output_dir}/motif/GTEx_GWAS/${name}_AS_inmotif_gwas_snps.txt
	awk -F '[\t+]' 'FNR==NR {x[$1"_"$4];next} ("chr"$12"_"$13 in x)' ${output_dir}/motif/${name}_AS_inmotif.txt ${basepath}/data/GWAS/gwas_all_associations.txt > ${output_dir}/motif/GTEx_GWAS/${name}_AS_inmotif_gwas_disease.txt

	awk -F '[\t|_]' 'FNR==NR {x[$1"_"$2];next} ($1"_"$4 in x)' ${gtex_file} ${output_dir}/motif/${name}_AS_inmotif.txt > ${output_dir}/motif/GTEx_GWAS/${name}_AS_inmotif_gtex_snps.txt
	awk -F '[\t|_]' 'FNR==NR {x[$1"_"$4];next} ($1"_"$2 in x)' ${output_dir}/motif/${name}_AS_inmotif.txt ${gtex_file} > ${output_dir}/motif/GTEx_GWAS/${name}_AS_inmotif_gtex_genepairs.txt
fi

echo 'End!'

