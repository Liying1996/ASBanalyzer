while getopts 'n:o:' arg 
do
        case $arg in
             n)
                name=$OPTARG
                ;;
             o)
                path=$OPTARG
                ;;
             ?)  
                echo "Unknown argument!"
                exit 1
                ;;
        esac
done

mkdir ${path}/delta_score/

awk -F '\t' 'FNR==NR {x[$1":"$2"-"$3];next} ($2 in x)' ${path}/${name}_AS_snps.bed ${path}/${name}_best.txt >  ${path}/delta_score/${name}_AS_tmp.txt

awk -F '\t' 'FNR==NR {x[$1"\t"$2"\t"$3"\t"$4"\t"$5];next} ($1"\t"$2"\t"$3"\t"$4"\t"$5 in x)' ${path}/delta_score/${name}_AS_tmp.txt  ${path}/${name}_ref_out.txt | sort > ${path}/delta_score/${name}_AS_ref_score.txt

awk -F '\t' 'FNR==NR {x[$1"\t"$2"\t"$3"\t"$4"\t"$5];next} ($1"\t"$2"\t"$3"\t"$4"\t"$5 in x)' ${path}/delta_score/${name}_AS_tmp.txt  ${path}/${name}_alt_out.txt | sort > ${path}delta_score/${name}_AS_alt_score.txt

paste -d '\t' <(cut -f2-9 ${path}/delta_score/${name}_AS_ref_score.txt) <(cut -f6-8 ${path}/delta_score/${name}_AS_alt_score.txt) > ${path}/delta_score/${name}_AS_scores.txt # get ref and alt scores

rm ${path}/delta_score/${name}_AS_tmp*.txt

# nonAS
awk -F '\t' 'FNR==NR {x[$1":"$2"-"$3];next} ($2 in x)' ${path}/${name}_nonAS_snps.bed ${path}/${name}_best.txt >  ${path}/delta_score/${name}_nonAS_tmp.txt

awk -F '\t' 'FNR==NR {x[$1"\t"$2"\t"$3"\t"$4"\t"$5];next} ($1"\t"$2"\t"$3"\t"$4"\t"$5 in x)' ${path}/delta_score/${name}_nonAS_tmp.txt  ${path}/${name}_ref_out.txt | sort > ${path}/delta_score/${name}_nonAS_ref_score.txt
awk -F '\t' 'FNR==NR {x[$1"\t"$2"\t"$3"\t"$4"\t"$5];next} ($1"\t"$2"\t"$3"\t"$4"\t"$5 in x)' ${path}/delta_score/${name}_nonAS_tmp.txt  ${path}/${name}_alt_out.txt | sort > ${path}/delta_score/${name}_nonAS_alt_score.txt

paste -d '\t' <(cut -f2-9 ${path}/delta_score/${name}_nonAS_ref_score.txt) <(cut -f6-8 ${path}/delta_score/${name}_nonAS_alt_score.txt) > ${path}/delta_score/${name}_nonAS_scores.txt

rm ${path}/delta_score/${name}_nonAS_tmp*.txt


awk '{print $0"\tAS"}' ${path}/delta_score/${name}_AS_scores.txt > ${path}/delta_score/${name}_all_scores.txt

awk '{print $0"\tnonAS"}' ${path}/delta_score/${name}_nonAS_scores.txt >> ${path}/delta_score/${name}_all_scores.txt

awk -F '\t' 'FNR==NR {x[$1":"$2-21"-"$2+20];next} ($1 in x)' ${path}/../../inpeak/${name}_counts_AS_0.1.txt ${path}/delta_score/${name}_all_scores.txt | sort >  ${path}/delta_score/${name}_AR_tmp1.txt

awk -F '\t' 'FNR==NR {x[$1];next} ($1":"$2-21"-"$2+20 in x)' ${path}/delta_score/${name}_all_scores.txt ${path}/../../inpeak/${name}_counts_AS_0.1.txt  | sort >  ${path}/delta_score/${name}_AR_tmp2.txt


paste -d '\t' <(cut -f1-9 ${path}/delta_score/${name}_AR_tmp2.txt) <(cat ${path}/delta_score/${name}_AR_tmp1.txt) > ${path}/delta_score/${name}_AR_score.txt


