import argparse
from collections import defaultdict

# python summary.py -n ${name} -dir ${output_dir} -o ${output_dir}/html_summary/

parser = argparse.ArgumentParser()
parser.add_argument("-n", help = "a biosample's name, eg. ENCFF000OBU", required = True)
parser.add_argument("-dir", help = "dir of results", required = True)
parser.add_argument("-o", help = "output dir", required = True)

args = parser.parse_args()

name = args.n
path = args.dir
output = args.o

d_ccre = {}
d_snpeff = {}
d_gtex = defaultdict(list)
d_gwas = defaultdict(list)
d_rsid = defaultdict(list)
list_motif = []
ccre_file = '{}/annotation/{}_ccre.txt'.format(path, name)
snf_file = '{}/annotation/{}.ann.txt'.format(path, name)
gtex_file = '{}/motif/GTEx_GWAS/{}_all_genepairs.txt'.format(path, name)
gwas_file = '{}/motif/GTEx_GWAS/{}_all_gwas_associations.txt'.format(path, name)
rsid_file =  '{}/html_summary/{}_rsID.txt'.format(path, name)
motif_file = '{}/motif/{}_AS_inmotif.txt'.format(path, name)

with open(ccre_file) as f:
	for line in f:
		line = line.strip()
		pos = line.split('\t')[0] + '\t' + line.split('\t')[2]
		d_ccre[pos] = line.split('\t')[-1] + '\t' + line.split('\t')[-2]

with open(snf_file) as f:
	for line in f:
		if line[0] != '#':
			line = line.strip()
			pos = '\t'.join(line.split('\t')[:2])
			d_snpeff[pos] = '_'.join(line.split('\t')[7].split('|')[1].split('_')[:-1])

with open(gtex_file) as f:
	for line in f:
		line = line.strip()
		pos = '\t'.join(line.split('_')[:2])
		d_gtex[pos].append(line.split('\t')[1])

with open(gwas_file) as f:
	for line in f:
		line = line.strip()
		pos = 'chr' + '\t'.join(line.split('\t')[11:13])
		d_gwas[pos].append(line.split('\t')[7].replace(' ', '_'))

with open(rsid_file) as f:
	for line in f:
		line = line.strip()
		pos = '\t'.join(line.split('\t')[2:4])
		d_rsid[pos].append(line.split('\t')[1])

with open(motif_file) as f:
	for line in f:
		line = line.strip()
		pos = line.split('\t')[0] + '\t' + line.split('\t')[3]
		list_motif.append(pos)


links_file =  '{}/html_summary/{}_screenshots_links.txt'.format(path, name)
output_file = open('{}/{}_all_summary.txt'.format(output, name), 'w')
with open(links_file) as f0:
	for line in f0:
		line = line.strip()
		pos = '\t'.join(line.split('\t')[:2])
		line1 = '\t'.join(line.split('\t')[2:-2])
		line2 = '\t'.join(line.split('\t')[-2:])
		
		if pos in d_ccre and d_ccre[pos] != 'Unclassified\t-':
			ccre = d_ccre[pos]
		else:
			ccre = '-\t-'
		
		if pos in d_snpeff:
			anno = d_snpeff[pos]
		else:
			anno = '-'

		if pos in d_gtex:
			gtex = ','.join(d_gtex[pos])
		else:
			gtex = '-'

		if pos in d_gwas:
			gwas =  ','.join(d_gwas[pos])
		else:
			gwas = '-'
		if pos in d_rsid:
			rsid = ','.join(d_rsid[pos])
		else:
			rsid = '-'
		if pos in list_motif:
			motif = 'YES'
		else:
			motif = '-'

		output_file.write(pos + '\t' + rsid + '\t' + line1 + '\t' + anno + '\t' + ccre + '\t' + motif + '\t' + gtex + '\t' + gwas + '\t' + line2 + '\n')
output_file.close()


