import argparse
from collections import defaultdict

parser = argparse.ArgumentParser()
parser.add_argument("-n", help = "a biosample's name, eg. ENCFF758RQJ", required = True)
parser.add_argument("-f", help = "ASB results file", required = True)
parser.add_argument("-p", help = "PFM (after cleaning)", required = True)
parser.add_argument("-m", help = "Motif dir", required = True)

args = parser.parse_args()

name = args.n 
ASB_file = args.f
ppm_file = args.p
motif_dir = args.m

asb_dict = {}
with open(ASB_file, 'r') as f0:
	for line in f0:
		line = line.strip()
		asb_dict['_'.join(line.split('\t')[:2])] = line


# PPM
ppm_list = []
with open(ppm_file) as f1:
	for line in f1:
		line = line.strip()
		ppm_list.append(line.split('\t'))

# A C G T
correspond_dict = {'A':0, 'C':1, 'G':2, 'T':3}

reverse_dict = {'A':'T', 'C':'G', 'G':'C', 'T':'A'}


motif_score_file = open(motif_dir + '/' + name + '_score.txt', 'w')
with open(motif_dir + '/temp_files/' + name + '_thresh.txt') as f2:
	for line in f2:
		line = line.strip()
		pos = line.split('\t')[1]
		strand = line.split('\t')[4]
		start = line.split('\t')[2]
		end = line.split('\t')[3]
		
		as_pos = int(pos.split('-')[-1]) - 20
		chrom = pos.split(':')[0]
		message = asb_dict[chrom + '_' + str(as_pos)]

		ref = message.split('\t')[2]
		alt = message.split('\t')[3]

		if strand == '-':
			disrupt_pos = int(end) - 21
			### Reverse Complement
			ref = reverse_dict[ref]
			alt = reverse_dict[alt]


		else:
			disrupt_pos =  21 - int(start)


		ref_score = ppm_list[disrupt_pos][correspond_dict[ref]]
		alt_score = ppm_list[disrupt_pos][correspond_dict[alt]]

		motif_score_file.write(message + '\t' + ref_score + '\t' + alt_score + '\n')
motif_score_file.close()





