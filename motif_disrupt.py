# Obtain SNPs' positions of motifs and calculate the frequency of pos of disruption
import argparse
from collections import defaultdict

parser = argparse.ArgumentParser()
parser.add_argument("-n", help = "a biosample's name, eg. ENCFF758RQJ", required = True)
parser.add_argument("-m", help = "Motif dir", required = True)
parser.add_argument("-l", help = "Motif Length", required = True)
parser.add_argument("-o", help = "output dir", required = True)
args = parser.parse_args()

name = args.n 
motif_dir = args.m
length = int(args.l)
output = args.o


dic = defaultdict(list)  # Pos of Motif disruption
with open(motif_dir + '/temp_files/' + name + '_thresh.txt') as f1:
	for line in f1:
		line = line.strip()
		pos = line.split('\t')[1]
		strand = line.split('\t')[4]
		start = line.split('\t')[2]
		end = line.split('\t')[3]
		if strand == '-':
			disrupt_pos = int(end) - 21
		else:
			disrupt_pos =  21 - int(start)
		dic[pos].append(line)
		dic[pos].append(disrupt_pos)

AS_inmotif_file = motif_dir + '/' + name + '_AS_inmotif.txt'
nonAS_inmotif_file = motif_dir + '/' + name + '_nonAS_inmotif.txt'


def write_file(file, flag, freqfile):
	freqfile = open(freqfile, 'a')
	pos_dict = {x:0 for x in range(1, length+1)} # 1-maxLen
	sum_count = 0
	
	with open(file) as f2:
		for line in f2:
			line = line.strip()
			chrom = line.split('\t')[0]
			pos1 = line.split('\t')[1]
			pos2 = line.split('\t')[2]
			pos = chrom + ':' + pos1 + '-' + pos2

			pos_dict[dic[pos][1]+1] += 1 # calculate frequency
			sum_count += 1

	for p in pos_dict:
		freqfile.write(flag +  '\t' + str(p) + '\t' + str(pos_dict[p]) + '\t' + str(sum_count) + '\t' + str(pos_dict[p]/sum_count) + '\n')


	freqfile.close()

write_file(AS_inmotif_file, 'AS', output + '/' + name + '_pos_freq.txt')
write_file(nonAS_inmotif_file, 'nonAS', output + '/' + name + '_pos_freq.txt')


