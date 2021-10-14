## python to get alt fasta
import argparse, os
parser = argparse.ArgumentParser()
parser.add_argument("-name", help = "output name, eg. ENCFF758RQJ", required = True)
parser.add_argument("-input_file", help = "please give an AS bed file's path", required = True)
parser.add_argument("-fa", help = "please give the path of fastq file", required = True)
parser.add_argument("-output", help = "please give the output dir", required = True)

args = parser.parse_args()
biosample = args.name 
input_file = args.input_file
fasta = args.fa
output = args.output

seq = []
name = []

with open(fasta) as f:
    for line in f:
        line = line.strip()
        if line[0] != '>':
            seq.append(line)
        else:
            name.append(line)

seq2 = []
n = 0

alt_fasta = open(output + '/' + biosample + '_alt.fasta', 'w')
with open(input_file) as f2:
    for line in f2:
        line = line.strip()
        ref = line.split('\t')[-2]
        alt = line.split('\t')[-1]
        if ref == seq[n][20]:
            seq2.append(seq[n][:20] + alt + seq[n][21:])
        n += 1

for i in range(len(seq)):
    alt_fasta.write(name[i] + '\n')
    alt_fasta.write(seq2[i] + '\n')

alt_fasta.close()
