import argparse
from collections import defaultdict

parser = argparse.ArgumentParser()
parser.add_argument("-n", help = "please give an output name, eg. ENCFF758RQJ", required = True)
parser.add_argument("-input_file", help = "fimo out results", required = True)
parser.add_argument("-output", help = "output dir", required = True)
args = parser.parse_args()

biosample = args.n 
input_file = args.input_file
output = args.output

dic = defaultdict(list)

with open(input_file) as f:
    for line in f:
        line = line.strip()
        pos = line.split('\t')[1]
        dic[pos].append(line)

new = open(output + '/' + biosample + '_best.txt', 'w')
for pos in dic:
    max_score = float(dic[pos][0].split('\t')[5])
    max_mess = dic[pos][0]
    for message in dic[pos]:
        score = float(message.split('\t')[5])
        if score > max_score:
            max_score = score
            max_mess = message
    new.write(max_mess + '\n')

new.close()
