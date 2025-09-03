# seqstats.py by Kassidy

import argparse
import readfasta

parser = argparse.ArgumentParser(description='FASTA file statistics')
parser.add_argument('file', type=str, help='FASTA file')
parser.add_argument('-c', '--codon', action='store_true', 
	help='report codon frequencies')
arg = parser.parse_args()

seq_lens = []
nt_dict = {}
if arg.codon:
	codon_ct = {}
	codon = []

# read through file and keep track of counts
reader = readfasta.faopen(arg.file)
while True:
	record = readfasta.faread(reader)
	if record == None: break
	defline, seq = record
	for nt in seq:
		if nt not in nt_dict.keys(): nt_dict[nt] = 0
		nt_dict[nt] += 1
		if not arg.codon: continue
		# rest of this block only for counting codons
		if len(codon) < 3:
			codon.append(nt)
		else:
			codonstr = ''.join(codon)
			if codonstr not in codon_ct.keys(): 
				codon_ct[codonstr] = 0
			codon_ct[codonstr] += 1
			# includes reverse comp codons
			comp = str.maketrans('ACGTRYMKWSBDHVacgtrymkwsbdhv',
				'TGCAYRKMWSVHDBtgcayrkmwsvhdb') # copied from mcb185.py
			anti = codonstr.translate(comp)[::-1]
			if anti not in codon_ct.keys(): codon_ct[anti] = 0
			codon_ct[anti] += 1
			codon = codon[1:]
			codon.append(nt)
	seq_lens.append(len(seq))

# calculate statistics
seq_lens.sort(reverse=True)
num = len(seq_lens) # total number of sequences
max = seq_lens[0]
min = seq_lens[-1]
sum = sum(seq_lens)
mean = sum / num

# median calculation depends on even or odd number of sequences
if num % 2 == 1:
	median = seq_lens[num // 2]
else:
	median = (seq_lens[int(num / 2)] + seq_lens[int(num/2 - 1)]) / 2

running_sum = 0
for i in range(len(seq_lens)):
	if running_sum < sum / 2:
		running_sum += seq_lens[i]
	else: # use previous index
		n50 = seq_lens[i-1]
		break

nt_freq = {}
for nt, count in nt_dict.items():
	nt_freq[nt] = count / sum # take counts of each nucleotide over total letters

# similar to calculating nucleotide frequencies
if arg.codon:
	codon_freq = {}
	tot_codons = 0
	for ct in codon_ct.values():
		tot_codons += ct
	for codon, count in codon_ct.items():
		codon_freq[codon] = count / tot_codons

# report the statistics
print('There are', num, 'sequences with a combined total of', sum, 'letters')
print('The minimum sequence length is', min, 'and the maximum is', max)
print(f'The mean is {mean:.2f}')
print('The median is', median)
print('The N50 is', n50)

print('These are the relative frequencies of each letter')
for nt, freq in nt_freq.items():
	print(nt, f'{freq:.4f}')

if arg.codon:
	print('These are the relative frequencies of each codon')
	for codon, freq in codon_freq.items():
		print(codon, f'{freq:.4f}')