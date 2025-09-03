# randomseq.py by Kassidy

import argparse
import random

# specify number and length of sequences
parser = argparse.ArgumentParser(description='random sequence generator')
parser.add_argument('type', type=str, choices=['nt', 'aa'],
	help='nucleotide (nt) or amino acid (aa)')
parser.add_argument('num', type=int, help='number of sequences')
parser.add_argument('length', type=int, help='sequence length')
parser.add_argument('-l', '--line', type=int, default=80, 
	help='line length [%(default)i]')
parser.add_argument('-f', '--freq', type=str, default='0.25,0.25,0.25,0.25',
	help='frequencies of ACGT respectively, comma separated [%(default)s]')
parser.add_argument('-q', '--fastq', action='store_true', 
	help='output FASTQ with quality values')
	# user quality values with FASTQ? otherwise random?
parser.add_argument('-p', '--pfreq', type=str, default='ecoli', 
	choices=['ecoli', 'flat5', 'celegans'], 
	help='protein frequencies [%(default)s]')
arg = parser.parse_args()

# amino acid frequency list
aa_list = ['F', 'L', 'Y', 'H', 'Q', 'I', 'M', 'N', 'K', 'V', 
		   'D', 'E', 'S', 'C', 'W', 'P', 'R', 'T', 'A', 'G']
# based on GenScript Codon Usage Frequency Table(chart) Tool, freq/thousand
ecoli_freq = [38.1, 102, 29.7, 21.8, 43, 60.3, 26.4, 42, 47.7, 70.1,
              51.9, 57.8, 62, 11.3, 13.9, 42.4, 55.1, 55.3, 92.6, 73.4]
celegans_freq = [47.9, 86.3, 31.3, 23.3, 41.4, 60.9, 25.9, 48.7, 63.9, 61.8,
                 52.6, 65, 80.4, 20.3, 11, 48.9, 52.6, 58.1, 62.7, 53.7]

if arg.type == 'nt':
	probs = []
	for val in arg.freq.split(','):
		probs.append(float(val))
	for i in range(arg.num): # generate this many sequences
		seq = random.choices(['A', 'C', 'G', 'T'], weights=probs, k=int(arg.length))
		if arg.fastq:
			print('@', end='')
		else:
			print('>', end='')
		# title of sequence
		print('randomntseq', i+1, 'of', arg.num, sep='')
		if len(seq) <= arg.line:
			print(''.join(seq))
		else:
			for j in range(len(seq) // arg.line):
				print(''.join(seq[j:j+arg.line]))
			if len(seq) % arg.line != 0: # leftover seq
				print(''.join(seq[j*arg.line + arg.line:])) 
		if arg.fastq:
			print('+')
			print('insert quality values')
			# need to implement user quality values
else:
	for i in range(arg.num):
		# generate random sequences, weights based on chosen proteome / flat5
		if   arg.pfreq == 'flat5': # no weights, all equal
			seq = random.choices(aa_list, k=int(arg.length))
		elif arg.pfreq == 'ecoli':
			seq = random.choices(aa_list, weights=ecoli_freq, k=int(arg.length))
		else:
			seq = random.choices(aa_list, weights=celegans_freq, k=int(arg.length))
		print('>randomaaseq', i+1, 'of', arg.num, '-', arg.pfreq, sep='')
		if len(seq) <= arg.line:
			print(''.join(seq))
		else:
			for j in range(len(seq) // arg.line):
				print(''.join(seq[j:j+arg.line]))
			if len(seq) % arg.line != 0: # leftover seq
				print(''.join(seq[j*arg.line + arg.line:]))