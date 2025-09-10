# dust.py by Kassidy

from collections import defaultdict
import argparse
import readfasta
import math

parser = argparse.ArgumentParser(
	description='low-complexity filter for nucleotide sequences')
parser.add_argument('file', type=str, help='FASTA file of nt sequence(s)')
parser.add_argument('-s', '--size', type=int, default=20, # used for calc
	help='window size [%(default)i]')
parser.add_argument('-e', '--entropy', type=float, default=1.4, 
	help='entropy threshold [%(default).3f]')
parser.add_argument('-l', '--lower', action='store_true', 
	help='soft mask (lowercase masking)')
parser.add_argument('-g', '--gff', action='store_true', 
	help='output low-complexity regions as GFF file')
arg = parser.parse_args()

reader = readfasta.faopen(arg.file)
masked = []
# initialize dictionary, use defline as keys, list of start as values
if arg.gff: gdict = defaultdict(list)

while True:
	record = readfasta.faread(reader)
	if record == None: break
	defline, seq = record
	output = []
	# start off 
	a = seq[:arg.size].count('A')
	c = seq[:arg.size].count('C')
	g = seq[:arg.size].count('G')
	t = seq[:arg.size].count('T')
	nts = [a, c, g, t]
	e_calc = 0
	for val in nts:
		if val == 0: continue
		e_calc -= (val/arg.size) * math.log2(val)
	if e_calc < arg.entropy:
		if arg.lower: output[:arg.size] = output[:arg.size].lower()
		else:         output[:arg.size] = ['N'] * arg.size
		if arg.gff: gdict[defline].append(1) # start is 1 for 1-base in gff
	else:
		for i in range(arg.size):
			output.append(seq[i])

	for i in range(len(seq) - arg.size):
		if   seq[i] == 'A': nts[0] -= 1
		elif seq[i] == 'C': nts[1] -= 1
		elif seq[i] == 'G': nts[2] -= 1
		elif seq[i] == 'T': nts[3] -= 1
		if   seq[i + arg.size] == 'A': nts[0] += 1
		elif seq[i + arg.size] == 'C': nts[1] += 1
		elif seq[i + arg.size] == 'G': nts[2] += 1
		elif seq[i + arg.size] == 'T': nts[3] += 1
		e_calc = 0 # reset entropy calculation
		for val in nts:
			if val == 0: continue
			e_calc -= (val/arg.size) * math.log2(val)
		# figure out whether to mask or not
		if e_calc >= arg.entropy:
			output.append(seq[i + arg.size])
		elif output[-1] in 'Nacgt': # accounts for either masking
			if arg.lower:
				output.append(seq[i + arg.size].lower())
			else:
				output.append('N')
			# only when gff is needed, account for logic of indexing
			if arg.gff: gdict[defline].append(i + 2)
		else:
			if arg.gff: gdict[defline].append(i + 2) # same as above
			if arg.lower: 
				output[i : i+arg.size] = output[i : i+arg.size].lower()
			else:
				output[i : i+arg.size] = ['N'] * arg.size
	masked.append((defline, output))

with open('dust_output.fa', 'w') as outfile:
	for tup in masked:
		print('>', tup[0], sep='', file=outfile)
		for i in range(0, len(tup[1]) -59, 60):
			print(''.join(tup[1][i:i + 60]), file=outfile)
		if len(tup[1]) % 60 != 0:
			print(''.join(tup[1][i+60 : i+60 + len(tup[1])%60]), file=outfile)
outfile.close()

if arg.gff:
	with open('low_complexity.gff', 'w') as gff:
		for defline, slist in gdict.items():
			for i in range(len(slist)):
				# gff format output
				print(defline, 'dust', 'low_entropy_region',     # col 1-3
				slist[i], slist[i]+arg.size, '.', '?', '.',      # col 4-8
				f'ID=region{i};entropy_threshold={arg.entropy}', # col 9
				sep='\t', file=gff)
gff.close()