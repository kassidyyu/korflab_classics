# skewer.py by Kassidy

import argparse
import readfasta

parser = argparse.ArgumentParser(description='GC skew and composition')
parser.add_argument('file', type=str, help='FASTA file of genome sequence')
parser.add_argument('-w', '--window', type=int, default=100, 
	help='size of window to comput GC skew')
arg = parser.parse_args()

beddict = {}
reader = readfasta.faopen(arg.file)
while True:
	record = readfasta.faread(reader)
	if record == None: break
	defline, seq = record
	bedlist = []                             # initialize list of tuples
	g = seq[:arg.window].count('G')
	c = seq[:arg.window].count('C')
	if g + c == 0: skew = 0
	else:          skew = (g-c) / (g+c)
	bedlist.append((skew, (g+c)/arg.window)) # tuple with skew & composition
	for i in range(len(seq) - arg.window):
		if   seq[i] == 'G':
			g -= 1
		elif seq[i] == 'C':
			c -= 1
		if   seq[i + arg.window] == 'G':
			g += 1
		elif seq[i + arg.window] == 'C':
			c += 1
		if g + c == 0: skew = 0
		else:          skew = (g-c) / (g+c)
		bedlist.append((skew, (g+c) / arg.window))
	# use deflines as keys (might change) to save all calculations
	beddict[defline] = bedlist

# write BED file