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
	comp = (g+c)/arg.window
	# tuple with skew & composition rounded
	bedlist.append((f'{skew:.4f}', f'{comp:.4f}')) 
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
		comp = (g+c)/arg.window
		bedlist.append((f'{skew:.4f}', f'{comp:.4f}'))
	# use deflines as keys (might change) to save all calculations
	beddict[defline] = bedlist

# write BED file - GC comp & skew in name column
with open('gcoutput.bed', 'w') as outfile:
	# title
	print('track name=gcCalcs description="GC skew and composition from FA"',
		file=outfile)
	# column titles
	print('chrom', 'start', 'end', 'name_gc_skew_comp', sep='\t', 
	   file=outfile) # s##_c##
	for defline, tlist in beddict.items():
		defline = defline.split()
		chrom = defline[0]
		for i in range(len(tlist)):
			print(chrom, i, i+arg.window, sep='\t', end='\t', file=outfile)
			sc = tlist[i]
			print('s', sc[0], '_c', sc[1], sep='', file=outfile)
outfile.close()