# readfasta.py by Kassidy

import sys
import gzip

filename = sys.argv[1]

def faopen(fasta):
	if   fasta == '-':          return sys.stdin
	elif fasta.endswith('.gz'): return gzip.open(fasta, 'rt')
	else:                       return open(fasta)

def faread(reader):
	defline = ''
	seq = []
	while True:                  # loop breaks with return
		start = reader.tell()	
		line = reader.readline().strip()
		if not line and not seq: # if line & seq are empty, file is done
			reader.close()       # close the file
			return None
		elif not line:           # line empty means last record in file
			return defline, ''.join(seq)
		elif defline and line.startswith('>'):
			# already found a defline, off to a new seq's defline 
			reader.seek(start)   # go back to beginning of line
			return defline, ''.join(seq)
		else:
			if line.startswith('>'):
				defline = line[1:]
			else:
				seq.append(line)

reader = faopen(filename)
while True:
	record = faread(reader)
	if record == None: break
	defline, seq = record