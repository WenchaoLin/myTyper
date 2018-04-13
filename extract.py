import sys
from Bio import SeqIO
from Bio.Seq import Seq


dicF = {}
dicR = {}

for record in SeqIO.parse(sys.argv[1],'fasta'):
	try:
		dicF[str(record.seq[-16:].reverse_complement())] += 1
		#dicR[str(record.seq[-16:])] += 1
	except:
		dicF[str(record.seq[-16:].reverse_complement())] = 1
		#dicR[str(record.seq[-16:])] = 1

for keyf,valuef in sorted(dicF.items(),key = lambda item:item[1]):
	print("%s\t%s" % (keyf,valuef))

'''
for keyr,valuer in sorted(dicR.items(),key = lambda item:item[1]):
	print("%s\t%s" % (keyr,valuer))
'''