import os
import sys
import re
#This script take as input the *_final_parse.txt of each LTR-LC combination and generates multiple files: one for every single IS (chr and pos), which will contain every read ID corresponding to the IS written in the file name.

baseDir = os.path.dirname(sys.argv[1])
baseFile = os.path.basename(sys.argv[1])
outDir = os.path.join(baseDir,"ISout")

#oytDir = outDir = os.path.join(baseDir,"ISoutTest")

if not os.path.exists(outDir):
    os.makedirs(outDir)

myfile=open(sys.argv[1])
mydict={}
for l in myfile:
	l=l.rstrip()
	l=l.split('\t')
	tmp= l[4],l[5]
	try:
		mydict[tmp].append(l[1])
	except:
		mydict[tmp]=[l[1]]
#
for chr, pos in mydict:
	print chr, pos, len(mydict[chr, pos])
	out_f=os.path.join(outDir,baseFile[:-25]+'_'+baseFile[-12:-4]+'_'+str(chr)+'_'+str(pos)+'.txt')
	my_out=open(out_f, 'w')
	n=0
	for read in mydict[chr, pos]:
		if n<=len(mydict[chr, pos]):
			my_out.write(str(read)+'\n')
			n+=1
		else:
			my_out.close()
#
myfile.close()

