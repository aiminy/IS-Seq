import sys
import operator
import re
import os
import subprocess
import numpy as np
from itertools import groupby
import fnmatch

#This script take as input the "final_parse" file of each sample and put in a list all the files corresponding to the same sample that contain the reads id per intergration site (IS). As output it prints a list of IS (chr, pos) with the strand +/- considering the LTR direction registered and the number of reads per IS.
def unique(seq):
	# order preserving
	noDupes = []
	[noDupes.append(i) for i in seq if not noDupes.count(i)]
	return noDupes

#
myblast=open(sys.argv[1])
my_legend=open(sys.argv[2])
dwdt=sys.argv[3]
os.chdir(dwdt)

baseDir = os.path.dirname(sys.argv[1])
baseFile = os.path.basename(sys.argv[1])
outDir = os.path.join(baseDir,"DiversityOut")

if not os.path.exists(outDir):
    os.makedirs(outDir)

my_report=open(os.path.join(outDir,baseFile[:-4]+'_report.txt'),'w')
grouped_IS=open(os.path.join(outDir,baseFile[:-4]+'_grouped_IS.txt'),'w')
NonGrouped=open(os.path.join(outDir,baseFile[:-4]+'_NonGrouped.txt'),'w')

list_pos=[]
list_strand=[]
list_diff=[]
list_item1=[]
list_item3=[]
final_list=[]
mydict={}
myleg={}
myData={}
my_pos={}
dist_cutoff = 7
#List all the file containing the reads per IS
for f in os.listdir(dwdt):
	if fnmatch.fnmatch(f,'*LTR*LC*filter*_chr*.txt'):
		list_pos.append(f)
#Read the input file (final_parse) and store the infos about strand, chr and pos
for l in myblast.readlines():
	l= l.rstrip()
	if l[0]=='#':
		continue
	l= l.split('\t')
	tmp= l[4], l[5]
	try:
		mydict[tmp].append(l[0])
	except:
		mydict[tmp]=[l[0]]
#
for l in my_legend.readlines():
	l= l.rstrip()
	if l[0]=='#':
		continue
	l= l.split('\t')
	try:
		myleg[l[1]].append(l[0])
	except:
		myleg[l[1]]=[l[0]]
#Defines the corrispondence between the input file and the files into the list
sys.argv[1]=re.split('\.|\_', sys.argv[1])
print sys.argv[1][2]+'\t'+sys.argv[1][4]+'\t'+sys.argv[1][7]+'\n'
for element in list_pos:
	element=re.split('\.|\_', element)
#	print element[2]+'\n'
	#Finds the correct LTR and LC barcode
	if element[2]==sys.argv[1][2] and element[4]==sys.argv[1][4]:
		#Finds the correct filter: filterNo, filter30, filter45, filter60
		if element[5]==sys.argv[1][7]:	
			print element[2]+'\t'+element[4]+'\t'+element[5]+'__'+sys.argv[1][2]+'\t'+sys.argv[1][4]+'\t'+sys.argv[1][7]+'\n'
			for chr,pos in mydict:
				#if the IS falls into a standard chr
				if len(element)==9:
					if chr == element[6] and pos == element[7]:
						if str(set(mydict[chr,pos]))[6:-3] == 'R1_rev':
							for name in myleg:
								if element[6]==name:
									for code in myleg[name]:
										#print chr, pos, code, '-', str(len(mydict[chr,pos]))
										list_strand.append([chr, pos, code, '-', str(len(mydict[chr,pos]))])
						else:
							for name in myleg:
								if element[6]==name:
									for code in myleg[name]:
										#print chr, pos, code, '-', str(len(mydict[chr,pos]))
										list_strand.append([chr, pos, code, '+', str(len(mydict[chr,pos]))])
				#if the IS falls into a chr_something
				if len(element)==10:
					if chr == (element[6]+'_'+element[7]) and pos == element[8]:
						if str(set(mydict[chr,pos]))[6:-3] == 'R1_rev':
							for name in myleg:
								if (element[6]+'_'+element[7])==name:
									for code in myleg[name]:
										#print chr, pos, code, '-', str(len(mydict[chr,pos]))
										list_strand.append([chr, pos, code, '-', str(len(mydict[chr,pos]))])
						else:
							for name in myleg:
								if (element[6]+'_'+element[7])==name:
									for code in myleg[name]:
										#print chr, pos, code, '+', str(len(mydict[chr,pos]))
										list_strand.append([chr, pos, code, '+', str(len(mydict[chr,pos]))])
				#if the IS falls into a chr_something_something
				if len(element)==11:
					if chr == (element[6]+'_'+element[7]+'_'+element[8]) and pos == element[9]:
						if str(set(mydict[chr,pos]))[6:-3] == 'R1_rev':
							for name in myleg:
								if (element[6]+'_'+element[7]+'_'+element[8])==name:
									for code in myleg[name]:
										#print chr, pos, code, '-', str(len(mydict[chr,pos]))
										list_strand.append([chr, pos, code, '-', str(len(mydict[chr,pos]))])
						else:
							for name in myleg:
								if (element[6]+'_'+element[7]+'_'+element[8])==name:
									for code in myleg[name]:
										#print chr, pos, code, '+', str(len(mydict[chr,pos]))
										list_strand.append([chr, pos, code, '+', str(len(mydict[chr,pos]))])
#Sorts the list by chr and pos
for item in sorted(list_strand, key=lambda x: (int(x[2]), int(x[1]))):
	try:
		myData[item[2]].append(int(item[1]))
	except:
		myData[item[2]]=[int(item[1])]
#
for chr in myData:
   if len(myData[chr])>=2:
	   count=0
	   for count in range(0,(len(myData[chr])-1)):
			next_count=count+1
			if count==next_count:
				continue
			else:
				#compare the i element of this list with the next one
				myDiff= int(sorted(myData[chr])[count]) - int(sorted(myData[chr])[next_count])
				if abs(myDiff)<= dist_cutoff:
					#if the distance between two element of this list is lower than the cut_off, make another list
					list_diff.append([chr, myData[chr][count], myData[chr][next_count]])
#
for item in sorted(list_strand, key=lambda x: (int(x[2]), int(x[1]))):
	control= True
	for r in sorted(list_diff):
		if int(item[2])==int(r[0]) and int(item[1])==int(r[1]):
			list_item1.append([item[0], item[1], item[2], item[3], int(item[4])])
			control= False
		if int(item[2])==int(r[0]) and int(item[1])==int(r[2]):
			list_item1.append([item[0], item[1], item[2], item[3], int(item[4])])
			control= False
	if control:
		list_item3.append([item[0], item[1], item[2], item[3], int(item[4])])
#
unique_list=unique(list_item1)
for i in unique_list:
	tmp=i[0],i[1],i[2],i[3],i[4]
	try:
		my_pos[int(i[2])].append(tmp)
	except:
		my_pos[int(i[2])]=[tmp]
for a in my_pos:
	b=my_pos[a]
	#print d
	d=sorted(b, key=lambda x: (int(x[2]), int(x[1])))
	diff = [int(d[i+1][1])-int(d[i][1]) for i in range(len(d)-1)]
	avg = sum(diff) / len(diff)
	c=[[d[0][0]]]
	m=[[d[0][1]]]
	n=[[d[0][4]]]
	strand=[[d[0][3]]]
	for i in range(0,len(d)-1):
		if int(d[i+1][1])-int(d[i][1]) <= dist_cutoff and d[i+1][3]==d[i][3]:
			c[-1].append(d[i+1][0])
			m[-1].append(d[i+1][1])
			n[-1].append(d[i+1][4])
			strand[-1].append(d[i+1][3])
		else:
			c.append([d[i+1][0]])
			m.append([d[i+1][1]])
			n.append([d[i+1][4]])
			strand.append([d[i+1][3]])
	c=np.array(c,dtype=object)
	m=np.array(m,dtype=object)
	n=np.array(n,dtype=object)
	strand=np.array(strand,dtype=object)
	#Check if this script is trying to merge different integration sites together.
	#Opposite strand will not be part of the same integration site and you will get the warnining message.
	for segno in strand:
		for contatore in range(0,len(segno)-1):
			if segno[contatore+1]==segno[contatore]:
				continue
			else:
				print 'Warninig: The merged sites contain opposite strands.'
				print 'Please check the report and modify the grouped_IS file, if necessary.'
	print_array= np.dstack((c,m,n,strand))
	my_array= np.vstack(([c],[m],[n],[strand]))
#Print a report in which the reads of the same integration sites are grouped together
	for i in my_array[[2],:]:
		count=0
		for x in i:
			x= np.array(x,dtype=np.int64)
			reads_count= sum(x)
			l= x.argmax()
			final_list.append([str(my_array[0][count][l]), str(my_array[1][count][l]), str(d[0][2]), str(my_array[3][count][l]), str(reads_count)])
			if len(my_array[0][count])>=2:
				print '>'+str(my_array[0][count][l])+'_'+str(my_array[1][count][l])+'_'+str(d[0][2])+'_'+str(my_array[3][count][l])+'_'+str(reads_count)
				my_report.write('>'+str(my_array[0][count][l])+'_'+str(my_array[1][count][l])+'_'+str(d[0][2])+'_'+str(my_array[3][count][l])+'_'+str(reads_count)+'\n')
				for valori in range(len(my_array[0][count])):
					my_report.write(str(my_array[0][count][valori])+'\t'+str(my_array[1][count][valori])+'\t'+str(d[0][2])+'\t'+str(my_array[3][count][valori])+'\t'+str(my_array[2][count][valori])+'\n')
					print my_array[0][count][valori], my_array[1][count][valori], d[0][2], my_array[3][count][valori], my_array[2][count][valori]
			count+=1
#
for c in list_item3:
	NonGrouped.write(str(c[0])+'\t'+str(c[1])+'\t'+str(c[2])+'\t'+str(c[3])+'\t'+str(c[4])+'\n')
	final_list.append([str(c[0]), str(c[1]), str(c[2]), str(c[3]), str(c[4])])
for y in sorted(final_list, key=lambda x: (int(x[2]), int(x[1]))):
	grouped_IS.write('\t'.join(y)+'\n')
#
myblast.close()
my_legend.close()
my_report.close()
#
