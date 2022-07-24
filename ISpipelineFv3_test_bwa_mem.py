#! /usr/bin/python2

import sys, getopt
import os
import subprocess
import select
from subprocess import call
from threading import Thread
import numbers
import re
import fnmatch
import os.path
import numpy as np
from itertools import groupby
import shutil
from shutil import copyfile
from more_itertools import locate
from Bio import SeqIO
from collections import OrderedDict

tmpFolder = ""

def usage():
	print ('Usage: python '+sys.argv[0]+' -1 <r1.fq.gz> -2 <r2.fq.gz> -s <sampleName> -o <outputFolder> -t <suffix> -r <researchFile> -u <referenceDataDir> -p <utilsProgramDir> -a <analysisType> -c <previousRun> -q <repRegQual>')
	print 'Example:'
	print ('python '+sys.argv[0]+' -1 /home/ayan/Aimin/ispipe/data/PL0431_S1_L001_R1_001.fastq.gz -2 /home/ayan/Aimin/ispipe/data/PL0431_S1_L001_R2_001.fastq.gz -s POOL-UCL-CPL-Re -o /home/ayan/Aimin/UploadToEgnyte/Aimin/ISseqOutput -t Mar04 -r /home/ayan/Aimin/ispipe/sample_research/Association_pool_CPL_Nov19Fix.csv -u /home/ayan/Aimin/ispipe/utilsRefData -p /home/ayan/Aimin/ispipe/utils -a read -c nothing -q 30')

def unique(seq):
	# order preserving
	noDupes = []
	[noDupes.append(i) for i in seq if not noDupes.count(i)]
	return noDupes

def getInputDataFromSampleResearch(sampleResearch,sampleName):

	dataSOfSampleResearch=[]

	mylegend=open(sampleResearch)

	for l in mylegend.readlines():
		l= l.rstrip()
		l= l.split(',')
		dataOfSampleResearch=OrderedDict()
		if (l[9]==sampleName):
			dataOfSampleResearch['LAM_PCR_ID']=l[0]
			dataOfSampleResearch['PT']=l[1]
			dataOfSampleResearch['Transduction_ID']=l[2]
			dataOfSampleResearch['Source']=l[3]
			dataOfSampleResearch['Sample_Type']=l[4]
			dataOfSampleResearch['Research_Clinic']=l[5]
			dataOfSampleResearch['TimePoint']=l[6]
			dataOfSampleResearch['Library']=l[9]
			dataOfSampleResearch['Organism']=l[11]
			dataOfSampleResearch['VectorType']=l[12]
			dataOfSampleResearch['vectorBed']=l[13]
			dataOfSampleResearch['Linker_Cassette']=l[14]
			dataOfSampleResearch['Vector']=l[15]
			dataSOfSampleResearch.append(dataOfSampleResearch)

	for i in dataSOfSampleResearch:
		print(i)

	return dataSOfSampleResearch

def reNameFile(DiversityOut,outputDir,inputPattern,sampleResearch):

	####### Rename files ###################
	#/lustre2/scratch/dpellin/collision/Lenti_Human/filter60/MLD01

	#dwdcFilter60=dwdc+"filter60/"+sampleName+"/"
	dwdcFilter60=os.path.join(outputDir,inputPattern+"/db/")
	print "collision Folder "+inputPattern+':' +dwdcFilter60

	if not os.path.exists(dwdcFilter60):
		os.makedirs(dwdcFilter60)
		#print(l[9]+"\t"+l[7]+"\t"+l[8]+"\n")

	inputFiles = []
	for f in fnmatch.filter(os.listdir(DiversityOut),'*'+inputPattern+'_grouped_IS.txt'):
		inputFiles.append(os.path.join(DiversityOut,f))

	outputFiles = []
	for f in fnmatch.filter(os.listdir(dwdcFilter60),'*_grouped_IS'):
		outputFiles.append(os.path.join(dwdcFilter60,f))

	#print inputFiles

	#print outputFiles

	check = compareFilesTime(inputFiles,outputFiles)

	if check:
		my_file=[]
		for filename in os.listdir(DiversityOut):
			if fnmatch.fnmatch(filename,'*'+inputPattern+'_grouped_IS.txt'):
				filename = os.path.join(DiversityOut, filename)
				my_file.append(filename)
				print(filename+"\n")

		#Open a legend were the real name are stored
		mylegend=open(sampleResearch)

		for l in mylegend.readlines():
			l= l.rstrip()
			if l[0]=='#':
				continue
			l= l.split(',')
			for a in my_file:
				sourceFile=a
				temp = os.path.basename(a)
				a= temp
				a=a.split('_')
				#print(a[0]+"\t"+l[9]+"\n"+a[1]+"\t"+l[7]+"\n"+a[2]+"\t"+l[8]+"\n")
				if ((a[0]==l[9]) and (a[1]==l[7])  and (a[2]==l[8])):
					destFile=l[9]+'''_'''+l[1]+'''_'''+l[2]+'''_'''+l[3]+'''_'''+l[4]+'''_'''+l[6]+'''_grouped_IS'''
					copyfile(sourceFile,dwdcFilter60+destFile)

def sumToTable(utilsDir,dwdcFilterNo,suffix,PreviousGroupedISfolder):

	inputFiles = []
	for f in fnmatch.filter(os.listdir(dwdcFilterNo),'*_grouped_IS'):
		inputFiles.append(os.path.join(dwdcFilterNo,f))

	dwdcFilterNoSuffix=os.path.join(dwdcFilterNo,suffix)

	if not os.path.exists(dwdcFilterNoSuffix):
		os.makedirs(dwdcFilterNoSuffix)

	outputFiles = []
	for f in fnmatch.filter(os.listdir(dwdcFilterNoSuffix),'*'):
		outputFiles.append(os.path.join(dwdcFilterNoSuffix,f))

	check = compareFilesTime(inputFiles,outputFiles)

	if check:

		def call_script_R(wd,suffix):
			mycmd='''Rscript --vanilla ''' +os.path.join(utilsDir,'''collisionTable.R''')+''' '''+wd+''' '''+suffix+''' '''+PreviousGroupedISfolder
			print mycmd
			subprocess.call(mycmd,shell=True)

		t1 = Thread(target=call_script_R, args=(dwdcFilterNo,suffix))
		t1.start()
		t1.join()

def annotateISite(dwdcFilter60,sortedKnownGene,genomeSorted,NT,utilsRef,utilsDir,vectorBed,VectorMask,suffix):

		threads = []
		def call_script_filterMap(filename):
			mycmd='''bedtools closest -a '''+filename+ ''' -b  '''+sortedKnownGene+''' -D "ref" -g '''+genomeSorted+''' > '''+filename+'''_closest_knownGenesV19.txt'''
			print mycmd
			subprocess.call(mycmd,shell=True)

		for filename in os.listdir(dwdcFilter60):
			if fnmatch.fnmatch(filename,"*_BedFormat"):
				t1 = Thread(target=call_script_filterMap, args=(os.path.join(dwdcFilter60,filename),))
				threads.append(t1)

		for i in xrange(0, len(threads), NT):
			# Start all threads
			for x in threads[i:i+NT]:
				x.start()
	 		# Wait for all of them to finish
			for x in threads[i:i+NT]:
				x.join()

		if (vectorBed!=""):
			threads = []
			def call_script_filterMap(filename):
				mycmd='''bedtools intersect -v -a '''+filename+''' -b '''+os.path.join(utilsRef,vectorBed)+''' -wa > '''+filename[:-4]+'''_VectMask.txt'''
				print mycmd
				subprocess.call(mycmd,shell=True)

			my_file=[]
			for filename in os.listdir(dwdcFilter60):
				if fnmatch.fnmatch(filename,"*_closest_knownGenes*.txt"):
					t1 = Thread(target=call_script_filterMap, args=(os.path.join(dwdcFilter60,filename),))
					threads.append(t1)

			for i in xrange(0, len(threads), NT):
				# Start all threads
				for x in threads[i:i+NT]:
					x.start()
		 		# Wait for all of them to finish
				for x in threads[i:i+NT]:
					x.join()

		if (VectorMask==""):
			for filename in os.listdir(dwdcFilter60):
				if fnmatch.fnmatch(filename,"*_closest_knownGenes*.txt"):
					shutil.copy(os.path.join(dwdcFilter60,filename),os.path.join(dwdcFilter60,filename[:-4]+'_VectMask.txt'))


			threads = []
			def call_script_filterMap(filename):
				mycmd='''cut -f 1,2,4,5,6,10,11,12,13 '''+filename+''' | sort -u > '''+filename[:-4]+'''_uniq.txt'''
				print mycmd
				subprocess.call(mycmd,shell=True)

			for filename in os.listdir(dwdcFilter60):
				if fnmatch.fnmatch(filename,"*_closest_knownGenes*_VectMask.txt"):
					t1 = Thread(target=call_script_filterMap, args=(os.path.join(dwdcFilter60,filename),))
					threads.append(t1)

			for i in xrange(0, len(threads), NT):
				# Start all threads
				for x in threads[i:i+NT]:
					x.start()
	 			# Wait for all of them to finish
				for x in threads[i:i+NT]:
					x.join()

		def call_script_R(wd,suffix):
			mycmd='''Rscript --vanilla ''' +os.path.join(utilsDir,'''ISAnnotation.R''')+''' '''+wd+''' '''+suffix
			print mycmd
			subprocess.call(mycmd,shell=True)

		t1 = Thread(target=call_script_R, args=(dwdcFilter60,suffix))
		t1.start()
		t1.join()

def main():

	if len(sys.argv) <= 1:
		usage()
		sys.exit(1)

	repRegQual=30

	try:
		opts, args = getopt.getopt(sys.argv[1:],"1:2:s:o:t:r:u:p:a:c:q:h")
	except getopt.GetoptError as err:
		print(err)
		usage()
		sys.exit(2)

	for o, a in opts:
		if o == '-h':
			usage()
			sys.exit()
		elif o == '-1':
			r1=a
		elif o == '-2':
			r2=a
		elif o == '-s':
			sampleName=str(a)
		elif o == '-o':
			outputFolder=str(a)
		elif o == '-t':
			suffix=str(a)
		elif o == '-r':
			researchFile=str(a)
		elif o == '-a':
		 	analysisType=str(a)
		elif o == '-u':
			utilsRefData=str(a)
		elif o == '-p':
			utils=str(a)
		elif o == '-c':
			PreviousCollision=str(a)
		elif o == '-q':
			repRegQual=int(a)

	#if not 'tmpFolder' in globals():
	#	tmpFolder=sampleName

	NT=3
	#repRegQual=30
	#repRegQual=5
	#repRegQual=0

	print repRegQual

	baseFolder=os.path.dirname(sys.argv[0])
	print baseFolder

	inputFolder= os.path.dirname(r1)
	outputFolder= os.path.join(outputFolder,suffix)

	if not os.path.exists(outputFolder):
		os.makedirs(outputFolder)

	if os.path.isdir(os.path.dirname(researchFile)):
		sampleResearch=researchFile
	else:
		sampleResearch=os.path.join(baseFolder,"sample_research",researchFile)

	print sampleResearch

	DataFromSampleResearch = getInputDataFromSampleResearch(sampleResearch,sampleName)

	print DataFromSampleResearch[-1]['Organism']

	Organism = DataFromSampleResearch[-1]['Organism']
	VectorType = DataFromSampleResearch[-1]['VectorType']
	Linker_Cassette = DataFromSampleResearch[-1]['Linker_Cassette']
	Library=DataFromSampleResearch[-1]['Library']
	vectorBed=DataFromSampleResearch[-1]['vectorBed']
	Vector=DataFromSampleResearch[-1]['Vector']

	if (Organism=="Human" and VectorType=="SIN-LV"):
		collisionFolder="Lenti_Human"
	if (Organism=="hg38" and VectorType=="SIN-LV"):
		collisionFolder="Lenti_Human"
	if (Organism=="Mouse" and VectorType=="SIN-LV"):
		collisionFolder="Lenti_Mouse"
	if (Organism=="Human" and VectorType=="RV"):
		collisionFolder="Retro_Human"
	if (Organism=="Mouse" and VectorType=="RV"):
		collisionFolder="Retro_Mouse"

	if (VectorType=="SIN-LV"):
		LTRFileFa="LTR_lentiviral.fa"
	if (VectorType=="RV"):
		LTRFileFa="LTR_retrovirus.fa"

	print collisionFolder +' ' +LTRFileFa

	randBrc=0
	if (Linker_Cassette=="LCbrcd"):
		randBrc=1

	if (VectorType=="SIN-LV"):
		LCFileFa="LC_completo.fa"
	if (VectorType=="RV"):
		LCFileFa="LC_completo.fa"

	print collisionFolder +' ' +LTRFileFa + ' '+LCFileFa

	utilsDir=utils
	utilsRef=utilsRefData

	dwdc=os.path.join(outputFolder,"collision",collisionFolder)

	print utilsDir
	print dwdc


	if (Organism=="Human"):
		maskFile=os.path.join(utilsRef,'hg19',"repeatMaskerHg19BED")
		chrList=os.path.join(utilsRef,'hg19',"hg19_chr_list_num.txt")
		dirGenome= os.path.join(utilsRef,'hg19',"hg19chrOnly.fa")
		sortedKnownGene=os.path.join(utilsRef,'hg19','wgEncodeGencodeAttrsV19_genesKNOWN_sorted.bed')
		genomeSorted=os.path.join(utilsRef,'hg19','hg19_genome_sorted.txt')
	if (Organism=="hg38"):
		maskFile=os.path.join(utilsRef,'hg38','repeatMaskerhg38BED')
		chrList=os.path.join(utilsRef,'hg38','hg38.genome.num.txt')
		dirGenome= os.path.join(utilsRef,'hg38','hg38ChrOnly.fa')
		sortedKnownGene=os.path.join(utilsRef,'hg38','hg38_genesKNOWN_sorted.bed')
		genomeSorted=os.path.join(utilsRef,'hg38','hg38.genome.sorted.txt')
	if (Organism=="hg18"):
		maskFile=os.path.join(utilsRef,'hg18','repeatMaskerhg38BED')
		chrList=os.path.join(utilsRef,'hg18','hg18.genome.num.txt')
		dirGenome= os.path.join(utilsRef,'h18','hg18ChrOnly.fa')
		sortedKnownGene=os.path.join(utilsRef,'hg18','hg18_genesKNOWN_sorted.bed')
		genomeSorted=os.path.join(utilsRef,'hg18','hg18.genome.sorted.txt')
	if (Organism=="Mouse"):
		maskFile=os.path.join(utilsRef,'mm10',"repeatMaskerMM10BED")
		chrList=os.path.join(utilsRef,'mm10',"mm10_chr_list_num.txt")
		dirGenome=os.path.join(utilsRef,'mm10',"mm10ChrOnly.fa")
		sortedKnownGene=os.path.join(utilsRef,'mm10','mm10_genesKNOWN_sorted.bed')
		genomeSorted=os.path.join(utilsRef,'mm10','mm10_genome_sorted.txt')

	seqPlat="MiSeq"
	VectorMask=""

	dwdt=outputFolder

	print ("looking for file in:"+dwdt)
	print utilsDir
	print LCFileFa
	print LTRFileFa

	lcFa=open(os.path.join(utilsRef,LCFileFa))
	lcFaLines=lcFa.readlines()
	lcFaSeq=lcFaLines[1]
	lcFaSeq=lcFaSeq.rstrip()
	lcminMatch=round(len(lcFaSeq)*0.9)
	lcFaSeqLength=len(lcFaSeq)

	ltrFa=open(os.path.join(utilsRef,LTRFileFa))
	ltrFaLines=ltrFa.readlines()
	ltrFaSeq=ltrFaLines[1]
	ltrFaSeq=ltrFaSeq.rstrip()
	ltrminMatch=round(len(ltrFaSeq)*0.9)


	# Display input and output file name passed as the args
	print ("R1 file : %s and r2 file: %s" % (r1,r2) )
	print ("sampleName: %s  Library: %s  Organism: %s   collisionFolder: %s  LTRFileFa: %s    LCFileFa: %s    randBrc: %s   seqPlat: %s" % (sampleName,Library,Organism,collisionFolder,LTRFileFa,LCFileFa,randBrc,seqPlat))
	#print tmpFolder


	fqGzinput = inputFolder
	inputFiles = []
	for f in fnmatch.filter(os.listdir(fqGzinput),'*.fastq.gz'):
		inputFiles.append(os.path.join(fqGzinput,f))

	outputFiles = []
	for f in fnmatch.filter(os.listdir(dwdt),'*_MatchBlastLtrLc'):
		outputFiles.append(os.path.join(dwdt,f))

	check = compareFilesTime(inputFiles,outputFiles)

	print inputFiles
	print outputFiles
	print check
	dwd=inputFolder

	if check:
		fq2MatchBlastltrLc(r1,r2,dwd,dwdt,utilsRef,LTRFileFa,LCFileFa,ltrminMatch,lcminMatch,lcFaSeqLength,check,NT)


	#####DEMULTIPLEXING ####################
	inputFiles = []
	for f in fnmatch.filter(os.listdir(dwdt),'*_trim12nt_qcTrimmed_MatchBlastLtrLc'):
		inputFiles.append(os.path.join(dwdt,f))

	outputFiles = []
	dir_4_DEMULTIPLEXING = os.path.join(dwdt,"R1_fastq_trim12nt_qcTrimmed_MatchBlastLtrLc"+"DEMULTIPLEXING")
	if not os.path.exists(dir_4_DEMULTIPLEXING):
		os.mkdir(dir_4_DEMULTIPLEXING)
	else:
		for f in fnmatch.filter(os.listdir(dir_4_DEMULTIPLEXING),'*_Barcode_*'):
			outputFiles.append(os.path.join(dir_4_DEMULTIPLEXING,f))

	dir_4_DEMULTIPLEXING = os.path.join(dwdt,"R2_fastq_trim12nt_qcTrimmed_MatchBlastLtrLc"+"DEMULTIPLEXING")
	if not os.path.exists(dir_4_DEMULTIPLEXING):
		os.mkdir(dir_4_DEMULTIPLEXING)
	else:
		for f in fnmatch.filter(os.listdir(dir_4_DEMULTIPLEXING),'*_Barcode_*'):
			outputFiles.append(os.path.join(dir_4_DEMULTIPLEXING,f))

	#print inputFiles
	#print outputFiles
	check = compareFilesTime(inputFiles,outputFiles)
	print check

	if check:
		demultiplexingWithBarcode(VectorType,check,dwdt,utilsRef)
		trimwithCutAdapt(dwdt,utilsRef,LTRFileFa,LCFileFa)

	########Random Barcode removal start ! ################################
	inputFiles = []
	temp = os.path.join(dwdt,"R2_fastq_trim12nt_qcTrimmed_MatchBlastLtrLcDEMULTIPLEXINGTofq")
	for f in fnmatch.filter(os.listdir(temp),'*'):
		inputFiles.append(os.path.join(temp,f))

	outputFiles = []
	RandomBarcodRemovalOut=os.path.join(dwdt,"RandomBarcodRemovalOutPut4R2CuReRun")
	if not os.path.exists(RandomBarcodRemovalOut):
		os.makedirs(RandomBarcodRemovalOut)

	for f in fnmatch.filter(os.listdir(RandomBarcodRemovalOut),'*'):
		outputFiles.append(os.path.join(RandomBarcodRemovalOut,f))
	check = compareFilesTime(inputFiles,outputFiles)
	print check
	if check:
		if (randBrc==1):
			dwdtBrc=os.path.join(dwdt,"LCbrcdCu")
			if not os.path.exists(dwdtBrc):
				os.makedirs(dwdtBrc)

			os.chdir(dwdtBrc)
			inputDir = os.path.join(dwdt,"R2_fastq_trim12nt_qcTrimmed_MatchBlastLtrLcDEMULTIPLEXINGTofq")
			inputPattern = 'R2_fastq*.fq_trimwithCutAdapt'
			seqFile = os.path.join(utilsRef,"LC1brcd_ancora.fa")
			seqFile1 = os.path.join(utilsRef,"LC1brcd_last.fa")
			outputDir = dwdtBrc
			RandomBarcodRemoval(inputDir,inputPattern,outputDir,seqFile,seqFile1,RandomBarcodRemovalOut,NT)
			os.chdir(dwdt)

	if analysisType == "read":
		print "read based analysis start"
		R1Out = os.path.join(dwdt,"R1_fastq_trim12nt_qcTrimmed_MatchBlastLtrLcDEMULTIPLEXINGTofq")
		R2Out = os.path.join(dwdt,"RandomBarcodRemovalOutPut4R2CuReRun")
		outputDir = os.path.join(dwdt,"CutAdapt")
		if not os.path.exists(outputDir):
			os.makedirs(outputDir)
		getCollisionTable(seqPlat,R1Out,R2Out,outputDir,sampleResearch,utilsDir,utilsRef,chrList,dirGenome,sampleName,sortedKnownGene,genomeSorted,vectorBed,suffix,NT,maskFile,VectorMask,repRegQual,PreviousCollision)

	if analysisType == "align2Vector":
		print "align vector start"
		outputDir = os.path.join(dwdt,"vector")

		if not os.path.exists(outputDir):
			os.makedirs(outputDir)

		R1Out = os.path.join(dwdt,"R1_fastq_trim12nt_qcTrimmed_MatchBlastLtrLcDEMULTIPLEXINGTofq")
		R2Out = os.path.join(dwdt,"RandomBarcodRemovalOutPut4R2CuReRun")

		if Vector == 'MND-GFP':
			dirGenome= os.path.join(utilsRef,'vector','pCDY.MND.GFP-KAN-BB-Aldevron-5LTR-3SIN-LTR.fa')
		if Vector == 'GlucoCco':
			dirGenome= os.path.join(utilsRef,'vector','GlucoCco','pCDY.EFS.GlucoCco_KAN_BB_Aldevron_Azadeh_5_LTR_3_SIN_LTR_R_U5-3_SIN_LTR.fa')
		if Vector == 'IUPF-CTNS':
			dirGenome= os.path.join(utilsRef,'vector','IUPF-CTNS','IUVPF_CTNS_LTRtoLTR.fa')
		if Vector == 'pCDY-EFS-hGLAco':
			dirGenome= os.path.join(utilsRef,'vector','pCDY-EFS-hGLAco','pCDY-EFS-hGLAco_LTRtoLTR.fa')


		align2vector(seqPlat,R1Out,R2Out,outputDir,sampleResearch,dirGenome,sampleName)


	if analysisType == "missingIS":
		print "missingIS analysis start"
		outputDirMissIS = os.path.join(dwdt,"CutAdapt","missingIS")
		if not os.path.exists(outputDirMissIS):
			os.makedirs(outputDirMissIS)

		threads = []
		for filename in os.listdir(os.path.join(dwdt,"CutAdapt","BAMSorted")):
			if fnmatch.fnmatch(filename,'*_aligned_mem_sort_inMask.bam'):
				rDir = utilsDir
				rScript = "getMissIS.R"
				refFasta = dirGenome
				IdentityThreshold=95
				input = os.path.join(dwdt,"CutAdapt","BAMSorted",filename)
				x = sampleName+"_"+filename[14:-28]+"_final_parse_filterNo_grouped_IS.txt"
				output =  os.path.join(outputDirMissIS,x)

				t1 = Thread(target=callGetMissIS,args=(rDir,rScript,refFasta,IdentityThreshold,input,output))
				threads.append(t1)

		for i in xrange(0, len(threads), NT):
			# Start all threads
			for x in threads[i:i+NT]:
				x.start()
			# Wait for all of them to finish
			for x in threads[i:i+NT]:
				x.join()

		getMissIsTable(outputDirMissIS,sampleResearch,utilsDir,utilsRef,chrList,dirGenome,sampleName,sortedKnownGene,genomeSorted,vectorBed,suffix,NT,maskFile,VectorMask)

	if analysisType == "VectorCount":
		outputDirMissIS = os.path.join(dwdt,"CutAdapt","missingIS","vectorAlign")
		if not os.path.exists(outputDirMissIS):
			os.makedirs(outputDirMissIS)

		threads = []
		for filename in os.listdir(os.path.join(dwdt,"CutAdapt","BAMSorted")):
			if fnmatch.fnmatch(filename,'*_aligned_mem_sort_inMask.bam'):
				rDir = utilsDir
				rScript = "align2vector.R"

				if Vector == 'MND-GFP':
					refFasta= os.path.join(utilsRef,'vector','pCDY.MND.GFP-KAN-BB-Aldevron-5LTR-3SIN-LTR.fa')
				if Vector == 'GlucoCco':
					refFasta= os.path.join(utilsRef,'vector','GlucoCco','pCDY.EFS.GlucoCco_KAN_BB_Aldevron_Azadeh_5_LTR_3_SIN_LTR_R_U5-3_SIN_LTR.fa')
				if Vector == 'IUPF-CTNS':
					dirGenome= os.path.join(utilsRef,'vector','IUPF-CTNS','IUVPF_CTNS_LTRtoLTR.fa')
				if Vector == 'pCDY-EFS-hGLAco':
					dirGenome= os.path.join(utilsRef,'vector','pCDY-EFS-hGLAco','pCDY-EFS-hGLAco_LTRtoLTR.fa')

				IdentityThreshold=95
				input = os.path.join(dwdt,"CutAdapt","BAMSorted",filename)
				x = sampleName+"_"+filename[14:-28]+"_final_parse_filterNo_grouped_IS.txt"
				output =  os.path.join(outputDirMissIS,x)

				t1 = Thread(target=callGetMissIS,args=(rDir,rScript,refFasta,IdentityThreshold,input,output))
				threads.append(t1)

		for i in xrange(0, len(threads), NT):
			# Start all threads
			for x in threads[i:i+NT]:
				x.start()
			# Wait for all of them to finish
			for x in threads[i:i+NT]:
				x.join()

		if Vector == 'MND-GFP':
			chrList=os.path.join(utilsRef,'vector','vector_chr_list_num.txt')
			sortedKnownGene=os.path.join(utilsRef,'vector','pCDY.MND.GFP-KAN-BB-Aldevron-5LTR-3SIN-LTR_sorted.bed')
			genomeSorted=os.path.join(utilsRef,'vector','vector_genome_sorted.txt')

		if Vector == 'GlucoCco':
			chrList=os.path.join(utilsRef,'vector','GlucoCco','vector_chr_list_num.txt')
			sortedKnownGene=os.path.join(utilsRef,'vector','GlucoCco','pCDY.EFS.GlucoCco_KAN_BB_Aldevron_Azadeh_5_LTR-3_SIN_LTR_R_U5-3_SIN_LTR_sorted.bed')
			genomeSorted=os.path.join(utilsRef,'vector','GlucoCco','vector_genome_sorted.txt')

		if Vector == 'IUPF-CTNS':
			chrList=os.path.join(utilsRef,'vector','IUPF-CTNS','vector_chr_list_num.txt')
			sortedKnownGene=os.path.join(utilsRef,'vector','IUPF-CTNS','IUVPF_CTNS_LTRtoLTR_sorted.bed')
			genomeSorted=os.path.join(utilsRef,'vector','IUPF-CTNS','vector_genome_sorted.txt')

		if Vector == 'pCDY-EFS-hGLAco':
			chrList=os.path.join(utilsRef,'vector','pCDY-EFS-hGLAco','vector_chr_list_num.txt')
			sortedKnownGene=os.path.join(utilsRef,'vector','pCDY-EFS-hGLAco','pCDY-EFS-hGLAco_LTRtoLTR_sorted.bed')
			genomeSorted=os.path.join(utilsRef,'vector','pCDY-EFS-hGLAco','vector_genome_sorted.txt')

		getVectorIsTable(outputDirMissIS,sampleResearch,utilsDir,utilsRef,chrList,sampleName,sortedKnownGene,genomeSorted,vectorBed,suffix,NT,VectorMask)

	if analysisType == "umi":

		print "umi based analysis start"

		inputDir = os.path.join(dwdt,"LCbrcdCu")

		#getUmiOut= os.path.join(dwdt,"GroupBasedUmi")
		getUmiOut= os.path.join(dwdt,"UmiBased")
		if not os.path.exists(getUmiOut):
			os.makedirs(getUmiOut)

		outDir=os.path.join(getUmiOut,"DiversityOut")
		if not os.path.exists(outDir):
			os.makedirs(outDir)

		outClusterUmi=os.path.join(getUmiOut,"UmiCluster")
		if not os.path.exists(outClusterUmi):
			os.makedirs(outClusterUmi)

		inputFiles = []
		for f in fnmatch.filter(os.listdir(inputDir),'*.fq_trimwithCutAdapt_reads_withAncora'):
			inputFiles.append(os.path.join(inputDir,f))

		outputFiles = []
		for f in fnmatch.filter(os.listdir(outDir),'*.txt'):
			outputFiles.append(os.path.join(outDir,f))

		check = compareFilesTime(inputFiles,outputFiles)
		print check

		if check:
			threads = []
			for filename in os.listdir(inputDir):
				if fnmatch.fnmatch(filename,'*.fq_trimwithCutAdapt_reads_withAncora'):
					t1 = Thread(target=extractUmiOnly, args=(os.path.join(inputDir,filename),getUmiOut))
					threads.append(t1)

			for i in xrange(0, len(threads), NT):
				# Start all threads
				for x in threads[i:i+NT]:
					x.start()
	 			# Wait for all of them to finish
				for x in threads[i:i+NT]:
					x.join()

			umi=getReadNameUmi(getUmiOut)

			threads = []
			finalParseFileDir=os.path.join(dwdt,"CutAdapt")
			chrListFile=chrList
			ISoutDir=os.path.join(dwdt,"CutAdapt","ISout")

			for filename in os.listdir(finalParseFileDir):
				if fnmatch.fnmatch(filename,'*_final_parse_*.txt'):
					t1 = Thread(target=calculateDiversity, args=(os.path.join(finalParseFileDir,filename),umi,chrListFile,ISoutDir,outDir,analysisType))
					threads.append(t1)

			for i in xrange(0, len(threads), NT):
				# Start all threads
				for x in threads[i:i+NT]:
					x.start()
				# Wait for all of them to finish
				for x in threads[i:i+NT]:
					x.join()

		outputDir=os.path.join(getUmiOut,"collision",collisionFolder)

		if not os.path.exists(outputDir):
			os.makedirs(outputDir)

		outputCollFiles = []
		for f in fnmatch.filter(os.listdir(outputDir),'*'):
			outputCollFiles.append(os.path.join(outputDir,f))

		check = compareFilesTime(outputFiles,outputCollFiles)

		check=True
		print check

		if check:

			inputPattern="filter60"
			reNameFile(outDir,outputDir,inputPattern,sampleResearch)
			dwdcFilter=os.path.join(outputDir,inputPattern+"/db")

			if PreviousCollision == "nothing":
				PreviousGroupedISfolder="nothing"
			else:
				PreviousGroupedISfolder=os.path.join(PreviousCollision,"UmiBased","collision","Lenti_Human",inputPattern+"/db")

			sumToTable(utilsDir,dwdcFilter,suffix,PreviousGroupedISfolder)
			dwdcFilterSuffix=os.path.join(dwdcFilter,suffix)
			annotateISite(dwdcFilterSuffix,sortedKnownGene,genomeSorted,NT,utilsRef,utilsDir,vectorBed,VectorMask,suffix)

			inputPattern="filterNo"
			reNameFile(outDir,outputDir,inputPattern,sampleResearch)
			dwdcFilter=os.path.join(outputDir,inputPattern+"/db")

			if PreviousCollision == "nothing":
				PreviousGroupedISfolder="nothing"
			else:
				PreviousGroupedISfolder=os.path.join(PreviousCollision,"UmiBased","collision","Lenti_Human",inputPattern+"/db")

			sumToTable(utilsDir,dwdcFilter,suffix,PreviousGroupedISfolder)
 			dwdcFilterSuffix=os.path.join(dwdcFilter,suffix)
			annotateISite(dwdcFilterSuffix,sortedKnownGene,genomeSorted,NT,utilsRef,utilsDir,vectorBed,VectorMask,suffix)

		inputExtractUmiFiles = []
		for f in fnmatch.filter(os.listdir(getUmiOut),'*_umi_extract'):
			inputExtractUmiFiles.append(os.path.join(getUmiOut,f))

		outputUmiClusterFiles = []
		for f in fnmatch.filter(os.listdir(outClusterUmi),'*'):
			outputCollFiles.append(os.path.join(outClusterUmi,f))

		check = compareFilesTime(inputExtractUmiFiles,outputUmiClusterFiles)
		print check

		if check:

			threads = []
			for filename in os.listdir(getUmiOut):
				if fnmatch.fnmatch(filename,'*_umi_extract'):
					t1 = Thread(target=getUmi,args=(os.path.join(getUmiOut,filename),outClusterUmi))
					threads.append(t1)

			for i in xrange(0, len(threads), NT):
				# Start all threads
				for x in threads[i:i+NT]:
					x.start()
				# Wait for all of them to finish
				for x in threads[i:i+NT]:
					x.join()

			mycmd='''cat ''' +os.path.join(outClusterUmi,'''*_umi_extract_umi > ''') +os.path.join(outClusterUmi,'''R2_fastq_trim12nt_qcTrimmed_MatchBlastLtrLc_Barcode_allFiles.fq_trimwithCutAdapt_reads_withAncora_umi_extract_umi''')
			subprocess.call(mycmd,shell=True)

			threads = []
			for filename in os.listdir(outClusterUmi):
				if fnmatch.fnmatch(filename,'*_umi_extract_umi'):
					t1 = Thread(target=clusterUmi,args=(os.path.join(outClusterUmi,filename),))
					threads.append(t1)

			for i in xrange(0, len(threads), NT):
				# Start all threads
				for x in threads[i:i+NT]:
					x.start()
				# Wait for all of them to finish
				for x in threads[i:i+NT]:
					x.join()

			callR(utilsDir,"drawBarPlot4Umi.R",outClusterUmi,os.path.join(outClusterUmi,"UmiBarPlots"))

	if analysisType == "fragment":

		print "fragment based analysis start"

		fragmentOut= os.path.join(dwdt,"FragmentBased2")
		if not os.path.exists(fragmentOut):
			os.makedirs(fragmentOut)

		bamSortedOut=os.path.join(dwdt,"CutAdapt","BAMSorted")

		getFragmentLength(utilsDir,sampleName,bamSortedOut,fragmentOut)
		filter60(fragmentOut)
		getISsite(utilsDir,fragmentOut,NT)

		threads = []
		finalParseFileDir=fragmentOut
		chrListFile=chrList
		ISoutDir=os.path.join(fragmentOut,"ISout")

		outDir=os.path.join(fragmentOut,"DiversityOut")
		if not os.path.exists(outDir):
			os.makedirs(outDir)

		for filename in os.listdir(finalParseFileDir):
			if fnmatch.fnmatch(filename,'*_final_parse_*.txt'):
				t1 = Thread(target=getCount4FragmentBased, args=(os.path.join(finalParseFileDir,filename),chrListFile,ISoutDir,outDir,analysisType))
				threads.append(t1)

		for i in xrange(0, len(threads), NT):
			# Start all threads
			for x in threads[i:i+NT]:
				x.start()
			# Wait for all of them to finish
			for x in threads[i:i+NT]:
				x.join()

		outputFiles = []
		for f in fnmatch.filter(os.listdir(outDir),'*.txt'):
			outputFiles.append(os.path.join(outDir,f))

		outputDir=os.path.join(fragmentOut,"collision",collisionFolder)

		if not os.path.exists(outputDir):
			os.makedirs(outputDir)

		outputCollFiles = []
		for f in fnmatch.filter(os.listdir(outputDir),'*'):
			outputCollFiles.append(os.path.join(outputDir,f))

		check = compareFilesTime(outputFiles,outputCollFiles)
		check=True
		print check

		if check:

			inputPattern="filter60"
			reNameFile(outDir,outputDir,inputPattern,sampleResearch)
			dwdcFilter=os.path.join(outputDir,inputPattern+"/db")

			if PreviousCollision == "nothing":
				PreviousGroupedISfolder="nothing"
			else:
				PreviousGroupedISfolder=os.path.join(PreviousCollision,"FragmentBased","collision","Lenti_Human",inputPattern+"/db")

			sumToTable(utilsDir,dwdcFilter,suffix,PreviousGroupedISfolder)
			dwdcFilterSuffix=os.path.join(dwdcFilter,suffix)
			annotateISite(dwdcFilterSuffix,sortedKnownGene,genomeSorted,NT,utilsRef,utilsDir,vectorBed,VectorMask,suffix)

			inputPattern="filterNo"
			reNameFile(outDir,outputDir,inputPattern,sampleResearch)
			dwdcFilter=os.path.join(outputDir,inputPattern+"/db")

			if PreviousCollision == "nothing":
				PreviousGroupedISfolder="nothing"
			else:
				PreviousGroupedISfolder=os.path.join(PreviousCollision,"FragmentBased","collision","Lenti_Human",inputPattern+"/db")

			sumToTable(utilsDir,dwdcFilter,suffix,PreviousGroupedISfolder)
 			dwdcFilterSuffix=os.path.join(dwdcFilter,suffix)
			annotateISite(dwdcFilterSuffix,sortedKnownGene,genomeSorted,NT,utilsRef,utilsDir,vectorBed,VectorMask,suffix)

def filter60(outputDir):
	###### filter alignment 60 bp ###########
	for filename in os.listdir(outputDir):
		if fnmatch.fnmatch(filename,'*_final_parse_filterNo.txt'):
			myblast=open(os.path.join(outputDir,filename))
			myblastout=os.path.join(outputDir,filename[:-13]+"_filter60.txt")
			with open(myblastout, 'w') as outfile:
				for l in myblast.readlines():
					l= l.rstrip()
					if l[0]=='#':
						continue
					l= l.split('\t')
					if int(l[3])>=60:
						outfile.write('\t'.join(l))
						outfile.write('\n')

def getISsite(utilsDir,outputDir,NT):
	###### IS ###########
	threads = []
	def call_script_filterMap(filename):
		mycmd='''python2 '''+os.path.join(utilsDir,'''try.py ''')+filename
		print mycmd
		subprocess.call(mycmd,shell=True)

	for filename in os.listdir(outputDir):
		if fnmatch.fnmatch(filename,'*_final_parse_*.txt'):
			t1 = Thread(target=call_script_filterMap, args=(os.path.join(outputDir,filename),))
			threads.append(t1)

	for i in xrange(0, len(threads), NT):
		# Start all threads
		for x in threads[i:i+NT]:
			x.start()
	 	# Wait for all of them to finish
		for x in threads[i:i+NT]:
			x.join()

def getFragmentLength(utilsDir,sampleName,bamSortedOut,outputDir):

	for filename in os.listdir(bamSortedOut):
		if fnmatch.fnmatch(filename,'*_allFilter_rehead_exact3nt_nonSupplementary.bam'):
			inputFile = os.path.join(bamSortedOut,filename)
			mycmd='''python '''+os.path.join(utilsDir,'''try_pysam_with_fragment_length.py ''')+inputFile+ ''' '''+ sampleName +''' '''+outputDir
			subprocess.call(mycmd,shell=True)

def getCount4FragmentBased(finalParseFile,chrListFile,ISoutDir,outDir,analysisType):

	myblast=open(finalParseFile)
	my_legend=open(chrListFile)
	dwdt=ISoutDir

	os.chdir(dwdt)

	baseDir = os.path.dirname(finalParseFile)
	baseFile = os.path.basename(finalParseFile)
	outDir = outDir

	if not os.path.exists(outDir):
		os.makedirs(outDir)

	my_report=open(os.path.join(outDir,baseFile[:-4]+'_report.txt'),'w')
	NonGrouped=open(os.path.join(outDir,baseFile[:-4]+'_NonGrouped.txt'),'w')
	grouped_IS=open(os.path.join(outDir,baseFile[:-4]+'_grouped_IS.txt'),'w')

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
	mydictUmi={}

	#List all the file containing the reads per IS
	for f in os.listdir(dwdt):
		if fnmatch.fnmatch(f,'*LTR*LC*filter*_chr*.txt'):
			list_pos.append(os.path.join(dwdt,f))


	#Read the input file (final_parse) and store the infos about strand, chr and pos
	for l in myblast.readlines():
		l= l.rstrip()
		if l[0]=='#':
			continue
		l= l.split('\t')
		tmp= l[4], l[5]
		#print "tmp:"+str(tmp)
		try:
			mydict[tmp].append(l[0])
			mydictUmi[tmp].append(l[6])
		except:
			mydict[tmp]=[l[0]]
			mydictUmi[tmp]=[l[6]]

	#print "mydictUmi"
	for chr,pos in mydictUmi:
		#print chr,pos,mydictUmi[chr,pos]
		umiList=[]
		for x in mydictUmi[chr,pos]:
				umiList.append(str(x).strip("['']"))
		umiCount = len(umiList)
		umiUniqCount = len(set(umiList))
		mydictUmi[chr,pos]=umiUniqCount

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
	finalParseFileTmp=re.split('\.|\_', baseFile)
	print finalParseFileTmp[2]+'\t'+finalParseFileTmp[4]+'\t'+finalParseFileTmp[7]+'\n'

	for element in list_pos:
		element=re.split('\.|\_', element)
		#Finds the correct LTR and LC barcode
		if element[2]==finalParseFileTmp[2] and element[4]==finalParseFileTmp[4]:
			#Finds the correct filter: filterNo, filter30, filter45, filter60
			if element[5]==finalParseFileTmp[7]:
				print element[2]+'\t'+element[4]+'\t'+element[5]+'__'+finalParseFileTmp[2]+'\t'+finalParseFileTmp[4]+'\t'+finalParseFileTmp[7]+'\n'
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
											#print chr, pos, code, '-', str(len(mydict[chr,pos])
											list_strand.append([chr, pos, code, '-', str(len(mydict[chr,pos]))])
							else:
								for name in myleg:
									if (element[6]+'_'+element[7])==name:
										for code in myleg[name]:
											#print chr, pos, code, '+', str(len(mydict[chr,pos]))
											list_strand.append([chr, pos, code, '-', str(len(mydict[chr,pos]))])
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
											list_strand.append([chr, pos, code, '-',str(len(mydict[chr,pos]))])


	if analysisType == 'fragment':

		for x in list_strand:
			chr = x[0].strip("'")
			pos = x[1].strip("'")
			code = x[2]
			strand = x[3]
			for y in mydictUmi:
				if  y[0] == chr and y[1]==pos:
					umiCount=str(mydictUmi[y])
			x[4]=umiCount

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
		#print dcondition
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

		#print ("my_array:\n ", my_array)

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
						print str(my_array[0][count][valori])+'\t'+str(my_array[1][count][valori])+'\t'+str(d[0][2])+'\t'+str(my_array[3][count][valori])+'\t'+str(my_array[2][count][valori])
						print '\n'
				count+=1

	#print len(list_item3)
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

def callR(rDir,rScript,input,output):

	mycmd = '''Rscript ''' +os.path.join(rDir,rScript)+''' -i '''+input+''' -o '''+output
	subprocess.call(mycmd,shell=True)

def callGetMissIS(rDir,rScript,refFasta,IdentityThreshold,input,output):

	check = compareFileTime(input,output)

	if check:
		mycmd = '''Rscript ''' +os.path.join(rDir,rScript)+''' -i '''+input+''' -r '''+refFasta+''' -q '''+str(IdentityThreshold)+''' -o '''+output
		print mycmd
		subprocess.call(mycmd,shell=True)

def getUmi(file,outputDir):
	goodID=[]

	temp = os.path.basename(file)
	outputFile = os.path.join(outputDir,temp+ "_umi")

	check = compareFileTime(file,outputFile)

	if check:
		with open(file,'r') as f:
			lines=f.readlines()

		readName=[item[:-1] for item in lines[::4]] #get rid of '\n'

		with open(outputFile, "w") as umi:
			for listitem in readName:

				listitem=re.split(' ',listitem)
				rUmi = listitem[0][-18:]

				umi.write('%s\n' % rUmi)
		f.close()
		umi.close()

def getReadNameUmi(UmiFileDir):

	mydictUmi={}
	for filename in os.listdir(UmiFileDir):
		if fnmatch.fnmatch(filename,"*_umi_extract"):
			with open(os.path.join(UmiFileDir,filename),'r') as f:
				lines=f.readlines()
			f.close()
			readName=[item[:-1] for item in lines[::4]]
			for listitem in readName:
				listitem=re.split(' ',listitem)
				rName = listitem[0][1:-19]
				rUmi = listitem[0][-18:]
				try:
					mydictUmi[rName].append(rUmi)
				except:
					mydictUmi[rName]=[rUmi]

	return mydictUmi

def compareFileTime(inputFile,outputFile):

	if not os.path.exists(outputFile):
		check =True
	else:
		InputTime = os.path.getmtime(inputFile)
		OutTime = os.path.getmtime(outputFile)
		check = True if InputTime > OutTime else False
	return check

def compareFilesTime(inputFiles,outputFiles):

	InputTime  = []
	for x in inputFiles:
		InputTime.append(os.path.getmtime(x))

	OutTime  = []
	for x in outputFiles:
		OutTime.append(os.path.getmtime(x))

	a = np.array(InputTime)
	b = np.array(OutTime)

	check = []
	for x in a:
		for y in b:
			z = x > y
			check.append(z)

	check = all(check)

	return check

def fq2MatchBlastltrLc(r1,r2,dwd,dwdt,utilsRef,LTRFileFa,LCFileFa,ltrminMatch,lcminMatch,lcFaSeqLength,check,NT):

	#### UNZIP START ####################
	fq = os.path.join(dwdt,"fq")
	if not os.path.exists(fq):
		os.makedirs(fq)

	fqSplit = os.path.join(dwdt,"fqSplit")
	if not os.path.exists(fqSplit):
		os.makedirs(fqSplit)

	fqSplitSum = os.path.join(dwdt,"fqSplitSum")
	if not os.path.exists(fqSplitSum):
		os.makedirs(fqSplitSum)

	def call_script(FileRx,Rx):
		inputFile= os.path.join(dwd,FileRx)
		temp=os.path.join(dwdt,"fq")
		if not os.path.exists(temp):
			os.makedirs(temp)
		outputFile = os.path.join(temp,Rx+"_fastq")
		check = compareFileTime(inputFile,outputFile)
		mycmd='''gzip -cd '''+ inputFile+ ''' > ''' +outputFile
		if check:
			print mycmd
			subprocess.call(mycmd,shell=True)

	threads = []
	t1 = Thread(target=call_script, args=(r1,"R1"))
	t2 = Thread(target=call_script, args=(r2,"R2"))
	t1.start()
 	t2.start()
 	t1.join()
	t2.join()
	print ("Unzip ended!\n")

	###### split fa files for blast######################

	if check:
		def call_script(FileRx):
			inputFile= os.path.join(dwdt,"fq",FileRx)
			temp=os.path.join(dwdt,"fqSplit")
			if not os.path.exists(temp):
				os.makedirs(temp)
			outputFile = os.path.join(temp,FileRx+'''FqSplitaaaaaaaaaa''')
			check = compareFileTime(inputFile,outputFile)
			mycmd='''split -l10000000 -a10 ''' +inputFile+ ''' '''  +os.path.join(temp,FileRx+'''FqSplit''')
			if check:
				print mycmd
				subprocess.call(mycmd,shell=True)

		t1 = Thread(target=call_script, args=("R1_fastq",))
		t2 = Thread(target=call_script, args=("R2_fastq",))
		t1.start()
		t2.start()
		t1.join()
		t2.join()
		print ("Splitted Fq\n")

		# #### TRIM FIRST 12 NT START ####################
	if check:
		if not 'ftrimmer' in globals():
			ftrimmer=13
		if not 'Qtrimmer' in globals():
			Qtrimmer=33

		threads = []
		def call_script(FileRx):
			inputFile= FileRx
			outputFile = FileRx+'''_trim12nt'''
			check = compareFileTime(inputFile,outputFile)
			mycmd='''fastx_trimmer -f '''+str(ftrimmer)+''' -Q '''+str(Qtrimmer)+''' -i '''+inputFile+''' -o '''+outputFile
			if check:
				print mycmd
				subprocess.call(mycmd,shell=True)

		temp = os.path.join(dwdt,"fqSplit")
		for filename in os.listdir(temp):
			if fnmatch.fnmatch(filename,'*FqSplit??????????'):
				t1 = Thread(target=call_script, args=(os.path.join(temp,filename),))
				threads.append(t1)

		for i in xrange(0, len(threads), NT):
			# Start all threads
			for x in threads[i:i+NT]:
				x.start()
	 		# Wait for all of them to finish
			for x in threads[i:i+NT]:
				x.join()

		print ("Trimming ended!\n")

			##### TRIM BASED ON QUALITY START ####################
	if check:
		qualityTrimmer_t=20
		qualityTrimmer_l=100
		qualityTrimmer_Q=33

		threads = []
		def call_script(FileRx):
			inputFile= FileRx
			outputFile = FileRx+'''_qcTrimmed'''
			check = compareFileTime(inputFile,outputFile)

			mycmd='''fastq_quality_trimmer -t '''+str(qualityTrimmer_t)+ ''' -l '''+str(qualityTrimmer_l)+''' -Q '''+str(qualityTrimmer_Q)+''' -i '''+inputFile+''' -o '''+outputFile
			if check:
				print mycmd
				subprocess.call(mycmd,shell=True)

		temp = os.path.join(dwdt,"fqSplit")
		for filename in os.listdir(temp):
			if fnmatch.fnmatch(filename,'*_trim12nt'):
				t1 = Thread(target=call_script, args=(os.path.join(temp,filename),))
				threads.append(t1)

		for i in xrange(0, len(threads), NT):
			# Start all threads
			for x in threads[i:i+NT]:
				x.start()
	 		# Wait for all of them to finish
			for x in threads[i:i+NT]:
				x.join()

		print ("QC Trimming ended! (parameters  t= %s l=: %s ) \n" % (qualityTrimmer_t,qualityTrimmer_l))

		##### CONVERT FASTQ2FASTA ####################

		if check:
			threads = []
			def call_script(FileRx):
				inputFile= FileRx
				outputFile = FileRx+'''.fa'''
				check = compareFileTime(inputFile,outputFile)
				mycmd='''awk 'BEGIN{P=1}{if(P==1||P==2){gsub(/^[@]/,">");print}; if(P==4)P=0; P++}' ''' +inputFile+ ''' > '''  +outputFile
				if check:
					print mycmd
					subprocess.call(mycmd,shell=True)

			temp = os.path.join(dwdt,"fqSplit")
			for filename in os.listdir(temp):
				if fnmatch.fnmatch(filename,'*_qcTrimmed'):
					t1 = Thread(target=call_script, args=(os.path.join(temp,filename),))
					threads.append(t1)

			for i in xrange(0, len(threads), NT):
				# Start all threads
				for x in threads[i:i+NT]:
					x.start()
	 			# Wait for all of them to finish
				for x in threads[i:i+NT]:
					x.join()

			print ("Converted to Fa\n")



		##### LTR in R1 ####################
		if check:
			threads = []
			def call_script_filterMap(filename):

				inputFile= filename

				outputFile = filename+ '''_blast.txt'''

				check = compareFileTime(inputFile,outputFile)

				mycmd='''blastn -task blastn -word_size 7 -query ''' +inputFile+ '''  -subject '''+os.path.join(utilsRef,LTRFileFa)+''' -outfmt "6 std btop slen" -max_target_seqs 1 > '''  +outputFile

				if check:
					print mycmd
					subprocess.call(mycmd,shell=True)

			temp = os.path.join(dwdt,"fqSplit")
			for filename in os.listdir(temp):
				if fnmatch.fnmatch(filename,'R1_*_qcTrimmed.fa'):
					t1 = Thread(target=call_script_filterMap, args=(os.path.join(temp,filename),))
					threads.append(t1)

			for i in xrange(0, len(threads), NT):
				# Start all threads
				for x in threads[i:i+NT]:
					x.start()
				# Wait for all of them to finish
				for x in threads[i:i+NT]:
					x.join()
			print ("LTR sequence: %s searched on R1" % (LTRFileFa))

			##### LC in R2 ####################

			if check:
				threads = []
				def call_script_filterMap(filename):

					inputFile= filename

					outputFile = filename+ '''_blast.txt'''

					mycmd='''blastn -task blastn -word_size 7 -query ''' +inputFile+ '''  -subject '''+os.path.join(utilsRef,LCFileFa)+''' -outfmt "6 std btop slen" -max_target_seqs 1 > '''  +outputFile

					check = compareFileTime(inputFile,outputFile)

					if check:
						print mycmd
						subprocess.call(mycmd,shell=True)

				temp = os.path.join(dwdt,"fqSplit")
				for filename in os.listdir(temp):
					if fnmatch.fnmatch(filename,'R2_*_qcTrimmed.fa'):
						t1 = Thread(target=call_script_filterMap, args=(os.path.join(temp,filename),))
						threads.append(t1)

				for i in xrange(0, len(threads), NT):
					# Start all threads
					for x in threads[i:i+NT]:
						x.start()
	 				# Wait for all of them to finish
					for x in threads[i:i+NT]:
						x.join()
				print ("LC sequence: %s searched on R2" % (LCFileFa))

				#####Merge blast output#################################
				if check:
					def call_script(FileRx):
						temp= os.path.join(dwdt,"fqSplit")
						mycmd='''cat ''' +os.path.join(temp,FileRx[:-27])+'''*_blast.txt > '''  +os.path.join(fqSplitSum,FileRx+ '''_blast.txt''')
						print mycmd
						subprocess.call(mycmd,shell=True)

					t1 = Thread(target=call_script, args=("R1_fastq_trim12nt_qcTrimmed.fa",))
					t2 = Thread(target=call_script, args=("R2_fastq_trim12nt_qcTrimmed.fa",))
					t1.start()
					t2.start()
					t1.join()
					t2.join()
					print ("Merging blast outputs ")

				# ### FILTER LTR#####################

				if check:
					myblast=open(os.path.join(fqSplitSum,"R1_fastq_trim12nt_qcTrimmed.fa_blast.txt"))
					goodID=[]

					with open(os.path.join(fqSplitSum,"R1_fastq_trim12nt_qcTrimmed.fa_blast_filteredLTR"), "a") as myblastOut:
						for blastLine in myblast.readlines():
							blastLine= blastLine.rstrip()
							if blastLine[0]=='#':
								continue
							blastLine= blastLine.split('\t')
							if int(blastLine[8])<int(blastLine[9]):
								if int(blastLine[3])>=ltrminMatch and int(blastLine[9])>=ltrminMatch:
									try:
										if int(blastLine[12][-2]) in range(1,9) and int(blastLine[4])<5:
											myblastOut.write('\t'.join(blastLine))
											myblastOut.write('\n')
											goodID+=[blastLine[0]]
									except:
										if int(blastLine[12][-1])>=5 and int(blastLine[4])<5:
											myblastOut.write('\t'.join(blastLine))
											myblastOut.write('\n')
											goodID+=[blastLine[0]]
							else:
								if int(blastLine[3])>=ltrminMatch and int(blastLine[8])>=ltrminMatch:
									try:
										if int(blastLine[12][-2]) in range(1,9) and int(blastLine[4])<5:
											myblastOut.write('\t'.join(blastLine))
											myblastOut.write('\n')
											goodID+=[blastLine[0]]
									except:
										if int(blastLine[12][-1])>=5 and int(blastLine[4])<5:
											myblastOut.write('\t'.join(blastLine))
											myblastOut.write('\n')
											goodID+=[blastLine[0]]

					myblast.close()
					myblastOut.close()

					goodIDset = set(goodID)
					with open(os.path.join(fqSplitSum,"R1_fastq_trim12nt_qcTrimmed.fa_blast_filteredLTR_UniqID"), "w") as myblastOutID:
						myblastOutID.write('\n'.join(goodIDset))
					print ("Filter LTR blast ouput")

					### FILTER LC#####################

					myblast=open(os.path.join(fqSplitSum,"R2_fastq_trim12nt_qcTrimmed.fa_blast.txt"))
					goodID=[]
					with open(os.path.join(fqSplitSum,"R2_fastq_trim12nt_qcTrimmed.fa_blast_filteredLC"), "a") as myblastOut:
						for blastLine in myblast.readlines():
							blastLine= blastLine.rstrip()
							if blastLine[0]=='#':
								continue
							blastLine= blastLine.split('\t')
							if int(blastLine[4])<5:
								if int(blastLine[8])< int(blastLine[9]):
									if int(blastLine[9])>=lcminMatch and int(blastLine[3])==lcFaSeqLength:
										myblastOut.write('\t'.join(blastLine))
										myblastOut.write('\n')
										goodID+=[blastLine[0]]
								else:
									if int(blastLine[8])>=lcminMatch and int(blastLine[3])==lcFaSeqLength:
										myblastOut.write('\t'.join(blastLine))
										myblastOut.write('\n')
										goodID+=[blastLine[0]]
					myblast.close()
					myblastOut.close()

					goodIDset = set(goodID)
					with open(os.path.join(fqSplitSum,"R2_fastq_trim12nt_qcTrimmed.fa_blast_filteredLC_UniqID"), "w") as myblastOutID:
						myblastOutID.write('\n'.join(goodIDset))
					print ("Filter LC blast ouput")


					####### MERGE TRIMMED SEQUENCES ######################

					if check:
						def call_script(FileRx):
							temp= os.path.join(dwdt,"fqSplit")
							mycmd='''cat ''' +os.path.join(temp,FileRx+"FqSplit*_trim12nt_qcTrimmed") + ''' > '''  +os.path.join(fqSplitSum,FileRx+'''_trim12nt_qcTrimmed''')
							print mycmd
							subprocess.call(mycmd,shell=True)

						t1 = Thread(target=call_script, args=("R1_fastq",))
						t2 = Thread(target=call_script, args=("R2_fastq",))

						t1.start()
						t2.start()

						t1.join()
						t2.join()

					# ###### REMOVE TMP SPLITTED TRIMMED SEQUENCES ###################### NB. By means of blast algortihm we identify a select those seqs with a good match to LTR and LC seqs

					if check:
						def call_script(FileRx):
							temp = os.path.join(dwdt,"fqSplit")
							mycmd='''rm ''' +os.path.join(temp,FileRx+'''FqSplit*''')
							print mycmd
							subprocess.call(mycmd,shell=True)

						t1 = Thread(target=call_script, args=("R1_fastq",))
						t2 = Thread(target=call_script, args=("R2_fastq",))
						t1.start()
						t2.start()
						t1.join()
						t2.join()


					###SELECT SEQUENCES WITH MATCH IN BOTH LTR AND LC,By means of blast algortihm we identify a select those seqs with a good match to LTR and LC seqs

					if check:
						def call_script(FileRx,FileRxID):
							mycmd='''seqtk subseq ''' +os.path.join(fqSplitSum,FileRx)+ '''  ''' +os.path.join(fqSplitSum,FileRxID)+ ''' > '''  +os.path.join(dwdt,FileRx+ '''_MatchBlastLtrLc''')
							print mycmd
							inputFile= os.path.join(fqSplitSum,FileRx)
							outputFile = os.path.join(dwdt,FileRx+ '''_MatchBlastLtrLc''')
							check = compareFileTime(inputFile,outputFile)
							if check:
								subprocess.call(mycmd,shell=True)

						t1 = Thread(target=call_script, args=("R1_fastq_trim12nt_qcTrimmed","R1_fastq_trim12nt_qcTrimmed.fa_blast_filteredLTR_UniqID"))
						t2 = Thread(target=call_script, args=("R2_fastq_trim12nt_qcTrimmed","R2_fastq_trim12nt_qcTrimmed.fa_blast_filteredLC_UniqID"))

						t1.start()
						t2.start()

						t1.join()
						t2.join()

						print ("Select subseq matching LTR and LC (completo)")

def demultiplexingWithBarcode(VectorType,check,dwdt,utilsRef):
	if check:
		if (VectorType=="SIN-LV"):
			barcodeLTRseq="barcode_LTR.fil"
			barcodeLCseq="barcode_LC.fil"

		if (VectorType=="RV"):
			barcodeLTRseq="barcode_LTR_Retrovirus.fil"
			barcodeLCseq="barcode_LC_Retrovirus.fil"

		def call_script(FileRx,FileRxBarcode):
			dir_4_DEMULTIPLEXING = os.path.join(dwdt,FileRx+"DEMULTIPLEXING")
			check = True
			outputFiles = []
			if not os.path.exists(dir_4_DEMULTIPLEXING):
				os.mkdir(dir_4_DEMULTIPLEXING)
			else:
				for f in fnmatch.filter(os.listdir(dir_4_DEMULTIPLEXING),'*_Barcode_*'):
					outputFiles.append(os.path.join(dir_4_DEMULTIPLEXING,f))

			inputFiles = [os.path.join(dwdt,FileRx)]

			check = compareFilesTime(inputFiles,outputFiles)

			file_4_DEMULTIPLEXING = os.path.join(dir_4_DEMULTIPLEXING,FileRx+'''_Barcode_%''')

			mycmd='''fastq-multx -b ''' +os.path.join(utilsRef,FileRxBarcode)+ ''' ''' +os.path.join(dwdt,FileRx)+ ''' -o '''  + file_4_DEMULTIPLEXING

			if check:
				print mycmd
				subprocess.call(mycmd,shell=True)

		t1 = Thread(target=call_script, args=("R1_fastq_trim12nt_qcTrimmed_MatchBlastLtrLc",barcodeLTRseq))
		t2 = Thread(target=call_script, args=("R2_fastq_trim12nt_qcTrimmed_MatchBlastLtrLc",barcodeLCseq))

		t1.start()
		t2.start()
		t1.join()
		t2.join()
		print ("DEMULTIPLEXING\n")

def trimwithCutAdapt(dwdt,utilsRef,LTRFileFa,LCFileFa):

	# #### trim LTR and LC withcutadapt ####################
	def call_script(FileRx,virusFa):
		temp = os.path.basename(os.path.dirname(FileRx))
		dir_4_DEMULTIPLEXING = os.path.join(dwdt,temp+"Tofq")
		if not os.path.exists(dir_4_DEMULTIPLEXING):
			os.mkdir(dir_4_DEMULTIPLEXING)
		inputFile= FileRx
		outputFile = os.path.join(dir_4_DEMULTIPLEXING,os.path.basename(FileRx)+ '''.fq''')
		check = compareFileTime(inputFile,outputFile)
		if check:
			mycmd='''cp ''' +inputFile+''' '''+outputFile
			print mycmd
			subprocess.call(mycmd,shell=True)
			outputFile1 = outputFile + "_trimwithCutAdapt"
			mycmd='''cutadapt -b file:'''+os.path.join(utilsRef,virusFa)+ ''' ''' +outputFile+ ''' -o '''  +outputFile1
			print mycmd
			subprocess.call(mycmd,shell=True)

	dir_4_DEMULTIPLEXING = os.path.join(dwdt,"R1_fastq_trim12nt_qcTrimmed_MatchBlastLtrLcDEMULTIPLEXING")
	for filename in os.listdir(dir_4_DEMULTIPLEXING):
		if fnmatch.fnmatch(filename,'R1*_Barcode_*'):
			t1 = Thread(target=call_script, args=(os.path.join(dir_4_DEMULTIPLEXING,filename),LTRFileFa))
			t1.start()
			t1.join()

	dir_4_DEMULTIPLEXING = os.path.join(dwdt,"R2_fastq_trim12nt_qcTrimmed_MatchBlastLtrLcDEMULTIPLEXING")
	for filename in os.listdir(dir_4_DEMULTIPLEXING):
		if fnmatch.fnmatch(filename,'R2*_Barcode_*'):
			t1 = Thread(target=call_script, args=(os.path.join(dir_4_DEMULTIPLEXING,filename),LCFileFa))
			t1.start()
			t1.join()

def trimwithMCF(dwdt):

	##### trim LTR and LC withMCF ####################
	inputFiles = []
	for f in fnmatch.filter(os.listdir(dwdt),'*_MatchBlastLtrLc'):
		inputFiles.append(os.path.join(dwdt,f))

	outputFiles = []
	temp = dwdt + "R1_fastq_trim12nt_qcTrimmed_MatchBlastLtrLcDEMULTIPLEXINGtrimwithMCFtrimwithFlexbar"

	if not os.path.exists(temp):
		os.mkdir(temp)

	for f in fnmatch.filter(os.listdir(temp),'*'):
		outputFiles.append(os.path.join(temp,f))

	temp = dwdt + "R2_fastq_trim12nt_qcTrimmed_MatchBlastLtrLcDEMULTIPLEXINGtrimwithMCFtrimwithFlexbar"

	if not os.path.exists(temp):
		os.mkdir(temp)

	for f in fnmatch.filter(os.listdir(temp),'*'):
		outputFiles.append(os.path.join(temp,f))

	check = compareFilesTime(inputFiles,outputFiles)

	if check:
		def call_script(FileRx,virusFa):
			temp = os.path.basename(os.path.dirname(FileRx))
			dir_4_DEMULTIPLEXING = os.path.join(dwdt,temp+"trimwithMCF")

			if not os.path.exists(dir_4_DEMULTIPLEXING):
				os.mkdir(dir_4_DEMULTIPLEXING)


				inputFile= FileRx

				outputFile = os.path.join(dir_4_DEMULTIPLEXING,os.path.basename(FileRx)+ '''_trimwithMCF''')

				check = compareFileTime(inputFile,outputFile)

				if check:
					mycmd='''fastq-mcf -S ''' +os.path.join(utilsRef,virusFa)+ ''' ''' +inputFile+ ''' -o '''  +outputFile
					print mycmd
					subprocess.call(mycmd,shell=True)

		dir_4_DEMULTIPLEXING = os.path.join(dwdt,"R1_fastq_trim12nt_qcTrimmed_MatchBlastLtrLcDEMULTIPLEXING")
		for filename in os.listdir(dir_4_DEMULTIPLEXING):
			if fnmatch.fnmatch(filename,'R1*_Barcode_*'):
				t1 = Thread(target=call_script, args=(os.path.join(dir_4_DEMULTIPLEXING,filename),LTRFileFa))
				t1.start()
				t1.join()

		dir_4_DEMULTIPLEXING = os.path.join(dwdt,"R2_fastq_trim12nt_qcTrimmed_MatchBlastLtrLcDEMULTIPLEXING")
		for filename in os.listdir(dir_4_DEMULTIPLEXING):
			if fnmatch.fnmatch(filename,'R2*_Barcode_*'):
				t1 = Thread(target=call_script, args=(os.path.join(dir_4_DEMULTIPLEXING,filename),LCFileFa))
				t1.start()
				t1.join()


		#### trim skipped with FLEXBAR ####################
		def call_script(FileRx,virusFa):

			temp = os.path.basename(os.path.dirname(FileRx))
			dir_4_DEMULTIPLEXING = os.path.join(dwdt,temp+"skip2fq")

			if not os.path.exists(dir_4_DEMULTIPLEXING):
				os.mkdir(dir_4_DEMULTIPLEXING)


			inputFile= FileRx

			outputFile = os.path.join(dir_4_DEMULTIPLEXING,os.path.basename(FileRx)+ '''.fq''')

			check = compareFileTime(inputFile,outputFile)

			if check:
				mycmd='''cp ''' +inputFile+''' '''+outputFile
				print mycmd
				subprocess.call(mycmd,shell=True)

				outputFile1 = outputFile + "_trimwithFlexbar"
				mycmd='''flexbar -r ''' +outputFile+ ''' -as "''' +virusFa+'''" -ae 0.1 -at LEFT -t ''' +outputFile1
				print mycmd
				subprocess.call(mycmd,shell=True)

		dir_4_DEMULTIPLEXING = os.path.join(dwdt,"R1_fastq_trim12nt_qcTrimmed_MatchBlastLtrLcDEMULTIPLEXINGtrimwithMCF")
		for filename in os.listdir(dir_4_DEMULTIPLEXING):
			if fnmatch.fnmatch(filename,'R1*_Barcode_*skip'):
				t1 = Thread(target=call_script, args=(os.path.join(dir_4_DEMULTIPLEXING,filename),ltrFaSeq))
				t1.start()
				t1.join()

		dir_4_DEMULTIPLEXING = os.path.join(dwdt,"R2_fastq_trim12nt_qcTrimmed_MatchBlastLtrLcDEMULTIPLEXINGtrimwithMCF")
		for filename in os.listdir(dir_4_DEMULTIPLEXING):
			if fnmatch.fnmatch(filename,'R2*_Barcode_*skip'):
				t1 = Thread(target=call_script, args=(os.path.join(dir_4_DEMULTIPLEXING,filename),lcFaSeq))
				t1.start()
				t1.join()

		def mergeMCFandFLEXBAR(temp):
			for filename in os.listdir(temp):
				if fnmatch.fnmatch(filename,'*_fastq_trim12nt_qcTrimmed_MatchBlastLtrLc_Barcode_*_trimwithMCF'):
					infilename = os.path.join(temp,filename)
					flexName=os.path.join(temp+"skip2fq",filename+".skip.fq_trimwithFlexbar.fastq")

					temp1 = temp+"trimwithFlexbar"
					if not os.path.exists(temp1):
						os.mkdir(temp1)


					outName=os.path.join(temp1,filename+"_trimwithFlexbar")

					with open(outName, 'w') as outfile:
						with open(infilename) as infile:
							for line in infile:
								outfile.write(line)
						infile.close()
						if os.path.isfile(flexName):
							with open(flexName) as infile:
								for line in infile:
									outfile.write(line)
							infile.close()
					outfile.close()

		temp = os.path.join(dwdt,"R1_fastq_trim12nt_qcTrimmed_MatchBlastLtrLcDEMULTIPLEXINGtrimwithMCF")
		mergeMCFandFLEXBAR(temp)

		temp = os.path.join(dwdt,"R2_fastq_trim12nt_qcTrimmed_MatchBlastLtrLcDEMULTIPLEXINGtrimwithMCF")
		mergeMCFandFLEXBAR(temp)

def RandomBarcodRemoval(inputDir,inputPattern,outputDir,seqFile,seqFile1,RandomBarcodRemovalOut,NT):

	########Random Barcode removal start ! ##############
	def call_scriptFqFa(filename):
		inputFile = filename
		temp = os.path.basename(filename)

		outputFile = os.path.join(outputDir,temp+".fa")
		check = compareFileTime(inputFile,outputFile)

		if check:
			mycmd='''awk 'BEGIN{P=1}{if(P==1||P==2){gsub(/^[@]/,">");print}; if(P==4)P=0; P++}' ''' +inputFile+ ''' > '''  +outputFile
			print mycmd
			subprocess.call(mycmd,shell=True)

	threads = []
	for filename in os.listdir(inputDir):
		if fnmatch.fnmatch(filename,inputPattern):
			t1 = Thread(target=call_scriptFqFa, args=(os.path.join(inputDir,filename),))
			threads.append(t1)

	for i in xrange(0, len(threads), NT):
		for x in threads[i:i+NT]:
			x.start()
		for x in threads[i:i+NT]:
			x.join()

	## b
	threads = []
	def call_script_Blast(filename,seqFile):
		inputFile = filename
		outputFile = os.path.join(outputDir,filename[:-3]+ '''_blastAnchor.txt''')
		check = compareFileTime(inputFile,outputFile)
		if check:
			mycmd='''blastn -task blastn-short -word_size 4 -query ''' +inputFile+ ''' -subject '''+seqFile+''' -outfmt "6 std btop slen" -max_target_seqs 1 > '''  +outputFile
			print mycmd
			subprocess.call(mycmd,shell=True)

	for filename in os.listdir(outputDir):
		if fnmatch.fnmatch(filename,'*.fq_trimwithCutAdapt.fa'):
			t1 = Thread(target=call_script_Blast, args=(os.path.join(outputDir,filename),seqFile))
			threads.append(t1)

	for i in xrange(0, len(threads), NT):
		for x in threads[i:i+NT]:
			x.start()
		for x in threads[i:i+NT]:
			x.join()

	for filename in os.listdir(outputDir):
		if fnmatch.fnmatch(filename,'*_blastAnchor.txt'):

			inputFile = os.path.join(outputDir,filename)
			outputFile = os.path.join(outputDir,filename[:-4]+ '''_7.txt''')

			check = compareFileTime(inputFile,outputFile)

			if check:
				goodID=[]
				myblast=open(inputFile)
				myout_7=open(outputFile,'w')

				for l in myblast.readlines():
					l= l.rstrip()
					if l[0]=='#':
						continue
					l= l.split('\t')
					if int(l[4])==0 and int(l[5])==0:
						if int(l[6])==7 and int(l[7])==12:
							myout_7.write('\t'.join(l)+'\n')
							goodID+=[l[0]]
					#if int(l[6])==6 and int(l[7])==11:
					#    myout_6.write('\t'.join(l)+'\n')
				goodIDset = set(goodID)
				with open(os.path.join(outputDir,filename[:-4]+'_ID_blast_7.lst'), "w") as myblastOutID:
					myblastOutID.write('\n'.join(goodIDset))
				myblast.close()
				myblastOutID.close()


	def call_script(filename):

		inputFile = filename
		temp = os.path.basename(filename)

		inputFile1 = os.path.join(outputDir,temp+"_blastAnchor_ID_blast_7.lst")

		outputFile = os.path.join(outputDir,temp+"_reads_withAncora")

		check = compareFileTime(inputFile,outputFile)

		if check:
			mycmd='''seqtk subseq '''+inputFile+ ''' ''' +inputFile1 + ''' > '''+outputFile
			print mycmd
			subprocess.call(mycmd,shell=True)


	threads = []
	for filename in os.listdir(inputDir,):
		if fnmatch.fnmatch(filename,inputPattern):
			t1 = Thread(target=call_script, args=(os.path.join(inputDir,filename),))
			threads.append(t1)

	for i in xrange(0, len(threads), NT):
		# Start all threads
		for x in threads[i:i+NT]:
			x.start()
		# Wait for all of them to finish
		for x in threads[i:i+NT]:
			x.join()

	def call_script(filename):
		inputFile = filename
		temp = os.path.basename(filename)

		outputFile = os.path.join(outputDir,temp+"_out_trim18")
		check = compareFileTime(inputFile,outputFile)
		if check:
			mycmd='''fastx_trimmer -f 19 -Q 33 -i '''+inputFile+''' -o '''+outputFile
			print mycmd
			subprocess.call(mycmd,shell=True)

	threads = []
	for filename in os.listdir(outputDir):
		if fnmatch.fnmatch(filename,'*_reads_withAncora'):
			t1 = Thread(target=call_script, args=(os.path.join(outputDir,filename),))
			threads.append(t1)

	for i in xrange(0, len(threads), NT):
		# Start all threads
		for x in threads[i:i+NT]:
			x.start()
			# Wait for all of them to finish
		for x in threads[i:i+NT]:
			x.join()

	threads = []
	def call_scriptFqFa(filename):
		inputFile = filename
		print inputFile
		temp = os.path.basename(filename)
		outputFile = os.path.join(outputDir,temp+".fa")

		check = compareFileTime(inputFile,outputFile)

		if check:
			mycmd='''awk 'BEGIN{P=1}{if(P==1||P==2){gsub(/^[@]/,">");print}; if(P==4)P=0; P++}' ''' +inputFile+ ''' > '''  +outputFile
			print mycmd
			subprocess.call(mycmd,shell=True)

	for filename in os.listdir(outputDir):
		if fnmatch.fnmatch(filename,'*withAncora_out_trim18'):
			t1 = Thread(target=call_scriptFqFa, args=(os.path.join(outputDir,filename),))
			threads.append(t1)

	for i in xrange(0, len(threads), NT):
		for x in threads[i:i+NT]:
			x.start()
		for x in threads[i:i+NT]:
			x.join()

	threads = []
	def call_script_Blast(filename,seqFile1):
		inputFile = filename
		temp = os.path.basename(filename)
		temp = temp[:-3]+ '''_blastLast.txt'''

		outputFile = os.path.join(outputDir,temp)
		check = compareFileTime(inputFile,outputFile)

		if check:
			mycmd='''blastn -task blastn -word_size 7 -query ''' +inputFile+ ''' -subject '''+seqFile1+''' -outfmt "6 std btop slen" -max_target_seqs 1 > '''  +outputFile
			print mycmd
			subprocess.call(mycmd,shell=True)

	for filename in os.listdir(outputDir):
		if fnmatch.fnmatch(filename,'*withAncora_out_trim18.fa'):
			t1 = Thread(target=call_script_Blast, args=(os.path.join(outputDir,filename),seqFile1))
			threads.append(t1)

	for i in xrange(0, len(threads), NT):
		for x in threads[i:i+NT]:
			x.start()
		for x in threads[i:i+NT]:
			x.join()

	for filename in os.listdir(outputDir):
		if fnmatch.fnmatch(filename,'*_blastLast.txt'):
			inputFile = os.path.join(outputDir,filename)
			outputFile = os.path.join(outputDir,filename[:-14]+'_ID_blast_last.lst')
			check = compareFileTime(inputFile,outputFile)

			if check:
				goodID=[]
				myblast=open(inputFile)
				for l in myblast.readlines():
					l= l.rstrip()
					if l[0]=='#':
						continue
					l= l.split('\t')
					if int(l[4])<=3:
						if int(l[6])<=2 and int(l[7])>=19:
							goodID+=[l[0]]
				myblast.close()
				goodIDset = set(goodID)
				with open(outputFile, "w") as myblastOutID:
					myblastOutID.write('\n'.join(goodIDset))


	def call_script(filename):
		inputFile = filename
		temp = os.path.basename(filename)

		inputFile1 = os.path.join(outputDir,temp+"_ID_blast_last.lst")

		outputFile = os.path.join(outputDir,temp+ "_reads_withLast")

		check = compareFileTime(inputFile,outputFile)

		if check:
			mycmd='''seqtk subseq '''  +inputFile+ ''' ''' +inputFile1+''' > '''  +outputFile
			print mycmd
			subprocess.call(mycmd,shell=True)

	threads = []
	for filename in os.listdir(outputDir):
		if fnmatch.fnmatch(filename,'R2_fastq*_out_trim18'):
			t1 = Thread(target=call_script, args=(os.path.join(outputDir,filename),))
			threads.append(t1)

	for i in xrange(0, len(threads), NT):
		# Start all threads
		for x in threads[i:i+NT]:
			x.start()
		# Wait for all of them to finish
		for x in threads[i:i+NT]:
			x.join()

	def call_script(filename):
		inputFile = filename

		temp1 = os.path.basename(filename)

		outputFile = os.path.join(RandomBarcodRemovalOut,temp1[:-43])

		check = compareFileTime(inputFile,outputFile)

		if check:
			mycmd='''fastx_trimmer -f 21 -Q 33 -i '''+inputFile+''' -o '''+outputFile
			print mycmd
			subprocess.call(mycmd,shell=True)

	threads = []
	for filename in os.listdir(outputDir):
		if fnmatch.fnmatch(filename,'*_reads_withLast'):
			t1 = Thread(target=call_script, args=(os.path.join(outputDir,filename),))
			threads.append(t1)

	for i in xrange(0, len(threads), NT):
		# Start all threads
		for x in threads[i:i+NT]:
			x.start()
		 # Wait for all of them to finish
		for x in threads[i:i+NT]:
			x.join()

def clusterUmi(filename):

	inputFile = filename

	temp1 = os.path.basename(filename)
	outDir = os.path.dirname(filename)

	outputFile = os.path.join(outDir,temp1+"_cluster")

	check = compareFileTime(inputFile,outputFile)

	if check:
		mycmd='''starcode --seq-id -d0  -r3 -t1 '''+inputFile+''' -o '''+outputFile
		print mycmd
		subprocess.call(mycmd,shell=True)

def extractUmiOnly(filename,outDir):

	########extract UMI from reads ! ####

	temp = os.path.basename(filename)

	outDir1=outDir

	outputFile = os.path.join(outDir1,temp+ "_umi_extract")

	check = compareFileTime(filename,outputFile)

	if check:

		mycmd='''eval "$(conda shell.bash hook)" && conda activate py3.7 && umi_tools extract --extract-method=string -p NNNNNNNNNNNNNNNNNN -I '''+filename+''' -S '''+outputFile
		print mycmd
		subprocess.call(mycmd,shell=True)

def getCollisionTable(seqPlat,R1Out,R2Out,outputDir,sampleResearch,utilsDir,utilsRef,chrList,dirGenome,sampleName,sortedKnownGene,genomeSorted,vectorBed,suffix,NT,maskFile,VectorMask,repRegQual,PreviousCollision):

	### getCollisionTable function #########
	alignOut = os.path.join(outputDir,"align")
	if not os.path.exists(alignOut):
		os.makedirs(alignOut)

	inputFiles = []
	for f in fnmatch.filter(os.listdir(R1Out),'*.fq_trimwithCutAdapt'):
		inputFiles.append(os.path.join(R1Out,f))

	for f in fnmatch.filter(os.listdir(R2Out),'*.fq_trimwithCutAdapt'):
		inputFiles.append(os.path.join(R2Out,f))

	outputFiles = []

	for f in fnmatch.filter(os.listdir(alignOut),'*_aligned_mem.sam'):
		outputFiles.append(os.path.join(alignOut,f))

	print inputFiles

	print outputFiles

	check = compareFilesTime(inputFiles,outputFiles)

	print check

	if check:
		if (seqPlat=="MiSeq"):
			pattern = re.compile('^@M')
		if (seqPlat=="HiSeq"):
			pattern = re.compile('^@HISEQ')

		with open(sampleResearch) as barCodeFile:
			for barCodeline in barCodeFile:
				print barCodeline
				barCodeline=barCodeline.rstrip()
				lineSplitBrc=barCodeline.split(',')

				print lineSplitBrc[9]
				print "OK"
				print sampleName

				if (lineSplitBrc[9]==sampleName):

					print "matched"

					print lineSplitBrc[7]

					print lineSplitBrc[8]

					r1filename=os.path.join(R1Out,"R1_fastq_trim12nt_qcTrimmed_MatchBlastLtrLc_Barcode_"+lineSplitBrc[7]+".fq_trimwithCutAdapt")
					r2filename=os.path.join(R2Out,"R2_fastq_trim12nt_qcTrimmed_MatchBlastLtrLc_Barcode_"+lineSplitBrc[8]+".fq_trimwithCutAdapt")
					outName=os.path.join(outputDir,"R1_R2_Barcode_"+lineSplitBrc[7]+"_"+lineSplitBrc[8]+"_trimmedID")

					print r1filename
					print r2filename
					print outName

					r1ID=[]
					r2ID=[]
					r12Int=[]
					if (os.path.isfile(r1filename) and os.path.isfile(r2filename)):
						print r1filename
						print r2filename
						with open(r1filename) as r1infile:
							print r1infile
							for r1line in r1infile:
								if pattern.search(r1line):
									print r1line
									lineSplit=r1line.split(' ')
									r1ID.append(lineSplit[0][1:])
									print lineSplit[0][1:]
						with open(r2filename) as r2infile:
							for r2line in r2infile:
								if pattern.search(r2line):
									lineSplit=r2line.split(' ')
									r2ID.append(lineSplit[0][1:])
									print lineSplit[0][1:]

						r12Int=set(r1ID).intersection(r2ID)

						print r12Int
						with open(outName, 'w') as outfile:
							outfile.write("\n".join(r12Int))

						mycmd='''seqtk subseq '''  +r1filename+ ''' ''' +outName+ ''' > '''  +r1filename+ '''_ReadyToAlign'''
						print mycmd
						process = subprocess.Popen(mycmd,shell=True)
						process.wait()

						mycmd='''fastqutils sort ''' +r1filename+ '''_ReadyToAlign > ''' +r1filename+ '''_ReadyToAlignSort'''
						print mycmd
						process = subprocess.Popen(mycmd,shell=True)
						process.wait()

						mycmd='''seqtk subseq '''  +r2filename+ ''' ''' +outName+ ''' > '''  +r2filename+ '''_ReadyToAlign'''
						print mycmd
						process = subprocess.Popen(mycmd,shell=True)
						process.wait()

						mycmd='''fastqutils sort ''' +r2filename+ '''_ReadyToAlign > ''' +r2filename+ '''_ReadyToAlignSort'''
						print mycmd
						process = subprocess.Popen(mycmd,shell=True)
						process.wait()

						mycmd='''bwa-mem mem -t 8 '''+dirGenome+''' '''+r1filename+ '''_ReadyToAlignSort ''' +r2filename+ '''_ReadyToAlignSort > ''' +os.path.join(alignOut,"R1_R2_Barcode_"+lineSplitBrc[7]+"_"+lineSplitBrc[8]+"_aligned_mem.sam")
						print mycmd
						process = subprocess.Popen(mycmd,shell=True)
						process.wait()


	rScript='getNumOfRead.R'
	output = os.path.join(os.path.dirname(outputDir),'TotalReads','totalReads.rds')
	callR(utilsDir,rScript,alignOut,output)

	### BAM -> SAM #########
	bamOut = os.path.join(outputDir,"BAM")
	if not os.path.exists(bamOut):
		os.makedirs(bamOut)

	inputFiles = []
	for f in fnmatch.filter(os.listdir(alignOut),'*_aligned_mem.sam'):
		inputFiles.append(os.path.join(alignOut,f))

	outputFiles = []
	for f in fnmatch.filter(os.listdir(bamOut),'*'):
		outputFiles.append(os.path.join(bamOut,f))

	print inputFiles

	print outputFiles

	check = compareFilesTime(inputFiles,outputFiles)

	print check

	if check:
		threads = []

		def call_script_bam(filename):
			temp1 = os.path.basename(filename)
			outputFile = os.path.join(bamOut,temp1[:-4]+'''.bam ''')
			mycmd='''samtools view -bS ''' +filename+ ''' > '''  +outputFile
			print mycmd
			subprocess.call(mycmd,shell=True)

		for filename in os.listdir(alignOut):
			if fnmatch.fnmatch(filename,'*_aligned_mem.sam'):
				t1 = Thread(target=call_script_bam, args=(os.path.join(alignOut,filename),))
				threads.append(t1)


		for i in xrange(0, len(threads), NT):
			# Start all threads
			for x in threads[i:i+NT]:
				x.start()
	 		# Wait for all of them to finish
			for x in threads[i:i+NT]:
				x.join()

		###### filter out unmapped ###########

		threads = []
		def call_script_filterMap(filename):

			temp1 = os.path.basename(filename)
			outputFile = os.path.join(bamOut,temp1[:-4]+'''_mapped.bam ''')

			mycmd='''samtools view -b -F 4 '''+filename+''' > '''+outputFile
			print mycmd
			subprocess.call(mycmd,shell=True)


		for filename in os.listdir(bamOut):
			if fnmatch.fnmatch(filename,'*.bam'):
				t1 = Thread(target=call_script_filterMap, args=(os.path.join(bamOut,filename),))
				threads.append(t1)


		for i in xrange(0, len(threads), NT):
			# Start all threads
			for x in threads[i:i+NT]:
				x.start()
	 		# Wait for all of them to finish
			for x in threads[i:i+NT]:
				x.join()

		##### filter out non primary mapping ############

		threads = []
		def call_script_filterMap(filename):

			temp1 = os.path.basename(filename)
			outputFile = os.path.join(bamOut,temp1[:-11]+'''_mapped_primary.bam ''')

			mycmd='''samtools view -b -F 256 '''+filename+''' > '''+outputFile
			print mycmd
			subprocess.call(mycmd,shell=True)


		for filename in os.listdir(bamOut):
			if fnmatch.fnmatch(filename,'*_mapped.bam'):
				t1 = Thread(target=call_script_filterMap, args=(os.path.join(bamOut,filename),))
				threads.append(t1)


		for i in xrange(0, len(threads), NT):
			# Start all threads
			for x in threads[i:i+NT]:
				x.start()
			# Wait for all of them to finish
			for x in threads[i:i+NT]:
				x.join()


	####### sort bam file ##########

	bamSortedOut = os.path.join(outputDir,"BAMSorted")
	if not os.path.exists(bamSortedOut):
		os.makedirs(bamSortedOut)

	inputFiles = []
	for f in fnmatch.filter(os.listdir(bamOut),'*'):
		inputFiles.append(os.path.join(bamOut,f))

	outputFiles = []
	for f in fnmatch.filter(os.listdir(bamSortedOut),'*'):
		outputFiles.append(os.path.join(bamSortedOut,f))

	print inputFiles

	print outputFiles

	check = compareFilesTime(inputFiles,outputFiles)

	print check

	if check:
		threads = []

		def call_script_filterMap(filename):
			temp1 = os.path.basename(filename)
			outputFile = os.path.join(bamSortedOut,temp1[:-19]+"_mapped_primary_sort.bam")
			mycmd='''samtools sort '''+filename+''' -o '''+outputFile
			print(mycmd)
			subprocess.call(mycmd,shell=True)

		for filename in os.listdir(bamOut):
			if fnmatch.fnmatch(filename,'*_mapped_primary.bam'):
				t1 = Thread(target=call_script_filterMap, args=(os.path.join(bamOut,filename),))
				threads.append(t1)

		for i in xrange(0, len(threads), NT):
			# Start all threads
			for x in threads[i:i+NT]:
				x.start()
	 		# Wait for all of them to finish
			for x in threads[i:i+NT]:
				x.join()

		###### indexing bam file ###########

		threads = []
		def call_script_filterMap(filename):
			mycmd='''samtools index '''+filename
			print mycmd
			subprocess.call(mycmd,shell=True)


		for filename in os.listdir(bamSortedOut):
			if fnmatch.fnmatch(filename,'*_mapped_primary_sort.bam'):
				t1 = Thread(target=call_script_filterMap, args=(os.path.join(bamSortedOut,filename),))
				threads.append(t1)

		for i in xrange(0, len(threads), NT):
			# Start all threads
			for x in threads[i:i+NT]:
				x.start()
	 		# Wait for all of them to finish
			for x in threads[i:i+NT]:
				x.join()

		#######select reads in repeated ##############

		threads = []
		def call_script_filterMap(filename):
			temp1 = os.path.basename(filename)
			outputFile = os.path.join(bamSortedOut,temp1[:-24]+"_sort_inMask.bam")
			#mycmd='''samtools view -b -L '''+maskFile+''' '''+filename+''' > '''+outputFile
			mycmd='''bedtools intersect -abam '''+filename+''' -b '''+maskFile+ ''' > '''  +outputFile
			print mycmd
			subprocess.call(mycmd,shell=True)


		for filename in os.listdir(bamSortedOut):
			if fnmatch.fnmatch(filename,'*_mapped_primary_sort.bam'):
				t1 = Thread(target=call_script_filterMap, args=(os.path.join(bamSortedOut,filename),))
				threads.append(t1)


		for i in xrange(0, len(threads), NT):
			# Start all threads
			for x in threads[i:i+NT]:
				x.start()
	 		# Wait for all of them to finish
			for x in threads[i:i+NT]:
				x.join()


		#######quality filter on reads falling in repeated param: _repRegQual_ ##############

		threads = []
		def call_script_filterMap(filename):
			temp1 = os.path.basename(filename)
			outputFile = os.path.join(bamSortedOut,temp1[:-16]+"_sort_inMask_qual.bam")
			mycmd='''samtools view -bq '''+ str(repRegQual)+''' '''+filename+''' > '''  +outputFile
			print mycmd
			subprocess.call(mycmd,shell=True)


		for filename in os.listdir(bamSortedOut):
			if fnmatch.fnmatch(filename,'*_sort_inMask.bam'):
				t1 = Thread(target=call_script_filterMap, args=(os.path.join(bamSortedOut,filename),))
				threads.append(t1)


		for i in xrange(0, len(threads), NT):
			# Start all threads
			for x in threads[i:i+NT]:
				x.start()
	 		# Wait for all of them to finish
			for x in threads[i:i+NT]:
				x.join()

		#######extract reads in non-repeated ##############

		threads = []
		def call_script_filterMap(filename):
			temp1 = os.path.basename(filename)
			outputFile = os.path.join(bamSortedOut,temp1[:-24]+'''_sort_nonMask.bam''')
			mycmd='''bedtools intersect -v -abam '''+filename+''' -b '''+maskFile+ ''' > '''  +outputFile
			print mycmd
			subprocess.call(mycmd,shell=True)

		for filename in os.listdir(bamSortedOut):
			if fnmatch.fnmatch(filename,'*_mapped_primary_sort.bam'):
				t1 = Thread(target=call_script_filterMap, args=(os.path.join(bamSortedOut,filename),))
				threads.append(t1)

		for i in xrange(0, len(threads), NT):
			# Start all threads
			for x in threads[i:i+NT]:
				x.start()
	 		# Wait for all of them to finish
			for x in threads[i:i+NT]:
				x.join()

		#######merge non-repeated and high quality repeated ##############

		threads = []
		def call_script_filterMap(filename):
			temp1 = os.path.basename(filename)
			outputFile = os.path.join(bamSortedOut,temp1[:-17]+'''_allFilter.bam''')
			mycmd='''samtools merge '''+outputFile+''' '''+filename[:-17]+'''_sort_nonMask.bam '''+filename[:-17]+'''_sort_inMask_qual.bam'''
			print mycmd
			subprocess.call(mycmd,shell=True)


		for filename in os.listdir(bamSortedOut):
			if fnmatch.fnmatch(filename,'*_sort_nonMask.bam'):
				t1 = Thread(target=call_script_filterMap, args=(os.path.join(bamSortedOut,filename),))
				threads.append(t1)

		for i in xrange(0, len(threads), NT):
			# Start all threads
			for x in threads[i:i+NT]:
				x.start()
	 		# Wait for all of them to finish
			for x in threads[i:i+NT]:
				x.join()

	####### rehead  allFilter.bam ##############
	inputFiles = []
	for f in fnmatch.filter(os.listdir(bamSortedOut),'*_allFilter.bam'):
		inputFiles.append(os.path.join(bamSortedOut,f))

	outputFiles = []
	for f in fnmatch.filter(os.listdir(bamSortedOut),'*_allFilter_rehead.bam'):
		outputFiles.append(os.path.join(bamSortedOut,f))

	print inputFiles

	print outputFiles

	check = compareFilesTime(inputFiles,outputFiles)

	if check:
		threads = []
		def call_script_filterMap(filename):
			temp1 = os.path.basename(filename)
			outputFile = os.path.join(bamSortedOut,temp1[:-14]+'''_allFilter_rehead.bam''')
			mycmd='''PicardCommandLine AddOrReplaceReadGroups I='''+filename+''' O='''+outputFile+''' SORT_ORDER=coordinate RGID=foo RGLB=bar RGPL=illumina RGPU=shoot RGSM=DePristo'''
			print mycmd
			subprocess.call(mycmd,shell=True)


		for filename in os.listdir(bamSortedOut):
			if fnmatch.fnmatch(filename,'*_allFilter.bam'):
				t1 = Thread(target=call_script_filterMap, args=(os.path.join(bamSortedOut,filename),))
				threads.append(t1)


		for i in xrange(0, len(threads), NT):
			# Start all threads
			for x in threads[i:i+NT]:
				x.start()
	 		# Wait for all of them to finish
			for x in threads[i:i+NT]:
				x.join()

		###### indexing bam file ###########

		threads = []
		def call_script_filterMap(filename):
			mycmd='''samtools index '''+filename
			print mycmd
			subprocess.call(mycmd,shell=True)


		for filename in os.listdir(bamSortedOut):
			if fnmatch.fnmatch(filename,'*_allFilter_rehead.bam'):
				t1 = Thread(target=call_script_filterMap, args=(os.path.join(bamSortedOut,filename),))
				threads.append(t1)


		for i in xrange(0, len(threads), NT):
			# Start all threads
			for x in threads[i:i+NT]:
				x.start()
	 		# Wait for all of them to finish
			for x in threads[i:i+NT]:
				x.join()

	###### check R1 first 3 nt ###########
	inputFiles = []
	for f in fnmatch.filter(os.listdir(bamSortedOut),'*_allFilter_rehead.bam'):
		inputFiles.append(os.path.join(bamSortedOut,f))

	outputFiles = []
	for f in fnmatch.filter(os.listdir(bamSortedOut),'*_allFilter_rehead_exact3nt.bam'):
		outputFiles.append(os.path.join(bamSortedOut,f))

	print inputFiles

	print outputFiles

	check = compareFilesTime(inputFiles,outputFiles)

	print check

	if check:
		threads = []
		def call_script_filterMap(filename):
			mycmd='''python2 '''+os.path.join(utilsDir,'''pysam_parse.py ''') +filename
			print mycmd
			subprocess.call(mycmd,shell=True)

		for filename in os.listdir(bamSortedOut):
			if fnmatch.fnmatch(filename,'*_allFilter_rehead.bam'):
				t1 = Thread(target=call_script_filterMap, args=(os.path.join(bamSortedOut,filename),))
				threads.append(t1)

		for i in xrange(0, len(threads), NT):
			# Start all threads
			for x in threads[i:i+NT]:
				x.start()
	 		# Wait for all of them to finish
			for x in threads[i:i+NT]:
				x.join()

	###### count unmapped ###########
	inputFiles = []
	for f in fnmatch.filter(os.listdir(bamSortedOut),'*_allFilter_rehead_exact3nt.bam'):
		inputFiles.append(os.path.join(bamSortedOut,f))

	outputFiles = []
	for f in fnmatch.filter(os.listdir(bamSortedOut),'*_allFilter_rehead_exact3nt_nonSupplementary.bam'):
		outputFiles.append(os.path.join(bamSortedOut,f))

	print inputFiles

	print outputFiles

	check = compareFilesTime(inputFiles,outputFiles)

	print check

	if check:
		threads = []
		def call_script_filterMap(filename):
			temp1 = os.path.basename(filename)
			outputFile = os.path.join(bamSortedOut,temp1[:-30]+'''_allFilter_rehead_exact3nt_unmapped.bam''')
			mycmd='''samtools view -b -f 4 '''+filename+''' > ''' + outputFile
			print mycmd
			subprocess.call(mycmd,shell=True)


		for filename in os.listdir(bamSortedOut):
			if fnmatch.fnmatch(filename,'*_allFilter_rehead_exact3nt.bam'):
				t1 = Thread(target=call_script_filterMap, args=(os.path.join(bamSortedOut,filename),))
				threads.append(t1)

		for i in xrange(0, len(threads), NT):
			# Start all threads
			for x in threads[i:i+NT]:
				x.start()
	 		# Wait for all of them to finish
			for x in threads[i:i+NT]:
				x.join()


		###### count supplementary ###########

		threads = []
		def call_script_filterMap(filename):
			temp1 = os.path.basename(filename)
			bamSortedOut= os.path.dirname(filename)
			outputFile = os.path.join(bamSortedOut,temp1[:-30]+'''_allFilter_rehead_exact3nt_supplementary.bam''')
			mycmd='''samtools view -b -f 2048 '''+filename+''' > '''  +outputFile
			print mycmd
			subprocess.call(mycmd,shell=True)


		for filename in os.listdir(bamSortedOut):
			if fnmatch.fnmatch(filename,'*_allFilter_rehead_exact3nt.bam'):
				t1 = Thread(target=call_script_filterMap, args=(os.path.join(bamSortedOut,filename),))
				threads.append(t1)

		for i in xrange(0, len(threads), NT):
			# Start all threads
			for x in threads[i:i+NT]:
				x.start()
	 		# Wait for all of them to finish
			for x in threads[i:i+NT]:
				x.join()

		###### select non-supplementary ###########

		threads = []
		def call_script_filterMap(filename):
			temp1 = os.path.basename(filename)
			bamSortedOut = os.path.dirname(filename)
			outputFile = os.path.join(bamSortedOut,temp1[:-30]+'''_allFilter_rehead_exact3nt_nonSupplementary.bam''')
			mycmd='''samtools view -b -F 2048 '''+filename+''' > '''  +outputFile
			print mycmd
			subprocess.call(mycmd,shell=True)

		for filename in os.listdir(bamSortedOut):
			if fnmatch.fnmatch(filename,'*_allFilter_rehead_exact3nt.bam'):
				t1 = Thread(target=call_script_filterMap, args=(os.path.join(bamSortedOut,filename),))
				threads.append(t1)

		for i in xrange(0, len(threads), NT):
			# Start all threads
			for x in threads[i:i+NT]:
				x.start()
	 		# Wait for all of them to finish
			for x in threads[i:i+NT]:
				x.join()


		###### indexing non-supplementary bam file ###########

		threads = []
		def call_script_filterMap(filename):
			mycmd='''samtools index '''+filename
			print mycmd
			subprocess.call(mycmd,shell=True)

		for filename in os.listdir(bamSortedOut):
			if fnmatch.fnmatch(filename,'*_allFilter_rehead_exact3nt_nonSupplementary.bam'):
				t1 = Thread(target=call_script_filterMap, args=(os.path.join(bamSortedOut,filename),))
				threads.append(t1)

		for i in xrange(0, len(threads), NT):
			# Start all threads
			for x in threads[i:i+NT]:
				x.start()
	 		# Wait for all of them to finish
			for x in threads[i:i+NT]:
				x.join()


	###### export infos ###########
	inputFiles = []
	for f in fnmatch.filter(os.listdir(bamSortedOut),'*_allFilter_rehead_exact3nt_nonSupplementary.bam'):
		inputFiles.append(os.path.join(bamSortedOut,f))

	outputFiles = []
	outIS = os.path.join(outputDir,"ISout")
	if not os.path.exists(outIS):
		os.makedirs(outIS)

	for f in fnmatch.filter(os.listdir(outIS),'*.txt'):
		outputFiles.append(os.path.join(outIS,f))

	check = compareFilesTime(inputFiles,outputFiles)

	print check

	if check:
		threads = []
		def call_script_filterMap(filename):
			mycmd='''python '''+os.path.join(utilsDir,'''try_pysam.py ''')+filename+ ''' '''+ sampleName
			print mycmd
			subprocess.call(mycmd,shell=True)

		for filename in os.listdir(bamSortedOut):
			if fnmatch.fnmatch(filename,'*_allFilter_rehead_exact3nt_nonSupplementary.bam'):
				print filename
				t1 = Thread(target=call_script_filterMap, args=(os.path.join(bamSortedOut,filename),))
				threads.append(t1)


		for i in xrange(0, len(threads), NT):
			# Start all threads
			for x in threads[i:i+NT]:
				x.start()
	 		# Wait for all of them to finish
			for x in threads[i:i+NT]:
				x.join()

		###### filter alignment 60 bp ###########

		for filename in os.listdir(outputDir):
			if fnmatch.fnmatch(filename,'*_final_parse_filterNo.txt'):
				myblast=open(os.path.join(outputDir,filename))
				myblastout=os.path.join(outputDir,filename[:-13]+"_filter60.txt")
				with open(myblastout, 'w') as outfile:
					for l in myblast.readlines():
						l= l.rstrip()
						if l[0]=='#':
							continue
						l= l.split('\t')
						#n = int(l[3]) if l[3].is_integer() else int(float(l[3]))
						print l
						if int(l[3])>=60:
							outfile.write('\t'.join(l))
							outfile.write('\n')

		###### IS ###########

		threads = []
		def call_script_filterMap(filename):
			mycmd='''python2 '''+os.path.join(utilsDir,'''try.py ''')+filename
			print mycmd
			subprocess.call(mycmd,shell=True)

		for filename in os.listdir(outputDir):
			if fnmatch.fnmatch(filename,'*_final_parse_*.txt'):
				t1 = Thread(target=call_script_filterMap, args=(os.path.join(outputDir,filename),))
				threads.append(t1)


		for i in xrange(0, len(threads), NT):
			# Start all threads
			for x in threads[i:i+NT]:
				x.start()
	 		# Wait for all of them to finish
			for x in threads[i:i+NT]:
				x.join()


		###### Calculate diversity ###########
	inputFiles = []
	for f in fnmatch.filter(os.listdir(outputDir),'*_final_parse_*.txt'):
		inputFiles.append(os.path.join(outputDir,f))

	outputFiles = []

	DiversityOut = os.path.join(outputDir,"DiversityOut")

	if not os.path.exists(DiversityOut):
		os.makedirs(DiversityOut)

	for f in fnmatch.filter(os.listdir(DiversityOut),'*.txt'):
		outputFiles.append(os.path.join(DiversityOut,f))

	check = compareFilesTime(inputFiles,outputFiles)

	print check

	if check:
		threads = []
		def call_script_filterMap(filename):
			mycmd='''python '''+os.path.join(utilsDir,'''calculate_diversity.py ''') +filename+ ''' '''+chrList+''' '''+outIS
			print mycmd
			subprocess.call(mycmd,shell=True)

		for filename in os.listdir(outputDir):
			if fnmatch.fnmatch(filename,'*_final_parse_*.txt'):
				t1 = Thread(target=call_script_filterMap, args=(os.path.join(outputDir,filename),))
				threads.append(t1)

		for i in xrange(0, len(threads), NT):
			# Start all threads
			for x in threads[i:i+NT]:
				x.start()
	 		# Wait for all of them to finish
			for x in threads[i:i+NT]:
				x.join()


	#inputFiles = []
	#for f in fnmatch.filter(os.listdir(DiversityOut),'*_final_parse_*.txt'):
	#	inputFiles.append(os.path.join(DiversityOut,f))

	#inputPattern="filter60"
	#dwdcFilter=os.path.join(outputDir,inputPattern+"/db/")
	#print "collision Folder "+inputPattern+':' +dwdcFilter60

	#if not os.path.exists(dwdcFilter):
	#	os.makedirs(dwdcFilter)

	#outputCollFiles = []
	#for f in fnmatch.filter(os.listdir(dwdcFilter),'*'):
	#	outputCollFiles.append(os.path.join(dwdcFilter,f))

	#check = compareFilesTime(inputFiles,outputCollFiles)
	#check = True

	#print check


	inputPattern="filter60"
	reNameFile(DiversityOut,outputDir,inputPattern,sampleResearch)
	dwdcFilter=os.path.join(outputDir,inputPattern+"/db")

	if PreviousCollision == "nothing":
		PreviousGroupedISfolder="nothing"
	else:
		PreviousGroupedISfolder=os.path.join(PreviousCollision,"CutAdapt",inputPattern+"/db")

	sumToTable(utilsDir,dwdcFilter,suffix,PreviousGroupedISfolder)
	dwdcFilterSuffix=os.path.join(dwdcFilter,suffix)
	annotateISite(dwdcFilterSuffix,sortedKnownGene,genomeSorted,NT,utilsRef,utilsDir,vectorBed,VectorMask,suffix)

	inputPattern="filterNo"
	reNameFile(DiversityOut,outputDir,inputPattern,sampleResearch)
	dwdcFilter=os.path.join(outputDir,inputPattern+"/db")

	if PreviousCollision == "nothing":
		PreviousGroupedISfolder="nothing"
	else:
		PreviousGroupedISfolder=os.path.join(PreviousCollision,"CutAdapt",inputPattern+"/db")

	sumToTable(utilsDir,dwdcFilter,suffix,PreviousGroupedISfolder)
	dwdcFilterSuffix=os.path.join(dwdcFilter,suffix)
	annotateISite(dwdcFilterSuffix,sortedKnownGene,genomeSorted,NT,utilsRef,utilsDir,vectorBed,VectorMask,suffix)

def getMissIsTable(outputDir,sampleResearch,utilsDir,utilsRef,chrList,dirGenome,sampleName,sortedKnownGene,genomeSorted,vectorBed,suffix,NT,maskFile,VectorMask):

	####### Rename files ###################
	dwdcFilterNo=os.path.join(outputDir,"filterNo","db")
	print "collision Folder filterNo:" +dwdcFilterNo

	if not os.path.exists(dwdcFilterNo):
		os.makedirs(dwdcFilterNo)

	inputFiles = []
	for f in fnmatch.filter(os.listdir(outputDir),'*filterNo_grouped_IS.txt'):
		inputFiles.append(os.path.join(outputDir,f))

	outputFiles = []
	for f in fnmatch.filter(os.listdir(dwdcFilterNo),'*_grouped_IS'):
		outputFiles.append(os.path.join(dwdcFilterNo,f))

	print inputFiles
	print outputFiles

	check = compareFilesTime(inputFiles,outputFiles)

	print check

	if check:
		my_file=[]
		for filename in os.listdir(outputDir):
			if fnmatch.fnmatch(filename,'*filterNo_grouped_IS.txt'):
				my_file.append(os.path.join(outputDir,filename))

			#Open a legend were the real name are stored
		mylegend=open(sampleResearch)

		for l in mylegend.readlines():
			l= l.rstrip()
			if l[0]=='#':
				continue
			l= l.split(',')
			for a in my_file:
				sourceFile=a
				temp = os.path.basename(a)
				a = temp
				a=a.split('_')
				if ((a[0]==l[9]) and (a[1]==l[7])  and (a[2]==l[8])):
					destFile=l[9]+'''_'''+l[1]+'''_'''+l[2]+'''_'''+l[3]+'''_'''+l[4]+'''_'''+l[6]+'''_grouped_IS'''
					copyfile(sourceFile,os.path.join(dwdcFilterNo,destFile))


	inputFiles = []
	for f in fnmatch.filter(os.listdir(dwdcFilterNo),'*_grouped_IS'):
		inputFiles.append(os.path.join(dwdcFilterNo,f))

	dwdcFilterNo=os.path.join(outputDir,"filterNo","db",suffix)

	if not os.path.exists(dwdcFilterNo):
		os.makedirs(dwdcFilterNo)

	outputFiles = []
	for f in fnmatch.filter(os.listdir(dwdcFilterNo),'*'):
		outputFiles.append(os.path.join(dwdcFilterNo,f))

	print inputFiles

	print outputFiles

	check = compareFilesTime(inputFiles,outputFiles)

	print check

	if check:

		def call_script_R(wd,suffix):
			mycmd='''Rscript --vanilla ''' +os.path.join(utilsDir,'''collisionTable.R ''')+wd+''' '''+suffix
			print mycmd
			subprocess.call(mycmd,shell=True)

		t1 = Thread(target=call_script_R, args=(os.path.join(outputDir,"filterNo","db"),suffix))

		t1.start()
		t1.join()

		my_legend=open(chrList)
	#
		myleg={}
	#
		threads = []

		if (VectorMask==""):

			threads = []
			def call_script_filterMap(filename):
				mycmd='''bedtools closest -a '''+filename+ ''' -b  '''+sortedKnownGene+''' -D "ref" -g '''+genomeSorted+''' > '''+filename+'''_closest_knownGenesV19.txt'''
				print mycmd
				subprocess.call(mycmd,shell=True)

			for filename in os.listdir(dwdcFilterNo):
				if fnmatch.fnmatch(filename,"*_BedFormat"):
					t1 = Thread(target=call_script_filterMap, args=(os.path.join(dwdcFilterNo,filename),))
					threads.append(t1)


			for i in xrange(0, len(threads), NT):
				# Start all threads
				for x in threads[i:i+NT]:
					x.start()
	 			# Wait for all of them to finish
				for x in threads[i:i+NT]:
					x.join()

		if (vectorBed!=""):
			threads = []
			def call_script_filterMap(filename):
				mycmd='''bedtools intersect -v -a '''+filename+''' -b '''+os.path.join(utilsRef,vectorBed)+''' -wa > '''+filename[:-4]+'''_VectMask.txt'''
				print mycmd
				subprocess.call(mycmd,shell=True)

			my_file=[]
			for filename in os.listdir(dwdcFilterNo):
				if fnmatch.fnmatch(filename,"*_closest_knownGenes*.txt"):
					t1 = Thread(target=call_script_filterMap, args=(os.path.join(dwdcFilterNo,filename),))
					threads.append(t1)

			for i in xrange(0, len(threads), NT):
				# Start all threads
				for x in threads[i:i+NT]:
					x.start()
		 		# Wait for all of them to finish
				for x in threads[i:i+NT]:
					x.join()

		if (VectorMask==''):
			for filename in os.listdir(dwdcFilterNo):
				if fnmatch.fnmatch(filename,"*_closest_knownGenes*.txt"):
					shutil.copy(os.path.join(dwdcFilterNo,filename),os.path.join(dwdcFilterNo,filename[:-4]+'_VectMask.txt'))

			threads = []
			def call_script_filterMap(filename):
				mycmd='''cut -f 1,2,4,5,6,10,11,12,13 '''+filename+''' | sort -u > '''+filename[:-4]+'''_uniq.txt'''
				print mycmd
				subprocess.call(mycmd,shell=True)

			my_file=[]
			for filename in os.listdir(dwdcFilterNo):
				if fnmatch.fnmatch(filename,"*_closest_knownGenes*_VectMask.txt"):
					t1 = Thread(target=call_script_filterMap, args=(os.path.join(dwdcFilterNo,filename),))
					threads.append(t1)

			for i in xrange(0, len(threads), NT):
				# Start all threads
				for x in threads[i:i+NT]:
					x.start()
	 			# Wait for all of them to finish
				for x in threads[i:i+NT]:
					x.join()

		def call_script_R(wd,suffix):
			mycmd='''Rscript --vanilla ''' +os.path.join(utilsDir,'''ISAnnotation.R ''')+wd+''' '''+suffix
			print mycmd
			subprocess.call(mycmd,shell=True)

		t1 = Thread(target=call_script_R, args=(dwdcFilterNo,suffix))

		t1.start()
		t1.join()

def getVectorIsTable(outputDir,sampleResearch,utilsDir,utilsRef,chrList,sampleName,sortedKnownGene,genomeSorted,vectorBed,suffix,NT,VectorMask):

	####### Rename files ###################
	dwdcFilterNo=os.path.join(outputDir,"filterNo","db")
	print "collision Folder filterNo:" +dwdcFilterNo

	if not os.path.exists(dwdcFilterNo):
		os.makedirs(dwdcFilterNo)

	inputFiles = []
	for f in fnmatch.filter(os.listdir(outputDir),'*filterNo_grouped_IS.txt'):
		inputFiles.append(os.path.join(outputDir,f))

	outputFiles = []
	for f in fnmatch.filter(os.listdir(dwdcFilterNo),'*_grouped_IS'):
		outputFiles.append(os.path.join(dwdcFilterNo,f))

	print inputFiles
	print outputFiles

	check = compareFilesTime(inputFiles,outputFiles)

	print check

	if check:
		my_file=[]
		for filename in os.listdir(outputDir):
			if fnmatch.fnmatch(filename,'*filterNo_grouped_IS.txt'):
				my_file.append(os.path.join(outputDir,filename))

			#Open a legend were the real name are stored
		mylegend=open(sampleResearch)

		for l in mylegend.readlines():
			l= l.rstrip()
			if l[0]=='#':
				continue
			l= l.split(',')
			for a in my_file:
				sourceFile=a
				temp = os.path.basename(a)
				a = temp
				a=a.split('_')
				if ((a[0]==l[9]) and (a[1]==l[7])  and (a[2]==l[8])):
					destFile=l[9]+'''_'''+l[1]+'''_'''+l[2]+'''_'''+l[3]+'''_'''+l[4]+'''_'''+l[6]+'''_grouped_IS'''
					copyfile(sourceFile,os.path.join(dwdcFilterNo,destFile))


	inputFiles = []
	for f in fnmatch.filter(os.listdir(dwdcFilterNo),'*_grouped_IS'):
		inputFiles.append(os.path.join(dwdcFilterNo,f))

	dwdcFilterNo=os.path.join(outputDir,"filterNo","db",suffix)

	if not os.path.exists(dwdcFilterNo):
		os.makedirs(dwdcFilterNo)

	outputFiles = []
	for f in fnmatch.filter(os.listdir(dwdcFilterNo),'*'):
		outputFiles.append(os.path.join(dwdcFilterNo,f))

	print inputFiles

	print outputFiles

	check = compareFilesTime(inputFiles,outputFiles)

	print check

	if check:

		def call_script_R(wd,suffix):
			mycmd='''Rscript --vanilla ''' +os.path.join(utilsDir,'''collisionTableVector.R ''')+wd+''' '''+suffix
			print mycmd
			subprocess.call(mycmd,shell=True)

		t1 = Thread(target=call_script_R, args=(os.path.join(outputDir,"filterNo","db"),suffix))

		t1.start()
		t1.join()

		my_legend=open(chrList)
	#
		myleg={}
	#
		threads = []

		if (VectorMask==""):

			threads = []
			def call_script_filterMap(filename):
				mycmd='''bedtools closest -a '''+filename+ ''' -b  '''+sortedKnownGene+''' -D "ref" -g '''+genomeSorted+''' > '''+filename+'''_closest_knownGenesV19.txt'''
				print mycmd
				subprocess.call(mycmd,shell=True)

			for filename in os.listdir(dwdcFilterNo):
				if fnmatch.fnmatch(filename,"*_BedFormat"):
					t1 = Thread(target=call_script_filterMap, args=(os.path.join(dwdcFilterNo,filename),))
					threads.append(t1)


			for i in xrange(0, len(threads), NT):
				# Start all threads
				for x in threads[i:i+NT]:
					x.start()
	 			# Wait for all of them to finish
				for x in threads[i:i+NT]:
					x.join()

		if (vectorBed!=""):
			threads = []
			def call_script_filterMap(filename):
				mycmd='''bedtools intersect -v -a '''+filename+''' -b '''+os.path.join(utilsRef,vectorBed)+''' -wa > '''+filename[:-4]+'''_VectMask.txt'''
				print mycmd
				subprocess.call(mycmd,shell=True)

			my_file=[]
			for filename in os.listdir(dwdcFilterNo):
				if fnmatch.fnmatch(filename,"*_closest_knownGenes*.txt"):
					t1 = Thread(target=call_script_filterMap, args=(os.path.join(dwdcFilterNo,filename),))
					threads.append(t1)

			for i in xrange(0, len(threads), NT):
				# Start all threads
				for x in threads[i:i+NT]:
					x.start()
		 		# Wait for all of them to finish
				for x in threads[i:i+NT]:
					x.join()

		if (VectorMask==''):
			for filename in os.listdir(dwdcFilterNo):
				if fnmatch.fnmatch(filename,"*_closest_knownGenes*.txt"):
					shutil.copy(os.path.join(dwdcFilterNo,filename),os.path.join(dwdcFilterNo,filename[:-4]+'_VectMask.txt'))

			threads = []
			def call_script_filterMap(filename):
				mycmd='''cut -f 1,2,4,5,6,10,11,12,13 '''+filename+''' | sort -u > '''+filename[:-4]+'''_uniq.txt'''
				print mycmd
				subprocess.call(mycmd,shell=True)

			my_file=[]
			for filename in os.listdir(dwdcFilterNo):
				if fnmatch.fnmatch(filename,"*_closest_knownGenes*_VectMask.txt"):
					t1 = Thread(target=call_script_filterMap, args=(os.path.join(dwdcFilterNo,filename),))
					threads.append(t1)

			for i in xrange(0, len(threads), NT):
				# Start all threads
				for x in threads[i:i+NT]:
					x.start()
	 			# Wait for all of them to finish
				for x in threads[i:i+NT]:
					x.join()

		def call_script_R(wd,suffix):
			mycmd='''Rscript --vanilla ''' +os.path.join(utilsDir,'''ISAnnotation.R ''')+wd+''' '''+suffix
			print mycmd
			subprocess.call(mycmd,shell=True)

		t1 = Thread(target=call_script_R, args=(dwdcFilterNo,suffix))

		t1.start()
		t1.join()

def calculateDiversity(finalParseFile,umi,chrListFile,ISoutDir,outDir,analysisType):

	myblast=open(finalParseFile)
	my_legend=open(chrListFile)
	dwdt=ISoutDir

	os.chdir(dwdt)

	baseDir = os.path.dirname(finalParseFile)
	baseFile = os.path.basename(finalParseFile)
	outDir = outDir

	if not os.path.exists(outDir):
		os.makedirs(outDir)

	my_report=open(os.path.join(outDir,baseFile[:-4]+'_report.txt'),'w')
	NonGrouped=open(os.path.join(outDir,baseFile[:-4]+'_NonGrouped.txt'),'w')
	#final_IS= open(os.path.join(outDir,baseFile[:-4]+'_final_IS.txt'),'w')
	grouped_IS=open(os.path.join(outDir,baseFile[:-4]+'_grouped_IS.txt'),'w')

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
	mydictUmi={}

	#List all the file containing the reads per IS
	for f in os.listdir(dwdt):
		if fnmatch.fnmatch(f,'*LTR*LC*filter*_chr*.txt'):
			list_pos.append(os.path.join(dwdt,f))

	#ISfiles=open(os.path.join(outDir,'ISfiles.txt'),'w')

	#for f in list_pos:
	#	ISfiles.write(f+'\n')
	#ISfiles.close()

	#Read the input file (final_parse) and store the infos about strand, chr and pos
	for l in myblast.readlines():
		l= l.rstrip()
		if l[0]=='#':
			continue
		l= l.split('\t')
		tmp= l[4], l[5]
		#print "tmp:"+str(tmp)
		try:
			mydict[tmp].append(l[0])
			mydictUmi[tmp].append(l[1])
		except:
			mydict[tmp]=[l[0]]
			mydictUmi[tmp]=[l[1]]

	#print "myDictLen:"+str(len(mydict))

	#for chr,pos in mydict:
	#	print chr,pos,mydict[chr,pos]
	#	print str(set(mydict[chr,pos]))[6:-3]
	#	print str(len(mydict[chr,pos]))

	#print "mydictUmi"
	for chr,pos in mydictUmi:
		#print chr,pos,mydictUmi[chr,pos]
		umiList=[]
		for x in mydictUmi[chr,pos]:
			if x in umi:
				umiList.append(str(umi[x]).strip("['']"))
		umiCount = len(umiList)
		umiUniqCount = len(set(umiList))
		#print umiList
		#print umiCount
		#print set(umiList)
		#print umiUniqCount
		mydictUmi[chr,pos]=umiUniqCount

	#for chr,pos in mydictUmi:
	#	print chr,pos,mydictUmi[chr,pos]

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
	finalParseFileTmp=re.split('\.|\_', baseFile)
	print finalParseFileTmp[2]+'\t'+finalParseFileTmp[4]+'\t'+finalParseFileTmp[7]+'\n'

	for element in list_pos:
		element=re.split('\.|\_', element)
		#print 'elementLen:'+str(len(element))
		#print element[2]+'\n'
		#Finds the correct LTR and LC barcode
		if element[2]==finalParseFileTmp[2] and element[4]==finalParseFileTmp[4]:
			#Finds the correct filter: filterNo, filter30, filter45, filter60
			if element[5]==finalParseFileTmp[7]:
				print element[2]+'\t'+element[4]+'\t'+element[5]+'__'+finalParseFileTmp[2]+'\t'+finalParseFileTmp[4]+'\t'+finalParseFileTmp[7]+'\n'
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
											#print chr, pos, code, '-', str(len(mydict[chr,pos])
											list_strand.append([chr, pos, code, '-', str(len(mydict[chr,pos]))])
							else:
								for name in myleg:
									if (element[6]+'_'+element[7])==name:
										for code in myleg[name]:
											#print chr, pos, code, '+', str(len(mydict[chr,pos]))
											list_strand.append([chr, pos, code, '-', str(len(mydict[chr,pos]))])
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
											list_strand.append([chr, pos, code, '-',str(len(mydict[chr,pos]))])



	#final_IS= open(os.path.join(outDir,baseFile[:-4]+'_list_strand.txt'),'w')

	#print(*list_strand, sep = "\n")
	#if analysis == read:
		#print "read_based"
	#	fin('\n'.join(map(str, list_strand)))

	if analysisType == 'umi':

		for x in list_strand:
			chr = x[0].strip("'")
			pos = x[1].strip("'")
			code = x[2]
			strand = x[3]
			for y in mydictUmi:
				if  y[0] == chr and y[1]==pos:
					umiCount=str(mydictUmi[y])
			x[4]=umiCount
		#print 'umi_based'

	#final_IS.write('\n'.join(map(str, list_strand)))
	#final_IS.close()

	#Sorts the list by chr and pos
	for item in sorted(list_strand, key=lambda x: (int(x[2]), int(x[1]))):
		try:
			myData[item[2]].append(int(item[1]))
		except:
			myData[item[2]]=[int(item[1])]

	#print "myData"
	#for x, y in myData.items():
  	#	print(x, y)
	#	print str(len(myData[x]))

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

	#print('\n'.join(map(str, list_diff)))

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

	#print "my_pos"
	#for x,y in my_pos.items():
	#	print(x,y)

	for a in my_pos:
		b=my_pos[a]
		#print dcondition
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

		#print ("my_array:\n ", my_array)

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
						print str(my_array[0][count][valori])+'\t'+str(my_array[1][count][valori])+'\t'+str(d[0][2])+'\t'+str(my_array[3][count][valori])+'\t'+str(my_array[2][count][valori])
						print '\n'
				count+=1

	#print len(list_item3)
						#
	for c in list_item3:
		NonGrouped.write(str(c[0])+'\t'+str(c[1])+'\t'+str(c[2])+'\t'+str(c[3])+'\t'+str(c[4])+'\n')
		final_list.append([str(c[0]), str(c[1]), str(c[2]), str(c[3]), str(c[4])])

	#print len(final_list)
	#for y in final_list:
	#	final_IS.write('\t'.join(y)+'\n')

	for y in sorted(final_list, key=lambda x: (int(x[2]), int(x[1]))):
		grouped_IS.write('\t'.join(y)+'\n')

	#
	myblast.close()
	my_legend.close()
	my_report.close()
	#

def align2vector(seqPlat,R1Out,R2Out,outputDir,sampleResearch,dirGenome,sampleName):

	### getCollisionTable function #########
	alignOut = os.path.join(outputDir,"align")
	if not os.path.exists(alignOut):
		os.makedirs(alignOut)

	OutFqCount = os.path.join(outputDir,"countFromFq")
	if not os.path.exists(OutFqCount):
		os.makedirs(OutFqCount)

	inputFiles = []
	for f in fnmatch.filter(os.listdir(R1Out),'*.fq_trimwithCutAdapt'):
		inputFiles.append(os.path.join(R1Out,f))

	for f in fnmatch.filter(os.listdir(R2Out),'*.fq_trimwithCutAdapt'):
		inputFiles.append(os.path.join(R2Out,f))

	outputFiles = []

	for f in fnmatch.filter(os.listdir(alignOut),'*_aligned_mem.sam'):
		outputFiles.append(os.path.join(alignOut,f))

	print inputFiles

	print outputFiles

	check = compareFilesTime(inputFiles,outputFiles)

	print check

	if check:
		if (seqPlat=="MiSeq"):
			pattern = re.compile('^@M')
		if (seqPlat=="HiSeq"):
			pattern = re.compile('^@HISEQ')

		with open(sampleResearch) as barCodeFile:
			for barCodeline in barCodeFile:
				barCodeline=barCodeline.rstrip()
				lineSplitBrc=barCodeline.split(',')

				#print lineSplitBrc[9]
				#print sampleName

				if (lineSplitBrc[9]==sampleName):
					r1filename=os.path.join(R1Out,"R1_fastq_trim12nt_qcTrimmed_MatchBlastLtrLc_Barcode_"+lineSplitBrc[7]+".fq_trimwithCutAdapt")
					r2filename=os.path.join(R2Out,"R2_fastq_trim12nt_qcTrimmed_MatchBlastLtrLc_Barcode_"+lineSplitBrc[8]+".fq_trimwithCutAdapt")
					outName=os.path.join(outputDir,"R1_R2_Barcode_"+lineSplitBrc[7]+"_"+lineSplitBrc[8]+"_trimmedID")

					r1ID=[]
					r2ID=[]
					r12Int=[]
					if (os.path.isfile(r1filename) and os.path.isfile(r2filename)):
						with open(r1filename) as r1infile:
							for r1line in r1infile:
								if pattern.search(r1line):
									lineSplit=r1line.split(' ')
									r1ID.append(lineSplit[0][1:])
									#print lineSplit[0][1:]
						with open(r2filename) as r2infile:
							for r2line in r2infile:
								if pattern.search(r2line):
									lineSplit=r2line.split(' ')
									r2ID.append(lineSplit[0][1:])

						r12Int=set(r1ID).intersection(r2ID)
						with open(outName, 'w') as outfile:
							outfile.write("\n".join(r12Int))

						mycmd='''seqtk subseq '''  +r1filename+ ''' ''' +outName+ ''' > '''  +r1filename+ '''_ReadyToAlign'''
						print mycmd
						process = subprocess.Popen(mycmd,shell=True)
						process.wait()

						mycmd='''fastqutils sort ''' +r1filename+ '''_ReadyToAlign > ''' +r1filename+ '''_ReadyToAlignSort'''
						print mycmd
						process = subprocess.Popen(mycmd,shell=True)
						process.wait()

						mycmd='''seqtk subseq '''  +r2filename+ ''' ''' +outName+ ''' > '''  +r2filename+ '''_ReadyToAlign'''
						print mycmd
						process = subprocess.Popen(mycmd,shell=True)
						process.wait()

						mycmd='''fastqutils sort ''' +r2filename+ '''_ReadyToAlign > ''' +r2filename+ '''_ReadyToAlignSort'''
						print mycmd
						process = subprocess.Popen(mycmd,shell=True)
						process.wait()

						mycmd='''wc -l '''+r1filename+ '''_ReadyToAlignSort > '''+os.path.join(OutFqCount,"R1_R2_Barcode_"+lineSplitBrc[7]+"_"+lineSplitBrc[8]+"r1_FqCount.txt")
						print mycmd
						process = subprocess.Popen(mycmd,shell=True)
						process.wait()

						mycmd='''wc -l '''+r2filename+ '''_ReadyToAlignSort > '''+os.path.join(OutFqCount,"R1_R2_Barcode_"+lineSplitBrc[7]+"_"+lineSplitBrc[8]+"r2_FqCount.txt")
						print mycmd
						process = subprocess.Popen(mycmd,shell=True)
						process.wait()

						mycmd='''bwa-mem mem -t 8 '''+dirGenome+''' '''+r1filename+ '''_ReadyToAlignSort ''' +r2filename+ '''_ReadyToAlignSort > ''' +os.path.join(alignOut,"R1_R2_Barcode_"+lineSplitBrc[7]+"_"+lineSplitBrc[8]+"_aligned_mem.sam")
						print mycmd
						process = subprocess.Popen(mycmd,shell=True)
						process.wait()


if __name__ == "__main__":
	main()
