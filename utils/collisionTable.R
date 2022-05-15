
args = commandArgs(trailingOnly=TRUE)
if (length(args)<2) {
  stop("Two argument first working directory and suffix name for file(date)", call.=FALSE)
} else{
  # default output file
  wd=args[1]
  suffixCollFile=args[2]
  previous_grouped_IS_folder=args[3]
}

#wd="/media/danilopellin/6d51b13b-052e-4544-a57e-554a1394026e/pipeDfci/collision/Lenti_Human/filter60/dbTMP"
#suffixCollFile="31Oct2018"

#length_chr_file <-'/home/ayan/Aimin/ispipe/utilsRefData/hg38/hg38.genome.sorted.txt'

#length_chr_file_data <- read.table(length_chr_file,header = FALSE)
#colnames(length_chr_file_data) <- c('chr_n','length')

library(plyr)
#wd='/home/danilopellin/projects/pipeDfci/collision/Retro_Human/filter60/db'
#suffixCollFile='20Oct2017'
library(reshape)
#setwd(wd)
#setwd("/home/danilopellin/projects/pipeDfci/collision/Lenti_Human/filter60/BOS-LV1")
#setwd("/home/danilopellin/HSR/MLD07")
#files=dir(pattern = "_grouped_IS$")
if(previous_grouped_IS_folder=='nothing'){files=list.files(path=wd,pattern = "_grouped_IS$",full.names = TRUE)}else{
  files=list.files(path=wd,pattern = "_grouped_IS$",full.names = TRUE)
  previous.files=list.files(path=previous_grouped_IS_folder,pattern = "_grouped_IS$",full.names = TRUE)
  files <- c(previous.files,files)
}

#suffixCollFile="10Feb2017"
# length_chr=data.frame(
#   chr_int=c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24),
#   length=c(249250621,243199373,198022430,191154276,180915260,171115067,159138663,146364022,141213431,135534747,135006516,133851895,115169878,107349540,102531392,90354753,81195210,78077248,59128983,63025520,48129895,51304566,155270560,59373566),
#   chr_n=c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX","chrY"),
#   chr_per_hist=c("chr01","chr02","chr03","chr04","chr05","chr06","chr07","chr08","chr09","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX","chrY")
# )
mergeCutoff=7

#length_chr_backup <- length_chr

#length_chr$length <- length_chr_file_data[match(length_chr$chr_n,length_chr_file_data$chr_n),]$length

for(f in 1:length(files)){
  if(file.info(files[f])$size>0){
    so=read.table(files[f],sep = "\t")
    ann=data.frame(t(unlist(strsplit(basename(files[f]),"_"))))
    anndf <- ann[rep(row.names(ann), dim(so)[1]), 1:6]
    so_ann=(cbind(so,anndf))
    if (!(exists("allsoB"))){allsoB=so_ann}  else {allsoB=rbind(allsoB,so_ann)}
  }
}
colnames(allsoB)=c('chr','pos','chrInt','strand','readsCount','library','ptDonor','trasd','source','sampleType','time')
allsoB=sort_df(allsoB,vars = c('chr','pos'))

trasdOverallReadsC=data.frame(trasd=unique(allsoB$trasd),expr=NA)
for(ll in 1:dim(trasdOverallReadsC)[1]){
  trasdOverallReadsC$expr[ll]=sum(allsoB$readsCount[allsoB$trasd==trasdOverallReadsC$trasd[ll]])
}
allsoB$relReadsCount=NA
for(i in 1:dim(allsoB)[1]){
  allsoB$relReadsCount[i]=allsoB$readsCount[i]/ trasdOverallReadsC$expr[trasdOverallReadsC$trasd==allsoB$trasd[i]]
}

#allsoB=allsoB[allsoB$chr=="chr6" & allsoB$pos>76000000 & allsoB$pos<100000000,]
#chr8	8866487	8	+
  
for(St in c("+","-")){
  allso=allsoB[allsoB$strand==St,]
  allso$diff=abs(diff(c(0,allso$pos)))
  min7=which(allso$diff<=mergeCutoff)
  toInv=sort(c(min7,min7-1))
  toInv=unique(toInv)

  allsoToVer=allso[toInv,]

  allsoOK=allso[-toInv,]

  ind=1
  isEv=c(ind)
  while (ind < dim(allsoToVer)[1]){
    if(abs(allsoToVer$pos[ind+1]-allsoToVer$pos[ind])<=mergeCutoff){
      ind=ind+1
      isEv=c(isEv,ind)
    } else{
      allsoToVerisEv=allsoToVer[isEv,]
      allsoToVerisEv$pos=allsoToVerisEv$pos[1]
      trasdReadsC=data.frame(trasd=unique(allsoToVerisEv$trasd),expr=NA)
      for(ll in 1:dim(trasdReadsC)[1]){
        trasdReadsC$expr[ll]=sum(allsoToVerisEv$relReadsCount[allsoToVerisEv$trasd==trasdReadsC$trasd[ll]])
      }
      trasdReadsC=trasdReadsC[order(trasdReadsC$expr,decreasing = T),]
      if(dim(trasdReadsC)[1]>=2){
        if( ((allsoToVerisEv$chr[1]=="chr19") & (allsoToVerisEv$pos[1]==49461738) & (St=='-') ) ){ # force this positive  to assignedISAfterCollision table
          allsoOK=rbind(allsoOK,allsoToVerisEv)
          write.table(allsoToVerisEv,file = file.path(wd,paste("assignedISAfterCollision",suffixCollFile,sep="")),sep = "\t",row.names = F,quote = F,col.names = F,append = T)
        }else if ((trasdReadsC$expr[1]/trasdReadsC$expr[2])>=10){
          allsoOK=rbind(allsoOK,allsoToVerisEv[allsoToVerisEv$trasd==trasdReadsC$trasd[1],])
          write.table(allsoToVerisEv,file = file.path(wd,paste("assignedISAfterCollision",suffixCollFile,sep="")),sep = "\t",row.names = F,quote = F,col.names = F,append = T)
        }else{
          write.table(allsoToVerisEv,file = file.path(wd,paste("deletedISAfterCollision",suffixCollFile,sep="")),sep = "\t",row.names = F,quote = F,col.names = F,append = T)
        }
      }
      if(dim(trasdReadsC)[1]==1){
        allsoOK=rbind(allsoOK,allsoToVerisEv[allsoToVerisEv$trasd==trasdReadsC$trasd[1],])
        write.table(allsoToVerisEv,file = file.path(wd,paste("assignedISAfterCollision",suffixCollFile,sep="")),sep = "\t",row.names = F,quote = F,col.names = F,append = T)
      }

      ind=ind+1
      isEv=c(ind)

    }
  }
  if(!(exists("allsoBOK"))){allsoBOK=allsoOK}
  else{allsoBOK=rbind(allsoBOK,allsoOK)  }
}

allsoBOK$label=apply(allsoBOK,1,function(x) {paste(as.character(x[1]),as.numeric(x[2],as.character(x[4])),sep='_',collapse='_')})

#setwd(paste(wd,"/",suffixCollFile,sep=""))

out = file.path(wd,suffixCollFile)
if(!dir.exists(out)){dir.create(out,recursive = TRUE)}
           
trasdAll=unique(allsoBOK$trasd)
for(ll in 1:length(trasdAll)){
  trasdSpec=allsoBOK[allsoBOK$trasd==trasdAll[ll],]
  unik <- !duplicated(trasdSpec$label)
  allsoUnique=trasdSpec[unik ,c(14,1,2,3,4)]
  trasdAllTimes=sort(unique(trasdSpec$time))
  trasdAllSource=sort(unique(trasdSpec$source))
  trasdAllSampleType=sort(unique(trasdSpec$sampleType))
  for(Time in trasdAllTimes){
    for(Source in trasdAllSource){
      for(SampleType in trasdAllSampleType){
        trasdSpecCurr=trasdSpec[(trasdSpec$time==Time & trasdSpec$source==Source &  trasdSpec$sampleType==SampleType),c(14,5)]
        if(dim(trasdSpecCurr)[1]>0){
          trasdSpecCurr=aggregate(trasdSpecCurr['readsCount'], by=trasdSpecCurr['label'], sum)
          if(length(trasdSpecCurr$label)>0){
            coln=c( colnames(allsoUnique),paste(SampleType,Source,Time,sep="_"))
            allsoUniqueM=merge(allsoUnique, trasdSpecCurr, by.x ='label', by.y = 'label', all.x = T, all.y = F)
            colnames(allsoUniqueM)=coln
            allsoUnique=allsoUniqueM
           # print(dim(allsoUnique))
          }
        }
      }
    }
  }
  allsoUnique=sort_df(allsoUnique,vars = c('chr','pos'))
  allsoUnique= allsoUnique[,-1]
  allsoUnique[is.na(allsoUnique)]=0
  write.table(allsoUnique,file = file.path(out,paste(trasdSpec$library[1],trasdSpec$ptDonor[1],trasdSpec$trasd[1],suffixCollFile,"CollisionClean",sep="_")),sep = "\t",row.names = F,quote = F,col.names = T)
  allsoUniqueBed=data.frame(chr=allsoUnique[,1],start=allsoUnique[,2],end=(allsoUnique[,2]+1),chrInt=allsoUnique[,3],c2=0,strand=allsoUnique[,4])
  allsoUniqueBed=sort_df(allsoUniqueBed,vars = c('chrInt','start'))

  write.table(allsoUniqueBed,file = file.path(out,paste(trasdSpec$library[1],trasdSpec$ptDonor[1],trasdSpec$trasd[1],suffixCollFile,"CollisionClean_BedFormat",sep="_")),sep = "\t",row.names = F,quote = F,col.names = F)

}