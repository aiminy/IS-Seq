
args = commandArgs(trailingOnly=TRUE)
if (length(args)<2) {
  stop("Two argument first working directory and suffix name for file(date)", call.=FALSE)
} else{
  # default output file
  wd=args[1]
  suffixCollFile=args[2]
  previous_grouped_IS_folder=args[3]
}

library(plyr)
library(reshape)

if(previous_grouped_IS_folder=='nothing'){files=list.files(path=wd,pattern = "_grouped_IS$",full.names = TRUE)}else{
  files=list.files(path=wd,pattern = "_grouped_IS$",full.names = TRUE)
  previous.files=list.files(path=previous_grouped_IS_folder,pattern = "_grouped_IS$",full.names = TRUE)
  files <- c(previous.files,files)
}

mergeCutoff=7

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


GetIs4OneIs <- function(allsoB,wd,suffixCollFile) {
  trasdOverallReadsC=data.frame(trasd=unique(allsoB$trasd),expr=NA)
  for(ll in 1:dim(trasdOverallReadsC)[1]){
    trasdOverallReadsC$expr[ll]=sum(allsoB$readsCount[allsoB$trasd==trasdOverallReadsC$trasd[ll]])
  }
  allsoB$relReadsCount=NA
  for(i in 1:dim(allsoB)[1]){
    allsoB$relReadsCount[i]=allsoB$readsCount[i]/ trasdOverallReadsC$expr[trasdOverallReadsC$trasd==allsoB$trasd[i]]
  }
  
  allsoOK=allsoB
  
  if(!(exists("allsoBOK"))){allsoBOK=allsoOK}else{allsoBOK=rbind(allsoBOK,allsoOK)}
  
  allsoBOK$label=apply(allsoBOK,1,function(x) {paste(as.character(x[1]),as.numeric(x[2],as.character(x[4])),sep='_',collapse='_')})
  
  out = file.path(wd,suffixCollFile)
  if(!dir.exists(out)){dir.create(out,recursive = TRUE)}
  
  
  trasdAll=unique(allsoBOK$trasd)
  trasdSpec=allsoBOK
  allsoUnique=trasdSpec[,c(13,1,2,3,4)]
  
  
  trasdAllTimes=sort(unique(trasdSpec$time))
  trasdAllSource=sort(unique(trasdSpec$source))
  trasdAllSampleType=sort(unique(trasdSpec$sampleType))
  for(Time in trasdAllTimes){
    for(Source in trasdAllSource){
      for(SampleType in trasdAllSampleType){
        trasdSpecCurr=trasdSpec[(trasdSpec$time==Time & trasdSpec$source==Source &  trasdSpec$sampleType==SampleType),c(13,5)]
        if(dim(trasdSpecCurr)[1]>0){
          trasdSpecCurr=aggregate(trasdSpecCurr['readsCount'], by=trasdSpecCurr['label'], sum)
          if(length(trasdSpecCurr$label)>0){
            coln=c( colnames(allsoUnique),paste(SampleType,Source,Time,sep="_"))
            allsoUniqueM=merge(allsoUnique, trasdSpecCurr, by.x ='label', by.y = 'label', all.x = T, all.y = F)
            colnames(allsoUniqueM)=coln
            allsoUnique=allsoUniqueM
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

if(length(unique(allsoB$trasd))==1){
  GetIs4OneIs(allsoB,wd,suffixCollFile)
}else{
  
  trasdOverallReadsC=data.frame(trasd=unique(allsoB$trasd),expr=NA)
  for(ll in 1:dim(trasdOverallReadsC)[1]){
    trasdOverallReadsC$expr[ll]=sum(allsoB$readsCount[allsoB$trasd==trasdOverallReadsC$trasd[ll]])
  }
  allsoB$relReadsCount=NA
  for(i in 1:dim(allsoB)[1]){
    allsoB$relReadsCount[i]=allsoB$readsCount[i]/ trasdOverallReadsC$expr[trasdOverallReadsC$trasd==allsoB$trasd[i]]
  }
  
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
        print(trasdReadsC)
        if(dim(trasdReadsC)[1]>=2&&(!is.na(trasdReadsC$expr[1])&&(!is.na(trasdReadsC$expr[2])))){
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
  
}