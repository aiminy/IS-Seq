# Examples:

# Rscript ~/Aimin/ispipe/R/GetINSPFragMLE0.R /Volumes/share/ISseqOutput/Dec282021/IsaByINSPIIRED/fa /Volumes/share/ISseqOutput/IsSeq_paper_results/HL60/INSP_rev95 95 

# Rscript ~/Aimin/ispipe/R/GetINSPFragMLE0.R /Volumes/share/ISseqOutput/Dec282021/IsaByINSPIIRED/fa /Volumes/share/ISseqOutput/IsSeq_paper_results/HL60/INSP_rev80 80

# Rscript ~/IS-Seq/R/GetFragMLE.R ~/SHARE/Aimin/INSPIIRED_test_output/clone1-1/IS_0/allSites.rds ~/SHARE/Aimin/INSPIIRED_test_output/clone1-1/IS_0

args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} else if (length(args)>=1) {
  # default output file
  input_file=args[1]
  output.dir=args[2]
  sequence.identity.threshold=args[3]
}
if(!dir.exists(output.dir)){dir.create(output.dir,recursive = TRUE)}

#allSites <- readRDS("/Volumes/share/ISseqOutput/INSPIIRED_test_run/pool3-1/allSites.rds")
#sites.final <- readRDS("/Volumes/share/ISseqOutput/INSPIIRED_test_run/pool3-1/sites.final.rds")
#multihitData <- readRDS("/Volumes/share/ISseqOutput/INSPIIRED_test_run/pool3-1/multihitData.rds")
library(GenomicRanges)
GetReadCount <- function(sites) {
  sites.reduced <- flank(sites, -1, start=TRUE)
  sites.reduced <- reduce(sites.reduced, min.gapwidth = 0L, with.revmap=TRUE)
  sites.reduced$counts <- sapply(sites.reduced$revmap, length)
  sites.reduced
}

#GetReadCount(allSites)

GetFragLen <- function(sites) {
  
  sites$FragLen <- width(sites)
  
  sites.reduced <- flank(sites, -1, start=TRUE)
  
  sites.reduced$posid <- paste0(as.character(seqnames(sites.reduced)),as.character(strand(sites.reduced)),start(flank(sites.reduced, width=-1, start=TRUE)))
  
  sites.reduced
}

#allSites.FragLen <- GetFragLen(allSites)

GetFragLenCount <- function(sites) {
  
  posid <- unique(sites$posid)
  
  counLen <- lapply(1:length(posid), function(u,sites){
    
    n <- length(unique(sites[which(sites$posid==posid[u]),]$FragLen))
    
  },sites)
  
  names(counLen) <- posid
  
  res <- data.frame(pos=posid,FragLenCount=do.call(rbind,counLen))
  res
}

#GetFragLenCount(allSites.FragLen)

GetFragCountMLE <- function(allSites.FragLen) {
  estAbund.df <- data.frame(Chromosome=seqnames(allSites.FragLen),Position=start(allSites.FragLen),Ort=strand(allSites.FragLen),length=allSites.FragLen$FragLen)
  estAbund.df <- unique(estAbund.df)
  id.estAbund.df <- with(estAbund.df, paste(Chromosome, Position, Ort))
  library(sonicLength)
  fit.real.INSP <- estAbund(id.estAbund.df, estAbund.df$length,min.length=min(estAbund.df$length))
  res <- fit.real.INSP$theta[unique(id.estAbund.df)]
  res <- as.data.frame(res)
  colnames(res) <- 'MLE'
  res
}

#FragMLE <- GetFragCountMLE(allSites.FragLen)

#input_dir <- '/Volumes/share/ISseqOutput/Dec282021/IsaByINSPIIRED/fa/rev'

#input_dir <- '/Volumes/share/ISseqOutput/Dec282021/IsaByINSPIIRED/fa'

input_dir <- '/home/ubuntu/SHARE/Aimin/INSPIIRED_test_output/clone1-1/IS_0'

input_dir <- '~/SHARE/ISseqOutput/INSPIIRED_test_run_using_IsSeqHg18/clone1-2'
input_dir <- '/home/ubuntu/SHARE/ISseqOutput/INSPIIRED_test_run/clone1-1'


allSites.files <- list.files(path=input_dir,pattern = 'allSites.rds',full.names = T,recursive = T)

if(sequence.identity.threshold==95){
  allSites.files.0 <- allSites.files[-grep('rev80|rev0',allSites.files)]
  allSites.files.1 <- allSites.files.0[grep('rev',allSites.files.0)]
  allSites.files.2 <- allSites.files.1[-1]
  allSites.files <- allSites.files.2
  SampleName <- basename(dirname(allSites.files))
}

if(sequence.identity.threshold==80){
  allSites.files <- allSites.files[grep('rev80',allSites.files)]
  SampleName <- basename(dirname(dirname(allSites.files)))
}

if(sequence.identity.threshold==0){
  allSites.files <- allSites.files[grep('rev0',allSites.files)]
  SampleName <- basename(dirname(dirname(allSites.files)))
}

print(allSites.files)

res.allSites<- lapply(1:length(allSites.files), function(u){
  allSites <- readRDS(allSites.files[u])
})


names(res.allSites) <- SampleName

res.count <- lapply(1:length(allSites.files), function(u){
  allSites <- readRDS(allSites.files[u])
  
  if(length(allSites)!=0){
    allSites.FragLen <- GetFragLen(allSites)
    z <- GetFragLenCount(allSites.FragLen)
    #estAbund.df <- data.frame(Chromosome=seqnames(allSites.FragLen),Position=start(allSites.FragLen),Ort=strand(allSites.FragLen),length=allSites.FragLen$FragLen)
    #z <- data.frame(estAbund.df,replicates=rep(u,dim(estAbund.df)[1]))
  }else{
    z <- NULL
  }
  z
})

#names(res.count) <- basename(dirname(allSites.files))
names(res.count) <- SampleName

res.count.1 <- res.count[!sapply(res.count,is.null)]

res.count <- do.call(rbind,res.count.1)

#print(res.1)

saveRDS(res.count,file.path(output.dir,'FragCount.rds'))

res.mle <- lapply(1:length(allSites.files), function(u){

  # u <- 1
  allSites <- readRDS(allSites.files[u])
  # allSites.FragLen <- GetFragLen(allSites)
  # estAbund.df <- data.frame(Chromosome=seqnames(allSites.FragLen),Position=start(allSites.FragLen),Ort=strand(allSites.FragLen),length=allSites.FragLen$FragLen)
  # estAbund.df <- unique(estAbund.df)
  # id.estAbund.df <- with(estAbund.df, paste(Chromosome, Position, Ort))
  # library(sonicLength)
  # fit.real.INSP <- estAbund(id.estAbund.df, estAbund.df$length,min.length=min(estAbund.df$length))
  
  #z <- GetFragCountMLE(allSites.FragLen)
  
  if(length(allSites)!=0){
    allSites.FragLen <- GetFragLen(allSites)
    z <- GetFragCountMLE(allSites.FragLen)
    # estAbund.df <- data.frame(Chromosome=seqnames(allSites.FragLen),Position=start(allSites.FragLen),Ort=strand(allSites.FragLen),length=allSites.FragLen$FragLen)
    # z <- data.frame(estAbund.df,replicates=rep(u,dim(estAbund.df)[1]))
  }else{
    z <- NULL
  }
  z
})

names(res.mle) <- SampleName

res.mle.1 <- res.mle[!sapply(res.mle,is.null)]

res.mle <- do.call(rbind,res.mle.1)

saveRDS(res.mle,file.path(output.dir,'FragMLE.rds'))

# estAbund.df <- res.1
# id.estAbund.df <- with(estAbund.df, paste(Chromosome, Position, Ort))
# library(sonicLength)
# fit.real.INSP <- estAbund(id.estAbund.df, estAbund.df$length,min.length=min(estAbund.df$length),replicates=estAbund.df$replicates)
# res <- fit.real.INSP$theta[unique(id.estAbund.df)]

res.mle.reformat <- lapply(1:length(res.mle.1), function(u){
  
  z <- data.frame(location=gsub(' ','_',row.names(res.mle.1[[u]])),sampleName=names(res.mle.1)[u],estAbund=round(res.mle.1[[u]]$MLE))
  z
})

res.mle.reformat.1 <- do.call(rbind,res.mle.reformat)
library(reshape2)
data_wide <- dcast(res.mle.reformat.1, location~sampleName, value.var='estAbund',sum)

#print(data_wide)

saveRDS(res.mle.reformat.1,file.path(output.dir,"res.mle.reformat.rds"))
saveRDS(data_wide,file.path(output.dir,"data_wide.rds"))

MakeGrLongGr <- function(res.mle.reformat.1) {
  is.gr.long <- makeGRangesFromDataFrame(data.frame(chr=unlist(lapply(strsplit(res.mle.reformat.1$id,split='_'),'[[',1)),start=unlist(lapply(strsplit(res.mle.reformat.1$id,split='_'),'[[',2)),end=unlist(lapply(strsplit(res.mle.reformat.1$id,split='_'),'[[',2)),strand=unlist(lapply(strsplit(res.mle.reformat.1$id,split='_'),'[[',3)),res.mle.reformat.1[,-1]),keep.extra.column=T)
}

MakeGr <- function(res.mle.reformat.1) {
  data_wide <- dcast(res.mle.reformat.1, id~sample, value.var='est',sum)
  is.gr <- makeGRangesFromDataFrame(data.frame(chr=unlist(lapply(strsplit(data_wide$id,split='_'),'[[',1)),start=unlist(lapply(strsplit(data_wide$id,split='_'),'[[',2)),end=unlist(lapply(strsplit(data_wide$id,split='_'),'[[',2)),strand=unlist(lapply(strsplit(data_wide$id,split='_'),'[[',3)),data_wide[,-1]),keep.extra.column=T)
}

INSP.mle <- res.mle.reformat.1

colnames(INSP.mle) <- c('id','sample','est')

INSP.mle.long.gr <- MakeGrLongGr(INSP.mle)
INSP.mle.is.gr <- MakeGr(INSP.mle)

saveRDS(INSP.mle.long.gr,file.path(output.dir,'INSP.mle.long.gr.rds'))
saveRDS(INSP.mle.is.gr,file.path(output.dir,'INSP.mle.is.gr.rds'))

source('/Users/aiminyan/Aimin/ispipe/R/GroupIsFunction.R')

source('~/ispipe/R/GroupIsFunction.R')

IS.grouped <- GroupIs(INSP.mle.long.gr,INSP.mle.is.gr)

saveRDS(IS.grouped,file.path(output.dir,'IS.grouped.rds'))

IS.grouped.gr <- makeGRangesFromDataFrame(data.frame(chr=IS.grouped$seqnames,start=IS.grouped$start,end=IS.grouped$start,strand=IS.grouped$strand,IS.grouped[,-c(1:3)]),keep.extra.columns = T)

getRelativeAbundance <- function(data.collision) {
  A=S4Vectors::as.matrix(data.collision)
  B <- apply(A,2,function(u){round((u/sum(u))*100,3)})
  B
}

z <- getRelativeAbundance(mcols(IS.grouped.gr))

IS.grouped.gr.relative <- cbind(as.data.frame(granges(IS.grouped.gr)),z)

data.ordered.1 <- data.frame(Is_location=paste0(IS.grouped.gr.relative$seqnames,'_',IS.grouped.gr.relative$start,'_',IS.grouped.gr.relative$strand),IS.grouped.gr.relative[,-c(1:5)])


GetSortedIsByRelAbun <- function(data.ordered.1) {
  trt <- colnames(data.ordered.1)[-1]
  
  data.ordered.L <- lapply(1:length(trt), function(u){
    
    #u <- 1
    index <- which(colnames(data.ordered.1)==trt[u])
    
    y <- data.ordered.1[,c(1,index)]
    
    z <- y[y[,2]!=0,]
    
    w <- data.frame(trt=rep(colnames(z)[2],dim(z)[1]),location=z[,1],value=z[,2])
    
    w <- w[order(w$value,decreasing = T),]
    
    n <- dim(w)[1]
    
    w$value <- as.numeric(w$value)
    
    
    if (n<=10){
      t <- data.frame(trt=unique(w[1:n,]$trt),location='Others',value= sum(w[-c(1:n),]$value))
      ww <- rbind(w[1:n,],t)
    }else
    {
      t <- data.frame(trt=unique(w[1:10,]$trt),location='Others',value= sum(w[-c(1:10),]$value))
      ww <- rbind(w[1:10,],t)
    }
    
    ww
  })
  
  names(data.ordered.L) <- trt
  data.ordered.L
}

#data.ordered.1 <- Is.Rel.gr[[1]]

INSP.top.10.IS <- GetSortedIsByRelAbun(data.ordered.1)

saveRDS(INSP.top.10.IS,file.path(output.dir,'INSP.top.10.IS.rds'))

GetMSE <- function(IsSeq.top.10.IS) {
  
  MSE <- lapply(1:length(IsSeq.top.10.IS), function(u,IsSeq.top.10.IS){
    
    if(names(IsSeq.top.10.IS)[u]=='CL6'){
      #u <- 1
      mse <- sum((IsSeq.top.10.IS[[u]][1,]$value-c(100))^2)/1
    }else if(names(IsSeq.top.10.IS)[u]=='CLH4'){
      mse <- sum((IsSeq.top.10.IS[[u]][1:5,]$value-c(20,20,20,20,20))^2)/5
    }else if(names(IsSeq.top.10.IS)[u]=='MOI30CLB6'){
      mse <- sum((IsSeq.top.10.IS[[u]][1:3,]$value-c(33.33333,33.33333,33.33333))^2)/3
    }else{
      mse <- sum((IsSeq.top.10.IS[[u]][1:2,]$value-c(50,50))^2)/2
    }
    
    mse <- round(mse,2)
    mse
  },IsSeq.top.10.IS)
  
  names(MSE) <- names(IsSeq.top.10.IS)
  MSE <- t(as.data.frame(MSE))
  MSE
}

names(INSP.top.10.IS) <- gsub('\\.','',names(INSP.top.10.IS))
INSP.mse <- GetMSE(INSP.top.10.IS)

print(INSP.mse)

saveRDS(INSP.mse,file.path(output.dir,'INSP.mse.rds'))
save.image(file.path(output.dir,"Results.RData"))

