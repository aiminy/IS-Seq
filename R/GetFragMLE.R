# Examples:

# Rscript ~/IS-Seq/R/GetFragMLE.R ~/SHARE/Aimin/INSPIIRED_test_output/clone1-1/IS_0/allSites.rds clone1-1 ~/SHARE/Aimin/INSPIIRED_test_output/clone1-1/IS_0

initial.options <- commandArgs(trailingOnly = FALSE)
print(initial.options)

file.arg.name <- "--file="
script.name <- sub(file.arg.name, "", initial.options[grep(file.arg.name, initial.options)])
script.dirname <- dirname(script.name)
print(script.dirname)

args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} else if (length(args)>=1) {
  # default output file
  input_file=args[1]
  SampleName=args[2]
  output.dir=args[3]
}

if(!dir.exists(output.dir)){dir.create(output.dir,recursive = TRUE)}

library(GenomicRanges)
GetReadCount <- function(sites) {
  sites.reduced <- flank(sites, -1, start=TRUE)
  sites.reduced <- reduce(sites.reduced, min.gapwidth = 0L, with.revmap=TRUE)
  sites.reduced$counts <- sapply(sites.reduced$revmap, length)
  sites.reduced
}

GetFragLen <- function(sites) {
  
  sites$FragLen <- width(sites)
  sites.reduced <- flank(sites, -1, start=TRUE)
  sites.reduced$posid <- paste0(as.character(seqnames(sites.reduced)),as.character(strand(sites.reduced)),start(flank(sites.reduced, width=-1, start=TRUE)))
  sites.reduced
}

GetFragLenCount <- function(sites) {
  
  posid <- unique(sites$posid)
  
  counLen <- lapply(1:length(posid), function(u,sites){
    
    n <- length(unique(sites[which(sites$posid==posid[u]),]$FragLen))
    
  },sites)
  
  names(counLen) <- posid
  
  res <- data.frame(pos=posid,FragLenCount=do.call(rbind,counLen))
  res
}

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

allSites.files <- input_file
  
print(allSites.files)

res.allSites<- readRDS(allSites.files)
print(res.allSites)
names(res.allSites) <- SampleName

res.count <- lapply(1:length(allSites.files), function(u){
  allSites <- readRDS(allSites.files[u])
  if(length(allSites)!=0){
    allSites.FragLen <- GetFragLen(allSites)
    z <- GetFragLenCount(allSites.FragLen)
  }else{
    z <- NULL
  }
  z
})

names(res.count) <- SampleName
res.count.1 <- res.count[!sapply(res.count,is.null)]
res.count <- do.call(rbind,res.count.1)

saveRDS(res.count,file.path(output.dir,'FragCount.rds'))

res.mle <- lapply(1:length(allSites.files), function(u){
  allSites <- readRDS(allSites.files[u])
  if(length(allSites)!=0){
    allSites.FragLen <- GetFragLen(allSites)
    z <- GetFragCountMLE(allSites.FragLen)
  }else{
    z <- NULL
  }
  z
})

names(res.mle) <- SampleName

res.mle.1 <- res.mle[!sapply(res.mle,is.null)]

res.mle <- do.call(rbind,res.mle.1)

print(res.mle)

saveRDS(res.mle,file.path(output.dir,'FragMLE.rds'))

res.mle.reformat <- lapply(1:length(res.mle.1), function(u){
  z <- data.frame(location=gsub(' ','_',row.names(res.mle.1[[u]])),sampleName=names(res.mle.1)[u],estAbund=round(res.mle.1[[u]]$MLE))
  z
})

res.mle.reformat.1 <- do.call(rbind,res.mle.reformat)

print(res.mle.reformat.1)

library(reshape2)
data_wide <- dcast(res.mle.reformat.1, location~sampleName, value.var='estAbund',sum)

print(data_wide)

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

source(file.path(script.dirname,'GroupIsFunction.R'))

IS.grouped <- GroupIs(INSP.mle.long.gr,INSP.mle.is.gr)

print(IS.grouped)

saveRDS(IS.grouped,file.path(output.dir,'IS.grouped.rds'))

save.image(file.path(output.dir,"Results.RData"))
