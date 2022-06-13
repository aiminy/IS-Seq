GroupIs <- function(is.gr.long,is.gr) {
  
  GetGroupedIs <- function(is.gr) {
    
    is.gr.reduced <- reduce(is.gr,min.gapwidth=7L,with.revmap=T)
    
    TotalCount <- lapply(is.gr.reduced$revmap,function(u){
      x <- sum(mcols(is.gr[u])$est)
    })
    
    is.gr.reduced$TotalCount <-  TotalCount
    
    AveCount <- lapply(is.gr.reduced$revmap,function(u){
      x <- round(mean(mcols(is.gr[u])$est))
    })
    
    is.gr.reduced$AveCount <-  AveCount
    
    MaxCount <- lapply(is.gr.reduced$revmap,function(u){
      x <- max((mcols(is.gr[u])$est))
    })
    
    is.gr.reduced$MaxCount <-  MaxCount
    
    MinCount <- lapply(is.gr.reduced$revmap,function(u){
      x <- min((mcols(is.gr[u])$est))
    })
    
    is.gr.reduced$MinCount <-  MinCount
    
    pos.max <- lapply(is.gr.reduced$revmap,function(u){
      x <-  is.gr[u]
      y <- x[which.max(mcols(x)$est)]
    })
    
    is.gr.reduced$pos.max <- pos.max
    
    pos.min <- lapply(is.gr.reduced$revmap,function(u){
      x <-  is.gr[u]
      y <- x[which.min(mcols(x)$est)]
    })
    
    is.gr.reduced$pos.min <- pos.min
    
    is.gr.reduced
  }
  
  is.gr.long.reduced <- GetGroupedIs(is.gr.long)
  
  GetGroupIs1 <- function(is.gr){
    
    z.redu <- lapply(1:dim(mcols(is.gr))[2],function(u){
      
      z <- is.gr[,u]
      z.reduced <- reduce(z,min.gapwidth=7L,with.revmap=T)
      
      count <- lapply(z.reduced$revmap,function(u){
        x <- sum(mcols(z[u])[,1])
      })
      
      z.reduced$count <- count
      
      pos.max <- lapply(z.reduced$revmap,function(u){
        x <- z[u]
        y <- x[which.max(mcols(x)[,1])]
      })
      
      z.reduced$pos.max <- pos.max
      
      z.reduced
      
    })
    
    names(z.redu) <- colnames(mcols(is.gr))
    z.redu
  }
  
  is.gr.1 <- GetGroupIs1(is.gr)
  
  is.gr.long.reduced.1 <- reduce(do.call(c,is.gr.long.reduced$pos.max))
  
  is.gr.2 <- lapply(1:length(is.gr.1), function(u,is.gr.long.reduced.1){
    z <- is.gr.long.reduced.1
    hits <- findOverlaps(is.gr.long.reduced.1,is.gr.1[[u]])
    z$count <- unlist(is.gr.1[[u]][subjectHits(hits)]$count)
    z$SampleName <- names(is.gr.1)[u]
    z <- as.data.frame(z)
    z
  },is.gr.long.reduced.1)
  
  is.gr.3 <- do.call(rbind,is.gr.2)
  #paste0(seqnames,'_',start,'_',strand)
  
  library(reshape2)
  is.gr.4 <- dcast(is.gr.3[,c(1:2,5:7)], seqnames+start+strand~SampleName, value.var='count',sum)
  
}

