library(rtracklayer)

y <- "hg38"
mySession = browserSession("UCSC")
genome(mySession) <- y

query <- ucscTableQuery(mySession, "rmsk")

rmsk.table <- track(query)

x <- '/Users/aiminyan'
w <- file.path(x,paste0('repeatMasker',y,'BED'))

con <- file(w)
export(rmsk.table,con,format = "bed",trackLine = FALSE)