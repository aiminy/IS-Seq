load("~/intsitecaller/testCases/intSiteValidation/completeMetadata.RData")
completeMetadata.new <- completeMetadata[,-c(18:19)]

output.dir <- '~/IS-Seq/utilsRefData/INSPIIRED'

if(!dir.exists(output.dir)){dir.create(output.dir,recursive = TRUE)}

con <- file.path(output.dir,'completeMetadata.RData')

save(completeMetadata.new,file = con)
