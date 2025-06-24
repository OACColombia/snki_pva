# Formating the original data set for parameter estimation:
# We will assume that the 4,5 and 6 counts constitute each three replicates
# The output of this code consists of 
# a) three matrices, each with 6 + 1 rows (one per pop. + 'other') and
#    as many columns as years the SNKI monitoring has lasted.
# b) three accompannying matrices that have a number one or a zero in each cell
#    indicating whether for that cell (a combination of year and population)
#    an "NA" is in the corresponding count matrix.  These matrices will be 
#    useeful for the JAGS estimation procedure

# Reading the data

rawdat <- read.csv("snailkite counts 1997_2025.csv")
View(rawdat)
names.rawdat <- names(rawdat)
deme.names <- c("EAST","EVER","KRV","OKEE","PP","SJM","OTHER" )
cols.with.counts <- match(deme.names, names.rawdat)


years <- unique(rawdat$year)
nyears <- length(years)
year.col <- which(names.rawdat=="year",arr.in=TRUE)


# I will use a simple time scale from 0 to 28, we can use the 'years' vector
# above to do the plots
tt <- years - years[1]

# List the survey numbers:
surv.num <- unique(rawdat$survey_num)

# Create a list with three empty matrices, one per survey (Surveys 4,5 and 6)

S4mat <- matrix(0, nrow=7, ncol=nyears)
S5mat <- S4mat
S6mat <- S4mat

Counts.list <- list(S4mat = S4mat, S5mat=S5mat, S6mat=S6mat)
NAS.list  <- list(NA4mat = S4mat, NA5mat=S5mat, NA6mat=S6mat)
short.inds.list <- list()

for(i in 1:3){
  
  SN <- surv.num[(i+3)]
  
  ith.rows   <- which(rawdat$survey_num==SN, arr.ind=TRUE)
  
  #ith.datmat <- rawdat[ith.rows, c(year.col,cols.with.counts)]
  ith.datmat <- t(rawdat[ith.rows, cols.with.counts])
  colnames(ith.datmat) <- tt
  
  Counts.list[[i]] <- ith.datmat
  
  ith.isnamat <- 1-is.na(ith.datmat)
  NAS.list[[i]] <- ith.isnamat
  
  ith.shortind.list <- list()
  for(j in 1:7){
    
    ith.shortind.list[[j]] <- which(ith.isnamat[j,]==1,arr.ind=TRUE)
    
  }
  names(ith.shortind.list) <- deme.names
  short.inds.list[[i]] <- ith.shortind.list 
    
}
names(short.inds.list) <- c("short.inds4","short.inds5", "short.inds6")


save.image("SNKIdata.RData")
