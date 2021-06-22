library(plotrix)
library(smacof)
#library(unfolding)
args <- commandArgs(TRUE)
fname <- args[1]
data <- read.table(fname, header = TRUE, sep = ",", check.names = FALSE)
wm <- data.matrix(read.table(args[3],sep = ",",check.names = FALSE), rownames.force = NA)
dim <- as.numeric(args[4])
loo <- as.numeric(args[5])

out <- smacofRect(data,ndim=dim,weightmat=wm,verbose=FALSE,itmax=5000)
stress1 <- round(as.numeric(out[5]),digits=2)
stress1

virus <- paste(fname,'.virus', sep='')
serum <- paste(fname,'.serum', sep='')
error <- paste(fname,'.error', sep='')
stressfile <- paste(fname, '.stress', sep='')
kerror <-  paste(fname, '.kerror', sep='')

write.table(out$confdiss, file = args[2], append = FALSE, quote = FALSE, sep = ",",
            eol = "\n", na = "NA", dec = ".", row.names = TRUE,
            col.names = TRUE, qmethod = c("escape", "double"),
            fileEncoding = "")
if (dim == 2)
{
colnames(out$conf.col) <- c('x','y')
colnames (out$conf.row) <- c('x','y')
}
if (dim == 3)
{
colnames(out$conf.col) <- c('x','y','z')
colnames (out$conf.row) <- c('x','y','z')
}
if (dim == 2 || dim == 3)
{
write.table(out$conf.col, file = virus, append = FALSE, quote = FALSE, sep = ",",
            eol = "\n", na = "NA", dec = ".", row.names = TRUE,
            col.names = TRUE, qmethod = c("escape", "double"),
            fileEncoding = "")
write.table(out$conf.row, file = serum, append = FALSE, quote = FALSE, sep = ",",
            eol = "\n", na = "NA", dec = ".", row.names = TRUE,
            col.names = TRUE, qmethod = c("escape", "double"),
            fileEncoding = "")
}
write(out$stress,file = stressfile, append = FALSE, sep=",")
kruskalerror <- 100*abs(out$obsdiss-out$confdiss)/out$obsdiss
write.table(kruskalerror, file = error, append = FALSE, quote = FALSE, sep = ",",
            eol = "\n", na = "NA", dec = ".", row.names = TRUE,
            col.names = TRUE, qmethod = c("escape", "double"),
            fileEncoding = "")
if(any(is.infinite(kruskalerror))) {
  kruskalerror[is.infinite(kruskalerror)] = NA
}
write(mean(kruskalerror, na.rm=TRUE), file = kerror, append = FALSE, sep = ",")
rl=dim(data)[1]
cl=dim(data)[2]
newmatrix = data
neterror<-c()
if (loo == 1)
{
  print("Computing Loo erros. Please wait ...")
  for (i in 1:rl)
  {
    for (j in 1:cl)
    {
      newmatrix[i,j]=NaN
      #print (newmatrix)
      out <- smacofRect(newmatrix,ndim=dim,weightmat=wm,verbose=FALSE,itmax=5000)
      error<- 100*abs(data[i,j] - out$confdiss[i,j])/data[i,j]
      neterror<-c(neterror,error)
      stress1 <- round(as.numeric(out[5]),digits=2)
      newmatrix[i,j]=data[i,j]
    }
  }
  if(any(is.infinite(neterror))) {
    neterror[is.infinite(neterror)] = NA
  }
loofile <- paste(fname, '.loo', sep='')
write(mean(neterror, na.rm=TRUE), file = loofile, append = FALSE, sep = ",")
  
}

#plot(out, "stressplot")
