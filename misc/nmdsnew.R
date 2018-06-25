#library(plotrix)
library(smacof)
#library(unfolding)
args <- commandArgs(TRUE)
fname <- 'data_dist.csv'
data <- read.table(fname, header = TRUE, sep = ",", check.names = FALSE)
wm <- data.matrix(read.table('data_wm.csv',sep = ",",check.names = FALSE), rownames.force = NA)
rl=dim(data)[1]
cl=dim(data)[2]
#wm
out <- smacofRect(data,ndim=2,verbose=TRUE,itmax=100000)
stress1 <- round(as.numeric(out[5]),digits=2)
stress1
#plot (out, type="n", xlab="", ylab="", xaxt="n", yaxt="n", asp=1)
plotname <- paste(fname,".pdf",sep="")
#pdf(plotname)
plot(out, plot.type="resplot")
plot(out, type="p",plot.type = "confplot",main="Neutralization Map of HCV",xlab="Residual Infectivity",ylab="Residual Infectivity",col.rows=4,col.columns=2,pch=1,joint=TRUE,
     label.conf.columns = list(label = TRUE,pos =1, col = 2, cex = 0.75),
     label.conf.rows = list(label = TRUE,pos =1, col = 4, cex = 0.75),cex=1)
#plot(out, type='')
usr <- par("usr")
abline(h=-10:10, v=-10:10, col="gray", lty=3)
legend(usr[2]-1.75,usr[3]+1,inset=0.02,c("Antigen","Antiserum"),cex=0.75,col=c(2,4),pch=1,text.font=3,box.lty=0,bty="n",x.intersp=0.25,y.intersp =0.5,xjust=0)
text(usr[2]-3,usr[4]-0.15,pos=1,labels = paste("Stress ",stress1,sep=": "),cex=1,col=3)
newmatrix = data
neterror<-c()
for (i in 1:rl)
  {
  for (j in 1:cl)
    {
    newmatrix[i,j]=NaN
    #print (newmatrix)
    out <- smacofRect(newmatrix,ndim=2,verbose=FALSE)
    error<- 100*abs(data[i,j] - out$confdiss[i,j])/data[i,j]
    neterror<-c(neterror,error)
    stress1 <- round(as.numeric(out[5]),digits=2)
    print (stress1)
    print(error)
    newmatrix[i,j]=data[i,j]
    }
}

print(mean(neterror, na.rm = TRUE))
#write.table(out[2], file = args[2], append = FALSE, quote = TRUE, sep = "\t",
#            eol = "\n", na = "NA", dec = ".", row.names = TRUE,
#            col.names = TRUE, qmethod = c("escape", "double"),
#            fileEncoding = "")
#plot(out, "stressplot")

