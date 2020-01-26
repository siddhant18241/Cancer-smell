#Histograms
#remove first 2 lines from the input file
#Header should have names of ORs
#Once done follow this code.

file<-read.csv("Log_transformed_TPM_expression.csv")
x<-file$OR10H1
x<-na.omit(x)
h<-hist(x, breaks=10, col="red", xlab="Log2(TPM+1)",main="OR10H1",ylim=c(0,2000))
xfit<-seq(min(x),max(x),length=40)
yfit<-dnorm(xfit,mean=mean(x),sd=sd(x))
yfit <- yfit*diff(h$mids[1:2])*length(x)
lines(xfit, yfit, col="blue", lwd=5)

x<-file$OR5AS1
x<-na.omit(x)
h<-hist(x, breaks=10, col="red", xlab="Log2(TPM+1)",main="OR5AS1",ylim=c(0,80))
xfit<-seq(min(x),max(x),length=40)
yfit<-dnorm(xfit,mean=mean(x),sd=sd(x))
yfit <- yfit*diff(h$mids[1:2])*length(x)
lines(xfit, yfit, col="blue", lwd=5)

x<-file$OR7D2
x<-na.omit(x)
h<-hist(x, breaks=10, col="red", xlab="Log2(TPM+1)",main="OR7D2",ylim=c(0,2000))
xfit<-seq(min(x),max(x),length=40)
yfit<-dnorm(xfit,mean=mean(x),sd=sd(x))
yfit <- yfit*diff(h$mids[1:2])*length(x)
lines(xfit, yfit, col="blue", lwd=5)

x<-file$OR4F17
x<-na.omit(x)
h<-hist(x, breaks=10, col="red", xlab="Log2(TPM+1)",main="OR4F17",ylim=c(0,250))
xfit<-seq(min(x),max(x),length=40)
yfit<-dnorm(xfit,mean=mean(x),sd=sd(x))
yfit <- yfit*diff(h$mids[1:2])*length(x)
lines(xfit, yfit, col="blue", lwd=5)

x<-file$OR2M3
x<-na.omit(x)
h<-hist(x, breaks=10, col="red", xlab="Log2(TPM+1)",main="OR2M3",ylim=c(0,30))
xfit<-seq(min(x),max(x),length=40)
yfit<-dnorm(xfit,mean=mean(x),sd=sd(x))
yfit <- yfit*diff(h$mids[1:2])*length(x)
lines(xfit, yfit, col="blue", lwd=5)

x<-file$OR6C75
x<-na.omit(x)
h<-hist(x, breaks=10, col="red", xlab="Log2(TPM+1)",main="OR6C75",ylim=c(0,30))
xfit<-seq(min(x),max(x),length=40)
yfit<-dnorm(xfit,mean=mean(x),sd=sd(x))
yfit <- yfit*diff(h$mids[1:2])*length(x)
lines(xfit, yfit, col="blue", lwd=5)

x<-file$OR1A1
x<-na.omit(x)
h<-hist(x, breaks=10, col="red", xlab="Log2(TPM+1)",main="OR1A1",ylim=c(0,60))
xfit<-seq(min(x),max(x),length=40)
yfit<-dnorm(xfit,mean=mean(x),sd=sd(x))
yfit <- yfit*diff(h$mids[1:2])*length(x)
lines(xfit, yfit, col="blue", lwd=5)