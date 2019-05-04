library(fda)
setwd("/Users/hamzaqureshi/Google Drive/MSC/AARMS 2018 (UPEI)/Functional Analysis for Big Data")

my_data <- read.delim(file.choose())

F2 = my_data[1:34,]

knots=seq(1,10,by=1)
basiss = create.bspline.basis(rangeval=c(0,10),breaks=knots,norder=6)
aaa = fdPar(basiss,Lfdobj = int2Lfd(2),lambda = 1e-5)
precfd = smooth.basis(1:10,t(F2[,-1]),aaa)
quartz()
plot(precfd$fd,xlab = 'Observations',ylab = 'Values',cex.lab=1.5,cex.axis=1.5,col="red")

# Do functional principal component analysis on the 35 temperature curves
# We choose 4 FPCs
temppca = pca.fd(precfd$fd,nharm=4)

#temppca$values are the eigenvalues
quartz()
par(mfrow=c(1,1),mar = c(8, 8, 4, 2))
plot(temppca$values[1:8],xlab='component',ylab='variance',col="blue",
     cex.lab=1.5,cex.axis=1.5,cex=2)

# plot the cumulative percentage explained total variations
quartz()
par(mfrow=c(1,1),mar = c(8, 8, 4, 2))
plot(cumsum(temppca$values[1:10])/sum(temppca$values),xlab='component',
     ylab='cumulative variance explained',col=4,cex.lab=3.5,
     cex.axis=3.5,cex=4)
abline(h=0.99)

# Show the mean curves - temppca$meanfd
quartz()
par(mfrow=c(1,1),mar = c(8, 8, 4, 2))
plot(precfd$fd,xlab='Observations',ylab='Values',cex.lab=1.5,cex.axis=1.5,col=4)
lines(temppca$meanfd,lwd=2.5,col=2)

# functional principal components
harmfd = temppca$harmonics
harmvals = eval.fd(1:10,harmfd)
dim(harmvals) # The top 4 FPCs

# plot the second FPC
quartz()
par(mfrow=c(1,1),mar = c(8, 8, 4, 2))
plot(1:10,harmvals[,2],xlab='Observations',ylab='PCs',
     lwd=4,lty=1,cex.lab=2,cex.axis=2,type='l')

# plot all 4 FPCs
quartz()
par(mfrow=c(1,1),mar = c(8, 8, 4, 2))
matplot(1:10,harmvals,xlab='Observations',ylab='PCs',
        lwd=4,lty=1,cex.lab=2.5,cex.axis=2.5,type='l')
legend(0,-0.07,c('PC1','PC2','PC3','PC4'),col=1:4,lty=1,lwd=5)
title('Principle Component Functions')

#Smoothing with smoothing parameter, lambda = 10^7
harmLfd = vec2Lfd(c(0,(2*pi/10)^2,0), c(0, 10))
tempfdPar2 = fdPar(basiss,harmLfd,1e-2)
tempfdPar3 = fdPar(basiss,harmLfd,1e7)
tempfd2 = smooth.basis(1:10,t(F2[,-1]),tempfdPar2)
ptemppca = pca.fd(precfd$fd,nharm=4,harmfdPar=tempfdPar3)
pharmfd = ptemppca$harmonics
pharmvals = eval.fd(1:10,pharmfd)
quartz()
par(mfrow=c(1,1),mar = c(8, 8, 4, 2))
matplot(1:10,pharmvals,xlab='Observations',ylab='PCs',
        lwd=2,lty=1,cex.lab=1.5,cex.axis=1.5,type='l')
legend(0,-0.07,c('PC1','PC2','PC3','PC4'),col=1:4,lty=1,lwd=2)
title('Smooth Principle Component Functions')

# plot the first FPC scores vs. the second FPC scores 
quartz()
par(mfrow=c(1,1),mar = c(8, 8, 4, 2))
plot(temppca$scores[,1:2],xlab='PC Score 1',ylab='PC Score 2',col=4,
     cex.lab=1.5,cex.axis=1.5,cex=1)
text(temppca$scores[,1],temppca$scores[,2],cex=1)

#Clustering (K means, 6 clusters)
install.packages("cluster")
library(cluster)
install.packages("shape")
library(shape)
set.seed(150)
kmean<-kmeans(temppca$scores[,c(2,1)],6,iter.max=300,nstart=10)
clusplot(temppca$scores[,c(2,1)],kmean$cluster,line=0,shade=F,color=TRUE,labels=4,plotchar=TRUE,span=TRUE,main=paste('Clusters'))