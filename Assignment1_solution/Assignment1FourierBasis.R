library('fda')

y = melanoma

for(i in c(3,7,15,23,37)){
  
  yearbasiss37 <- create.fourier.basis(c(0, 37),i)
  bvalss = eval.basis(1:37,yearbasiss37)
  
  # Evaluate the first derivative of the basis functions at 1,2,...,37
  d1bvalss = eval.basis(1:37,yearbasiss37,Lfdobj=1)
  
  # Evaluate the second derivative of the basis functions at 1,2,...,37
  d2bvalss = eval.basis(1:37,yearbasiss37,Lfdobj=2)
  
  Q = bvalss
  
  yhatt = Q%*%solve(t(Q)%*%Q)%*%(t(Q)%*%as.matrix(y))
  
  # Estimate the first derivative of f(t)
  d1yhatt = d1bvalss%*%solve(t(Q)%*%Q)%*%(t(Q)%*%as.matrix(y))
  
  # Estimate the second derivative of f(t)
  d2yhatt = d2bvalss%*%solve(t(Q)%*%Q)%*%(t(Q)%*%as.matrix(y))
  
  quartz()
  par(mfrow=c(1,1),mar = c(8, 8, 4, 2))
  title = paste(i,' Basis Functions',sep='')
  plot(y,col=2,cex=1.5,xlab='year',ylab='incidence',cex.lab=1.5,
       main=title(),cex.axis=1.5)
  lines(yhatt,lwd=2,col=4)
  quartz()
  par(mfrow=c(1,1),mar = c(8, 8, 4, 2))
  plot(d1yhatt,col=4,type='l',lwd=2,cex.lab=1.5,cex.axis=1.5,
       ylab = 'D1 incidence',xlab='year',main=title())
  lines(yhatt,lwd=2,col=4)
  quartz()
  par(mfrow=c(1,1),mar = c(8, 8, 4, 2))
  plot(d2yhatt,col=4,type='l',lwd=2,cex.lab=1.5,cex.axis=1.5,
       ylab = 'D2 incidence',xlab='year',main=title())
  readline(prompt="Press [enter] to continue")
}

