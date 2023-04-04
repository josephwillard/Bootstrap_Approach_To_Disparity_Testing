WD<-c(0.0292,0.0482,0.1116,0.2610,0.5141,0.7494,0.8510)
WOD<-c(0.0888,0.1434,0.2886,0.5084,0.7247,0.8650,0.9163)
per<-c(0.05,0.10,0.25,0.50,0.75,0.90,0.95)
rnames<-c("WD","WOD")
cnames<-c("5%","10%","25%","50%","75%","90%","95%")
matrix.boot<-matrix(c(WD,WOD),byrow=TRUE,nrow=2,ncol=7)
dimnames=list(rnames,cnames)
barplot(matrix.boot,beside=TRUE,main="Percentiles with and without disparity",
xlab="Percentages",ylab="Percentiles",col=c("red","green"),
names.arg=cnames,cex.main="1.0")
legend("topleft",c("with disparity","without disparity"),
fill=c("red","green"),bg="yellow")
