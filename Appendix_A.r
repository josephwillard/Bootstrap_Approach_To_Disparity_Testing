#Bootstrap Approach to Disparity Testing with Source Uncertainty in Data
#(I)####################
cat("GENERATE 2^n-2 CONFIGURATIONS FOR DATA SET OF LENGTH n\n")
cat("INPUT n ON NEXT LINE\n")
n<-12
mm<-rep(list(c(-1,1)),n)
con<-expand.grid(mm)
con[1,]
con[2^n,]
config<-con[-c(1,2^n),]
config<-as.matrix(config)
dim(config)
config[1,]
config[(2^n)-2,]
config.tr<-t(config)
dim(config.tr)
#(II)####################
#INSERT SIMULATED DATA AND po FOR 'total' RECORDS
#use set.seed(xxx) if sample draw is to be repeated
#set.seed(300)
#reca records for pop A ~ N(muA,sd=sig); recb for pop B ~ N(muB,sd=sig)
#po are the probabilities that the obs came from A
#po draws for A ~ beta(aA,bA); for B ~ beta(aB,bB)
#If X ~ beta(a,b), E(X)=a/(a+b), var(X)=ab/[(a+b+1)(a+b)^2]
reca=1000;recb=1000;mua=4;mub=4;sig=1
aA<-58.1;bA<-24.9;aB<-24.9;bB<-58.1
total<-reca+recb
xa<-rnorm(reca,mean=mua,sd=sig)
xb<-rnorm(recb,mean=mub,sd=sig)
poa<-rbeta(reca,aA,bA)
pob<-rbeta(recb,aB,bB)
data<-c(xa,xb)
po<-c(poa,pob)
#print first five values of data and po
round(cbind(data[1:5],po[1:5]),4)
#check on the means and std devs of xa,xb,poa and pob
checkdataa<-c(mean(data[1:reca]),sd(data[1:reca]))
inc<-reca+1
checkdatab<-c(mean(data[inc:total]),sd(data[inc:total]))
checkpoa<-c(mean(po[1:reca]),sd(po[1:reca]))
checkpob<-c(mean(po[inc:total]),sd(po[inc:total]))
round(cbind(checkdataa,checkdatab,checkpoa,checkpob),4)
#(III)####################
message("DRAW 'd' RANDOM SAMPLES OF SIZE ",n," FROM (data,po)")
#match the data and po to the random sample numbers
message("data1 and po1 are values for the random sample of size ", n)
pmean1<-NULL; pstd.dev1<-NULL;per<-NULL
perct_05<-NULL;perct_10<-NULL;perct_25<-NULL;perct_50<-NULL
perct_75<-NULL;perct_90<-NULL;perct_95<-NULL
#set d for the number of random samples (data1,po1)
cat("INPUT d ON NEXT LINE\n")
d<-1000
perct<-NULL
for (m in 1:d) {
data1<-NULL;po1<-NULL;sam<-NULL
sam<-sample(1:total,n,replace=TRUE)
for (k in 1:n) {
	data1[k]<-data[sam[k]]
	po1[k]<-po[sam[k]]
}
#(IV)####################
#COMPUTE THE VECTOR OF PROBABILITIES FOR THE (2^n)-2 CONFIGURATIONS
#AND FOR THE ADJUSTED PROBABILITIES
#po1 are the probabilities that the obs came from A
#probs is the vector of probabilities for the (2^n)-2 configurations
p<-NULL
lim<-(2^n)-2
probs<-c(rep(1),lim)
for (j in 1:lim) {
  for (i in 1:n){
    if (config.tr[i,j]==-1) p[i]<-1-po1[i]
    if (config.tr[i,j]==1) p[i]<-po1[i]}
  probs[j]<-prod(p)
}
#adjprobs is the vector of probabilties for the (2^n)-2 configurations
#normalized to sum to one
adjprobs<-NULL
for (j in 1:lim) {adjprobs[j]<-probs[j]/sum(probs)}
#(V)####################
#COMPUTE THE SUM OF data1 AND SAMPLE SIZES FOR EACH OF THE CONFIGS
#sumA (sumB) is sum of values assigned to category A (B)
#sumN (sumM) is the number of values assigned to cat A (B)
#for each of the (2^n)-2 allocations sumN+sumM=n and
#sumA+sumB=total, the sum of the observations
sumA<-c(rep(0,lim)); sumN<-c(rep(0,lim));
sumB<-c(rep(0,lim)); sumM<-c(rep(0,lim))
for (i in 1:lim){
for (j in 1:n){
if (config.tr[j,i]==1){
  sumA[i]=sumA[i]+data1[j]
  sumN[i]=sumN[i]+config.tr[j,i]}
if (config.tr[j,i]==-1){
  sumB[i]=sumB[i]+data1[j]
  sumM[i]=sumM[i]-config.tr[j,i]}
	}
}
#(VI)####################
#COMPUTE TEST STATISTIC AND p-VALUES
#test Ho: meanA = meanB vs Ha: meanA > meanB
#assume normal populations with common known std deviation
#input common known std. dev. value, sigma, based on assumed data model
sigma<-1
z<-NULL; pval<-NULL; diff<-NULL; sddiff<-NULL
for (i in 1:lim){
	diff[i] <- (sumA[i]/sumN[i])-(sumB[i]/sumM[i])
	sddiff[i] <- sigma*sqrt((1/sumN[i])+(1/sumM[i]))
	z[i] <- diff[i]/sddiff[i]
	pval[i] <- 1-pnorm(z[i],mean=0,sd=1)
	}
pval.ordered<-NULL
pval.ordered<-sort(pval)
diff.pval.ord<-NULL
for (j in 2:lim) {diff.pval.ord[j-1]<-pval.ordered[j]-pval.ordered[j-1]}
#
#construct data frame for configs, pvals, adjprobs, cdf
#sort data frame by pval
#calculate cdf of pvals
run<-c(seq(from=1,to=lim))
test.df<-data.frame(run,pval,adjprobs)
test.df<-test.df[order(pval),]
cdf.pval<-NULL
cdf.pval<-cumsum(test.df$adjprobs)
test.df["cdf.pval"]<-cdf.pval
#(VII)####################
#CALCULATE THE MEAN, STANDARD DEV, AND PERCENTILES OF pval
pmean<-c(rep(0),m); pvar<-c(rep(0),m); pstd.dev<-c(rep(0),m)
pmean<-sum(test.df$pval*test.df$adjprobs)
pvar<-sum(((test.df$pval-pmean)^2)*test.df$adjprobs)
pstd.dev<-sqrt(pvar)
#calculate the percentiles of pval
#using linear interpolation with two bordering pval for a given percentile
per<-c(5,10,25,50,75,90,95)
lengthp<-length(per)
lim1<-lim-1
for (j in 1:lengthp) {
	for (i in 1:lim1) {
		if ((cdf.pval[i]<=per[j]/100)&(cdf.pval[i+1]>per[j]/100))
		{x1<-cdf.pval[i]
		x2<-cdf.pval[i+1]
		y1<-pval.ordered[i]
		y2<-pval.ordered[i+1]
		x<-per[j]/100
		perct[j]<-y1+((x-x1)*(y2-y1))/(x2-x1)
		}
	}
}
pmean1[m]<-pmean
pstd.dev1[m]<-pstd.dev
perct_05[m]<-perct[1]
perct_10[m]<-perct[2]
perct_25[m]<-perct[3]
perct_50[m]<-perct[4]
perct_75[m]<-perct[5]
perct_90[m]<-perct[6]
perct_95[m]<-perct[7]
#end of sampling from (data,po)
}
#(VIII)####################
#DISPLAY RESULTS
message("These results are for two normal populations with means ",mua ," and " ,mub)
message("and common standard deviation equal to ",sig)
message("The parameters for the beta distribution for pop A are ",aA ," and ",bA)
message("and for population B ",aB," and ",bB)
message("The sample sizes for the two populations are ",reca," and ",recb)
message("The number of bootstrap samples is ",d, " and each is of size ",n)
cbind(sumA[1:5],sumN[1:5],sumB[1:5],sumM[1:5])
#print the mean & std.dev of p-values for the first 5 samples of size n
round(rbind(pmean1[1:5],pstd.dev1[1:5]),4)
message("percentiles(5,10,25,50,75,90,95)")
#print percentiles of first five samples of size n
message("If NA appears, it means for one or more of the ",d," bootstrap samples
the smallest positive cdf.pval, cdf.pval[1], is larger than per/100")
round(cbind(perct_05[1:5],perct_10[1:5],perct_25[1:5],
perct_50[1:5],perct_75[1:5],perct_90[1:5],perct_95[1:5]),4)
results<-data.frame(perct_05,perct_10,perct_25,perct_50,perct_75,perct_90,perct_95)
round(apply(results,2,mean),4)
round(apply(results,2,sd),4)
#calculate std dev of percentile sample means
sd.meanp1<-sd(perct_05)/sqrt(d)
sd.meanp2<-sd(perct_10)/sqrt(d)
sd.meanp3<-sd(perct_25)/sqrt(d)
sd.meanp4<-sd(perct_50)/sqrt(d)
sd.meanp5<-sd(perct_75)/sqrt(d)
sd.meanp6<-sd(perct_90)/sqrt(d)
sd.meanp7<-sd(perct_95)/sqrt(d)
sd.meanp<-c(sd.meanp1,sd.meanp2,sd.meanp3,sd.meanp4,sd.meanp5,sd.meanp6,sd.meanp7)
round(sd.meanp,4)
#(IX)####################
cat("A SUMMARY OF RESULTS \n")
cat("ESTIMATES OF POPULATION p-VALUE PERCENTILES ARE MEANS OF BOOTSTRAP p-VALUE PERCENTILES \n")
cat("STANDARD ERRORS OF THESE ESTIMATES ARE THE STANDARD DEVIATIONS OF THE BOOTSTRAP \n")
cat("PERCENTILES DIVIDED BY THE SQUARE ROOT OF THE NUMBER OF THE BOOTSTRAP SAMPLES, d \n")
results1<-data.frame(perct_05,perct_10,perct_25,perct_50,perct_75,
perct_90,perct_95)
percentile.estimates<-apply(results1,2,mean)
std.errors<-apply(results1,2,sd)/sqrt(d)
analysis<-data.frame(percentile.estimates,std.errors)
message("A SUMMARY OF THE ANALYSIS RESULTS FOLLOWS...")
round(analysis,4)
message("The averages of the ",d," p-value means and std. devs. are ",
	round(mean(pmean1),4)," and ", round(mean(pstd.dev1),4))
message("The standard deviations of the ",d," p-value means and std.devs. are ",
	round(sd(pmean1),4)," and ", round(sd(pstd.dev1),4))
message("The standard errors of the avgs. of the ",d," p-value means and std.devs. are ",
	round(sd(pmean1)/sqrt(d),4)," and ", round(sd(pstd.dev1)/sqrt(d),4))
#CHOOSE THE LEVEL OF CONFIDENCE FOR THE HYPOTHESIS TEST
#90, 95, or 99 percent
levcon<-99
if (levcon==90){zalpha<-1.28}
if (levcon==95){zalpha<-1.645}
if (levcon==99){zalpha<-2.33}
hyp<-c(levcon,zalpha)
hyp
pmean2<-mean(pmean1)
semean2<-sd(pmean1)/sqrt(d)
UCL<-pmean2+zalpha*semean2
UCL<-round(UCL,4)
if (UCL<=0.5) {message("The UCL is ",UCL ," <= 0.5, so
	reject the null hypothesis in favor of disparity
	at confidence level ",levcon ,"%")
	}
if (UCL>0.5) {message("The UCL is ",UCL ," > 0.5, so
	do not reject the null hypothesis of no disparity
	at confidence level ",levcon ,"%")
	}
#(X)####################
#message("CALCULATE p-VALUE FOR RECORDS 'DATA' ASSUMING NO UNCERTAINTY IN DATA SOURCE")
p.val<-NULL
z.test<-function(asum,bsum,asize,bsize,sig){
	differ<-(asum/asize)-(bsum/bsize)
	diffsd<-sig*sqrt((1/asize)+(1/bsize))
	z.test<-differ/diffsd
	p.val<-1-pnorm(z.test,mean=0,sd=1)
	return(p.val)
}
asum<-sum(data[1:reca])
bsum<-sum(data[reca+1:recb])
asize<-reca
bsize<-recb
u<-z.test(asum,bsum,asize,bsize,sig)
message("The p-value for the data assuming no uncertainty in the data source is ", u)
#(XI)####################
message("ASSESSMENT FOR THE LAST BOOTSTRAP SAMPLE, ", d)
message("FOR EACH BOOTSTRAP SAMPLE, FOR EVERY z VALUE THERE IS A CORRESPONDING -z VALUE
	THUS, FOR THOSE TWO z VALUES THE SUM OF THEIR p-VALUES IS EQUAL TO 1
	SO THE SUM OF THE p-VALUES FOR THE (2^n-2) CONFIGURATIONS IS (2^n-2)/2 = 2^(n-1)-1
	AND THE AVERAGE OF THE p-VALUES IS 0.5")
sum.pval<-sum(pval)
avg.pval<-sum.pval/(2^n-2)
message("CHECK: SUM OF p-VALUES FOR THIS CONFIGURATION IS ",sum.pval,
	" AND THE AVERAGE IS ",avg.pval)
new<-data.frame(z,pval,adjprobs)
dim(new)
revise<-new[order(pval),]
message("The plot is for mua = ",mua,"; mub = ",mub,"; n = ",n,"; d = ",d,
	"; reca = ",reca,"; recb = ",recb)
color1<-c("red","blue","purple","black","green","orange","brown")
color2<-c("brown","orange","green","black","purple","blue","red")
tag<-c("95%","90%","75%","50%","25%","10%","5%")
par(mfrow=c(1,2))
plot(pval,adjprobs,main="VERTICAL LINES CUT SEVEN p-VALUE PERCENTILES
	FOR THE LAST BOOTSTRAP SAMPLE",xlab="p-VALUE",ylab="PROBABILITY")
abline(v=perct,col=color1)
legend(x="topright",box.col="brown",bg="yellow",
	legend=tag,fill=color2)
plot(test.df$pval,test.df$cdf.pval,main="CDF OF p-VALUE FOR LAST BOOTSTRAP SAMPLE
	horizontal lines for cdf values; vertical lines for percentiles",ylab="CDF of p-VALUES",
	xlab="p-VALUE")
abline(h=per/100,col=color1)
abline(v=perct,col=color1)
legend(x="topleft",box.col="brown",bg="yellow",
	legend=tag,fill=color2)
perct
pmean
pstd.dev
cdf.pval[1:10]
