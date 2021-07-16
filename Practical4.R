##practical 4: Population structure

#Set working directory
setwd(" C:/Users/brookse/Desktop/UWSISG_UW_2021")

library(gaston)
#library(devtools)
#install_github("jgx65/JGTeach")
library(hierfstat)
#library("JGTeach")


##Drift
dum<-drift(nind=100,p0=0.5,nrep=1000,ngen=100,PlotIt = TRUE)


x<-drift(nind=100,ngen=1000)
gens<-c(1,2,3,5,10,20,50,100,1000)

par(mfrow=c(3,3))
for (i in gens)
  hist(x[i,],breaks=0:50/50,main=paste("generation:", i),xlab="Freq",ylab="")
par(mfrow=c(1,1))


x1<-drift(nind=50,ngen=1000,nrep=10000)
par(mfrow=c(3,3))
for (i in gens)
  hist(x1[i,],breaks=seq(0,1,0.01),main=paste("generation:", i),xlab="Freq",ylab="")
par(mfrow=c(1,1))


# distribution of fixation times
dist.tfix<-apply(x1,2,function(y) which(y>=1.0)[1])
# fixation probability
(pfix<-sum(!is.na(dist.tfix))/10000)
# mean time to fixation
(tfix<-mean(dist.tfix,na.rm=T)) 


hist(dist.tfix, main="Distribution of time to fixation, p0=0.50, N=50",xlab="Generations",ylab="")
abline(v=tfix, col="red",lwd=2)


# 50 loci with 20 alleles each, 
dat<-sim.genot.t(size=100,nbal=20,N=c(100,1000,10000),
                 nbloc=50,mig=0.0,t=50)

# the betas estimates population specific Fst from a
# fstat format data frame assuming random mating
betas(dat,nboot=100)$ci


##Individual populations Betas
# 50 loci with 20 alleles each, 
dat<-sim.genot.t(size=100,nbal=20,N=c(100,1000,10000),
                 nbloc=50,mig=0.0,t=200)

# the betas estimates population specific Fst from a
# fsat format data frame assuming random mating
betas(dat,nboot=100)$ci
#expectation
200/c(200,2000,20000)


# assumes mu=1e-8, r=1e-8, N0=1e+5, nbp=1e+5, 
# simulate 10 replicates, 100 diploid individuals from each pop
#  
# pop 4 will be the ancestral population with N4=N0=1e+5
# Only pops 1-3 are sampled, and there is no migration:
# -I 4 200 200 200 0 0 
# pop1 is a 1000th of anc pop, pop2 a 100th and pop 3 a 10th:
# -n 1 0.001 -n 2 0.01 -n 3 0.1
# split occured (backward) 50 generation ago
# in unit of 4* Anc pop size ->50/4/100000=0.000125:
# -ej 0.000125 1 4 -ej 0.000125 2 4 -ej 0.000125 3 4


bed<-ms2bed("3popsdrift.txt")
# expectation is (50/2/c(100,1000,10000)), 
# but works well for t<0.2N only
fst.dosage(bed,pop=rep(1:3,each=100)) 
50/c(200,2000,20000)


##2 populations system
# reproduces fig 1 in Weir and Goudet (2017)
#top row
first<-thet.bet.2pops(mu=0.0,m2=0.0,m1=0.0,n1=10000,n2=100,ngen=10000)
#middle row
second<-thet.bet.2pops(mu=0.001,m2=0.0,m1=0.0,n1=10000,n2=100,ngen=10000)
#bottom row
third<-thet.bet.2pops(mu=0.001,m2=0.0,m1=0.01,n1=10000,n2=100,ngen=10000)


#library(hierfstat)
gens<-1:10*100
ngen<-gens[length(gens)]
#generate data
mig.mat2<-matrix(c(.99,0.01,0,1),nrow=2,byrow=TRUE)
mig.mat2
mutrate<-10^{-3} #10^{-6} initial value

# As sim.genot.metapop.t only outputs genotypes from the 
# last generation, necessary to run the function
# with several end time points. This is the essence of 
# using lapply to a here

#  x<-lapply(gens,function(y){ 
#  sim.genot.metapop.t(t=y,nbal=20,nbpop=2,N=c(10000,100),
#                      mig=mig.mat2,nbloc=1000,mut=mutrate)})

# estimate betas for data sets
# beta.x<-lapply(x,betas,nboot=1000)


load("fig3.RData")
betas.ci<-lapply(beta.x,function(x) x$ci)
#expected values

etb1<-thet.bet.2pops(mu=mutrate,n1=10000,n2=100,m1=0.01,m2=0,ngen=ngen,plotit=FALSE)

gens<-1:10*100
ngen<-gens[length(gens)]

#create the plot
plot(1:ngen,etb1$Be[,1],type="l",ylim=range(etb1$Be,betas.ci),col="red",lwd=2,
     xlab="Generations",ylab=expression(beta))
lines(1:ngen,etb1$Be[,2],col="blue",lwd=2)
lines(1:ngen,etb1$Be[,3],lwd=2)
abline(h=0)
segments(gens,unlist(lapply(betas.ci,function(x) x[1,2])),gens,
         unlist(lapply(betas.ci,function(x) x[2,2])),lwd=2,col="blue")
segments(gens,unlist(lapply(betas.ci,function(x) x[1,1])),gens,
         unlist(lapply(betas.ci,function(x) x[2,1])),lwd=2,col="red")
points(gens,unlist(lapply(beta.x,function(x) x$betaiovl[1])),
       col="red",cex=1.5,pch=16)
points(gens,unlist(lapply(beta.x,function(x) x$betaiovl[2])),
       col="blue",cex=1.5,pch=16)
#beta_w added too
points(gens,unlist(lapply(beta.x,function(x) mean(x$betaiovl))),
       col="black",cex=1.5,pch=16)
title("2 populations with migration model. \n        N1=10'000; N2=100; m1=0.01; m2=0.0")


##Population specific Fsts from the 1000 genomes
ch22<-read.VCF("chr22_Mb0_20.recode.vcf.gz")
samp.desc.url<-"https://www2.unil.ch/popgen/teaching/SISG18/integrated_call_samples_v3.20130502.ALL.panel"
samp.desc<-read.table(samp.desc.url,header=TRUE,stringsAsFactors = TRUE)


with(samp.desc,table(pop,super_pop))


all.equal(ch22@ped$id,as.character(samp.desc$sample),check.attributes = FALSE)

# if not (in different order and or some missing) :
# use match(ch22@ped$id,samp.desc$sample) 


# From dosage data:
# fst.ch22.pop<-fst.dosage(ch22,pop=samp.desc$pop)
# fst.ch22.cont<-fst.dosage(ch22,pop=samp.desc$super_pop)

# More efficiently from matching data
# Matching.ch22<-matching(ch22)
# takes some time, so saved in matching.ch22.RDS available from # website

Matching.ch22<-readRDS("matching.ch22.RDS")
fs.ch22.pop<-fs.dosage(Matching.ch22,pop=samp.desc$pop,matching=TRUE)
fs.ch22.cont<-fs.dosage(Matching.ch22,pop=samp.desc$super_pop,matching=TRUE)
fs.ch22.pop$Fs[2,]
fs.ch22.cont$Fs[2,]


#extract popnames in the correct order
pnames<-unlist(lapply(strsplit(with(samp.desc,
                                    levels(factor(paste(super_pop,pop,sep=":")))),":"),function(x) x[[2]]))
#forces factor to take  this order                            
pop<-with(samp.desc,factor(pop,levels=pnames))
fs.ch22.pop<-fs.dosage(Matching.ch22,pop=pop,matching=TRUE)
coul<-rep(c("orange","gold","red","green","blue","purple"),c(2,5,4,5,5,5))
plot(fs.ch22.pop,las=2,cex.axis=0.5,col=coul)


##PCA
np<-5
dat5<-sim.genot(nbal=2,nbloc=1000,nbpop=np,mig=0.003)


dos<-biall2dos(dat5[,-1])


ps<-colMeans(dos)/2
sdp<-(ps*(1-ps))^(0.5)
#first, filter out the fixed loci
nfixed<-which(ps>0 & ps <1)
dosnf<-dos[,nfixed]

X<-scale(dosnf)

# or
# X<-sweep(dosnf,2,ps[nfixed],FUN="-")
# X<-sweep(X,2,sdp[nfixed],FUN="/")


XXT<-tcrossprod(X)


eigx<-eigen(XXT)


ind.coord<-sweep(eigx$vectors,2,eigx$values^0.5,FUN="*")


prx<-prcomp(X)
#check equality of the absolute values of the  PCs

all.equal(abs(matrix(prx$x)),abs(matrix(ind.coord))) 


par(mfrow=c(2,2))
plot((eigx$values/sum(eigx$values))[1:20],type="h",
     main="screeplot",xlab="Eigen values",ylab="")
plot(prx$x[,1:2],col=rep(1:np,each=50),pch=16)
plot(prx$x[,3:4],col=rep(1:np,each=50),pch=16)
plot(prx$x[,5:6],col=rep(1:np,each=50),pch=16)

par(mfrow=c(1,1))


colpca<-rep(1:np,each=50)
ns<-5
plot(1:ns,prx$x[1,1:ns],col=colpca[1],type="l",ylim=range(prx$x[,1:ns]),
     xlab="Axis",ylab="coord.",main="par. coord. plot")
for (i in 2:nrow(prx$x)) lines(1:ns,prx$x[i,1:ns],col=colpca[i])


X<-scale(dosnf,scale=FALSE)
XXT<-tcrossprod(X)
eigx<-eigen(XXT)
ind.coord<-sweep(eigx$vectors,2,eigx$values^0.5,FUN="*")

prx<-prcomp(X)

all.equal(abs(matrix(prx$x)),abs(matrix(ind.coord))) 


par(mfrow=c(2,2))
plot((eigx$values/sum(eigx$values))[1:20],type="h",
     main="screeplot",xlab="Eigen values",ylab="")
plot(prx$x[,1:2],col=rep(1:np,each=50),pch=16)
plot(prx$x[,3:4],col=rep(1:np,each=50),pch=16)
plot(prx$x[,5:6],col=rep(1:np,each=50),pch=16)

par(mfrow=c(1,1))


Kas<-beta.dosage(dos)
eigKas<-eigen(Kas)

ind.coord<-sweep(eigKas$vectors,2,eigKas$values^0.5,FUN="*")

par(mfrow=c(2,2))
plot((eigKas$values/sum(eigKas$values))[1:20],type="h",
     main="screeplot",xlab="Eigen values",ylab="")
plot(ind.coord[,1:2],col=rep(1:np,each=50),pch=16)
plot(ind.coord[,3:4],col=rep(1:np,each=50),pch=16)
plot(ind.coord[,5:6],col=rep(1:np,each=50),pch=16)

par(mfrow=c(1,1))


np<-3
dat3<-sim.genot(nbal=2,nbloc=1000,nbpop=np,mig=0.003)
dos<-biall2dos(dat3[,-1])
ps<-colMeans(dos)/2
sdp<-(ps*(1-ps))^(0.5)
#first, filter out the fixed loci
nfixed<-which(ps>0 & ps <1)
dosnf<-dos[,nfixed]

X<-scale(dosnf)
prx<-prcomp(X)
par(mfrow=c(2,2))
plot((prx$sdev^2/sum(prx$sdev^2))[1:20],type="h",
     main="screeplot",xlab="Eigen values",ylab="")
plot(prx$x[,1:2],col=rep(1:np,each=50),pch=16)
plot(prx$x[,3:4],col=rep(1:np,each=50),pch=16)
plot(prx$x[,5:6],col=rep(1:np,each=50),pch=16)

par(mfrow=c(1,1))


np<-5
#migration 3 times larger

dat5<-sim.genot(nbal=2,nbloc=1000,nbpop=np,mig=0.01)
dos<-biall2dos(dat5[,-1])
ps<-colMeans(dos)/2
sdp<-(ps*(1-ps))^(0.5)
#first, filter out the fixed loci
nfixed<-which(ps>0 & ps <1)
dosnf<-dos[,nfixed]

X<-scale(dosnf)
prx<-prcomp(X)
par(mfrow=c(2,2))
plot((prx$sdev^2/sum(prx$sdev^2))[1:20],type="h",
     main="screeplot",xlab="Eigen values",ylab="")
plot(prx$x[,1:2],col=rep(1:np,each=50),pch=16)
plot(prx$x[,3:4],col=rep(1:np,each=50),pch=16)
plot(prx$x[,5:6],col=rep(1:np,each=50),pch=16)

par(mfrow=c(1,1))

colpca<-rep(1:np,each=50)
ns<-5
plot(1:ns,prx$x[1,1:ns],col=colpca[1],type="l",
     ylim=range(prx$x[,1:ns]),xlab="Axis",
     ylab="coord.",main="par. coord. plot")
for (i in 2:nrow(prx$x)) lines(1:ns,prx$x[i,1:ns],col=colpca[i])


np<-8
nl<-1000
t<-1000
M<-matrix(0,ncol=np,nrow=np)
diag(M[-1,-np])<-0.05
diag(M[-np,-1])<-0.05
diag(M)<-1-rowSums(M,na.rm=TRUE)

ss1d<-sim.genot.metapop.t(nbal=2,nbpop=np,nbloc=nl,mig=M,t=t)
dos<-biall2dos(ss1d[,-1])

ps<-colMeans(dos)/2
sdp<-(ps*(1-ps))^(0.5)
#first, filter out the fixed loci
nfixed<-which(ps>0 & ps <1)
dosnf<-dos[,nfixed]

X<-scale(dosnf)
prx<-prcomp(X)
par(mfrow=c(2,2))
plot((prx$sdev^2/sum(prx$sdev^2))[1:20],type="h",
     main="screeplot",xlab="Eigen values",ylab="")
plot(prx$x[,1:2],col=rep(1:np,each=50),pch=16)
plot(prx$x[,3:4],col=rep(1:np,each=50),pch=16)
plot(prx$x[,5:6],col=rep(1:np,each=50),pch=16)

par(mfrow=c(1,1))