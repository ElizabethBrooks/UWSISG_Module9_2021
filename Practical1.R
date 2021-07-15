#install packages
install.packages("hierfstat")
install.packages("gaston")
install.packages("HardyWeinberg")
install.packages("adegenet")
install.packages("pegas")
install.packages("ape")
#install.packages("devtools")
#library(devtools)
#install_github("jgx65/JGTeach")


#Set working directory
setwd("C:/Users/brookse/Downloads/")


##Importing data into R
library(hierfstat)
ch22<-read.VCF("chr22_Mb0_20.recode.vcf.gz")


#ch22 matrix with nrows inds and ncol snps
dim(ch22)
str(ch22)


#library(gaston)
#chr22<-read.bed.matrix("chr22.1kg")


#library(SNPRelate)
#snpgdsVCF2GDS("chr22.vcf.gz","chr22.gds")
#f<-snpgdsOpen("chr22.gds")


samp.desc.fname<-"integrated_call_samples_v3.20130502.ALL.panel"
ftp.path<-"ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/release/20130502/"
samp.desc.url<-paste0(ftp.path,samp.desc.fname)
onekg.desc<-read.table(samp.desc.url,header=TRUE,stringsAsFactors = TRUE)
#checks that the order of samples in bed file and description file are the same
all.equal(as.character(ch22@ped$id),as.character(onekg.desc$sample))
boxplot(ch22@ped$hz~onekg.desc$super_pop)
#and per population
boxplot(ch22@ped$hz~onekg.desc$pop,las=2)
#same, sorted by continent
boxplot(ch22@ped$hz~with(onekg.desc,factor(super_pop:pop)),las=2)


##Simulating Genetic/Genomic data
library(hierfstat)
#precede the function name with ? to get help, e.g. ?sim.genot

dat<-sim.genot(nbpop=4,nbloc=100,nbal=2,size=50,N=1000,mig=0.001)


#dim to get the number of rows and columns
dim(dat)

#column names
head(names(dat))

#structure of the data set

str(dat)


dos<-biall2dos(dat[,-1])
str(dos)


##Allelic frequencies
#library(gaston)
#library(hierfstat)
#library(JGTeach)


pan<-ms2bed("pan.txt")
hist(pan@p,breaks=101)


ch22<-read.VCF("chr22_Mb0_20.recode.vcf.gz")
samp.desc.file<-"https://www2.unil.ch/popgen/teaching/SISG18/integrated_call_samples_v3.20130502.ALL.panel"
samp.desc<-read.table(samp.desc.file,header=TRUE)
names(samp.desc)
str(samp.desc)

AFR<-which(samp.desc$super_pop=="AFR")
EAS<-which(samp.desc$super_pop=="EAS")
par(mfrow=c(2,2))
#AFR hist
hist(ch22[AFR,]@p,breaks=101,main="AFR",xlab="Allele Freq.")
#EAS hist
hist(ch22[EAS,]@p,breaks=101,main="EAS",xlab="Allele Freq.")
#PAN ms hist
hist(pan@p,breaks=101,main="panmictic pop with ms",xlab="Allele Freq.")
#simulate data from panmictic pop with sim.genot
dat<-sim.genot(nbal=2,nbpop=1,size=100,nbloc=10000)
#convert to dosage
dos<-biall2dos(dat[,-1])
#colMeans will give twice the frequencies
#->colMean(dos)/2 gives the frequencies
hist(colMeans(dos)/2,breaks=101,main="panmictic pop with sim.genot",xlab="Allele Freq.")
par(mfrow=c(1,1))


tab<-matrix(c(100,200,100,150,100,150,200,
              0,200,4,72,324,22,36,342,40,0,360),ncol=3,byrow=TRUE)
pr<-(tab[,1]+tab[,2]/2)/400 #freq of reference
f<-1-tab[,2]/400/(2*pr*(1-pr)) # inbreeding coefficient f
v<-pr*(1-pr)*(1+f)/2/400 # variance from genotype counts
round(cbind(tab,pr,f,v),digits=5)


#read pan.txt into a bed object
pan<-ms2bed("pan.txt")
#plot snps hz against p
plot(pan@p,pan@snps$hz,col="red",pch=16,cex=0.3)

#add expected prop of heterozygotes given the frequency
x<-0:1000/1000
lines(x,2*x*(1-x),col="blue")


pan50<-pan[1:50,]
plot(pan50@p,pan50@snps$hz,col="red",pch=16,cex=0.3)
lines(x,2*x*(1-x),col="blue")