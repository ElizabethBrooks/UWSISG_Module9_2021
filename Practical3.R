##practical 3: Inbreeding and kinshhip
library(gaston)
library(hierfstat)

#writeLines('PATH="${RTOOLS40_HOME}\\usr\\bin;${PATH}"', con = "~/.Renviron")
#Sys.which("make")
#BiocManager::install("JGTeach", type = "source")
#library(JGTeach)


#Set working directory
setwd(" C:/Users/brookse/Desktop/UWSISG_UW_2021")


##Kinship and inbreeding from pedigrees
ped<-read.table("https://www2.unil.ch/popgen/teaching/SISGData/pedmono.txt",header=TRUE)
str(ped)
dim(ped)

# number of founders
# founders are individuals with missing sires and dams

sum(is.na(ped$sire) & is.na(ped$dam))

# families:  28:32; 59:62
ped[28:32,]
ped[59:62,]

#first inds of gen 1:5 :  c(21,47,77,117,159) 


ARM<-JGTeach::pedARM(ped$sire,ped$dam)

ARM[28:32,28:32]
ped[28:32,]
ARM[28:32,c(ped[28,2],ped[28,3])]
ARM[59:62,59:62]
ARM[59:62,c(ped[59,2],ped[59,3])]

#relatedness higher than 0.5 with their sibs
ARM[92:95,92:95]


#and their parents
ARM[92:95,c(ped[92,2],ped[92,3])]

#their parents are related
ARM[c(ped[92,2],ped[92,3]),c(ped[92,2],ped[92,3])]


hist(mat2vec(ARM))

ped.kin<-hierfstat::grm2kinship(ARM)
image(1:195,1:195,ped.kin)

#add vertical and horizontal lines between generations
genlim<-c(21,47,77,117,159)
abline(h=genlim-.5,v=genlim-.5)


##Kinship and inbreeding from markers, compared to pedigrees
genoped<-readRDS("geno_mono.RDS")
genoped.Kas<-hierfstat::beta.dosage(genoped,inb=FALSE,Mb=TRUE)
genoped.M<-with(genoped.Kas,betas*(1-MB)+MB)
genoped.Kc0<-JGTeach::Kc0(genoped.M,matching=TRUE)
bed.genoped<-gaston::as.bed.matrix(genoped)
genoped.std.GRM<-gaston::GRM(bed.genoped,autosome.only=FALSE)

par(mfrow=c(1,3))
plot(mat2vec(ARM)/2,mat2vec(genoped.Kas$betas),
     cex=0.5,col="red",main=expression(K[AS]),
     xlab="ped kinship",ylab="est. kinship")
abline(c(0,1))
plot(mat2vec(ARM)/2,mat2vec(genoped.Kc0),
     cex=0.5,col="red",main=expression(K[c0]),
     xlab="ped kinship",ylab="est. kinship")
abline(c(0,1))

# gaston::GRM reports relatedness rather than kinship, hence the /2

plot(mat2vec(ARM)/2,mat2vec(genoped.std.GRM)/2,
     cex=0.5,col="red",main="standard kinship",
     xlab="ped kinship",ylab="est. kinship")
abline(c(0,1))
par(mfrow=c(1,1))


#we divide by two because we want kinships rather than relatedness
mARM<-mean(mat2vec(ARM)/2)
ARMc<-(ARM/2-mARM)/(1-mARM)

par(mfrow=c(1,3))
plot(mat2vec(ARMc),mat2vec(genoped.Kas$betas),
     cex=0.5,col="red",main=expression(K[AS]),
     xlab="ped kinship",ylab="est. kinship")
abline(c(0,1))
plot(mat2vec(ARMc),mat2vec(genoped.Kc0),
     cex=0.5,col="red",main=expression(K[c0]),
     xlab="ped kinship",ylab="est. kinship")
abline(c(0,1))
plot(mat2vec(ARMc),mat2vec(genoped.std.GRM)/2,
     cex=0.5,col="red",main="standard kinship",
     xlab="ped kinship",ylab="est. kinship")
abline(c(0,1))
par(mfrow=c(1,1))


##kinship and inbreeding in the 1000 genomes
ch22<-read.VCF("chr22_Mb0_20.recode.vcf.gz")
ch22.M<-readRDS("matching.ch22.RDS")
Mb<-mean(mat2vec(ch22.M))

ch22.Kas<-(ch22.M-Mb)/(1-Mb)
ch22.Kc0<-JGTeach::Kc0(ch22.M,matching=TRUE)
ch22.std.GRM<-gaston::GRM(ch22)


#extract self kinship and convert it to inbreeding
ch22.Kas.inb<-diag(ch22.Kas)*2-1/2
diag(ch22.Kas)<-NA
ch22.Kc0.inb<-diag(ch22.Kc0)*2-1/2
diag(ch22.Kc0)<-NA
ch22.std.Inb<-diag(ch22.std.GRM)-1
diag(ch22.std.GRM)<-NA


samp.desc<-read.table("ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel",header=TRUE,stringsAsFactors = TRUE)

#to order samples
ospp<-with(samp.desc,order(super_pop,pop))
lsp<-cumsum(table(samp.desc$super_pop))

#png("1kg_ch22chunck_Kas.png")
image(1:2504,1:2504,ch22.Kas[ospp,ospp])
abline(h=lsp,v=lsp)
#dev.off()
#png("1kg_ch22chunck_Kc0.png")
image(1:2504,1:2504,ch22.Kc0[ospp,ospp])
abline(h=lsp,v=lsp)
#dev.off()
#png("1kg_ch22chunck_Kstdk.png")
image(1:2504,1:2504,ch22.std.GRM[ospp,ospp]/2)
abline(h=lsp,v=lsp)
#dev.off()