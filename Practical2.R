#Set working directory
setwd("C:/Users/brookse/Desktop/UWSISG_UW_2021")


#install.packages("devtools")
#library(devtools)
#install_github("jgx65/JGTeach")

library(gaston)
library(hierfstat)
#library(JGTeach)


# F=1-Ho/He
pan<-ms2bed("pan.txt")
inb.coeff<-1-pan@snps$hz/2/pan@p/(1-pan@p)
#nb inds
ni<-dim(pan)[1]
x2<-ni*inb.coeff^2
p.val.x2<-pchisq(x2,df=1,lower=FALSE)
nl<-dim(pan)[2]
#theo dist p.val under null
theo.pval<-1:nl/nl
plot(-log10(theo.pval),-log10(sort(p.val.x2)),
     col="red",cex=0.5,xlab="Theo p val dist",
     ylab="emp p-val dist x2");abline(c(0,1))


