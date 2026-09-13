library(MASS)
library(boot)

estim_betagb<-read.csv('beta.csv',h=T)
varcovgb<-read.csv('varcov.csv',h=F)

numbetagb<-c(1:3)
lbeta<-length(numbetagb)
Mugb<-estim_betagb[numbetagb,2]
Sigmagb<-varcovgb[numbetagb,numbetagb]
n<-10000
rp<-mvrnorm(n, Mugb, Sigmagb)


pop<-matrix(NA,nrow=n,ncol=1)
fem<-matrix(NA,nrow=n,ncol=1)
mal<-matrix(NA,nrow=n,ncol=1)
sr<-matrix(NA,nrow=n,ncol=1)

for(i in 1:n)
{
  pop[i,1] <- 127 / (1 - (1 - inv.logit(rp[i,3]))^3)
  fem[i,1] <- 38 / (1 - (1 - inv.logit(rp[i,3]))^3)
  mal[i,1] <- 89 / (1 - (1 - inv.logit(rp[i,3]))^3)
  sr[i,1] <- mal[i,1]/fem[i,1]
}

popSize <- quantile(pop, probs=c(0.025,0.5,0.975))
femSize <- quantile(fem, probs=c(0.025,0.5,0.975))
malSize <- quantile(mal, probs=c(0.025,0.5,0.975))
SRat <- quantile(sr, probs=c(0.025,0.5,0.975))

quartiles_pop <- data.frame(
    Variable = (c("population size","females", "males", "sex ratio")),
     Q_0.025 = as.numeric(c(popSize[1],femSize[1],malSize[1],SRat[1])),
     Q_0.5   = as.numeric(c(popSize[2],femSize[2],malSize[2],SRat[2])),
     Q_0.975 = as.numeric(c(popSize[3],femSize[3],malSize[3],SRat[3])))

write.csv(quartiles_pop, file="results.csv", row.names = FALSE)