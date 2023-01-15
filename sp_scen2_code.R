## Surrogate Paradox Measures ##
## Code for Scenario 2        ##

# Packages
library(mvtnorm)
library(mnormt)
library(bayesSurv)
library(MASS)
library(coda)

# Fill in the values of the following variables
# n = total number of subjects across all trials 
# ntrial = number of trials 

# Fill in the following variables with your data
data=matrix(0,n,5)
data[,1:2]=ST     # ST is an nx2 matrix with the surrogate and true endpoint for each patient
data[,3]=Z        # Z is a length n vector with the treatment assignment
data[,4]=X.01     # X.01 is a length n vector with the covariate value
data[,5]=trialID  # trialID is a length n vector with the trial ID 


n=length(data[,1])
Y=array(0,c(2,n))
Y[1,]=data[,1]
Y[2,]=data[,2]

# Fixed effect matrix 
X=array(rep(0,2*8*n),dim=c(2,8,n))
for(i in 1:length(data[,1])){
  X[1,1,i]=1
  X[1,2,i]=0
  X[1,3,i]=Z[i]
  X[1,4,i]=0
  X[1,5,i]=X.01[i]
  X[1,6,i]=0
  X[1,7,i]=X.01[i]*Z[i]
  X[1,8,i]=0
  
  X[2,1,i]=0
  X[2,2,i]=1
  X[2,3,i]=0
  X[2,4,i]=Z[i]
  X[2,5,i]=0
  X[2,6,i]=X.01[i]
  X[2,7,i]=0
  X[2,8,i]=X.01[i]*Z[i]
}

# Random effect matrix
W=array(rep(0,2*8*n),dim=c(2,8,n))
for(i in 1:length(data[,1])){
  W[1,1,i]=1
  W[1,2,i]=0
  W[1,3,i]=Z[i]
  W[1,4,i]=0
  W[1,5,i]=X.01[i]
  W[1,6,i]=0
  W[1,7,i]=X.01[i]*Z[i]
  W[1,8,i]=0
  
  W[2,1,i]=1
  W[2,2,i]=0
  W[2,3,i]=Z[i]
  W[2,4,i]=0
  W[2,5,i]=X.01[i]
  W[2,6,i]=0
  W[2,7,i]=X.01[i]*Z[i]
  W[2,8,i]=0
}

Ds<-matrix(rep(0,8*8),8,8)
diag(Ds) <- 1

Dr<-matrix(rep(0,8*8),8,8)
Dr[!(diag(Dr))] <- 0.3
diag(Dr) <- 1

D=Ds%*%Dr%*%Ds

Ss <- matrix(rep(0,2*2),2,2)
Ss[1,1] <- 1
Ss[2,2] <- 1
Ss[1,2] <- 0
Ss[2,1] <- 0

Sr <- matrix(rep(0,2*2),2,2)
Sr[1,1] <- 1
Sr[2,2] <- 1
Sr[1,2] <- 0.5
Sr[2,1] <- 0.5

S=Ss%*%Sr%*%Ss

# Starting values
beta=c(1, 1, 2, 1, 0, 0, -1, 1)
eta=rmnorm(ntrial, c(0,0,0,0,0,0,0,0), Ds%*%Dr%*%Ds)

# Store values 
holdbet<-matrix(rep(0,8*SIM),SIM,8)
holdeta<-array(rep(0,ntrial*8*SIM),dim=c(ntrial,8,SIM))
holdD<-array(rep(0,8*8*SIM),dim=c(8,8,SIM))
holdsig<-array(rep(0,2*2*SIM),dim=c(2,2,SIM))

vd <- 9 # 5 for scenario 1, 9 for scenario 2 
vsigma <- 3

F<-matrix(rep(0,8*8),8,8)
diag(F) <- 1/10 # 1/(q+2) where q = length of random effects (4 in scenario 1, 8 in scenario 2)

E<-matrix(rep(0,2*2),2,2)
diag(E) <- 1/4 # 1/(q+2) where q = length of residual variance effects (2)

sig<-matrix(rep(0,4*4),4,4)
diag(sig) <- 100

SIM<-10000
sim=1
while(sim<=SIM){
  
  # Dinv conditional
  x=0
  for(k in 1:ntrial){
    x1=t(t(eta[k,]))%*%t(eta[k,])
    x=x+x1
  }
  scale=solve((x+solve(F)))
  
  Dinv=stats::rWishart(1, ntrial+vd, scale)[,,1]
  D=solve(Dinv)
  
  # Sigmainv conditional 
  x=0
  mu=array(0,c(2,1))
  et=array(0,c(2,1))
  for(i in 1:n){
    mu[1,]=beta[1]+beta[3]*Z[i]+beta[5]*X.01[i]+beta[7]*Z[i]*X.01[i]
    mu[2,]=beta[2]+beta[4]*Z[i]+beta[6]*X.01[i]+beta[8]*Z[i]*X.01[i]
    for(k in 1:ntrial){
      if(trialID[i]==k){
        et[1,]=eta[k,1]+eta[k,3]*Z[i]+eta[k,5]*X.01[i]+eta[k,7]*Z[i]*X.01[i]
        et[2,]=eta[k,2]+eta[k,4]*Z[i]+eta[k,6]*X.01[i]+eta[k,8]*Z[i]*X.01[i]
      }}
    lik=t(t(Y[,i]-et-mu))%*%t(Y[,i]-et-mu)
    x=x+lik
  }
  
  scale=solve((x+solve(E)))
  
  Sinv=stats::rWishart(1, n+vsigma, scale)[,,1]
  S=solve(Sinv)
  
  
  x<-array(rep(0,8*8*ntrial),dim=c(8,8,ntrial))
  y<-array(rep(0,8*ntrial),dim=c(8,ntrial))
  
  # S2 inv 
  for(i in 1:n){
    for(k in 1:ntrial){
      if(trialID[i]==k){
        x1=t(X[,,i])%*%solve(S)%*%X[,,i]
        x[,,k]=x[,,k]+x1
      }}}
  
  for(k in 1:ntrial){
    x[,,k]=x[,,k]+solve(D)
  }
  
  # gamma_i conditional
  for(i in 1:n){
    for(k in 1:ntrial){
      if(trialID[i]==k){
        y1=t(X[,,i])%*%solve(S)%*%(Y[,i]-X[,,i]%*%t(t(beta)))
        y[,k]=y[,k]+y1
      }}}
  
  for(k in 1:ntrial){
    mean=solve(x[,,k])%*%y[,k]
    eta[k,]=mvrnorm(1,mean,solve(x[,,k]))
  }
  
  # S3
  x=0
  for(i in 1:n){
    x1=t(X[,,i])%*%solve(S)%*%X[,,i]
    x=x+x1
  }
  
  # mu conditional
  y=array(rep(0,8*ntrial),dim=c(8,ntrial))
  for(i in 1:n){
    for(k in 1:ntrial){
      if(trialID[i]==k){
        ynew=t(X[,,i])%*%solve(S)%*%(Y[,i]-W[,,i]%*%t(t(eta[k,])))
        y[,k]=y[,k]+ynew
      }}}
  
  m=array(0,c(8,1))
  for(i in 1:8){
    m[i,]=sum(y[i,])
  }
  
  beta=mvrnorm(1,solve(x)%*%m,solve(x))
  
  holdbet[sim,]<-beta
  holdeta[,,sim]<-eta
  holdD[,,sim]<-D
  holdsig[,,sim]<-S
  
  sim=sim+1
  
}

keep.psiSP13 <- matrix(0, SIM, 2)
keep.psiSP123 <- matrix(0, SIM, 2)
keep.s <- matrix(0, SIM, 2)

n0.0 <- sum(X.01==0) - sum(Z[which(X.01 == 0)])
n1.0 <- sum(Z[which(X.01 == 0)])
n0.1 <- sum(X.01==1) - sum(Z[which(X.01 == 1)])
n1.1 <- sum(Z[which(X.01 == 1)])


for (i in 1:SIM){
  D <- holdD[,,i]
  D_0 <- D[3:4,3:4]
  D_1 <- matrix(c(D[3,3] + D[7,7] + 2*D[3,7],
                  D[3,4] + D[3,8] + D[4,7] + D[7,8],
                  D[3,4] + D[3,8] + D[4,7] + D[7,8],
                  D[4,4] + D[8,8] + 2*D[4,8]), 
                nrow = 2, ncol = 2, byrow = T)
  
  ss <- sqrt(holdsig[1,1,i])
  daa.0=D_0[1,1]+ss*((1/n1.0)+(1/n0.0))
  daa.1=D_1[1,1]+ss*((1/n1.1)+(1/n0.1))
  
  keep.psiSP13[i,1] <- 1 - pmvnorm(lower = -Inf, upper = 0, mean = holdbet[i,3], sigma = D_0[1,1]) -
    pmvnorm(lower = -Inf, upper = 0, mean = holdbet[i,4], sigma = D_0[2,2]) +
    2*pmvnorm(lower = -Inf, upper = c(0,0), mean = holdbet[i,3:4], sigma = D_0)
  keep.psiSP13[i, 2] <- 1 - pmvnorm(lower = -Inf, upper = 0, mean = holdbet[i,3] + holdbet[i,7], sigma = D_1[1,1]) -
    pmvnorm(lower = -Inf, upper = 0, mean = holdbet[i,4] + holdbet[i,8], sigma = D_1[2,2]) +
    2*pmvnorm(lower = -Inf, upper = c(0,0), mean = holdbet[i,3:4] + holdbet[i,7:8] , sigma = D_1)
  keep.psiSP123[i, 1] <- 1 - pmvnorm(lower = -Inf, upper = 0, mean = holdbet[i,4], sigma = D_0[2,2]) +
    pmvnorm(lower = -Inf, upper = c(0,0), mean = holdbet[i,3:4], sigma = D_0)
  keep.psiSP123[i, 2] <- 1 - pmvnorm(lower = -Inf, upper = 0, mean = holdbet[i,4] + holdbet[i,8], sigma = D_1[2,2]) +
    pmvnorm(lower = -Inf, upper = c(0,0), mean = holdbet[i,3:4] + holdbet[i,7:8] , sigma = D_1)
  
  keep.s[i,1]=(holdbet[i,3])-(daa.0/D_0[1,2])*((holdbet[i,4])+qnorm(0.95)*sqrt(D_0[2,2]-(D_0[1,2]^2/daa.0)))
  keep.s[i,2]=(holdbet[i,3] + holdbet[i,7])-(daa.1/D_1[1,2])*((holdbet[i,4] + holdbet[i,8])+qnorm(0.95)*sqrt(D_1[2,2]-(D_1[1,2]^2/daa.1)))
  
}

sp13_x0 <- keep.psiSP13[8000:SIM,1]
sp13_x1 <- keep.psiSP13[8000:SIM,2]
sp123_x0 <- keep.psiSP123[8000:SIM,1]
sp123_x1 <- keep.psiSP123[8000:SIM,2]

s_x0 <- keep.s[8000:SIM,1]
s_x1 <- keep.s[8000:SIM,2]


estimates <- c(sp13_x0_mean = mean(sp13_x0),
               sp13_x0_se = sqrt(var(sp13_x0)),
               sp13_x0_lower = quantile(sp13_x0,  probs = 0.025),
               sp13_x0_upper =quantile(sp13_x0,  probs = 0.975),
               sp13_x1_mean = mean(sp13_x1),
               sp13_x1_se = sqrt(var(sp13_x1)),
               sp13_x1_lower = quantile(sp13_x1,  probs = 0.025),
               sp13_x1_upper =quantile(sp13_x1,  probs = 0.975),
               sp123_x0_mean = mean(sp123_x0),
               sp123_x0_se = sqrt(var(sp123_x0)),
               sp123_x0_lower = quantile(sp123_x0,  probs = 0.025),
               sp123_x0_upper =quantile(sp123_x0,  probs = 0.975),
               sp123_x1_mean = mean(sp123_x1),
               sp123_x1_se = sqrt(var(sp123_x1)),
               sp123_x1_lower = quantile(sp123_x1,  probs = 0.025),
               sp123_x1_upper =quantile(sp123_x1,  probs = 0.975),
               
               s_x0_mean = mean(s_x0),
               s_x0_se = sqrt(var(s_x0)),
               s_x0_lower = quantile(s_x0,  probs = 0.025),
               s_x0_upper =quantile(s_x0,  probs = 0.975),
               s_x1_mean = mean(s_x1),
               s_x1_se = sqrt(var(s_x1)),
               s_x1_lower = quantile(s_x1,  probs = 0.025),
               s_x1_upper = quantile(s_x1,  probs = 0.975))

