library('nimble')


##############
# Question 1 #
##############

# DATA PREPARATION

d<-c(0,62.5,125,250,500)
N<-c(282,225,290,261,141)
y<-c(67,34,193,250,141)
nr<-5

model.data <- list('y' = y)
model.constant <- list('d' = d, 'N' = N,'nr'=nr)

# DEFINE INITIAL VALUES

model.inits_1 <- list(alpha=-1,beta=0.010)
model.inits_2 <- list(alpha=1,beta=-0.010)
model.inits<-list(model.inits_1,model.inits_2)

# MODEL SPECIFICATION 

set.seed(20)
CBI_DR <- nimbleCode(
  {
    # Specification likelihood
    for (i in 1:nr)
    {
      p[i] <- expit(alpha+beta*d[i])
      y[i] ~ dbin(p[i],N[i])
      
    }
    # Prior specification
    alpha ~ dnorm(0,0.001)
    beta ~ dnorm(0,0.001)
    
  })


# SET UP MODEL OBJECT 

out<-nimbleMCMC(code=CBI_DR,data=model.data,constants=model.constant, thin=15,niter=100000,nburnin=2000,nchains=2,summary=TRUE,inits=model.inits)
out$summary

# CONVERGENCE
# plots

library(coda)
par(mfrow=c(1,2))

traceplot(as.mcmc(out$samples$chain1))
densplot(as.mcmc(out$samples$chain1))
acf(as.mcmc(out$samples$chain1))


traceplot(as.mcmc(out$samples$chain2))
densplot(as.mcmc(out$samples$chain2))
acf(as.mcmc(out$samples$chain2))

# Gelman and Rubin convergence diagnostic
mclist<-mcmc.list(as.mcmc(out$samples$chain1),as.mcmc(out$samples$chain2))
gelman.diag(mclist)

# SENSITIVITY ANALYSIS

CBI_DR_t <- nimbleCode(
  {
    # Specification likelihood
    for (i in 1:nr)
    {
      p[i] <- expit(alpha+beta*d[i])
      y[i] ~ dbin(p[i],N[i])
    }
    # Prior specification
    alpha ~ dnorm(0,0.0001)
    beta ~ dt(0,0.0001,4)
    
  })


out_t<-nimbleMCMC(code=CBI_DR_t,data=model.data,constants=model.constant, thin=15,niter=100000,nburnin=2000,nchains=2,summary=TRUE,inits=model.inits)
out_t$summary


traceplot(as.mcmc(out_t$samples$chain1))
densplot(as.mcmc(out_t$samples$chain1))
acf(as.mcmc(out_t$samples$chain1))


traceplot(as.mcmc(out_t$samples$chain2))
densplot(as.mcmc(out_t$samples$chain2))
acf(as.mcmc(out_t$samples$chain2))

# Gelman and Rubin convergence diagnostic
mclist_t<-mcmc.list(as.mcmc(out_t$samples$chain1),as.mcmc(out_t$samples$chain2))
gelman.diag(mclist_t)

# Combine samples
samples=rbind(out$samples$chain1,out$samples$chain2)
samples_t=rbind(out_t$samples$chain1,out_t$samples$chain2)
densplot(as.mcmc(samples))
densplot(as.mcmc(samples_t))


library(MCMCvis)
MCMCsummary(out$samples,round=2)

# traceplot and posterior distributions, Rhat,
# number of effective samples
MCMCtrace(out$samples,
          pdf=FALSE,
          ind=TRUE,
          Rhat=TRUE,
          n.eff=TRUE)


head(output$samples$chain1)
head(output$samples$chain2)

MCMCsummary(out_t$samples,round=2)

# traceplot and posterior distributions, Rhat,
# number of effective samples
MCMCtrace(out_t$samples,
          pdf=FALSE,
          ind=TRUE,
          Rhat=TRUE,
          n.eff=TRUE)


head(out_t$samples$chain1)
head(out_t$samples$chain2)


library(mcmcse)

# monte carlo error smaller than 5% of the posterior standard error
mcse(as.mcmc(out$samples$chain1[,1]))
mcse(as.mcmc(out$samples$chain1[,2]))
mcse(as.mcmc(out$samples$chain2[,1]))
mcse(as.mcmc(out$samples$chain2[,2]))
mcse(as.mcmc(samples[,1]))
mcse(as.mcmc(samples[,2]))

# monte carlo error smaller than 5% of the posterior standard error
mcse(as.mcmc(out_t$samples$chain1[,1]))
mcse(as.mcmc(out_t$samples$chain1[,2]))
mcse(as.mcmc(out_t$samples$chain2[,1]))
mcse(as.mcmc(out_t$samples$chain2[,2]))
mcse(as.mcmc(samples_t[,1]))
mcse(as.mcmc(samples_t[,2]))


# Posterior HPD
HPDinterval(as.mcmc(out$samples$chain1))
HPDinterval(as.mcmc(out$samples$chain2))
HPDinterval(as.mcmc(samples))

HPDinterval(as.mcmc(out_t$samples$chain1))
HPDinterval(as.mcmc(out_t$samples$chain2))
HPDinterval(as.mcmc(samples_t))





# Posterior dose-response relationship

alpha=out$summary$all.chains[1,1]
beta=out$summary$all.chains[2,1]

x=seq(0,500,0.1)
par(mfrow=c(1,1))
lines(x,expit(alpha+beta*x), xlab="dose",ylab="probability",col='grey')
points(d,y/N,pch=19)


alpha_t=out_t$summary$all.chains[1,1]
beta_t=out_t$summary$all.chains[2,1]

x=seq(0,500,0.1)
par(mfrow=c(1,1))
plot(x,expit(alpha_t+beta_t*x), xlab="dose",ylab="probability",col='grey')
points(d,y/N,pch=19)

# Safe level of exposure

q=0.01
alpha_samples=samples[,1]
beta_samples=samples[,2]
q_star=q*(1-expit(alpha_samples))+expit(alpha_samples)
BMD=(logit(q_star)-alpha_samples)/beta_samples
summary(BMD)
plot(density(BMD))
HPDinterval(as.mcmc(BMD))

q=0.01
alpha_samples_t=samples_t[,1]
beta_samples_t=samples_t[,2]
q_star=q*(1-expit(alpha_samples_t))+expit(alpha_samples_t)
BMD_t=(logit(q_star)-alpha_samples_t)/beta_samples_t
summary(BMD_t)
plot(density(BMD_t))
HPDinterval(as.mcmc(BMD_t))

# Number of malformations
y_pred=vector()
for (i in 1:13066){
  y_pred[i]=rbinom(1,240,expit(alpha_samples[i]+beta_samples[i]*100))
}
summary(y_pred)
plot(density(y_pred))
HPDinterval(as.mcmc(y_pred))

y_pred_t=vector()
for (i in 1:13066){
  y_pred_t[i]=rbinom(1,240,expit(alpha_samples_t[i]+beta_samples_t[i]*100))
}
summary(y_pred_t)
plot(density(y_pred_t))
HPDinterval(as.mcmc(y_pred_t))


##############
# Question 2 #
##############

# DATA PREPARATION

N<-c(272,87,322,176,94,387,279,194,65,110,266,397,152,231)
Z<-c(17,15,71,17,9,23,78,59,47,34,43,57,29,17)
type<-c(1,1,1,1,1,1,0,0,0,0,0,0,0,0)
nr<-14

model.data <- list('Z' = Z)
model.constant <- list('N' = N, 'type' = type,'nr'=nr)

# DEFINE INITIAL VALUES

model.inits_1 <- list(beta_0=-10,beta_1=10)
model.inits_2 <- list(beta_0=10,beta_1=-10)
model.inits<-list(model.inits_1,model.inits_2)

# MODEL SPECIFICATION 

set.seed(20)
CBI_AP <- nimbleCode(
  {
    # Specification likelihood
    for (i in 1:nr)
    {
      p[i] <- expit(beta_0+beta_1*type[i])
      Z[i] ~ dbin(p[i],N[i])
    }
    # Prior specification
    beta_0 ~ dnorm(0,0.001)
    beta_1 ~ dnorm(0,0.001)
  })

# SET UP MODEL OBJECT 
out<-nimbleMCMC(code=CBI_AP,data=model.data,constants=model.constant,niter=50000,nburnin=1000,thin=10,nchains=2,summary=TRUE,inits=model.inits)
out$summary


# Plots

library(coda)
par(mfrow=c(1,2))


traceplot(as.mcmc(out$samples$chain1))
densplot(as.mcmc(out$samples$chain1))
acf(as.mcmc(out$samples$chain1))


traceplot(as.mcmc(out$samples$chain2))
densplot(as.mcmc(out$samples$chain2))
acf(as.mcmc(out$samples$chain2))


MCMCtrace(out$samples,
          pdf=FALSE,
          ind=TRUE,
          Rhat=TRUE,
          n.eff=TRUE)

# Gelman and Rubin convergence diagnostic
mclist<-mcmc.list(as.mcmc(out$samples$chain1),as.mcmc(out$samples$chain2))
gelman.diag(mclist)

# Posterior summary statistics
samples=rbind(out$samples$chain1,out$samples$chain2)
HPDinterval(as.mcmc(out$samples$chain1))
HPDinterval(as.mcmc(out$samples$chain2))
HPDinterval(as.mcmc(samples))


# monte carlo error smaller than 5% of the posterior standard error
mcse(as.mcmc(out$samples$chain1[,1]))
mcse(as.mcmc(out$samples$chain1[,2]))
mcse(as.mcmc(out$samples$chain2[,1]))
mcse(as.mcmc(out$samples$chain2[,2]))
mcse(as.mcmc(samples[,1]))
mcse(as.mcmc(samples[,2]))


par(mfrow=c(1,2))
densplot(as.mcmc(samples))

# apparent prevalence
beta_0_sample=samples[,1]
beta_1_sample=samples[,2]
p0_sample=expit(beta_0_sample)
p1_sample=expit(beta_0_sample+beta_1_sample)
summary(p0_sample)
summary(p1_sample)
par(mfrow=c(1,1))
plot(density(p0_sample),xlim=c(0,0.3),ylim=c(0,50),col="blue", main='Apparent animal prevalence')
lines(density(p1_sample),col="red")
legend(x=0.01,y=40,legend=c("type 0","type 1"),col=c("blue","red"),lty=1)

difference=p1_sample-p0_sample
summary(difference)

HPDinterval(as.mcmc(difference))
HPDinterval(as.mcmc(p0_sample))
HPDinterval(as.mcmc(p1_sample))

# prior on specificity and sensitivity
# Calculating alpha and beta for beta distribution

# Sensitivity

func_sens<- function(alpha) pbeta(0.90,alpha,(0.15*alpha+0.7)/0.85)-pbeta(0.82,alpha,(0.15*alpha+0.7)/0.85)-0.95
alpha_sens <- uniroot(func_sens,c(0.01,1000), extendInt="yes")$root
beta_sens=(0.15*alpha_sens+0.7)/0.85
(alpha_sens-1)/(alpha_sens+beta_sens-2)
pbeta(0.90,alpha_sens,beta_sens)-pbeta(0.82,alpha_sens,beta_sens)


# Specificity
func_spec <- function(alpha) pbeta(0.97,alpha,(0.05*alpha+0.9)/0.95)-pbeta(0.90,alpha,(0.05*alpha+0.9)/0.95)-0.95
alpha_spec <- uniroot(func_spec,c(0.01,1000), extendInt="yes")$root
beta_spec=(0.05*alpha_spec+0.9)/0.95
(alpha_spec-1)/(alpha_spec+beta_spec-2)
pbeta(0.97,alpha_spec,beta_spec)-pbeta(0.90,alpha_spec,beta_spec)

x=seq(0.7,1,0.001)
plot(x,dbeta(x,alpha_sens,beta_sens),col="dark green",type='l',ylab='',main='Sensitivity and Specificity',ylim=c(0,30))
lines(x,dbeta(x,alpha_spec,beta_spec),col="dark blue",type='l')
legend(x=0.7,y=25,col=c("dark green", "dark blue"), legend=c("Sensitivity",'Specificity'),lty=1)


hpdbeta <- function(alpha,beta)
{
  p2 <- alpha
  q2 <- beta
  f <- function(x,p=p2,q=q2){
    b<-qbeta(pbeta(x,p,q)+0.95,p,q);(dbeta(x,p,q)-dbeta(b,p,q))^2}
  hpdmin <- optimize(f,lower=0,upper=qbeta(0.05,p2,q2),p=p2,q=q2)$minimum
  hpdmax <- qbeta(pbeta(hpdmin,p2,q2)+0.95,p2,q2)
  return(c(hpdmin,hpdmax))
}
hpdbeta(alpha_sens,beta_sens)
hpdbeta(alpha_spec,beta_spec)

mean_sens=alpha_sens/(alpha_sens+beta_sens)
mean_spec=alpha_spec/(alpha_spec+beta_spec)

qbeta(0.25,alpha_sens,beta_sens)
qbeta(0.75,alpha_sens,beta_sens)
qbeta(0.25,alpha_spec,beta_spec)
qbeta(0.75,alpha_spec,beta_spec)

# posterior true animal prevalence

Se<-rbeta(8000,alpha_sens,beta_sens)
Sp<-rbeta(8000,alpha_spec,beta_spec)

pi_0=(p0_sample+Sp-1)/(Se+Sp-1)
pi_1=(p1_sample+Sp-1)/(Se+Sp-1)
difference=pi_1-pi_0
mean(pi_0)
mean(pi_1)

plot(density(pi_0), xlim=c(0,0.30),ylim=c(0,30),col="blue",main='True animal prevalence')
lines(density(pi_1),col="red")
legend(x=0.01,y=30,legend=c("type 0","type 1"),col=c("blue","red"),lty=1)
summary(pi_0)
summary(pi_1)
summary(difference)
HPDinterval(as.mcmc(pi_0))
HPDinterval(as.mcmc(pi_1))
HPDinterval(as.mcmc(difference))

# alternative
CBI_AP_2 <- nimbleCode(
  {
    # Specification likelihood
    for (i in 1:nr)
    {
      p[i] <- expit(beta_0+beta_1*type[i])
      Z[i] ~ dbin(p[i],N[i])
    }
    # Prior specification
    beta_0 ~ dnorm(0,0.001)
    beta_1 ~ dnorm(0,0.001)
    Se ~ dbeta(387.61,69.23)
    Sp ~ dbeta(188.75,10.88)
  })

# SET UP MODEL OBJECT 
out2<-nimbleMCMC(code=CBI_AP_2,data=model.data,constants=model.constant,niter=5000,nburnin=1000,nchains=2,summary=TRUE,inits=model.inits)
out2$summary

samples2 <- rbind(out2$samples$chain1,out2$samples$chain2)

se_sample = samples2[,1]
sp_sample = samples2[,2]
beta_0_sample=samples2[,3]
beta_1_sample=samples2[,4]
p0_sample=expit(beta_0_sample)
p1_sample=expit(beta_0_sample+beta_1_sample)

Pi_0_sampled1 <- (p0_sample+ sp_sample - 1)/(se_sample+sp_sample-1)
Pi_1_sampled1 <- (p1_sample + sp_sample - 1)/(se_sample+sp_sample-1)

Pi_0_sampled1 <- (p0_sample+ sp_sample - 1)/(se_sample+sp_sample-1)
Pi_1_sampled1 <- (p1_sample + sp_sample - 1)/(se_sample+sp_sample-1)

mean(Pi_0_sampled1)
mean(Pi_1_sampled1)

par(mfrow=c(1,1))
plot(density(Pi_0_sampled1),xlim=c(0,0.3),ylim=c(0,30),col="blue",
     main = "Posterior distribution of Pi")  
lines(density(Pi_1_sampled1),col="red")
legend(x="topleft",
       cex = 0.6,
       box.col = 11,
       box.lty = 2,
       box.lwd = 1,
       legend = c("type 0","type 1"),
       title = "true animal prevalences Pi",
       title.col = c(1),
       title.adj = 0.5,
       lty = c(1,1),
       col = c("blue","red"),
       lwd = 1)
