# Load libraries
library(statsecol)
library(jagsUI)
library(MCMCvis)
library(ggplot2)
library(statsecol) 
# load data
data("voles")

# Specify model in BUGS language
sink("WVSSM.txt")
cat("
model{
  # Priors and constraints
  psi1 ~ dunif(0, 1)    
  beta0g ~ dnorm(0, 0.01)
  beta1g ~ dnorm(0, 0.01)
  beta0e ~ dnorm(0, 0.01)
  beta1e ~ dnorm(0, 0.01)
  mu_p ~ dnorm(0, 0.01)
  sig_p ~ dunif(0, 1)
  
  # Likelihood - State process
  for(i in 1:n.sites){
    z[i,1] ~ dbern(psi1) #initial
    for(t in 2:n.years){ #subsequent
      z[i,t] ~ dbern((1 - z[i,t-1])*gamma[i,t-1] + z[i,t-1]*(1-epsilon[i,t-1]))
    logit(gamma[i,t-1]) <- beta0g + beta1g*Connectivity[i]
    logit(epsilon[i,t-1]) <- beta0e + beta1e*Length[i]
    }
  }
  
  # Likelihood - Observation process
  for(i in 1:n.sites){
     q[i,1] ~ dnorm(mu_p,pow(sig_p,-2))
     logit(p[i,1]) <- q[i,1]   
     for(t in 2:n.years){
        q[i,t] ~ dnorm(mu_p,pow(sig_p,-2))
        logit(p[i,t]) <- q[i,t]
     }
     for(t in 1:n.years){
      y[i,t] ~ dbin(p[i,t]*z[i,t],n.visits[i,t])

      ysim[i,t] ~ dbin(p[i,t]*z[i,t],n.visits[i,t])
      yexp[i,t] <- p[i,t]*n.visits[i,t]*z[i,t] + 0.001
      x2.obs[i,t] <- pow((y[i,t] - yexp[i,t]),2)/(yexp[i,t])    # for observed data
      x2.sim[i,t] <- pow((ysim[i,t] - yexp[i,t]),2)/(yexp[i,t]) # for 'ideal' data
      ft.obs[i,t] <- pow(sqrt(y[i,t])   - sqrt(yexp[i,t]),2)    # for observed data
      ft.sim[i,t] <- pow(sqrt(ysim[i,t]) - sqrt(yexp[i,t]),2)   # for 'ideal' data
    }
  }
  Chi2.obs <- sum(x2.obs[,])
  Chi2.sim <- sum(x2.sim[,])
  Chi2.ratio <- x2.obs/x2.sim
  FT.obs <- sum(ft.obs[,])
  FT.sim <- sum(ft.sim[,])
  FT.ratio <- FT.obs/FT.sim
  
  # Derived quantities
  for(t in 1:n.years){
    propocc[t] <- sum(z[,t])/n.sites
  }
}

",fill = TRUE)
sink()

# Bundle data
WVdata <- list(y = matrix(c(voles$y1,voles$y2,voles$y3,voles$y4),114,4), 
               n.years = 4, n.sites = nrow(voles), 
               n.visits = matrix(c(voles$j1,voles$j2,voles$j3,voles$j4),114,4),
               Length = voles$Length, Connectivity = voles$Connectivity)

# Initial values
WVinits <- function(){
  list(psi1 = runif(1, 0, 1), 
       beta0g = runif(1, -2, 2),
       beta1g = runif(1, -2, 2),
       beta0e = runif(1, -2, 2),
       beta1e = runif(1, -2, 2),
       mu_p = runif(1, -2, 2),
       sig_p = runif(1, 0, 1),
       z = ifelse(matrix(c(voles$y1,voles$y2,voles$y3,voles$y4),114,4)>0,1,0))
}

# Parameters monitored
WVparms <- c("gamma","epsilon","beta0g","beta1g","beta0e","beta1e","mu_p","sig_p",
             "propocc","p","Chi2.obs", "Chi2.sim", "FT.obs", "FT.sim")

# MCMC settings
ni <- 10000
nt <- 4
nb <- 1000
nc <- 3

WVout <- jags(data = WVdata,
                  inits = WVinits,
                  parameters.to.save = WVparms,
                  model.file = "WVSSM.txt",
                  n.chains = nc,
                  n.iter = ni,
                  n.burnin = nb,
                  n.thin = nt) 

# Check convergence of gamma
MCMCtrace(WVout,                 #the fitted model
          params = WVparms[1],   #core model parameters
          iter = ni,                 #plot all iterations
          pdf = FALSE,           #DON'T write to a PDF
          ind=TRUE)                
# Check convergence of epsilon
MCMCtrace(WVout,                 #the fitted model
          params = WVparms[2],   #core model parameters
          iter = ni,                 #plot all iterations
          pdf = FALSE,               #DON'T write to a PDF
          ind=TRUE)   
# Check convergence of beta0g
MCMCtrace(WVout,                 #the fitted model
          params = WVparms[3],   #core model parameters
          iter = ni,                 #plot all iterations
          pdf = FALSE,               #DON'T write to a PDF
          ind=TRUE)       
# Check convergence of beta1g
MCMCtrace(WVout,                 #the fitted model
          params = WVparms[4],   #core model parameters
          iter = ni,                 #plot all iterations
          pdf = FALSE,               #DON'T write to a PDF
          ind=TRUE)   
# Check convergence of beta0e
MCMCtrace(WVout,                 #the fitted model
          params = WVparms[5],   #core model parameters
          iter = ni,                 #plot all iterations
          pdf = FALSE,               #DON'T write to a PDF
          ind=TRUE) 
# Check convergence of beta1e
MCMCtrace(WVout,                 #the fitted model
          params = WVparms[6],   #core model parameters
          iter = ni,                 #plot all iterations
          pdf = FALSE,               #DON'T write to a PDF
          ind=TRUE)    
# Check convergence of mu_p
MCMCtrace(WVout,                 #the fitted model
          params = WVparms[7],   #core model parameters
          iter = ni,                 #plot all iterations
          pdf = FALSE,               #DON'T write to a PDF
          ind=TRUE)  
# Check convergence of sig_p
MCMCtrace(WVout,                 #the fitted model
          params = WVparms[8],   #core model parameters
          iter = ni,                 #plot all iterations
          pdf = FALSE,               #DON'T write to a PDF
          ind=TRUE)   

# Check convergence of propocc
MCMCtrace(WVout,                 #the fitted model
          params = WVparms[10],   #core model parameters
          iter = ni,                 #plot all iterations
          pdf = FALSE,               #DON'T write to a PDF
          ind=TRUE)   

# Check Rhat and summaries of all parameters and make table
MCMC_Param_Summary<-MCMCsummary(WVout,
            params = WVparms[3:9]) #out parameters of interest

colnames(MCMC_Param_Summary) <- c("Mean", "Sd", "Lower", "Medium", "Upper",
                                  "Rhat", "n.eff")
MCMC_Param_Summary <- rownames_to_column(MCMC_Param_Summary, var = "Parameter")
MCMC_Param_Summary[,2:6] <- MCMC_Param_Summary[,2:6] %>%
  mutate(across(where(is.numeric), round, digits = 3))

formattable(MCMC_Param_Summary, 
            align =c("l","c","c","c","c", "c", "c", "c", "r"), 
            list('Parameter' = formatter(
              "span", style = ~ style(color = "grey",font.weight = "bold")) 
            ))

# Posterior predictive checking
obs.chi <- WVout$sims.list$Chi2.obs
sim.chi <- WVout$sims.list$Chi2.sim
obs.ft <- WVout$sims.list$FT.obs
sim.ft <- WVout$sims.list$FT.sim

par(mfrow=c(1,2), oma=c(0,0,0,0), mar=c(4,4,1,1))
plot(sim.chi~obs.chi,col=adjustcolor(ifelse(obs.chi>sim.chi, "maroon","darkgreen"),0.25), 
     pch=16, asp=1, xlab="Observed χ^2 distribution", ylab="Simulated χ^2 distribution")
abline(0,1, lwd=2)

plot(sim.ft~obs.ft,col=adjustcolor(ifelse(obs.ft>sim.ft, "maroon","darkgreen"),0.25), 
     pch=16, asp=1, xlab="Observed FT distribution", ylab="Simulated FT distribution")
abline(0,1, lwd=2)

# Bayesian p-values
mean(obs.chi>sim.chi)
mean(obs.ft>sim.ft)

# Plot extinction and colonisation probabilities against predictors
df <- data.frame(Connectivity = voles$Connectivity,
                     Mean = WVout$mean$gamma[,1],
                     Lower = WVout$q2.5$gamma[,1],
                     Upper = WVout$q97.5$gamma[,1],
                     Length = voles$Length,
                     Meane = WVout$mean$epsilon[,1],
                     Lowere = WVout$q2.5$epsilon[,1],
                     Uppere = WVout$q97.5$epsilon[,1] )

ggplot(data=df) + 
  geom_line(aes(x=Connectivity, y=Mean),size=2,colour="blue") +
  ylab("Colonisation Probability γ")+
  geom_ribbon(aes(x=Connectivity, y=Mean, ymin=Lower, ymax = Upper), 
              fill="black", alpha=0.2) 

ggplot(data=df) + 
  geom_line(aes(x=Length, y=Meane),size=2,colour="red") +
  ylab("Extinction Probability ε")+
  geom_ribbon(aes(x=Length, y=Meane, ymin=Lowere, ymax = Uppere), 
              fill="black", alpha=0.2) + xlab("Length (km)")
