##This R code is meant to analyze the feeding rate of uninfected, singly 
##infected, and uninfected animals across age groups. Based on feeding rates, 
##We also analyze host susceptibility across those categories

#Load libraries
library(here)
library(ggplot2)
library(tidyr)
require(plyr)
require(dplyr)
require(lubridate)
require(MASS)
require(fitdistrplus)

#Import feeding rate data
feeding <- read.csv("feeding_rate_data.csv")
feeding<-na.omit(feeding)
feeding$treatment <- as.factor(feeding$treatment)

#Import spore count data
spores <- read.csv('spore_counts.csv')
spores$coinf <- as.factor(spores$coinf)
spores$coinf_treat <- as.factor(spores$coinf_treat)
spores$status <- as.factor(spores$status)
spores$treatment <- as.factor(spores$treatment)

#Get rid of contaminated individuals
spores <- spores %>%
  filter(missing == 0) %>%
  filter(cont == 0)

#Merge feeding rate data with spores data, so we can relate feeding rate and whether
##animals got infected
feeding <- merge(feeding,spores,by=c("treatment","replicate"),all.x = TRUE)

##So that we can eventually add in only successfully infected individuals, merge
##With infection data


######1. feeding rate vs. age in healthy individuals

feeding.healthy <- feeding[ which(feeding$treatment=='1'), ]

#We don't trust middle ages

feeding.healthy <- feeding.healthy[ which(feeding.healthy$age.x < 12 | feeding.healthy$age.x > 18), ]


######2. feeding rate vs. age purely in the presence of metsch only or past only

feeding.mexp <- feeding[ which(feeding$treatment=='3' | feeding$treatment=='5' | feeding$treatment=='7' | feeding$treatment=='9' | feeding$treatment=='11' | feeding$treatment=='13'), ]
feeding.mexp <- feeding.mexp[ which(feeding.mexp$age.x=='4' | feeding.mexp$age.x=='20'), ]

feeding.pexp <- feeding[ which(feeding$treatment=='2' | feeding$treatment=='4' | feeding$treatment=='6' | feeding$treatment=='8' | feeding$treatment=='10' | feeding$treatment=='12'), ]
feeding.pexp <- feeding.pexp[ which(feeding.pexp$age.x=='4' | feeding.pexp$age.x=='20'), ]

##Add in those individuals who were exposed and successfully infected previously

feeding.mexp.pinf <- feeding[ which(feeding$treatment=='4' | feeding$treatment=='12'), ]
feeding.mexp.pinf <- feeding.mexp.pinf[ which(feeding.mexp.pinf$age.x=='7' | feeding.mexp.pinf$age.x=='23'), ]
feeding.mexp.pinf <- feeding.mexp.pinf[ which(feeding.mexp.pinf$average_past>0), ]

feeding.pexp.minf <- feeding[ which(feeding$treatment=='5' | feeding$treatment=='13'), ]
feeding.pexp.minf <- feeding.pexp.minf[ which(feeding.pexp.minf$age.x=='7' | feeding.pexp.minf$age.x=='23'), ]
feeding.pexp.minf <- feeding.pexp.minf[ which(feeding.pexp.minf$average_met>0), ]


healthy.mexp.pexp <- rbind(feeding.healthy,feeding.mexp,feeding.pexp,feeding.mexp.pinf,feeding.pexp.minf)

##OK, filter out negative feeding rates
healthy.mexp.pexp <- healthy.mexp.pexp[which(healthy.mexp.pexp$filtration > 0),]

##OK, now we measure feeding rates
#full model. No square term because only had two feeding dates (essentially)
##This is what we report in the paper
feeding_age_lm <- lm(filtration ~ age.x * status.x, data=healthy.mexp.pexp)
summary(feeding_age_lm)

#Model for significance- just compare animals at first exposure to unexposed animals
#This is to parameterize the model- doesn't actually make any difference
mexp.pexp.feed.param <- healthy.mexp.pexp[which(healthy.mexp.pexp$status.x != 'mets_after_past' & healthy.mexp.pexp$status.x != 'past_after_mets'),]
feeding_age_lm <- lm(filtration ~ age.x + age.x:status.x, data=mexp.pexp.feed.param)
summary(feeding_age_lm)


##Our glm results indicate that that individuals change feeding rate when 
##exposed to either pasteuria or metschnikowia, but that this dissapears if they
##have previously been infected. So for plotting, one color for uninfected, one
##color for not previously infected and exposed to metsch, one for not
##previously infected and exposed to past

feeding.plot <- healthy.mexp.pexp[which(healthy.mexp.pexp$status.x != "cont_after_mets_after_past" | healthy.mexp.pexp$status.x != "cont_after_past_after_mets"),]

feeding.plot$status.x <- as.factor(feeding.plot$status.x)

levels(feeding.plot$status.x) <- c("Unexposed", "Exposed to Fungus","Exposed to Fungus, Prior Infection",
                                   "Exposed to Bacteria","Exposed to Bacteria, Prior Infection")

feeding.plot.lim <- feeding.plot[which(feeding.plot$status.x == "Unexposed" | feeding.plot$status.x == "Exposed to Bacteria" | feeding.plot$status.x == "Exposed to Fungus"),]

tiff(filename = "feeding_figure.tiff", width = 2500, height = 1500, pointsize = 12,res=300)

ggplot(feeding.plot,aes(y=filtration,x=age.x,color=factor(status.x),fill=factor(status.x),shape=factor(status.x))) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(colour = "black", size=1)) + 
  geom_jitter(width=1.5,height=0,alpha = 1.0,size=2) +
  stat_smooth(data = feeding.plot.lim, method="lm") +
  stat_smooth(data = feeding.plot.lim, method="lm", se = FALSE) +
  #coord_cartesian(ylim = c(0,250)) +
  scale_colour_manual(values = c("#000000", "#E69F00","#009E73","#56B4E9","#F0E442")) +
  scale_fill_manual(values = c("#000000", "#E69F00","#009E73","#56B4E9","#F0E442")) +
  xlab("Age (Days)") + ylab("Filtration Rate (ml/hour)") + #ggtitle("Metsch fitness inc. both infectivity and spore yield") +
  theme(axis.title.x = element_text(face="bold", size=16),axis.text.x  = element_text(size=16)) + 
  theme(axis.title.y = element_text(face="bold", size=16),axis.text.y  = element_text(size=16)) +
  theme(legend.title=element_blank())

dev.off()

#OK, let's get spores consumed per age. 
#OF past alone, past after metsch infection, metsch alone, and metsch after past infection
#This is those four classes
mexp.pexp <- rbind(feeding.mexp,feeding.pexp,feeding.mexp.pinf,feeding.pexp.minf)
#put in row for 1-exp(-(filt/15ml)*24hour) then multiple by spore dose
mexp.pexp$prop.filt <- 1-exp(-(mexp.pexp$filtration/15)*24)

mexp.pexp$consumed <- 0
for(i in 1:length(mexp.pexp$treatment)){
  if(mexp.pexp$status.x[i]=='mets'){mexp.pexp$consumed[i] <- mexp.pexp$prop.filt[i] * 15 * 200}
  if(mexp.pexp$status.x[i]=='past'){mexp.pexp$consumed[i] <- mexp.pexp$prop.filt[i] * 15 * 1000}
  if(mexp.pexp$status.x[i]=='mets_after_past'){mexp.pexp$consumed[i] <- mexp.pexp$prop.filt[i] * 15 * 200}
  if(mexp.pexp$status.x[i]=='past_after_mets'){mexp.pexp$consumed[i] <- mexp.pexp$prop.filt[i] * 15 * 1000}
}

#Remove individuals with negative filtration rates
mexp.pexp <- mexp.pexp[ which(mexp.pexp$consumed>0), ]

##OK, now have number of spores consumed per individual
##How do we actually do it. 
##GLM, family binomial
##treatments- past, mets, mets_after_past, past_after_mets
##times- only young and old
##a glm for probability of invectivity will give us prob based on just age and status
##but we want to incorporate amount of spores eaten

#only pasteuria susceptibility
mexp.pexp.p <- mexp.pexp[ which(mexp.pexp$status.x == 'past' | mexp.pexp$status.x == 'past_after_mets'), ]
mexp.pexp.p <- na.omit(mexp.pexp.p)
for(i in 1:length(mexp.pexp.p$status.x)){
  if(mexp.pexp.p$status.x[i] == "past"){mexp.pexp.p$status.x[i] = 0}
  if(mexp.pexp.p$status.x[i] == "past_after_mets"){mexp.pexp.p$status.x[i] = 1}
  
}
mexp.pexp.p$status.x <- as.numeric(mexp.pexp.p$status.x)


##OK, I don't think linear models are working. What if we do likelihood 
##based estimates
##THat worked, but now we need to randomize start parameters
library(deSolve)
library(ggplot2)
library(bbmle)
library(MASS)
library(bbmle)



##Calculate neg log likelihood for a given parameter set
sibr_nll = function(infect, age.int, status.int, stat.age){
  
  # Setting initial conditions and parameters just like before
  ode_res = 1-exp(-mexp.pexp.p$consumed * (infect + age.int * mexp.pexp.p$age.x
                                           + status.int * mexp.pexp.p$status.x
                                           + stat.age * mexp.pexp.p$status.x * mexp.pexp.p$age.x))
  for(i in 1:length(ode_res)){
    if(ode_res[i]>0.999){ode_res[i]=0.999}
    if(ode_res[i]<0.001){ode_res[i]=0.001}
  }
  # Removing the initial condition
  nll = -1*sum(dbinom(mexp.pexp.p$past_inf, size = 1, prob = ode_res, log=TRUE))
  return(nll)
  
}
past.explore <- matrix(0,21^4,5)
index <- 1
for(w in -10:10){
  for(x in -10:10){
    for(y in -10:10){
      for(z in -10:10){
        past.explore[index,1]<-(w/50)^3
        past.explore[index,2]<-(x/50)^3
        past.explore[index,3]<-(y/50)^3
        past.explore[index,4]<-(z/50)^3
        past.explore[index,5]<-sibr_nll((w/50)^3,(x/50)^3,(y/50)^3,(z/50)^3)
        index <- index + 1
      }
    }
  }
}


min.ind <- which.min(past.explore[,5])


start_params = c(past.explore[min.ind,1],past.explore[min.ind,2],past.explore[min.ind,3], past.explore[min.ind,4])


# Get the mle estimates for my parameters
fit_all_p = mle2(sibr_nll, 
               start=list(infect=start_params[1],
                          age.int=start_params[2],
                          status.int=start_params[3],
                          stat.age=start_params[4]),
               method="Nelder-Mead")
               
summary(fit_all_p)               

##visualize




#############NOW FOR MESCH###############3
mexp.pexp.m <- mexp.pexp[ which(mexp.pexp$status.x == 'mets' | mexp.pexp$status.x == 'mets_after_past'), ]
mexp.pexp.m <- na.omit(mexp.pexp.m)
for(i in 1:length(mexp.pexp.m$status.x)){
  if(mexp.pexp.m$status.x[i] == "mets"){mexp.pexp.m$status.x[i] = 0}
  if(mexp.pexp.m$status.x[i] == "mets_after_past"){mexp.pexp.m$status.x[i] = 1}
  
}
mexp.pexp.m$status.x <- as.numeric(mexp.pexp.m$status.x)

##OK, I don't think linear models are working. What if we do likelihood 
##based estimates

##Calculate neg log likelihood for a given parameter set
sibr_nll = function(infect, age.int, status.int, stat.age){
  
  # Setting initial conditions and parameters just like before
  ode_res = 1-exp(-mexp.pexp.m$consumed * (infect + age.int * mexp.pexp.m$age.x
                                           + status.int * mexp.pexp.m$status.x
                                           + stat.age * mexp.pexp.m$status.x * mexp.pexp.m$age.x))
  for(i in 1:length(ode_res)){
    if(ode_res[i]>0.999){ode_res[i]=0.999}
    if(ode_res[i]<0.001){ode_res[i]=0.001}
  }
  # Removing the initial condition
  nll = -1*sum(dbinom(mexp.pexp.m$met_inf, size = 1, prob = ode_res, log=TRUE))
  return(nll)
  
}

twin <- c(-20:20)/2
mets.explore <- matrix(0,41^4,5)
index <- 1
for(w in twin){
  for(x in twin){
    for(y in twin){
      for(z in twin){
        mets.explore[index,1]<-(w/50)^3
        mets.explore[index,2]<-(x/50)^3
        mets.explore[index,3]<-(y/50)^3
        mets.explore[index,4]<-(z/50)^3
        mets.explore[index,5]<-sibr_nll((w/50)^3,(x/50)^3,(y/50)^3,(z/50)^3)
        index <- index + 1
      }
    }
  }
}


min.ind <- which.min(mets.explore[,5])


start_params = c(mets.explore[min.ind,1],mets.explore[min.ind,2],mets.explore[min.ind,3], mets.explore[min.ind,4])

# Get the mle estimates for my parameters
fit_all_m = mle2(sibr_nll, 
               start=list(infect=start_params[1],
                          age.int=start_params[2],
                          status.int=start_params[3],
                          stat.age=start_params[4]),
               method="Nelder-Mead")

summary(fit_all_m)



####Visualize
####first thing to do is to make a sheet with mean, min, max of infectivity
psus <- matrix(0,23,6)
for(i in 1:23){
  psus[i,1] <- i
  psus[i,2] <- coef(summary(fit_all_p))[1,1] + coef(summary(fit_all_p))[2,1]*i
  psus[i,3] <- coef(summary(fit_all_p))[1,1] - coef(summary(fit_all_p))[1,2] + 
    (coef(summary(fit_all_p))[2,1] - coef(summary(fit_all_p))[2,2])*i
  psus[i,4] <- coef(summary(fit_all_p))[1,1] + coef(summary(fit_all_p))[1,2] + 
    (coef(summary(fit_all_p))[2,1] + coef(summary(fit_all_p))[2,2])*i
  psus[i,5] <- "Uninfected"
  psus[i,6] <- "Susceptibility to Bacteria"
}

msus <- matrix(0,23,6)
for(i in 1:23){
  msus[i,1] <- i
  msus[i,2] <- coef(summary(fit_all_m))[1,1] + coef(summary(fit_all_m))[2,1]*i
  msus[i,3] <- coef(summary(fit_all_m))[1,1] - coef(summary(fit_all_m))[1,2] + 
    (coef(summary(fit_all_m))[2,1] - coef(summary(fit_all_m))[2,2])*i
  msus[i,4] <- coef(summary(fit_all_m))[1,1] + coef(summary(fit_all_m))[1,2] + 
    (coef(summary(fit_all_m))[2,1] + coef(summary(fit_all_m))[2,2])*i
  msus[i,5] <- "Uninfected"
  msus[i,6] <- "Susceptibility to Fungus"
}

prempsus <- matrix(0,23,6)
for(i in 1:23){
  prempsus[i,1] <- i
  prempsus[i,2] <- coef(summary(fit_all_p))[1,1] + coef(summary(fit_all_p))[3,1] + 
    (coef(summary(fit_all_p))[2,1] + coef(summary(fit_all_p))[4,1]) * i
  prempsus[i,3] <- coef(summary(fit_all_p))[1,1] - coef(summary(fit_all_p))[1,2] +
    coef(summary(fit_all_p))[3,1] - coef(summary(fit_all_p))[3,2] +
    (coef(summary(fit_all_p))[2,1] - coef(summary(fit_all_p))[2,2] +
       coef(summary(fit_all_p))[4,1] - coef(summary(fit_all_p))[4,2]) * i
  prempsus[i,4] <- coef(summary(fit_all_p))[1,1] + coef(summary(fit_all_p))[1,2] +
  coef(summary(fit_all_p))[3,1] + coef(summary(fit_all_p))[3,2] +
  (coef(summary(fit_all_p))[2,1] + coef(summary(fit_all_p))[2,2] +
      coef(summary(fit_all_p))[4,1] + coef(summary(fit_all_p))[4,2]) * i
  prempsus[i,5] <- "Prior Fungal Infection"
  prempsus[i,6] <- "Susceptibility to Bacteria"
}

prepmsus <- matrix(0,23,6)
for(i in 1:23){
  prepmsus[i,1] <- i
  prepmsus[i,2] <- coef(summary(fit_all_m))[1,1] + coef(summary(fit_all_m))[3,1] + 
    (coef(summary(fit_all_m))[2,1] + coef(summary(fit_all_m))[4,1]) * i
  prepmsus[i,3] <- coef(summary(fit_all_m))[1,1] - coef(summary(fit_all_m))[1,2] +
  coef(summary(fit_all_m))[3,1] - coef(summary(fit_all_m))[3,2] +
  (coef(summary(fit_all_m))[2,1] - coef(summary(fit_all_m))[2,2] +
      coef(summary(fit_all_m))[4,1] - coef(summary(fit_all_m))[4,2]) * i
  prepmsus[i,4] <- coef(summary(fit_all_m))[1,1] + coef(summary(fit_all_m))[1,2] +
  coef(summary(fit_all_m))[3,1] + coef(summary(fit_all_m))[3,2] +
  (coef(summary(fit_all_m))[2,1] + coef(summary(fit_all_m))[2,2] +
      coef(summary(fit_all_m))[4,1] + coef(summary(fit_all_m))[4,2]) * i
  prepmsus[i,5] <- "Prior Bacterial Infection"
  prepmsus[i,6] <- "Susceptibility to Fungus"
}

sus_vis <- rbind(psus,msus,prempsus,prepmsus)
sus_vis <- as.data.frame(sus_vis)
colnames(sus_vis) <- c("days","mean","mean_lower","mean_upper","status","Experiment")

sus_vis$days <- as.numeric(sus_vis$days)
sus_vis$mean <- as.numeric(sus_vis$mean)
sus_vis$mean_lower <- as.numeric(sus_vis$mean_lower)
sus_vis$mean_upper <- as.numeric(sus_vis$mean_upper)
sus_vis$status <- as.factor(sus_vis$status)
sus_vis$Experiment <- as.factor(sus_vis$Experiment)


tiff(filename = "susceptibility.tiff", width = 3000, height = 1200, pointsize = 12,res=300)

ggplot(sus_vis,aes(x=days,y=mean,color=factor(status),fill=factor(status))) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(colour = "black", size=1)) + 
  geom_line() +
  geom_ribbon(aes(ymin = mean_lower, ymax = mean_upper),alpha = 0.5) +
  geom_line(lwd=2) +
  #coord_cartesian(ylim = c(0.00002,0.00023)) +
  scale_colour_manual(values = c("#009E73", "#F0E442","black")) +
  scale_fill_manual(values = c("#009E73", "#F0E442","black")) +
  xlab("Age at Infection (Days)") + ylab("Per Spore Infection Probability") + #ggtitle("Metsch fitness inc. both infectivity and spore yield") +
  theme(axis.title.x = element_text(face="bold", size=16),axis.text.x  = element_text(size=16)) + 
  theme(axis.title.y = element_text(face="bold", size=16),axis.text.y  = element_text(size=16)) +
  theme(legend.title=element_blank()) +
  facet_wrap(~ Experiment)#, scales = "free")

dev.off()
