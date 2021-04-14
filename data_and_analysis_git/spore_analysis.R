

#### The following code will Screen the results of our age_coinfection experiment for whether our various age/coinfection treatments had any effect on pathogen or host fitness.

#### This is a rough first visualization, using linear glm's to estimate impact of treatments. Further steps involve looking for non-linearities in the impact of age on host and pathogen fitness.  


require(ggplot2)
require(plyr)
require(dplyr)
require(tidyr)
require(lubridate)
require(MASS)
require(fitdistrplus)

spores <- read.csv('spore_counts.csv')

spores$coinf <- as.factor(spores$coinf)
spores$coinf_treat <- as.factor(spores$coinf_treat)
spores$status <- as.factor(spores$status)
spores$treatment <- as.factor(spores$treatment)

#Look at second arriving individuals outcompeting first arrivers.

#get only past exposed
past_alone <- spores %>%
  filter(past_treat == 1 & met_treat == 0) 

#get past first with metsch success
past_first <- spores %>%
  filter(treatment == 4 | treatment == 8 | treatment == 12) %>%
  filter(met_inf == 1)

super_metsch <- rbind(past_alone,past_first)

summary(lm(data = super_metsch, past_inf ~ met_treat))

#get only metsch exposed
mets_alone <- spores %>%
  filter(past_treat == 0 & met_treat == 1) 

#get mets first with past success
mets_first <- spores %>%
  filter(treatment == 5 | treatment == 9 | treatment == 13) %>%
  filter(past_inf == 1)

super_past <- rbind(mets_alone,mets_first)

summary(lm(data = super_past, met_inf ~ past_treat))

#####################


response_spores_past <- spores %>%
  filter(missing == 0) %>%
  filter(treat_fail == 0) %>%
  filter(past_treat == 1) %>%
  filter(cont == 0)

response_spores_met <- spores %>%
  filter(missing == 0) %>%
  filter(treat_fail == 0) %>%
  filter(met_treat == 1) %>%
  filter(cont == 0)

response_infection <- spores %>%
  filter(missing == 0) %>%
  filter(cont == 0)

proportion_infection <- response_infection %>%
  group_by(age,coinf_treat,status) %>%
  summarise(n = n()) %>%
  mutate(freq = n / sum(n))

past_success  <- spores %>%
  filter(missing == 0) %>%
  filter(past_treat == 1) %>%
  filter(cont == 0)

met_success  <- spores %>%
  filter(missing == 0) %>%
  filter(met_treat == 1) %>%
  filter(cont == 0)


#### First, we will look at how age at first infection and whether a host is singly infected, coinfected with pasteuria exposure first, or coinfected with metschnikowia exposure first alters the total number of pasteuria spores released from successfully infected hosts.

#### A simple GLM shows that age has a significant effect on spore yield, and has a significant interaction with coinfection when pasteuria arrives first.

#### Visualizing this data, we see that age at first infection has a negative relationship with pasteuria spore yield, except in coinfected hosts when pasteuria arrives first, where age has a positive relationship with pasteuria spore yield.

response_spores_past$age <- (response_spores_past$age * 8) -4

response_spores_past$coinf_treat <- factor(response_spores_past$coinf_treat, levels = c("1","2","3"),
                                           labels = c("Single Inf.","Past First","Met First"))

##First, what distributions does this data fall under
##Visualize fit of data to different distribution
fit.norm<-fitdist(response_spores_past$average_past,"norm")  
fit.gamma<-fitdist(response_spores_past$average_past,"gamma",lower=c(0,0))
fit.weibull<-fitdist(response_spores_past$average_past,"weibull",lower=c(0,0))
fit.lnorm<-fitdist(response_spores_past$average_past,"lnorm")
fit.pois<-fitdist(round(response_spores_past$average_past),"pois")
fit.nbin<-fitdist(round(response_spores_past$average_past),"nbinom")

#OK, can't make a correlation matrix, but that's ok

summary(fit.norm)   # AIC = 1044.638
summary(fit.gamma)   # AIC = 864.4365
summary(fit.weibull)   # AIC = 869.8534
summary(fit.lnorm)   # AIC = 882.5094
summary(fit.pois)   # AIC = 11015.45
summary(fit.nbin)   # AIC = 867.3897

#ok, so could feasibly do gamma or nbinom


##want to test linear, exp, and saturating functions.
##There's a possibility that each infection treatment will have aa different function
pspore.poly <- glm((average_past*1000) ~ age + coinf_treat + coinf_treat:age, 
                   data = response_spores_past, family=Gamma(link = "log"),
                   maxit = 100)

summary(pspore.poly,dispersion=1) ##need dispersion =1 for gamma family



levels(response_spores_past$coinf_treat) <- c("Bacteria Only","Bacteria First","Fungi First")


#Ok, let's try this again, but this time simulating the trend a bunch of times

bac.spore.bac.only.sim <- matrix(0,26,1000)
for(i in 1:1000){
  bac.spore.bac.only.sim[1,i] <- rnorm(1,mean = coef(summary(pspore.poly))[1,1],sd=coef(summary(pspore.poly))[1,2])
  bac.spore.bac.only.sim[2,i] <- rnorm(1,mean = coef(summary(pspore.poly))[2,1],sd=coef(summary(pspore.poly))[2,2])
  for(j in 0:23){
    bac.spore.bac.only.sim[j+3,i] <- bac.spore.bac.only.sim[1,i] + bac.spore.bac.only.sim[2,i] * j
  }
  }

bac.spore.bac.only <- matrix(0,24,5)
for(i in 0:23){
  bac.spore.bac.only[i+1,1] <- i
  bac.spore.bac.only[i+1,2] <- mean(bac.spore.bac.only.sim[i+3,])
  bac.spore.bac.only[i+1,3] <- mean(bac.spore.bac.only.sim[i+3,]) + sd(bac.spore.bac.only.sim[i+3,])
  bac.spore.bac.only[i+1,4] <- mean(bac.spore.bac.only.sim[i+3,]) - sd(bac.spore.bac.only.sim[i+3,])
  bac.spore.bac.only[i+1,5] <- "Bacteria Only"
}

bac.spore.bac.first.sim <- matrix(0,26,1000)
for(i in 1:1000){
  bac.spore.bac.first.sim[1,i] <- rnorm(1,mean = coef(summary(pspore.poly))[1,1] + coef(summary(pspore.poly))[3,1],sd=coef(summary(pspore.poly))[1,2] + coef(summary(pspore.poly))[3,2])
  bac.spore.bac.first.sim[2,i] <- rnorm(1,mean = coef(summary(pspore.poly))[2,1] + coef(summary(pspore.poly))[5,1],sd=coef(summary(pspore.poly))[2,2] + coef(summary(pspore.poly))[5,2])
  for(j in 0:23){
    bac.spore.bac.first.sim[j+3,i] <- bac.spore.bac.first.sim[1,i] + bac.spore.bac.first.sim[2,i] * j
  }
}

bac.spore.bac.first <- matrix(0,24,5)
for(i in 0:23){
  bac.spore.bac.first[i+1,1] <- i
  bac.spore.bac.first[i+1,2] <- mean(bac.spore.bac.first.sim[i+3,])
  bac.spore.bac.first[i+1,3] <- mean(bac.spore.bac.first.sim[i+3,]) + sd(bac.spore.bac.first.sim[i+3,])
  bac.spore.bac.first[i+1,4] <- mean(bac.spore.bac.first.sim[i+3,]) - sd(bac.spore.bac.first.sim[i+3,])
  bac.spore.bac.first[i+1,5] <- "Bacteria First"
}

bac.spore.fun.first.sim <- matrix(0,26,1000)
for(i in 1:1000){
  bac.spore.fun.first.sim[1,i] <- rnorm(1,mean = coef(summary(pspore.poly))[1,1],sd=coef(summary(pspore.poly))[1,2])
  bac.spore.fun.first.sim[2,i] <- rnorm(1,mean = coef(summary(pspore.poly))[2,1] + coef(summary(pspore.poly))[6,1],sd=coef(summary(pspore.poly))[2,2] + coef(summary(pspore.poly))[6,2])
  for(j in 0:23){
    bac.spore.fun.first.sim[j+3,i] <- bac.spore.fun.first.sim[1,i] + bac.spore.fun.first.sim[2,i] * j
  }
}

bac.spore.fun.first <- matrix(0,24,5)
for(i in 0:23){
  bac.spore.fun.first[i+1,1] <- i
  bac.spore.fun.first[i+1,2] <- mean(bac.spore.fun.first.sim[i+3,])
  bac.spore.fun.first[i+1,3] <- mean(bac.spore.fun.first.sim[i+3,]) + sd(bac.spore.fun.first.sim[i+3,])
  bac.spore.fun.first[i+1,4] <- mean(bac.spore.fun.first.sim[i+3,]) - sd(bac.spore.fun.first.sim[i+3,])
  bac.spore.fun.first[i+1,5] <- "Fungi First"
}

bac.spore <- rbind(bac.spore.bac.only,bac.spore.bac.first,bac.spore.fun.first)
bac.spore <- as.data.frame(bac.spore)
colnames(bac.spore) <- c("age","mean","max","min","coinf_treat")
bac.spore$age <- as.numeric(bac.spore$age)
bac.spore$mean <- as.numeric(bac.spore$mean)
bac.spore$min <- as.numeric(bac.spore$min)
bac.spore$max <- as.numeric(bac.spore$max)

levels(bac.spore$coinf_treat) <- c("Bacteria Only","Bacteria First","Fungi First")


tiff(filename = "past_spore_load.tiff", width = 2000, height = 1500, pointsize = 12,res=300)

ggplot() +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(colour = "black", size=1)) + 
  #stat_smooth(method="lm",formula = y ~ x,se=TRUE) +
  #stat_smooth(method="lm",formula = y ~ x,se=FALSE) +
  coord_cartesian(ylim = c(5,13),xlim = c(3,21)) +
  geom_ribbon(data = bac.spore,aes(ymin=min,ymax=max,x=age,color=factor(coinf_treat),fill=factor(coinf_treat)),alpha = 0.25) + 
  geom_line(data = bac.spore,aes(y=mean,x=age,color=factor(coinf_treat)),lwd=1) +
  geom_jitter(data = response_spores_past,aes(y=(log(average_past*1000)),x=age,color=factor(coinf_treat),fill=factor(coinf_treat),shape=factor(coinf_treat)),width=1, size = 2) +
  scale_colour_manual(values = c("#009E73","#56B4E9","#F0E442")) +
  scale_fill_manual(values = c("#009E73","#56B4E9","#F0E442")) +
  xlab("Age at Infection (Days)") + ylab("Log Bacteria Spores per Host") + #ggtitle("Impact of age and coinfection on Pasteuria spore load") +
  theme(axis.title.x = element_text(face="bold", size=16),axis.text.x  = element_text(size=16)) + 
  theme(axis.title.y = element_text(face="bold", size=16),axis.text.y  = element_text(size=16)) +
  theme(legend.title=element_blank())

dev.off()
#scale_colour_manual(values = c("#000000", "#E69F00","#009E73","#56B4E9","#F0E442")) +
  


#### Now we will repeat the process for Metschnikowia spore yield 

#### A simple GLM shows neither age nor coinfection alter metschnikowia spore yield.

#### Visualizing this data, We see that indeed relationships between age and spore yield largely overlap, and don't have strong trends. 

response_spores_met$age <- (response_spores_met$age * 8) -4

response_spores_met$coinf_treat <- factor(response_spores_met$coinf_treat, levels = c("1","2","3"),
                                          labels = c("Single Inf.","Past First","Met First"))

##First, what distributions does this data fall under
##Visualize fit of data to different distribution
fit.norm<-fitdist(response_spores_met$average_met,"norm")  
fit.gamma<-fitdist(response_spores_met$average_met,"gamma",lower=c(0,0))
fit.weibull<-fitdist(response_spores_met$average_met,"weibull",lower=c(0,0))
fit.lnorm<-fitdist(response_spores_met$average_met,"lnorm")
fit.pois<-fitdist(round(response_spores_met$average_met),"pois")
fit.nbin<-fitdist(round(response_spores_met$average_met),"nbinom")

#OK, can't make a correlation matrix, but that's ok

summary(fit.norm)  ##AIC = 1040
summary(fit.gamma)  ##AIC = 1053
summary(fit.weibull)  ##AIC = 1043
summary(fit.lnorm)  ##AIC = 1105
summary(fit.pois)  ##AIC = 3112
summary(fit.nbin)  ##AIC = 1050


##want to test linear, exp, and saturating functions.
##There's a possibility that each infection treatment will have aa different function
mspore.poly <- lm(average_met ~ age + I(age^2) + coinf_treat + coinf_treat:age, data = response_spores_met)
summary(mspore.poly) ##need dispersion =1 for gamma family

#Redo without the polynomial
mspore.poly <- lm(average_met ~ age + coinf_treat + coinf_treat:age, data = response_spores_met)
summary(mspore.poly) #Still no significance

#To get actual parameter, run without age or treatment
mspore.poly <- lm(average_met ~ 1, data = response_spores_met)
summary(mspore.poly) #parameter is 63.93


##Now what redo analysis taking away age and treatment

met_spore_fig <- response_spores_met
levels(met_spore_fig$coinf_treat) <- c("Fungi Only","Bacteria First","Fungi First")

#Ok, let's try this again, but this time simulating the trend a bunch of times

fun.spore.fun.only.sim <- matrix(0,26,1000)
for(i in 1:1000){
  fun.spore.fun.only.sim[1,i] <- log(rnorm(1,mean = coef(summary(mspore.poly))[1,1],sd=2*coef(summary(mspore.poly))[1,2])*1000)
  fun.spore.fun.only.sim[2,i] <- 0
  for(j in 0:23){
    fun.spore.fun.only.sim[j+3,i] <- fun.spore.fun.only.sim[1,i] + fun.spore.fun.only.sim[2,i] * j
  }
}

fun.spore.fun.only <- matrix(0,24,5)
for(i in 0:23){
  fun.spore.fun.only[i+1,1] <- i
  fun.spore.fun.only[i+1,2] <- mean(fun.spore.fun.only.sim[i+3,])
  fun.spore.fun.only[i+1,3] <- mean(fun.spore.fun.only.sim[i+3,]) + sd(fun.spore.fun.only.sim[i+3,])
  fun.spore.fun.only[i+1,4] <- mean(fun.spore.fun.only.sim[i+3,]) - sd(fun.spore.fun.only.sim[i+3,])
  fun.spore.fun.only[i+1,5] <- "Fungi Only"
}


fun.spore <- fun.spore.fun.only
fun.spore <- as.data.frame(fun.spore)
colnames(fun.spore) <- c("age","mean","max","min","coinf_treat")
fun.spore$age <- as.numeric(fun.spore$age)
fun.spore$mean <- as.numeric(fun.spore$mean)
fun.spore$min <- as.numeric(fun.spore$min)
fun.spore$max <- as.numeric(fun.spore$max)


tiff(filename = "mets_spore_load.tiff", width = 2000, height = 1500, pointsize = 12,res=300)

ggplot() +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(colour = "black", size=1)) + 
  #stat_smooth(method="lm",se=TRUE,color="black") +
  coord_cartesian(xlim = c(3,21)) +
  geom_ribbon(data = fun.spore,aes(ymin=min,ymax=max,x=age),alpha = 0.5) + 
  geom_line(data = fun.spore,aes(y=mean,x=age),lwd=1) +
  geom_jitter(data = met_spore_fig,aes(y=(log(average_met*1000)),x=age,color=factor(coinf_treat),fill=factor(coinf_treat),shape=factor(coinf_treat)),width=1, size = 2) +
  scale_colour_manual(values = c("#E69F00", "#009E73", "#F0E442")) +
  scale_fill_manual(values = c("#E69F00", "#009E73", "#F0E442")) +
  xlab("Age at Infection (Days)") + ylab("Log Fungi Spores per Host") + #ggtitle("Impact of age and coinfection on Metschnikowia spore load") +
  theme(axis.title.x = element_text(face="bold", size=16),axis.text.x  = element_text(size=16)) + 
  theme(axis.title.y = element_text(face="bold", size=16),axis.text.y  = element_text(size=16)) +
  theme(legend.title=element_blank())

dev.off()


#### Spore yield is, of course, only one aspect of pathogen fitness. We will also examine the likelihood that pasteuria successfully infects its hosts, given age at first exposure and coexposure by metschnikowia

#### Both a binomial glm and data visualization show that increasing age at infection lowers likelihood of pasteuria infecting a host, and that prior exposure to metschnikowia further lowers this likelihood. 


past_success$age <- (past_success$age * 8) -4

past_success$coinf_treat <- factor(past_success$coinf_treat, levels = c("1","2","3"),
                                   labels = c("Single Inf.","Past First","Met First"))

#The data is binomial so just go with that

past_prevalence_glm <- glm(past_inf ~ age + I(age^2) + coinf_treat + coinf_treat:age, family=binomial(link="logit"), data=past_success)

summary(past_prevalence_glm)


ggplot(past_success,aes(y=past_inf,x=age,color=factor(coinf_treat),fill=factor(coinf_treat))) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(colour = "black", size=1)) + 
  geom_jitter(width=1.5,height=.05) +
  stat_smooth(method="lm",se=TRUE) +
  #coord_cartesian(ylim = c(0,150)) +
  scale_colour_manual(values = c("blue", "black", "green")) +
  scale_fill_manual(values = c("blue", "black", "green")) +
  xlab("Age at Infection (Days)") + ylab("Proportion infected by past") + ggtitle("Impact of age and coinfection on Pasteuria prevalence") +
  theme(axis.title.x = element_text(face="bold", size=16),axis.text.x  = element_text(size=16)) + 
  theme(axis.title.y = element_text(face="bold", size=16),axis.text.y  = element_text(size=16))



#### Examining the impact of age and pasteuria exposure on the likelihood that metschnikowia can successfully infect a host, we see that when hosts are previously exposed to pasteuria, age at first infection has a negative relationship with likelihood of metschnikowia infection, but that other wise age doesn't have a large impact on metschnikowia infection success. 


met_success$coinf_treat <- factor(met_success$coinf_treat, levels = c("1","2","3"),
                                  labels = c("Single Inf.","Past First","Met First"))

met_success$age <- (met_success$age * 8) -4

met_prevalence_glm <- glm(met_inf ~ age + I(age^2) + coinf_treat + coinf_treat:age, family=binomial(link="logit"), data=met_success)

summary(met_prevalence_glm)


ggplot(met_success,aes(y=met_inf,x=age,color=factor(coinf_treat),fill=factor(coinf_treat))) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(colour = "black", size=1)) + 
  geom_jitter(width=1.5,height=.05) +
  stat_smooth(method="lm",se=TRUE) +
  #coord_cartesian(ylim = c(0,150)) +
  scale_colour_manual(values = c("blue", "black", "green")) +
  scale_fill_manual(values = c("blue", "black", "green")) +
  xlab("Age at Infection (Days)") + ylab("Proportion Infected by metsch") + ggtitle("Impact of age and coinfection on Metschnikowia prevalence") +
  theme(axis.title.x = element_text(face="bold", size=16),axis.text.x  = element_text(size=16)) + 
  theme(axis.title.y = element_text(face="bold", size=16),axis.text.y  = element_text(size=16))




#### Finally, let's look at the effect of age and metschnikowia exposure on total pasteuria fitness, including both likelihood of infection and average spore yield. The response metric here will be total spore yield, but this includes those individuals that were not successfully infected. 
#### Thus, we can interpret this as the predicted number of infectious spores to come from a host, given that it was exposed to pasteuria. 

#### Now, pretty much everything has a significant interactive effect on pasteuria fitness, and the patterns we saw in the first figure become more prominent.




##First, what distributions does this data fall under
##Visualize fit of data to different distribution
fit.norm.1<-fitdist(past_success$average_past,"norm")  
fit.gamma.1<-fitdist(past_success$average_past,"gamma",lower=c(0,0))
fit.weibull.1<-fitdist(past_success$average_past,"weibull",lower=c(0,0))
fit.pois.1<-fitdist(round(past_success$average_past),"pois")
fit.nbin.1<-fitdist(round(past_success$average_past),"nbinom")


summary(fit.norm.1)## nbinom must be bounded by 1, so normal actually wins
summary(fit.pois.1)
summary(fit.nbin.1) #(Nbinom definitely wins, but glm doesn't do negative binom)


past_total_glm <- glm.nb(round(average_past) ~ age + I(age^2) + coinf_treat + coinf_treat:age, data=past_success)

summary(past_total_glm)


past_spore_fig <- past_success
levels(past_spore_fig$coinf_treat) <- c("Bacteria Only","Bacteria First","Fungi First")

#Ok, let's try this again, but this time simulating the trend a bunch of times

bac.spore.bac.only.sim <- matrix(0,26,1000)
for(i in 1:1000){
  bac.spore.bac.only.sim[1,i] <- rnorm(1,mean = coef(summary(past_total_glm))[1,1],sd=coef(summary(past_total_glm))[1,2])
  bac.spore.bac.only.sim[2,i] <- 0
  for(j in 0:23){
    bac.spore.bac.only.sim[j+3,i] <- bac.spore.bac.only.sim[1,i] + bac.spore.bac.only.sim[2,i] * j
  }
}

bac.spore.bac.only <- matrix(0,24,5)
for(i in 0:23){
  bac.spore.bac.only[i+1,1] <- i
  bac.spore.bac.only[i+1,2] <- mean(bac.spore.bac.only.sim[i+3,])
  bac.spore.bac.only[i+1,3] <- mean(bac.spore.bac.only.sim[i+3,]) + sd(bac.spore.bac.only.sim[i+3,])
  bac.spore.bac.only[i+1,4] <- mean(bac.spore.bac.only.sim[i+3,]) - sd(bac.spore.bac.only.sim[i+3,])
  bac.spore.bac.only[i+1,5] <- "Bacteria Only"
}

bac.spore.bac.first.sim <- matrix(0,26,1000)
for(i in 1:1000){
  bac.spore.bac.first.sim[1,i] <- rnorm(1,mean = coef(summary(past_total_glm))[1,1] + coef(summary(past_total_glm))[4,1],sd=coef(summary(past_total_glm))[1,2] + coef(summary(past_total_glm))[4,2])
  bac.spore.bac.first.sim[2,i] <- 0
  for(j in 0:23){
    bac.spore.bac.first.sim[j+3,i] <- bac.spore.bac.first.sim[1,i] + bac.spore.bac.first.sim[2,i] * j
  }
}

bac.spore.bac.first <- matrix(0,24,5)
for(i in 0:23){
  bac.spore.bac.first[i+1,1] <- i
  bac.spore.bac.first[i+1,2] <- mean(bac.spore.bac.first.sim[i+3,])
  bac.spore.bac.first[i+1,3] <- mean(bac.spore.bac.first.sim[i+3,]) + sd(bac.spore.bac.first.sim[i+3,])
  bac.spore.bac.first[i+1,4] <- mean(bac.spore.bac.first.sim[i+3,]) - sd(bac.spore.bac.first.sim[i+3,])
  bac.spore.bac.first[i+1,5] <- "Bacteria First"
}

bac.spore.fun.first.sim <- matrix(0,26,1000)
for(i in 1:1000){
  bac.spore.fun.first.sim[1,i] <- rnorm(1,mean = coef(summary(past_total_glm))[1,1],sd=coef(summary(past_total_glm))[1,2])
  bac.spore.fun.first.sim[2,i] <- rnorm(1,mean = coef(summary(past_total_glm))[7,1],sd=coef(summary(past_total_glm))[7,2])
  for(j in 0:23){
    bac.spore.fun.first.sim[j+3,i] <- bac.spore.fun.first.sim[1,i] + bac.spore.fun.first.sim[2,i] * j
  }
}

bac.spore.fun.first <- matrix(0,24,5)
for(i in 0:23){
  bac.spore.fun.first[i+1,1] <- i
  bac.spore.fun.first[i+1,2] <- mean(bac.spore.fun.first.sim[i+3,])
  bac.spore.fun.first[i+1,3] <- mean(bac.spore.fun.first.sim[i+3,]) + sd(bac.spore.fun.first.sim[i+3,])
  bac.spore.fun.first[i+1,4] <- mean(bac.spore.fun.first.sim[i+3,]) - sd(bac.spore.fun.first.sim[i+3,])
  bac.spore.fun.first[i+1,5] <- "Fungi First"
}

bac.spore <- rbind(bac.spore.bac.only,bac.spore.bac.first,bac.spore.fun.first)
bac.spore <- as.data.frame(bac.spore)
colnames(bac.spore) <- c("age","mean","max","min","coinf_treat")
bac.spore$age <- as.numeric(bac.spore$age)
bac.spore$mean <- as.numeric(bac.spore$mean)
bac.spore$min <- as.numeric(bac.spore$min)
bac.spore$max <- as.numeric(bac.spore$max)

levels(bac.spore$coinf_treat) <- c("Bacteria Only","Bacteria First","Fungi First")


tiff(filename = "past_all_spores.tiff", width = 2000, height = 1500, pointsize = 12,res=300)

ggplot() +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(colour = "black", size=1)) + 
  #stat_smooth(method="lm",formula = y ~ x,se=TRUE) +
  #stat_smooth(method="lm",formula = y ~ x,se=FALSE) +
  coord_cartesian(ylim = c(-0.5,14),xlim = c(3,21)) +
  geom_ribbon(data = bac.spore,aes(ymin=min,ymax=max,x=age,color=factor(coinf_treat),fill=factor(coinf_treat)),alpha = 0.25) + 
  geom_line(data = bac.spore,aes(y=mean,x=age,color=factor(coinf_treat)),lwd=1) +
  geom_jitter(data = past_spore_fig,aes(y=(log(average_past*1000+1)),x=age,color=factor(coinf_treat),fill=factor(coinf_treat),shape=factor(coinf_treat)),width=1,height=.5, size = 2) +
  scale_colour_manual(values = c("#009E73","#56B4E9","#F0E442")) +
  scale_fill_manual(values = c("#009E73","#56B4E9","#F0E442")) +
  xlab("Age at Infection (Days)") + ylab("Log Bacteria Spores per Host") + #ggtitle("Impact of age and coinfection on Pasteuria spore load") +
  theme(axis.title.x = element_text(face="bold", size=16),axis.text.x  = element_text(size=16)) + 
  theme(axis.title.y = element_text(face="bold", size=16),axis.text.y  = element_text(size=16)) +
  theme(legend.title=element_blank())

dev.off()



#### Repeating the process for Metschnikowia, however, metschnikowia still has a pretty robust fitness, regardless of treatment. 




fit.norm.2<-fitdist(met_success$average_met,"norm")  
fit.gamma.2<-fitdist(met_success$average_met,"gamma",lower=c(0,0))
fit.weibull.2<-fitdist(met_success$average_met,"weibull",lower=c(0,0))
fit.lnorm.2<-fitdist(met_success$average_met,"lnorm")
fit.pois.2<-fitdist(round(met_success$average_met),"pois")
fit.nbinom.2<-fitdist(round(met_success$average_met),"nbinom")


#OK, can't make a correlation matrix, but that's ok

summary(fit.norm.2)## nbinom must be bounded by 1
summary(fit.gamma.2)
summary(fit.weibull.2)
summary(fit.lnorm.2)
summary(fit.pois.2)
summary(fit.nbinom.2)



met_total_glm <- glm.nb(average_met ~ age + I(age^2) + coinf_treat + coinf_treat:age, data=met_success)

summary(met_total_glm)

met_spore_fig <- met_success
levels(met_spore_fig$coinf_treat) <- c("Fungi Only","Bacteria First","Fungi First")


#Ok, let's try this again, but this time simulating the trend a bunch of times

fun.spore.fun.only.sim <- matrix(0,27,1000)
for(i in 1:1000){
  fun.spore.fun.only.sim[1,i] <- rnorm(1,mean = coef(summary(met_total_glm))[1,1],sd=coef(summary(met_total_glm))[1,2])
  fun.spore.fun.only.sim[2,i] <- rnorm(1,mean = coef(summary(met_total_glm))[2,1],sd=coef(summary(met_total_glm))[2,2])
  fun.spore.fun.only.sim[3,i] <- rnorm(1,mean = coef(summary(met_total_glm))[3,1],sd=coef(summary(met_total_glm))[3,2])
  for(j in 0:23){
    fun.spore.fun.only.sim[j+4,i] <- fun.spore.fun.only.sim[1,i] + fun.spore.fun.only.sim[2,i] * j + fun.spore.fun.only.sim[3,i] * j^2
  }
}

fun.spore.fun.only <- matrix(0,24,5)
for(i in 0:23){
  fun.spore.fun.only[i+1,1] <- i
  fun.spore.fun.only[i+1,2] <- mean(fun.spore.fun.only.sim[i+4,])
  fun.spore.fun.only[i+1,3] <- mean(fun.spore.fun.only.sim[i+4,]) + sd(fun.spore.fun.only.sim[i+4,])
  fun.spore.fun.only[i+1,4] <- mean(fun.spore.fun.only.sim[i+4,]) - sd(fun.spore.fun.only.sim[i+4,])
  fun.spore.fun.only[i+1,5] <- "Fungi Only"
}

fun.spore.bac.first.sim <- matrix(0,27,1000)
for(i in 1:1000){
  fun.spore.bac.first.sim[1,i] <- rnorm(1,mean = coef(summary(met_total_glm))[1,1] + coef(summary(met_total_glm))[4,1],sd=coef(summary(met_total_glm))[1,2] + coef(summary(met_total_glm))[4,2])
  fun.spore.bac.first.sim[2,i] <- rnorm(1,mean = coef(summary(met_total_glm))[2,1] + coef(summary(met_total_glm))[6,1],sd=coef(summary(met_total_glm))[2,2] + coef(summary(met_total_glm))[6,2])
  fun.spore.bac.first.sim[3,i] <- rnorm(1,mean = coef(summary(met_total_glm))[3,1],sd=coef(summary(met_total_glm))[3,2])
  for(j in 0:23){
    fun.spore.bac.first.sim[j+4,i] <- fun.spore.bac.first.sim[1,i] + fun.spore.bac.first.sim[2,i] * j + fun.spore.bac.first.sim[3,i] * j^2
  }
}

fun.spore.bac.first <- matrix(0,24,5)
for(i in 0:23){
  fun.spore.bac.first[i+1,1] <- i
  fun.spore.bac.first[i+1,2] <- mean(fun.spore.bac.first.sim[i+4,])
  fun.spore.bac.first[i+1,3] <- mean(fun.spore.bac.first.sim[i+4,]) + sd(fun.spore.bac.first.sim[i+4,])
  fun.spore.bac.first[i+1,4] <- mean(fun.spore.bac.first.sim[i+4,]) - sd(fun.spore.bac.first.sim[i+4,])
  fun.spore.bac.first[i+1,5] <- "Bacteria First"
}

fun.spore.fun.first.sim <- matrix(0,27,1000)
for(i in 1:1000){
  fun.spore.fun.first.sim[1,i] <- rnorm(1,mean = coef(summary(met_total_glm))[1,1],sd=coef(summary(met_total_glm))[1,2])
  fun.spore.fun.first.sim[2,i] <- rnorm(1,mean = coef(summary(met_total_glm))[2,1] + coef(summary(met_total_glm))[7,1],sd=coef(summary(met_total_glm))[2,2] + coef(summary(met_total_glm))[7,2])
  fun.spore.fun.first.sim[3,i] <- rnorm(1,mean = coef(summary(met_total_glm))[3,1],sd=coef(summary(met_total_glm))[3,2])
  for(j in 0:23){
    fun.spore.fun.first.sim[j+4,i] <- fun.spore.fun.first.sim[1,i] + fun.spore.fun.first.sim[2,i] * j + fun.spore.fun.first.sim[3,i] * j^2
  }
}

fun.spore.fun.first <- matrix(0,24,5)
for(i in 0:23){
  fun.spore.fun.first[i+1,1] <- i
  fun.spore.fun.first[i+1,2] <- mean(fun.spore.fun.first.sim[i+4,])
  fun.spore.fun.first[i+1,3] <- mean(fun.spore.fun.first.sim[i+4,]) + sd(fun.spore.fun.first.sim[i+4,])
  fun.spore.fun.first[i+1,4] <- mean(fun.spore.fun.first.sim[i+4,]) - sd(fun.spore.fun.first.sim[i+4,])
  fun.spore.fun.first[i+1,5] <- "Fungi First"
}

fun.spore <- rbind(fun.spore.fun.only,fun.spore.bac.first,fun.spore.fun.first)
fun.spore <- as.data.frame(fun.spore)
colnames(fun.spore) <- c("age","mean","max","min","coinf_treat")
fun.spore$age <- as.numeric(fun.spore$age)
fun.spore$mean <- as.numeric(fun.spore$mean)
fun.spore$min <- as.numeric(fun.spore$min)
fun.spore$max <- as.numeric(fun.spore$max)

levels(fun.spore$coinf_treat) <- c("Fungi Only","Fungi First","Bacteria First")

tiff(filename = "mets_all_spores.tiff", width = 2000, height = 1500, pointsize = 12,res=300)

ggplot() +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(colour = "black", size=1)) + 
  #stat_smooth(method="lm",formula = y ~ x,se=TRUE) +
  #stat_smooth(method="lm",formula = y ~ x,se=FALSE) +
  coord_cartesian(xlim = c(3,21)) +
  geom_ribbon(data = fun.spore,aes(ymin=min,ymax=max,x=age,color=factor(coinf_treat),fill=factor(coinf_treat)),alpha = 0.25) + 
  geom_line(data = fun.spore,aes(y=mean,x=age,color=factor(coinf_treat)),lwd=1) +
  geom_jitter(data = met_spore_fig,aes(y=(log(average_met*1000+1)),x=age,color=factor(coinf_treat),fill=factor(coinf_treat),shape=factor(coinf_treat)),width=1, height = 0.5,size = 2) +
  scale_colour_manual(values = c("#009E73",  "#F0E442","#E69F00")) +
  scale_fill_manual(values = c("#009E73",  "#F0E442","#E69F00")) +
  xlab("Age at Infection (Days)") + ylab("Log Fungi Spores per Host") + #ggtitle("Impact of age and coinfection on Pasteuria spore load") +
  theme(axis.title.x = element_text(face="bold", size=16),axis.text.x  = element_text(size=16)) + 
  theme(axis.title.y = element_text(face="bold", size=16),axis.text.y  = element_text(size=16)) +
  theme(legend.title=element_blank())

dev.off()

tiff(filename = "mets_all_spores_violin.tiff", width = 2000, height = 1500, pointsize = 12,res=300)

ggplot() +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(colour = "black", size=1)) + 
  geom_violin(data = met_spore_fig,aes(y=((average_met*1000)),x=1,color=factor(coinf_treat),fill=factor(coinf_treat))) +
  geom_boxplot(data = met_spore_fig,aes(y=((average_met*1000)),x=1,fill=factor(coinf_treat)),color="black") +
  scale_colour_manual(values = c("#009E73",  "#F0E442","#E69F00")) +
  scale_fill_manual(values = c("#009E73",  "#F0E442","#E69F00")) +
  xlab("Age at Infection (Panels)") + ylab("Fungi Spores per Host") + #ggtitle("Impact of age and coinfection on Pasteuria spore load") +
  theme(axis.title.x = element_text(face="bold", size=16),axis.text.x  = element_text(size=16)) + 
  theme(axis.title.y = element_text(face="bold", size=16),axis.text.y  = element_text(size=16)) +
  theme(legend.title=element_blank()) +
  facet_wrap(met_spore_fig$age)

dev.off()








#OK, now we need to look at some host fitness data
#metrics we can lookat are 
#total lifespan
#time til death after infection
#lifetime fecundity
#offspring per day

#First thing we want to do is compare fitness of uninfected in this experiment and in kailash's experiment

#First, to do that, we want to have a sheet where we have the fitness metrics for every successfully infected animal

#So we can first just find the metrics for every animal, and then merge into existing sheet. 

life_table <- read.csv('life_table.csv')

drops <- c("dead","missing","notes","month","day","month_adjust")

life_table <- life_table[ , !(names(life_table) %in% drops)]


#to get max lifespan, and lifespan post infection
life_table_max_life <- aggregate(life_table, by=list(life_table$treatment,life_table$replicate),FUN=max)

life_table_max_life$post_metsch_life <- life_table_max_life$host_age - life_table_max_life$day_metsch_infect
life_table_max_life$post_past_life <- life_table_max_life$host_age - life_table_max_life$day_past_infect

#lifetime_fecundity
life_table_total_babies <- aggregate(life_table, by=list(life_table$treatment,life_table$replicate),FUN=sum)
life_table_total_babies$treatment <- life_table_total_babies$Group.1
life_table_total_babies$replicate <- life_table_total_babies$Group.2

#Ok, now get average babies
life_table_average <- aggregate(life_table, by=list(life_table$treatment,life_table$replicate),FUN=mean)

#average babies post infection
life_table_post_metsch <- life_table %>%
  filter(host_age > day_infect)

life_table_post_past <- life_table %>%
  filter(host_age > day_infect)

life_table_post_metsch_average <- aggregate(life_table_post_metsch, by=list(life_table_post_metsch$treatment,life_table_post_metsch$replicate),FUN=mean)

life_table_post_past_average <- aggregate(life_table_post_past,by=list(life_table_post_past$treatment,life_table_post_past$replicate),FUN=mean)

##OK, now need to merge

useable_reps <- spores %>%
  filter(missing == 0) %>%
  filter(treat_fail == 0) %>%
  filter(cont == 0)

life_table_max_life_merge <- merge(useable_reps,life_table_max_life,by=c("treatment","replicate"))
life_table_total_babies_merge <- merge(useable_reps,life_table_total_babies,by=c("treatment","replicate"))
life_table_average_merge <- merge(useable_reps,life_table_average,by=c("treatment","replicate"))
life_table_post_metsch_merge <- merge(useable_reps,life_table_post_metsch_average,by=c("treatment","replicate"))
life_table_post_past_merge <- merge(useable_reps,life_table_post_past_average,by=c("treatment","replicate"))

#OK, now make everything relative to the healthy

life_table_max_life_healthy <- life_table_max_life_merge %>%
  filter(treatment == 1)
life_table_max_life_infected <- life_table_max_life_merge %>%
  filter(treatment != 1)

life_table_total_babies_healthy <- life_table_total_babies_merge %>%
  filter(treatment == 1)
life_table_total_babies_infected <- life_table_total_babies_merge %>%
  filter(treatment != 1)

life_table_average_healthy <- life_table_average_merge %>%
  filter(treatment == 1)
life_table_average_infected <- life_table_average_merge %>%
  filter(treatment != 1)

life_table_post_metsch_healthy <- life_table_post_metsch_merge %>%
  filter(treatment == 1)
life_table_post_metsch_infected <- life_table_post_metsch_merge %>%
  filter(treatment != 1)

life_table_post_past_healthy <- life_table_post_past_merge %>%
  filter(treatment == 1)
life_table_post_past_infected <- life_table_post_past_merge %>%
  filter(treatment != 1)

life_table_max_life_infected$babies <- life_table_max_life_infected$babies - mean(life_table_max_life_healthy$babies)
life_table_max_life_infected$host_age <- life_table_max_life_infected$host_age - mean(life_table_max_life_healthy$host_age)
#Doesn't make sense to subtract healthy stats from infection to death
life_table_total_babies_infected$babies <- life_table_total_babies_infected$babies - mean(life_table_total_babies_healthy$babies)
life_table_average_infected$babies <- life_table_average_infected$babies - mean(life_table_average_healthy$babies)
life_table_post_metsch_healthy$babies <- life_table_post_metsch_healthy$babies - mean(life_table_post_metsch_healthy$babies)
life_table_post_past_infected$babies <- life_table_post_past_infected$babies - mean(life_table_post_past_healthy$babies)

life_table_max_life_infected$coinf_treat.y <- factor(life_table_max_life_infected$coinf_treat.y, levels = c("1","2","3","4"),
                                                     labels = c("Single Past","Single Met","Past First", "Met First"))

#summing values makes this a bit different
life_table_total_babies_infected$coinf_treat.y <- life_table_max_life_infected$coinf_treat.y
life_table_total_babies_infected$day_infect <- life_table_max_life_infected$day_infect

life_table_average_infected$coinf_treat.y <- factor(life_table_average_infected$coinf_treat.y, levels = c("1","2","3","4"),
                                                    labels = c("Single Past","Single Met","Past First", "Met First"))

life_table_post_metsch_infected$coinf_treat.y <- factor(life_table_post_metsch_infected$coinf_treat.y, levels = c("1","2","3","4"),
                                                        labels = c("Single Past","Single Met","Past First", "Met First"))

life_table_post_past_infected$coinf_treat.y <- factor(life_table_post_past_infected$coinf_treat.y, levels = c("1","2","3","4"),
                                                      labels = c("Single Past","Single Met","Past First", "Met First"))

#OK, so overall we should be looking at 6 indicators of health
#Also making comparisons with the old experiment, so let's copy that sheet over

old_data <- read.csv('coinfR.csv')
old_data <- old_data %>%
  filter(Parasite == 1) %>%
  filter(Treatment == 1)



#### Ok, now let's examine host fitness
#### First let's examine lifetime offspring for each treatment. We will display lifetime offspring minus average lifetime offspring from healthy hosts.
#### In all cases, these are indeed parasites, and all infected hosts generally have less offspring than healthy hosts. 
#### We see that for all treatments, the later on in their life they are infected, the more offspring they have, except that this trend is much flatter for those only infected by metsch.
#### In fact, for pasteuria infection first, and single pasteuria inifection, may actually increase host fitness.





#fit.norm.3<-fitdist(life_table_total_babies_infected$babies,"norm")  
#fit.gamma.3<-fitdist(life_table_total_babies_infected$babies,"gamma",lower=c(0,0))
#fit.weibull.3<-fitdist(life_table_total_babies_infected$babies,"weibull",lower=c(0,0))
#fit.lnorm.3<-fitdist(life_table_total_babies_infected$babies,"lnorm")
#fit.pois.3<-fitdist(round(life_table_total_babies_infected$babies),"pois")

#OK, can't make a correlation matrix, but that's ok

#summary(fit.norm.3)## nbinom must be bounded by 1
#summary(fit.gamma.3)
#summary(fit.weibull.3)
#summary(fit.lnorm.3)
#summary(fit.pois.3)


lt_fitness_glm <- lm(babies ~ day_infect + coinf_treat.y + coinf_treat.y:day_infect, data=life_table_total_babies_infected)


summary(lt_fitness_glm)

life_fig <- life_table_max_life_infected
levels(life_fig$coinf_treat.y) <- c("Bacteria Only","Fungi Only","Bacteria First","Fungi First")

tiff(filename = "life_table.tiff", width = 2000, height = 1500, pointsize = 12,res=300)


ggplot(life_fig,aes(y=babies,x=day_infect,color=factor(coinf_treat.y),fill=factor(coinf_treat.y),shape=factor(coinf_treat.y))) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(colour = "black", size=1)) + 
  geom_jitter(width=1,height=.05,size=2) +
  stat_smooth(method="lm",se=TRUE) +
  stat_smooth(method="lm",se=FALSE) +
  #coord_cartesian(ylim = c(0,250)) +
  scale_colour_manual(values = c("#56B4E9","#E69F00", "#009E73", "#F0E442")) +
  scale_fill_manual(values = c("#56B4E9","#E69F00", "#009E73", "#F0E442")) +
  xlab("Age at Infection (Days)") + ylab("Lifetime Offspring") + #ggtitle("Metsch fitness inc. both infectivity and spore yield") +
  theme(axis.title.x = element_text(face="bold", size=16),axis.text.x  = element_text(size=16)) + 
  theme(axis.title.y = element_text(face="bold", size=16),axis.text.y  = element_text(size=16)) +
  theme(legend.title=element_blank())

dev.off()



#### So what's causing these trends? Let's look at babies per day after the hosts are infected, once again compared to average babies per day had by healthy hosts. 
#### Interesting. It appears as though the biggest hit to daily reproduction is the single pasteuria infections, moreso than coinfections early on. 


##Jesus H christ, That birth thing had turned into a hot ass mess. 
##What we need is treatments after first infection date, 
##Then we need to test age, age^2, treat.age, treat.age^2

#coinf_treat is a variable, with 0 = uninf, 1 = pat infect, 2 = metsch infect, 3 = pm, 4= mp


#Births after the first day of infection
life_table_post_infect <- life_table
for(i in 1:length(life_table_post_infect$day_infect)){
  life_table_post_infect$day_infect[i] <- max(life_table_post_infect$day_past_infect[i],life_table_post_infect$day_metsch_infect[i])
}

life_table_post_infect <- life_table_post_infect[which(life_table_post_infect$host_age > life_table_post_infect$day_infect),]

life_table_post_infect$since_infect <- life_table_post_infect$host_age - life_table_post_infect$day_infect

for(i in 1:length(life_table_post_infect$since_infect)){
  if(life_table_post_infect$treatment[i] == 1){life_table_post_infect$since_infect[i] <- 0}
}

#Ok, need to go back and remove individuals who were not succesfully infected 

life_table_post_infect_merge <- merge(useable_reps,life_table_post_infect,by=c("treatment","replicate"))


#Ok, now a lm
bab_per_day_lm <- lm(babies ~ host_age + I(host_age^2) + since_infect:as.factor(coinf_treat.y), data=life_table_post_infect_merge)
summary(bab_per_day_lm)


#try to plot? OK, so what we will do is say that the data points are just a giant,
#unhelpful cloud (it is)
#So will just plot the trend lines, faceted across day of infection

healthy.young.sim <- matrix(0,43,1000)
for(i in 1:1000){
  healthy.young.sim[1,i] <- rnorm(1,mean = coef(summary(bab_per_day_lm))[1,1],sd=coef(summary(bab_per_day_lm))[1,2])
  healthy.young.sim[2,i] <- rnorm(1,mean = coef(summary(bab_per_day_lm))[2,1],sd=coef(summary(bab_per_day_lm))[2,2])
  healthy.young.sim[3,i] <- rnorm(1,mean = coef(summary(bab_per_day_lm))[3,1],sd=coef(summary(bab_per_day_lm))[3,2])
  for(j in 0:39){
    healthy.young.sim[j+4,i] <- healthy.young.sim[1,i] + healthy.young.sim[2,i] * j + healthy.young.sim[3,i] * j^2
  }
}

healthy.young <- matrix(0,40,6)
for(i in 0:39){
  healthy.young[i+1,1] <- i
  healthy.young[i+1,2] <- mean(healthy.young.sim[i+4,])
  healthy.young[i+1,3] <- mean(healthy.young.sim[i+4,]) + sd(healthy.young.sim[i+4,])
  healthy.young[i+1,4] <- mean(healthy.young.sim[i+4,]) - sd(healthy.young.sim[i+4,])
  healthy.young[i+1,5] <- "Not Exposed"
  healthy.young[i+1,6] <- "Young"
}

healthy.middle.sim <- matrix(0,43,1000)
for(i in 1:1000){
  healthy.middle.sim[1,i] <- rnorm(1,mean = coef(summary(bab_per_day_lm))[1,1],sd=coef(summary(bab_per_day_lm))[1,2])
  healthy.middle.sim[2,i] <- rnorm(1,mean = coef(summary(bab_per_day_lm))[2,1],sd=coef(summary(bab_per_day_lm))[2,2])
  healthy.middle.sim[3,i] <- rnorm(1,mean = coef(summary(bab_per_day_lm))[3,1],sd=coef(summary(bab_per_day_lm))[3,2])
  for(j in 0:39){
    healthy.middle.sim[j+4,i] <- healthy.middle.sim[1,i] + healthy.middle.sim[2,i] * j + healthy.middle.sim[3,i] * j^2
  }
}

healthy.middle <- matrix(0,40,6)
for(i in 0:39){
  healthy.middle[i+1,1] <- i
  healthy.middle[i+1,2] <- mean(healthy.middle.sim[i+4,])
  healthy.middle[i+1,3] <- mean(healthy.middle.sim[i+4,]) + sd(healthy.middle.sim[i+4,])
  healthy.middle[i+1,4] <- mean(healthy.middle.sim[i+4,]) - sd(healthy.middle.sim[i+4,])
  healthy.middle[i+1,5] <- "Not Exposed"
  healthy.middle[i+1,6] <- "Middle"
}

healthy.old.sim <- matrix(0,43,1000)
for(i in 1:1000){
  healthy.old.sim[1,i] <- rnorm(1,mean = coef(summary(bab_per_day_lm))[1,1],sd=coef(summary(bab_per_day_lm))[1,2])
  healthy.old.sim[2,i] <- rnorm(1,mean = coef(summary(bab_per_day_lm))[2,1],sd=coef(summary(bab_per_day_lm))[2,2])
  healthy.old.sim[3,i] <- rnorm(1,mean = coef(summary(bab_per_day_lm))[3,1],sd=coef(summary(bab_per_day_lm))[3,2])
  for(j in 0:39){
    healthy.old.sim[j+4,i] <- healthy.old.sim[1,i] + healthy.old.sim[2,i] * j + healthy.old.sim[3,i] * j^2
  }
}

healthy.old <- matrix(0,40,6)
for(i in 0:39){
  healthy.old[i+1,1] <- i
  healthy.old[i+1,2] <- mean(healthy.old.sim[i+4,])
  healthy.old[i+1,3] <- mean(healthy.old.sim[i+4,]) + sd(healthy.old.sim[i+4,])
  healthy.old[i+1,4] <- mean(healthy.old.sim[i+4,]) - sd(healthy.old.sim[i+4,])
  healthy.old[i+1,5] <- "Not Exposed"
  healthy.old[i+1,6] <- "Old"
}

bacteria.first.young.sim <- matrix(0,44,1000)
for(i in 1:1000){
  bacteria.first.young.sim[1,i] <- rnorm(1,mean = coef(summary(bab_per_day_lm))[1,1],sd=coef(summary(bab_per_day_lm))[1,2])
  bacteria.first.young.sim[2,i] <- rnorm(1,mean = coef(summary(bab_per_day_lm))[2,1],sd=coef(summary(bab_per_day_lm))[2,2])
  bacteria.first.young.sim[3,i] <- rnorm(1,mean = coef(summary(bab_per_day_lm))[3,1],sd=coef(summary(bab_per_day_lm))[3,2])
  bacteria.first.young.sim[4,i] <- rnorm(1,mean = coef(summary(bab_per_day_lm))[6,1],sd=coef(summary(bab_per_day_lm))[6,2])
  for(j in 0:39){
    if(j>4){
      bacteria.first.young.sim[j+5,i] <- bacteria.first.young.sim[1,i] + bacteria.first.young.sim[2,i] * j + bacteria.first.young.sim[3,i] * (j^2) + bacteria.first.young.sim[4,i] * (j-4)
    }
    else{bacteria.first.young.sim[j+5,i] <- healthy.young.sim[j+5,i]}
  }
}

bacteria.first.young <- matrix(0,40,6)
for(i in 0:39){
  bacteria.first.young[i+1,1] <- i
  bacteria.first.young[i+1,2] <- mean(bacteria.first.young.sim[i+5,])
  bacteria.first.young[i+1,3] <- mean(bacteria.first.young.sim[i+5,]) + sd(bacteria.first.young.sim[i+5,])
  bacteria.first.young[i+1,4] <- mean(bacteria.first.young.sim[i+5,]) - sd(bacteria.first.young.sim[i+5,])
  bacteria.first.young[i+1,5] <- "Bacteria First"
  bacteria.first.young[i+1,6] <- "Young"
}

bacteria.first.middle.sim <- matrix(0,44,1000)
for(i in 1:1000){
  bacteria.first.middle.sim[1,i] <- rnorm(1,mean = coef(summary(bab_per_day_lm))[1,1],sd=coef(summary(bab_per_day_lm))[1,2])
  bacteria.first.middle.sim[2,i] <- rnorm(1,mean = coef(summary(bab_per_day_lm))[2,1],sd=coef(summary(bab_per_day_lm))[2,2])
  bacteria.first.middle.sim[3,i] <- rnorm(1,mean = coef(summary(bab_per_day_lm))[3,1],sd=coef(summary(bab_per_day_lm))[3,2])
  bacteria.first.middle.sim[4,i] <- rnorm(1,mean = coef(summary(bab_per_day_lm))[6,1],sd=coef(summary(bab_per_day_lm))[6,2])
  for(j in 0:39){
    if(j>12){
      bacteria.first.middle.sim[j+5,i] <- bacteria.first.middle.sim[1,i] + bacteria.first.middle.sim[2,i] * j + bacteria.first.middle.sim[3,i] * (j^2) + bacteria.first.middle.sim[4,i] * (j-12)
    }
    else{bacteria.first.middle.sim[j+5,i] <- healthy.middle.sim[j+5,i]}
  }
}

bacteria.first.middle <- matrix(0,40,6)
for(i in 0:39){
  bacteria.first.middle[i+1,1] <- i
  bacteria.first.middle[i+1,2] <- mean(bacteria.first.middle.sim[i+5,])
  bacteria.first.middle[i+1,3] <- mean(bacteria.first.middle.sim[i+5,]) + sd(bacteria.first.middle.sim[i+5,])
  bacteria.first.middle[i+1,4] <- mean(bacteria.first.middle.sim[i+5,]) - sd(bacteria.first.middle.sim[i+5,])
  bacteria.first.middle[i+1,5] <- "Bacteria First"
  bacteria.first.middle[i+1,6] <- "Middle"
}

bacteria.first.old.sim <- matrix(0,44,1000)
for(i in 1:1000){
  bacteria.first.old.sim[1,i] <- rnorm(1,mean = coef(summary(bab_per_day_lm))[1,1],sd=coef(summary(bab_per_day_lm))[1,2])
  bacteria.first.old.sim[2,i] <- rnorm(1,mean = coef(summary(bab_per_day_lm))[2,1],sd=coef(summary(bab_per_day_lm))[2,2])
  bacteria.first.old.sim[3,i] <- rnorm(1,mean = coef(summary(bab_per_day_lm))[3,1],sd=coef(summary(bab_per_day_lm))[3,2])
  bacteria.first.old.sim[4,i] <- rnorm(1,mean = coef(summary(bab_per_day_lm))[6,1],sd=coef(summary(bab_per_day_lm))[6,2])
  for(j in 0:39){
    if(j>20){
      bacteria.first.old.sim[j+5,i] <- bacteria.first.old.sim[1,i] + bacteria.first.old.sim[2,i] * j + bacteria.first.old.sim[3,i] * (j^2) + bacteria.first.old.sim[4,i] * (j-20)
    }
    else{bacteria.first.old.sim[j+5,i] <- healthy.old.sim[j+5,i]}
  }
}

bacteria.first.old <- matrix(0,40,6)
for(i in 0:39){
  bacteria.first.old[i+1,1] <- i
  bacteria.first.old[i+1,2] <- mean(bacteria.first.old.sim[i+5,])
  bacteria.first.old[i+1,3] <- mean(bacteria.first.old.sim[i+5,]) + sd(bacteria.first.old.sim[i+5,])
  bacteria.first.old[i+1,4] <- mean(bacteria.first.old.sim[i+5,]) - sd(bacteria.first.old.sim[i+5,])
  bacteria.first.old[i+1,5] <- "Bacteria First"
  bacteria.first.old[i+1,6] <- "Old"
}

bacteria.alone.young.sim <- matrix(0,44,1000)
for(i in 1:1000){
  bacteria.alone.young.sim[1,i] <- rnorm(1,mean = coef(summary(bab_per_day_lm))[1,1],sd=coef(summary(bab_per_day_lm))[1,2])
  bacteria.alone.young.sim[2,i] <- rnorm(1,mean = coef(summary(bab_per_day_lm))[2,1],sd=coef(summary(bab_per_day_lm))[2,2])
  bacteria.alone.young.sim[3,i] <- rnorm(1,mean = coef(summary(bab_per_day_lm))[3,1],sd=coef(summary(bab_per_day_lm))[3,2])
  bacteria.alone.young.sim[4,i] <- rnorm(1,mean = coef(summary(bab_per_day_lm))[4,1],sd=coef(summary(bab_per_day_lm))[4,2])
  for(j in 0:39){
    if(j>4){
    bacteria.alone.young.sim[j+5,i] <- bacteria.alone.young.sim[1,i] + bacteria.alone.young.sim[2,i] * j + bacteria.alone.young.sim[3,i] * (j^2) + bacteria.alone.young.sim[4,i] * (j-4)
    }
    else{bacteria.alone.young.sim[j+5,i] <- healthy.young.sim[j+5,i]}
  }
}

bacteria.alone.young <- matrix(0,40,6)
for(i in 0:39){
  bacteria.alone.young[i+1,1] <- i
  bacteria.alone.young[i+1,2] <- mean(bacteria.alone.young.sim[i+5,])
  bacteria.alone.young[i+1,3] <- mean(bacteria.alone.young.sim[i+5,]) + sd(bacteria.alone.young.sim[i+5,])
  bacteria.alone.young[i+1,4] <- mean(bacteria.alone.young.sim[i+5,]) - sd(bacteria.alone.young.sim[i+5,])
  bacteria.alone.young[i+1,5] <- "Bacteria Only"
  bacteria.alone.young[i+1,6] <- "Young"
}

bacteria.alone.middle.sim <- matrix(0,44,1000)
for(i in 1:1000){
  bacteria.alone.middle.sim[1,i] <- rnorm(1,mean = coef(summary(bab_per_day_lm))[1,1],sd=coef(summary(bab_per_day_lm))[1,2])
  bacteria.alone.middle.sim[2,i] <- rnorm(1,mean = coef(summary(bab_per_day_lm))[2,1],sd=coef(summary(bab_per_day_lm))[2,2])
  bacteria.alone.middle.sim[3,i] <- rnorm(1,mean = coef(summary(bab_per_day_lm))[3,1],sd=coef(summary(bab_per_day_lm))[3,2])
  bacteria.alone.middle.sim[4,i] <- rnorm(1,mean = coef(summary(bab_per_day_lm))[4,1],sd=coef(summary(bab_per_day_lm))[4,2])
  for(j in 0:39){
    if(j>12){
      bacteria.alone.middle.sim[j+5,i] <- bacteria.alone.middle.sim[1,i] + bacteria.alone.middle.sim[2,i] * j + bacteria.alone.middle.sim[3,i] * (j^2) + bacteria.alone.middle.sim[4,i] * (j-12)
    }
    else{bacteria.alone.middle.sim[j+5,i] <- healthy.middle.sim[j+5,i]}
  }
}

bacteria.alone.middle <- matrix(0,40,6)
for(i in 0:39){
  bacteria.alone.middle[i+1,1] <- i
  bacteria.alone.middle[i+1,2] <- mean(bacteria.alone.middle.sim[i+5,])
  bacteria.alone.middle[i+1,3] <- mean(bacteria.alone.middle.sim[i+5,]) + sd(bacteria.alone.middle.sim[i+5,])
  bacteria.alone.middle[i+1,4] <- mean(bacteria.alone.middle.sim[i+5,]) - sd(bacteria.alone.middle.sim[i+5,])
  bacteria.alone.middle[i+1,5] <- "Bacteria Only"
  bacteria.alone.middle[i+1,6] <- "Middle"
}

bacteria.alone.old.sim <- matrix(0,44,1000)
for(i in 1:1000){
  bacteria.alone.old.sim[1,i] <- rnorm(1,mean = coef(summary(bab_per_day_lm))[1,1],sd=coef(summary(bab_per_day_lm))[1,2])
  bacteria.alone.old.sim[2,i] <- rnorm(1,mean = coef(summary(bab_per_day_lm))[2,1],sd=coef(summary(bab_per_day_lm))[2,2])
  bacteria.alone.old.sim[3,i] <- rnorm(1,mean = coef(summary(bab_per_day_lm))[3,1],sd=coef(summary(bab_per_day_lm))[3,2])
  bacteria.alone.old.sim[4,i] <- rnorm(1,mean = coef(summary(bab_per_day_lm))[4,1],sd=coef(summary(bab_per_day_lm))[4,2])
  for(j in 0:39){
    if(j>20){
      bacteria.alone.old.sim[j+5,i] <- bacteria.alone.old.sim[1,i] + bacteria.alone.old.sim[2,i] * j + bacteria.alone.old.sim[3,i] * (j^2) + bacteria.alone.old.sim[4,i] * (j-20)
    }
    else{bacteria.alone.old.sim[j+5,i] <- healthy.old.sim[j+5,i]}
  }
}

bacteria.alone.old <- matrix(0,40,6)
for(i in 0:39){
  bacteria.alone.old[i+1,1] <- i
  bacteria.alone.old[i+1,2] <- mean(bacteria.alone.old.sim[i+5,])
  bacteria.alone.old[i+1,3] <- mean(bacteria.alone.old.sim[i+5,]) + sd(bacteria.alone.old.sim[i+5,])
  bacteria.alone.old[i+1,4] <- mean(bacteria.alone.old.sim[i+5,]) - sd(bacteria.alone.old.sim[i+5,])
  bacteria.alone.old[i+1,5] <- "Bacteria Only"
  bacteria.alone.old[i+1,6] <- "Old"
}

fungi.first.young.sim <- matrix(0,44,1000)
for(i in 1:1000){
  fungi.first.young.sim[1,i] <- rnorm(1,mean = coef(summary(bab_per_day_lm))[1,1],sd=coef(summary(bab_per_day_lm))[1,2])
  fungi.first.young.sim[2,i] <- rnorm(1,mean = coef(summary(bab_per_day_lm))[2,1],sd=coef(summary(bab_per_day_lm))[2,2])
  fungi.first.young.sim[3,i] <- rnorm(1,mean = coef(summary(bab_per_day_lm))[3,1],sd=coef(summary(bab_per_day_lm))[3,2])
  fungi.first.young.sim[4,i] <- rnorm(1,mean = coef(summary(bab_per_day_lm))[7,1],sd=coef(summary(bab_per_day_lm))[7,2])
  for(j in 0:39){
    if(j>4){
      fungi.first.young.sim[j+5,i] <- fungi.first.young.sim[1,i] + fungi.first.young.sim[2,i] * j + fungi.first.young.sim[3,i] * (j^2) + fungi.first.young.sim[4,i] * (j-4)
    }
    else{fungi.first.young.sim[j+5,i] <- healthy.young.sim[j+5,i]}
  }
}

fungi.first.young <- matrix(0,40,6)
for(i in 0:39){
  fungi.first.young[i+1,1] <- i
  fungi.first.young[i+1,2] <- mean(fungi.first.young.sim[i+5,])
  fungi.first.young[i+1,3] <- mean(fungi.first.young.sim[i+5,]) + sd(fungi.first.young.sim[i+5,])
  fungi.first.young[i+1,4] <- mean(fungi.first.young.sim[i+5,]) - sd(fungi.first.young.sim[i+5,])
  fungi.first.young[i+1,5] <- "Fungi First"
  fungi.first.young[i+1,6] <- "Young"
}

fungi.first.middle.sim <- matrix(0,44,1000)
for(i in 1:1000){
  fungi.first.middle.sim[1,i] <- rnorm(1,mean = coef(summary(bab_per_day_lm))[1,1],sd=coef(summary(bab_per_day_lm))[1,2])
  fungi.first.middle.sim[2,i] <- rnorm(1,mean = coef(summary(bab_per_day_lm))[2,1],sd=coef(summary(bab_per_day_lm))[2,2])
  fungi.first.middle.sim[3,i] <- rnorm(1,mean = coef(summary(bab_per_day_lm))[3,1],sd=coef(summary(bab_per_day_lm))[3,2])
  fungi.first.middle.sim[4,i] <- rnorm(1,mean = coef(summary(bab_per_day_lm))[7,1],sd=coef(summary(bab_per_day_lm))[7,2])
  for(j in 0:39){
    if(j>12){
      fungi.first.middle.sim[j+5,i] <- fungi.first.middle.sim[1,i] + fungi.first.middle.sim[2,i] * j + fungi.first.middle.sim[3,i] * (j^2) + fungi.first.middle.sim[4,i] * (j-12)
    }
    else{fungi.first.middle.sim[j+5,i] <- healthy.middle.sim[j+5,i]}
  }
}

fungi.first.middle <- matrix(0,40,6)
for(i in 0:39){
  fungi.first.middle[i+1,1] <- i
  fungi.first.middle[i+1,2] <- mean(fungi.first.middle.sim[i+5,])
  fungi.first.middle[i+1,3] <- mean(fungi.first.middle.sim[i+5,]) + sd(fungi.first.middle.sim[i+5,])
  fungi.first.middle[i+1,4] <- mean(fungi.first.middle.sim[i+5,]) - sd(fungi.first.middle.sim[i+5,])
  fungi.first.middle[i+1,5] <- "Fungi First"
  fungi.first.middle[i+1,6] <- "Middle"
}

fungi.first.old.sim <- matrix(0,44,1000)
for(i in 1:1000){
  fungi.first.old.sim[1,i] <- rnorm(1,mean = coef(summary(bab_per_day_lm))[1,1],sd=coef(summary(bab_per_day_lm))[1,2])
  fungi.first.old.sim[2,i] <- rnorm(1,mean = coef(summary(bab_per_day_lm))[2,1],sd=coef(summary(bab_per_day_lm))[2,2])
  fungi.first.old.sim[3,i] <- rnorm(1,mean = coef(summary(bab_per_day_lm))[3,1],sd=coef(summary(bab_per_day_lm))[3,2])
  fungi.first.old.sim[4,i] <- rnorm(1,mean = coef(summary(bab_per_day_lm))[7,1],sd=coef(summary(bab_per_day_lm))[7,2])
  for(j in 0:39){
    if(j>20){
      fungi.first.old.sim[j+5,i] <- fungi.first.old.sim[1,i] + fungi.first.old.sim[2,i] * j + fungi.first.old.sim[3,i] * (j^2) + fungi.first.old.sim[4,i] * (j-20)
    }
    else{fungi.first.old.sim[j+5,i] <- healthy.old.sim[j+5,i]}
  }
}

fungi.first.old <- matrix(0,40,6)
for(i in 0:39){
  fungi.first.old[i+1,1] <- i
  fungi.first.old[i+1,2] <- mean(fungi.first.old.sim[i+5,])
  fungi.first.old[i+1,3] <- mean(fungi.first.old.sim[i+5,]) + sd(fungi.first.old.sim[i+5,])
  fungi.first.old[i+1,4] <- mean(fungi.first.old.sim[i+5,]) - sd(fungi.first.old.sim[i+5,])
  fungi.first.old[i+1,5] <- "Fungi First"
  fungi.first.old[i+1,6] <- "Old"
}

birth.regres <- rbind(healthy.young,healthy.middle,healthy.old,bacteria.alone.young,
                      bacteria.alone.middle,bacteria.alone.old,bacteria.first.young,
                      bacteria.first.middle,bacteria.first.old,fungi.first.young,
                      fungi.first.middle,fungi.first.old)
birth.regres <- as.data.frame(birth.regres)
colnames(birth.regres) <- c("age","mean","max","min","Infection_Treatment","Age_Treatment")
birth.regres$age <- as.numeric(birth.regres$age)
birth.regres$mean <- as.numeric(birth.regres$mean)
birth.regres$min <- as.numeric(birth.regres$min)
birth.regres$max <- as.numeric(birth.regres$max)
birth.regres$Infection_Treatment <- as.factor(birth.regres$Infection_Treatment)
birth.regres$Age_Treatment <- as.factor(birth.regres$Age_Treatment)

#levels(birth.regres$Age_Treatment) <- c("Young","Middle","Old")

birth.regres$Age_Treatment <- factor(birth.regres$Age_Treatment, levels = sort(unique(birth.regres$Age_Treatment)))

tiff(filename = "birth_regres.tiff", width = 3000, height = 1300, pointsize = 12,res=300)

ggplot(birth.regres,aes(y=mean,x=age,color=factor(Infection_Treatment),fill=factor(Infection_Treatment))) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(colour = "black", size=1)) + 
  geom_ribbon(aes(ymin=min,ymax=max,x=age),alpha = 0.25) + 
  geom_line(aes(y=mean,x=age),lwd=2) +
  coord_cartesian(ylim = c(0,1.5)) +
  scale_colour_manual(values = c("#56B4E9","#E69F00", "#009E73", "#F0E442","black")) +
  scale_fill_manual(values = c("#56B4E9","#E69F00", "#009E73", "#F0E442","black")) +
  xlab("Age (Days)") + ylab("Offspring/Day") + #ggtitle("Metsch fitness inc. both infectivity and spore yield") +
  theme(axis.title.x = element_text(face="bold", size=16),axis.text.x  = element_text(size=16)) + 
  theme(axis.title.y = element_text(face="bold", size=16),axis.text.y  = element_text(size=16)) +
  theme(legend.title=element_blank()) +
  facet_wrap(birth.regres$Age_Treatment)

dev.off()

#In summary, age, age^2, age at infect, and age at infect ^2 matter for everything
#mp lowers birth, metsch alone makes age_infect have a negative impact on bab per day
#Past and pm make host_age have a less positive impact on birth
#Can't remove anything for these numbers


#### OK, how about average lifespan of hosts? Once again compared to the lifespan of uninfected hosts. 
#### Most interesting thing: If Pasteuria infects hosts late in their life, it seems to increase their lifespan.



lt_fitness_glm <- lm(host_age ~ day_infect + I(day_infect^2) + coinf_treat.y + coinf_treat.y:day_infect, data=life_table_max_life_infected)

summary(lt_fitness_glm)

ggplot(life_table_max_life_infected,aes(y=host_age,x=day_infect,color=factor(coinf_treat.y),fill=factor(coinf_treat.y))) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(colour = "black", size=1)) + 
  geom_jitter(width=1.5,height=.05) +
  stat_smooth(method="lm",se=TRUE) +
  #coord_cartesian(ylim = c(0,250)) +
  scale_colour_manual(values = c("blue", "black", "green","pink")) +
  scale_fill_manual(values = c("blue", "black", "green", "pink")) +
  xlab("Age at Infection (Days)") + ylab("Lifespan") + #ggtitle("Metsch fitness inc. both infectivity and spore yield") +
  theme(axis.title.x = element_text(face="bold", size=16),axis.text.x  = element_text(size=16)) + 
  theme(axis.title.y = element_text(face="bold", size=16),axis.text.y  = element_text(size=16))


#### Finally, let's examine how time til death after infection depended on age at infection
#### First, let's examine time lived after pasteuria infection
#### So we're really focusing on the single pasteuria infections here, since in the others we expect lifespan to be largely driven by metsch
#### But it does look like pasteuria does shorten lifespan, because all hosts die roughly 30 days after infection. 



life_table_max_life_pasteuria <- life_table_max_life_infected %>%
  filter(day_past_infect > 0)

life_table_max_life_only_pasteuria <- life_table_max_life_pasteuria %>%
  filter(status == 2)

lt_fitness_glm <- lm(post_past_life ~ day_past_infect + I(day_past_infect^2), data=life_table_max_life_only_pasteuria)

summary(lt_fitness_glm)

ggplot(life_table_max_life_only_pasteuria,aes(y=post_past_life,x=day_past_infect)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(colour = "black", size=1)) + 
  geom_jitter(width=1.5,height=.05) +
  stat_smooth(method="lm",se=TRUE) +
  #coord_cartesian(ylim = c(0,250)) +
  scale_colour_manual(values = c("blue", "black", "green","pink")) +
  scale_fill_manual(values = c("blue", "black", "green", "pink")) +
  xlab("Age at Infection (Days)") + ylab("Lifespan") + #ggtitle("Metsch fitness inc. both infectivity and spore yield") +
  theme(axis.title.x = element_text(face="bold", size=16),axis.text.x  = element_text(size=16)) + 
  theme(axis.title.y = element_text(face="bold", size=16),axis.text.y  = element_text(size=16))


#### Now, let's examine time lived after metschnikowia infection
#### Looks like doesn't much effect time til death here



life_table_max_life_metschnikowia <- life_table_max_life_infected %>%
  filter(day_metsch_infect > 0)

lt_fitness_glm <- lm(post_metsch_life ~ day_metsch_infect + I(day_metsch_infect^2) + coinf_treat.y + coinf_treat.y:day_metsch_infect, data=life_table_max_life_metschnikowia)



summary(lt_fitness_glm)

ggplot(life_table_max_life_metschnikowia,aes(y=post_metsch_life,x=day_metsch_infect,color=factor(coinf_treat.y),fill=factor(coinf_treat.y))) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(colour = "black", size=1)) + 
  geom_jitter(width=1.5,height=.05) +
  stat_smooth(method="lm",se=TRUE) +
  #coord_cartesian(ylim = c(0,250)) +
  scale_colour_manual(values = c("blue", "black", "green","pink")) +
  scale_fill_manual(values = c("blue", "black", "green", "pink")) +
  xlab("Age at Infection (Days)") + ylab("Lifespan") + #ggtitle("Metsch fitness inc. both infectivity and spore yield") +
  theme(axis.title.x = element_text(face="bold", size=16),axis.text.x  = element_text(size=16)) + 
  theme(axis.title.y = element_text(face="bold", size=16),axis.text.y  = element_text(size=16))


#### Alright, As before, We want to get a general sense of whether conditions were different between this experiment and the AmNat experiment
#### The following figures show less resources going to reproduction, and thus less chance for pasteuria to facilitate metschnikowia via castration. 

#### Left is boxplot of lifespan of amnat daphnia, right is current daphnia, for healthy individuals
#### Current daphnia lived much longer

```{r , echo=FALSE}
t.test(old_data$dayslife,life_table_max_life_healthy$host_age)
boxplot(old_data$dayslife,life_table_max_life_healthy$host_age)

```

#### Left is boxplot of lifetime reproduction of amnat daphnia, right is current daphnia, for healthy individuals
#### Current daphnia had much fewer babies
#### Thus less reproductive energy for pasteuria to highjack. 

```{r , echo=FALSE}
t.test(old_data$totbab,life_table_total_babies_healthy$babies)
boxplot(old_data$totbab,life_table_total_babies_healthy$babies)

```

#### Left is boxplot of daily reproduction of amnat daphnia, right is current daphnia, for healthy individuals
#### Current daphnia had much fewer babies per day as well. 

```{r , echo=FALSE}
t.test(old_data$babperday,life_table_average_healthy$babies)
boxplot(old_data$babperday,life_table_average_healthy$babies)

```