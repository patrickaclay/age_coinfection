library(here)
library(ggplot2)
############Comparison#######

##NetLogo output is awful, need to skip first19 rows
coin_comp <- read.csv("degradation_exploration.csv",header = FALSE,skip = 24)

long <- ((length(coin_comp)-1)/4)*2000
coin_comp_long <- matrix(0,long,4)
for(i in 1:long){
  coin_comp_long[i,1] = coin_comp[ifelse(i-(floor(i/2000)*2000)==0,2000,i-(floor(i/2000)*2000)),ifelse(i-(floor(i/2000)*2000)==0,2+(floor(i/2000)-1)*4,2+floor(i/2000)*4)]
  coin_comp_long[i,2] = coin_comp[ifelse(i-(floor(i/2000)*2000)==0,2000,i-(floor(i/2000)*2000)),ifelse(i-(floor(i/2000)*2000)==0,3+(floor(i/2000)-1)*4,3+floor(i/2000)*4)]
  coin_comp_long[i,3] = coin_comp[ifelse(i-(floor(i/2000)*2000)==0,2000,i-(floor(i/2000)*2000)),ifelse(i-(floor(i/2000)*2000)==0,4+(floor(i/2000)-1)*4,4+floor(i/2000)*4)]
  coin_comp_long[i,4] = coin_comp[ifelse(i-(floor(i/2000)*2000)==0,2000,i-(floor(i/2000)*2000)),ifelse(i-(floor(i/2000)*2000)==0,5+(floor(i/2000)-1)*4,5+floor(i/2000)*4)]
}

##OK, now need to put in the variables
coin_comp_long <- as.data.frame(coin_comp_long)

#Name state var columns
colnames(coin_comp_long) <- c('S','P','M','C')

long <- length(coin_comp_long$S)
#days column
coin_comp_long$days <- 0
for(i in 1:long){
  coin_comp_long$days[i] = ifelse(i-(floor(i/1000)*1000)==0,1000,i-(floor(i/1000)*1000))
  
}

##run number
coin_comp_long$run <- 0
for(i in 1:long){
  coin_comp_long$run[i] = floor(((i-1)/1000))+1
}

#switches
coin_comp_long$feedback <- 1

for(i in 1:long){
  if((coin_comp_long$run[i] %% 2) == 1 ){coin_comp_long$feedback[i] <- 0}
}

#put in age
coin_comp_long$age <- "4 days"

for(i in 1:long){
  if(coin_comp_long$run[i] > 96 & coin_comp_long$run[i] < 193){coin_comp_long$age[i] <- "12 days"}
  if(coin_comp_long$run[i] > 192){coin_comp_long$age[i] <- "20 days"}
}

#put in Treatment
coin_comp_long$Treatment <- "Coinfection"
#so c c p p m m
for(i in 1:long){
  if((floor((coin_comp_long$run[i]+1)/2) %% 3) == 2){coin_comp_long$Treatment[i] <- "Only Bacteria"}
  if((floor((coin_comp_long$run[i]+1)/2) %% 3) == 0){coin_comp_long$Treatment[i] <- "Only Fungi"}
}

#Put in degradation rates

coin_comp_long$Pdeg <- 0.35
#in chunks of 12 * 2
for(i in 1:long){
  if((floor((coin_comp_long$run[i]+1)/24) %% 4) == 2){coin_comp_long$Pdeg[i] <- 0.45}
  if((floor((coin_comp_long$run[i]+1)/24) %% 4) == 3){coin_comp_long$Pdeg[i] <- 0.55}
  if((floor((coin_comp_long$run[i]+1)/24) %% 4) == 0){coin_comp_long$Pdeg[i] <- 0.65}
}

coin_comp_long$Mdeg <- 0.06
#in chunks of 3 * 2
for(i in 1:long){
  if((floor((coin_comp_long$run[i]+1)/6) %% 4) == 2){coin_comp_long$Mdeg[i] <- 0.11}
  if((floor((coin_comp_long$run[i]+1)/6) %% 4) == 3){coin_comp_long$Mdeg[i] <- 0.26}
  if((floor((coin_comp_long$run[i]+1)/6) %% 4) == 0){coin_comp_long$Mdeg[i] <- 0.6}
}

##now make N and Mprev and Pprev
coin_comp_long$N <- (coin_comp_long$S + coin_comp_long$M + coin_comp_long$P + coin_comp_long$C)
coin_comp_long$M_prev <- (coin_comp_long$M + coin_comp_long$C)/(coin_comp_long$S + coin_comp_long$M + coin_comp_long$P + coin_comp_long$C)
coin_comp_long$P_prev <- (coin_comp_long$P + coin_comp_long$C)/(coin_comp_long$S + coin_comp_long$M + coin_comp_long$P + coin_comp_long$C)

write.csv(coin_comp_long,"degradation_compare_formatted.csv")
coin_comp_long <- read.csv("degradation_compare_formatted.csv")

##OK, only focus on last 500 days
coin_comp_long <- coin_comp_long[which(coin_comp_long$days > 500),]

#OK, how do we want to do this. response = prevalence, input = model, color = Treatment
#panel == degradation

#MAke a column for age or feedback
coin_comp_long$model <- "Age Effects"
for(i in 1:length(coin_comp_long$run)){
  if(coin_comp_long$feedback[i] == 0){coin_comp_long$model[i] <- coin_comp_long$age[i]}
  }
coin_comp_long$model <- as.factor(coin_comp_long$model)
levels(coin_comp_long$model) 
coin_comp_long$model <- factor(coin_comp_long$model, levels = c("4 days","12 days","20 days","Age Effects"))
#OK that worked

#Make sure for feedback, only do if prior age was 4, so that doesn't start extinct

coin_comp_long <- coin_comp_long[which(coin_comp_long$age == "4 days" | coin_comp_long$feedback == 0),]

ggplot(coin_comp_long,aes(y=M_prev,x=model,color = Treatment)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(colour = "black", size=1)) + 
  geom_boxplot() +
  #coord_cartesian(ylim = c(0,250)) +
  scale_colour_manual(values = c("blue", "dark green", "green")) +
  xlab("Age (Days)") + ylab("Prevalence") + #ggtitle("Feeding Feedback") +
  theme(axis.title.x = element_text(face="bold", size=16),axis.text.x  = element_text(size=16)) + 
  theme(axis.title.y = element_text(face="bold", size=16),axis.text.y  = element_text(size=16)) +
  facet_grid(Mdeg ~ Pdeg)

ggplot(coin_comp_long,aes(y=P_prev,x=model,color = Treatment)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(colour = "black", size=1)) + 
  geom_boxplot() +
  #coord_cartesian(ylim = c(0,250)) +
  scale_colour_manual(values = c("blue", "dark green", "green")) +
  xlab("Age (Days)") + ylab("Prevalence") + #ggtitle("Feeding Feedback") +
  theme(axis.title.x = element_text(face="bold", size=16),axis.text.x  = element_text(size=16)) + 
  theme(axis.title.y = element_text(face="bold", size=16),axis.text.y  = element_text(size=16)) +
  facet_grid(Mdeg ~ Pdeg)

##OK, so now we need to make these figures.

coin_comp_long <- coin_comp_long[which(coin_comp_long$Mdeg == 0.06 & coin_comp_long$Pdeg == 0.35),]

tiff(filename = "fungi_prev.tiff", width = 2000, height = 1500, pointsize = 12,res=300)

ggplot(coin_comp_long,aes(y=M_prev,x=model,color = Treatment)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(colour = "black", size=1)) + 
  geom_boxplot() +
  #coord_cartesian(ylim = c(0,250)) +
  scale_colour_manual(values = c("black","#56B4E9","#E69F00")) +
  xlab("Model") + ylab("Fungi Prevalence") + #ggtitle("Feeding Feedback") +
  theme(axis.title.x = element_text(face="bold", size=16),axis.text.x  = element_text(size=16)) + 
  theme(axis.title.y = element_text(face="bold", size=16),axis.text.y  = element_text(size=16),
        legend.text = element_text( size=12), legend.title = element_text(face="bold", size=12)) 

dev.off()

tiff(filename = "bacteria_prev.tiff", width = 2000, height = 1500, pointsize = 12,res=300)

ggplot(coin_comp_long,aes(y=P_prev,x=model,color = Treatment)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(colour = "black", size=1)) + 
  geom_boxplot() +
  #coord_cartesian(ylim = c(0,250)) +
  scale_colour_manual(values = c("black","#56B4E9","#E69F00")) +
  xlab("Model") + ylab("Bacteria Prevalence") + #ggtitle("Feeding Feedback") +
  theme(axis.title.x = element_text(face="bold", size=16),axis.text.x  = element_text(size=16)) + 
  theme(axis.title.y = element_text(face="bold", size=16),axis.text.y  = element_text(size=16),
        legend.text = element_text( size=12), legend.title = element_text(face="bold", size=12))

dev.off()

#######################################################################
########## Do It Again But start from initial conditions ##############
#######################################################################


##NetLogo output is awful, need to skip first19 rows
coin_comp <- read.csv("degradation_exploration_low_start.csv",header = FALSE,skip = 25)

long <- ((length(coin_comp)-1)/4)*1000
coin_comp_long <- matrix(0,long,4)
for(i in 1:long){
  coin_comp_long[i,1] = coin_comp[ifelse(i-(floor(i/1000)*1000)==0,1000,i-(floor(i/1000)*1000)),ifelse(i-(floor(i/1000)*1000)==0,2+(floor(i/1000)-1)*4,2+floor(i/1000)*4)]
  coin_comp_long[i,2] = coin_comp[ifelse(i-(floor(i/1000)*1000)==0,1000,i-(floor(i/1000)*1000)),ifelse(i-(floor(i/1000)*1000)==0,3+(floor(i/1000)-1)*4,3+floor(i/1000)*4)]
  coin_comp_long[i,3] = coin_comp[ifelse(i-(floor(i/1000)*1000)==0,1000,i-(floor(i/1000)*1000)),ifelse(i-(floor(i/1000)*1000)==0,4+(floor(i/1000)-1)*4,4+floor(i/1000)*4)]
  coin_comp_long[i,4] = coin_comp[ifelse(i-(floor(i/1000)*1000)==0,1000,i-(floor(i/1000)*1000)),ifelse(i-(floor(i/1000)*1000)==0,5+(floor(i/1000)-1)*4,5+floor(i/1000)*4)]
}

##OK, now need to put in the variables
coin_comp_long <- as.data.frame(coin_comp_long)

#Name state var columns
colnames(coin_comp_long) <- c('S','P','M','C')

long <- length(coin_comp_long$S)
#days column
coin_comp_long$days <- 0
for(i in 1:long){
  coin_comp_long$days[i] = ifelse(i-(floor(i/1000)*1000)==0,1000,i-(floor(i/1000)*1000))
  
}

##run number
coin_comp_long$run <- 0
for(i in 1:long){
  coin_comp_long$run[i] = floor(((i-1)/1000))+1
}

#switches
coin_comp_long$feedback <- 1

for(i in 1:long){
  if(coin_comp_long$run[i] > 144 ){coin_comp_long$feedback[i] <- 0}
}

#put in age
coin_comp_long$age <- "4 days"

for(i in 1:long){
  if(coin_comp_long$run[i] > 48 & coin_comp_long$run[i] < 97){coin_comp_long$age[i] <- "12 days"}
  if(coin_comp_long$run[i] > 96 & coin_comp_long$run[i] < 145){coin_comp_long$age[i] <- "20 days"}
  if(coin_comp_long$run[i] > 192 & coin_comp_long$run[i] < 241){coin_comp_long$age[i] <- "12 days"}
  if(coin_comp_long$run[i] > 240){coin_comp_long$age[i] <- "20 days"}
}

#put in Treatment
coin_comp_long$Treatment <- "Coinfection"
#so c c p p m m
for(i in 1:long){
  if((floor((coin_comp_long$run[i])) %% 3) == 2){coin_comp_long$Treatment[i] <- "Only Bacteria"}
  if((floor((coin_comp_long$run[i])) %% 3) == 0){coin_comp_long$Treatment[i] <- "Only Fungi"}
}

#Put in degradation rates

coin_comp_long$Pdeg <- 0.35
#in chunks of 12 * 2
for(i in 1:long){
  if((floor((coin_comp_long$run[i]+1)/12) %% 4) == 2){coin_comp_long$Pdeg[i] <- 0.45}
  if((floor((coin_comp_long$run[i]+1)/12) %% 4) == 3){coin_comp_long$Pdeg[i] <- 0.55}
  if((floor((coin_comp_long$run[i]+1)/12) %% 4) == 0){coin_comp_long$Pdeg[i] <- 0.65}
}

coin_comp_long$Mdeg <- 0.06
#in chunks of 3 * 2
for(i in 1:long){
  if((floor((coin_comp_long$run[i]+1)/3) %% 4) == 2){coin_comp_long$Mdeg[i] <- 0.11}
  if((floor((coin_comp_long$run[i]+1)/3) %% 4) == 3){coin_comp_long$Mdeg[i] <- 0.26}
  if((floor((coin_comp_long$run[i]+1)/3) %% 4) == 0){coin_comp_long$Mdeg[i] <- 0.6}
}

##now make N and Mprev and Pprev
coin_comp_long$N <- (coin_comp_long$S + coin_comp_long$M + coin_comp_long$P + coin_comp_long$C)
coin_comp_long$M_prev <- (coin_comp_long$M + coin_comp_long$C)/(coin_comp_long$S + coin_comp_long$M + coin_comp_long$P + coin_comp_long$C)
coin_comp_long$P_prev <- (coin_comp_long$P + coin_comp_long$C)/(coin_comp_long$S + coin_comp_long$M + coin_comp_long$P + coin_comp_long$C)


##OK, only focus on last 500 days
coin_comp_long <- coin_comp_long[which(coin_comp_long$days > 500),]

#OK, how do we want to do this. response = prevalence, input = model, color = Treatment
#panel == degradation

#MAke a column for age or feedback
coin_comp_long$model <- "Age Effects"
for(i in 1:length(coin_comp_long$run)){
  if(coin_comp_long$feedback[i] == 0){coin_comp_long$model[i] <- coin_comp_long$age[i]}
}
coin_comp_long$model <- as.factor(coin_comp_long$model)
levels(coin_comp_long$model) 
coin_comp_long$model <- factor(coin_comp_long$model, levels = c("4 days","12 days","20 days","Age Effects"))
#OK that worked

#Make sure for feedback, only do if prior age was 4, so that doesn't start extinct

coin_comp_long <- coin_comp_long[which(coin_comp_long$age == "4 days" | coin_comp_long$feedback == 0),]

ggplot(coin_comp_long,aes(y=M_prev,x=model,color = Treatment)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(colour = "black", size=1)) + 
  geom_boxplot() +
  #coord_cartesian(ylim = c(0,250)) +
  scale_colour_manual(values = c("blue", "dark green", "green")) +
  xlab("Age (Days)") + ylab("Prevalence") + #ggtitle("Feeding Feedback") +
  theme(axis.title.x = element_text(face="bold", size=16),axis.text.x  = element_text(size=16)) + 
  theme(axis.title.y = element_text(face="bold", size=16),axis.text.y  = element_text(size=16)) +
  facet_grid(Mdeg ~ Pdeg)

ggplot(coin_comp_long,aes(y=P_prev,x=model,color = Treatment)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(colour = "black", size=1)) + 
  geom_boxplot() +
  #coord_cartesian(ylim = c(0,250)) +
  scale_colour_manual(values = c("blue", "dark green", "green")) +
  xlab("Age (Days)") + ylab("Prevalence") + #ggtitle("Feeding Feedback") +
  theme(axis.title.x = element_text(face="bold", size=16),axis.text.x  = element_text(size=16)) + 
  theme(axis.title.y = element_text(face="bold", size=16),axis.text.y  = element_text(size=16)) +
  facet_grid(Mdeg ~ Pdeg)

##OK, so now we need to make these figures.

coin_comp_long <- coin_comp_long[which(coin_comp_long$Mdeg == 0.06 & coin_comp_long$Pdeg == 0.35),]

tiff(filename = "fungi_prev.tiff", width = 2000, height = 1500, pointsize = 12,res=300)

ggplot(coin_comp_long,aes(y=M_prev,x=model,color = Treatment)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(colour = "black", size=1)) + 
  geom_boxplot() +
  #coord_cartesian(ylim = c(0,250)) +
  scale_colour_manual(values = c("black","#56B4E9","#E69F00")) +
  xlab("Model") + ylab("Fungi Prevalence") + #ggtitle("Feeding Feedback") +
  theme(axis.title.x = element_text(face="bold", size=16),axis.text.x  = element_text(size=16)) + 
  theme(axis.title.y = element_text(face="bold", size=16),axis.text.y  = element_text(size=16),
        legend.text = element_text( size=12), legend.title = element_text(face="bold", size=12)) 

dev.off()

tiff(filename = "bacteria_prev.tiff", width = 2000, height = 1500, pointsize = 12,res=300)

ggplot(coin_comp_long,aes(y=P_prev,x=model,color = Treatment)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(colour = "black", size=1)) + 
  geom_boxplot() +
  #coord_cartesian(ylim = c(0,250)) +
  scale_colour_manual(values = c("black","#56B4E9","#E69F00")) +
  xlab("Model") + ylab("Bacteria Prevalence") + #ggtitle("Feeding Feedback") +
  theme(axis.title.x = element_text(face="bold", size=16),axis.text.x  = element_text(size=16)) + 
  theme(axis.title.y = element_text(face="bold", size=16),axis.text.y  = element_text(size=16),
        legend.text = element_text( size=12), legend.title = element_text(face="bold", size=12))

dev.off()


















##################################################################
##OK, plot ten runs with or without feedback

##NetLogo output is awful, need to skip first19 rows
coin_comp <- read.csv("main_plot_ten_runs.csv",header = FALSE,skip = 25)

long <- ((length(coin_comp)-1)/4)*1000
coin_comp_long <- matrix(0,long,4)
for(i in 1:long){
  coin_comp_long[i,1] = coin_comp[ifelse(i-(floor(i/1000)*1000)==0,1000,i-(floor(i/1000)*1000)),ifelse(i-(floor(i/1000)*1000)==0,2+(floor(i/1000)-1)*4,2+floor(i/1000)*4)]
  coin_comp_long[i,2] = coin_comp[ifelse(i-(floor(i/1000)*1000)==0,1000,i-(floor(i/1000)*1000)),ifelse(i-(floor(i/1000)*1000)==0,3+(floor(i/1000)-1)*4,3+floor(i/1000)*4)]
  coin_comp_long[i,3] = coin_comp[ifelse(i-(floor(i/1000)*1000)==0,1000,i-(floor(i/1000)*1000)),ifelse(i-(floor(i/1000)*1000)==0,4+(floor(i/1000)-1)*4,4+floor(i/1000)*4)]
  coin_comp_long[i,4] = coin_comp[ifelse(i-(floor(i/1000)*1000)==0,1000,i-(floor(i/1000)*1000)),ifelse(i-(floor(i/1000)*1000)==0,5+(floor(i/1000)-1)*4,5+floor(i/1000)*4)]
}

##OK, now need to put in the variables
coin_comp_long <- as.data.frame(coin_comp_long)

#Name state var columns
colnames(coin_comp_long) <- c('S','P','M','C')

long <- length(coin_comp_long$S)
#days column
coin_comp_long$days <- 0
for(i in 1:long){
  coin_comp_long$days[i] = ifelse(i-(floor(i/1000)*1000)==0,1000,i-(floor(i/1000)*1000))
  
}

##run number
coin_comp_long$run <- 0
for(i in 1:long){
  coin_comp_long$run[i] = floor(((i-1)/1000))+1
}

#switches
coin_comp_long$feedback <- 1

for(i in 1:long){
  if(coin_comp_long$run[i]  > 90 ){coin_comp_long$feedback[i] <- 0}
}

#put in age
coin_comp_long$age <- "4 days"

for(i in 1:long){
  if(coin_comp_long$run[i] > 30 & coin_comp_long$run[i] < 61){coin_comp_long$age[i] <- "12 days"}
  if(coin_comp_long$run[i] > 60 & coin_comp_long$run[i] < 91){coin_comp_long$age[i] <- "20 days"}
  if(coin_comp_long$run[i] > 120 & coin_comp_long$run[i] < 151){coin_comp_long$age[i] <- "12 days"}
  if(coin_comp_long$run[i] > 150){coin_comp_long$age[i] <- "20 days"}
}

#put in Treatment
coin_comp_long$Treatment <- "Coinfection"
#so c c p p m m
for(i in 1:long){
  if((floor((coin_comp_long$run[i]+1)/10) %% 3) == 2){coin_comp_long$Treatment[i] <- "Only Bacteria"}
  if((floor((coin_comp_long$run[i]+1)/10) %% 3) == 0){coin_comp_long$Treatment[i] <- "Only Fungi"}
}

#Put in degradation rates



##now make N and Mprev and Pprev
coin_comp_long$N <- (coin_comp_long$S + coin_comp_long$M + coin_comp_long$P + coin_comp_long$C)
coin_comp_long$M_prev <- (coin_comp_long$M + coin_comp_long$C)/(coin_comp_long$S + coin_comp_long$M + coin_comp_long$P + coin_comp_long$C)
coin_comp_long$P_prev <- (coin_comp_long$P + coin_comp_long$C)/(coin_comp_long$S + coin_comp_long$M + coin_comp_long$P + coin_comp_long$C)

##OK, only focus on last 500 days
coin_comp_long <- coin_comp_long[which(coin_comp_long$days > 500),]

##Ok, plot


ggplot(coin_comp_long,aes(y=M_prev,x=as.factor(run),fill=as.factor(age))) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(colour = "black", size=1)) + 
  geom_boxplot() +
  #coord_cartesian(ylim = c(0,250)) +
  #scale_colour_manual(values = c("blue", "dark green", "green")) +
  xlab("Age (Days)") + ylab("Prevalence") + #ggtitle("Feeding Feedback") +
  theme(axis.title.x = element_text(face="bold", size=16),axis.text.x  = element_text(size=16)) + 
  theme(axis.title.y = element_text(face="bold", size=16),axis.text.y  = element_text(size=16)) 


ggplot(coin_comp_long,aes(y=P_prev,x=as.factor(run),fill=as.factor(age))) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(colour = "black", size=1)) + 
  geom_boxplot() +
  #coord_cartesian(ylim = c(0,250)) +
  #scale_colour_manual(values = c("blue", "dark green", "green")) +
  xlab("Age (Days)") + ylab("Prevalence") + #ggtitle("Feeding Feedback") +
  theme(axis.title.x = element_text(face="bold", size=16),axis.text.x  = element_text(size=16)) + 
  theme(axis.title.y = element_text(face="bold", size=16),axis.text.y  = element_text(size=16)) 


#######################################################################
##### Looking at stuff across a variety of initial conditions #########
#######################################################################


##NetLogo output is awful, need to skip first19 rows
coin_comp <- read.csv("Initial_conditions_M_sing.csv",header = FALSE,skip = 27)

long <- ((length(coin_comp)-1)/4)*1000
coin_comp_long <- matrix(0,long,4)
for(i in 1:long){
  coin_comp_long[i,1] = coin_comp[ifelse(i-(floor(i/1000)*1000)==0,1000,i-(floor(i/1000)*1000)),ifelse(i-(floor(i/1000)*1000)==0,2+(floor(i/1000)-1)*4,2+floor(i/1000)*4)]
  coin_comp_long[i,2] = coin_comp[ifelse(i-(floor(i/1000)*1000)==0,1000,i-(floor(i/1000)*1000)),ifelse(i-(floor(i/1000)*1000)==0,3+(floor(i/1000)-1)*4,3+floor(i/1000)*4)]
  coin_comp_long[i,3] = coin_comp[ifelse(i-(floor(i/1000)*1000)==0,1000,i-(floor(i/1000)*1000)),ifelse(i-(floor(i/1000)*1000)==0,4+(floor(i/1000)-1)*4,4+floor(i/1000)*4)]
  coin_comp_long[i,4] = coin_comp[ifelse(i-(floor(i/1000)*1000)==0,1000,i-(floor(i/1000)*1000)),ifelse(i-(floor(i/1000)*1000)==0,5+(floor(i/1000)-1)*4,5+floor(i/1000)*4)]
}

##OK, now need to put in the variables
coin_comp_long <- as.data.frame(coin_comp_long)

#Name state var columns
colnames(coin_comp_long) <- c('S','P','M','C')

long <- length(coin_comp_long$S)
#days column
coin_comp_long$days <- 0
for(i in 1:long){
  coin_comp_long$days[i] = ifelse(i-(floor(i/1000)*1000)==0,1000,i-(floor(i/1000)*1000))
  
}

##run number
coin_comp_long$run <- 0
for(i in 1:long){
  coin_comp_long$run[i] = floor(((i-1)/1000))+1
}

#switches
coin_comp_long$feedback <- 1

for(i in 1:long){
  if(coin_comp_long$run[i]  > (max(coin_comp_long$run)/2) ){coin_comp_long$feedback[i] <- 0}
}

#put in age
coin_comp_long$age <- "4 days"

for(i in 1:long){
  if(coin_comp_long$run[i] > 9 & coin_comp_long$run[i] < 19){coin_comp_long$age[i] <- "12 days"}
  if(coin_comp_long$run[i] > 18 & coin_comp_long$run[i] < 28){coin_comp_long$age[i] <- "20 days"}
  if(coin_comp_long$run[i] > 36 & coin_comp_long$run[i] < 46){coin_comp_long$age[i] <- "12 days"}
  if(coin_comp_long$run[i] > 45 & coin_comp_long$run[i] < 55){coin_comp_long$age[i] <- "20 days"}
}

#put in Treatment
coin_comp_long$Treatment <- "Single"

#Mstart
coin_comp_long$Mstart <- 10

for(i in 1:long){
  if((coin_comp_long$run[i] %% 9) == 2){coin_comp_long$Mstart[i] <- 100}
  if((coin_comp_long$run[i] %% 9) == 3){coin_comp_long$Mstart[i] <- 200}
  if((coin_comp_long$run[i] %% 9) == 4){coin_comp_long$Mstart[i] <- 300}
  if((coin_comp_long$run[i] %% 9) == 5){coin_comp_long$Mstart[i] <- 400}
  if((coin_comp_long$run[i] %% 9) == 6){coin_comp_long$Mstart[i] <- 500}
  if((coin_comp_long$run[i] %% 9) == 7){coin_comp_long$Mstart[i] <- 600}
  if((coin_comp_long$run[i] %% 9) == 8){coin_comp_long$Mstart[i] <- 700}
  if((coin_comp_long$run[i] %% 9) == 0){coin_comp_long$Mstart[i] <- 900}
}


##now make N and Mprev and Pprev
coin_comp_long$N <- (coin_comp_long$S + coin_comp_long$M + coin_comp_long$P + coin_comp_long$C)
coin_comp_long$M_prev <- (coin_comp_long$M + coin_comp_long$C)/(coin_comp_long$S + coin_comp_long$M + coin_comp_long$P + coin_comp_long$C)
coin_comp_long$P_prev <- (coin_comp_long$P + coin_comp_long$C)/(coin_comp_long$S + coin_comp_long$M + coin_comp_long$P + coin_comp_long$C)


##OK, only focus on last 500 days
coin_comp_long <- coin_comp_long[which(coin_comp_long$days > 500),]

#OK, how do we want to do this. response = prevalence, input = model, color = Treatment
#panel == degradation

#MAke a column for age or feedback
coin_comp_long$model <- "Age Effects"
for(i in 1:length(coin_comp_long$run)){
  if(coin_comp_long$feedback[i] == 0){coin_comp_long$model[i] <- coin_comp_long$age[i]}
}
coin_comp_long$model <- as.factor(coin_comp_long$model)
levels(coin_comp_long$model) 
coin_comp_long$model <- factor(coin_comp_long$model, levels = c("4 days","12 days","20 days","Age Effects"))
#OK that worked

#Make sure for feedback, only do if prior age was 4, so that doesn't start extinct

coin_comp_long <- coin_comp_long[which(coin_comp_long$age == "4 days" | coin_comp_long$feedback == 0),]

ggplot(coin_comp_long,aes(y=M_prev,x=as.factor(Mstart))) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(colour = "black", size=1)) + 
  geom_boxplot() +
  #coord_cartesian(ylim = c(0,250)) +
  scale_colour_manual(values = c("blue", "dark green", "green")) +
  xlab("Age (Days)") + ylab("Prevalence") + #ggtitle("Feeding Feedback") +
  theme(axis.title.x = element_text(face="bold", size=16),axis.text.x  = element_text(size=16)) + 
  theme(axis.title.y = element_text(face="bold", size=16),axis.text.y  = element_text(size=16)) +
  facet_wrap(coin_comp_long$model)


##NetLogo output is awful, need to skip first19 rows
coin_comp <- read.csv("Initial_conditions_M_co.csv",header = FALSE,skip = 27)

long <- ((length(coin_comp)-1)/4)*1000
coin_comp_long <- matrix(0,long,4)
for(i in 1:long){
  coin_comp_long[i,1] = coin_comp[ifelse(i-(floor(i/1000)*1000)==0,1000,i-(floor(i/1000)*1000)),ifelse(i-(floor(i/1000)*1000)==0,2+(floor(i/1000)-1)*4,2+floor(i/1000)*4)]
  coin_comp_long[i,2] = coin_comp[ifelse(i-(floor(i/1000)*1000)==0,1000,i-(floor(i/1000)*1000)),ifelse(i-(floor(i/1000)*1000)==0,3+(floor(i/1000)-1)*4,3+floor(i/1000)*4)]
  coin_comp_long[i,3] = coin_comp[ifelse(i-(floor(i/1000)*1000)==0,1000,i-(floor(i/1000)*1000)),ifelse(i-(floor(i/1000)*1000)==0,4+(floor(i/1000)-1)*4,4+floor(i/1000)*4)]
  coin_comp_long[i,4] = coin_comp[ifelse(i-(floor(i/1000)*1000)==0,1000,i-(floor(i/1000)*1000)),ifelse(i-(floor(i/1000)*1000)==0,5+(floor(i/1000)-1)*4,5+floor(i/1000)*4)]
}

##OK, now need to put in the variables
coin_comp_long <- as.data.frame(coin_comp_long)

#Name state var columns
colnames(coin_comp_long) <- c('S','P','M','C')

long <- length(coin_comp_long$S)
#days column
coin_comp_long$days <- 0
for(i in 1:long){
  coin_comp_long$days[i] = ifelse(i-(floor(i/1000)*1000)==0,1000,i-(floor(i/1000)*1000))
  
}

##run number
coin_comp_long$run <- 0
for(i in 1:long){
  coin_comp_long$run[i] = floor(((i-1)/1000))+1
}

#switches
coin_comp_long$feedback <- 1

for(i in 1:long){
  if(coin_comp_long$run[i]  > (max(coin_comp_long$run)/2) ){coin_comp_long$feedback[i] <- 0}
}

#put in age
coin_comp_long$age <- "4 days"

for(i in 1:long){
  if(coin_comp_long$run[i] > 9 & coin_comp_long$run[i] < 19){coin_comp_long$age[i] <- "12 days"}
  if(coin_comp_long$run[i] > 18 & coin_comp_long$run[i] < 28){coin_comp_long$age[i] <- "20 days"}
  if(coin_comp_long$run[i] > 36 & coin_comp_long$run[i] < 46){coin_comp_long$age[i] <- "12 days"}
  if(coin_comp_long$run[i] > 45 & coin_comp_long$run[i] < 55){coin_comp_long$age[i] <- "20 days"}
}

#put in Treatment
coin_comp_long$Treatment <- "Single"

#Mstart
coin_comp_long$Mstart <- 10

for(i in 1:long){
  if((coin_comp_long$run[i] %% 9) == 2){coin_comp_long$Mstart[i] <- 100}
  if((coin_comp_long$run[i] %% 9) == 3){coin_comp_long$Mstart[i] <- 200}
  if((coin_comp_long$run[i] %% 9) == 4){coin_comp_long$Mstart[i] <- 300}
  if((coin_comp_long$run[i] %% 9) == 5){coin_comp_long$Mstart[i] <- 400}
  if((coin_comp_long$run[i] %% 9) == 6){coin_comp_long$Mstart[i] <- 500}
  if((coin_comp_long$run[i] %% 9) == 7){coin_comp_long$Mstart[i] <- 600}
  if((coin_comp_long$run[i] %% 9) == 8){coin_comp_long$Mstart[i] <- 700}
  if((coin_comp_long$run[i] %% 9) == 0){coin_comp_long$Mstart[i] <- 900}
}


##now make N and Mprev and Pprev
coin_comp_long$N <- (coin_comp_long$S + coin_comp_long$M + coin_comp_long$P + coin_comp_long$C)
coin_comp_long$M_prev <- (coin_comp_long$M + coin_comp_long$C)/(coin_comp_long$S + coin_comp_long$M + coin_comp_long$P + coin_comp_long$C)
coin_comp_long$P_prev <- (coin_comp_long$P + coin_comp_long$C)/(coin_comp_long$S + coin_comp_long$M + coin_comp_long$P + coin_comp_long$C)


##OK, only focus on last 500 days
coin_comp_long <- coin_comp_long[which(coin_comp_long$days > 500),]

#OK, how do we want to do this. response = prevalence, input = model, color = Treatment
#panel == degradation

#MAke a column for age or feedback
coin_comp_long$model <- "Age Effects"
for(i in 1:length(coin_comp_long$run)){
  if(coin_comp_long$feedback[i] == 0){coin_comp_long$model[i] <- coin_comp_long$age[i]}
}
coin_comp_long$model <- as.factor(coin_comp_long$model)
levels(coin_comp_long$model) 
coin_comp_long$model <- factor(coin_comp_long$model, levels = c("4 days","12 days","20 days","Age Effects"))
#OK that worked

#Make sure for feedback, only do if prior age was 4, so that doesn't start extinct

coin_comp_long <- coin_comp_long[which(coin_comp_long$age == "4 days" | coin_comp_long$feedback == 0),]

ggplot(coin_comp_long,aes(y=M_prev,x=as.factor(Mstart))) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(colour = "black", size=1)) + 
  geom_boxplot() +
  #coord_cartesian(ylim = c(0,250)) +
  scale_colour_manual(values = c("blue", "dark green", "green")) +
  xlab("Age (Days)") + ylab("Prevalence") + #ggtitle("Feeding Feedback") +
  theme(axis.title.x = element_text(face="bold", size=16),axis.text.x  = element_text(size=16)) + 
  theme(axis.title.y = element_text(face="bold", size=16),axis.text.y  = element_text(size=16)) +
  facet_wrap(coin_comp_long$model)


##NetLogo output is awful, need to skip first19 rows
coin_comp <- read.csv("Initial_conditions_P_sing.csv",header = FALSE,skip = 27)

long <- ((length(coin_comp)-1)/4)*1000
coin_comp_long <- matrix(0,long,4)
for(i in 1:long){
  coin_comp_long[i,1] = coin_comp[ifelse(i-(floor(i/1000)*1000)==0,1000,i-(floor(i/1000)*1000)),ifelse(i-(floor(i/1000)*1000)==0,2+(floor(i/1000)-1)*4,2+floor(i/1000)*4)]
  coin_comp_long[i,2] = coin_comp[ifelse(i-(floor(i/1000)*1000)==0,1000,i-(floor(i/1000)*1000)),ifelse(i-(floor(i/1000)*1000)==0,3+(floor(i/1000)-1)*4,3+floor(i/1000)*4)]
  coin_comp_long[i,3] = coin_comp[ifelse(i-(floor(i/1000)*1000)==0,1000,i-(floor(i/1000)*1000)),ifelse(i-(floor(i/1000)*1000)==0,4+(floor(i/1000)-1)*4,4+floor(i/1000)*4)]
  coin_comp_long[i,4] = coin_comp[ifelse(i-(floor(i/1000)*1000)==0,1000,i-(floor(i/1000)*1000)),ifelse(i-(floor(i/1000)*1000)==0,5+(floor(i/1000)-1)*4,5+floor(i/1000)*4)]
}

##OK, now need to put in the variables
coin_comp_long <- as.data.frame(coin_comp_long)

#Name state var columns
colnames(coin_comp_long) <- c('S','P','M','C')

long <- length(coin_comp_long$S)
#days column
coin_comp_long$days <- 0
for(i in 1:long){
  coin_comp_long$days[i] = ifelse(i-(floor(i/1000)*1000)==0,1000,i-(floor(i/1000)*1000))
  
}

##run number
coin_comp_long$run <- 0
for(i in 1:long){
  coin_comp_long$run[i] = floor(((i-1)/1000))+1
}

#switches
coin_comp_long$feedback <- 1

for(i in 1:long){
  if(coin_comp_long$run[i]  > (max(coin_comp_long$run)/2) ){coin_comp_long$feedback[i] <- 0}
}

#put in age
coin_comp_long$age <- "4 days"

for(i in 1:long){
  if(coin_comp_long$run[i] > 9 & coin_comp_long$run[i] < 19){coin_comp_long$age[i] <- "12 days"}
  if(coin_comp_long$run[i] > 18 & coin_comp_long$run[i] < 28){coin_comp_long$age[i] <- "20 days"}
  if(coin_comp_long$run[i] > 36 & coin_comp_long$run[i] < 46){coin_comp_long$age[i] <- "12 days"}
  if(coin_comp_long$run[i] > 45 & coin_comp_long$run[i] < 55){coin_comp_long$age[i] <- "20 days"}
}

#put in Treatment
coin_comp_long$Treatment <- "Single"

#Mstart
coin_comp_long$Pstart <- 10

for(i in 1:long){
  if((coin_comp_long$run[i] %% 9) == 2){coin_comp_long$Pstart[i] <- 100}
  if((coin_comp_long$run[i] %% 9) == 3){coin_comp_long$Pstart[i] <- 200}
  if((coin_comp_long$run[i] %% 9) == 4){coin_comp_long$Pstart[i] <- 300}
  if((coin_comp_long$run[i] %% 9) == 5){coin_comp_long$Pstart[i] <- 400}
  if((coin_comp_long$run[i] %% 9) == 6){coin_comp_long$Pstart[i] <- 500}
  if((coin_comp_long$run[i] %% 9) == 7){coin_comp_long$Pstart[i] <- 600}
  if((coin_comp_long$run[i] %% 9) == 8){coin_comp_long$Pstart[i] <- 700}
  if((coin_comp_long$run[i] %% 9) == 0){coin_comp_long$Pstart[i] <- 900}
}


##now make N and Mprev and Pprev
coin_comp_long$N <- (coin_comp_long$S + coin_comp_long$M + coin_comp_long$P + coin_comp_long$C)
coin_comp_long$M_prev <- (coin_comp_long$M + coin_comp_long$C)/(coin_comp_long$S + coin_comp_long$M + coin_comp_long$P + coin_comp_long$C)
coin_comp_long$P_prev <- (coin_comp_long$P + coin_comp_long$C)/(coin_comp_long$S + coin_comp_long$M + coin_comp_long$P + coin_comp_long$C)


##OK, only focus on last 500 days
coin_comp_long <- coin_comp_long[which(coin_comp_long$days > 500),]

#OK, how do we want to do this. response = prevalence, input = model, color = Treatment
#panel == degradation

#MAke a column for age or feedback
coin_comp_long$model <- "Age Effects"
for(i in 1:length(coin_comp_long$run)){
  if(coin_comp_long$feedback[i] == 0){coin_comp_long$model[i] <- coin_comp_long$age[i]}
}
coin_comp_long$model <- as.factor(coin_comp_long$model)
levels(coin_comp_long$model) 
coin_comp_long$model <- factor(coin_comp_long$model, levels = c("4 days","12 days","20 days","Age Effects"))
#OK that worked

#Make sure for feedback, only do if prior age was 4, so that doesn't start extinct

coin_comp_long <- coin_comp_long[which(coin_comp_long$age == "4 days" | coin_comp_long$feedback == 0),]

ggplot(coin_comp_long,aes(y=P_prev,x=as.factor(Pstart))) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(colour = "black", size=1)) + 
  geom_boxplot() +
  #coord_cartesian(ylim = c(0,250)) +
  scale_colour_manual(values = c("blue", "dark green", "green")) +
  xlab("Age (Days)") + ylab("Prevalence") + #ggtitle("Feeding Feedback") +
  theme(axis.title.x = element_text(face="bold", size=16),axis.text.x  = element_text(size=16)) + 
  theme(axis.title.y = element_text(face="bold", size=16),axis.text.y  = element_text(size=16)) +
  facet_wrap(coin_comp_long$model)

##NetLogo output is awful, need to skip first19 rows
coin_comp <- read.csv("Initial_conditions_P_co.csv",header = FALSE,skip = 27)

long <- ((length(coin_comp)-1)/4)*1000
coin_comp_long <- matrix(0,long,4)
for(i in 1:long){
  coin_comp_long[i,1] = coin_comp[ifelse(i-(floor(i/1000)*1000)==0,1000,i-(floor(i/1000)*1000)),ifelse(i-(floor(i/1000)*1000)==0,2+(floor(i/1000)-1)*4,2+floor(i/1000)*4)]
  coin_comp_long[i,2] = coin_comp[ifelse(i-(floor(i/1000)*1000)==0,1000,i-(floor(i/1000)*1000)),ifelse(i-(floor(i/1000)*1000)==0,3+(floor(i/1000)-1)*4,3+floor(i/1000)*4)]
  coin_comp_long[i,3] = coin_comp[ifelse(i-(floor(i/1000)*1000)==0,1000,i-(floor(i/1000)*1000)),ifelse(i-(floor(i/1000)*1000)==0,4+(floor(i/1000)-1)*4,4+floor(i/1000)*4)]
  coin_comp_long[i,4] = coin_comp[ifelse(i-(floor(i/1000)*1000)==0,1000,i-(floor(i/1000)*1000)),ifelse(i-(floor(i/1000)*1000)==0,5+(floor(i/1000)-1)*4,5+floor(i/1000)*4)]
}

##OK, now need to put in the variables
coin_comp_long <- as.data.frame(coin_comp_long)

#Name state var columns
colnames(coin_comp_long) <- c('S','P','M','C')

long <- length(coin_comp_long$S)
#days column
coin_comp_long$days <- 0
for(i in 1:long){
  coin_comp_long$days[i] = ifelse(i-(floor(i/1000)*1000)==0,1000,i-(floor(i/1000)*1000))
  
}

##run number
coin_comp_long$run <- 0
for(i in 1:long){
  coin_comp_long$run[i] = floor(((i-1)/1000))+1
}

#switches
coin_comp_long$feedback <- 1

for(i in 1:long){
  if(coin_comp_long$run[i]  > (max(coin_comp_long$run)/2) ){coin_comp_long$feedback[i] <- 0}
}

#put in age
coin_comp_long$age <- "4 days"

for(i in 1:long){
  if(coin_comp_long$run[i] > 9 & coin_comp_long$run[i] < 19){coin_comp_long$age[i] <- "12 days"}
  if(coin_comp_long$run[i] > 18 & coin_comp_long$run[i] < 28){coin_comp_long$age[i] <- "20 days"}
  if(coin_comp_long$run[i] > 36 & coin_comp_long$run[i] < 46){coin_comp_long$age[i] <- "12 days"}
  if(coin_comp_long$run[i] > 45 & coin_comp_long$run[i] < 55){coin_comp_long$age[i] <- "20 days"}
}

#put in Treatment
coin_comp_long$Treatment <- "Single"

#Mstart
coin_comp_long$Pstart <- 10

for(i in 1:long){
  if((coin_comp_long$run[i] %% 9) == 2){coin_comp_long$Pstart[i] <- 100}
  if((coin_comp_long$run[i] %% 9) == 3){coin_comp_long$Pstart[i] <- 200}
  if((coin_comp_long$run[i] %% 9) == 4){coin_comp_long$Pstart[i] <- 300}
  if((coin_comp_long$run[i] %% 9) == 5){coin_comp_long$Pstart[i] <- 400}
  if((coin_comp_long$run[i] %% 9) == 6){coin_comp_long$Pstart[i] <- 500}
  if((coin_comp_long$run[i] %% 9) == 7){coin_comp_long$Pstart[i] <- 600}
  if((coin_comp_long$run[i] %% 9) == 8){coin_comp_long$Pstart[i] <- 700}
  if((coin_comp_long$run[i] %% 9) == 0){coin_comp_long$Pstart[i] <- 900}
}


##now make N and Mprev and Pprev
coin_comp_long$N <- (coin_comp_long$S + coin_comp_long$M + coin_comp_long$P + coin_comp_long$C)
coin_comp_long$M_prev <- (coin_comp_long$M + coin_comp_long$C)/(coin_comp_long$S + coin_comp_long$M + coin_comp_long$P + coin_comp_long$C)
coin_comp_long$P_prev <- (coin_comp_long$P + coin_comp_long$C)/(coin_comp_long$S + coin_comp_long$M + coin_comp_long$P + coin_comp_long$C)


##OK, only focus on last 500 days
coin_comp_long <- coin_comp_long[which(coin_comp_long$days > 500),]

#OK, how do we want to do this. response = prevalence, input = model, color = Treatment
#panel == degradation

#MAke a column for age or feedback
coin_comp_long$model <- "Age Effects"
for(i in 1:length(coin_comp_long$run)){
  if(coin_comp_long$feedback[i] == 0){coin_comp_long$model[i] <- coin_comp_long$age[i]}
}
coin_comp_long$model <- as.factor(coin_comp_long$model)
levels(coin_comp_long$model) 
coin_comp_long$model <- factor(coin_comp_long$model, levels = c("4 days","12 days","20 days","Age Effects"))
#OK that worked

#Make sure for feedback, only do if prior age was 4, so that doesn't start extinct

coin_comp_long <- coin_comp_long[which(coin_comp_long$age == "4 days" | coin_comp_long$feedback == 0),]

ggplot(coin_comp_long,aes(y=P_prev,x=as.factor(Pstart))) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(colour = "black", size=1)) + 
  geom_boxplot() +
  #coord_cartesian(ylim = c(0,250)) +
  scale_colour_manual(values = c("blue", "dark green", "green")) +
  xlab("Age (Days)") + ylab("Prevalence") + #ggtitle("Feeding Feedback") +
  theme(axis.title.x = element_text(face="bold", size=16),axis.text.x  = element_text(size=16)) + 
  theme(axis.title.y = element_text(face="bold", size=16),axis.text.y  = element_text(size=16)) +
  facet_wrap(coin_comp_long$model)

############################################################333#####
############################FEEDBACK EXPLORE FINAL##################
####################################################################

##NetLogo output is awful, need to skip first19 rows
coin_comp <- read.csv("feedback_explore_high_initial.csv",header = FALSE,skip = 27)

long <- ((length(coin_comp)-1)/4)*3000
coin_comp_long <- matrix(0,long,4)
for(i in 1:long){
  coin_comp_long[i,1] = coin_comp[ifelse(i-(floor(i/3000)*3000)==0,3000,i-(floor(i/3000)*3000)),ifelse(i-(floor(i/3000)*3000)==0,2+(floor(i/3000)-1)*4,2+floor(i/3000)*4)]
  coin_comp_long[i,2] = coin_comp[ifelse(i-(floor(i/3000)*3000)==0,3000,i-(floor(i/3000)*3000)),ifelse(i-(floor(i/3000)*3000)==0,3+(floor(i/3000)-1)*4,3+floor(i/3000)*4)]
  coin_comp_long[i,3] = coin_comp[ifelse(i-(floor(i/3000)*3000)==0,3000,i-(floor(i/3000)*3000)),ifelse(i-(floor(i/3000)*3000)==0,4+(floor(i/3000)-1)*4,4+floor(i/3000)*4)]
  coin_comp_long[i,4] = coin_comp[ifelse(i-(floor(i/3000)*3000)==0,3000,i-(floor(i/3000)*3000)),ifelse(i-(floor(i/3000)*3000)==0,5+(floor(i/3000)-1)*4,5+floor(i/3000)*4)]
}

##OK, now need to put in the variables
coin_comp_long <- as.data.frame(coin_comp_long)

#Name state var columns
colnames(coin_comp_long) <- c('S','P','M','C')

long <- length(coin_comp_long$S)
#days column
coin_comp_long$days <- 0
for(i in 1:long){
  coin_comp_long$days[i] = ifelse(i-(floor(i/3000)*3000)==0,3000,i-(floor(i/3000)*3000))
  
}

##run number
coin_comp_long$run <- 0
for(i in 1:long){
  coin_comp_long$run[i] = floor(((i-1)/3000))+1
}

#switches
coin_comp_long$infect_feedback <- 1
coin_comp_long$fecund_feedback <- 1
coin_comp_long$spore_feedback <- 1
coin_comp_long$feed_feedback <- 1

for(i in 1:long){
  if((floor((coin_comp_long$run[i]-1)/72) %% 2) == 1){coin_comp_long$infect_feedback[i] <- 0}
  if((floor((coin_comp_long$run[i]-1)/36) %% 2) == 1){coin_comp_long$fecund_feedback[i] <- 0}
  if((floor((coin_comp_long$run[i]-1)/18) %% 2) == 1){coin_comp_long$spore_feedback[i] <- 0}
  if((floor((coin_comp_long$run[i]-1)/9) %% 2) == 1){coin_comp_long$feed_feedback[i] <- 0}
}

#put in age
coin_comp_long$age <- "4 days"

for(i in 1:long){
  if((floor((coin_comp_long$run[i]-1)/3) %% 3) == 1){coin_comp_long$age[i] <- "12 days"}
  if((floor((coin_comp_long$run[i]-1)/3) %% 3) == 2){coin_comp_long$age[i] <- "20 days"}
}

#put in Treatment
coin_comp_long$Treatment <- "Coinfection"
for(i in 1:long){
  if((floor((coin_comp_long$run[i]+1)/1) %% 3) == 0){coin_comp_long$Treatment[i] <- "Only Bacteria"}
  if((floor((coin_comp_long$run[i]+1)/1) %% 3) == 1){coin_comp_long$Treatment[i] <- "Only Fungi"}
}

##now make N and Mprev and Pprev
coin_comp_long$N <- (coin_comp_long$S + coin_comp_long$M + coin_comp_long$P + coin_comp_long$C)
coin_comp_long$M_prev <- (coin_comp_long$M + coin_comp_long$C)/(coin_comp_long$S + coin_comp_long$M + coin_comp_long$P + coin_comp_long$C)
coin_comp_long$P_prev <- (coin_comp_long$P + coin_comp_long$C)/(coin_comp_long$S + coin_comp_long$M + coin_comp_long$P + coin_comp_long$C)

##OK, only focus on last 500 days
coin_comp_long <- coin_comp_long[which(coin_comp_long$days > 2500),]


#MAke a column for age or feedback
coin_comp_long$model <- "Age Effects"
for(i in 1:length(coin_comp_long$run)){
  if(coin_comp_long$infect_feedback[i] == 0 & coin_comp_long$fecund_feedback[i] == 0 & 
     coin_comp_long$feed_feedback[i] == 0 & coin_comp_long$spore_feedback[i] == 0)
    {coin_comp_long$model[i] <- coin_comp_long$age[i]}
}
coin_comp_long$model <- as.factor(coin_comp_long$model)
levels(coin_comp_long$model) 
coin_comp_long$model <- factor(coin_comp_long$model, levels = c("4 days","12 days","20 days","Age Effects"))
#OK that worked
coin_comp_long$age <- factor(coin_comp_long$age, levels = c("4 days","12 days","20 days"))

#Make sure for feedback, only do if prior age was 4, so that doesn't start extinct

coin_comp_long_main <- coin_comp_long[which((coin_comp_long$infect_feedback == 1 & coin_comp_long$fecund_feedback == 1 & 
                                               coin_comp_long$feed_feedback == 1 & coin_comp_long$spore_feedback == 1) | 
                                              (coin_comp_long$infect_feedback == 0 & coin_comp_long$fecund_feedback == 0 & 
                                                 coin_comp_long$feed_feedback == 0 & coin_comp_long$spore_feedback == 0)),]
#In order to make sure that the feedback model doesn't have three runs included, filter out 
coin_comp_long_main <- coin_comp_long_main[which(coin_comp_long_main$model != "Age Effects" | coin_comp_long_main$age == "20 days"),]
#I checked, they are all robust

tiff(filename = "fungi_prev.tiff", width = 2000, height = 1500, pointsize = 12,res=300)

ggplot(coin_comp_long_main,aes(y=M_prev,x=model,color = Treatment)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(colour = "black", size=1)) + 
  geom_boxplot() +
  coord_cartesian(ylim = c(0,1)) +
  scale_colour_manual(values = c("black","#56B4E9","#E69F00")) +
  xlab("Model") + ylab("Fungi Prevalence") + #ggtitle("Feeding Feedback") +
  theme(axis.title.x = element_text(face="bold", size=16),axis.text.x  = element_text(size=16)) + 
  theme(axis.title.y = element_text(face="bold", size=16),axis.text.y  = element_text(size=16),
        legend.text = element_text( size=12), legend.title = element_text(face="bold", size=12))

dev.off()

tiff(filename = "bacteria_prev.tiff", width = 2000, height = 1500, pointsize = 12,res=300)

ggplot(coin_comp_long_main,aes(y=P_prev,x=model,color = Treatment)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(colour = "black", size=1)) + 
  geom_boxplot() +
  coord_cartesian(ylim = c(0,1)) +
  scale_colour_manual(values = c("black","#56B4E9","#E69F00")) +
  xlab("Model") + ylab("Bacteria Prevalence") + #ggtitle("Feeding Feedback") +
  theme(axis.title.x = element_text(face="bold", size=16),axis.text.x  = element_text(size=16)) + 
  theme(axis.title.y = element_text(face="bold", size=16),axis.text.y  = element_text(size=16),
        legend.text = element_text( size=12), legend.title = element_text(face="bold", size=12))

dev.off()

##Only feed_feedback

coin_comp_long_main <- coin_comp_long[which((coin_comp_long$infect_feedback == 0 & coin_comp_long$fecund_feedback == 0 & 
                                               coin_comp_long$feed_feedback == 1 & coin_comp_long$spore_feedback == 0) | 
                                              (coin_comp_long$infect_feedback == 0 & coin_comp_long$fecund_feedback == 0 & 
                                                 coin_comp_long$feed_feedback == 0 & coin_comp_long$spore_feedback == 0)),]


tiff(filename = "fungi_feed_feedback.tiff", width = 2500, height = 1500, pointsize = 12,res=300)

ggplot(coin_comp_long_main,aes(y=M_prev,x=age,color = Treatment)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(colour = "black", size=1)) + 
  geom_boxplot() +
  #coord_cartesian(ylim = c(0,250)) +
  scale_colour_manual(values = c("black","#56B4E9","#E69F00")) +
  xlab("Model") + ylab("Fungi Prevalence") + #ggtitle("Feeding Feedback") +
  theme(axis.title.x = element_text(face="bold", size=16),axis.text.x  = element_text(size=16)) + 
  theme(axis.title.y = element_text(face="bold", size=16),axis.text.y  = element_text(size=16),
        legend.text = element_text( size=12), legend.title = element_text(face="bold", size=12)) +
  facet_wrap(coin_comp_long_main$feed_feedback)

dev.off()

tiff(filename = "bacteria_feed_feedback.tiff", width = 2500, height = 1500, pointsize = 12,res=300)

ggplot(coin_comp_long_main,aes(y=P_prev,x=age,color = Treatment)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(colour = "black", size=1)) + 
  geom_boxplot() +
  #coord_cartesian(ylim = c(0,250)) +
  scale_colour_manual(values = c("black","#56B4E9","#E69F00")) +
  xlab("Model") + ylab("Bacteria Prevalence") + #ggtitle("Feeding Feedback") +
  theme(axis.title.x = element_text(face="bold", size=16),axis.text.x  = element_text(size=16)) + 
  theme(axis.title.y = element_text(face="bold", size=16),axis.text.y  = element_text(size=16),
        legend.text = element_text( size=12), legend.title = element_text(face="bold", size=12))+
  facet_wrap(coin_comp_long_main$feed_feedback)

dev.off()


##Only infect_feedback

coin_comp_long_main <- coin_comp_long[which((coin_comp_long$infect_feedback == 1 & coin_comp_long$fecund_feedback == 0 & 
                                               coin_comp_long$feed_feedback == 0 & coin_comp_long$spore_feedback == 0) | 
                                              (coin_comp_long$infect_feedback == 0 & coin_comp_long$fecund_feedback == 0 & 
                                                 coin_comp_long$feed_feedback == 0 & coin_comp_long$spore_feedback == 0)),]


tiff(filename = "fungi_suscept_feedback.tiff", width = 2500, height = 1500, pointsize = 12,res=300)

ggplot(coin_comp_long_main,aes(y=M_prev,x=age,color = Treatment)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(colour = "black", size=1)) + 
  geom_boxplot() +
  #coord_cartesian(ylim = c(0,250)) +
  scale_colour_manual(values = c("black","#56B4E9","#E69F00")) +
  xlab("Model") + ylab("Fungi Prevalence") + #ggtitle("Susceptibility Feedback") +
  theme(axis.title.x = element_text(face="bold", size=16),axis.text.x  = element_text(size=16)) + 
  theme(axis.title.y = element_text(face="bold", size=16),axis.text.y  = element_text(size=16),
        legend.text = element_text( size=12), legend.title = element_text(face="bold", size=12)) +
  facet_wrap(coin_comp_long_main$infect_feedback)

dev.off()

tiff(filename = "bacteria_suscept_feedback.tiff", width = 2500, height = 1500, pointsize = 12,res=300)

ggplot(coin_comp_long_main,aes(y=P_prev,x=age,color = Treatment)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(colour = "black", size=1)) + 
  geom_boxplot() +
  #coord_cartesian(ylim = c(0,250)) +
  scale_colour_manual(values = c("black","#56B4E9","#E69F00")) +
  xlab("Model") + ylab("Bacteria Prevalence") + #ggtitle("Feeding Feedback") +
  theme(axis.title.x = element_text(face="bold", size=16),axis.text.x  = element_text(size=16)) + 
  theme(axis.title.y = element_text(face="bold", size=16),axis.text.y  = element_text(size=16),
        legend.text = element_text( size=12), legend.title = element_text(face="bold", size=12))+
  facet_wrap(coin_comp_long_main$infect_feedback)

dev.off()

##Only fecund_feedback

coin_comp_long_main <- coin_comp_long[which((coin_comp_long$infect_feedback == 0 & coin_comp_long$fecund_feedback == 1 & 
                                               coin_comp_long$feed_feedback == 0 & coin_comp_long$spore_feedback == 0) | 
                                              (coin_comp_long$infect_feedback == 0 & coin_comp_long$fecund_feedback == 0 & 
                                                 coin_comp_long$feed_feedback == 0 & coin_comp_long$spore_feedback == 0)),]


tiff(filename = "fungi_fecund_feedback.tiff", width = 2500, height = 1500, pointsize = 12,res=300)

ggplot(coin_comp_long_main,aes(y=M_prev,x=age,color = Treatment)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(colour = "black", size=1)) + 
  geom_boxplot() +
  #coord_cartesian(ylim = c(0,250)) +
  scale_colour_manual(values = c("black","#56B4E9","#E69F00")) +
  xlab("Model") + ylab("Fungi Prevalence") + #ggtitle("Feeding Feedback") +
  theme(axis.title.x = element_text(face="bold", size=16),axis.text.x  = element_text(size=16)) + 
  theme(axis.title.y = element_text(face="bold", size=16),axis.text.y  = element_text(size=16),
        legend.text = element_text( size=12), legend.title = element_text(face="bold", size=12)) +
  facet_wrap(coin_comp_long_main$fecund_feedback)

dev.off()

tiff(filename = "bacteria_fecund_feedback.tiff", width = 2500, height = 1500, pointsize = 12,res=300)

ggplot(coin_comp_long_main,aes(y=P_prev,x=age,color = Treatment)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(colour = "black", size=1)) + 
  geom_boxplot() +
  #coord_cartesian(ylim = c(0,250)) +
  scale_colour_manual(values = c("black","#56B4E9","#E69F00")) +
  xlab("Model") + ylab("Bacteria Prevalence") + #ggtitle("Feeding Feedback") +
  theme(axis.title.x = element_text(face="bold", size=16),axis.text.x  = element_text(size=16)) + 
  theme(axis.title.y = element_text(face="bold", size=16),axis.text.y  = element_text(size=16),
        legend.text = element_text( size=12), legend.title = element_text(face="bold", size=12))+
  facet_wrap(coin_comp_long_main$fecund_feedback)

dev.off()
##Only spore_feedback

coin_comp_long_main <- coin_comp_long[which((coin_comp_long$infect_feedback == 0 & coin_comp_long$fecund_feedback == 0 & 
                                               coin_comp_long$feed_feedback == 0 & coin_comp_long$spore_feedback == 1) | 
                                              (coin_comp_long$infect_feedback == 0 & coin_comp_long$fecund_feedback == 0 & 
                                                 coin_comp_long$feed_feedback == 0 & coin_comp_long$spore_feedback == 0)),]


tiff(filename = "fungi_spore_feedback.tiff", width = 2500, height = 1500, pointsize = 12,res=300)

ggplot(coin_comp_long_main,aes(y=M_prev,x=age,color = Treatment)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(colour = "black", size=1)) + 
  geom_boxplot() +
  #coord_cartesian(ylim = c(0,250)) +
  scale_colour_manual(values = c("black","#56B4E9","#E69F00")) +
  xlab("Model") + ylab("Fungi Prevalence") + #ggtitle("Feeding Feedback") +
  theme(axis.title.x = element_text(face="bold", size=16),axis.text.x  = element_text(size=16)) + 
  theme(axis.title.y = element_text(face="bold", size=16),axis.text.y  = element_text(size=16),
        legend.text = element_text( size=12), legend.title = element_text(face="bold", size=12)) +
  facet_wrap(coin_comp_long_main$spore_feedback)

dev.off()

tiff(filename = "bacteria_spore_feedback.tiff", width = 2500, height = 1500, pointsize = 12,res=300)

ggplot(coin_comp_long_main,aes(y=P_prev,x=age,color = Treatment)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(colour = "black", size=1)) + 
  geom_boxplot() +
  #coord_cartesian(ylim = c(0,250)) +
  scale_colour_manual(values = c("black","#56B4E9","#E69F00")) +
  xlab("Model") + ylab("Bacteria Prevalence") + #ggtitle("Feeding Feedback") +
  theme(axis.title.x = element_text(face="bold", size=16),axis.text.x  = element_text(size=16)) + 
  theme(axis.title.y = element_text(face="bold", size=16),axis.text.y  = element_text(size=16),
        legend.text = element_text( size=12), legend.title = element_text(face="bold", size=12))+
  facet_wrap(coin_comp_long_main$spore_feedback)

dev.off()

##No feed_feedback

coin_comp_long_main <- coin_comp_long[which((coin_comp_long$infect_feedback == 1 & coin_comp_long$fecund_feedback == 1 & 
                                               coin_comp_long$feed_feedback == 0 & coin_comp_long$spore_feedback == 1) | 
                                              (coin_comp_long$infect_feedback == 0 & coin_comp_long$fecund_feedback == 0 & 
                                                 coin_comp_long$feed_feedback == 0 & coin_comp_long$spore_feedback == 0)),]


ggplot(coin_comp_long_main,aes(y=M_prev,x=age,color = Treatment)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(colour = "black", size=1)) + 
  geom_boxplot() +
  #coord_cartesian(ylim = c(0,250)) +
  scale_colour_manual(values = c("blue", "dark green", "green")) +
  xlab("Age (Days)") + ylab("Prevalence") + ggtitle("No Feeding Feedback") +
  theme(axis.title.x = element_text(face="bold", size=16),axis.text.x  = element_text(size=16)) + 
  theme(axis.title.y = element_text(face="bold", size=16),axis.text.y  = element_text(size=16)) +
  facet_wrap(coin_comp_long_main$infect_feedback)

ggplot(coin_comp_long_main,aes(y=P_prev,x=age,color = Treatment)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(colour = "black", size=1)) + 
  geom_boxplot() +
  #coord_cartesian(ylim = c(0,250)) +
  scale_colour_manual(values = c("blue", "dark green", "green")) +
  xlab("Age (Days)") + ylab("Prevalence") + ggtitle("No Feeding Feedback") +
  theme(axis.title.x = element_text(face="bold", size=16),axis.text.x  = element_text(size=16)) + 
  theme(axis.title.y = element_text(face="bold", size=16),axis.text.y  = element_text(size=16)) +
  facet_wrap(coin_comp_long_main$infect_feedback)


##No infect_feedback

coin_comp_long_main <- coin_comp_long[which((coin_comp_long$infect_feedback == 0 & coin_comp_long$fecund_feedback == 1 & 
                                               coin_comp_long$feed_feedback == 1 & coin_comp_long$spore_feedback == 1) | 
                                              (coin_comp_long$infect_feedback == 0 & coin_comp_long$fecund_feedback == 0 & 
                                                 coin_comp_long$feed_feedback == 0 & coin_comp_long$spore_feedback == 0)),]


ggplot(coin_comp_long_main,aes(y=M_prev,x=age,color = Treatment)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(colour = "black", size=1)) + 
  geom_boxplot() +
  #coord_cartesian(ylim = c(0,250)) +
  scale_colour_manual(values = c("blue", "dark green", "green")) +
  xlab("Age (Days)") + ylab("Prevalence") + ggtitle("No Infect Feedback") +
  theme(axis.title.x = element_text(face="bold", size=16),axis.text.x  = element_text(size=16)) + 
  theme(axis.title.y = element_text(face="bold", size=16),axis.text.y  = element_text(size=16)) +
  facet_wrap(coin_comp_long_main$feed_feedback)

ggplot(coin_comp_long_main,aes(y=P_prev,x=age,color = Treatment)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(colour = "black", size=1)) + 
  geom_boxplot() +
  #coord_cartesian(ylim = c(0,250)) +
  scale_colour_manual(values = c("blue", "dark green", "green")) +
  xlab("Age (Days)") + ylab("Prevalence") + ggtitle("Infect Feedback") +
  theme(axis.title.x = element_text(face="bold", size=16),axis.text.x  = element_text(size=16)) + 
  theme(axis.title.y = element_text(face="bold", size=16),axis.text.y  = element_text(size=16)) +
  facet_wrap(coin_comp_long_main$feed_feedback)


##No fecund_feedback

coin_comp_long_main <- coin_comp_long[which((coin_comp_long$infect_feedback == 1 & coin_comp_long$fecund_feedback == 0 & 
                                               coin_comp_long$feed_feedback == 1 & coin_comp_long$spore_feedback == 1) | 
                                              (coin_comp_long$infect_feedback == 0 & coin_comp_long$fecund_feedback == 0 & 
                                                 coin_comp_long$feed_feedback == 0 & coin_comp_long$spore_feedback == 0)),]


ggplot(coin_comp_long_main,aes(y=M_prev,x=age,color = Treatment)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(colour = "black", size=1)) + 
  geom_boxplot() +
  #coord_cartesian(ylim = c(0,250)) +
  scale_colour_manual(values = c("blue", "dark green", "green")) +
  xlab("Age (Days)") + ylab("Prevalence") + ggtitle("No Fecund Feedback") +
  theme(axis.title.x = element_text(face="bold", size=16),axis.text.x  = element_text(size=16)) + 
  theme(axis.title.y = element_text(face="bold", size=16),axis.text.y  = element_text(size=16)) +
  facet_wrap(coin_comp_long_main$feed_feedback)

ggplot(coin_comp_long_main,aes(y=P_prev,x=age,color = Treatment)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(colour = "black", size=1)) + 
  geom_boxplot() +
  #coord_cartesian(ylim = c(0,250)) +
  scale_colour_manual(values = c("blue", "dark green", "green")) +
  xlab("Age (Days)") + ylab("Prevalence") + ggtitle("No Fecund Feedback") +
  theme(axis.title.x = element_text(face="bold", size=16),axis.text.x  = element_text(size=16)) + 
  theme(axis.title.y = element_text(face="bold", size=16),axis.text.y  = element_text(size=16)) +
  facet_wrap(coin_comp_long_main$feed_feedback)


##No spore_feedback

coin_comp_long_main <- coin_comp_long[which((coin_comp_long$infect_feedback == 1 & coin_comp_long$fecund_feedback == 1 & 
                                               coin_comp_long$feed_feedback == 1 & coin_comp_long$spore_feedback == 0) | 
                                              (coin_comp_long$infect_feedback == 0 & coin_comp_long$fecund_feedback == 0 & 
                                                 coin_comp_long$feed_feedback == 0 & coin_comp_long$spore_feedback == 0)),]


ggplot(coin_comp_long_main,aes(y=M_prev,x=age,color = Treatment)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(colour = "black", size=1)) + 
  geom_boxplot() +
  #coord_cartesian(ylim = c(0,250)) +
  scale_colour_manual(values = c("blue", "dark green", "green")) +
  xlab("Age (Days)") + ylab("Prevalence") + ggtitle("No Spore Feedback") +
  theme(axis.title.x = element_text(face="bold", size=16),axis.text.x  = element_text(size=16)) + 
  theme(axis.title.y = element_text(face="bold", size=16),axis.text.y  = element_text(size=16)) +
  facet_wrap(coin_comp_long_main$feed_feedback)

ggplot(coin_comp_long_main,aes(y=P_prev,x=age,color = Treatment)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(colour = "black", size=1)) + 
  geom_boxplot() +
  #coord_cartesian(ylim = c(0,250)) +
  scale_colour_manual(values = c("blue", "dark green", "green")) +
  xlab("Age (Days)") + ylab("Prevalence") + ggtitle("No Spore Feedback") +
  theme(axis.title.x = element_text(face="bold", size=16),axis.text.x  = element_text(size=16)) + 
  theme(axis.title.y = element_text(face="bold", size=16),axis.text.y  = element_text(size=16)) +
  facet_wrap(coin_comp_long_main$feed_feedback)


###################OK, Timing of fungal prevalence over time

coin_comp_long_main <- coin_comp_long[which((coin_comp_long$infect_feedback == 1 & coin_comp_long$fecund_feedback == 1 & 
                                               coin_comp_long$feed_feedback == 1 & coin_comp_long$spore_feedback == 1) | 
                                              (coin_comp_long$infect_feedback == 0 & coin_comp_long$fecund_feedback == 0 & 
                                                 coin_comp_long$feed_feedback == 0 & coin_comp_long$spore_feedback == 0)),]
#In order to make sure that the feedback model doesn't have three runs included, filter out 
coin_comp_long_main <- coin_comp_long_main[which(coin_comp_long_main$model != "Age Effects" | coin_comp_long_main$age == "20 days"),]
#I checked, they are all robust

ggplot(coin_comp_long_main,aes(y=M_prev,x=days,color = Treatment)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(colour = "black", size=1)) + 
  geom_line() +
  #coord_cartesian(ylim = c(0,250)) +
  scale_colour_manual(values = c("blue", "dark green", "green")) +
  xlab("Age (Days)") + ylab("Prevalence") + ggtitle("No Spore Feedback") +
  theme(axis.title.x = element_text(face="bold", size=16),axis.text.x  = element_text(size=16)) + 
  theme(axis.title.y = element_text(face="bold", size=16),axis.text.y  = element_text(size=16)) +
  facet_wrap(coin_comp_long_main$model)

OK, looks like that may be temporary.

