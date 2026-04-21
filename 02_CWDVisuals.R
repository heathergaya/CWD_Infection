##### Start Header ####
## Code for: Host-pathogen dynamics and the hidden hazards of chronic wasting disease
## Heather E. Gaya, Marcelo H. Jorge, Michael J. Chamberlain, Amy V. Nalls, 
## Nathaniel D. Denkers, Candace K. Mathiason, Mark G. Ruder, 
## Gino J. D’Angelo, and Richard B. Chandler
## Code written by Heather Gaya 
## Last modified April 21, 2026
################# End Header ####

#### Load libraries ####
library(coda)
library(nimble)
library(dplyr)
library(MCMCvis)
library(ggplot2)
library(ggpubr)

### Grab Output ####
nim.out <- readRDS('UpdatedSurvival_wSex_Jan25.rds')
nim.out2 <- as.mcmc.list(nim.out)[,c('B0', 'B1[1]','B1[2]', 'B2', 'B3', 'B4', 'B5', 'B6', 'p[1]', 'p[2]', 'p[3]', 'p1[1]', 'p1[2]', 'p2'),] #just the parameters we want
fawns2 <- read.csv('Predicted_Hazard_Survival.csv') #fawn survival from Jorge 2024 dissertation 

### Plot Figure 4 #####
### **Figure 4 b ####
surv_betas <- as.matrix(nim.out2[,c('B0', 'B2', 'B3', 'B4', 'B6'),])
surv_outs <- array(NA, c(nrow(surv_betas), 18))
surv_neg_F <- surv_neg_M <- surv_pos_F <- surv_pos_M <- array(NA, c(nrow(surv_betas), 3285)) #from .5 to 9.5 years old (would be 2920 if not using fawns)
sa <- (seq(1.5, 2.5, by = 1/365)[-366]-4.2473)/1.989805
a1 <- (seq(2.5, 3.5, by = 1/365)[-366]-4.2473)/1.989805
a2 <- (seq(3.5, 4.5, by = 1/365)[-366]-4.2473)/1.989805
a3 <- (seq(4.5, 5.5, by = 1/365)[-366]-4.2473)/1.989805
a4 <- (seq(5.5, 6.5, by = 1/365)[-366]-4.2473)/1.989805
a5 <- (seq(6.5, 7.5, by = 1/365)[-366]-4.2473)/1.989805

#avg roughness is 0.5475648 for the deer in our study
di <- seq(1, 365, by = 1)/100
#for long curve, they get infected on day 365
for(k in 1:nrow(surv_betas)){
  surv_outs[k,1] <- prod(plogis(surv_betas[k,'B0'] + surv_betas[k, 'B2']*sa + surv_betas[k,'B3']*0 + surv_betas[k,'B6']*0 + surv_betas[k,'B4']*0.5475648)) #neg F SA
  surv_outs[k,2] <- prod(plogis(surv_betas[k,1] + surv_betas[k, 'B2']*a1 + surv_betas[k,'B3']*0 + surv_betas[k,'B6']*0+ surv_betas[k,'B4']*0.5475648)) #neg F A
  surv_outs[k,3] <- prod(plogis(surv_betas[k,1] + surv_betas[k, 'B2']*sa + surv_betas[k,'B3']*0 + surv_betas[k,'B6']*di+ surv_betas[k,'B4']*0.5475648)) #Pos F SA 1yr infect
  surv_outs[k,4] <- prod(plogis(surv_betas[k,1] + surv_betas[k, 'B2']*a1 + surv_betas[k,'B3']*0 + surv_betas[k,'B6']*(di+ 3.65)+ surv_betas[k,'B4']*0.5475648)) #Pos F A 2yr infect
  surv_outs[k,5] <- prod(plogis(surv_betas[k,1] + surv_betas[k, 'B2']*sa + surv_betas[k,'B3']*1 + surv_betas[k,'B6']*0+ surv_betas[k,'B4']*0.5475648)) #neg M SA
  surv_outs[k,6] <- prod(plogis(surv_betas[k,1] + surv_betas[k, 'B2']*a1 + surv_betas[k,'B3']*1 + surv_betas[k,'B6']*0+ surv_betas[k,'B4']*0.5475648)) #neg M A
  surv_outs[k,7] <- prod(plogis(surv_betas[k,1] + surv_betas[k, 'B2']*sa + surv_betas[k,'B3']*1 + surv_betas[k,'B6']*di+ surv_betas[k,'B4']*0.5475648)) #Pos M SA 1yr infect
  surv_outs[k,8] <- prod(plogis(surv_betas[k,1] + surv_betas[k, 'B2']*a1 + surv_betas[k,'B3']*1 + surv_betas[k,'B6']*(di+3.65)+ surv_betas[k,'B4']*0.5475648)) #Pos M A 2yr infect
  
  surv_outs[k,9] <- prod(plogis(surv_betas[k,1] + surv_betas[k, 'B2']*a2 + surv_betas[k,'B3']*0 + surv_betas[k,'B6']*(di+2*3.65)+ surv_betas[k,'B4']*0.5475648)) #Pos F A 3 yr infect
  surv_outs[k,10] <- prod(plogis(surv_betas[k,1] + surv_betas[k, 'B2']*a2 + surv_betas[k,'B3']*1 + surv_betas[k,'B6']*(di+2*3.65)+ surv_betas[k,'B4']*0.5475648)) #Pos M A 3 yr infect
  
  surv_outs[k,13] <- prod(plogis(surv_betas[k,1] + surv_betas[k, 'B2']*a3 + surv_betas[k,'B3']*0 + surv_betas[k,'B6']*(di+3*3.65)+ surv_betas[k,'B4']*0.5475648)) #Pos F A 4 yr infect
  surv_outs[k,14] <- prod(plogis(surv_betas[k,1] + surv_betas[k, 'B2']*a3 + surv_betas[k,'B3']*1 + surv_betas[k,'B6']*(di+3*3.65)+ surv_betas[k,'B4']*0.5475648)) #Pos M A 4 yr infect
  
  surv_outs[k,15] <- prod(plogis(surv_betas[k,1] + surv_betas[k, 'B2']*a4 + surv_betas[k,'B3']*0 + surv_betas[k,'B6']*(di+4*3.65)+ surv_betas[k,'B4']*0.5475648)) #Pos F A 5 yr infect
  surv_outs[k,16] <- prod(plogis(surv_betas[k,1] + surv_betas[k, 'B2']*a4 + surv_betas[k,'B3']*1 + surv_betas[k,'B6']*(di+4*3.65)+ surv_betas[k,'B4']*0.5475648)) #Pos M A 5 yr infect
  
  surv_outs[k,17] <- prod(plogis(surv_betas[k,1] + surv_betas[k, 'B2']*a5 + surv_betas[k,'B3']*0 + surv_betas[k,'B6']*(di+5*3.65)+ surv_betas[k,'B4']*0.5475648)) #Pos F A 6 yr infect
  surv_outs[k,18] <- prod(plogis(surv_betas[k,1] + surv_betas[k, 'B2']*a5 + surv_betas[k,'B3']*1 + surv_betas[k,'B6']*(di+5*3.65)+ surv_betas[k,'B4']*0.5475648)) #Pos M A 6 yr infect
  
  
  #survival for subadults, all 4 options, just day 1.
  #just the mean for fawns right now, might not work for graphing:
  for(tt in 1:365){
    s <- rnorm(1, mean = fawns2[tt, 2], sd = (fawns2[tt, 4]-fawns2[tt, 3])/3.92)
    s <- min(s, 1)
    surv_neg_F[k,tt] <- s
    surv_neg_M[k,tt] <- s
    surv_pos_F[k,tt] <- s
    surv_pos_M[k,tt] <- s
  }
  
  surv_outs[k,11] <- surv_neg_F[k,365] #fawns
  surv_outs[k,12] <- surv_neg_F[k,365] #fawns
  
  surv_neg_F[k,366] <-  plogis(surv_betas[k,1]  + surv_betas[k,'B3']*0 + surv_betas[k,'B6']*0+ surv_betas[k,'B4']*0.5475648)[1]*surv_neg_F[k,365] #+ surv_betas[k, 'B2']*sa
  
  surv_neg_M[k,366] <-  plogis(surv_betas[k,1] + surv_betas[k,'B3']*1 + surv_betas[k,'B6']*0+ surv_betas[k,'B4']*0.5475648)[1]*surv_neg_M[k,365]
  surv_pos_F[k,366] <-  plogis(surv_betas[k,1] + surv_betas[k,'B3']*0 + surv_betas[k,'B6']*1/100 + surv_betas[k,'B4']*0.5475648)[1]*surv_pos_F[k,365]
  surv_pos_M[k,366] <-  plogis(surv_betas[k,1] + surv_betas[k,'B3']*1 + surv_betas[k,'B6']*1/100 + surv_betas[k,'B4']*0.5475648)[1]*surv_pos_M[k,365]
}

survs <- apply(surv_outs, 2, FUN = function(x){quantile(x, c(0.025, .5, .975))})

surv_gg <- data.frame(Median = survs[2,],
                      LCI = survs[1,],
                      UCI = survs[3,],
                      Stage = c('- F SA', '- F A', '+ F SA', '+ F A',
                                '- M SA', '- M A', '+ M SA', '+ M A',
                                '+ F 2A', '+ M 2A', 'F Fawn','M Fawn', 
                                '+ F 3A', '+ M 3A','+ F 4A', '+ M 4A',
                                '+ F 5A', '+ M 5A'),
                      Status = c(rep(rep(c('Neg', 'Pos'), each = 2), 2), '2 yrs Pos','2 yrs Pos', 
                                 'Fawn', 'Fawn',  '3 yrs Pos','3 yrs Pos','4 yrs Pos','4 yrs Pos', '5 yrs Pos','5 yrs Pos'),
                      Sex = c(rep(c('Female', 'Male'), each = 4), 'Female', 
                              'Male', 'Female', 'Male', 'Female', 
                              'Male','Female', 'Male','Female', 'Male'))
surv_gg$Stage <- factor(surv_gg$Stage, 
                        levels = c('F Fawn', '- F SA', '+ F SA', '- F A', '+ F A',
                                   '+ F 2A','+ F 3A','+ F 4A','+ F 5A',
                                   'M Fawn', '- M SA', '+ M SA', '- M A', '+ M A',
                                   '+ M 2A','+ M 3A','+ M 4A','+ M 5A'))
surv_gg2 <- subset(surv_gg, surv_gg$Stage != '- M A' & surv_gg$Stage != '- F A')
surv_gg2$Stage2 <- c(rep(c('Negative', '+ 1', '+ 2'), 2), rep(c('+ 3',  '0', '+ 4', '+ 5', '+ 6'), each = 2))

aa2 <- ggplot(surv_gg2, aes(x = Stage, y = Median))+
  geom_point(aes(col = Stage2), alpha= .8, cex= 3)+
  geom_errorbar(aes(ymin = LCI, ymax = UCI), width = .2)+
  theme_linedraw()+
  facet_wrap(~Sex, scales = 'free_x')+
  scale_color_manual(values = c( '#ea8367','#d56b56', '#bf5345', '#a93c34', '#942223', '#900003', 'grey','grey40', 'skyblue'))+
  scale_x_discrete(
    labels = c(
      'F Fawn' = 'Fawn',
      '- F SA' = '0',
      '+ F SA' = ' 1',
      '+ F A'  = ' 2',
      '+ F 2A' = ' 3',
      '+ F 3A' = ' 4',
      '+ F 4A' = ' 5',
      '+ F 5A' = ' 6',
      
      'M Fawn' = 'Fawn',
      '- M SA' = '0',
      '+ M SA' = ' 1',
      '+ M A'  = ' 2',
      '+ M 2A' = ' 3',
      '+ M 3A' = ' 4',
      '+ M 4A' = ' 5',
      '+ M 5A' = ' 6'
    )
  ) +
  ylab('Annual survival')+
  xlab('Years since infection')+
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 22),
        legend.text = element_text(size = 18),
        legend.title = element_text(size = 20),
        strip.text = element_text(size = 25), 
        strip.background = element_rect(fill = 'grey50'),
        panel.grid.minor = element_blank(),
        legend.position = 'none',
        legend.position.inside =  c(.9, .8))

aa2

### **Figure 4 a ####
### Survival curves. These loops take several minutes to run. 
## add in fawns 
Cp_PosM <- Cp_PosF <- Cp_NegF <- Cp_NegM <- Cp_LatePosM <- Cp_LatePosF <- array(NA, c(nrow(surv_betas), 33)) #for CIs for conditional lifetimes
surv_latePos_F <- surv_latePos_M <- array(NA, dim(surv_neg_F))
for(k in 1:nrow(surv_betas)){
  for(m in 2:2920){ #to age 9.5
    surv_neg_F[k,(m+365)] <-  plogis(surv_betas[k,1] + surv_betas[k, 'B2']*((1.5+m/365)-4.271508)/1.993448 + surv_betas[k,'B3']*0 + surv_betas[k,'B6']*0+ surv_betas[k,'B4']*0.5475648)*surv_neg_F[k,(m+365)-1]
    surv_neg_M[k,(m+365)] <-  plogis(surv_betas[k,1] + surv_betas[k, 'B2']*((1.5+m/365)-4.271508)/1.993448 + surv_betas[k,'B3']*1 + surv_betas[k,'B6']*0+ surv_betas[k,'B4']*0.5475648)*surv_neg_M[k,(m+365)-1]
    surv_pos_F[k,(m+365)] <-  plogis(surv_betas[k,1] + surv_betas[k, 'B2']*((1.5+m/365)-4.271508)/1.993448 + surv_betas[k,'B3']*0 + surv_betas[k,'B6']*(1+m)/100 + surv_betas[k,'B4']*0.5475648)*surv_pos_F[k,(m+365)-1]
    surv_pos_M[k,(m+365)] <-  plogis(surv_betas[k,1] + surv_betas[k, 'B2']*((1.5+m/365)-4.271508)/1.993448 + surv_betas[k,'B3']*1 + surv_betas[k,'B6']*(1+m)/100 + surv_betas[k,'B4']*0.5475648)*surv_pos_M[k,(m+365)-1]
  }
  
  ## What about a F and a M that don't get infected until age 3? 
  surv_latePos_F[k,1:1095] <- surv_neg_F[k,1:1095]
  surv_latePos_M[k,1:1095] <- surv_neg_M[k,1:1095]
  for(m in 731:2920){
    surv_latePos_F[k,(m+365)] <-  plogis(surv_betas[k,1] + surv_betas[k, 'B2']*((1.5+m/365)-4.271508)/1.993448 + 
                                           surv_betas[k,'B3']*0 + surv_betas[k,'B6']*(1+(m-365))/100 + 
                                           surv_betas[k,'B4']*0.5475648)*surv_latePos_F[k,(m+365)-1]
    surv_latePos_M[k,(m+365)] <-  plogis(surv_betas[k,1] + surv_betas[k, 'B2']*((1.5+m/365)-4.271508)/1.993448 + 
                                           surv_betas[k,'B3']*1 + surv_betas[k,'B6']*(1+(m-365))/100 + 
                                           surv_betas[k,'B4']*0.5475648)*surv_latePos_M[k,(m+365)-1]
  }
  
  
  
  Cp_NegF[k,] <- surv_neg_F[k,seq(365, 3285, by = 365/4)]/sum(surv_neg_F[k,seq(365, 3285, by = 365/4)])
  Cp_PosF[k,] <- surv_pos_F[k,seq(365, 3285, by = 365/4)]/sum(surv_pos_F[k,seq(365, 3285, by = 365/4)])
  Cp_NegM[k,] <- surv_neg_M[k,seq(365, 3285, by = 365/4)]/sum(surv_neg_M[k,seq(365, 3285, by = 365/4)])
  Cp_PosM[k,] <- surv_pos_M[k,seq(365, 3285, by = 365/4)]/sum(surv_pos_M[k,seq(365, 3285, by = 365/4)])
  Cp_LatePosM[k,] <- surv_latePos_M[k,seq(365, 3285, by = 365/4)]/sum(surv_latePos_M[k,seq(365, 3285, by = 365/4)])
  Cp_LatePosF[k,] <- surv_latePos_F[k,seq(365, 3285, by = 365/4)]/sum(surv_latePos_F[k,seq(365, 3285, by = 365/4)])
}

survs_NF <- apply(surv_neg_F, 2, FUN = function(x){quantile(x, c(0.025, .5, .975))})
survs_NM <- apply(surv_neg_M, 2, FUN = function(x){quantile(x, c(0.025, .5, .975))})
survs_PF <- apply(surv_pos_F, 2, FUN = function(x){quantile(x, c(0.025, .5, .975))})
survs_PM <- apply(surv_pos_M, 2, FUN = function(x){quantile(x, c(0.025, .5, .975))})

survs_LPM <- apply(surv_latePos_M, 2, FUN = function(x){quantile(x, c(0.025, .5, .975))})
survs_LPF <- apply(surv_latePos_F, 2, FUN = function(x){quantile(x, c(0.025, .5, .975))})

survs_cpNF <- apply(Cp_NegF, 2, FUN = function(x){quantile(x, c(0.025, .5, .975))})
survs_cpNM <- apply(Cp_NegM, 2, FUN = function(x){quantile(x, c(0.025, .5, .975))})
survs_cpPF <- apply(Cp_PosF, 2, FUN = function(x){quantile(x, c(0.025, .5, .975))})
survs_cpPM <- apply(Cp_PosM, 2, FUN = function(x){quantile(x, c(0.025, .5, .975))})
survs_cpLPM <- apply(Cp_LatePosM, 2, FUN = function(x){quantile(x, c(0.025, .5, .975))})
survs_cpLPF <- apply(Cp_LatePosF, 2, FUN = function(x){quantile(x, c(0.025, .5, .975))})

surv_lines <- data.frame(Median = c(survs_NF[2,], survs_NM[2,],survs_PF[2,], survs_PM[2,], survs_LPF[2,], survs_LPM[2,]),
                         LCI = c(survs_NF[1,], survs_NM[1,],survs_PF[1,], survs_PM[1,], survs_LPF[1,], survs_LPM[1,]),
                         UCI = c(survs_NF[3,], survs_NM[3,],survs_PF[3,], survs_PM[3,], survs_LPF[3,], survs_LPM[3,]),
                         Stage = rep(c("Neg F", "Neg M", "Pos F", "Pos M", 'Late Pos F', 'Late Pos M'), each = ncol(survs_NF)),
                         Status = c(rep(c('CWD-', 'CWD+ @1yr'), each = ncol(survs_NF)*2), rep('CWD+ @3yr', ncol(survs_NF)*2)),
                         Sex = rep(c('Female', 'Male'), each = ncol(survs_NF)),
                         Age = rep(c(c(183:547)/365, 1.5+c(1:2920)/365), 6))

surv_4bb2 <- surv_lines[surv_lines$Age <= 7.5,]
surv_4bb2$lnUCI <- -log(surv_4bb2$UCI)
surv_4bb2$lnMed <- -log(surv_4bb2$Median)
surv_4bb2$lnLCI <- -log(surv_4bb2$LCI)

bb2 <- 
  ggplot(surv_4bb2, aes(x = Age, y = lnMed))+
  geom_line(aes(col = Status))+
  geom_ribbon(aes(ymin = lnLCI, ymax = lnUCI, fill = Status), alpha = .2)+
  theme_linedraw()+
  facet_wrap(~Sex, ncol = 2)+
  xlim(.5, 7.5)+
  scale_y_continuous(
    labels = scales::number_format(accuracy = 0.1))+
  ylab('Cummulative hazard ')+
  xlab('Age (years)')+
  scale_fill_manual(values = c('skyblue', 'red','purple'))+
  scale_color_manual(values = c('navyblue', 'firebrick3','purple'))+
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 22),
        legend.text = element_text(size = 18),
        legend.title = element_text(size = 20),
        strip.text = element_text(size = 25), 
        strip.background = element_rect(fill = 'grey50'),
        panel.grid.minor = element_blank(),
        legend.position = 'none',
        legend.position.inside =  c(.15, .25))

ggarrange(bb2, aa2, nrow = 2)

### **Figure 4 c ####
### Thinking about expected lifetimes. Assume they make it past the fawn stage:

ggLifeLine <- data.frame(Median_cp = c(survs_cpNF[2,],survs_cpPF[2,], survs_cpNM[2,],survs_cpPM[2,], survs_cpLPF[2,], survs_cpLPM[2,]),
                         LCI_cp = c(survs_cpNF[1,],survs_cpPF[1,], survs_cpNM[1,],survs_cpPM[1,], survs_cpLPF[1,], survs_cpLPM[1,]),
                         UCI_cp = c(survs_cpNF[3,],survs_cpPF[3,], survs_cpNM[3,],survs_cpPM[3,], survs_cpLPF[3,], survs_cpLPM[3,]),
                         Sex = c(rep(c('Female', 'Male'), each = ncol(survs_cpNF)*2), rep(c('Female', 'Male'), each = ncol(survs_cpNF))),
                         Age = seq(1.5, 9.5, by = .25),
                         Status = c(rep(rep(c('CWD-', 'CWD+ @1yr'), each = ncol(survs_cpNF)), 2), rep('CWD+ @3yr', ncol(survs_cpNF)*2)))

cc <- ggplot(ggLifeLine, aes(x = Age, y = Median_cp, fill = Status, col = Status))+
  facet_wrap(~Sex)+
  geom_ribbon(aes(ymin = LCI_cp, ymax = UCI_cp), alpha = .2, lwd = 0)+
  geom_line()+
  geom_vline(xintercept = 1.5, lty = 2)+
  scale_fill_manual(values = c('skyblue', 'red', 'purple'))+
  scale_color_manual(values = c('navyblue', 'firebrick3', 'purple'))+
  theme_linedraw()+
  xlab('Lifetime')+
  ylab('Probabiltiy density')+
  xlim(.5, 7.5)+
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 22),
        legend.text = element_text(size = 18),
        legend.title = element_text(size = 20),
        strip.text = element_text(size = 25), 
        strip.background = element_rect(fill = 'grey50'), 
        panel.grid.minor = element_blank(),
        legend.position = 'inside',
        legend.position.inside =  c(.85, .7))

cc
### **Figure 4 full ####
ggarrange(bb2, aa2,  cc, nrow = 3, labels = 'AUTO')

### Plot Figure 2 ####
#hazard[i,t] <- p2*(1-exp(-p1[sex]*Prion_dens[i,t-1]))/365 #this is the hazard equation in the model
nim_mat <- as.matrix(as.mcmc.list(nim.out)[,c('p1[1]', 'p1[2]', 'p2')]) 
infs_now <- seq(0, 4, by = .05) #total infected deer density in time t-1
infs_before <- c(0, 2, 4) #total infected deer density in time t-2
infect_pm <- infect_pf <- array(NA, c(nrow(nim_mat), length(infs_now), length(infs_before)))
for(j in 1:nrow(nim_mat)){
  for(q in 1:length(infs_before)){
    for(pp in 1:length(infs_now)){
      infect_pm[j,pp,q] <- nim_mat[j,'p2']*(1-exp(-(nim_mat[j,'p1[2]']*(.75*infs_now[pp]+.25*infs_before[q]))))
      infect_pf[j,pp,q] <- nim_mat[j,'p2']*(1-exp(-(nim_mat[j,'p1[1]']*(.75*infs_now[pp]+.25*infs_before[q]))))
    }
  }}

inf_outf <- apply(infect_pf, c(2,3), function(x){quantile(x, c(0.025,.5, .975))})
inf_outm <- apply(infect_pm, c(2,3), function(x){quantile(x, c(0.025,.5, .975))})

gg_df <- data.frame(
  median = c(inf_outm[2,,], inf_outf[2,,]),
  LCI = c(inf_outm[1,,], inf_outf[1,,]),
  UCI = c(inf_outm[3,,], inf_outf[3,,]),
  DeerDensity = rep(rep(infs_now, length(infs_before)), 2),
  OldDens = paste0("I[i*','*t-1] == ", 
                   rep(as.character(rep(infs_before, each = length(infs_now))), 2)),
  Sex = rep(c('Male', 'Female'), each = length(inf_outm[2,,]))
)

ggplot(gg_df, aes(x = DeerDensity))+ 
  geom_ribbon(aes(ymin = LCI, ymax = UCI), lty = 2, fill = 'grey90')+
  geom_line(aes(y = median), lwd = 1.2)+
  theme_linedraw()+
  xlab(expression("Infected deer density (" * I[i*","*t] * ")"))+
  ylab('Annual probability of infection')+
  facet_grid(rows = vars(Sex), cols = vars(OldDens), labeller = label_parsed)+
  scale_shape_manual(values = c(23, 24,22))+
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 22),
        legend.text = element_text(size = 18),
        legend.title = element_text(size = 20),
        strip.text = element_text(size = 25), 
        strip.background = element_rect(fill = 'grey50'),
        panel.grid = element_blank(),
        legend.direction = 'horizontal', legend.position = 'bottom')

## for average infection probs for males/females:
stuff <- rbind(
  gg_df[gg_df$Sex == 'Female' & gg_df$DeerDensity == 0.75,][2,], #from 2 to 0.75 
  gg_df[gg_df$Sex == 'Male' & gg_df$DeerDensity == 0.75,][2,], #from 2 to 0.75 
  gg_df[gg_df$Sex == 'Female' & gg_df$DeerDensity == 2.5,][3,], #from 4 to 2.5
  gg_df[gg_df$Sex == 'Male' & gg_df$DeerDensity == 2.5,][3,] #from 4 to 2.5 
)
kableExtra::kable(stuff, digits = 2)

#for abstract:
stuff2 <- rbind(
  gg_df[gg_df$Sex == 'Female' & gg_df$DeerDensity == 1,][3,], 
  gg_df[gg_df$Sex == 'Female' & gg_df$DeerDensity == 4,][3,], 
  gg_df[gg_df$Sex == 'Male' & gg_df$DeerDensity == 1,][3,], 
  gg_df[gg_df$Sex == 'Male' & gg_df$DeerDensity == 4,][3,] 
)
kableExtra::kable(stuff2, digits = 2)

# Plot Figure 3 ####
BigOut <- readRDS('UpdatedSurvival_LargePosterior.rds') #large posterior where we monitored extra params
mydays <- paste0('daysinfected[133, ', 1:697, ']')
mydays <- mydays[-c(240:330)]
mysurvs <- paste0('me.phi[133, ', 1:697, ']')
mysurvs <- mysurvs[-c(240:330)]
plotsurvs <- MCMCsummary(as.mcmc.list(BigOut)[,mysurvs,])
plotsurvs$Day <- c(1:606)

plotme <- as.matrix(as.mcmc.list(BigOut)[,mydays,])
firstinfect <- array(NA, nrow(plotme))
for(j in 1:nrow(plotme)){
  firstinfect[j] <- which(plotme[j,] == 1)
}
infectday <- quantile(firstinfect, c(0.025, .5, .975))
infectday #when was this animal infected

library(grid)
rect1 <- rectGrob(
  gp = gpar(fill = '#800080', alpha = .05, col = NA)
)
rect2 <- rectGrob(
  gp = gpar(fill = "#a93c34", alpha = 0.2, col = NA)
)

### **Figure 3 A.1 ####
(g1 <- ggplot()+
   geom_line(data = plotsurvs[plotsurvs$Day < 479,],aes(x = Day, y  = `50%`), lwd = .75, col = 'darkblue')+
   geom_ribbon(data = plotsurvs[plotsurvs$Day < 479,],aes(x = Day, ymin = `2.5%`, ymax = `97.5%`), col = '#009AEE', alpha = .2, fill = 'skyblue', lwd = .2)+
   geom_line(data = plotsurvs[plotsurvs$Day > 477 & plotsurvs$Day < 563,],aes(x = Day, y  = `50%`), lwd = .75, col = '#800080')+
   geom_ribbon(data = plotsurvs[plotsurvs$Day > 477 & plotsurvs$Day < 563,],aes(x = Day, ymin = `2.5%`, ymax = `97.5%`), col = '#800080', alpha = .05, fill = '#800080', lwd = .2)+
   geom_line(data = plotsurvs[plotsurvs$Day > 561,],aes(x = Day, y  = `50%`), lwd = .75, col = '#a93c34')+
   geom_ribbon(data = plotsurvs[plotsurvs$Day > 561,],aes(x = Day, ymin = `2.5%`, ymax = `97.5%`), col = '#a93c34', alpha = .2, fill = '#a93c34', lwd = .2)+
   geom_vline(xintercept = infectday[3], lty = 2, col = 'red', lwd = .5)+ 
   geom_vline(xintercept = 331, lty =2, col = 'blue', lwd = .5)+
   geom_vline(xintercept = 141, lty =2, col = 'blue', lwd = .5)+
   theme_bw()+
   ylab('Daily survival probability')+
   xlab('Day')+
   ylim(0.967, 1)+
   scale_x_continuous(limits = c(0, 606), breaks = c(0, 100, 300, 500), expand = c(.005, 0))+
   theme(axis.title = element_text(size = 15),
         axis.text = element_text(size = 15),
         plot.title = element_text(size  = 20))+
   ggtitle('Deer 133')
)

### **Figure 3 A.2 ####
plotOuts <- MCMCchains(as.mcmc.list(BigOut), params = mysurvs, exact= T, ISB = F)
cumprod_surv_all <- array(NA, c(nrow(plotOuts), length(mysurvs)))
for(j in 1:nrow(plotOuts)){
  cumprod_surv_all[j,] <- cumprod(plotOuts[j,])
}

cumprod_CIs <- apply(cumprod_surv_all, 2, function(x){quantile(x, c(0.025, .5, .975))})
cumProd_plotsurvs <- data.frame(Median = cumprod_CIs[2,],
                                LCI = cumprod_CIs[1,],
                                UCI = cumprod_CIs[3,],
                                Days= plotsurvs$Day)
(g2 <- ggplot(cumProd_plotsurvs, aes(x = Days))+
    theme_bw()+
    geom_line(data = cumProd_plotsurvs[cumProd_plotsurvs$Day < 479,],aes(x = Days, y  = Median), lwd = .75, col = 'darkblue')+
    geom_ribbon(data = cumProd_plotsurvs[cumProd_plotsurvs$Day < 479,],aes(x = Days, ymin = LCI, ymax = UCI), col = '#009AEE', alpha = .2, fill = 'skyblue', lwd = .2)+
    geom_line(data = cumProd_plotsurvs[cumProd_plotsurvs$Days > 477 & cumProd_plotsurvs$Days < 563,],aes(x = Days, y  = Median), lwd = .75, col = '#800080')+
    geom_ribbon(data = cumProd_plotsurvs[cumProd_plotsurvs$Days > 477 & cumProd_plotsurvs$Days < 563,],aes(x = Days, ymin = LCI, ymax = UCI), col = '#800080', alpha = .05, fill = '#800080', lwd = .2)+
    geom_line(data = cumProd_plotsurvs[cumProd_plotsurvs$Days > 561,],aes(x = Days, y  = Median), lwd = .75, col = '#a93c34')+
    geom_ribbon(data = cumProd_plotsurvs[cumProd_plotsurvs$Days > 561,],aes(x = Days, ymin = LCI, ymax = UCI),  col = '#a93c34', alpha = .2, fill = '#a93c34', lwd = .2)+
    geom_vline(xintercept = infectday[3], lty = 2, col = 'red', lwd = .5)+ 
    geom_vline(xintercept = 331, lty =2, col = 'blue', lwd = .5)+
    geom_vline(xintercept = 141, lty =2, col = 'blue', lwd = .5)+
    ylab('Survivorship')+
    xlab('Day')+
    ylim(0, 1)+
    theme(axis.title = element_text(size = 15),
          axis.text = element_text(size = 15))+
    scale_x_continuous(limits = c(0, 606), breaks = c(0, 100, 300, 500), expand = c(.005, 0))
)

### **Figure 3 A.3 ####

Temp_133 <-  data.frame(Temperature_s = nimstuff$consts$TAvg[133,c(1:697)[-c(240:330)]],
                        Days= plotsurvs$Day)
#mean 59.94255; sd 16.0104
#scaled*sd(x)+mean(x) = x
Temp_133$Temp <- Temp_133$Temperature_s*16.0104 + 59.94255
Temp_133$TempC <- (Temp_133$Temp - 32)*(5/9)

(g3 <- ggplot(Temp_133, aes(x = Days))+
    annotation_custom(rect1,
                      xmin = infectday[1], xmax = infectday[3], ymin = -Inf, ymax = Inf) +
    annotation_custom(rect2,
                      xmin = infectday[3], xmax = 606, ymin = -Inf, ymax = Inf) +
    geom_line(aes(y  = TempC), lwd = .75)+
    geom_vline(xintercept = infectday[3], lty = 2, col = 'red', lwd = .5)+ 
    geom_vline(xintercept = 331, lty =2, col = 'blue', lwd = .5)+
    geom_vline(xintercept = 141, lty =2, col = 'blue', lwd = .5)+
    theme_bw()+
    ylab('Temperature (°C)')+
    xlab('Day')+
    ylim(-12.3, 31.6)+
    scale_x_continuous(limits = c(0, 606), breaks = c(0, 100, 300, 500), expand = c(.005, 0))+
    theme(axis.title = element_text(size = 15),
          axis.text = element_text(size = 15))
)


ggarrange(g1, g2, g3, nrow = 3)

### **Figure 3 B.1 ####
mydays2 <- paste0('daysinfected[31, ', 1:768, ']')
mysurvs2 <- paste0('me.phi[31, ', 1:768, ']')
plotsurvs2 <- MCMCsummary(as.mcmc.list(BigOut)[,mysurvs2,])
plotsurvs2$Day <- 1:768

plotme2 <- as.matrix(as.mcmc.list(BigOut)[,mydays2,])
(gg1 <- ggplot(plotsurvs2, aes(x = Day))+
    geom_line(aes(y  = `50%`), lwd = .75, col ='#900003')+
    geom_ribbon(aes(ymin = `2.5%`, ymax = `97.5%`), col = '#a93c34', alpha = .2, fill = '#a93c34', lwd = .1)+
    geom_vline(xintercept = 1, lwd = .5,lty =2, col = 'red')+
    geom_vline(xintercept = 356, lwd = .5, lty =2, col = 'red')+
    geom_vline(xintercept = 769, lwd = .5, lty =2, col = 'red')+
    theme_bw()+
    ylab('')+
    xlab('Day')+
    ylim(0.967, 1)+
    ggtitle('Deer 31')+
    scale_x_continuous(limits = c(0, 769), expand = c(.005, 0))+
    theme(axis.title = element_text(size = 15),
          plot.title = element_text(size  = 20),
          axis.text = element_text(size = 15),
          axis.text.y = element_blank())
)

### **Figure 3 B.2 ####
plotOuts2 <- MCMCchains(as.mcmc.list(BigOut), params = mysurvs2, exact= T, ISB = F)
cumprod_surv_all2 <- array(NA, c(nrow(plotOuts2), length(mysurvs2)))
for(j in 1:nrow(plotOuts2)){
  cumprod_surv_all2[j,] <- cumprod(plotOuts2[j,])
}

cumprod_CIs2 <- apply(cumprod_surv_all2, 2, function(x){quantile(x, c(0.025, .5, .975))})
cumProd_plotsurvs2 <- data.frame(Median = cumprod_CIs2[2,],
                                 LCI = cumprod_CIs2[1,],
                                 UCI = cumprod_CIs2[3,],
                                 Days= plotsurvs2$Day)
(gg2 <- ggplot(cumProd_plotsurvs2, aes(x = Days))+
    geom_line(aes(y  = Median), lwd = .75, col ='#900003')+
    geom_ribbon(aes(ymin = LCI, ymax = UCI), col = '#a93c34', alpha = .2, fill = '#a93c34', lwd = .2)+
    geom_vline(xintercept = 1, lty =2, lwd = .5, col = 'red')+
    geom_vline(xintercept = 356, lty =2, lwd = .5,col = 'red')+
    geom_vline(xintercept = 769, lty =2, lwd = .5,col = 'red')+
    theme_bw()+
    xlab('Day')+
    ylim(0, 1)+
    ylab('')+
    scale_x_continuous(limits = c(0, 769), expand = c(.005, 0))+
    theme(axis.title = element_text(size = 15),
          axis.text = element_text(size = 15),
          axis.text.y = element_blank())
)

### **Figure 3 B.3 ####
Temp_31 <-  data.frame(Temperature_s = nimstuff$consts$TAvg[31,c(1:768)],
                       Days= plotsurvs2$Day)
#mean 59.94255; sd 16.0104
#scaled*sd(x)+mean(x) = x
Temp_31$Temp <- Temp_31$Temperature_s*16.0104 + 59.94255
Temp_31$TempC <- (Temp_31$Temp - 32)*(5/9)

(gg3 <- ggplot(Temp_31, aes(x = Days))+
    geom_line(aes(y = TempC), lwd = .75)+
    geom_vline(xintercept = 1, lty =2,lwd = .5, col = 'red')+
    geom_vline(xintercept = 356, lty =2,lwd = .5, col = 'red')+
    geom_vline(xintercept = 769, lty =2, lwd = .5,col = 'red')+
    theme_bw()+
    #ylab('Temperature (°C)')+
    xlab('Day')+
    ylim(-12.3, 31.6)+
    theme(axis.title = element_text(size = 15),
          axis.text = element_text(size = 15),
          axis.text.y = element_blank())+
    ylab('')+
    scale_x_continuous(limits = c(0, 769), expand = c(.005, 0))
)


ggarrange(gg1, gg2, gg3, nrow = 3)

### **Figure 3 full ####
ggarrange(ggarrange(g1, g2, g3, nrow = 3), ggarrange(gg1, gg2, gg3, nrow = 3), nrow = 1, labels = 'AUTO')


#### Plot Figure 1 ####
quick_prev <- read.csv('Prevalence_Data.csv')
table(table(quick_prev$Sample.ID)) #note this data also includes deer without GPS collars
quick_prev$age_nn <- ifelse(quick_prev$Age %in% c('8 mo', '8mo', '8mon', '0.5', '6mon'), 1.5,
                            ifelse(quick_prev$Age %in% c('6.5+', '6.5', '5.5', '4.5',"Ad 7.5","Ad 5.5","Ad 4.5","Ad 8.5","Ad 9.5"), 4.5,
                                   ifelse(quick_prev$Age == "Ad 2.5", 2.5, 
                                          ifelse(quick_prev$Age == "Ad 3.5", 3.5,
                                                 as.numeric(quick_prev$Age)))))
pre_age <- quick_prev %>% group_by(Sex, age_nn, Status, year) %>% summarize(Prevalence = n() )

deer_prop_s <- pre_age %>%
  group_by(Sex, age_nn,year) %>%
  summarise(
    total_deer = sum(Prevalence),  # Total number of deer for each combo
    pos_deer = sum(Prevalence[Status == 'Pos'], na.rm = T),  # Number of 'Pos' deer
    proportion_pos = pos_deer / total_deer  # prop of 'Pos'
  ) %>%
  ungroup()

deer_prop_A1 <- pre_age %>%
  group_by(Sex, age_nn) %>%
  summarise(
    total_deer = sum(Prevalence),  # Total number of deer for each combo
    pos_deer = sum(Prevalence[Status == 'Pos'], na.rm = T),  # Number of 'Pos' deer
    proportion_pos = pos_deer / total_deer  # prop of 'Pos'
  ) %>%
  ungroup() %>%
  mutate(year = 1111)

deer_prop_A2 <- pre_age %>%
  group_by(age_nn) %>%
  summarise(
    total_deer = sum(Prevalence),  # Total number of deer for each combo
    pos_deer = sum(Prevalence[Status == 'Pos'], na.rm = T),  # Number of 'Pos' deer
    proportion_pos = pos_deer / total_deer  # prop of 'Pos'
  ) %>%
  ungroup()%>%
  mutate(Sex = "All", year = 1111)

deer_prop_all <- pre_age %>%
  group_by(age_nn, year) %>%
  summarise(
    total_deer = sum(Prevalence),
    pos_deer = sum(Prevalence[Status == "Pos"], na.rm = T),
    proportion_pos = pos_deer / total_deer,
    .groups = "drop"
  ) %>%
  mutate(Sex = "All")

deer_prop <- bind_rows(deer_prop_s, deer_prop_all,deer_prop_A1,deer_prop_A2)


deer_prop$Sex <- ifelse(deer_prop$Sex == 'F', 'Female', 
                        ifelse(deer_prop$Sex == 'M', 
                               'Male', 'All'))
deer_prop$Year <- as.character(deer_prop$year)
deer_prop$Year <- ifelse(deer_prop$Year == "1111", 'Mean', deer_prop$Year)
d <- ggplot(deer_prop, aes(x = age_nn, y = proportion_pos, pch = Year, group = year))+
  geom_line(data = deer_prop[deer_prop$Year == 'Mean',], lwd = .5, lty = 2)+
  geom_point(cex = 3)+
  ylab("CWD sample prevalence")+
  xlab('Age')+
  #geom_label(aes(label = total_deer,  y = proportion_pos- .05), col = 'black', size = 4, position = position_dodge(width = .5))+
  theme_linedraw()+
  facet_wrap(~Sex, nrow = 1)+
  ylim(0,1)+
  theme(legend.position = 'bottom', legend.position.inside = c(.8, .1), 
        axis.text = element_text(size = 22),
        axis.title = element_text(size = 20),
        legend.text = element_text(size = 20),
        legend.title = element_text(size  =20),
        strip.text = element_text(size = 20),
        strip.background = element_rect(fill = 'grey50'),
        aspect.ratio = 1)+
  guides(shape = guide_legend(nrow = 1))+
  scale_x_continuous(breaks = c(1.5, 2.5, 3.5, 4.5),
                     labels = c('1', '2', '3', '4+')
  )
d


#### NRR ####
### Lifetime Repro (for females):
repro_NF <- surv_neg_F
repro_PF <- surv_pos_F
for(k in 1:nrow(surv_betas)){
  for(m in 2:2920){
    if (m %% 365 != 0) {
      # if m/365 is NOT an integer
      repro_NF[k, m + 365] <- repro_NF[k, m + 364]
      repro_PF[k, m + 365] <- repro_PF[k, m + 364]
    } else {
      if(m > 365){ #adults
        repro_NF[k, m + 365] <- surv_neg_F[k, m + 365] * 1.7 + repro_NF[k, m + 364]
        repro_PF[k, m + 365] <- surv_pos_F[k, m + 365] * 1.7 + repro_PF[k, m + 364]
      } else{
        repro_NF[k, m + 365] <- surv_neg_F[k, m + 365] * .3 + repro_NF[k, m + 364] #bio data says ~ 30% lactation for yearlings; only one fawn usually
        repro_PF[k, m + 365] <- surv_pos_F[k, m + 365] * .3 + repro_PF[k, m + 364]
      }
    }
  }
}

repro_NF2 <- apply(repro_NF, 2, FUN = function(x){quantile(x, c(0.025, .5, .975))})
repro_PF2 <- apply(repro_PF, 2, FUN = function(x){quantile(x, c(0.025, .5, .975))})

repro_lines <- data.frame(Median = c(repro_NF2[2,], repro_PF2[2,]),
                          LCI = c(repro_NF2[1,], repro_PF2[1,]),
                          UCI = c(repro_NF2[3,], repro_PF2[3,]),
                          Stage = rep(c("Neg F", "Pos F"), each = ncol(repro_NF2)),
                          Status = c(rep(c('Neg', 'Pos'), each = ncol(repro_NF2))),
                          Age = rep(c(c(183:547)/365, 1.5+c(1:2920)/365), 2))

## figure not in manuscript:
ggplot(repro_lines, aes(x = Age, y = Median))+
  geom_line(aes(col = Status), lwd = 1)+
  geom_ribbon(aes(ymin = LCI, ymax = UCI, fill = Status), alpha = .5)+
  theme_bw()+
  xlim(1.5, 9.5)+
  ylab('Lifetime reproductive output')+
  scale_fill_manual(values = c('cadetblue2', 'firebrick3'))+
  scale_color_manual(values = c('navyblue', 'firebrick3'))+
  theme(axis.text = element_text(size = 15), axis.title = element_text(size = 20), legend.text = element_text(size = 12), legend.title = element_text(size = 15))

repro_lines[repro_lines$Age == 8.5, ]
