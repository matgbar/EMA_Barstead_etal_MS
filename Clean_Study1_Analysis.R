############################################################################################################
#Authors:     Barstead, M.G., DeYoung, K.D., Anderson, A. S., & Shackman, A. J.

#Title:       The moment-to-moment affective experience of dispositionally negative individuals 

#Contact:     barstead@umd.edu

#Contents:    Bayesian regression analyses - Study 1
############################################################################################################
############################################################################################################
#This script represents a final version of analyses agreed upon following a meeting with AJS and JH 
#Two major changes to our analytic approach were agreed upon
#Both changes were designed to better map our models to the underlying theoretical framework we were guided by

#Change 1 - dichotomize the the momentary ratings of negative events and best events 
# The decision was to make a variable that represents frequency (which is what our model suggests is linked to DN)

#Change 2 - split negative events/positive events into two separate models to reduce collinearity and look at overall effects 
############################################################################################################

library(brms)
library(rstan)
library(rstanarm)
library(bayesplot)
library(pan)
library(mitml)
library(mice)
library(parallel)
library(RColorBrewer)
library(ggridges)
library(riverplot)
library(ggalluvial)
library(ggpubr)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

###########################################################################################################
user<-ifelse(Sys.getenv("USERPROFILE")=="", "~", Sys.getenv("USERPROFILE"))
wd<-paste0(user, '/Dropbox/UMD/Shackman Lab/EMA_MS')
data.folder<-paste0(wd, '/Data')
study1.out<-paste0(wd, '/Study 1 output')
study1.graphics<-paste0(study1.out, '/Graphics')
study1.model<-paste0(study1.out, '/Model summaries')
stan.code<-paste0(wd, '/Stan_code')
EDA.folder<-paste0(study1.out, '/EDA')

#Loading Study 1 Data (from Emotion MS - Shackman et al. 2017)
load(paste0(data.folder, '/Emotion MS environment.RData'))

#If data scripts have already been run then
load(paste0(data.folder, '/Study_1_Data_Environment.RData'))
###########################################################################################################
#Dichotomizing Momentary Negative Affect 
#Extracting Variables to be used in imputation model: 
###########################################################################################################

#Getting vector of IDs (required for group mean centering values)
IDs<-unique(dat$subid)

ID<-vector()
Worst<-vector()
Best<-vector()
c.DN<-vector()
NegAff<-vector()
PosAff<-vector()

#Selecting additional dispositional variables to consider in the imputation. 
Gender<-vector()
BFI_E<-vector()
BFI_A<-vector()
BFI_C<-vector()
BFI_O<-vector()
GenDep<-vector()

#Bringing in covariates for imputation - 
#Want to include gender, BFI_E, BFI_C, BFI_O, BFI_A, GenDep

for(i in 1:length(IDs)){
  dat.temp<-dat[dat$subid==IDs[i],]
  for(j in 1:length(dat.temp[,1])){
    ID<-c(ID, IDs[i])
    Worst<-c(Worst, dat.temp$WorstEvent_Neg[j])
    Best<-c(Best, dat.temp$BestEvent_Pos[j])
    c.DN<-c(c.DN, dat.temp$ZAP_Both[j])
    NegAff<-c(NegAff, dat.temp$MAFS_NA[j])
    PosAff<-c(PosAff, dat.temp$MAFS_PA[j])
    Gender<-c(Gender, dat.temp$Gender[j])
    BFI_E<-c(BFI_E, dat.temp$Battery_BFI_E[j])
    BFI_A<-c(BFI_A, dat.temp$Battery_BFI_A[j])
    BFI_C<-c(BFI_C, dat.temp$Battery_BFI_C[j])
    BFI_O<-c(BFI_O, dat.temp$Battery_BFI_O[j])
    GenDep<-c(GenDep, dat.temp$Battery_IDASII64_GenDep[j])
  }
}

#Dichotomizing on the 75% percentile. 
#In practice this was really dichomotizing at the 83rd percentile and the 77th percentile 

Worst_dich<-ifelse(Worst>quantile(Worst, .75, na.rm=TRUE), 1, 0)
Best_dich<-ifelse(Best>quantile(Best, .75, na.rm = TRUE), 1, 0)

dat.study1<-data.frame(ID, 
                       Worst, 
                       Best,
                       Worst_dich, 
                       Best_dich, 
                       c.DN, 
                       NegAff,
                       PosAff,
                       Gender, 
                       BFI_A, 
                       BFI_C, 
                       BFI_E, 
                       BFI_O, 
                       GenDep)                
psych::describe(dat.study1)

dat.study1$Gender<-dat.study1$Gender-1  #Putting this in 0-1 coding
#Males = 1, Females = 0
dat.study1$c.DN<-as.numeric(scale(dat.study1$c.DN)) #Re-standardizing to give variable unit variance in models

###########################################################################################################
#Missingness Evaluation
###########################################################################################################
#inspecting missing patterns in the data: 
md.pattern(dat.study1)

#Exploring differences in missing vs. non-missing distributions of scores: 
#Missinginess is almost entirely clustered - will explore using single variable
dat.study1$miss<-ifelse(!is.na(dat.study1$PosAff), 'Complete', 'Missing')

g1<-ggplot(data=dat.study1, 
           aes(x=BFI_O, y=miss, fill=miss))+
  geom_density_ridges(alpha=.5)+
  ylab('')+
  theme_bw()+
  xlab('BFI Openess Scores (1-5)')+
  guides(fill = guide_legend(""))+
  scale_fill_brewer(palette = 'Set1')+
  geom_vline(xintercept = mean(dat.study1$BFI_O[dat.study1$miss=='Missing']), color=brewer.pal(3,'Set1')[2])+
  geom_vline(xintercept = mean(dat.study1$BFI_O[dat.study1$miss=='Complete']), color=brewer.pal(3,'Set1')[1])+
  annotate(geom = 'text', label=paste('Mean Complete =', 
                                      round(mean(dat.study1$BFI_O[dat.study1$miss=='Complete']), digits = 2)), 
           color=brewer.pal(3,'Set1')[1], 
           x=2.9, 
           y=.5)+
  annotate(geom = 'text', label=paste('Mean Missing =', 
                                      round(mean(dat.study1$BFI_O[dat.study1$miss=='Missing']), digits = 2)), 
           color=brewer.pal(3,'Set1')[2], 
           x=4.25, 
           y=.5)

g2<-ggplot(data=dat.study1, 
           aes(x=BFI_C, y=miss, fill=miss))+
  geom_density_ridges(alpha=.5)+
  ylab('')+
  theme_bw()+
  xlab('BFI Conscientiousness Scores (1-5)')+
  guides(fill = guide_legend(""))+
  scale_fill_brewer(palette = 'Set1')+
  geom_vline(xintercept = mean(dat.study1$BFI_C[dat.study1$miss=='Missing']), color=brewer.pal(3,'Set1')[2])+
  geom_vline(xintercept = mean(dat.study1$BFI_C[dat.study1$miss=='Complete']), color=brewer.pal(3,'Set1')[1])+
  annotate(geom = 'text', label=paste('Mean Complete =', 
                                      round(mean(dat.study1$BFI_C[dat.study1$miss=='Complete']), digits = 2)), 
           color=brewer.pal(3,'Set1')[1], 
           x=4.35, 
           y=.5)+
  annotate(geom = 'text', label=paste('Mean Missing =', 
                                      round(mean(dat.study1$BFI_C[dat.study1$miss=='Missing']), digits = 2)), 
           color=brewer.pal(3,'Set1')[2], 
           x=2.75, 
           y=.5)

g3<-ggplot(data=dat.study1, 
           aes(x=BFI_E, y=miss, fill=miss))+
  geom_density_ridges(alpha=.5)+
  ylab('')+
  theme_bw()+
  xlab('BFI Extraversion Scores (1-5)')+
  guides(fill = guide_legend(""))+
  scale_fill_brewer(palette = 'Set1')+
  geom_vline(xintercept = mean(dat.study1$BFI_E[dat.study1$miss=='Missing']), color=brewer.pal(3,'Set1')[2])+
  geom_vline(xintercept = mean(dat.study1$BFI_E[dat.study1$miss=='Complete']), color=brewer.pal(3,'Set1')[1])+
  annotate(geom = 'text', label=paste('Mean Complete =', 
                                      round(mean(dat.study1$BFI_E[dat.study1$miss=='Complete']), digits = 2)), 
           color=brewer.pal(3,'Set1')[1], 
           x=2.2, 
           y=.5)+
  annotate(geom = 'text', label=paste('Mean Missing =', 
                                      round(mean(dat.study1$BFI_E[dat.study1$miss=='Missing']), digits = 2)), 
           color=brewer.pal(3,'Set1')[2], 
           x=4, 
           y=.5)


g4<-ggplot(data=dat.study1, 
           aes(x=BFI_A, y=miss, fill=miss))+
  geom_density_ridges(alpha=.5)+
  ylab('')+
  theme_bw()+
  xlab('BFI Agreeableness Scores (1-5)')+
  guides(fill = guide_legend(""))+
  scale_fill_brewer(palette = 'Set1')+
  geom_vline(xintercept = mean(dat.study1$BFI_A[dat.study1$miss=='Missing']), color=brewer.pal(3,'Set1')[2])+
  geom_vline(xintercept = mean(dat.study1$BFI_A[dat.study1$miss=='Complete']), color=brewer.pal(3,'Set1')[1])+
  annotate(geom = 'text', label=paste('Mean Complete =', 
                                      round(mean(dat.study1$BFI_A[dat.study1$miss=='Complete']), digits = 2)), 
           color=brewer.pal(3,'Set1')[1], 
           x=4.5, 
           y=.5)+
  annotate(geom = 'text', label=paste('Mean Missing =', 
                                      round(mean(dat.study1$BFI_A[dat.study1$miss=='Missing']), digits = 2)), 
           color=brewer.pal(3,'Set1')[2], 
           x=3.15, 
           y=.5)

g5<-ggplot(data=dat.study1, 
           aes(x=c.DN, y=miss, fill=miss))+
  geom_density_ridges(alpha=.5)+
  ylab('')+
  theme_bw()+
  xlab('Centered DN scores')+
  guides(fill = guide_legend(""))+
  scale_fill_brewer(palette = 'Set1')+
  geom_vline(xintercept = mean(dat.study1$c.DN[dat.study1$miss=='Missing']), color=brewer.pal(3,'Set1')[2])+
  geom_vline(xintercept = mean(dat.study1$c.DN[dat.study1$miss=='Complete']), color=brewer.pal(3,'Set1')[1])+
  annotate(geom = 'text', label=paste('Mean Complete =', 
                                      round(mean(dat.study1$c.DN[dat.study1$miss=='Complete']), digits = 2)), 
           color=brewer.pal(3,'Set1')[1], 
           x=-.95, 
           y=.5)+
  annotate(geom = 'text', label=paste('Mean Missing =', 
                                      round(mean(dat.study1$c.DN[dat.study1$miss=='Missing']), digits = 2)), 
           color=brewer.pal(3,'Set1')[2], 
           x=.95, 
           y=.5)

g6<-ggplot(data=dat.study1, 
           aes(x=GenDep, y=miss, fill=miss))+
  geom_density_ridges(alpha=.5)+
  ylab('')+
  theme_bw()+
  xlab('IDAS Depression Scores')+
  guides(fill = guide_legend(""))+
  scale_fill_brewer(palette = 'Set1')+
  geom_vline(xintercept = mean(dat.study1$GenDep[dat.study1$miss=='Missing']), color=brewer.pal(3,'Set1')[2])+
  geom_vline(xintercept = mean(dat.study1$GenDep[dat.study1$miss=='Complete']), color=brewer.pal(3,'Set1')[1])+
  annotate(geom = 'text', label=paste('Mean Complete =', 
                                      round(mean(dat.study1$GenDep[dat.study1$miss=='Complete']), digits = 2)), 
           color=brewer.pal(3,'Set1')[1], 
           x=30, 
           y=.5)+
  annotate(geom = 'text', label=paste('Mean Missing =', 
                                      round(mean(dat.study1$GenDep[dat.study1$miss=='Missing']), digits = 2)), 
           color=brewer.pal(3,'Set1')[2], 
           x=55, 
           y=.5)
png(paste0(EDA.folder, '/Missing_Data_Exploration.png'), 
    res=600, 
    units='in', 
    height=10, 
    width =16)
cowplot::plot_grid(g1,g2,g3,g4,g5,g6, nrow = 3)
dev.off()

#Nothing obvious in this relative unsophisticated way of looking at the data
#Will Model missingness using a random effects approach. 
dat.study1$miss.dich<-ifelse(dat.study1$miss=='Missing', 1, 0)
dat.study1$zGenDep<-as.numeric(scale(dat.study1$GenDep))

fit.miss<-lme4::glmer(miss.dich~1+c.DN+Gender+zGenDep+BFI_A+BFI_C+BFI_E+BFI_O+
                        (1|ID), 
                      data=dat.study1, 
                      family='binomial',
                      control=glmerControl(optimizer="bobyqa",
                                           optCtrl=list(maxfun=2e5)))
summary(fit.miss)
#Suggests that the following variables are related to missingness: 
#   1. BFI Agreeableness (less likely to have missing data)
#   2. BFI Conscientiouness (less likely to have missing data)
#   3. BFI Openness (more likely to have missing data... interesting)

#=========================================================================================================
#Imputation Model
#=========================================================================================================
fml<- Worst + Best + NegAff + PosAff  ~ 
  1 + c.DN + BFI_A + BFI_C + BFI_O + (1|ID)

#Set number of data sets to impute 
M<-20
imp<-panImpute(dat.study1, 
               formula=fml, 
               n.burn=10000, 
               n.iter = 5000, 
               m=M, 
               seed = 0716)

#plotting results - assessing convergence on final distributions
dat.imp<-mitmlComplete(imp)

dat.long<-data.frame()
for(i in 1:M){
  dat.temp<-dat.imp[[i]]
  dat.temp$IMP<-rep(i, length(dat.temp[,1]))
  dat.long<-rbind(dat.long, dat.temp)
}

dat.study1$IMP<-rep(0, length(dat.study1[,1]))
dat.long<-rbind(dat.study1, dat.long)
dat.long$Orig<-ifelse(dat.long$IMP<1, 'Original', 'Imputed')

#Plotting Negative Affect Imputation
g1<-ggplot()+
  geom_freqpoly(data = dat.long, aes(x=NegAff, group=IMP, color=Orig), stat = 'density')+
  theme(legend.title = element_blank())+
  ggtitle('Imputed Distributions of Momentary Negative Affect')+
  ylab('Density')+
  xlab('Momentary Negative Affect')
g1

png(paste0(EDA.folder, '/S1_Imputation_NegativeAffect.png'), res=300, 
    units='in', height=5, width = 8)
g1
dev.off()

#Plotting Positive Affect Imputation
g1<-ggplot()+
  geom_freqpoly(data = dat.long, aes(x=PosAff, group=IMP, color=Orig), stat = 'density')+
  theme(legend.title = element_blank())+
  ggtitle('Imputed Distributions of Momentary Positive Affect')+
  ylab('Density')+
  xlab('Momentary Positive Affect')
g1

png(paste0(EDA.folder, '/S1_Imputation_PositiveAffect.png'), res=300, 
    units='in', height=5, width = 8)
g1
dev.off()

#Plotting Deviations from Mean Worst Ratings
g1<-ggplot()+
  geom_freqpoly(data = dat.long, aes(x=Worst, group=IMP, color=Orig), stat = 'density')+
  theme(legend.title = element_blank())+
  ggtitle('Imputed Distributions of Worst Event')+
  ylab('Density')+
  xlab('Worst Event Ratings')
g1

png(paste0(EDA.folder, '/S1_Imputation_WorstEvent.png'), res=300, 
    units='in', height=5, width = 8)
g1
dev.off()

#Plotting Deviations from Mean Best Ratings
g1<-ggplot()+
  geom_freqpoly(data = dat.long, aes(x=Best, group=IMP, color=Orig), stat = 'density')+
  theme(legend.title = element_blank())+
  ggtitle('Imputed Distributions of Best Event')+
  ylab('Density')+
  xlab('Best Event Ratings')
g1

png(paste0(EDA.folder, '/S1_Imputation_BestEvent.png'), res=300, 
    units='in', height=5, width = 8)
g1
dev.off()

#summarizing imputations - checking for convergence
sink(paste0(study1.out, '/Imputation_Summary.txt'))
summary(imp)
psych::describeBy(dat.long, group='IMP')
sink()    #warnings have to do with factor variables for the most part

###########################################################################################################
#Converting Imputed Data for Use
#Requires following steps
#1. Restoring variable properties (i.e., truncating DVs at 1 and 5)
#1. Dichomotizing using 75% percentile as cut point 
#2. Identifying individual mean proportions for each event type 
#3. Group-mean centering each of the dichotomized variables - ensures full decomposition 
###########################################################################################################

for(i in 1:M){
  dat.imp[[i]]$T.NegAff<-ifelse(dat.imp[[i]]$NegAff>5, 5,dat.imp[[i]]$NegAff)
  dat.imp[[i]]$T.NegAff<-ifelse(dat.imp[[i]]$T.NegAff<1, 1, dat.imp[[i]]$T.NegAff)
  dat.imp[[i]]$T.PosAff<-ifelse(dat.imp[[i]]$PosAff>5, 5, dat.imp[[i]]$PosAff)
  dat.imp[[i]]$T.PosAff<-ifelse(dat.imp[[i]]$T.PosAff<1, 1, dat.imp[[i]]$T.PosAff)
  dat.imp[[i]]$T.Worst<-ifelse(dat.imp[[i]]$Worst>5, 5, dat.imp[[i]]$Worst)
  dat.imp[[i]]$T.Worst<-ifelse(dat.imp[[i]]$T.Worst<1, 1, dat.imp[[i]]$T.Worst)
  dat.imp[[i]]$T.Best<-ifelse(dat.imp[[i]]$Best>5, 5, dat.imp[[i]]$Best)
  dat.imp[[i]]$T.Best<-ifelse(dat.imp[[i]]$T.Best<1, 1, dat.imp[[i]]$T.Best)
  dat.imp[[i]]$Worst_dich<-ifelse(dat.imp[[i]]$T.Worst>3, 1, 0)
  dat.imp[[i]]$Best_dich<-ifelse(dat.imp[[i]]$T.Best>3, 1, 0)
}

#Quirky formatting problems - getting around it this way... 
#   The issue is that brm_multiple requires a list of datasets
#   It cannot recognize an mitml.list object though
#   So I made a list out of list
dat.study1.list<-list()

for(i in 1:M){
  dat.study1.list[[i]]<-dat.imp[[i]]
}

#Obtaining individual mean ratings for Best and Worst events
#Individually mean centering predictors to decompose between and within effects

for(i in 1:M){
  IDs<-unique(dat.study1$ID)
  DF.temp<-dat.study1.list[[i]]
  for(j in 1:length(IDs)){
    DF.temp$c.Worst[DF.temp$ID==IDs[j]]<-
      DF.temp$Worst_dich[DF.temp$ID==IDs[j]] - 
      mean(DF.temp$Worst_dich[DF.temp$ID==IDs[j]])
    
    DF.temp$prop_Worst[DF.temp$ID==IDs[j]]<-rep(mean(DF.temp$Worst_dich[DF.temp$ID==IDs[j]]))
    
    DF.temp$c.Best[DF.temp$ID==IDs[j]]<-
      DF.temp$Best_dich[DF.temp$ID==IDs[j]] - 
      mean(DF.temp$Best_dich[DF.temp$ID==IDs[j]])
    
    DF.temp$prop_Best[DF.temp$ID==IDs[j]]<-rep(mean(DF.temp$Best_dich[DF.temp$ID==IDs[j]]))
    
    #For Eventual Inclusion in variance decomp 
    DF.temp$DNxWorst<-DF.temp$c.DN*DF.temp$c.Worst
    DF.temp$DNxBest <-DF.temp$c.DN*DF.temp$c.Best
  
  }
  dat.study1.list[[i]]<-DF.temp
}

##########################################################################################################
#*** IMPORTANT!!!!!!! ***
#Change stan code to map onto to each model model (i.e.,. change file name) 
#*** IMPORTANT!!!!!!! ***
#*** IMPORTANT!!!!!!! ***
#*** IMPORTANT!!!!!!! ***
#*** IMPORTANT!!!!!!! ***
##########################################################################################################


###########################################################################################################
#Negative Affect Models
###########################################################################################################
mean(log(dat.study1$NegAff), na.rm=TRUE)
sd(log(dat.study1$NegAff), na.rm = TRUE)

Neg_ucm<-brms::brm_multiple(T.NegAff~1+(1|ID), 
                            data = dat.study1.list, 
                            family = 'lognormal',
                            prior = c(set_prior('normal(0,5)', class='sd'), 
                                      set_prior('normal(0,1)', class='Intercept')),
                            warmup = 2000, 
                            iter = 3000,
                            chains = 3,
                            control = list(adapt_delta=.99, 
                                           max_treedepth=15), 
                            save_model = paste0(stan.code, '/S1_ucm.stan'))
pp_check(Neg_ucm, nsamples = 50)
round(Neg_ucm$rhats[,1:3], 3)

#----------------------------------------------------------------------------------------------------------
#Worst Event Models
#----------------------------------------------------------------------------------------------------------

#--
Neg_Worst<-brms::brm_multiple(T.NegAff~1+c.Worst+(1+c.Worst|ID), 
                            data = dat.study1.list, 
                            family = 'lognormal',
                            prior = c(set_prior('normal(0,5)', class='sd'), 
                                      set_prior('normal(0,1)', class='Intercept'),
                                      set_prior('lkj(2)', class = 'cor'), 
                                      set_prior('normal(0,5)', class='b')),
                            warmup = 3000, 
                            iter = 4000,
                            chains = 2,
                            control = list(adapt_delta=.95, 
                                           max_treedepth=15), 
                            save_model = paste0(stan.code, '/S1_ucm.stan'))

#--
Neg_Worst_c.DN<-brms::brm_multiple(T.NegAff~1+c.Worst+
                                 c.DN+
                                 (1+c.Worst|ID), 
                                data = dat.study1.list, 
                                family = 'lognormal',
                                prior = c(set_prior('normal(0,5)', class='sd'), 
                                          set_prior('normal(0,1)', class='Intercept'),
                                          set_prior('lkj(2)', class = 'cor'), 
                                          set_prior('normal(0,5)', class='b')),
                                warmup = 3000, 
                                iter = 4000,
                                chains = 2,
                                control = list(adapt_delta=.95, 
                                               max_treedepth=15), 
                                save_model = paste0(stan.code, '/S1_ucm.stan'))

#--
Neg_Worst_prop_Worst<-brms::brm_multiple(T.NegAff~1+c.Worst+
                                 prop_Worst+
                                 (1+c.Worst|ID), 
                               data = dat.study1.list, 
                               family = 'lognormal',
                               prior = c(set_prior('normal(0,5)', class='sd'), 
                                         set_prior('normal(0,1)', class='Intercept'),
                                         set_prior('lkj(2)', class = 'cor'), 
                                         set_prior('normal(0,5)', class='b')),
                               warmup = 3000, 
                               iter = 4000,
                               chains = 2,
                               control = list(adapt_delta=.95, 
                                              max_treedepth=15), 
                               save_model = paste0(stan.code, '/S1_ucm.stan'))

#--
Neg_Worst_All<-brms::brm_multiple(T.NegAff~1+c.Worst+
                                         c.DN+prop_Worst+
                                         (1+c.Worst|ID), 
                                       data = dat.study1.list, 
                                       family = 'lognormal',
                                       prior = c(set_prior('normal(0,5)', class='sd'), 
                                                 set_prior('normal(0,1)', class='Intercept'),
                                                 set_prior('lkj(2)', class = 'cor'), 
                                                 set_prior('normal(0,5)', class='b')),
                                       warmup = 3000, 
                                       iter = 4000,
                                       chains = 2,
                                       control = list(adapt_delta=.95, 
                                                      max_treedepth=15), 
                                       save_model = paste0(stan.code, '/S1_ucm.stan'))

#--
Neg_Worst_cross<-brms::brm_multiple(T.NegAff~1+c.Worst+
                                c.DN+prop_Worst+
                                c.Worst:c.DN+
                                  (1+c.Worst|ID), 
                                data = dat.study1.list, 
                                family = 'lognormal',
                                prior = c(set_prior('normal(0,5)', class='sd'), 
                                          set_prior('normal(0,1)', class='Intercept'),
                                          set_prior('lkj(2)', class = 'cor'), 
                                          set_prior('normal(0,5)', class='b')),
                                warmup = 3000, 
                                iter = 4000,
                                chains = 2,
                                control = list(adapt_delta=.95, 
                                               max_treedepth=15), 
                                save_model = paste0(stan.code, '/S1_ucm.stan'))

#----------------------------------------------------------------------------------------------------------
#Best Event Models
#----------------------------------------------------------------------------------------------------------

#--
Neg_Best<-brms::brm_multiple(T.NegAff~1+c.Best+(1+c.Best|ID), 
                              data = dat.study1.list, 
                              family = 'lognormal',
                              prior = c(set_prior('normal(0,5)', class='sd'), 
                                        set_prior('normal(0,1)', class='Intercept'),
                                        set_prior('lkj(2)', class = 'cor'), 
                                        set_prior('normal(0,5)', class='b')),
                              warmup = 3000, 
                              iter = 4000,
                              chains = 2,
                              control = list(adapt_delta=.95, 
                                             max_treedepth=15), 
                              save_model = paste0(stan.code, '/S1_ucm.stan'))

#--
Neg_Best_c.DN<-brms::brm_multiple(T.NegAff~1+c.Best+
                                     c.DN+
                                     (1+c.Best|ID), 
                                   data = dat.study1.list, 
                                   family = 'lognormal',
                                   prior = c(set_prior('normal(0,5)', class='sd'), 
                                             set_prior('normal(0,1)', class='Intercept'),
                                             set_prior('lkj(2)', class = 'cor'), 
                                             set_prior('normal(0,5)', class='b')),
                                   warmup = 3000, 
                                   iter = 4000,
                                   chains = 2,
                                   control = list(adapt_delta=.95, 
                                                  max_treedepth=15), 
                                   save_model = paste0(stan.code, '/S1_ucm.stan'))

#--
Neg_Best_prop_Best<-brms::brm_multiple(T.NegAff~1+c.Best+
                                           prop_Best+
                                           (1+c.Best|ID), 
                                         data = dat.study1.list, 
                                         family = 'lognormal',
                                         prior = c(set_prior('normal(0,5)', class='sd'), 
                                                   set_prior('normal(0,1)', class='Intercept'),
                                                   set_prior('lkj(2)', class = 'cor'), 
                                                   set_prior('normal(0,5)', class='b')),
                                         warmup = 3000, 
                                         iter = 4000,
                                         chains = 2,
                                         control = list(adapt_delta=.95, 
                                                        max_treedepth=15), 
                                         save_model = paste0(stan.code, '/S1_ucm.stan'))

#--
Neg_Best_All<-brms::brm_multiple(T.NegAff~1+c.Best+
                                    c.DN+prop_Best+
                                    (1+c.Best|ID), 
                                  data = dat.study1.list, 
                                  family = 'lognormal',
                                  prior = c(set_prior('normal(0,5)', class='sd'), 
                                            set_prior('normal(0,1)', class='Intercept'),
                                            set_prior('lkj(2)', class = 'cor'), 
                                            set_prior('normal(0,5)', class='b')),
                                  warmup = 3000, 
                                  iter = 4000,
                                  chains = 2,
                                  control = list(adapt_delta=.95, 
                                                 max_treedepth=15), 
                                  save_model = paste0(stan.code, '/S1_ucm.stan'))

#--
Neg_Best_cross<-brms::brm_multiple(T.NegAff~1+c.Best+
                                      c.DN+prop_Best+
                                      c.Best:c.DN+
                                      (1+c.Best|ID), 
                                    data = dat.study1.list, 
                                    family = 'lognormal',
                                    prior = c(set_prior('normal(0,5)', class='sd'), 
                                              set_prior('normal(0,1)', class='Intercept'),
                                              set_prior('lkj(2)', class = 'cor'), 
                                              set_prior('normal(0,5)', class='b')),
                                    warmup = 3000, 
                                    iter = 4000,
                                    chains = 2,
                                    control = list(adapt_delta=.95, 
                                                   max_treedepth=15), 
                                    save_model = paste0(stan.code, '/S1_ucm.stan'))

###########################################################################################################
#Plotting Results - Riverplot creation
###########################################################################################################
#Rights and Sterba (2018) function: 
r2MLM <- function(data,within_covs,between_covs,random_covs,
                  gamma_w,gamma_b,Tau,sigma2,has_intercept=T,clustermeancentered=T){
  #browser()
  if(has_intercept==T){
    if(length(gamma_b)>1) gamma <- c(1,gamma_w,gamma_b[2:length(gamma_b)])
    if(length(gamma_b)==1) gamma <- c(1,gamma_w)
    if(is.null(within_covs)==T) gamma_w <- 0
  }
  if(has_intercept==F){
    gamma <- c(gamma_w,gamma_b)
    if(is.null(within_covs)==T) gamma_w <- 0
    if(is.null(between_covs)==T) gamma_b <- 0
  }
  if(is.null(gamma)) gamma <- 0
  ##compute phi
  phi <- var(cbind(1,data[,c(within_covs)],data[,c(between_covs)]),na.rm=T)
  if(has_intercept==F) phi <- var(cbind(data[,c(within_covs)],data[,c(between_covs)]),na.rm=T)
  if(is.null(within_covs)==T & is.null(within_covs)==T & has_intercept==F) phi <- 0
  phi_w <- var(data[,within_covs],na.rm=T)
  if(is.null(within_covs)==T) phi_w <- 0
  phi_b <- var(cbind(1,data[,between_covs]),na.rm=T)
  if(is.null(between_covs)==T) phi_b <- 0
  ##compute psi and kappa
  var_randomcovs <- var(cbind(1,data[,c(random_covs)]),na.rm=T)
  if(length(Tau)>1) psi <- matrix(c(diag(Tau)),ncol=1)
  if(length(Tau)==1) psi <- Tau
  if(length(Tau)>1) kappa <- matrix(c(Tau[lower.tri(Tau)==TRUE]),ncol=1)
  if(length(Tau)==1) kappa <- 0
  v <- matrix(c(diag(var_randomcovs)),ncol=1)
  r <- matrix(c(var_randomcovs[lower.tri(var_randomcovs)==TRUE]),ncol=1)
  if(is.null(random_covs)==TRUE){
    v <- 0
    r <- 0
    m <- matrix(1,ncol=1)
  }
  if(length(random_covs)>0) m <- matrix(c(colMeans(cbind(1,data[,c(random_covs)]),na.rm=T)),ncol=1)
  ##total variance
  totalvar_notdecomp <- t(v)%*%psi + 2*(t(r)%*%kappa) + t(gamma)%*%phi%*%gamma + t(m)%*%Tau%*%m + sigma2
  totalwithinvar <- (t(gamma_w)%*%phi_w%*%gamma_w) + (t(v)%*%psi + 2*(t(r)%*%kappa)) + sigma2
  totalbetweenvar <- (t(gamma_b)%*%phi_b%*%gamma_b) + Tau[1]
  totalvar <- totalwithinvar + totalbetweenvar
  ##total decomp
  decomp_fixed_notdecomp <- (t(gamma)%*%phi%*%gamma) / totalvar
  decomp_fixed_within <- (t(gamma_w)%*%phi_w%*%gamma_w) / totalvar
  decomp_fixed_between <- (t(gamma_b)%*%phi_b%*%gamma_b) / totalvar
  decomp_fixed <- decomp_fixed_within + decomp_fixed_between
  decomp_varslopes <- (t(v)%*%psi + 2*(t(r)%*%kappa)) / totalvar
  decomp_varmeans <- (t(m)%*%Tau%*%m) / totalvar
  decomp_sigma <- sigma2/totalvar
  ##within decomp
  decomp_fixed_within_w <- (t(gamma_w)%*%phi_w%*%gamma_w) / totalwithinvar
  decomp_varslopes_w <- (t(v)%*%psi + 2*(t(r)%*%kappa)) / totalwithinvar
  decomp_sigma_w <- sigma2/totalwithinvar
  ##between decomp
  decomp_fixed_between_b <- (t(gamma_b)%*%phi_b%*%gamma_b) / totalbetweenvar
  decomp_varmeans_b <- Tau[1] / totalbetweenvar
  #NEW measures
  if (clustermeancentered==TRUE){
    R2_f <- decomp_fixed
    R2_f1 <- decomp_fixed_within
    R2_f2 <- decomp_fixed_between
    R2_fv <- decomp_fixed + decomp_varslopes
    R2_fvm <- decomp_fixed + decomp_varslopes + decomp_varmeans
    R2_v <- decomp_varslopes
    R2_m <- decomp_varmeans
    R2_f_w <- decomp_fixed_within_w
    R2_f_b <- decomp_fixed_between_b
    R2_fv_w <- decomp_fixed_within_w + decomp_varslopes_w
    R2_v_w <- decomp_varslopes_w
    R2_m_b <- decomp_varmeans_b
  }
  if (clustermeancentered==FALSE){
    R2_f <- decomp_fixed_notdecomp
    R2_fv <- decomp_fixed_notdecomp + decomp_varslopes
    R2_fvm <- decomp_fixed_notdecomp + decomp_varslopes + decomp_varmeans
    R2_v <- decomp_varslopes
    R2_m <- decomp_varmeans
  }
  if(clustermeancentered==TRUE){
    decomp_table <- matrix(c(decomp_fixed_within,decomp_fixed_between,decomp_varslopes,decomp_varmeans,decomp_sigma,
                             decomp_fixed_within_w,"NA",decomp_varslopes_w,"NA",decomp_sigma_w,
                             "NA",decomp_fixed_between_b,"NA",decomp_varmeans_b,"NA"),ncol=3)
    rownames(decomp_table) <- c("fixed, within","fixed, between","slope variation","mean variation","sigma2")
    colnames(decomp_table) <- c("total","within","between")
    R2_table <- matrix(c(R2_f1,R2_f2,R2_v,R2_m,R2_f,R2_fv,R2_fvm,
                         R2_f_w,"NA",R2_v_w,"NA","NA",R2_fv_w,"NA",
                         "NA",R2_f_b,"NA",R2_m_b,"NA","NA","NA")
                       ,ncol=3)
    rownames(R2_table) <- c("f1","f2","v","m","f","fv","fvm")
    colnames(R2_table) <- c("total","within","between")
  }
  ##barchart
  if(clustermeancentered==TRUE){
    contributions_stacked <- matrix(c(decomp_fixed_within,decomp_fixed_between,decomp_varslopes,decomp_varmeans,decomp_sigma,
                                      decomp_fixed_within_w,0,decomp_varslopes_w,0,decomp_sigma_w,
                                      0,decomp_fixed_between_b,0,decomp_varmeans_b,0),5,3)
    colnames(contributions_stacked) <- c("total","within","between")
    rownames(contributions_stacked) <- c("fixed slopes (within)",
                                         "fixed slopes (between)",
                                         "slope variation (within)",
                                         "intercept variation (between)",
                                         "residual (within)")
    #barplot(contributions_stacked, main="Decomposition", horiz=FALSE,
    #        ylim=c(0,1),col=c("darkred","steelblue","darkred","midnightblue","white"),ylab="proportion of variance",
    #        density=c(NA,NA,30,40,NA),angle=c(0,45,0,135,0),xlim=c(0,1),width=c(.3,.3))
    #legend(.30,-.1,legend=rownames(contributions_stacked),fill=c("darkred","steelblue","darkred","midnightblue","white"),
    #       cex=.7, pt.cex = 1,xpd=T,density=c(NA,NA,30,40,NA),angle=c(0,45,0,135,0))
  }
  if(clustermeancentered==FALSE){
    decomp_table <- matrix(c(decomp_fixed_notdecomp,decomp_varslopes,decomp_varmeans,decomp_sigma),ncol=1)
    rownames(decomp_table) <- c("fixed","slope variation","mean variation","sigma2")
    colnames(decomp_table) <- c("total")
    R2_table <- matrix(c(R2_f,R2_v,R2_m,R2_fv,R2_fvm),ncol=1)
    rownames(R2_table) <- c("f","v","m","fv","fvm")
    colnames(R2_table) <- c("total")
    ##barchar
    contributions_stacked <- matrix(c(decomp_fixed_notdecomp,decomp_varslopes,decomp_varmeans,decomp_sigma),4,1)
    colnames(contributions_stacked) <- c("total")
    rownames(contributions_stacked) <- c("fixed slopes",
                                         "slope variation",
                                         "intercept variation",
                                         "residual")
    #barplot(contributions_stacked, main="Decomposition", horiz=FALSE,
    #        ylim=c(0,1),col=c("darkblue","darkblue","darkblue","white"),ylab="proportion of variance",
    #        density=c(NA,30,40,NA),angle=c(0,0,135,0),xlim=c(0,1),width=c(.6))
    #legend(.30,-.1,legend=rownames(contributions_stacked),fill=c("darkblue","darkblue","darkblue","white"),
    #       cex=.7, pt.cex = 1,xpd=TRUE,density=c(NA,30,40,NA),angle=c(0,0,135,0))
  }
  Output <- list(noquote(decomp_table),noquote(R2_table))
  names(Output) <- c("Decompositions","R2s")
  return(Output)
}

#==========================================================================================================
#NEGATIVE MOOD MODELS - Variance Decomposition
# :: Negative Mood - Worst Event Models ::
#==========================================================================================================

#==========================================================================================================
#----------------------------------------------------------------------------------------------------------
#Negative Mood, Negative Events, and Dispositional Negativity: 
#==========================================================================================================
#Will need the grand intercept from the null model (on the link scale)
beta_00<-fixef(Neg_ucm)[1,1]
post_iter<-2000   #No. of chains per data set x No. post-warmup iterations per chain
                  #for loops below require this value is correctly specified

#Simple Variance Estimate From Level 1 Equation Only 
#=-=-=-=-=-=-=-=-=-=-=-=-=-=
within_cov<-c(22)                    #Columns with group-mean centered predictors
random_cov<-c(22)                    #Level 1 predictors with random slopes

#Getting distributions of fixed effects
Intercept = posterior_samples(Neg_Worst, pars = 'b_Intercept') 
Worst = posterior_samples(Neg_Worst, pars = 'b_c.Worst') 

#Getting between-subjects variance terms for random effects
Int_var = posterior_samples(Neg_Worst, pars = 'sd_ID__Intercept')^2
Worst_var = posterior_samples(Neg_Worst, pars = 'sd_ID__c.Worst')^2 

#Getting covariances of random effects (i.e., tau matrix)
cov_Int_Worst = posterior_samples(Neg_Worst, pars = 'sd_ID__Intercept')*
  posterior_samples(Neg_Worst, pars = 'sd_ID__c.Worst')*
  posterior_samples(Neg_Worst, pars = 'cor_ID__Intercept__c.Worst')

#Getting level 1 error variance
sigma<-log(1+posterior_samples(Neg_Worst, pars = 'sigma')/beta_00)

post_samples<-data.frame(Intercept,
                         Worst, 
                         Int_var, 
                         Worst_var, 
                         cov_Int_Worst, 
                         sigma)

colnames(post_samples)<-c('Intercept', 
                          'Worst', 
                          'Int_var', 
                          'Worst_var', 
                          'cov_Int_Worst', 
                          'sigma')

#Aggregate across imputed datasets
#Currently going to take 1000 draws from posterior distributions
#Will then apply across all 20 data sets

Neg_Worst_var<-list(between_var=vector(), 
                         within_var=vector(), 
                         between_All_tot=vector(),
                         between_All_btw=vector(),
                         between_res_btw=vector(),
                         within_fix_wthn=vector(),
                         within_fix_tot=vector(),
                         within_slope_var_wthn=vector(),
                         within_res_wthn=vector(), 
                         within_unmod_tot=vector())

for(i in 1:length(dat.study1.list)){
  print(paste('Starting w/ imputed set', i, 'out of', M))
  samp_frame<-1:post_iter+(i-1)*post_iter
  D<-sample(samp_frame, size = 400, replace = FALSE)
  for(d in 1:length(D)){
    #browser()
    Gamma_w<-as.vector(
      as.matrix(
        post_samples[D[d], c('Worst')]))       #Make sure the effects line up - in order of within_cov
    Gamma_b<-as.vector(
      as.matrix(
        post_samples[D[d], c('Intercept')]))    #level fixed intercept and any level 2 fixed effects
    tau<-c(post_samples[D[d], c('Int_var', 'cov_Int_Worst')], 
           post_samples[D[d], c('cov_Int_Worst', 'Worst_var')])
    
    tau<-matrix(unlist(tau), 
                byrow=TRUE, 
                ncol=2)
    
    #Matrix rows/columns ordered starting with tau.00 (intercept variance)
    #then add columns/rows for random effects of slopes in order of within_cov
    
    Sigma2<-post_samples[D[d], 'sigma']
    
    r2mlm_DN<-r2MLM(data=dat.study1.list[[i]], 
                    within_covs = within_cov,
                    between_covs = NULL,
                    random_covs = random_cov, 
                    gamma_w = Gamma_w, 
                    gamma_b = Gamma_b, 
                    Tau = tau, 
                    sigma2 = Sigma2, 
                    has_intercept = TRUE, 
                    clustermeancentered = TRUE)
    
    #Extracting relevant values across imputed datasets
    Neg_Worst_var$between_var<-c(Neg_Worst_var$between_var, as.numeric(r2mlm_DN$Decompositions['fixed, between', 'total'])+
                                        as.numeric(r2mlm_DN$Decompositions['mean variation', 'total']))
    Neg_Worst_var$within_var<-c(Neg_Worst_var$within_var, as.numeric(r2mlm_DN$Decompositions['sigma2', 'total'])+
                                       as.numeric(r2mlm_DN$Decompositions['slope variation', 'total'])+
                                       as.numeric(r2mlm_DN$Decompositions['fixed, within', 'total']))
    Neg_Worst_var$between_All_tot<-c(Neg_Worst_var$between_All_tot, as.numeric(r2mlm_DN$R2s['f2', 'total']))
    Neg_Worst_var$between_All_btw<-c(Neg_Worst_var$between_All_btw, as.numeric(r2mlm_DN$R2s['f2', 'between']))
    Neg_Worst_var$between_res_btw<-c(Neg_Worst_var$between_res_btw, as.numeric(r2mlm_DN$R2s['m', 'between']))
    Neg_Worst_var$within_fix_wthn<-c(Neg_Worst_var$within_fix_wthn, as.numeric(r2mlm_DN$Decompositions['fixed, within', 'within']))
    Neg_Worst_var$within_fix_tot<-c(Neg_Worst_var$within_fix_tot, as.numeric(r2mlm_DN$Decompositions['fixed, within', 'total']))
    Neg_Worst_var$within_slope_var_wthn<-c(Neg_Worst_var$within_slope_var_wthn, as.numeric(r2mlm_DN$Decompositions['slope variation', 'within']))
    Neg_Worst_var$within_res_wthn<-c(Neg_Worst_var$within_res_wthn, as.numeric(r2mlm_DN$Decompositions['sigma2', 'within']))
    Neg_Worst_var$within_unmod_tot<-c(Neg_Worst_var$within_unmod_tot, as.numeric(r2mlm_DN$Decompositions['sigma2', 'total'])+
                                             as.numeric(r2mlm_DN$Decompositions['slope variation', 'total']))
    
    print(paste("Sampling", d, "out of", length(D), '- Dataset:', i))
  }
  print(paste("*** Finshed imputed data set", i, 'out of', M, "***"))
}
hist(Neg_Worst_var$between_var)
hist(Neg_Worst_var$within_var)

gc()

#MODEL 1 - LV2 Predictor = DN, LV1 = individually mean-centered negative event indicator
#=-=-=-=-=-=-=-=-=-=-=-=-=-=
within_cov<-c(22)                    #Columns with group-mean centered predictors
between_cov<-c(8)                       #Columns with between-subject predictors
random_cov<-c(22)                    #Level 1 predictors with random slopes

#Getting distributions of fixed effects
Intercept = posterior_samples(Neg_Worst_c.DN, pars = 'b_Intercept') 
DN = posterior_samples(Neg_Worst_c.DN, pars = 'b_c.DN') 
Worst = posterior_samples(Neg_Worst_c.DN, pars = 'b_c.Worst') 

#Getting between-subjects variance terms for random effects
Int_var = posterior_samples(Neg_Worst_c.DN, pars = 'sd_ID__Intercept')^2
Worst_var = posterior_samples(Neg_Worst_c.DN, pars = 'sd_ID__c.Worst')^2 

#Getting covariances of random effects (i.e., tau matrix)
cov_Int_Worst = posterior_samples(Neg_Worst_c.DN, pars = 'sd_ID__Intercept')*
  posterior_samples(Neg_Worst_c.DN, pars = 'sd_ID__c.Worst')*
  posterior_samples(Neg_Worst_c.DN, pars = 'cor_ID__Intercept__c.Worst')

#Getting level 1 error variance
sigma<-log(1+posterior_samples(Neg_Worst_c.DN, pars = 'sigma')/beta_00)

post_samples<-data.frame(Intercept,
                         DN, 
                         Worst, 
                         Int_var, 
                         Worst_var, 
                         cov_Int_Worst, 
                         sigma)

colnames(post_samples)<-c('Intercept', 
                          'DN', 
                          'Worst', 
                          'Int_var', 
                          'Worst_var', 
                          'cov_Int_Worst', 
                          'sigma')

#Aggregate across imputed datasets
#Currently going to take 1000 draws from posterior distributions
#Will then apply across all 20 data sets

Neg_Worst_c.DN_var<-list(between_var=vector(), 
               within_var=vector(), 
               between_All_tot=vector(),
               between_All_btw=vector(),
               between_res_btw=vector(),
               within_fix_wthn=vector(),
               within_fix_tot=vector(),
               within_slope_var_wthn=vector(),
               within_res_wthn=vector(), 
               within_unmod_tot=vector())

for(i in 1:length(dat.study1.list)){
  print(paste('Starting w/ imputed set', i, 'out of', M))
  samp_frame<-1:post_iter+(i-1)*post_iter
  D<-sample(samp_frame, size = 400, replace = FALSE)
  for(d in 1:length(D)){
    #browser()
    Gamma_w<-as.vector(
      as.matrix(
        post_samples[D[d], c('Worst')]))       #Make sure the effects line up - in order of within_cov
    Gamma_b<-as.vector(
      as.matrix(
        post_samples[D[d], c('Intercept', 'DN')]))    #level fixed intercept and any level 2 fixed effects
    tau<-c(post_samples[D[d], c('Int_var', 'cov_Int_Worst')], 
           post_samples[D[d], c('cov_Int_Worst', 'Worst_var')])
    
    tau<-matrix(unlist(tau), 
                byrow=TRUE, 
                ncol=2)
    
    #Matrix rows/columns ordered starting with tau.00 (intercept variance)
    #then add columns/rows for random effects of slopes in order of within_cov
    
    Sigma2<-post_samples[D[d], 'sigma']
    
    r2mlm_DN<-r2MLM(data=dat.study1.list[[i]], 
                    within_covs = within_cov, 
                    between_covs = between_cov, 
                    random_covs = random_cov, 
                    gamma_w = Gamma_w, 
                    gamma_b = Gamma_b, 
                    Tau = tau, 
                    sigma2 = Sigma2, 
                    has_intercept = TRUE, 
                    clustermeancentered = TRUE)
    
    #Extracting relevant values across imputed datasets
    Neg_Worst_c.DN_var$between_var<-c(Neg_Worst_c.DN_var$between_var, as.numeric(r2mlm_DN$Decompositions['fixed, between', 'total'])+
                              as.numeric(r2mlm_DN$Decompositions['mean variation', 'total']))
    Neg_Worst_c.DN_var$within_var<-c(Neg_Worst_c.DN_var$within_var, as.numeric(r2mlm_DN$Decompositions['sigma2', 'total'])+
                             as.numeric(r2mlm_DN$Decompositions['slope variation', 'total'])+
                             as.numeric(r2mlm_DN$Decompositions['fixed, within', 'total']))
    Neg_Worst_c.DN_var$between_All_tot<-c(Neg_Worst_c.DN_var$between_All_tot, as.numeric(r2mlm_DN$R2s['f2', 'total']))
    Neg_Worst_c.DN_var$between_All_btw<-c(Neg_Worst_c.DN_var$between_All_btw, as.numeric(r2mlm_DN$R2s['f2', 'between']))
    Neg_Worst_c.DN_var$between_res_btw<-c(Neg_Worst_c.DN_var$between_res_btw, as.numeric(r2mlm_DN$R2s['m', 'between']))
    Neg_Worst_c.DN_var$within_fix_wthn<-c(Neg_Worst_c.DN_var$within_fix_wthn, as.numeric(r2mlm_DN$Decompositions['fixed, within', 'within']))
    Neg_Worst_c.DN_var$within_fix_tot<-c(Neg_Worst_c.DN_var$within_fix_tot, as.numeric(r2mlm_DN$Decompositions['fixed, within', 'total']))
    Neg_Worst_c.DN_var$within_slope_var_wthn<-c(Neg_Worst_c.DN_var$within_slope_var_wthn, as.numeric(r2mlm_DN$Decompositions['slope variation', 'within']))
    Neg_Worst_c.DN_var$within_res_wthn<-c(Neg_Worst_c.DN_var$within_res_wthn, as.numeric(r2mlm_DN$Decompositions['sigma2', 'within']))
    Neg_Worst_c.DN_var$within_unmod_tot<-c(Neg_Worst_c.DN_var$within_unmod_tot, as.numeric(r2mlm_DN$Decompositions['sigma2', 'total'])+
                                   as.numeric(r2mlm_DN$Decompositions['slope variation', 'total']))
    
    print(paste("Sampling", d, "out of", length(D), '- Dataset:', i))
  }
  print(paste("*** Finshed imputed data set", i, 'out of', M, "***"))
}
hist(Neg_Worst_c.DN_var$between_var)
hist(Neg_Worst_c.DN_var$within_var)

gc()

#MODEL 2 - LV2 Predictor = Mean Negative Events, LV1 = individually mean-centered negative event indicator
#=-=-=-=-=-=-=-=-=-=-=-=-=-=
within_cov<-c(22)                    #Columns with group-mean centered predictors
between_cov<-c(23)                       #Columns with between-subject predictors
random_cov<-c(22)                    #Level 1 predictors with random slopes

#Getting distributions of fixed effects
Intercept = posterior_samples(Neg_Worst_prop_Worst, pars = 'b_Intercept') 
prop_Worst = posterior_samples(Neg_Worst_prop_Worst, pars = 'b_prop_Worst') 
Worst = posterior_samples(Neg_Worst_prop_Worst, pars = 'b_c.Worst') 

#Getting between-subjects variance terms for random effects
Int_var = posterior_samples(Neg_Worst_prop_Worst, pars = 'sd_ID__Intercept')^2
Worst_var = posterior_samples(Neg_Worst_prop_Worst, pars = 'sd_ID__c.Worst')^2 

#Getting covariances of random effects (i.e., tau matrix)
cov_Int_Worst = posterior_samples(Neg_Worst_prop_Worst, pars = 'sd_ID__Intercept')*
  posterior_samples(Neg_Worst_prop_Worst, pars = 'sd_ID__c.Worst')*
  posterior_samples(Neg_Worst_prop_Worst, pars = 'cor_ID__Intercept__c.Worst')

#Getting level 1 error variance
sigma<-log(1+posterior_samples(Neg_Worst_prop_Worst, pars = 'sigma')/beta_00)

post_samples<-data.frame(Intercept,
                         prop_Worst, 
                         Worst, 
                         Int_var, 
                         Worst_var, 
                         cov_Int_Worst, 
                         sigma)

colnames(post_samples)<-c('Intercept', 
                          'prop_Worst', 
                          'Worst', 
                          'Int_var', 
                          'Worst_var', 
                          'cov_Int_Worst', 
                          'sigma')

#Aggregate across imputed datasets
#Currently going to take 1000 draws from posterior distributions
#Will then apply across all 20 data sets

Neg_Worst_prop_Worst_var<-list(between_var=vector(), 
                         within_var=vector(), 
                         between_All_tot=vector(),
                         between_All_btw=vector(),
                         between_res_btw=vector(),
                         within_fix_wthn=vector(),
                         within_fix_tot=vector(),
                         within_slope_var_wthn=vector(),
                         within_res_wthn=vector(), 
                         within_unmod_tot=vector())

for(i in 1:length(dat.study1.list)){
  print(paste('Starting w/ imputed set', i, 'out of', M))
  samp_frame<-1:post_iter+(i-1)*post_iter
  D<-sample(samp_frame, size = 400, replace = FALSE)
  for(d in 1:length(D)){
    #browser()
    Gamma_w<-as.vector(
      as.matrix(
        post_samples[D[d], c('Worst')]))       #Make sure the effects line up - in order of within_cov
    Gamma_b<-as.vector(
      as.matrix(
        post_samples[D[d], c('Intercept', 'prop_Worst')]))    #level fixed intercept and any level 2 fixed effects
    tau<-c(post_samples[D[d], c('Int_var', 'cov_Int_Worst')], 
           post_samples[D[d], c('cov_Int_Worst', 'Worst_var')])
    
    tau<-matrix(unlist(tau), 
                byrow=TRUE, 
                ncol=2)
    
    #Matrix rows/columns ordered starting with tau.00 (intercept variance)
    #then add columns/rows for random effects of slopes in order of within_cov
    
    Sigma2<-post_samples[D[d], 'sigma']
    
    r2mlm_prop_Worst<-r2MLM(data=dat.study1.list[[i]], 
                    within_covs = within_cov, 
                    between_covs = between_cov, 
                    random_covs = random_cov, 
                    gamma_w = Gamma_w, 
                    gamma_b = Gamma_b, 
                    Tau = tau, 
                    sigma2 = Sigma2, 
                    has_intercept = TRUE, 
                    clustermeancentered = TRUE)
    
    #Extracting relevant values across imputed datasets
    Neg_Worst_prop_Worst_var$between_var<-c(Neg_Worst_prop_Worst_var$between_var, as.numeric(r2mlm_prop_Worst$Decompositions['fixed, between', 'total'])+
                                        as.numeric(r2mlm_prop_Worst$Decompositions['mean variation', 'total']))
    Neg_Worst_prop_Worst_var$within_var<-c(Neg_Worst_prop_Worst_var$within_var, as.numeric(r2mlm_prop_Worst$Decompositions['sigma2', 'total'])+
                                       as.numeric(r2mlm_prop_Worst$Decompositions['slope variation', 'total'])+
                                       as.numeric(r2mlm_prop_Worst$Decompositions['fixed, within', 'total']))
    Neg_Worst_prop_Worst_var$between_All_tot<-c(Neg_Worst_prop_Worst_var$between_All_tot, as.numeric(r2mlm_prop_Worst$R2s['f2', 'total']))
    Neg_Worst_prop_Worst_var$between_All_btw<-c(Neg_Worst_prop_Worst_var$between_All_btw, as.numeric(r2mlm_prop_Worst$R2s['f2', 'between']))
    Neg_Worst_prop_Worst_var$between_res_btw<-c(Neg_Worst_prop_Worst_var$between_res_btw, as.numeric(r2mlm_prop_Worst$R2s['m', 'between']))
    Neg_Worst_prop_Worst_var$within_fix_wthn<-c(Neg_Worst_prop_Worst_var$within_fix_wthn, as.numeric(r2mlm_prop_Worst$Decompositions['fixed, within', 'within']))
    Neg_Worst_prop_Worst_var$within_fix_tot<-c(Neg_Worst_prop_Worst_var$within_fix_tot, as.numeric(r2mlm_prop_Worst$Decompositions['fixed, within', 'total']))
    Neg_Worst_prop_Worst_var$within_slope_var_wthn<-c(Neg_Worst_prop_Worst_var$within_slope_var_wthn, as.numeric(r2mlm_prop_Worst$Decompositions['slope variation', 'within']))
    Neg_Worst_prop_Worst_var$within_res_wthn<-c(Neg_Worst_prop_Worst_var$within_res_wthn, as.numeric(r2mlm_prop_Worst$Decompositions['sigma2', 'within']))
    Neg_Worst_prop_Worst_var$within_unmod_tot<-c(Neg_Worst_prop_Worst_var$within_unmod_tot, as.numeric(r2mlm_prop_Worst$Decompositions['sigma2', 'total'])+
                                             as.numeric(r2mlm_prop_Worst$Decompositions['slope variation', 'total']))
    
    print(paste("Sampling", d, "out of", length(D), '- Dataset:', i))
  }
  print(paste("*** Finshed imputed data set", i, 'out of', M, "***"))
}

hist(Neg_Worst_prop_Worst_var$between_var)
hist(Neg_Worst_prop_Worst_var$within_var)

#MODEL 3 - LV2 Predictors = DN and Mean Negative Events, LV1 = individually mean-centered negative event indicator
#=-=-=-=-=-=-=-=-=-=-=-=-=-=
within_cov<-c(22)                    #Columns with group-mean centered predictors
between_cov<-c(8,23)                       #Columns with between-subject predictors
random_cov<-c(22)                    #Level 1 predictors with random slopes

#Getting distributions of fixed effects
Intercept = posterior_samples(Neg_Worst_All, pars = 'b_Intercept') 
DN = posterior_samples(Neg_Worst_All, pars = 'b_c.DN') 
prop_Worst = posterior_samples(Neg_Worst_All, pars = 'b_prop_Worst') 
Worst = posterior_samples(Neg_Worst_All, pars = 'b_c.Worst') 

#Getting between-subjects variance terms for random effects
Int_var = posterior_samples(Neg_Worst_All, pars = 'sd_ID__Intercept')^2
Worst_var = posterior_samples(Neg_Worst_All, pars = 'sd_ID__c.Worst')^2 

#Getting covariances of random effects (i.e., tau matrix)
cov_Int_Worst = posterior_samples(Neg_Worst_All, pars = 'sd_ID__Intercept')*
  posterior_samples(Neg_Worst_All, pars = 'sd_ID__c.Worst')*
  posterior_samples(Neg_Worst_All, pars = 'cor_ID__Intercept__c.Worst')

#Getting level 1 error variance
sigma<-log(1+posterior_samples(Neg_Worst_All, pars = 'sigma')/beta_00)

post_samples<-data.frame(Intercept,
                         DN,
                         prop_Worst, 
                         Worst, 
                         Int_var, 
                         Worst_var, 
                         cov_Int_Worst, 
                         sigma)

colnames(post_samples)<-c('Intercept', 
                          'DN',
                          'prop_Worst', 
                          'Worst', 
                          'Int_var', 
                          'Worst_var', 
                          'cov_Int_Worst', 
                          'sigma')

#Aggregate across imputed datasets
#Currently going to take 1000 draws from posterior distributions
#Will then apply across all 20 data sets

Neg_Worst_All_var<-list(between_var=vector(), 
                               within_var=vector(), 
                               between_All_tot=vector(),
                               between_All_btw=vector(),
                               between_res_btw=vector(),
                               within_fix_wthn=vector(),
                               within_fix_tot=vector(),
                               within_slope_var_wthn=vector(),
                               within_res_wthn=vector(), 
                               within_unmod_tot=vector())

for(i in 1:length(dat.study1.list)){
  print(paste('Starting w/ imputed set', i, 'out of', M))
  samp_frame<-1:post_iter+(i-1)*post_iter
  D<-sample(samp_frame, size = 400, replace = FALSE)
  for(d in 1:length(D)){
    #browser()
    Gamma_w<-as.vector(
      as.matrix(
        post_samples[D[d], c('Worst')]))       #Make sure the effects line up - in order of within_cov
    Gamma_b<-as.vector(
      as.matrix(
        post_samples[D[d], c('Intercept', 'DN', 'prop_Worst')]))    #level fixed intercept and any level 2 fixed effects
    tau<-c(post_samples[D[d], c('Int_var', 'cov_Int_Worst')], 
           post_samples[D[d], c('cov_Int_Worst', 'Worst_var')])
    
    tau<-matrix(unlist(tau), 
                byrow=TRUE, 
                ncol=2)
    
    #Matrix rows/columns ordered starting with tau.00 (intercept variance)
    #then add columns/rows for random effects of slopes in order of within_cov
    
    Sigma2<-post_samples[D[d], 'sigma']
    
    r2mlm_prop_Worst<-r2MLM(data=dat.study1.list[[i]], 
                            within_covs = within_cov, 
                            between_covs = between_cov, 
                            random_covs = random_cov, 
                            gamma_w = Gamma_w, 
                            gamma_b = Gamma_b, 
                            Tau = tau, 
                            sigma2 = Sigma2, 
                            has_intercept = TRUE, 
                            clustermeancentered = TRUE)
    
    #Extracting relevant values across imputed datasets
    Neg_Worst_All_var$between_var<-c(Neg_Worst_All_var$between_var, as.numeric(r2mlm_prop_Worst$Decompositions['fixed, between', 'total'])+
                                              as.numeric(r2mlm_prop_Worst$Decompositions['mean variation', 'total']))
    Neg_Worst_All_var$within_var<-c(Neg_Worst_All_var$within_var, as.numeric(r2mlm_prop_Worst$Decompositions['sigma2', 'total'])+
                                             as.numeric(r2mlm_prop_Worst$Decompositions['slope variation', 'total'])+
                                             as.numeric(r2mlm_prop_Worst$Decompositions['fixed, within', 'total']))
    Neg_Worst_All_var$between_All_tot<-c(Neg_Worst_All_var$between_All_tot, as.numeric(r2mlm_prop_Worst$R2s['f2', 'total']))
    Neg_Worst_All_var$between_All_btw<-c(Neg_Worst_All_var$between_All_btw, as.numeric(r2mlm_prop_Worst$R2s['f2', 'between']))
    Neg_Worst_All_var$between_res_btw<-c(Neg_Worst_All_var$between_res_btw, as.numeric(r2mlm_prop_Worst$R2s['m', 'between']))
    Neg_Worst_All_var$within_fix_wthn<-c(Neg_Worst_All_var$within_fix_wthn, as.numeric(r2mlm_prop_Worst$Decompositions['fixed, within', 'within']))
    Neg_Worst_All_var$within_fix_tot<-c(Neg_Worst_All_var$within_fix_tot, as.numeric(r2mlm_prop_Worst$Decompositions['fixed, within', 'total']))
    Neg_Worst_All_var$within_slope_var_wthn<-c(Neg_Worst_All_var$within_slope_var_wthn, as.numeric(r2mlm_prop_Worst$Decompositions['slope variation', 'within']))
    Neg_Worst_All_var$within_res_wthn<-c(Neg_Worst_All_var$within_res_wthn, as.numeric(r2mlm_prop_Worst$Decompositions['sigma2', 'within']))
    Neg_Worst_All_var$within_unmod_tot<-c(Neg_Worst_All_var$within_unmod_tot, as.numeric(r2mlm_prop_Worst$Decompositions['sigma2', 'total'])+
                                                   as.numeric(r2mlm_prop_Worst$Decompositions['slope variation', 'total']))
    
    print(paste("Sampling", d, "out of", length(D), '- Dataset:', i))
  }
  print(paste("*** Finshed imputed data set", i, 'out of', M, "***"))
}
hist(Neg_Worst_All_var$between_var)
hist(Neg_Worst_All_var$within_var)


#MODEL 4 - Full cross-level interaction with 
#=-=-=-=-=-=-=-=-=-=-=-=-=-=
within_cov<-c(22,26)                    #Columns with group-mean centered predictors
between_cov<-c(8,23)                       #Columns with between-subject predictors
random_cov<-c(22)                    #Level 1 predictors with random slopes

#Getting distributions of fixed effects
Intercept = posterior_samples(Neg_Worst_cross, pars = 'b_Intercept') 
DN = posterior_samples(Neg_Worst_cross, pars = 'b_c.DN') 
prop_Worst = posterior_samples(Neg_Worst_cross, pars = 'b_prop_Worst') 
Worst = posterior_samples(Neg_Worst_cross, pars = 'b_c.Worst')[,1]  #Note posterior_samples() searches for text match
                                                                    #With cross-level interaction returns two values 
                                                                    #Need to select the one with the fixed effect for lv-1
DNxWorst = posterior_samples(Neg_Worst_cross, pars = 'b_c.Worst:c.DN')

#Getting between-subjects variance terms for random effects
Int_var = posterior_samples(Neg_Worst_cross, pars = 'sd_ID__Intercept')^2
Worst_var = posterior_samples(Neg_Worst_cross, pars = 'sd_ID__c.Worst')^2 

#Getting covariances of random effects (i.e., tau matrix)
cov_Int_Worst = posterior_samples(Neg_Worst_cross, pars = 'sd_ID__Intercept')*
  posterior_samples(Neg_Worst_cross, pars = 'sd_ID__c.Worst')*
  posterior_samples(Neg_Worst_cross, pars = 'cor_ID__Intercept__c.Worst')

#Getting level 1 error variance
sigma<-log(1+posterior_samples(Neg_Worst_cross, pars = 'sigma')/beta_00)

post_samples<-data.frame(Intercept,
                         DN,
                         prop_Worst, 
                         Worst,
                         DNxWorst,
                         Int_var, 
                         Worst_var, 
                         cov_Int_Worst, 
                         sigma)

colnames(post_samples)<-c('Intercept', 
                          'DN',
                          'prop_Worst', 
                          'Worst', 
                          'DNxWorst',
                          'Int_var', 
                          'Worst_var', 
                          'cov_Int_Worst', 
                          'sigma')

#Aggregate across imputed datasets
#Currently going to take 1000 draws from posterior distributions
#Will then apply across all 20 data sets

Neg_Worst_cross_var<-list(between_var=vector(), 
                        within_var=vector(), 
                        between_All_tot=vector(),
                        between_All_btw=vector(),
                        between_res_btw=vector(),
                        within_fix_wthn=vector(),
                        within_fix_tot=vector(),
                        within_slope_var_wthn=vector(),
                        within_res_wthn=vector(), 
                        within_unmod_tot=vector())

for(i in 1:length(dat.study1.list)){
  print(paste('Starting w/ imputed set', i, 'out of', M))
  samp_frame<-1:post_iter+(i-1)*post_iter
  D<-sample(samp_frame, size = 400, replace = FALSE)
  for(d in 1:length(D)){
    #browser()
    Gamma_w<-as.vector(
      as.matrix(
        post_samples[D[d], c('Worst', 'DNxWorst')]))       #Make sure the effects line up - in order of within_cov
    Gamma_b<-as.vector(
      as.matrix(
        post_samples[D[d], c('Intercept', 'DN', 'prop_Worst')]))    #level fixed intercept and any level 2 fixed effects
    tau<-c(post_samples[D[d], c('Int_var', 'cov_Int_Worst')], 
           post_samples[D[d], c('cov_Int_Worst', 'Worst_var')])
    
    tau<-matrix(unlist(tau), 
                byrow=TRUE, 
                ncol=2)
    
    #Matrix rows/columns ordered starting with tau.00 (intercept variance)
    #then add columns/rows for random effects of slopes in order of within_cov
    
    Sigma2<-post_samples[D[d], 'sigma']
    
    r2mlm_prop_Worst<-r2MLM(data=dat.study1.list[[i]], 
                            within_covs = within_cov, 
                            between_covs = between_cov, 
                            random_covs = random_cov, 
                            gamma_w = Gamma_w, 
                            gamma_b = Gamma_b, 
                            Tau = tau, 
                            sigma2 = Sigma2, 
                            has_intercept = TRUE, 
                            clustermeancentered = TRUE)
    
    #Extracting relevant values across imputed datasets
    Neg_Worst_cross_var$between_var<-c(Neg_Worst_cross_var$between_var, as.numeric(r2mlm_prop_Worst$Decompositions['fixed, between', 'total'])+
                                       as.numeric(r2mlm_prop_Worst$Decompositions['mean variation', 'total']))
    Neg_Worst_cross_var$within_var<-c(Neg_Worst_cross_var$within_var, as.numeric(r2mlm_prop_Worst$Decompositions['sigma2', 'total'])+
                                      as.numeric(r2mlm_prop_Worst$Decompositions['slope variation', 'total'])+
                                      as.numeric(r2mlm_prop_Worst$Decompositions['fixed, within', 'total']))
    Neg_Worst_cross_var$between_All_tot<-c(Neg_Worst_cross_var$between_All_tot, as.numeric(r2mlm_prop_Worst$R2s['f2', 'total']))
    Neg_Worst_cross_var$between_All_btw<-c(Neg_Worst_cross_var$between_All_btw, as.numeric(r2mlm_prop_Worst$R2s['f2', 'between']))
    Neg_Worst_cross_var$between_res_btw<-c(Neg_Worst_cross_var$between_res_btw, as.numeric(r2mlm_prop_Worst$R2s['m', 'between']))
    Neg_Worst_cross_var$within_fix_wthn<-c(Neg_Worst_cross_var$within_fix_wthn, as.numeric(r2mlm_prop_Worst$Decompositions['fixed, within', 'within']))
    Neg_Worst_cross_var$within_fix_tot<-c(Neg_Worst_cross_var$within_fix_tot, as.numeric(r2mlm_prop_Worst$Decompositions['fixed, within', 'total']))
    Neg_Worst_cross_var$within_slope_var_wthn<-c(Neg_Worst_cross_var$within_slope_var_wthn, as.numeric(r2mlm_prop_Worst$Decompositions['slope variation', 'within']))
    Neg_Worst_cross_var$within_res_wthn<-c(Neg_Worst_cross_var$within_res_wthn, as.numeric(r2mlm_prop_Worst$Decompositions['sigma2', 'within']))
    Neg_Worst_cross_var$within_unmod_tot<-c(Neg_Worst_cross_var$within_unmod_tot, as.numeric(r2mlm_prop_Worst$Decompositions['sigma2', 'total'])+
                                            as.numeric(r2mlm_prop_Worst$Decompositions['slope variation', 'total']))
    
    print(paste("Sampling", d, "out of", length(D), '- Dataset:', i))
  }
  print(paste("*** Finshed imputed data set", i, 'out of', M, "***"))
}
hist(Neg_Worst_cross_var$between_var)
hist(Neg_Worst_cross_var$within_var)

###########################################################################################################
#Preparing values for River plots :: Negative Mood - Worst Event Models ::
###########################################################################################################
#Checking across models - between subjects variance
mean(Neg_Worst_c.DN_var$between_var)
mean(Neg_Worst_prop_Worst_var$between_var)
mean(Neg_Worst_All_var$between_var)
mean(Neg_Worst_cross_var$between_var)

#Differemce between R2's due to sampling variability - is only 0.0008066132
#Difference is <.1% in total variance
mean(Neg_Worst_cross_var$between_var)-mean(Neg_Worst_c.DN_var$between_var)

#Between Variance Model Terms: 
Neg_Worst_btw<-mean(Neg_Worst_var$between_var)
Neg_Worst_btw_All<-mean(Neg_Worst_cross_var$between_All_tot)
Neg_Worst_btw_umod<-Neg_Worst_btw-Neg_Worst_btw_All

Neg_Worst_btw_DN_unique<-(mean(Neg_Worst_All_var$between_All_btw)-mean(Neg_Worst_prop_Worst_var$between_All_btw))*Neg_Worst_btw
Neg_Worst_btw_Exp_unique<-(mean(Neg_Worst_All_var$between_All_btw)-mean(Neg_Worst_c.DN_var$between_All_btw))*Neg_Worst_btw

Neg_Worst_btw_Shared<-Neg_Worst_btw_All-(Neg_Worst_btw_DN_unique+Neg_Worst_btw_Exp_unique)

#Within Variance Model Terms: 
Neg_Worst_wthn<-mean(Neg_Worst_var$within_var)
Neg_Worst_wthn_Worst<-mean(Neg_Worst_All_var$within_fix_tot)
Neg_Worst_wthn_Worst_DN<-(mean(Neg_Worst_All_var$within_res_wthn)-mean(Neg_Worst_cross_var$within_res_wthn))*Neg_Worst_wthn
Neg_Worst_wthn_unmod<-Neg_Worst_wthn-Neg_Worst_wthn_Worst-Neg_Worst_wthn_Worst_DN
#Note final value produced above is a measure of total variance attributable to cross-level interaction

#Creating Total DN Effect for separate river plot decomposition: 
Total_DN<-Neg_Worst_btw_DN_unique+Neg_Worst_btw_Shared+Neg_Worst_wthn_Worst_DN

###########################################################################################################
#Main Variance Decompisition River Plot :: Negative Mood - Worst Event Models :: 
###########################################################################################################
Neg_Worst_River_DF<-data.frame(N1 = c('DN',
                                     'DN <--> NDE \n Exposure',
                                     'NDE \n Exposure', 
                                     'Unmodeled \n Between', 
                                     'Momentary \n NDE', 
                                     'DN x NDE \n Reactivity',
                                     'Unmodeled \n Within', 
                                     'Total \n Between', 
                                     'Total \n Within'), 
                              N2 = c('Total \n Between', 
                                     'Total \n Between', 
                                     'Total \n Between', 
                                     'Total \n Between', 
                                     'Total \n Within', 
                                     'Total \n Within',
                                     'Total \n Within', 
                                     'Total \n Variance', 
                                     'Total \n Variance'), 
                              Value = c(Neg_Worst_btw_DN_unique, 
                                        Neg_Worst_btw_Shared,
                                        Neg_Worst_btw_Exp_unique, 
                                        Neg_Worst_btw_umod,
                                        Neg_Worst_wthn_Worst, 
                                        Neg_Worst_wthn_Worst_DN,
                                        Neg_Worst_wthn_unmod,
                                        Neg_Worst_btw, 
                                        Neg_Worst_wthn), 
                              ID = 1:9)

Neg_Worst_River_DF$N1<-paste(Neg_Worst_River_DF$N1, 
                            '\n', 
                            paste0(round(Neg_Worst_River_DF$Value*100, 
                                         digits = 2), '%'
                            )
)

Neg_Worst_River_DF$N2<-c(rep(Neg_Worst_River_DF$N1[8], 4),
                        rep(Neg_Worst_River_DF$N1[9], 3), 
                        rep('Total Variance', 2)
)

nodes<-data.frame(ID = c(Neg_Worst_River_DF$N1, 
                         'Total Variance'), 
                  x = c(1,1,1,1,1,1,1,2,2,3), 
                  y = c(0,2,4,6,8,10,12,3,9,6))

palette = c(paste0(brewer.pal(9, "Blues"), 90)[c(6, 6, 6, 6)], 
            paste0(brewer.pal(9, "Reds"), 90)[c(6,6,6)], 
            paste0(brewer.pal(9, "Blues"), 90)[9], 
            paste0(brewer.pal(9, "Reds"), 90)[9], 
            paste0(brewer.pal(9, 'Purples'), 90)[9])

styles = lapply(nodes$y, function(n) {
  list(col = palette[n], lty = 0, textcol = "black")
})

#Not sure why color is not mapping correctly - trying to force it with this fuction... 
for(i in 1:length(palette)){
  styles[[i]]$col<-palette[i]
}

names(styles) = nodes$ID

Neg_Worst_riv<-makeRiver(nodes = nodes, 
                        edges =  Neg_Worst_River_DF, 
                        node_styles = styles)

png(paste0(study1.graphics, '/S1_Neg_Worst_river_overall.png'), 
    units = 'in', 
    res = 1200, 
    height = 14, 
    width = 8)
riverplot(Neg_Worst_riv, 
          nodewidth = 2, 
          plot_area = .95)
title(ylab = 'Study 1 - Riverplot of Total Variance Decomposition Estimated from Negative Mood Models: Negative Daily Events')
dev.off()

#----------------------------------------------------------------------------------------------------------
#DN - Specific Plot 
# :: Negative Mood - Worst Event Models :: 
#----------------------------------------------------------------------------------------------------------
Neg_Worst_DN_River_DF<-data.frame(N1 = c('DN',
                                        'DN <--> NDE \n Exposure',
                                        'DN x NDE \n Reactivity'
), 
N2 = rep(paste0("Combined DN Effect ",
                round(Total_DN*100, digits = 2),
                '%'), 
         3), 
Value = c(Neg_Worst_btw_DN_unique, 
          Neg_Worst_btw_Shared, 
          Neg_Worst_wthn_Worst_DN), 
ID = 1:3)

Neg_Worst_DN_River_DF$N1<-paste(Neg_Worst_DN_River_DF$N1, 
                               '\n', 
                               paste0(round(Neg_Worst_DN_River_DF$Value*100, 
                                            digits = 2), '%'
                               )
)

nodes<-data.frame(ID = c(Neg_Worst_DN_River_DF$N1, 
                         paste0("Combined DN Effect ",
                                round(Total_DN*100, digits = 2),
                                '%')), 
                  x = c(1,1,1,2), 
                  y = c(0,1,2,1))

palette = c(paste0(brewer.pal(9, "Blues"),90)[c(6,9)],
            paste0(brewer.pal(9, "Reds"),90)[6], 
            paste0(brewer.pal(9, "Purples"),90)[9])

styles = lapply(nodes$y, function(n) {
  list(col = palette[n], lty = 0, textcol = "black")
})

#Not sure why color is not mapping correctly - trying to force it with this fuction... 
for(i in 1:length(palette)){
  styles[[i]]$col<-palette[i]
}

names(styles) = nodes$ID

Neg_Worst_DN_riv<-makeRiver(nodes = nodes, 
                           edges =  Neg_Worst_DN_River_DF, 
                           node_styles = styles)

png(paste0(study1.graphics, '/S1_Neg_Worst_river_DN.png'), 
    units = 'in', 
    res = 900, 
    height = 10, 
    width = 10)
riverplot(Neg_Worst_DN_riv, 
          nodewidth = 2, 
          plot_area = .95)
title(ylab = 'Study 1 - Riverplot of DN Effect Decomposition Derived from Negative Mood Model: Negative Daily Events')
dev.off()

#Come back to this area to produce publication graphics (or arrange as graphics in ggplot)
#png(paste0(study1.graphics, '/Combined_Neg_Worst_Rivers.png'), 
png(paste0(getwd(), '/Combined_Neg_Worst_Rivers.png'),
    units = 'in', 
    res = 900, 
    height = 10, 
    width = 6)
layout(matrix(c(1,2,2), nrow=3, ncol=1))
par(mar = c(0,0,0,0))

riverplot(Neg_Worst_DN_riv, 
          nodewidth = 2, 
          plot_area = .85)
title(ylab = 'Study 1 - Riverplot of DN Effect Decomposition Derived from Negative Mood Model: Negative Daily Events')

riverplot(Neg_Worst_riv, 
          nodewidth = 2, 
          plot_area = .85)
title(ylab = 'Study 1 - Riverplot of Total Variance Decomposition Estimated from Negative Mood Models: Negative Daily Events')
dev.off()

#==========================================================================================================
#NEGATIVE MOOD MODELS - Variance Decomposition
# :: Negative Mood - Best Event Models ::
#==========================================================================================================

#==========================================================================================================
#----------------------------------------------------------------------------------------------------------
#Negative Mood, Negative Events, and Dispositional Negativity: 
#==========================================================================================================
#Will need the grand intercept from the null model (on the link scale)
beta_00<-fixef(Neg_ucm)[1,1]
post_iter<-2000   #No. of chains per data set x No. post-warmup iterations per chain
#for loops below require this value is correctly specified

#NULL model - Using simple level 1 to ground between/within estimates
#=-=-=-=-=-=-=-=-=-=-=-=-=-=
within_cov<-c(24)                    #Columns with group-mean centered predictors
random_cov<-c(24)                    #Level 1 predictors with random slopes

#Getting distributions of fixed effects
Intercept = posterior_samples(Neg_Best, pars = 'b_Intercept') 
Best = posterior_samples(Neg_Best, pars = 'b_c.Best') 

#Getting between-subjects variance terms for random effects
Int_var = posterior_samples(Neg_Best, pars = 'sd_ID__Intercept')^2
Best_var = posterior_samples(Neg_Best, pars = 'sd_ID__c.Best')^2 

#Getting covariances of random effects (i.e., tau matrix)
cov_Int_Best = posterior_samples(Neg_Best, pars = 'sd_ID__Intercept')*
  posterior_samples(Neg_Best, pars = 'sd_ID__c.Best')*
  posterior_samples(Neg_Best, pars = 'cor_ID__Intercept__c.Best')

#Getting level 1 error variance
sigma<-log(1+posterior_samples(Neg_Best, pars = 'sigma')/beta_00)

post_samples<-data.frame(Intercept,
                         Best, 
                         Int_var, 
                         Best_var, 
                         cov_Int_Best, 
                         sigma)

colnames(post_samples)<-c('Intercept', 
                          'Best', 
                          'Int_var', 
                          'Best_var', 
                          'cov_Int_Best', 
                          'sigma')

#Aggregate across imputed datasets
#Currently going to take 1000 draws from posterior distributions
#Will then apply across all 20 data sets

Neg_Best_var<-list(between_var=vector(), 
                        within_var=vector(), 
                        between_All_tot=vector(),
                        between_All_btw=vector(),
                        between_res_btw=vector(),
                        within_fix_wthn=vector(),
                        within_fix_tot=vector(),
                        within_slope_var_wthn=vector(),
                        within_res_wthn=vector(), 
                        within_unmod_tot=vector())

for(i in 1:length(dat.study1.list)){
  print(paste('Starting w/ imputed set', i, 'out of', M))
  samp_frame<-1:post_iter+(i-1)*post_iter
  D<-sample(samp_frame, size = 400, replace = FALSE)
  for(d in 1:length(D)){
    #browser()
    Gamma_w<-as.vector(
      as.matrix(
        post_samples[D[d], c('Best')]))       #Make sure the effects line up - in order of within_cov
    Gamma_b<-as.vector(
      as.matrix(
        post_samples[D[d], c('Intercept')]))    #level fixed intercept and any level 2 fixed effects
    tau<-c(post_samples[D[d], c('Int_var', 'cov_Int_Best')], 
           post_samples[D[d], c('cov_Int_Best', 'Best_var')])
    
    tau<-matrix(unlist(tau), 
                byrow=TRUE, 
                ncol=2)
    
    #Matrix rows/columns ordered starting with tau.00 (intercept variance)
    #then add columns/rows for random effects of slopes in order of within_cov
    
    Sigma2<-post_samples[D[d], 'sigma']
    
    r2mlm_DN<-r2MLM(data=dat.study1.list[[i]], 
                    within_covs = within_cov, 
                    between_covs = NULL, 
                    random_covs = random_cov, 
                    gamma_w = Gamma_w, 
                    gamma_b = Gamma_b, 
                    Tau = tau, 
                    sigma2 = Sigma2, 
                    has_intercept = TRUE, 
                    clustermeancentered = TRUE)
    
    #Extracting relevant values across imputed datasets
    Neg_Best_var$between_var<-c(Neg_Best_var$between_var, as.numeric(r2mlm_DN$Decompositions['fixed, between', 'total'])+
                                       as.numeric(r2mlm_DN$Decompositions['mean variation', 'total']))
    Neg_Best_var$within_var<-c(Neg_Best_var$within_var, as.numeric(r2mlm_DN$Decompositions['sigma2', 'total'])+
                                      as.numeric(r2mlm_DN$Decompositions['slope variation', 'total'])+
                                      as.numeric(r2mlm_DN$Decompositions['fixed, within', 'total']))
    Neg_Best_var$between_All_tot<-c(Neg_Best_var$between_All_tot, as.numeric(r2mlm_DN$R2s['f2', 'total']))
    Neg_Best_var$between_All_btw<-c(Neg_Best_var$between_All_btw, as.numeric(r2mlm_DN$R2s['f2', 'between']))
    Neg_Best_var$between_res_btw<-c(Neg_Best_var$between_res_btw, as.numeric(r2mlm_DN$R2s['m', 'between']))
    Neg_Best_var$within_fix_wthn<-c(Neg_Best_var$within_fix_wthn, as.numeric(r2mlm_DN$Decompositions['fixed, within', 'within']))
    Neg_Best_var$within_fix_tot<-c(Neg_Best_var$within_fix_tot, as.numeric(r2mlm_DN$Decompositions['fixed, within', 'total']))
    Neg_Best_var$within_slope_var_wthn<-c(Neg_Best_var$within_slope_var_wthn, as.numeric(r2mlm_DN$Decompositions['slope variation', 'within']))
    Neg_Best_var$within_res_wthn<-c(Neg_Best_var$within_res_wthn, as.numeric(r2mlm_DN$Decompositions['sigma2', 'within']))
    Neg_Best_var$within_unmod_tot<-c(Neg_Best_var$within_unmod_tot, as.numeric(r2mlm_DN$Decompositions['sigma2', 'total'])+
                                            as.numeric(r2mlm_DN$Decompositions['slope variation', 'total']))
    
    print(paste("Sampling", d, "out of", length(D), '- Dataset:', i))
  }
  print(paste("*** Finshed imputed data set", i, 'out of', M, "***"))
}

gc()

hist(Neg_Best_var$between_var)
hist(Neg_Best_var$within_var)


#MODEL 1 - LV2 Predictor = DN, LV1 = individually mean-centered negative event indicator
#=-=-=-=-=-=-=-=-=-=-=-=-=-=
within_cov<-c(24)                    #Columns with group-mean centered predictors
between_cov<-c(8)                       #Columns with between-subject predictors
random_cov<-c(24)                    #Level 1 predictors with random slopes

#Getting distributions of fixed effects
Intercept = posterior_samples(Neg_Best_c.DN, pars = 'b_Intercept') 
DN = posterior_samples(Neg_Best_c.DN, pars = 'b_c.DN') 
Best = posterior_samples(Neg_Best_c.DN, pars = 'b_c.Best') 

#Getting between-subjects variance terms for random effects
Int_var = posterior_samples(Neg_Best_c.DN, pars = 'sd_ID__Intercept')^2
Best_var = posterior_samples(Neg_Best_c.DN, pars = 'sd_ID__c.Best')^2 

#Getting covariances of random effects (i.e., tau matrix)
cov_Int_Best = posterior_samples(Neg_Best_c.DN, pars = 'sd_ID__Intercept')*
  posterior_samples(Neg_Best_c.DN, pars = 'sd_ID__c.Best')*
  posterior_samples(Neg_Best_c.DN, pars = 'cor_ID__Intercept__c.Best')

#Getting level 1 error variance
sigma<-log(1+posterior_samples(Neg_Best_c.DN, pars = 'sigma')/beta_00)

post_samples<-data.frame(Intercept,
                         DN, 
                         Best, 
                         Int_var, 
                         Best_var, 
                         cov_Int_Best, 
                         sigma)

colnames(post_samples)<-c('Intercept', 
                          'DN', 
                          'Best', 
                          'Int_var', 
                          'Best_var', 
                          'cov_Int_Best', 
                          'sigma')

#Aggregate across imputed datasets
#Currently going to take 1000 draws from posterior distributions
#Will then apply across all 20 data sets

Neg_Best_c.DN_var<-list(between_var=vector(), 
                         within_var=vector(), 
                         between_All_tot=vector(),
                         between_All_btw=vector(),
                         between_res_btw=vector(),
                         within_fix_wthn=vector(),
                         within_fix_tot=vector(),
                         within_slope_var_wthn=vector(),
                         within_res_wthn=vector(), 
                         within_unmod_tot=vector())

for(i in 1:length(dat.study1.list)){
  print(paste('Starting w/ imputed set', i, 'out of', M))
  samp_frame<-1:post_iter+(i-1)*post_iter
  D<-sample(samp_frame, size = 400, replace = FALSE)
  for(d in 1:length(D)){
    #browser()
    Gamma_w<-as.vector(
      as.matrix(
        post_samples[D[d], c('Best')]))       #Make sure the effects line up - in order of within_cov
    Gamma_b<-as.vector(
      as.matrix(
        post_samples[D[d], c('Intercept', 'DN')]))    #level fixed intercept and any level 2 fixed effects
    tau<-c(post_samples[D[d], c('Int_var', 'cov_Int_Best')], 
           post_samples[D[d], c('cov_Int_Best', 'Best_var')])
    
    tau<-matrix(unlist(tau), 
                byrow=TRUE, 
                ncol=2)
    
    #Matrix rows/columns ordered starting with tau.00 (intercept variance)
    #then add columns/rows for random effects of slopes in order of within_cov
    
    Sigma2<-post_samples[D[d], 'sigma']
    
    r2mlm_DN<-r2MLM(data=dat.study1.list[[i]], 
                    within_covs = within_cov, 
                    between_covs = between_cov, 
                    random_covs = random_cov, 
                    gamma_w = Gamma_w, 
                    gamma_b = Gamma_b, 
                    Tau = tau, 
                    sigma2 = Sigma2, 
                    has_intercept = TRUE, 
                    clustermeancentered = TRUE)
    
    #Extracting relevant values across imputed datasets
    Neg_Best_c.DN_var$between_var<-c(Neg_Best_c.DN_var$between_var, as.numeric(r2mlm_DN$Decompositions['fixed, between', 'total'])+
                                        as.numeric(r2mlm_DN$Decompositions['mean variation', 'total']))
    Neg_Best_c.DN_var$within_var<-c(Neg_Best_c.DN_var$within_var, as.numeric(r2mlm_DN$Decompositions['sigma2', 'total'])+
                                       as.numeric(r2mlm_DN$Decompositions['slope variation', 'total'])+
                                       as.numeric(r2mlm_DN$Decompositions['fixed, within', 'total']))
    Neg_Best_c.DN_var$between_All_tot<-c(Neg_Best_c.DN_var$between_All_tot, as.numeric(r2mlm_DN$R2s['f2', 'total']))
    Neg_Best_c.DN_var$between_All_btw<-c(Neg_Best_c.DN_var$between_All_btw, as.numeric(r2mlm_DN$R2s['f2', 'between']))
    Neg_Best_c.DN_var$between_res_btw<-c(Neg_Best_c.DN_var$between_res_btw, as.numeric(r2mlm_DN$R2s['m', 'between']))
    Neg_Best_c.DN_var$within_fix_wthn<-c(Neg_Best_c.DN_var$within_fix_wthn, as.numeric(r2mlm_DN$Decompositions['fixed, within', 'within']))
    Neg_Best_c.DN_var$within_fix_tot<-c(Neg_Best_c.DN_var$within_fix_tot, as.numeric(r2mlm_DN$Decompositions['fixed, within', 'total']))
    Neg_Best_c.DN_var$within_slope_var_wthn<-c(Neg_Best_c.DN_var$within_slope_var_wthn, as.numeric(r2mlm_DN$Decompositions['slope variation', 'within']))
    Neg_Best_c.DN_var$within_res_wthn<-c(Neg_Best_c.DN_var$within_res_wthn, as.numeric(r2mlm_DN$Decompositions['sigma2', 'within']))
    Neg_Best_c.DN_var$within_unmod_tot<-c(Neg_Best_c.DN_var$within_unmod_tot, as.numeric(r2mlm_DN$Decompositions['sigma2', 'total'])+
                                             as.numeric(r2mlm_DN$Decompositions['slope variation', 'total']))
    
    print(paste("Sampling", d, "out of", length(D), '- Dataset:', i))
  }
  print(paste("*** Finshed imputed data set", i, 'out of', M, "***"))
}

gc()

hist(Neg_Best_c.DN_var$between_var)
hist(Neg_Best_c.DN_var$within_var)

#MODEL 2 - LV2 Predictor = Mean Negative Events, LV1 = individually mean-centered negative event indicator
#=-=-=-=-=-=-=-=-=-=-=-=-=-=
within_cov<-c(24)                    #Columns with group-mean centered predictors
between_cov<-c(25)                       #Columns with between-subject predictors
random_cov<-c(24)                    #Level 1 predictors with random slopes

#Getting distributions of fixed effects
Intercept = posterior_samples(Neg_Best_prop_Best, pars = 'b_Intercept') 
prop_Best = posterior_samples(Neg_Best_prop_Best, pars = 'b_prop_Best') 
Best = posterior_samples(Neg_Best_prop_Best, pars = 'b_c.Best') 

#Getting between-subjects variance terms for random effects
Int_var = posterior_samples(Neg_Best_prop_Best, pars = 'sd_ID__Intercept')^2
Best_var = posterior_samples(Neg_Best_prop_Best, pars = 'sd_ID__c.Best')^2 

#Getting covariances of random effects (i.e., tau matrix)
cov_Int_Best = posterior_samples(Neg_Best_prop_Best, pars = 'sd_ID__Intercept')*
  posterior_samples(Neg_Best_prop_Best, pars = 'sd_ID__c.Best')*
  posterior_samples(Neg_Best_prop_Best, pars = 'cor_ID__Intercept__c.Best')

#Getting level 1 error variance
sigma<-log(1+posterior_samples(Neg_Best_prop_Best, pars = 'sigma')/beta_00)

post_samples<-data.frame(Intercept,
                         prop_Best, 
                         Best, 
                         Int_var, 
                         Best_var, 
                         cov_Int_Best, 
                         sigma)

colnames(post_samples)<-c('Intercept', 
                          'prop_Best', 
                          'Best', 
                          'Int_var', 
                          'Best_var', 
                          'cov_Int_Best', 
                          'sigma')

#Aggregate across imputed datasets
#Currently going to take 1000 draws from posterior distributions
#Will then apply across all 20 data sets

Neg_Best_prop_Best_var<-list(between_var=vector(), 
                               within_var=vector(), 
                               between_All_tot=vector(),
                               between_All_btw=vector(),
                               between_res_btw=vector(),
                               within_fix_wthn=vector(),
                               within_fix_tot=vector(),
                               within_slope_var_wthn=vector(),
                               within_res_wthn=vector(), 
                               within_unmod_tot=vector())

for(i in 1:length(dat.study1.list)){
  print(paste('Starting w/ imputed set', i, 'out of', M))
  samp_frame<-1:post_iter+(i-1)*post_iter
  D<-sample(samp_frame, size = 400, replace = FALSE)
  for(d in 1:length(D)){
    #browser()
    Gamma_w<-as.vector(
      as.matrix(
        post_samples[D[d], c('Best')]))       #Make sure the effects line up - in order of within_cov
    Gamma_b<-as.vector(
      as.matrix(
        post_samples[D[d], c('Intercept', 'prop_Best')]))    #level fixed intercept and any level 2 fixed effects
    tau<-c(post_samples[D[d], c('Int_var', 'cov_Int_Best')], 
           post_samples[D[d], c('cov_Int_Best', 'Best_var')])
    
    tau<-matrix(unlist(tau), 
                byrow=TRUE, 
                ncol=2)
    
    #Matrix rows/columns ordered starting with tau.00 (intercept variance)
    #then add columns/rows for random effects of slopes in order of within_cov
    
    Sigma2<-post_samples[D[d], 'sigma']
    
    r2mlm_prop_Best<-r2MLM(data=dat.study1.list[[i]], 
                            within_covs = within_cov, 
                            between_covs = between_cov, 
                            random_covs = random_cov, 
                            gamma_w = Gamma_w, 
                            gamma_b = Gamma_b, 
                            Tau = tau, 
                            sigma2 = Sigma2, 
                            has_intercept = TRUE, 
                            clustermeancentered = TRUE)
    
    #Extracting relevant values across imputed datasets
    Neg_Best_prop_Best_var$between_var<-c(Neg_Best_prop_Best_var$between_var, as.numeric(r2mlm_prop_Best$Decompositions['fixed, between', 'total'])+
                                              as.numeric(r2mlm_prop_Best$Decompositions['mean variation', 'total']))
    Neg_Best_prop_Best_var$within_var<-c(Neg_Best_prop_Best_var$within_var, as.numeric(r2mlm_prop_Best$Decompositions['sigma2', 'total'])+
                                             as.numeric(r2mlm_prop_Best$Decompositions['slope variation', 'total'])+
                                             as.numeric(r2mlm_prop_Best$Decompositions['fixed, within', 'total']))
    Neg_Best_prop_Best_var$between_All_tot<-c(Neg_Best_prop_Best_var$between_All_tot, as.numeric(r2mlm_prop_Best$R2s['f2', 'total']))
    Neg_Best_prop_Best_var$between_All_btw<-c(Neg_Best_prop_Best_var$between_All_btw, as.numeric(r2mlm_prop_Best$R2s['f2', 'between']))
    Neg_Best_prop_Best_var$between_res_btw<-c(Neg_Best_prop_Best_var$between_res_btw, as.numeric(r2mlm_prop_Best$R2s['m', 'between']))
    Neg_Best_prop_Best_var$within_fix_wthn<-c(Neg_Best_prop_Best_var$within_fix_wthn, as.numeric(r2mlm_prop_Best$Decompositions['fixed, within', 'within']))
    Neg_Best_prop_Best_var$within_fix_tot<-c(Neg_Best_prop_Best_var$within_fix_tot, as.numeric(r2mlm_prop_Best$Decompositions['fixed, within', 'total']))
    Neg_Best_prop_Best_var$within_slope_var_wthn<-c(Neg_Best_prop_Best_var$within_slope_var_wthn, as.numeric(r2mlm_prop_Best$Decompositions['slope variation', 'within']))
    Neg_Best_prop_Best_var$within_res_wthn<-c(Neg_Best_prop_Best_var$within_res_wthn, as.numeric(r2mlm_prop_Best$Decompositions['sigma2', 'within']))
    Neg_Best_prop_Best_var$within_unmod_tot<-c(Neg_Best_prop_Best_var$within_unmod_tot, as.numeric(r2mlm_prop_Best$Decompositions['sigma2', 'total'])+
                                                   as.numeric(r2mlm_prop_Best$Decompositions['slope variation', 'total']))
    
    print(paste("Sampling", d, "out of", length(D), '- Dataset:', i))
  }
  print(paste("*** Finshed imputed data set", i, 'out of', M, "***"))
}
hist(Neg_Best_prop_Best_var$between_var)
mean(Neg_Best_prop_Best_var$between_var)
hist(Neg_Best_prop_Best_var$within_var)

#MODEL 3 - LV2 Predictors = DN and Mean Negative Events, LV1 = individually mean-centered negative event indicator
#=-=-=-=-=-=-=-=-=-=-=-=-=-=
within_cov<-c(24)                    #Columns with group-mean centered predictors
between_cov<-c(8,25)                       #Columns with between-subject predictors
random_cov<-c(24)                    #Level 1 predictors with random slopes

#Getting distributions of fixed effects
Intercept = posterior_samples(Neg_Best_All, pars = 'b_Intercept') 
DN = posterior_samples(Neg_Best_All, pars = 'b_c.DN') 
prop_Best = posterior_samples(Neg_Best_All, pars = 'b_prop_Best') 
Best = posterior_samples(Neg_Best_All, pars = 'b_c.Best') 

#Getting between-subjects variance terms for random effects
Int_var = posterior_samples(Neg_Best_All, pars = 'sd_ID__Intercept')^2
Best_var = posterior_samples(Neg_Best_All, pars = 'sd_ID__c.Best')^2 

#Getting covariances of random effects (i.e., tau matrix)
cov_Int_Best = posterior_samples(Neg_Best_All, pars = 'sd_ID__Intercept')*
  posterior_samples(Neg_Best_All, pars = 'sd_ID__c.Best')*
  posterior_samples(Neg_Best_All, pars = 'cor_ID__Intercept__c.Best')

#Getting level 1 error variance
sigma<-log(1+posterior_samples(Neg_Best_All, pars = 'sigma')/beta_00)

post_samples<-data.frame(Intercept,
                         DN,
                         prop_Best, 
                         Best, 
                         Int_var, 
                         Best_var, 
                         cov_Int_Best, 
                         sigma)

colnames(post_samples)<-c('Intercept', 
                          'DN',
                          'prop_Best', 
                          'Best', 
                          'Int_var', 
                          'Best_var', 
                          'cov_Int_Best', 
                          'sigma')

# ***** IMPORTANT *****
# REMOVE ALL SECTIONS LIKE THE ONE BELOW
# ***** IMPORTANT *****
#Aggregate across imputed datasets
#Currently going to take 1000 draws from posterior distributions
#Will then apply across all 20 data sets

Neg_Best_All_var<-list(between_var=vector(), 
                        within_var=vector(), 
                        between_All_tot=vector(),
                        between_All_btw=vector(),
                        between_res_btw=vector(),
                        within_fix_wthn=vector(),
                        within_fix_tot=vector(),
                        within_slope_var_wthn=vector(),
                        within_res_wthn=vector(), 
                        within_unmod_tot=vector())

for(i in 1:length(dat.study1.list)){
  print(paste('Starting w/ imputed set', i, 'out of', M))
  samp_frame<-1:post_iter+(i-1)*post_iter
  D<-sample(samp_frame, size = 400, replace = FALSE)
  for(d in 1:length(D)){
    #browser()
    Gamma_w<-as.vector(
      as.matrix(
        post_samples[D[d], c('Best')]))       #Make sure the effects line up - in order of within_cov
    Gamma_b<-as.vector(
      as.matrix(
        post_samples[D[d], c('Intercept', 'DN', 'prop_Best')]))    #level fixed intercept and any level 2 fixed effects
    tau<-c(post_samples[D[d], c('Int_var', 'cov_Int_Best')], 
           post_samples[D[d], c('cov_Int_Best', 'Best_var')])
    
    tau<-matrix(unlist(tau), 
                byrow=TRUE, 
                ncol=2)
    
    #Matrix rows/columns ordered starting with tau.00 (intercept variance)
    #then add columns/rows for random effects of slopes in order of within_cov
    
    Sigma2<-post_samples[D[d], 'sigma']
    
    r2mlm_prop_Best<-r2MLM(data=dat.study1.list[[i]], 
                            within_covs = within_cov, 
                            between_covs = between_cov, 
                            random_covs = random_cov, 
                            gamma_w = Gamma_w, 
                            gamma_b = Gamma_b, 
                            Tau = tau, 
                            sigma2 = Sigma2, 
                            has_intercept = TRUE, 
                            clustermeancentered = TRUE)
    
    #Extracting relevant values across imputed datasets
    Neg_Best_All_var$between_var<-c(Neg_Best_All_var$between_var, as.numeric(r2mlm_prop_Best$Decompositions['fixed, between', 'total'])+
                                       as.numeric(r2mlm_prop_Best$Decompositions['mean variation', 'total']))
    Neg_Best_All_var$within_var<-c(Neg_Best_All_var$within_var, as.numeric(r2mlm_prop_Best$Decompositions['sigma2', 'total'])+
                                      as.numeric(r2mlm_prop_Best$Decompositions['slope variation', 'total'])+
                                      as.numeric(r2mlm_prop_Best$Decompositions['fixed, within', 'total']))
    Neg_Best_All_var$between_All_tot<-c(Neg_Best_All_var$between_All_tot, as.numeric(r2mlm_prop_Best$R2s['f2', 'total']))
    Neg_Best_All_var$between_All_btw<-c(Neg_Best_All_var$between_All_btw, as.numeric(r2mlm_prop_Best$R2s['f2', 'between']))
    Neg_Best_All_var$between_res_btw<-c(Neg_Best_All_var$between_res_btw, as.numeric(r2mlm_prop_Best$R2s['m', 'between']))
    Neg_Best_All_var$within_fix_wthn<-c(Neg_Best_All_var$within_fix_wthn, as.numeric(r2mlm_prop_Best$Decompositions['fixed, within', 'within']))
    Neg_Best_All_var$within_fix_tot<-c(Neg_Best_All_var$within_fix_tot, as.numeric(r2mlm_prop_Best$Decompositions['fixed, within', 'total']))
    Neg_Best_All_var$within_slope_var_wthn<-c(Neg_Best_All_var$within_slope_var_wthn, as.numeric(r2mlm_prop_Best$Decompositions['slope variation', 'within']))
    Neg_Best_All_var$within_res_wthn<-c(Neg_Best_All_var$within_res_wthn, as.numeric(r2mlm_prop_Best$Decompositions['sigma2', 'within']))
    Neg_Best_All_var$within_unmod_tot<-c(Neg_Best_All_var$within_unmod_tot, as.numeric(r2mlm_prop_Best$Decompositions['sigma2', 'total'])+
                                            as.numeric(r2mlm_prop_Best$Decompositions['slope variation', 'total']))
    
    print(paste("Sampling", d, "out of", length(D), '- Dataset:', i))
  }
  print(paste("*** Finshed imputed data set", i, 'out of', M, "***"))
}
hist(Neg_Best_All_var$between_var)
mean(Neg_Best_All_var$between_var)
hist(Neg_Best_All_var$within_var)

#MODEL 4 - Full cross-level interaction with 
#=-=-=-=-=-=-=-=-=-=-=-=-=-=
within_cov<-c(24,27)                    #Columns with group-mean centered predictors
between_cov<-c(8,25)                       #Columns with between-subject predictors
random_cov<-c(24)                    #Level 1 predictors with random slopes

#Getting distributions of fixed effects
Intercept = posterior_samples(Neg_Best_cross, pars = 'b_Intercept') 
DN = posterior_samples(Neg_Best_cross, pars = 'b_c.DN') 
prop_Best = posterior_samples(Neg_Best_cross, pars = 'b_prop_Best') 
Best = posterior_samples(Neg_Best_cross, pars = 'b_c.Best')[,1]  #Note posterior_samples() searches for text match
#With cross-level interaction returns two values 
#Need to select the one with the fixed effect for lv-1
DNxBest = posterior_samples(Neg_Best_cross, pars = 'b_c.Best:c.DN')

#Getting between-subjects variance terms for random effects
Int_var = posterior_samples(Neg_Best_cross, pars = 'sd_ID__Intercept')^2
Best_var = posterior_samples(Neg_Best_cross, pars = 'sd_ID__c.Best')^2 

#Getting covariances of random effects (i.e., tau matrix)
cov_Int_Best = posterior_samples(Neg_Best_cross, pars = 'sd_ID__Intercept')*
  posterior_samples(Neg_Best_cross, pars = 'sd_ID__c.Best')*
  posterior_samples(Neg_Best_cross, pars = 'cor_ID__Intercept__c.Best')

#Getting level 1 error variance
sigma<-log(1+posterior_samples(Neg_Best_cross, pars = 'sigma')/beta_00)

post_samples<-data.frame(Intercept,
                         DN,
                         prop_Best, 
                         Best,
                         DNxBest,
                         Int_var, 
                         Best_var, 
                         cov_Int_Best, 
                         sigma)

colnames(post_samples)<-c('Intercept', 
                          'DN',
                          'prop_Best', 
                          'Best', 
                          'DNxBest',
                          'Int_var', 
                          'Best_var', 
                          'cov_Int_Best', 
                          'sigma')

#Aggregate across imputed datasets
#Currently going to take 1000 draws from posterior distributions
#Will then apply across all 20 data sets

Neg_Best_cross_var<-list(between_var=vector(), 
                          within_var=vector(), 
                          between_All_tot=vector(),
                          between_All_btw=vector(),
                          between_res_btw=vector(),
                          within_fix_wthn=vector(),
                          within_fix_tot=vector(),
                          within_slope_var_wthn=vector(),
                          within_res_wthn=vector(), 
                          within_unmod_tot=vector())

for(i in 1:length(dat.study1.list)){
  print(paste('Starting w/ imputed set', i, 'out of', M))
  samp_frame<-1:post_iter+(i-1)*post_iter
  D<-sample(samp_frame, size = 400, replace = FALSE)
  for(d in 1:length(D)){
    #browser()
    Gamma_w<-as.vector(
      as.matrix(
        post_samples[D[d], c('Best', 'DNxBest')]))       #Make sure the effects line up - in order of within_cov
    Gamma_b<-as.vector(
      as.matrix(
        post_samples[D[d], c('Intercept', 'DN', 'prop_Best')]))    #level fixed intercept and any level 2 fixed effects
    tau<-c(post_samples[D[d], c('Int_var', 'cov_Int_Best')], 
           post_samples[D[d], c('cov_Int_Best', 'Best_var')])
    
    tau<-matrix(unlist(tau), 
                byrow=TRUE, 
                ncol=2)
    
    #Matrix rows/columns ordered starting with tau.00 (intercept variance)
    #then add columns/rows for random effects of slopes in order of within_cov
    
    Sigma2<-post_samples[D[d], 'sigma']
    
    r2mlm_prop_Best<-r2MLM(data=dat.study1.list[[i]], 
                            within_covs = within_cov, 
                            between_covs = between_cov, 
                            random_covs = random_cov, 
                            gamma_w = Gamma_w, 
                            gamma_b = Gamma_b, 
                            Tau = tau, 
                            sigma2 = Sigma2, 
                            has_intercept = TRUE, 
                            clustermeancentered = TRUE)
    
    #Extracting relevant values across imputed datasets
    Neg_Best_cross_var$between_var<-c(Neg_Best_cross_var$between_var, as.numeric(r2mlm_prop_Best$Decompositions['fixed, between', 'total'])+
                                         as.numeric(r2mlm_prop_Best$Decompositions['mean variation', 'total']))
    Neg_Best_cross_var$within_var<-c(Neg_Best_cross_var$within_var, as.numeric(r2mlm_prop_Best$Decompositions['sigma2', 'total'])+
                                        as.numeric(r2mlm_prop_Best$Decompositions['slope variation', 'total'])+
                                        as.numeric(r2mlm_prop_Best$Decompositions['fixed, within', 'total']))
    Neg_Best_cross_var$between_All_tot<-c(Neg_Best_cross_var$between_All_tot, as.numeric(r2mlm_prop_Best$R2s['f2', 'total']))
    Neg_Best_cross_var$between_All_btw<-c(Neg_Best_cross_var$between_All_btw, as.numeric(r2mlm_prop_Best$R2s['f2', 'between']))
    Neg_Best_cross_var$between_res_btw<-c(Neg_Best_cross_var$between_res_btw, as.numeric(r2mlm_prop_Best$R2s['m', 'between']))
    Neg_Best_cross_var$within_fix_wthn<-c(Neg_Best_cross_var$within_fix_wthn, as.numeric(r2mlm_prop_Best$Decompositions['fixed, within', 'within']))
    Neg_Best_cross_var$within_fix_tot<-c(Neg_Best_cross_var$within_fix_tot, as.numeric(r2mlm_prop_Best$Decompositions['fixed, within', 'total']))
    Neg_Best_cross_var$within_slope_var_wthn<-c(Neg_Best_cross_var$within_slope_var_wthn, as.numeric(r2mlm_prop_Best$Decompositions['slope variation', 'within']))
    Neg_Best_cross_var$within_res_wthn<-c(Neg_Best_cross_var$within_res_wthn, as.numeric(r2mlm_prop_Best$Decompositions['sigma2', 'within']))
    Neg_Best_cross_var$within_unmod_tot<-c(Neg_Best_cross_var$within_unmod_tot, as.numeric(r2mlm_prop_Best$Decompositions['sigma2', 'total'])+
                                              as.numeric(r2mlm_prop_Best$Decompositions['slope variation', 'total']))
    
    print(paste("Sampling", d, "out of", length(D), '- Dataset:', i))
  }
  print(paste("*** Finshed imputed data set", i, 'out of', M, "***"))
}
hist(Neg_Best_cross_var$between_var)
hist(Neg_Best_cross_var$within_var)

###########################################################################################################
#Preparing values for River plots :: Negative Mood - Best Event Models ::
###########################################################################################################
#Checking across models - between subjects variance
mean(Neg_Best_c.DN_var$between_var)
mean(Neg_Best_prop_Best_var$between_var)
mean(Neg_Best_All_var$between_var)
mean(Neg_Best_cross_var$between_var)

#Differemce between R2's due to sampling variability - is only 0.0008066132
#Difference is <.1% in total variance
mean(Neg_Best_cross_var$between_var)-mean(Neg_Best_c.DN_var$between_var)

#Between Variance Model Terms: 
Neg_Best_btw<-mean(Neg_Best_var$between_var)
Neg_Best_btw_All<-mean(Neg_Best_cross_var$between_All_tot)
Neg_Best_btw_umod<-Neg_Best_btw-Neg_Best_btw_All

Neg_Best_btw_DN_unique<-(mean(Neg_Best_All_var$between_All_btw)-mean(Neg_Best_prop_Best_var$between_All_btw))*Neg_Best_btw
Neg_Best_btw_Exp_unique<-(mean(Neg_Best_All_var$between_All_btw)-mean(Neg_Best_c.DN_var$between_All_btw))*Neg_Best_btw

Neg_Best_btw_Shared<-Neg_Best_btw_All-(Neg_Best_btw_DN_unique+Neg_Best_btw_Exp_unique)

#Within Variance Model Terms: 
Neg_Best_wthn<-mean(Neg_Best_var$within_var)
Neg_Best_wthn_Best<-mean(Neg_Best_All_var$within_fix_tot)
Neg_Best_wthn_Best_DN<-(mean(Neg_Best_All_var$within_res_wthn)-mean(Neg_Best_cross_var$within_res_wthn))*Neg_Best_wthn
Neg_Best_wthn_unmod<-Neg_Best_wthn-Neg_Best_wthn_Best-Neg_Best_wthn_Best_DN
#Note final value produced above is a measure of total variance attributable to cross-level interaction

#Creating Total DN Effect for separate river plot decomposition: 
Total_DN<-Neg_Best_btw_DN_unique+Neg_Best_btw_Shared+Neg_Best_wthn_Best_DN

###########################################################################################################
#Main Variance Decompisition River Plot :: Negative Mood - Best Event Models :: 
###########################################################################################################
Neg_Best_River_DF<-data.frame(N1 = c('DN',
                                      'DN <--> PDE \n Exposure',
                                      'PDE \n Exposure', 
                                      'Unmodeled \n Between', 
                                      'Momentary \n PDE', 
                                      'DN x PDE \n Reactivity',
                                      'Unmodeled \n Within', 
                                      'Total \n Between', 
                                      'Total \n Within'), 
                               N2 = c('Total \n Between', 
                                      'Total \n Between', 
                                      'Total \n Between', 
                                      'Total \n Between', 
                                      'Total \n Within', 
                                      'Total \n Within',
                                      'Total \n Within', 
                                      'Total \n Variance', 
                                      'Total \n Variance'), 
                               Value = c(Neg_Best_btw_DN_unique, 
                                         Neg_Best_btw_Shared,
                                         Neg_Best_btw_Exp_unique, 
                                         Neg_Best_btw_umod,
                                         Neg_Best_wthn_Best, 
                                         Neg_Best_wthn_Best_DN,
                                         Neg_Best_wthn_unmod,
                                         Neg_Best_btw, 
                                         Neg_Best_wthn), 
                               ID = 1:9)

Neg_Best_River_DF$N1<-paste(Neg_Best_River_DF$N1, 
                             '\n', 
                             paste0(round(Neg_Best_River_DF$Value*100, 
                                          digits = 2), '%'
                             )
)

Neg_Best_River_DF$N2<-c(rep(Neg_Best_River_DF$N1[8], 4),
                         rep(Neg_Best_River_DF$N1[9], 3), 
                         rep('Total Variance', 2)
)

nodes<-data.frame(ID = c(Neg_Best_River_DF$N1, 
                         'Total Variance'), 
                  x = c(1,1,1,1,1,1,1,2,2,3), 
                  y = c(0,2,4,6,8,10,12,3,9,6))

palette = c(paste0(brewer.pal(9, "Blues"), 90)[c(6,6,6,6)], 
            paste0(brewer.pal(9, "Reds"), 90)[c(6,6,6)], 
            paste0(brewer.pal(9, "Blues"), 90)[9], 
            paste0(brewer.pal(9, "Reds"), 90)[9], 
            paste0(brewer.pal(9, 'Purples'), 90)[9])

styles = lapply(nodes$y, function(n) {
  list(col = palette[n], lty = 0, textcol = "black")
})

#Not sure why color is not mapping correctly - trying to force it with this fuction... 
for(i in 1:length(palette)){
  styles[[i]]$col<-palette[i]
}

names(styles) = nodes$ID

Neg_Best_riv<-makeRiver(nodes = nodes, 
                         edges =  Neg_Best_River_DF, 
                         node_styles = styles)

png(paste0(study1.graphics, '/S1_Neg_Best_river_overall.png'), 
    units = 'in', 
    res = 1200, 
    height = 14, 
    width = 8)
riverplot(Neg_Best_riv, 
          nodewidth = 3, 
          plot_area = 1)
title(ylab = 'Study 1 - Riverplot of Total Variance Decomposition Estimated from Negative Mood Models: Positive Daily Events')
dev.off()

#----------------------------------------------------------------------------------------------------------
#DN - Specific Plot 
# :: Negative Mood - Best Event Models :: 
#----------------------------------------------------------------------------------------------------------
Neg_Best_DN_River_DF<-data.frame(N1 = c('DN',
                                         'DN <--> PDE Exposure',
                                         'DN x PDE Reactivity'
), 
N2 = rep(paste0("Combined DN Effect ",
                round(Total_DN*100, digits = 2),
                '%'), 
         3), 
Value = c(Neg_Best_btw_DN_unique, 
          Neg_Best_btw_Shared, 
          Neg_Best_wthn_Best_DN), 
ID = 1:3)

Neg_Best_DN_River_DF$N1<-paste(Neg_Best_DN_River_DF$N1, 
                                '\n', 
                                paste0(round(Neg_Best_DN_River_DF$Value*100, 
                                             digits = 2), '%'
                                )
)

nodes<-data.frame(ID = c(Neg_Best_DN_River_DF$N1, 
                         paste0("Combined DN Effect ",
                                round(Total_DN*100, digits = 2),
                                '%')), 
                  x = c(1,1,1,2), 
                  y = c(0,1,2,1))

palette = c(paste0(brewer.pal(9, "Blues"),90)[c(7,9)],
            paste0(brewer.pal(9, "Reds"),90)[7], 
            paste0(brewer.pal(9, "Purples"),90)[9])

styles = lapply(nodes$y, function(n) {
  list(col = palette[n], lty = 0, textcol = "black")
})

#Not sure why color is not mapping correctly - trying to force it with this fuction... 
for(i in 1:length(palette)){
  styles[[i]]$col<-palette[i]
}

names(styles) = nodes$ID

Neg_Best_DN_riv<-makeRiver(nodes = nodes, 
                            edges =  Neg_Best_DN_River_DF, 
                            node_styles = styles)

png(paste0(study1.graphics, '/S1_Neg_Best_river_DN.png'), 
    units = 'in', 
    res = 900, 
    height = 10, 
    width = 10)
riverplot(Neg_Best_DN_riv, 
          nodewidth = 3, 
          plot_area = .95)
title(ylab = 'Study 1 - Riverplot of DN Effect Decomposition Derived from Negative Mood Model: Positive Daily Events')
dev.off()

#
#
#
# Below 
# Are 
# The 
# Models
# Attuned 
# Not to negative but to positive mood
#
#
#




###########################################################################################################
#Positive Affect Models
###########################################################################################################
mean(dat.study1$PosAff, na.rm=TRUE)

Pos_ucm<-brms::brm_multiple(T.PosAff~1+(1|ID), 
                            data = dat.study1.list, 
                            family = 'normal',
                            prior = c(set_prior('normal(0,5)', class='sd'), 
                                      set_prior('normal(3.10,1)', class='Intercept')),
                            warmup = 5000, 
                            iter = 6000,
                            chains = 3,
                            control = list(adapt_delta=.99, 
                                           max_treedepth=15), 
                            save_model = paste0(stan.code, '/S1_ucm.stan'))
pp_check(Pos_ucm, nsamples = 100)
pairs(Pos_ucm)
round(Pos_ucm$rhats[,1:3], 3)
summary(Pos_ucm)

#----------------------------------------------------------------------------------------------------------
#Worst Event Models
#----------------------------------------------------------------------------------------------------------

#--
Pos_Worst<-brms::brm_multiple(T.PosAff~1+c.Worst+(1+c.Worst|ID), 
                              data = dat.study1.list, 
                              family = 'normal',
                              prior = c(set_prior('normal(0,5)', class='sd'), 
                                        set_prior('normal(3.10,1)', class='Intercept'),
                                        set_prior('lkj(2)', class = 'cor'), 
                                        set_prior('normal(0,5)', class='b')),
                              warmup = 5000, 
                              iter = 6000,
                              chains = 2,
                              control = list(adapt_delta=.95, 
                                             max_treedepth=15), 
                              save_model = paste0(stan.code, '/S1_Worst.stan'))

#--
Pos_Worst_c.DN<-brms::brm_multiple(T.PosAff~1+c.Worst+
                                     c.DN+
                                     (1+c.Worst|ID), 
                                   data = dat.study1.list, 
                                   family = 'normal',
                                   prior = c(set_prior('normal(0,5)', class='sd'), 
                                             set_prior('normal(3.10,1)', class='Intercept'),
                                             set_prior('lkj(2)', class = 'cor'), 
                                             set_prior('normal(0,5)', class='b')),
                                   warmup = 5000, 
                                   iter = 6000,
                                   chains = 2,
                                   control = list(adapt_delta=.95, 
                                                  max_treedepth=15), 
                                   save_model = paste0(stan.code, '/S1_Worst_DN.stan'))

#--
Pos_Worst_prop_Worst<-brms::brm_multiple(T.PosAff~1+c.Worst+
                                           prop_Worst+
                                           (1+c.Worst|ID), 
                                         data = dat.study1.list, 
                                         family = 'normal',
                                         prior = c(set_prior('normal(0,5)', class='sd'), 
                                                   set_prior('normal(3.10,1)', class='Intercept'),
                                                   set_prior('lkj(2)', class = 'cor'), 
                                                   set_prior('normal(0,5)', class='b')),
                                         warmup = 5000, 
                                         iter = 6000,
                                         chains = 2,
                                         control = list(adapt_delta=.95, 
                                                        max_treedepth=15), 
                                         save_model = paste0(stan.code, '/S1_Worst_prop_Worst.stan'))

#--
Pos_Worst_All<-brms::brm_multiple(T.PosAff~1+c.Worst+
                                    c.DN+prop_Worst+
                                    (1+c.Worst|ID), 
                                  data = dat.study1.list, 
                                  family = 'normal',
                                  prior = c(set_prior('normal(0,5)', class='sd'), 
                                            set_prior('normal(3.10,1)', class='Intercept'),
                                            set_prior('lkj(2)', class = 'cor'), 
                                            set_prior('normal(0,5)', class='b')),
                                  warmup = 5000, 
                                  iter = 6000,
                                  chains = 2,
                                  control = list(adapt_delta=.95, 
                                                 max_treedepth=15), 
                                  save_model = paste0(stan.code, '/S1_Worst_All.stan'))

#--
Pos_Worst_cross<-brms::brm_multiple(T.PosAff~1+c.Worst+
                                      c.DN+prop_Worst+
                                      c.Worst:c.DN+
                                      (1+c.Worst|ID), 
                                    data = dat.study1.list, 
                                    family = 'normal',
                                    prior = c(set_prior('normal(0,5)', class='sd'), 
                                              set_prior('normal(3.10,1)', class='Intercept'),
                                              set_prior('lkj(2)', class = 'cor'), 
                                              set_prior('normal(0,5)', class='b')),
                                    warmup = 5000, 
                                    iter = 6000,
                                    chains = 2,
                                    control = list(adapt_delta=.95, 
                                                   max_treedepth=15), 
                                    save_model = paste0(stan.code, '/S1_Worst_cross.stan'))

#----------------------------------------------------------------------------------------------------------
#Best Event Models
#----------------------------------------------------------------------------------------------------------

#--
Pos_Best<-brms::brm_multiple(T.PosAff~1+c.Best+(1+c.Best|ID), 
                             data = dat.study1.list, 
                             family = 'normal',
                             prior = c(set_prior('normal(0,5)', class='sd'), 
                                       set_prior('normal(3.10,1)', class='Intercept'),
                                       set_prior('lkj(2)', class = 'cor'), 
                                       set_prior('normal(0,5)', class='b')),
                             warmup = 5000, 
                             iter = 6000,
                             chains = 2,
                             control = list(adapt_delta=.95, 
                                            max_treedepth=15))

#--
Pos_Best_c.DN<-brms::brm_multiple(T.PosAff~1+c.Best+
                                    c.DN+
                                    (1+c.Best|ID), 
                                  data = dat.study1.list, 
                                  family = 'normal',
                                  prior = c(set_prior('normal(0,5)', class='sd'), 
                                            set_prior('normal(3.10,1)', class='Intercept'),
                                            set_prior('lkj(2)', class = 'cor'), 
                                            set_prior('normal(0,5)', class='b')),
                                  warmup = 5000, 
                                  iter = 6000,
                                  chains = 2,
                                  control = list(adapt_delta=.95, 
                                                 max_treedepth=15))

#--
Pos_Best_prop_Best<-brms::brm_multiple(T.PosAff~1+c.Best+
                                         prop_Best+
                                         (1+c.Best|ID), 
                                       data = dat.study1.list, 
                                       family = 'normal',
                                       prior = c(set_prior('normal(0,5)', class='sd'), 
                                                 set_prior('normal(3.10,1)', class='Intercept'),
                                                 set_prior('lkj(2)', class = 'cor'), 
                                                 set_prior('normal(0,5)', class='b')),
                                       warmup = 5000, 
                                       iter = 6000,
                                       chains = 2,
                                       control = list(adapt_delta=.95, 
                                                      max_treedepth=15))

#--
Pos_Best_All<-brms::brm_multiple(T.PosAff~1+c.Best+
                                   c.DN+prop_Best+
                                   (1+c.Best|ID), 
                                 data = dat.study1.list, 
                                 family = 'normal',
                                 prior = c(set_prior('normal(0,5)', class='sd'), 
                                           set_prior('normal(3.10,1)', class='Intercept'),
                                           set_prior('lkj(2)', class = 'cor'), 
                                           set_prior('normal(0,5)', class='b')),
                                 warmup = 5000, 
                                 iter = 6000,
                                 chains = 2,
                                 control = list(adapt_delta=.95, 
                                                max_treedepth=15))

#--
Pos_Best_cross<-brms::brm_multiple(T.PosAff~1+c.Best+
                                     c.DN+prop_Best+
                                     c.Best:c.DN+
                                     (1+c.Best|ID), 
                                   data = dat.study1.list, 
                                   family = 'normal',
                                   prior = c(set_prior('normal(0,5)', class='sd'), 
                                             set_prior('normal(3.10,1)', class='Intercept'),
                                             set_prior('lkj(2)', class = 'cor'), 
                                             set_prior('normal(0,5)', class='b')),
                                   warmup = 5000, 
                                   iter = 6000,
                                   chains = 2,
                                   control = list(adapt_delta=.95, 
                                                  max_treedepth=15))

#==========================================================================================================
#POSITIVE MOOD MODELS - Variance Decomposition
# :: Positive Mood - Worst Event Models ::
#==========================================================================================================

#==========================================================================================================
#----------------------------------------------------------------------------------------------------------
#Positive Mood, Positive Events, and Dispositional Posativity: 
#==========================================================================================================
post_iter<-2000   #No. of chains per data set x No. post-warmup iterations per chain
#for loops below require this value is correctly specified

#MODEL 1 - LV2 Predictor = DN, LV1 = individually mean-centered negative event indicator
#=-=-=-=-=-=-=-=-=-=-=-=-=-=
within_cov<-c(22)                    #Columns with group-mean centered predictors
random_cov<-c(22)                    #Level 1 predictors with random slopes

#Getting distributions of fixed effects
Intercept = posterior_samples(Pos_Worst, pars = 'b_Intercept') 
Worst = posterior_samples(Pos_Worst, pars = 'b_c.Worst') 

#Getting between-subjects variance terms for random effects
Int_var = posterior_samples(Pos_Worst, pars = 'sd_ID__Intercept')^2
Worst_var = posterior_samples(Pos_Worst, pars = 'sd_ID__c.Worst')^2 

#Getting covariances of random effects (i.e., tau matrix)
cov_Int_Worst = posterior_samples(Pos_Worst, pars = 'sd_ID__Intercept')*
  posterior_samples(Pos_Worst, pars = 'sd_ID__c.Worst')*
  posterior_samples(Pos_Worst, pars = 'cor_ID__Intercept__c.Worst')

#Getting level 1 error variance
sigma<-posterior_samples(Pos_Worst, pars = 'sigma')

post_samples<-data.frame(Intercept,
                         Worst, 
                         Int_var, 
                         Worst_var, 
                         cov_Int_Worst, 
                         sigma)

colnames(post_samples)<-c('Intercept', 
                          'Worst', 
                          'Int_var', 
                          'Worst_var', 
                          'cov_Int_Worst', 
                          'sigma')

#Aggregate across imputed datasets
#Currently going to take 1000 draws from posterior distributions
#Will then apply across all 20 data sets

Pos_Worst_var<-list(between_var=vector(), 
                         within_var=vector(), 
                         between_All_tot=vector(),
                         between_All_btw=vector(),
                         between_res_btw=vector(),
                         within_fix_wthn=vector(),
                         within_fix_tot=vector(),
                         within_slope_var_wthn=vector(),
                         within_res_wthn=vector(), 
                         within_unmod_tot=vector())

for(i in 1:length(dat.study1.list)){
  print(paste('Starting w/ imputed set', i, 'out of', M))
  samp_frame<-1:post_iter+(i-1)*post_iter
  D<-sample(samp_frame, size = 400, replace = FALSE)
  for(d in 1:length(D)){
    #browser()
    Gamma_w<-as.vector(
      as.matrix(
        post_samples[D[d], c('Worst')]))       #Make sure the effects line up - in order of within_cov
    Gamma_b<-as.vector(
      as.matrix(
        post_samples[D[d], c('Intercept')]))    #level fixed intercept and any level 2 fixed effects
    tau<-c(post_samples[D[d], c('Int_var', 'cov_Int_Worst')], 
           post_samples[D[d], c('cov_Int_Worst', 'Worst_var')])
    
    tau<-matrix(unlist(tau), 
                byrow=TRUE, 
                ncol=2)
    
    #Matrix rows/columns ordered starting with tau.00 (intercept variance)
    #then add columns/rows for random effects of slopes in order of within_cov
    
    Sigma2<-post_samples[D[d], 'sigma']
    
    r2mlm_DN<-r2MLM(data=dat.study1.list[[i]], 
                    within_covs = within_cov, 
                    between_covs = NULL, 
                    random_covs = random_cov, 
                    gamma_w = Gamma_w, 
                    gamma_b = Gamma_b, 
                    Tau = tau, 
                    sigma2 = Sigma2, 
                    has_intercept = TRUE, 
                    clustermeancentered = TRUE)
    
    #Extracting relevant values across imputed datasets
    Pos_Worst_var$between_var<-c(Pos_Worst_var$between_var, as.numeric(r2mlm_DN$Decompositions['fixed, between', 'total'])+
                                        as.numeric(r2mlm_DN$Decompositions['mean variation', 'total']))
    Pos_Worst_var$within_var<-c(Pos_Worst_var$within_var, as.numeric(r2mlm_DN$Decompositions['sigma2', 'total'])+
                                       as.numeric(r2mlm_DN$Decompositions['slope variation', 'total'])+
                                       as.numeric(r2mlm_DN$Decompositions['fixed, within', 'total']))
    Pos_Worst_var$between_All_tot<-c(Pos_Worst_var$between_All_tot, as.numeric(r2mlm_DN$R2s['f2', 'total']))
    Pos_Worst_var$between_All_btw<-c(Pos_Worst_var$between_All_btw, as.numeric(r2mlm_DN$R2s['f2', 'between']))
    Pos_Worst_var$between_res_btw<-c(Pos_Worst_var$between_res_btw, as.numeric(r2mlm_DN$R2s['m', 'between']))
    Pos_Worst_var$within_fix_wthn<-c(Pos_Worst_var$within_fix_wthn, as.numeric(r2mlm_DN$Decompositions['fixed, within', 'within']))
    Pos_Worst_var$within_fix_tot<-c(Pos_Worst_var$within_fix_tot, as.numeric(r2mlm_DN$Decompositions['fixed, within', 'total']))
    Pos_Worst_var$within_slope_var_wthn<-c(Pos_Worst_var$within_slope_var_wthn, as.numeric(r2mlm_DN$Decompositions['slope variation', 'within']))
    Pos_Worst_var$within_res_wthn<-c(Pos_Worst_var$within_res_wthn, as.numeric(r2mlm_DN$Decompositions['sigma2', 'within']))
    Pos_Worst_var$within_unmod_tot<-c(Pos_Worst_var$within_unmod_tot, as.numeric(r2mlm_DN$Decompositions['sigma2', 'total'])+
                                             as.numeric(r2mlm_DN$Decompositions['slope variation', 'total']))
    
    print(paste("Sampling", d, "out of", length(D), '- Dataset:', i))
  }
  print(paste("*** Finshed imputed data set", i, 'out of', M, "***"))
}
hist(Pos_Worst_var$between_var)
hist(Pos_Worst_var$within_var)

gc()

#MODEL 1 - LV2 Predictor = DN, LV1 = individually mean-centered negative event indicator
#=-=-=-=-=-=-=-=-=-=-=-=-=-=
within_cov<-c(22)                    #Columns with group-mean centered predictors
between_cov<-c(8)                       #Columns with between-subject predictors
random_cov<-c(22)                    #Level 1 predictors with random slopes

#Getting distributions of fixed effects
Intercept = posterior_samples(Pos_Worst_c.DN, pars = 'b_Intercept') 
DN = posterior_samples(Pos_Worst_c.DN, pars = 'b_c.DN') 
Worst = posterior_samples(Pos_Worst_c.DN, pars = 'b_c.Worst') 

#Getting between-subjects variance terms for random effects
Int_var = posterior_samples(Pos_Worst_c.DN, pars = 'sd_ID__Intercept')^2
Worst_var = posterior_samples(Pos_Worst_c.DN, pars = 'sd_ID__c.Worst')^2 

#Getting covariances of random effects (i.e., tau matrix)
cov_Int_Worst = posterior_samples(Pos_Worst_c.DN, pars = 'sd_ID__Intercept')*
  posterior_samples(Pos_Worst_c.DN, pars = 'sd_ID__c.Worst')*
  posterior_samples(Pos_Worst_c.DN, pars = 'cor_ID__Intercept__c.Worst')

#Getting level 1 error variance
sigma<-posterior_samples(Pos_Worst_c.DN, pars = 'sigma')

post_samples<-data.frame(Intercept,
                         DN, 
                         Worst, 
                         Int_var, 
                         Worst_var, 
                         cov_Int_Worst, 
                         sigma)

colnames(post_samples)<-c('Intercept', 
                          'DN', 
                          'Worst', 
                          'Int_var', 
                          'Worst_var', 
                          'cov_Int_Worst', 
                          'sigma')

#Aggregate across imputed datasets
#Currently going to take 1000 draws from posterior distributions
#Will then apply across all 20 data sets

Pos_Worst_c.DN_var<-list(between_var=vector(), 
                         within_var=vector(), 
                         between_All_tot=vector(),
                         between_All_btw=vector(),
                         between_res_btw=vector(),
                         within_fix_wthn=vector(),
                         within_fix_tot=vector(),
                         within_slope_var_wthn=vector(),
                         within_res_wthn=vector(), 
                         within_unmod_tot=vector())

for(i in 1:length(dat.study1.list)){
  print(paste('Starting w/ imputed set', i, 'out of', M))
  samp_frame<-1:post_iter+(i-1)*post_iter
  D<-sample(samp_frame, size = 400, replace = FALSE)
  for(d in 1:length(D)){
    #browser()
    Gamma_w<-as.vector(
      as.matrix(
        post_samples[D[d], c('Worst')]))       #Make sure the effects line up - in order of within_cov
    Gamma_b<-as.vector(
      as.matrix(
        post_samples[D[d], c('Intercept', 'DN')]))    #level fixed intercept and any level 2 fixed effects
    tau<-c(post_samples[D[d], c('Int_var', 'cov_Int_Worst')], 
           post_samples[D[d], c('cov_Int_Worst', 'Worst_var')])
    
    tau<-matrix(unlist(tau), 
                byrow=TRUE, 
                ncol=2)
    
    #Matrix rows/columns ordered starting with tau.00 (intercept variance)
    #then add columns/rows for random effects of slopes in order of within_cov
    
    Sigma2<-post_samples[D[d], 'sigma']
    
    r2mlm_DN<-r2MLM(data=dat.study1.list[[i]], 
                    within_covs = within_cov, 
                    between_covs = between_cov, 
                    random_covs = random_cov, 
                    gamma_w = Gamma_w, 
                    gamma_b = Gamma_b, 
                    Tau = tau, 
                    sigma2 = Sigma2, 
                    has_intercept = TRUE, 
                    clustermeancentered = TRUE)
    
    #Extracting relevant values across imputed datasets
    Pos_Worst_c.DN_var$between_var<-c(Pos_Worst_c.DN_var$between_var, as.numeric(r2mlm_DN$Decompositions['fixed, between', 'total'])+
                                        as.numeric(r2mlm_DN$Decompositions['mean variation', 'total']))
    Pos_Worst_c.DN_var$within_var<-c(Pos_Worst_c.DN_var$within_var, as.numeric(r2mlm_DN$Decompositions['sigma2', 'total'])+
                                       as.numeric(r2mlm_DN$Decompositions['slope variation', 'total'])+
                                       as.numeric(r2mlm_DN$Decompositions['fixed, within', 'total']))
    Pos_Worst_c.DN_var$between_All_tot<-c(Pos_Worst_c.DN_var$between_All_tot, as.numeric(r2mlm_DN$R2s['f2', 'total']))
    Pos_Worst_c.DN_var$between_All_btw<-c(Pos_Worst_c.DN_var$between_All_btw, as.numeric(r2mlm_DN$R2s['f2', 'between']))
    Pos_Worst_c.DN_var$between_res_btw<-c(Pos_Worst_c.DN_var$between_res_btw, as.numeric(r2mlm_DN$R2s['m', 'between']))
    Pos_Worst_c.DN_var$within_fix_wthn<-c(Pos_Worst_c.DN_var$within_fix_wthn, as.numeric(r2mlm_DN$Decompositions['fixed, within', 'within']))
    Pos_Worst_c.DN_var$within_fix_tot<-c(Pos_Worst_c.DN_var$within_fix_tot, as.numeric(r2mlm_DN$Decompositions['fixed, within', 'total']))
    Pos_Worst_c.DN_var$within_slope_var_wthn<-c(Pos_Worst_c.DN_var$within_slope_var_wthn, as.numeric(r2mlm_DN$Decompositions['slope variation', 'within']))
    Pos_Worst_c.DN_var$within_res_wthn<-c(Pos_Worst_c.DN_var$within_res_wthn, as.numeric(r2mlm_DN$Decompositions['sigma2', 'within']))
    Pos_Worst_c.DN_var$within_unmod_tot<-c(Pos_Worst_c.DN_var$within_unmod_tot, as.numeric(r2mlm_DN$Decompositions['sigma2', 'total'])+
                                             as.numeric(r2mlm_DN$Decompositions['slope variation', 'total']))
    
    print(paste("Sampling", d, "out of", length(D), '- Dataset:', i))
  }
  print(paste("*** Finshed imputed data set", i, 'out of', M, "***"))
}
hist(Pos_Worst_c.DN_var$between_var)
hist(Pos_Worst_c.DN_var$within_var)

gc()

#MODEL 2 - LV2 Predictor = Mean Positive Events, LV1 = individually mean-centered negative event indicator
#=-=-=-=-=-=-=-=-=-=-=-=-=-=
within_cov<-c(22)                    #Columns with group-mean centered predictors
between_cov<-c(23)                       #Columns with between-subject predictors
random_cov<-c(22)                    #Level 1 predictors with random slopes

#Getting distributions of fixed effects
Intercept = posterior_samples(Pos_Worst_prop_Worst, pars = 'b_Intercept') 
prop_Worst = posterior_samples(Pos_Worst_prop_Worst, pars = 'b_prop_Worst') 
Worst = posterior_samples(Pos_Worst_prop_Worst, pars = 'b_c.Worst') 

#Getting between-subjects variance terms for random effects
Int_var = posterior_samples(Pos_Worst_prop_Worst, pars = 'sd_ID__Intercept')^2
Worst_var = posterior_samples(Pos_Worst_prop_Worst, pars = 'sd_ID__c.Worst')^2 

#Getting covariances of random effects (i.e., tau matrix)
cov_Int_Worst = posterior_samples(Pos_Worst_prop_Worst, pars = 'sd_ID__Intercept')*
  posterior_samples(Pos_Worst_prop_Worst, pars = 'sd_ID__c.Worst')*
  posterior_samples(Pos_Worst_prop_Worst, pars = 'cor_ID__Intercept__c.Worst')

#Getting level 1 error variance
sigma<-posterior_samples(Pos_Worst_prop_Worst, pars = 'sigma')

post_samples<-data.frame(Intercept,
                         prop_Worst, 
                         Worst, 
                         Int_var, 
                         Worst_var, 
                         cov_Int_Worst, 
                         sigma)

colnames(post_samples)<-c('Intercept', 
                          'prop_Worst', 
                          'Worst', 
                          'Int_var', 
                          'Worst_var', 
                          'cov_Int_Worst', 
                          'sigma')

#Aggregate across imputed datasets
#Currently going to take 1000 draws from posterior distributions
#Will then apply across all 20 data sets

Pos_Worst_prop_Worst_var<-list(between_var=vector(), 
                               within_var=vector(), 
                               between_All_tot=vector(),
                               between_All_btw=vector(),
                               between_res_btw=vector(),
                               within_fix_wthn=vector(),
                               within_fix_tot=vector(),
                               within_slope_var_wthn=vector(),
                               within_res_wthn=vector(), 
                               within_unmod_tot=vector())

for(i in 1:length(dat.study1.list)){
  print(paste('Starting w/ imputed set', i, 'out of', M))
  samp_frame<-1:post_iter+(i-1)*post_iter
  D<-sample(samp_frame, size = 400, replace = FALSE)
  for(d in 1:length(D)){
    #browser()
    Gamma_w<-as.vector(
      as.matrix(
        post_samples[D[d], c('Worst')]))       #Make sure the effects line up - in order of within_cov
    Gamma_b<-as.vector(
      as.matrix(
        post_samples[D[d], c('Intercept', 'prop_Worst')]))    #level fixed intercept and any level 2 fixed effects
    tau<-c(post_samples[D[d], c('Int_var', 'cov_Int_Worst')], 
           post_samples[D[d], c('cov_Int_Worst', 'Worst_var')])
    
    tau<-matrix(unlist(tau), 
                byrow=TRUE, 
                ncol=2)
    
    #Matrix rows/columns ordered starting with tau.00 (intercept variance)
    #then add columns/rows for random effects of slopes in order of within_cov
    
    Sigma2<-post_samples[D[d], 'sigma']
    
    r2mlm_prop_Worst<-r2MLM(data=dat.study1.list[[i]], 
                            within_covs = within_cov, 
                            between_covs = between_cov, 
                            random_covs = random_cov, 
                            gamma_w = Gamma_w, 
                            gamma_b = Gamma_b, 
                            Tau = tau, 
                            sigma2 = Sigma2, 
                            has_intercept = TRUE, 
                            clustermeancentered = TRUE)
    
    #Extracting relevant values across imputed datasets
    Pos_Worst_prop_Worst_var$between_var<-c(Pos_Worst_prop_Worst_var$between_var, as.numeric(r2mlm_prop_Worst$Decompositions['fixed, between', 'total'])+
                                              as.numeric(r2mlm_prop_Worst$Decompositions['mean variation', 'total']))
    Pos_Worst_prop_Worst_var$within_var<-c(Pos_Worst_prop_Worst_var$within_var, as.numeric(r2mlm_prop_Worst$Decompositions['sigma2', 'total'])+
                                             as.numeric(r2mlm_prop_Worst$Decompositions['slope variation', 'total'])+
                                             as.numeric(r2mlm_prop_Worst$Decompositions['fixed, within', 'total']))
    Pos_Worst_prop_Worst_var$between_All_tot<-c(Pos_Worst_prop_Worst_var$between_All_tot, as.numeric(r2mlm_prop_Worst$R2s['f2', 'total']))
    Pos_Worst_prop_Worst_var$between_All_btw<-c(Pos_Worst_prop_Worst_var$between_All_btw, as.numeric(r2mlm_prop_Worst$R2s['f2', 'between']))
    Pos_Worst_prop_Worst_var$between_res_btw<-c(Pos_Worst_prop_Worst_var$between_res_btw, as.numeric(r2mlm_prop_Worst$R2s['m', 'between']))
    Pos_Worst_prop_Worst_var$within_fix_wthn<-c(Pos_Worst_prop_Worst_var$within_fix_wthn, as.numeric(r2mlm_prop_Worst$Decompositions['fixed, within', 'within']))
    Pos_Worst_prop_Worst_var$within_fix_tot<-c(Pos_Worst_prop_Worst_var$within_fix_tot, as.numeric(r2mlm_prop_Worst$Decompositions['fixed, within', 'total']))
    Pos_Worst_prop_Worst_var$within_slope_var_wthn<-c(Pos_Worst_prop_Worst_var$within_slope_var_wthn, as.numeric(r2mlm_prop_Worst$Decompositions['slope variation', 'within']))
    Pos_Worst_prop_Worst_var$within_res_wthn<-c(Pos_Worst_prop_Worst_var$within_res_wthn, as.numeric(r2mlm_prop_Worst$Decompositions['sigma2', 'within']))
    Pos_Worst_prop_Worst_var$within_unmod_tot<-c(Pos_Worst_prop_Worst_var$within_unmod_tot, as.numeric(r2mlm_prop_Worst$Decompositions['sigma2', 'total'])+
                                                   as.numeric(r2mlm_prop_Worst$Decompositions['slope variation', 'total']))
    
    print(paste("Sampling", d, "out of", length(D), '- Dataset:', i))
  }
  print(paste("*** Finshed imputed data set", i, 'out of', M, "***"))
}

hist(Pos_Worst_prop_Worst_var$between_var)
hist(Pos_Worst_prop_Worst_var$within_var)

#MODEL 3 - LV2 Predictors = DN and Mean Positive Events, LV1 = individually mean-centered negative event indicator
#=-=-=-=-=-=-=-=-=-=-=-=-=-=
within_cov<-c(22)                    #Columns with group-mean centered predictors
between_cov<-c(8,23)                       #Columns with between-subject predictors
random_cov<-c(22)                    #Level 1 predictors with random slopes

#Getting distributions of fixed effects
Intercept = posterior_samples(Pos_Worst_All, pars = 'b_Intercept') 
DN = posterior_samples(Pos_Worst_All, pars = 'b_c.DN') 
prop_Worst = posterior_samples(Pos_Worst_All, pars = 'b_prop_Worst') 
Worst = posterior_samples(Pos_Worst_All, pars = 'b_c.Worst') 

#Getting between-subjects variance terms for random effects
Int_var = posterior_samples(Pos_Worst_All, pars = 'sd_ID__Intercept')^2
Worst_var = posterior_samples(Pos_Worst_All, pars = 'sd_ID__c.Worst')^2 

#Getting covariances of random effects (i.e., tau matrix)
cov_Int_Worst = posterior_samples(Pos_Worst_All, pars = 'sd_ID__Intercept')*
  posterior_samples(Pos_Worst_All, pars = 'sd_ID__c.Worst')*
  posterior_samples(Pos_Worst_All, pars = 'cor_ID__Intercept__c.Worst')

#Getting level 1 error variance
sigma<-posterior_samples(Pos_Worst_All, pars = 'sigma')

post_samples<-data.frame(Intercept,
                         DN,
                         prop_Worst, 
                         Worst, 
                         Int_var, 
                         Worst_var, 
                         cov_Int_Worst, 
                         sigma)

colnames(post_samples)<-c('Intercept', 
                          'DN',
                          'prop_Worst', 
                          'Worst', 
                          'Int_var', 
                          'Worst_var', 
                          'cov_Int_Worst', 
                          'sigma')

#Aggregate across imputed datasets
#Currently going to take 1000 draws from posterior distributions
#Will then apply across all 20 data sets

Pos_Worst_All_var<-list(between_var=vector(), 
                        within_var=vector(), 
                        between_All_tot=vector(),
                        between_All_btw=vector(),
                        between_res_btw=vector(),
                        within_fix_wthn=vector(),
                        within_fix_tot=vector(),
                        within_slope_var_wthn=vector(),
                        within_res_wthn=vector(), 
                        within_unmod_tot=vector())

for(i in 1:length(dat.study1.list)){
  print(paste('Starting w/ imputed set', i, 'out of', M))
  samp_frame<-1:post_iter+(i-1)*post_iter
  D<-sample(samp_frame, size = 400, replace = FALSE)
  for(d in 1:length(D)){
    #browser()
    Gamma_w<-as.vector(
      as.matrix(
        post_samples[D[d], c('Worst')]))       #Make sure the effects line up - in order of within_cov
    Gamma_b<-as.vector(
      as.matrix(
        post_samples[D[d], c('Intercept', 'DN', 'prop_Worst')]))    #level fixed intercept and any level 2 fixed effects
    tau<-c(post_samples[D[d], c('Int_var', 'cov_Int_Worst')], 
           post_samples[D[d], c('cov_Int_Worst', 'Worst_var')])
    
    tau<-matrix(unlist(tau), 
                byrow=TRUE, 
                ncol=2)
    
    #Matrix rows/columns ordered starting with tau.00 (intercept variance)
    #then add columns/rows for random effects of slopes in order of within_cov
    
    Sigma2<-post_samples[D[d], 'sigma']
    
    r2mlm_prop_Worst<-r2MLM(data=dat.study1.list[[i]], 
                            within_covs = within_cov, 
                            between_covs = between_cov, 
                            random_covs = random_cov, 
                            gamma_w = Gamma_w, 
                            gamma_b = Gamma_b, 
                            Tau = tau, 
                            sigma2 = Sigma2, 
                            has_intercept = TRUE, 
                            clustermeancentered = TRUE)
    
    #Extracting relevant values across imputed datasets
    Pos_Worst_All_var$between_var<-c(Pos_Worst_All_var$between_var, as.numeric(r2mlm_prop_Worst$Decompositions['fixed, between', 'total'])+
                                       as.numeric(r2mlm_prop_Worst$Decompositions['mean variation', 'total']))
    Pos_Worst_All_var$within_var<-c(Pos_Worst_All_var$within_var, as.numeric(r2mlm_prop_Worst$Decompositions['sigma2', 'total'])+
                                      as.numeric(r2mlm_prop_Worst$Decompositions['slope variation', 'total'])+
                                      as.numeric(r2mlm_prop_Worst$Decompositions['fixed, within', 'total']))
    Pos_Worst_All_var$between_All_tot<-c(Pos_Worst_All_var$between_All_tot, as.numeric(r2mlm_prop_Worst$R2s['f2', 'total']))
    Pos_Worst_All_var$between_All_btw<-c(Pos_Worst_All_var$between_All_btw, as.numeric(r2mlm_prop_Worst$R2s['f2', 'between']))
    Pos_Worst_All_var$between_res_btw<-c(Pos_Worst_All_var$between_res_btw, as.numeric(r2mlm_prop_Worst$R2s['m', 'between']))
    Pos_Worst_All_var$within_fix_wthn<-c(Pos_Worst_All_var$within_fix_wthn, as.numeric(r2mlm_prop_Worst$Decompositions['fixed, within', 'within']))
    Pos_Worst_All_var$within_fix_tot<-c(Pos_Worst_All_var$within_fix_tot, as.numeric(r2mlm_prop_Worst$Decompositions['fixed, within', 'total']))
    Pos_Worst_All_var$within_slope_var_wthn<-c(Pos_Worst_All_var$within_slope_var_wthn, as.numeric(r2mlm_prop_Worst$Decompositions['slope variation', 'within']))
    Pos_Worst_All_var$within_res_wthn<-c(Pos_Worst_All_var$within_res_wthn, as.numeric(r2mlm_prop_Worst$Decompositions['sigma2', 'within']))
    Pos_Worst_All_var$within_unmod_tot<-c(Pos_Worst_All_var$within_unmod_tot, as.numeric(r2mlm_prop_Worst$Decompositions['sigma2', 'total'])+
                                            as.numeric(r2mlm_prop_Worst$Decompositions['slope variation', 'total']))
    
    print(paste("Sampling", d, "out of", length(D), '- Dataset:', i))
  }
  print(paste("*** Finshed imputed data set", i, 'out of', M, "***"))
}
hist(Pos_Worst_All_var$between_var)
hist(Pos_Worst_All_var$within_var)


#MODEL 4 - Full cross-level interaction with 
#=-=-=-=-=-=-=-=-=-=-=-=-=-=
within_cov<-c(22,26)                    #Columns with group-mean centered predictors
between_cov<-c(8,23)                       #Columns with between-subject predictors
random_cov<-c(22)                    #Level 1 predictors with random slopes

#Getting distributions of fixed effects
Intercept = posterior_samples(Pos_Worst_cross, pars = 'b_Intercept') 
DN = posterior_samples(Pos_Worst_cross, pars = 'b_c.DN') 
prop_Worst = posterior_samples(Pos_Worst_cross, pars = 'b_prop_Worst') 
Worst = posterior_samples(Pos_Worst_cross, pars = 'b_c.Worst')[,1]  #Note posterior_samples() searches for text match
#With cross-level interaction returns two values 
#Need to select the one with the fixed effect for lv-1
DNxWorst = posterior_samples(Pos_Worst_cross, pars = 'b_c.Worst:c.DN')

#Getting between-subjects variance terms for random effects
Int_var = posterior_samples(Pos_Worst_cross, pars = 'sd_ID__Intercept')^2
Worst_var = posterior_samples(Pos_Worst_cross, pars = 'sd_ID__c.Worst')^2 

#Getting covariances of random effects (i.e., tau matrix)
cov_Int_Worst = posterior_samples(Pos_Worst_cross, pars = 'sd_ID__Intercept')*
  posterior_samples(Pos_Worst_cross, pars = 'sd_ID__c.Worst')*
  posterior_samples(Pos_Worst_cross, pars = 'cor_ID__Intercept__c.Worst')

#Getting level 1 error variance
sigma<-posterior_samples(Pos_Worst_cross, pars = 'sigma')

post_samples<-data.frame(Intercept,
                         DN,
                         prop_Worst, 
                         Worst,
                         DNxWorst,
                         Int_var, 
                         Worst_var, 
                         cov_Int_Worst, 
                         sigma)

colnames(post_samples)<-c('Intercept', 
                          'DN',
                          'prop_Worst', 
                          'Worst', 
                          'DNxWorst',
                          'Int_var', 
                          'Worst_var', 
                          'cov_Int_Worst', 
                          'sigma')

#Aggregate across imputed datasets
#Currently going to take 1000 draws from posterior distributions
#Will then apply across all 20 data sets

Pos_Worst_cross_var<-list(between_var=vector(), 
                          within_var=vector(), 
                          between_All_tot=vector(),
                          between_All_btw=vector(),
                          between_res_btw=vector(),
                          within_fix_wthn=vector(),
                          within_fix_tot=vector(),
                          within_slope_var_wthn=vector(),
                          within_res_wthn=vector(), 
                          within_unmod_tot=vector())

for(i in 1:length(dat.study1.list)){
  print(paste('Starting w/ imputed set', i, 'out of', M))
  samp_frame<-1:post_iter+(i-1)*post_iter
  D<-sample(samp_frame, size = 400, replace = FALSE)
  for(d in 1:length(D)){
    #browser()
    Gamma_w<-as.vector(
      as.matrix(
        post_samples[D[d], c('Worst', 'DNxWorst')]))       #Make sure the effects line up - in order of within_cov
    Gamma_b<-as.vector(
      as.matrix(
        post_samples[D[d], c('Intercept', 'DN', 'prop_Worst')]))    #level fixed intercept and any level 2 fixed effects
    tau<-c(post_samples[D[d], c('Int_var', 'cov_Int_Worst')], 
           post_samples[D[d], c('cov_Int_Worst', 'Worst_var')])
    
    tau<-matrix(unlist(tau), 
                byrow=TRUE, 
                ncol=2)
    
    #Matrix rows/columns ordered starting with tau.00 (intercept variance)
    #then add columns/rows for random effects of slopes in order of within_cov
    
    Sigma2<-post_samples[D[d], 'sigma']
    
    r2mlm_prop_Worst<-r2MLM(data=dat.study1.list[[i]], 
                            within_covs = within_cov, 
                            between_covs = between_cov, 
                            random_covs = random_cov, 
                            gamma_w = Gamma_w, 
                            gamma_b = Gamma_b, 
                            Tau = tau, 
                            sigma2 = Sigma2, 
                            has_intercept = TRUE, 
                            clustermeancentered = TRUE)
    
    #Extracting relevant values across imputed datasets
    Pos_Worst_cross_var$between_var<-c(Pos_Worst_cross_var$between_var, as.numeric(r2mlm_prop_Worst$Decompositions['fixed, between', 'total'])+
                                         as.numeric(r2mlm_prop_Worst$Decompositions['mean variation', 'total']))
    Pos_Worst_cross_var$within_var<-c(Pos_Worst_cross_var$within_var, as.numeric(r2mlm_prop_Worst$Decompositions['sigma2', 'total'])+
                                        as.numeric(r2mlm_prop_Worst$Decompositions['slope variation', 'total'])+
                                        as.numeric(r2mlm_prop_Worst$Decompositions['fixed, within', 'total']))
    Pos_Worst_cross_var$between_All_tot<-c(Pos_Worst_cross_var$between_All_tot, as.numeric(r2mlm_prop_Worst$R2s['f2', 'total']))
    Pos_Worst_cross_var$between_All_btw<-c(Pos_Worst_cross_var$between_All_btw, as.numeric(r2mlm_prop_Worst$R2s['f2', 'between']))
    Pos_Worst_cross_var$between_res_btw<-c(Pos_Worst_cross_var$between_res_btw, as.numeric(r2mlm_prop_Worst$R2s['m', 'between']))
    Pos_Worst_cross_var$within_fix_wthn<-c(Pos_Worst_cross_var$within_fix_wthn, as.numeric(r2mlm_prop_Worst$Decompositions['fixed, within', 'within']))
    Pos_Worst_cross_var$within_fix_tot<-c(Pos_Worst_cross_var$within_fix_tot, as.numeric(r2mlm_prop_Worst$Decompositions['fixed, within', 'total']))
    Pos_Worst_cross_var$within_slope_var_wthn<-c(Pos_Worst_cross_var$within_slope_var_wthn, as.numeric(r2mlm_prop_Worst$Decompositions['slope variation', 'within']))
    Pos_Worst_cross_var$within_res_wthn<-c(Pos_Worst_cross_var$within_res_wthn, as.numeric(r2mlm_prop_Worst$Decompositions['sigma2', 'within']))
    Pos_Worst_cross_var$within_unmod_tot<-c(Pos_Worst_cross_var$within_unmod_tot, as.numeric(r2mlm_prop_Worst$Decompositions['sigma2', 'total'])+
                                              as.numeric(r2mlm_prop_Worst$Decompositions['slope variation', 'total']))
    
    print(paste("Sampling", d, "out of", length(D), '- Dataset:', i))
  }
  print(paste("*** Finshed imputed data set", i, 'out of', M, "***"))
}
hist(Pos_Worst_cross_var$between_var)
hist(Pos_Worst_cross_var$within_var)

###########################################################################################################
#Preparing values for River plots :: Positive Mood - Worst Event Models ::
###########################################################################################################
#Checking across models - between subjects variance
mean(Pos_Worst_c.DN_var$between_var)
mean(Pos_Worst_prop_Worst_var$between_var)
mean(Pos_Worst_All_var$between_var)
mean(Pos_Worst_cross_var$between_var)

#Differemce between R2's due to Bayesian sampling variability - is only 0.001621475
#Difference is <.2% in total variance
mean(Pos_Worst_cross_var$between_var)-mean(Pos_Worst_c.DN_var$between_var)

#Between Variance Model Terms: 
Pos_Worst_btw<-mean(Pos_Worst_cross_var$between_var)
Pos_Worst_btw_All<-mean(Pos_Worst_cross_var$between_All_tot)
Pos_Worst_btw_umod<-Pos_Worst_btw-Pos_Worst_btw_All

Pos_Worst_btw_DN_unique<-(mean(Pos_Worst_All_var$between_All_btw)-mean(Pos_Worst_prop_Worst_var$between_All_btw))*Pos_Worst_btw
Pos_Worst_btw_Exp_unique<-(mean(Pos_Worst_All_var$between_All_btw)-mean(Pos_Worst_c.DN_var$between_All_btw))*Pos_Worst_btw

Pos_Worst_btw_Shared<-Pos_Worst_btw_All-(Pos_Worst_btw_DN_unique+Pos_Worst_btw_Exp_unique)

#Within Variance Model Terms: 
Pos_Worst_wthn<-mean(Pos_Worst_cross_var$within_var)
Pos_Worst_wthn_Worst<-mean(Pos_Worst_All_var$within_fix_tot)
Pos_Worst_wthn_Worst_DN<-(mean(Pos_Worst_All_var$within_res_wthn)-mean(Pos_Worst_cross_var$within_res_wthn))*Pos_Worst_wthn
Pos_Worst_wthn_unmod<-Pos_Worst_wthn-Pos_Worst_wthn_Worst-Pos_Worst_wthn_Worst_DN
#Note final value produced above is a measure of total variance attributable to cross-level interaction

#Creating Total DN Effect for separate river plot decomposition: 
Total_DN<-Pos_Worst_btw_DN_unique+Pos_Worst_btw_Shared+Pos_Worst_wthn_Worst_DN

###########################################################################################################
#Main Variance Decompisition River Plot :: Positive Mood - Worst Event Models :: 
###########################################################################################################
Pos_Worst_River_DF<-data.frame(N1 = c('DN',
                                     'DN <--> NDE Exposure',
                                     'NDE Exposure', 
                                     'Unmodeled Between', 
                                     'Momentary NDE', 
                                     'DN x NDE Reactivity',
                                     'Unmodeled Within', 
                                     'Total Between', 
                                     'Total Within'), 
                              N2 = c('Total Between', 
                                     'Total Between', 
                                     'Total Between', 
                                     'Total Between', 
                                     'Total Within', 
                                     'Total Within',
                                     'Total Within', 
                                     'Total Variance', 
                                     'Total Variance'), 
                              Value = c(Pos_Worst_btw_DN_unique, 
                                        Pos_Worst_btw_Shared,
                                        Pos_Worst_btw_Exp_unique, 
                                        Pos_Worst_btw_umod,
                                        Pos_Worst_wthn_Worst, 
                                        Pos_Worst_wthn_Worst_DN,
                                        Pos_Worst_wthn_unmod,
                                        Pos_Worst_btw, 
                                        Pos_Worst_wthn), 
                              ID = 1:9)

Pos_Worst_River_DF$N1<-paste(Pos_Worst_River_DF$N1, 
                            '\n', 
                            paste0(round(Pos_Worst_River_DF$Value*100, 
                                         digits = 2), '%'
                            )
)

Pos_Worst_River_DF$N2<-c(rep(Pos_Worst_River_DF$N1[8], 4),
                        rep(Pos_Worst_River_DF$N1[9], 3), 
                        rep('Total Variance', 2)
)

nodes<-data.frame(ID = c(Pos_Worst_River_DF$N1, 
                         'Total Variance'), 
                  x = c(1,1,1,1,1,1,1,2,2,3), 
                  y = c(0,2,4,6,8,10,12,3,9,6))

palette = c(paste0(brewer.pal(9, "Blues"), 90)[c(2, 4, 6, 8)], 
            paste0(brewer.pal(9, "Reds"), 90)[c(4,6,8)], 
            paste0(brewer.pal(9, "Blues"), 90)[9], 
            paste0(brewer.pal(9, "Reds"), 90)[9], 
            paste0(brewer.pal(9, 'Purples'), 90)[9])

styles = lapply(nodes$y, function(n) {
  list(col = palette[n], lty = 0, textcol = "black")
})

#Not sure why color is not mapping correctly - trying to force it with this fuction... 
for(i in 1:length(palette)){
  styles[[i]]$col<-palette[i]
}

names(styles) = nodes$ID

Pos_Worst_riv<-makeRiver(nodes = nodes, 
                        edges =  Pos_Worst_River_DF, 
                        node_styles = styles)

png(paste0(study1.graphics, '/S1_Pos_Worst_river_overall.png'), 
    units = 'in', 
    res = 1200, 
    height = 14, 
    width = 8)
riverplot(Pos_Worst_riv, 
          nodewidth = 3, 
          plot_area = 1)
title(ylab = 'Study 1 - Riverplot of Total Variance Decomposition Estimated from Positive Mood Models: Negative Daily Events')
dev.off()

#----------------------------------------------------------------------------------------------------------
#DN - Specific Plot 
# :: Positive Mood - Worst Event Models :: 
#----------------------------------------------------------------------------------------------------------
Pos_Worst_DN_River_DF<-data.frame(N1 = c('DN',
                                        'DN <--> NDE Exposure',
                                        'DN x NDE Reactivity'
), 
N2 = rep(paste0("Combined DN Effect ",
                round(Total_DN*100, digits = 2),
                '%'), 
         3), 
Value = c(Pos_Worst_btw_DN_unique, 
          Pos_Worst_btw_Shared, 
          Pos_Worst_wthn_Worst_DN), 
ID = 1:3)

Pos_Worst_DN_River_DF$N1<-paste(Pos_Worst_DN_River_DF$N1, 
                               '\n', 
                               paste0(round(Pos_Worst_DN_River_DF$Value*100, 
                                            digits = 2), '%'
                               )
)

nodes<-data.frame(ID = c(Pos_Worst_DN_River_DF$N1, 
                         paste0("Combined DN Effect ",
                                round(Total_DN*100, digits = 2),
                                '%')), 
                  x = c(1,1,1,2), 
                  y = c(0,1,2,1))

palette = c(paste0(brewer.pal(9, "Blues"),90)[c(7,9)],
            paste0(brewer.pal(9, "Reds"),90)[7], 
            paste0(brewer.pal(9, "Purples"),90)[9])

styles = lapply(nodes$y, function(n) {
  list(col = palette[n], lty = 0, textcol = "black")
})

#Not sure why color is not mapping correctly - trying to force it with this fuction... 
for(i in 1:length(palette)){
  styles[[i]]$col<-palette[i]
}

names(styles) = nodes$ID

Pos_Worst_DN_riv<-makeRiver(nodes = nodes, 
                           edges =  Pos_Worst_DN_River_DF, 
                           node_styles = styles)

png(paste0(study1.graphics, '/S1_Pos_Worst_river_DN.png'), 
    units = 'in', 
    res = 900, 
    height = 10, 
    width = 10)
riverplot(Pos_Worst_DN_riv, 
          nodewidth = 3, 
          plot_area = .95)
title(ylab = 'Study 1 - Riverplot of DN Effect Decomposition Derived from Positive Mood Model: Positive Daily Events')
dev.off()
#Come back to this area to produce publication graphics (or arrange as graphics in ggplot)
#png(paste0(study1.graphics, '/Combined_Pos_Worst_Rivers.png'), 
#    units = 'in', 
#    res = 900, 
#    height = 20, 
#    width = 10)
#layout(matrix(c(1,2,2), nrow=3, ncol=1))
#par(mar = c(0,0,0,0))

#riverplot(Pos_Worst_DN_riv, 
#          nodewidth = 3, 
#          plot_area = .95)
#title(ylab = 'Study 1 - Riverplot of DN Effect Decomposition Derived from Positive Mood Model: Positive Daily Events')

#riverplot(Pos_Worst_riv, 
#          nodewidth = 3, 
#          plot_area = 1)
#title(ylab = 'Study 1 - Riverplot of Total Variance Decomposition Estimated from Positive Mood Models: Positive Daily Events')
#dev.off()

#==========================================================================================================
#POSITIVE MOOD MODELS - Variance Decomposition
# :: Positive Mood - Best Event Models ::
#==========================================================================================================

#==========================================================================================================
#----------------------------------------------------------------------------------------------------------
#Positive Mood, Positive Events, and Dispositional Negativity: 
#==========================================================================================================
post_iter<-2000   #No. of chains per data set x No. post-warmup iterations per chain
#for loops below require this value is correctly specified

#MODEL 1 - LV2 Predictor = DN, LV1 = individually mean-centered negative event indicator
#=-=-=-=-=-=-=-=-=-=-=-=-=-=
within_cov<-c(24)                    #Columns with group-mean centered predictors
between_cov<-c(8)                       #Columns with between-subject predictors
random_cov<-c(24)                    #Level 1 predictors with random slopes

#Getting distributions of fixed effects
Intercept = posterior_samples(Pos_Best_c.DN, pars = 'b_Intercept') 
DN = posterior_samples(Pos_Best_c.DN, pars = 'b_c.DN') 
Best = posterior_samples(Pos_Best_c.DN, pars = 'b_c.Best') 

#Getting between-subjects variance terms for random effects
Int_var = posterior_samples(Pos_Best_c.DN, pars = 'sd_ID__Intercept')^2
Best_var = posterior_samples(Pos_Best_c.DN, pars = 'sd_ID__c.Best')^2 

#Getting covariances of random effects (i.e., tau matrix)
cov_Int_Best = posterior_samples(Pos_Best_c.DN, pars = 'sd_ID__Intercept')*
  posterior_samples(Pos_Best_c.DN, pars = 'sd_ID__c.Best')*
  posterior_samples(Pos_Best_c.DN, pars = 'cor_ID__Intercept__c.Best')

#Getting level 1 error variance
sigma<-posterior_samples(Pos_Best_c.DN, pars = 'sigma')

post_samples<-data.frame(Intercept,
                         DN, 
                         Best, 
                         Int_var, 
                         Best_var, 
                         cov_Int_Best, 
                         sigma)

colnames(post_samples)<-c('Intercept', 
                          'DN', 
                          'Best', 
                          'Int_var', 
                          'Best_var', 
                          'cov_Int_Best', 
                          'sigma')

#Aggregate across imputed datasets
#Currently going to take 1000 draws from posterior distributions
#Will then apply across all 20 data sets

Pos_Best_c.DN_var<-list(between_var=vector(), 
                        within_var=vector(), 
                        between_All_tot=vector(),
                        between_All_btw=vector(),
                        between_res_btw=vector(),
                        within_fix_wthn=vector(),
                        within_fix_tot=vector(),
                        within_slope_var_wthn=vector(),
                        within_res_wthn=vector(), 
                        within_unmod_tot=vector())

for(i in 1:length(dat.study1.list)){
  print(paste('Starting w/ imputed set', i, 'out of', M))
  samp_frame<-1:post_iter+(i-1)*post_iter
  D<-sample(samp_frame, size = 400, replace = FALSE)
  for(d in 1:length(D)){
    #browser()
    Gamma_w<-as.vector(
      as.matrix(
        post_samples[D[d], c('Best')]))       #Make sure the effects line up - in order of within_cov
    Gamma_b<-as.vector(
      as.matrix(
        post_samples[D[d], c('Intercept', 'DN')]))    #level fixed intercept and any level 2 fixed effects
    tau<-c(post_samples[D[d], c('Int_var', 'cov_Int_Best')], 
           post_samples[D[d], c('cov_Int_Best', 'Best_var')])
    
    tau<-matrix(unlist(tau), 
                byrow=TRUE, 
                ncol=2)
    
    #Matrix rows/columns ordered starting with tau.00 (intercept variance)
    #then add columns/rows for random effects of slopes in order of within_cov
    
    Sigma2<-post_samples[D[d], 'sigma']
    
    r2mlm_DN<-r2MLM(data=dat.study1.list[[i]], 
                    within_covs = within_cov, 
                    between_covs = between_cov, 
                    random_covs = random_cov, 
                    gamma_w = Gamma_w, 
                    gamma_b = Gamma_b, 
                    Tau = tau, 
                    sigma2 = Sigma2, 
                    has_intercept = TRUE, 
                    clustermeancentered = TRUE)
    
    #Extracting relevant values across imputed datasets
    Pos_Best_c.DN_var$between_var<-c(Pos_Best_c.DN_var$between_var, as.numeric(r2mlm_DN$Decompositions['fixed, between', 'total'])+
                                       as.numeric(r2mlm_DN$Decompositions['mean variation', 'total']))
    Pos_Best_c.DN_var$within_var<-c(Pos_Best_c.DN_var$within_var, as.numeric(r2mlm_DN$Decompositions['sigma2', 'total'])+
                                      as.numeric(r2mlm_DN$Decompositions['slope variation', 'total'])+
                                      as.numeric(r2mlm_DN$Decompositions['fixed, within', 'total']))
    Pos_Best_c.DN_var$between_All_tot<-c(Pos_Best_c.DN_var$between_All_tot, as.numeric(r2mlm_DN$R2s['f2', 'total']))
    Pos_Best_c.DN_var$between_All_btw<-c(Pos_Best_c.DN_var$between_All_btw, as.numeric(r2mlm_DN$R2s['f2', 'between']))
    Pos_Best_c.DN_var$between_res_btw<-c(Pos_Best_c.DN_var$between_res_btw, as.numeric(r2mlm_DN$R2s['m', 'between']))
    Pos_Best_c.DN_var$within_fix_wthn<-c(Pos_Best_c.DN_var$within_fix_wthn, as.numeric(r2mlm_DN$Decompositions['fixed, within', 'within']))
    Pos_Best_c.DN_var$within_fix_tot<-c(Pos_Best_c.DN_var$within_fix_tot, as.numeric(r2mlm_DN$Decompositions['fixed, within', 'total']))
    Pos_Best_c.DN_var$within_slope_var_wthn<-c(Pos_Best_c.DN_var$within_slope_var_wthn, as.numeric(r2mlm_DN$Decompositions['slope variation', 'within']))
    Pos_Best_c.DN_var$within_res_wthn<-c(Pos_Best_c.DN_var$within_res_wthn, as.numeric(r2mlm_DN$Decompositions['sigma2', 'within']))
    Pos_Best_c.DN_var$within_unmod_tot<-c(Pos_Best_c.DN_var$within_unmod_tot, as.numeric(r2mlm_DN$Decompositions['sigma2', 'total'])+
                                            as.numeric(r2mlm_DN$Decompositions['slope variation', 'total']))
    
    print(paste("Sampling", d, "out of", length(D), '- Dataset:', i))
  }
  print(paste("*** Finshed imputed data set", i, 'out of', M, "***"))
}
hist(Pos_Best_c.DN_var$between_var)
mean(Pos_Best_c.DN_var$between_var)

hist(Pos_Best_c.DN_var$within_var)

#MODEL 2 - LV2 Predictor = Mean Positive Events, LV1 = individually mean-centered negative event indicator
#=-=-=-=-=-=-=-=-=-=-=-=-=-=
within_cov<-c(24)                    #Columns with group-mean centered predictors
between_cov<-c(25)                       #Columns with between-subject predictors
random_cov<-c(24)                    #Level 1 predictors with random slopes

#Getting distributions of fixed effects
Intercept = posterior_samples(Pos_Best_prop_Best, pars = 'b_Intercept') 
prop_Best = posterior_samples(Pos_Best_prop_Best, pars = 'b_prop_Best') 
Best = posterior_samples(Pos_Best_prop_Best, pars = 'b_c.Best') 

#Getting between-subjects variance terms for random effects
Int_var = posterior_samples(Pos_Best_prop_Best, pars = 'sd_ID__Intercept')^2
Best_var = posterior_samples(Pos_Best_prop_Best, pars = 'sd_ID__c.Best')^2 

#Getting covariances of random effects (i.e., tau matrix)
cov_Int_Best = posterior_samples(Pos_Best_prop_Best, pars = 'sd_ID__Intercept')*
  posterior_samples(Pos_Best_prop_Best, pars = 'sd_ID__c.Best')*
  posterior_samples(Pos_Best_prop_Best, pars = 'cor_ID__Intercept__c.Best')

#Getting level 1 error variance
sigma<-posterior_samples(Pos_Best_prop_Best, pars = 'sigma')

post_samples<-data.frame(Intercept,
                         prop_Best, 
                         Best, 
                         Int_var, 
                         Best_var, 
                         cov_Int_Best, 
                         sigma)

colnames(post_samples)<-c('Intercept', 
                          'prop_Best', 
                          'Best', 
                          'Int_var', 
                          'Best_var', 
                          'cov_Int_Best', 
                          'sigma')

#Aggregate across imputed datasets
#Currently going to take 1000 draws from posterior distributions
#Will then apply across all 20 data sets

Pos_Best_prop_Best_var<-list(between_var=vector(), 
                             within_var=vector(), 
                             between_All_tot=vector(),
                             between_All_btw=vector(),
                             between_res_btw=vector(),
                             within_fix_wthn=vector(),
                             within_fix_tot=vector(),
                             within_slope_var_wthn=vector(),
                             within_res_wthn=vector(), 
                             within_unmod_tot=vector())

for(i in 1:length(dat.study1.list)){
  print(paste('Starting w/ imputed set', i, 'out of', M))
  samp_frame<-1:post_iter+(i-1)*post_iter
  D<-sample(samp_frame, size = 400, replace = FALSE)
  for(d in 1:length(D)){
    #browser()
    Gamma_w<-as.vector(
      as.matrix(
        post_samples[D[d], c('Best')]))       #Make sure the effects line up - in order of within_cov
    Gamma_b<-as.vector(
      as.matrix(
        post_samples[D[d], c('Intercept', 'prop_Best')]))    #level fixed intercept and any level 2 fixed effects
    tau<-c(post_samples[D[d], c('Int_var', 'cov_Int_Best')], 
           post_samples[D[d], c('cov_Int_Best', 'Best_var')])
    
    tau<-matrix(unlist(tau), 
                byrow=TRUE, 
                ncol=2)
    
    #Matrix rows/columns ordered starting with tau.00 (intercept variance)
    #then add columns/rows for random effects of slopes in order of within_cov
    
    Sigma2<-post_samples[D[d], 'sigma']
    
    r2mlm_prop_Best<-r2MLM(data=dat.study1.list[[i]], 
                           within_covs = within_cov, 
                           between_covs = between_cov, 
                           random_covs = random_cov, 
                           gamma_w = Gamma_w, 
                           gamma_b = Gamma_b, 
                           Tau = tau, 
                           sigma2 = Sigma2, 
                           has_intercept = TRUE, 
                           clustermeancentered = TRUE)
    
    #Extracting relevant values across imputed datasets
    Pos_Best_prop_Best_var$between_var<-c(Pos_Best_prop_Best_var$between_var, as.numeric(r2mlm_prop_Best$Decompositions['fixed, between', 'total'])+
                                            as.numeric(r2mlm_prop_Best$Decompositions['mean variation', 'total']))
    Pos_Best_prop_Best_var$within_var<-c(Pos_Best_prop_Best_var$within_var, as.numeric(r2mlm_prop_Best$Decompositions['sigma2', 'total'])+
                                           as.numeric(r2mlm_prop_Best$Decompositions['slope variation', 'total'])+
                                           as.numeric(r2mlm_prop_Best$Decompositions['fixed, within', 'total']))
    Pos_Best_prop_Best_var$between_All_tot<-c(Pos_Best_prop_Best_var$between_All_tot, as.numeric(r2mlm_prop_Best$R2s['f2', 'total']))
    Pos_Best_prop_Best_var$between_All_btw<-c(Pos_Best_prop_Best_var$between_All_btw, as.numeric(r2mlm_prop_Best$R2s['f2', 'between']))
    Pos_Best_prop_Best_var$between_res_btw<-c(Pos_Best_prop_Best_var$between_res_btw, as.numeric(r2mlm_prop_Best$R2s['m', 'between']))
    Pos_Best_prop_Best_var$within_fix_wthn<-c(Pos_Best_prop_Best_var$within_fix_wthn, as.numeric(r2mlm_prop_Best$Decompositions['fixed, within', 'within']))
    Pos_Best_prop_Best_var$within_fix_tot<-c(Pos_Best_prop_Best_var$within_fix_tot, as.numeric(r2mlm_prop_Best$Decompositions['fixed, within', 'total']))
    Pos_Best_prop_Best_var$within_slope_var_wthn<-c(Pos_Best_prop_Best_var$within_slope_var_wthn, as.numeric(r2mlm_prop_Best$Decompositions['slope variation', 'within']))
    Pos_Best_prop_Best_var$within_res_wthn<-c(Pos_Best_prop_Best_var$within_res_wthn, as.numeric(r2mlm_prop_Best$Decompositions['sigma2', 'within']))
    Pos_Best_prop_Best_var$within_unmod_tot<-c(Pos_Best_prop_Best_var$within_unmod_tot, as.numeric(r2mlm_prop_Best$Decompositions['sigma2', 'total'])+
                                                 as.numeric(r2mlm_prop_Best$Decompositions['slope variation', 'total']))
    
    print(paste("Sampling", d, "out of", length(D), '- Dataset:', i))
  }
  print(paste("*** Finshed imputed data set", i, 'out of', M, "***"))
}
hist(Pos_Best_prop_Best_var$between_var)
mean(Pos_Best_prop_Best_var$between_var)
hist(Pos_Best_prop_Best_var$within_var)

#MODEL 3 - LV2 Predictors = DN and Mean Positive Events, LV1 = individually mean-centered negative event indicator
#=-=-=-=-=-=-=-=-=-=-=-=-=-=
within_cov<-c(24)                    #Columns with group-mean centered predictors
between_cov<-c(8,25)                       #Columns with between-subject predictors
random_cov<-c(24)                    #Level 1 predictors with random slopes

#Getting distributions of fixed effects
Intercept = posterior_samples(Pos_Best_All, pars = 'b_Intercept') 
DN = posterior_samples(Pos_Best_All, pars = 'b_c.DN') 
prop_Best = posterior_samples(Pos_Best_All, pars = 'b_prop_Best') 
Best = posterior_samples(Pos_Best_All, pars = 'b_c.Best') 

#Getting between-subjects variance terms for random effects
Int_var = posterior_samples(Pos_Best_All, pars = 'sd_ID__Intercept')^2
Best_var = posterior_samples(Pos_Best_All, pars = 'sd_ID__c.Best')^2 

#Getting covariances of random effects (i.e., tau matrix)
cov_Int_Best = posterior_samples(Pos_Best_All, pars = 'sd_ID__Intercept')*
  posterior_samples(Pos_Best_All, pars = 'sd_ID__c.Best')*
  posterior_samples(Pos_Best_All, pars = 'cor_ID__Intercept__c.Best')

#Getting level 1 error variance
sigma<-posterior_samples(Pos_Best_All, pars = 'sigma')

post_samples<-data.frame(Intercept,
                         DN,
                         prop_Best, 
                         Best, 
                         Int_var, 
                         Best_var, 
                         cov_Int_Best, 
                         sigma)

colnames(post_samples)<-c('Intercept', 
                          'DN',
                          'prop_Best', 
                          'Best', 
                          'Int_var', 
                          'Best_var', 
                          'cov_Int_Best', 
                          'sigma')

# ***** IMPORTANT *****
# REMOVE ALL SECTIONS LIKE THE ONE BELOW
# ***** IMPORTANT *****
#Aggregate across imputed datasets
#Currently going to take 1000 draws from posterior distributions
#Will then apply across all 20 data sets

Pos_Best_All_var<-list(between_var=vector(), 
                       within_var=vector(), 
                       between_All_tot=vector(),
                       between_All_btw=vector(),
                       between_res_btw=vector(),
                       within_fix_wthn=vector(),
                       within_fix_tot=vector(),
                       within_slope_var_wthn=vector(),
                       within_res_wthn=vector(), 
                       within_unmod_tot=vector())

for(i in 1:length(dat.study1.list)){
  print(paste('Starting w/ imputed set', i, 'out of', M))
  samp_frame<-1:post_iter+(i-1)*post_iter
  D<-sample(samp_frame, size = 400, replace = FALSE)
  for(d in 1:length(D)){
    #browser()
    Gamma_w<-as.vector(
      as.matrix(
        post_samples[D[d], c('Best')]))       #Make sure the effects line up - in order of within_cov
    Gamma_b<-as.vector(
      as.matrix(
        post_samples[D[d], c('Intercept', 'DN', 'prop_Best')]))    #level fixed intercept and any level 2 fixed effects
    tau<-c(post_samples[D[d], c('Int_var', 'cov_Int_Best')], 
           post_samples[D[d], c('cov_Int_Best', 'Best_var')])
    
    tau<-matrix(unlist(tau), 
                byrow=TRUE, 
                ncol=2)
    
    #Matrix rows/columns ordered starting with tau.00 (intercept variance)
    #then add columns/rows for random effects of slopes in order of within_cov
    
    Sigma2<-post_samples[D[d], 'sigma']
    
    r2mlm_prop_Best<-r2MLM(data=dat.study1.list[[i]], 
                           within_covs = within_cov, 
                           between_covs = between_cov, 
                           random_covs = random_cov, 
                           gamma_w = Gamma_w, 
                           gamma_b = Gamma_b, 
                           Tau = tau, 
                           sigma2 = Sigma2, 
                           has_intercept = TRUE, 
                           clustermeancentered = TRUE)
    
    #Extracting relevant values across imputed datasets
    Pos_Best_All_var$between_var<-c(Pos_Best_All_var$between_var, as.numeric(r2mlm_prop_Best$Decompositions['fixed, between', 'total'])+
                                      as.numeric(r2mlm_prop_Best$Decompositions['mean variation', 'total']))
    Pos_Best_All_var$within_var<-c(Pos_Best_All_var$within_var, as.numeric(r2mlm_prop_Best$Decompositions['sigma2', 'total'])+
                                     as.numeric(r2mlm_prop_Best$Decompositions['slope variation', 'total'])+
                                     as.numeric(r2mlm_prop_Best$Decompositions['fixed, within', 'total']))
    Pos_Best_All_var$between_All_tot<-c(Pos_Best_All_var$between_All_tot, as.numeric(r2mlm_prop_Best$R2s['f2', 'total']))
    Pos_Best_All_var$between_All_btw<-c(Pos_Best_All_var$between_All_btw, as.numeric(r2mlm_prop_Best$R2s['f2', 'between']))
    Pos_Best_All_var$between_res_btw<-c(Pos_Best_All_var$between_res_btw, as.numeric(r2mlm_prop_Best$R2s['m', 'between']))
    Pos_Best_All_var$within_fix_wthn<-c(Pos_Best_All_var$within_fix_wthn, as.numeric(r2mlm_prop_Best$Decompositions['fixed, within', 'within']))
    Pos_Best_All_var$within_fix_tot<-c(Pos_Best_All_var$within_fix_tot, as.numeric(r2mlm_prop_Best$Decompositions['fixed, within', 'total']))
    Pos_Best_All_var$within_slope_var_wthn<-c(Pos_Best_All_var$within_slope_var_wthn, as.numeric(r2mlm_prop_Best$Decompositions['slope variation', 'within']))
    Pos_Best_All_var$within_res_wthn<-c(Pos_Best_All_var$within_res_wthn, as.numeric(r2mlm_prop_Best$Decompositions['sigma2', 'within']))
    Pos_Best_All_var$within_unmod_tot<-c(Pos_Best_All_var$within_unmod_tot, as.numeric(r2mlm_prop_Best$Decompositions['sigma2', 'total'])+
                                           as.numeric(r2mlm_prop_Best$Decompositions['slope variation', 'total']))
    
    print(paste("Sampling", d, "out of", length(D), '- Dataset:', i))
  }
  print(paste("*** Finshed imputed data set", i, 'out of', M, "***"))
}
hist(Pos_Best_All_var$between_var)
mean(Pos_Best_All_var$between_var)
hist(Pos_Best_All_var$within_var)

#MODEL 4 - Full cross-level interaction with 
#=-=-=-=-=-=-=-=-=-=-=-=-=-=
within_cov<-c(24,27)                    #Columns with group-mean centered predictors
between_cov<-c(8,25)                       #Columns with between-subject predictors
random_cov<-c(24)                    #Level 1 predictors with random slopes

#Getting distributions of fixed effects
Intercept = posterior_samples(Pos_Best_cross, pars = 'b_Intercept') 
DN = posterior_samples(Pos_Best_cross, pars = 'b_c.DN') 
prop_Best = posterior_samples(Pos_Best_cross, pars = 'b_prop_Best') 
Best = posterior_samples(Pos_Best_cross, pars = 'b_c.Best')[,1]  #Note posterior_samples() searches for text match
#With cross-level interaction returns two values 
#Need to select the one with the fixed effect for lv-1
DNxBest = posterior_samples(Pos_Best_cross, pars = 'b_c.Best:c.DN')

#Getting between-subjects variance terms for random effects
Int_var = posterior_samples(Pos_Best_cross, pars = 'sd_ID__Intercept')^2
Best_var = posterior_samples(Pos_Best_cross, pars = 'sd_ID__c.Best')^2 

#Getting covariances of random effects (i.e., tau matrix)
cov_Int_Best = posterior_samples(Pos_Best_cross, pars = 'sd_ID__Intercept')*
  posterior_samples(Pos_Best_cross, pars = 'sd_ID__c.Best')*
  posterior_samples(Pos_Best_cross, pars = 'cor_ID__Intercept__c.Best')

#Getting level 1 error variance
sigma<-posterior_samples(Pos_Best_cross, pars = 'sigma')

post_samples<-data.frame(Intercept,
                         DN,
                         prop_Best, 
                         Best,
                         DNxBest,
                         Int_var, 
                         Best_var, 
                         cov_Int_Best, 
                         sigma)

colnames(post_samples)<-c('Intercept', 
                          'DN',
                          'prop_Best', 
                          'Best', 
                          'DNxBest',
                          'Int_var', 
                          'Best_var', 
                          'cov_Int_Best', 
                          'sigma')

#Aggregate across imputed datasets
#Currently going to take 1000 draws from posterior distributions
#Will then apply across all 20 data sets

Pos_Best_cross_var<-list(between_var=vector(), 
                         within_var=vector(), 
                         between_All_tot=vector(),
                         between_All_btw=vector(),
                         between_res_btw=vector(),
                         within_fix_wthn=vector(),
                         within_fix_tot=vector(),
                         within_slope_var_wthn=vector(),
                         within_res_wthn=vector(), 
                         within_unmod_tot=vector())

for(i in 1:length(dat.study1.list)){
  print(paste('Starting w/ imputed set', i, 'out of', M))
  samp_frame<-1:post_iter+(i-1)*post_iter
  D<-sample(samp_frame, size = 400, replace = FALSE)
  for(d in 1:length(D)){
    #browser()
    Gamma_w<-as.vector(
      as.matrix(
        post_samples[D[d], c('Best', 'DNxBest')]))       #Make sure the effects line up - in order of within_cov
    Gamma_b<-as.vector(
      as.matrix(
        post_samples[D[d], c('Intercept', 'DN', 'prop_Best')]))    #level fixed intercept and any level 2 fixed effects
    tau<-c(post_samples[D[d], c('Int_var', 'cov_Int_Best')], 
           post_samples[D[d], c('cov_Int_Best', 'Best_var')])
    
    tau<-matrix(unlist(tau), 
                byrow=TRUE, 
                ncol=2)
    
    #Matrix rows/columns ordered starting with tau.00 (intercept variance)
    #then add columns/rows for random effects of slopes in order of within_cov
    
    Sigma2<-post_samples[D[d], 'sigma']
    
    r2mlm_prop_Best<-r2MLM(data=dat.study1.list[[i]], 
                           within_covs = within_cov, 
                           between_covs = between_cov, 
                           random_covs = random_cov, 
                           gamma_w = Gamma_w, 
                           gamma_b = Gamma_b, 
                           Tau = tau, 
                           sigma2 = Sigma2, 
                           has_intercept = TRUE, 
                           clustermeancentered = TRUE)
    
    #Extracting relevant values across imputed datasets
    Pos_Best_cross_var$between_var<-c(Pos_Best_cross_var$between_var, as.numeric(r2mlm_prop_Best$Decompositions['fixed, between', 'total'])+
                                        as.numeric(r2mlm_prop_Best$Decompositions['mean variation', 'total']))
    Pos_Best_cross_var$within_var<-c(Pos_Best_cross_var$within_var, as.numeric(r2mlm_prop_Best$Decompositions['sigma2', 'total'])+
                                       as.numeric(r2mlm_prop_Best$Decompositions['slope variation', 'total'])+
                                       as.numeric(r2mlm_prop_Best$Decompositions['fixed, within', 'total']))
    Pos_Best_cross_var$between_All_tot<-c(Pos_Best_cross_var$between_All_tot, as.numeric(r2mlm_prop_Best$R2s['f2', 'total']))
    Pos_Best_cross_var$between_All_btw<-c(Pos_Best_cross_var$between_All_btw, as.numeric(r2mlm_prop_Best$R2s['f2', 'between']))
    Pos_Best_cross_var$between_res_btw<-c(Pos_Best_cross_var$between_res_btw, as.numeric(r2mlm_prop_Best$R2s['m', 'between']))
    Pos_Best_cross_var$within_fix_wthn<-c(Pos_Best_cross_var$within_fix_wthn, as.numeric(r2mlm_prop_Best$Decompositions['fixed, within', 'within']))
    Pos_Best_cross_var$within_fix_tot<-c(Pos_Best_cross_var$within_fix_tot, as.numeric(r2mlm_prop_Best$Decompositions['fixed, within', 'total']))
    Pos_Best_cross_var$within_slope_var_wthn<-c(Pos_Best_cross_var$within_slope_var_wthn, as.numeric(r2mlm_prop_Best$Decompositions['slope variation', 'within']))
    Pos_Best_cross_var$within_res_wthn<-c(Pos_Best_cross_var$within_res_wthn, as.numeric(r2mlm_prop_Best$Decompositions['sigma2', 'within']))
    Pos_Best_cross_var$within_unmod_tot<-c(Pos_Best_cross_var$within_unmod_tot, as.numeric(r2mlm_prop_Best$Decompositions['sigma2', 'total'])+
                                             as.numeric(r2mlm_prop_Best$Decompositions['slope variation', 'total']))
    
    print(paste("Sampling", d, "out of", length(D), '- Dataset:', i))
  }
  print(paste("*** Finshed imputed data set", i, 'out of', M, "***"))
}
hist(Pos_Best_cross_var$between_var)
hist(Pos_Best_cross_var$within_var)

###########################################################################################################
#Preparing values for River plots :: Positive Mood - Best Event Models ::
###########################################################################################################
#Checking across models - between subjects variance
mean(Pos_Best_c.DN_var$between_var)
mean(Pos_Best_prop_Best_var$between_var)
mean(Pos_Best_All_var$between_var)
mean(Pos_Best_cross_var$between_var)

#Differemce between R2's due to sampling variability - is only 0.0008066132
#Difference is <.1% in total variance
mean(Pos_Best_cross_var$between_var)-mean(Pos_Best_c.DN_var$between_var)

#Between Variance Model Terms: 
Pos_Best_btw<-mean(Pos_Best_cross_var$between_var)
Pos_Best_btw_All<-mean(Pos_Best_cross_var$between_All_tot)
Pos_Best_btw_umod<-Pos_Best_btw-Pos_Best_btw_All

Pos_Best_btw_DN_unique<-(mean(Pos_Best_All_var$between_All_btw)-mean(Pos_Best_prop_Best_var$between_All_btw))*Pos_Best_btw
Pos_Best_btw_Exp_unique<-(mean(Pos_Best_All_var$between_All_btw)-mean(Pos_Best_c.DN_var$between_All_btw))*Pos_Best_btw

Pos_Best_btw_Shared<-Pos_Best_btw_All-(Pos_Best_btw_DN_unique+Pos_Best_btw_Exp_unique)

#Within Variance Model Terms: 
Pos_Best_wthn<-mean(Pos_Best_cross_var$within_var)
Pos_Best_wthn_Best<-mean(Pos_Best_All_var$within_fix_tot)
Pos_Best_wthn_Best_DN<-(mean(Pos_Best_All_var$within_res_wthn)-mean(Pos_Best_cross_var$within_res_wthn))*Pos_Best_wthn
Pos_Best_wthn_unmod<-Pos_Best_wthn-Pos_Best_wthn_Best-Pos_Best_wthn_Best_DN
#Note final value produced above is a measure of total variance attributable to cross-level interaction

#Creating Total DN Effect for separate river plot decomposition: 
Total_DN<-Pos_Best_btw_DN_unique+Pos_Best_btw_Shared+Pos_Best_wthn_Best_DN

###########################################################################################################
#Main Variance Decompisition River Plot :: Positive Mood - Best Event Models :: 
###########################################################################################################
Pos_Best_River_DF<-data.frame(N1 = c('DN',
                                     'DN <--> PDE Exposure',
                                     'PDE Exposure', 
                                     'Unmodeled Between', 
                                     'Momentary PDE', 
                                     'DN x PDE Reactivity',
                                     'Unmodeled Within', 
                                     'Total Between', 
                                     'Total Within'), 
                              N2 = c('Total Between', 
                                     'Total Between', 
                                     'Total Between', 
                                     'Total Between', 
                                     'Total Within', 
                                     'Total Within',
                                     'Total Within', 
                                     'Total Variance', 
                                     'Total Variance'), 
                              Value = c(Pos_Best_btw_DN_unique, 
                                        Pos_Best_btw_Shared,
                                        Pos_Best_btw_Exp_unique, 
                                        Pos_Best_btw_umod,
                                        Pos_Best_wthn_Best, 
                                        Pos_Best_wthn_Best_DN,
                                        Pos_Best_wthn_unmod,
                                        Pos_Best_btw, 
                                        Pos_Best_wthn), 
                              ID = 1:9)

Pos_Best_River_DF$N1<-paste(Pos_Best_River_DF$N1, 
                            '\n', 
                            paste0(round(Pos_Best_River_DF$Value*100, 
                                         digits = 2), '%'
                            )
)

Pos_Best_River_DF$N2<-c(rep(Pos_Best_River_DF$N1[8], 4),
                        rep(Pos_Best_River_DF$N1[9], 3), 
                        rep('Total Variance', 2)
)

nodes<-data.frame(ID = c(Pos_Best_River_DF$N1, 
                         'Total Variance'), 
                  x = c(1,1,1,1,1,1,1,2,2,3), 
                  y = c(0,2,4,6,8,10,12,3,9,6))

palette = c(paste0(brewer.pal(9, "Blues"), 90)[c(2, 4, 6, 8)], 
            paste0(brewer.pal(9, "Reds"), 90)[c(4,6,8)], 
            paste0(brewer.pal(9, "Blues"), 90)[9], 
            paste0(brewer.pal(9, "Reds"), 90)[9], 
            paste0(brewer.pal(9, 'Purples'), 90)[9])

styles = lapply(nodes$y, function(n) {
  list(col = palette[n], lty = 0, textcol = "black")
})

#Not sure why color is not mapping correctly - trying to force it with this fuction... 
for(i in 1:length(palette)){
  styles[[i]]$col<-palette[i]
}

names(styles) = nodes$ID

Pos_Best_riv<-makeRiver(nodes = nodes, 
                        edges =  Pos_Best_River_DF, 
                        node_styles = styles)

#Creating Dataset to add in text - Will have to manipulate the y-values to get text to line up 

#Creating Dataset to add in text - Will have to manipulate the y-values to get text to line up 

png(paste0(study1.graphics, '/S1_Pos_Best_river_overall.png'), 
    units = 'in', 
    res = 1200, 
    height = 14, 
    width = 8)
riverplot(Pos_Best_riv, 
          nodewidth = 3, 
          plot_area = 1)
title(ylab = 'Study 1 - Riverplot of Total Variance Decomposition Estimated from Positive Mood Models: Positive Daily Events')
dev.off()

#----------------------------------------------------------------------------------------------------------
#DN - Specific Plot 
# :: Positive Mood - Best Event Models :: 
#----------------------------------------------------------------------------------------------------------
Pos_Best_DN_River_DF<-data.frame(N1 = c('DN',
                                        'DN <--> PDE Exposure',
                                        'DN x PDE Reactivity'
), 
N2 = rep(paste0("Combined DN Effect ",
                round(Total_DN*100, digits = 2),
                '%'), 
         3), 
Value = c(Pos_Best_btw_DN_unique, 
          Pos_Best_btw_Shared, 
          Pos_Best_wthn_Best_DN), 
ID = 1:3)

Pos_Best_DN_River_DF$N1<-paste(Pos_Best_DN_River_DF$N1, 
                               '\n', 
                               paste0(round(Pos_Best_DN_River_DF$Value*100, 
                                            digits = 2), '%'
                               )
)

nodes<-data.frame(ID = c(Pos_Best_DN_River_DF$N1, 
                         paste0("Combined DN Effect ",
                                round(Total_DN*100, digits = 2),
                                '%')), 
                  x = c(1,1,1,2), 
                  y = c(0,1,2,1))

palette = c(paste0(brewer.pal(9, "Blues"),90)[c(7,9)],
            paste0(brewer.pal(9, "Reds"),90)[7], 
            paste0(brewer.pal(9, "Purples"),90)[9])

styles = lapply(nodes$y, function(n) {
  list(col = palette[n], lty = 0, textcol = "black")
})

#Not sure why color is not mapping correctly - trying to force it with this fuction... 
for(i in 1:length(palette)){
  styles[[i]]$col<-palette[i]
}

names(styles) = nodes$ID

Pos_Best_DN_riv<-makeRiver(nodes = nodes, 
                           edges =  Pos_Best_DN_River_DF, 
                           node_styles = styles)

png(paste0(study1.graphics, '/S1_Pos_Best_river_DN.png'), 
    units = 'in', 
    res = 900, 
    height = 10, 
    width = 10)
riverplot(Pos_Best_DN_riv, 
          nodewidth = 3, 
          plot_area = .95)
title(ylab = 'Study 1 - Riverplot of DN Effect Decomposition Derived from Positive Mood Model: Positive Daily Events')
dev.off()


###########################################################################################################
#Simple Exposure Question 
###########################################################################################################

Worst_Evnt<-brm_multiple(Worst_dich~1+c.DN + (1|ID),
                          data = dat.study1.list, 
                         family = 'bernoulli',
                         prior = c(set_prior('student_t(3,0,10)', class='sd'), 
                                   set_prior('normal(0,5)', class='Intercept'), 
                                   set_prior('normal(0,5)', class='b')),
                         warmup = 3000, 
                         iter = 4000,
                         chains = 2,
                         control = list(adapt_delta=.95, 
                                        max_treedepth=15))
summary(Worst_Evnt)

Best_Evnt<-brm_multiple(Best_dich~1+c.DN + (1|ID),
                         data = dat.study1.list, 
                         family = 'bernoulli',
                         prior = c(set_prior('student_t(3,0,10)', class='sd'), 
                                   set_prior('normal(0,5)', class='Intercept'), 
                                   set_prior('normal(0,5)', class='b')),
                         warmup = 3000, 
                         iter = 4000,
                         chains = 2,
                         control = list(adapt_delta=.95, 
                                        max_treedepth=15))
summary(Best_Evnt)

