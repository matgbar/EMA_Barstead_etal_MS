############################################################################################################
#Authors:     Barstead, M.G., DeYoung, K.D., Anderson, A. S., & Shackman, A. J.

#Title:       The moment-to-moment affective experience of dispositionally negative individuals 

#Contact:     barstead@umd.edu

#Contents:    Bayesian regression analyses - Study 1
############################################################################################################
############################################################################################################
library(brms)
library(rstan)
library(bayesplot)
library(pan)
library(mitml)
library(mice)
library(parallel)

#Stan options
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
############################################################################################################
wd<-'/home/mbarsted/Dropbox/EMA_MS_drafts'
data.folder<-paste0(wd, '/Model_Data')
study1.out<-paste0(wd, '/Study1')
study1.graphics<-paste0(study1.out, '/Graphics')
study1.model<-paste0(study1.out, '/Model_Summaries')
stan.code<-paste0(wd, '/Stan_code')

#Loading Study 1 Data (from Emotion MS - Shackman et al. 2017)
load(paste0(data.folder, '/Emotion MS environment.RData'))
dat$c.AP_both<-dat$AP_Both-mean(dat$AP_Both, na.rm=T) #mean centering level 2 pred

############################################################################################################
#Creating a Function to Finalize Imputed Model Output: 
bayes.to.txt<-function(model=NULL, 
                       out.folder=NULL, 
                       file=NULL, 
                       DF=NULL, 
                       tot.pars=3, 
                       impute=TRUE){
  sink(paste0(out.folder, '/', file, '.txt'))
  cat('System Information:')
  cat('\n=====================================================================================')
  cat(paste0('\nProcessor:', '\t\t\t', benchmarkme::get_cpu()$model_name))
  cat(paste0('\nNumber of Threads:', '\t\t', parallel::detectCores(logical=T)))
  cat(paste0('\nRAM:', '\t\t\t\t', paste(round(benchmarkme::get_ram()/1073741824), 'GB')))
  cat('\n=====================================================================================\n')
  
  cat('\n\nModel Information:')
  cat('\n=====================================================================================')
  cat('\n\nFormula (lme4 syntax):')
  cat('\n-------------------------------------------------------------------------------------\n')
  print(model$formula)
  cat('\n\nPriors:')
  cat('\n-------------------------------------------------------------------------------------\n')
  print(prior_summary(model))
  cat('\n\nStan Code:')
  cat('\n-------------------------------------------------------------------------------------\n')
  print(make_stancode(model$formula, data=DF, family=model$family$family))
  cat('\n=====================================================================================\n')
  
  cat('\n\nStan Arguments:')
  cat('\n=====================================================================================')
  cat(paste0('\nAdapt Delta:', '\t\t\t', model$fit@stan_args[[1]]$control$adapt_delta))
  cat(paste0('\nMaximum Tree Depth:', '\t\t', model$fit@stan_args[[1]]$control$max_treedepth))
  cat('\n=====================================================================================\n')
  
  cat('\n\nVariance Explained:')
  cat('\n=====================================================================================')
  cat('\n\nResidual Variance:')
  cat('\n-------------------------------------------------------------------------------------\n')
  print(sjstats::icc(model, posterior = T), prob =.95, digits=3)
  cat('\n\nBayesian R-squared (overall variance explained):')
  cat('\n-------------------------------------------------------------------------------------\n')
  print(bayes_R2(model))
  cat('\n=====================================================================================\n')
  
  cat('\n\nModel Summary:')
  cat('\n=====================================================================================\n')
  print(summary(model))
  
  if(impute==0){
    cat('WARNING! If using multiply imputed datasets ignore R=hat')
    cat('\nUse the code below to obtain interpretable R=hat values for each data set')
    cat('\nround(modelname$rhats[,1:tot.pars], 3)')  
    cat('\nwhere "tot.pars" is the total number of model parameters to return (usually fixed effects)')
  }
  else{
    cat('\n\nPotential Scale Reduction Factor for Each Imputed Dataset:')
    cat('\n=====================================================================================\n')
    print(round(model$rhats[,1:tot.pars], 3))
    cat('\n=====================================================================================\n')
  }
  sink()
}

############################################################################################################
#Getting vector of IDs (required for group mean centering values)
IDs<-unique(dat$subid)

ID<-vector()
c.Worst<-vector()
c.Best<-vector()
c.DN<-vector()
NegAff<-vector()
PosAff<-vector()
BFI_O<-vector()
BFI_C<-vector()
BFI_E<-vector()
BFI_A<-vector()
Dep<-vector()
Gender<-vector()

#Bringing in covariates for imputation - 
#Want to include gender, BFI_E, BFI_C, BFI_O, BFI_A, GenDep

for(i in 1:length(IDs)){
  dat.temp<-dat[dat$subid==IDs[i],]
  mean.Worst.temp<-mean(dat.temp$WorstEvent_Neg, na.rm=T)
  mean.Best.temp<-mean(dat.temp$BestEvent_Pos, na.rm=T)
  for(j in 1:length(dat.temp[,1])){
    ID<-c(ID, IDs[i])
    c.Worst<-c(c.Worst, dat.temp$WorstEvent_Neg[j]-mean.Worst.temp)
    c.Best<-c(c.Best, dat.temp$BestEvent_Pos[j]-mean.Best.temp)
    c.DN<-c(c.DN, dat.temp$c.AP_both[j])
    NegAff<-c(NegAff, dat.temp$MAFS_NA[j])
    PosAff<-c(PosAff, dat.temp$MAFS_PA[j])
    BFI_O<-c(BFI_O, dat.temp$Battery_BFI_O[j])
    BFI_C<-c(BFI_C, dat.temp$Battery_BFI_C[j])
    BFI_E<-c(BFI_E, dat.temp$Battery_BFI_E[j])
    BFI_A<-c(BFI_A, dat.temp$Battery_BFI_A[j])
    Dep<-c(Dep, dat.temp$Battery_IDASII64_GenDep[j])
    Gender<-c(Gender, dat.temp$SexNum[j])
  }
}

dat.study1<-data.frame(ID, 
                       c.Worst, 
                       c.Best, 
                       c.DN, 
                       NegAff,
                       PosAff, 
                       BFI_O, 
                       BFI_C, 
                       BFI_E, 
                       BFI_A,
                       Dep, 
                       Gender=Gender-1) #Note that I subtracted 1 
                                        #Be sure to confirm male/female codes
                                      

##########################################################################################################
##########################################################################################################
#Creating a brms model for Study 1 - using as a discovery sample
##########################################################################################################
##########################################################################################################

#Multilevel imputation will take an incredibly long time (even when running on eleven cores)

#inspecting missing patterns in the data: 
md.pattern(dat.study1)

#Approximately 20% of the data is missing - will impute using a multilevel approach

#=========================================================================================================
#Create initial matrix for convenience
fml<- c.Worst + c.Best + NegAff + PosAff  ~ 
  1 + c.DN + Gender + (1|ID)

imp<-panImpute(dat.study1, formula=fml, n.burn=10000, 
               n.iter = 5000, m=20, seed = 0716)

#some of the imputed outcomes are below 1 and above 5 - will recode these variables later
dat<-mitmlComplete(imp)

#Inspecting imputation quality properties - especially interested in outcome measures  
dat.long<-data.frame()
for(i in 1:20){
  dat.temp<-dat[[i]]
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

png(paste0(study1.graphics, '/S1_Imputation_NegativeAffect.png'), res=300, 
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

png(paste0(study1.graphics, '/S1_Imputation_PositiveAffect.png'), res=300, 
    units='in', height=5, width = 8)
g1
dev.off()

#Plotting Deviations from Mean Worst Ratings
g1<-ggplot()+
  geom_freqpoly(data = dat.long, aes(x=c.Worst, group=IMP, color=Orig), stat = 'density')+
  theme(legend.title = element_blank())+
  ggtitle('Imputed Distributions of Worst Event (Individually Mean-Centered)')+
  ylab('Density')+
  xlab('Deviations from Mean Worst Event Ratings')
g1

png(paste0(study1.graphics, '/S1_Imputation_WorstEvent.png'), res=300, 
    units='in', height=5, width = 8)
g1
dev.off()

#Plotting Deviations from Mean Best Ratings
g1<-ggplot()+
  geom_freqpoly(data = dat.long, aes(x=c.Best, group=IMP, color=Orig), stat = 'density')+
  theme(legend.title = element_blank())+
  ggtitle('Imputed Distributions of Best Event (Individually Mean-Centered)')+
  ylab('Density')+
  xlab('Deviations from Mean Best Event Ratings')
g1

png(paste0(study1.graphics, '/S1_Imputation_BestEvent.png'), res=300, 
    units='in', height=5, width = 8)
g1
dev.off()

#summarizing imputations - checking for convergence
sink(paste0(study1.out, '/Imputation_Summary.txt'))
summary(imp)
psych::describeBy(dat.long, group='IMP')
sink()    #warnings have to do with factor variables for the most part


#Looping through and constraining values to be 1 <= x <= 5.  
#Note that 20 is the number of imputed data sets:
for(i in 1:20){
  dat[[i]]$T.NegAff<-ifelse(dat[[i]]$NegAff>5, 5, dat[[i]]$NegAff)
  dat[[i]]$T.NegAff<-ifelse(dat[[i]]$T.NegAff<1, 1, dat[[i]]$T.NegAff)
  dat[[i]]$T.PosAff<-ifelse(dat[[i]]$PosAff>5, 5, dat[[i]]$PosAff)
  dat[[i]]$T.PosAff<-ifelse(dat[[i]]$T.PosAff<1, 1, dat[[i]]$T.PosAff)
}

#Quirky formatting problems - getting around it this way... 
#   The issue is that brm_multiple requires a list of datasets
#   It cannot recognize an mitml.list object though
#   So I made a list out of list
dat.study1.list<-list()

for(i in 1:20){
  dat.study1.list[[i]]<-dat[[i]]
}
#---------------------------------------------------------------------------------------------------------
#NEGATIVE AFFECT MODEL
#---------------------------------------------------------------------------------------------------------

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#Negative unconditional model (i.e., random intercepts only)
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Neg_ucm<-brms::brm_multiple(T.NegAff~1+(1|ID), 
                            data = dat.study1.list, 
                            family = 'lognormal',
                            prior = c(set_prior('student_t(3,0,1)', class='sd'), 
                                      set_prior('normal(0,2)', class='Intercept')),
                            warmup = 3000, 
                            iter = 6000,
                            thin = 3, 
                            chains = 4,
                            control = list(adapt_delta=.90, 
                                           max_treedepth=10), 
                            save_model = paste0(stan.code, '/S1_ucm.stan'))

#Examine Convergence in each imputed data set
round(Neg_ucm$rhats[,1:3], 3)   #Can get random effects as well 
summary(Neg_ucm)                #note that overall R-hat "looks" bad
#that is due to the fact that each model is working with different data
#need to look at the R-hats for the submodels (see code above)

#plot convergence
png(paste0(study1.graphics, '/S1_NegAff_ucm_traceplots.png'), 
    units='in', width = 8, height=8, res=300)
plot(Neg_ucm, pars = c("b_Intercept", 'sigma', 'sd'))
dev.off()

#plot recovery of original distribution
png(paste0(study1.graphics, '/S1_NegAff_ucm_pp_dens.png'), 
    units='in', width = 8, height=8, res=300)
pp_check(Neg_ucm, nsamples=100)
dev.off()

#Storing most important aspects of Summary
bayes.to.txt(model = Neg_ucm, 
             out.folder = study1.model,
             file = paste0('S1_NegAff_ucm_Summary_', Sys.Date()),
             DF = dat.study1.list[[1]], 
             tot.pars = 3, 
             impute = TRUE)

#----------------------------------------------------------------------------------------
Neg_lv1<-brms::brm_multiple(T.NegAff~1+
                              c.Worst+c.Best+
                              (1+c.Worst+c.Best|ID), 
                            data = dat.study1.list, 
                            family = 'lognormal',
                            prior = c(set_prior('student_t(3,0,1)', class='sd'), 
                                      set_prior('normal(0,2)', class='Intercept'), 
                                      set_prior('normal(0,2)', class='b'), 
                                      set_prior('lkj(2)', class='cor')),
                            warmup = 3000, 
                            iter = 6000,
                            thin = 3, 
                            chains = 4,
                            control = list(adapt_delta=.90, 
                                           max_treedepth=10), 
                            save_model = paste0(stan.code, '/S1_lv1.stan'))

#Examine Convergence in each imputed data set
round(Neg_lv1$rhats[,1:7], 3)   #Can get random effects as well 
summary(Neg_lv1)                #note that overall R-hat "looks" bad
#that is due to the fact that each model is working with different data
#need to look at the R-hats for the submodels (see code above)

#plot convergence
png(paste0(study1.graphics, '/S1_NegAff_lv1_traceplots.png'), 
    units='in', width = 8, height=8, res=300)
plot(Neg_lv1, pars = c("^b_", 'sigma', '^sd_'), N=7)
dev.off()

#plot recovery of original distribution
png(paste0(study1.graphics, '/S1_NegAff_lv1_pp_dens.png'), 
    units='in', width = 8, height=8, res=300)
pp_check(Neg_lv1, nsamples=100)
dev.off()

#Storing most important aspects of Summary
bayes.to.txt(model = Neg_lv1, 
             out.folder = study1.model,
             file = paste0('S1_NegAff_lv1_Summary_', Sys.Date()),
             DF = dat.study1.list[[1]], 
             tot.pars = 7, 
             impute = TRUE)

#----------------------------------------------------------------------------------------
Neg_lv2<-brms::brm_multiple(T.NegAff~1+c.DN+
                              (1|ID), 
                            data = dat.study1.list, 
                            family = 'lognormal',
                            prior = c(set_prior('student_t(3,0,1)', class='sd'), 
                                      set_prior('normal(0,2)', class='Intercept'), 
                                      set_prior('normal(0,2)', class='b')),
                            warmup = 3000, 
                            iter = 6000,
                            thin = 3, 
                            chains = 4,
                            control = list(adapt_delta=.90, 
                                           max_treedepth=10), 
                            save_model = paste0(stan.code, '/S1_lv2.stan'))

#Examine Convergence in each imputed data set
round(Neg_lv2$rhats[,1:4], 3)   #Can get random effects as well 
summary(Neg_lv2)                #note that overall R-hat "looks" bad
#that is due to the fact that each model is working with different data
#need to look at the R-hats for the submodels (see code above)

#plot convergence
png(paste0(study1.graphics, '/S1_NegAff_lv2_traceplots.png'), 
    units='in', width = 8, height=8, res=300)
plot(Neg_lv2, pars = c("^b_", 'sigma', '^sd_'), N=4)
dev.off()

#plot recovery of original distribution
png(paste0(study1.graphics, '/S1_NegAff_lv2_pp_dens.png'), 
    units='in', width = 8, height=8, res=300)
pp_check(Neg_lv2, nsamples=100)
dev.off()

#Storing most important aspects of Summary
bayes.to.txt(model = Neg_lv2, 
             out.folder = study1.model,
             file = paste0('S1_NegAff_lv2_Summary_', Sys.Date()),
             DF = dat.study1.list[[1]], 
             tot.pars = 4, 
             impute = TRUE)

#----------------------------------------------------------------------------------------
Neg_cross<-brms::brm_multiple(T.NegAff~1+c.DN*c.Worst+c.DN*c.Best+
                              (1+c.Worst+c.Best|ID), 
                            data = dat.study1.list, 
                            family = 'lognormal',
                            prior = c(set_prior('student_t(3,0,1)', class='sd'), 
                                      set_prior('normal(0,2)', class='Intercept'), 
                                      set_prior('normal(0,2)', class='b')),
                            warmup = 3000, 
                            iter = 6000,
                            thin = 3, 
                            chains = 4,
                            control = list(adapt_delta=.90, 
                                           max_treedepth=10), 
                            save_model = paste0(stan.code, '/S1_cross.stan'))

#Examine Convergence in each imputed data set
round(Neg_cross$rhats[,1:13], 3)   #Can get random effects as well 
summary(Neg_cross)                #note that overall R-hat "looks" bad
#that is due to the fact that each model is working with different data
#need to look at the R-hats for the submodels (see code above)

#plot convergence
png(paste0(study1.graphics, '/S1_NegAff_cross_traceplots.png'), 
    units='in', width = 8, height=14, res=300)
plot(Neg_cross, pars = c("^b_", 'sigma', '^sd_'), N=13)
dev.off()

#plot recovery of original distribution
png(paste0(study1.graphics, '/S1_NegAff_cross_pp_dens.png'), 
    units='in', width = 8, height=8, res=300)
pp_check(Neg_cross, nsamples=100)
dev.off()

#Storing most important aspects of Summary
bayes.to.txt(model = Neg_cross, 
             out.folder = study1.model,
             file = paste0('S1_NegAff_cross_Summary_', Sys.Date()),
             DF = dat.study1.list[[1]], 
             tot.pars = 13, 
             impute = TRUE)

#---------------------------------------------------------------------------------------------------------
#POSITIVE AFFECT MODEL
#---------------------------------------------------------------------------------------------------------
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#Positive unconditional model (i.e., random intercepts only)
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Pos_ucm<-brms::brm_multiple(T.PosAff~1+(1|ID), 
                            data = dat.study1.list, 
                            family = 'normal',
                            prior = c(set_prior('student_t(3,0,1)', class='sd'), 
                                      set_prior('normal(0,2)', class='Intercept')),
                            warmup = 3000, 
                            iter = 6000,
                            thin = 3, 
                            chains = 4,
                            control = list(adapt_delta=.90, 
                                           max_treedepth=10), 
                            save_model = paste0(stan.code, '/S1_ucm.stan'))

#Examine Convergence in each imputed data set
round(Pos_ucm$rhats[,1:3], 3)   #Can get random effects as well 
summary(Pos_ucm)                #note that overall R-hat "looks" bad
#that is due to the fact that each model is working with different data
#need to look at the R-hats for the submodels (see code above)

#plot convergence
png(paste0(study1.graphics, '/S1_PosAff_ucm_traceplots.png'), 
    units='in', width = 8, height=8, res=300)
plot(Pos_ucm, pars = c("b_Intercept", 'sigma', 'sd'))
dev.off()

#plot recovery of original distribution
png(paste0(study1.graphics, '/S1_PosAff_ucm_pp_dens.png'), 
    units='in', width = 8, height=8, res=300)
pp_check(Pos_ucm, nsamples=100)
dev.off()

#Storing most important aspects of Summary
bayes.to.txt(model = Pos_ucm, 
             out.folder = study1.model,
             file = paste0('S1_PosAff_ucm_Summary_', Sys.Date()),
             DF = dat.study1.list[[1]], 
             tot.pars = 3, 
             impute = TRUE)

#----------------------------------------------------------------------------------------
Pos_lv1<-brms::brm_multiple(T.PosAff~1+
                              c.Worst+c.Best+
                              (1+c.Worst+c.Best|ID), 
                            data = dat.study1.list, 
                            family = 'normal',
                            prior = c(set_prior('student_t(3,0,1)', class='sd'), 
                                      set_prior('normal(0,2)', class='Intercept'), 
                                      set_prior('normal(0,2)', class='b'), 
                                      set_prior('lkj(2)', class='cor')),
                            warmup = 3000, 
                            iter = 6000,
                            thin = 3, 
                            chains = 4,
                            control = list(adapt_delta=.90, 
                                           max_treedepth=10), 
                            save_model = paste0(stan.code, '/S1_lv1.stan'))

#Examine Convergence in each imputed data set
round(Pos_lv1$rhats[,1:7], 3)   #Can get random effects as well 
summary(Pos_lv1)                #note that overall R-hat "looks" bad
#that is due to the fact that each model is working with different data
#need to look at the R-hats for the submodels (see code above)

#plot convergence
png(paste0(study1.graphics, '/S1_PosAff_lv1_traceplots.png'), 
    units='in', width = 8, height=8, res=300)
plot(Pos_lv1, pars = c("^b_", 'sigma', '^sd_'), N=7)
dev.off()

#plot recovery of original distribution
png(paste0(study1.graphics, '/S1_PosAff_lv1_pp_dens.png'), 
    units='in', width = 8, height=8, res=300)
pp_check(Pos_lv1, nsamples=100)
dev.off()

#Storing most important aspects of Summary
bayes.to.txt(model = Pos_lv1, 
             out.folder = study1.model,
             file = paste0('S1_PosAff_lv1_Summary_', Sys.Date()),
             DF = dat.study1.list[[1]], 
             tot.pars = 7, 
             impute = TRUE)

#----------------------------------------------------------------------------------------
Pos_lv2<-brms::brm_multiple(T.PosAff~1+c.DN+
                              (1|ID), 
                            data = dat.study1.list, 
                            family = 'normal',
                            prior = c(set_prior('student_t(3,0,1)', class='sd'), 
                                      set_prior('normal(0,2)', class='Intercept'), 
                                      set_prior('normal(0,2)', class='b')),
                            warmup = 3000, 
                            iter = 6000,
                            thin = 3, 
                            chains = 4,
                            control = list(adapt_delta=.90, 
                                           max_treedepth=10), 
                            save_model = paste0(stan.code, '/S1_lv2.stan'))

#Examine Convergence in each imputed data set
round(Pos_lv2$rhats[,1:4], 3)   #Can get random effects as well 
summary(Pos_lv2)                #note that overall R-hat "looks" bad
#that is due to the fact that each model is working with different data
#need to look at the R-hats for the submodels (see code above)

#plot convergence
png(paste0(study1.graphics, '/S1_PosAff_lv2_traceplots.png'), 
    units='in', width = 8, height=8, res=300)
plot(Pos_lv2, pars = c("^b_", 'sigma', '^sd_'), N=4)
dev.off()

#plot recovery of original distribution
png(paste0(study1.graphics, '/S1_PosAff_lv2_pp_dens.png'), 
    units='in', width = 8, height=8, res=300)
pp_check(Pos_lv2, nsamples=100)
dev.off()

#Storing most important aspects of Summary
bayes.to.txt(model = Pos_lv2, 
             out.folder = study1.model,
             file = paste0('S1_PosAff_lv2_Summary_', Sys.Date()),
             DF = dat.study1.list[[1]], 
             tot.pars = 4, 
             impute = TRUE)

#----------------------------------------------------------------------------------------
Pos_cross<-brms::brm_multiple(T.PosAff~1+c.DN*c.Worst+c.DN*c.Best+
                                (1+c.Worst+c.Best|ID), 
                              data = dat.study1.list, 
                              family = 'normal',
                              prior = c(set_prior('student_t(3,0,1)', class='sd'), 
                                        set_prior('normal(0,2)', class='Intercept'), 
                                        set_prior('normal(0,2)', class='b')),
                              warmup = 3000, 
                              iter = 6000,
                              thin = 3, 
                              chains = 4,
                              control = list(adapt_delta=.90, 
                                             max_treedepth=10), 
                              save_model = paste0(stan.code, '/S1_cross.stan'))

#Examine Convergence in each imputed data set
round(Pos_cross$rhats[,1:13], 3)   #Can get random effects as well 
summary(Pos_cross)                #note that overall R-hat "looks" bad
#that is due to the fact that each model is working with different data
#need to look at the R-hats for the submodels (see code above)

#plot convergence
png(paste0(study1.graphics, '/S1_PosAff_cross_traceplots.png'), 
    units='in', width = 8, height=14, res=300)
plot(Pos_cross, pars = c("^b_", 'sigma', '^sd_'), N=13)
dev.off()

#plot recovery of original distribution
png(paste0(study1.graphics, '/S1_PosAff_cross_pp_dens.png'), 
    units='in', width = 8, height=8, res=300)
pp_check(Pos_cross, nsamples=100)
dev.off()

#Storing most important aspects of Summary
bayes.to.txt(model = Pos_cross, 
             out.folder = study1.model,
             file = paste0('S1_PosAff_cross_Summary_', Sys.Date()),
             DF = dat.study1.list[[1]], 
             tot.pars = 13, 
             impute = TRUE)