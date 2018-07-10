############################################################################################################
#Authors:     Barstead, M.G., DeYoung, K.D., Anderson, A. S., & Shackman, A. J.

#Title:       The moment-to-moment affective experience of dispositionally negative individuals 

#Contact:     barstead@umd.edu

#Contents:    Bayesian Regression Analyses - Study 2
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
study2.out<-paste0(wd, '/Study2')
study2.graphics<-paste0(study2.out, '/Graphics')
study2.model<-paste0(study2.out, '/Model_Summaries')
stan.code<-paste0(wd, '/Stan_code')

dat.lv1<-read.csv(paste0(data.folder, '/MB_PAX_T1_Baseline_EMA_FINAL MASTER_070318.csv'), 
                  header = TRUE, stringsAsFactors = FALSE)

dat.lv1$ID.int<-as.numeric(substr(dat.lv1$ID, start = 4, stop=6))
dat.lv1<-dat.lv1[dat.lv1$STUDY.YEAR=='Y1',]   #selecting only year 1

dat.lv1.use<-dat.lv1[dat.lv1$AvailableForAnalyses==1,]

joy.cols<-5:7
calm.cols<-8:10
anx.cols<-11:13
ang.cols<-14:15
tired.cols<-c(16,18)
dep.cols<-c(17,19)

dat.lv1.trunc<-data.frame(ID=dat.lv1.use$ID, 
                          JOY=rowMeans(dat.lv1.use[,joy.cols])+1,
                          CALM=rowMeans(dat.lv1.use[,calm.cols])+1,
                          ANX=rowMeans(dat.lv1.use[,anx.cols])+1,
                          ANG=rowMeans(dat.lv1.use[,ang.cols])+1,
                          TIRED=rowMeans(dat.lv1.use[,tired.cols])+1,
                          DEP=rowMeans(dat.lv1.use[,dep.cols])+1, 
                          POS=rowMeans(dat.lv1.use[,c(joy.cols, calm.cols)])+1, 
                          NEG=rowMeans(dat.lv1.use[,c(anx.cols, dep.cols)])+1,
                          PosEvnt = dat.lv1.use$EMA_12_PosEvent, 
                          NegEvnt = dat.lv1.use$EMA_13_NegEvent, 
                          PosEvnt.Amnt=dat.lv1.use$EMA_12a_PosAmt, 
                          NegEvnt.Amnt=dat.lv1.use$EMA_13a_NegAmt,
                          PosEvnt.Imp=dat.lv1.use$EMA_12b_PosImport, 
                          NegEvnt.Imp=dat.lv1.use$EMA_13b_NegImport,
                          DAY_Corrected = dat.lv1.use$DAY_Corrected, 
                          SIG_Corrected = dat.lv1.use$SIG_Corrected,
                          stringsAsFactors = FALSE)

psych::pairs.panels(dat.lv1.trunc[,2:9])

dat.lv1.trunc$MB_id<-paste0(dat.lv1.trunc$ID, dat.lv1.trunc$DAY_Corrected, dat.lv1.trunc$SIG_Corrected)

dat.lv1.temp<-data.frame(ID = rep(unique(dat.lv1.trunc$ID), each = 56), 
                         DAY = rep(rep(1:7, each=8), times=length(unique(dat.lv1.trunc$ID))),
                         SIG = rep(101:108), 
                         stringsAsFactors = F)

dat.lv1.temp$MB_id<-paste0(dat.lv1.temp$ID, dat.lv1.temp$DAY, dat.lv1.temp$SIG)

dat.lv1.trunc<-merge(dat.lv1.trunc, dat.lv1.temp, by='MB_id', all=TRUE)
#A little bit of cleanup - trimming down unncessary columns made in process
dat.lv1.trunc<-dat.lv1.trunc[,-17:-19]
colnames(dat.lv1.trunc)[2]<-'ID'

dat.lv1.trunc$ID<-ifelse(!is.na(dat.lv1.trunc$ID), dat.lv1.trunc$ID, substr(dat.lv1.trunc$MB_id, start=1, stop=6))
dat.lv2<-read.csv(paste0(data.folder, '/PAX_yr1_Level2_dataset.csv'), header=TRUE, 
                  stringsAsFactors = FALSE)

colnames(dat.lv2)[1]<-'ID'
dat.lv2.trunc<-dat.lv2[,c(1,3,72)]

dat.study2<-merge(dat.lv1.trunc, dat.lv2.trunc, by='ID')
table(dat.study2$ID)  #Everyone should be at 56 observations now - have added missingness... 

############################################################################################################
#Creating a Function to Finalize Imputed Model Output: 
bayes.to.txt<-function(model=NULL, 
                       out.folder=NULL, 
                       file=NULL, 
                       DF=NULL, 
                       tot.pars=3, 
                       impute=TRUE){
  require(parallel)
  require(benchmarkme)
  sink(paste0(out.folder, '/', file, '.txt'))
  cat('System Information:')
  cat('\n=====================================================================================')
  cat(paste0('\nProcessor:', '\t\t\t', benchmarkme::get_cpu()$model_name))
  cat(paste0('\nNumber of Threads:', '\t\t', detectCores(logical=T)))
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
#inspecting missing patterns in the data: 
md.pattern(dat.study2)  

#Will only be imputing for the Potential Dependent variables
#Create initial matrix for convenience
fml<- JOY + CALM + ANX + ANG + TIRED + DEP + POS + NEG + PosEvnt + NegEvnt ~ 
  1 + DN_comb + Sex + (1|ID)

imp<-panImpute(dat.study2, formula=fml, n.burn=10000, 
               n.iter = 5000, m=20, seed = 0320)

dat<-mitmlComplete(imp, force.list = TRUE)

#Inspecting imputation properties - especially interested in outcome measures  
dat.long<-data.frame()
for(i in 1:20){
  dat.temp<-dat[[i]]
  dat.temp$IMP<-rep(i, length(dat.temp[,1]))
  dat.long<-rbind(dat.long, dat.temp)
}

#Creating variables for plotting imputed vs. original distributions:
dat.study2$IMP<-rep(0, length(dat.study2[,1]))
dat.long<-rbind(dat.study2, dat.long)

dat.long$Orig<-ifelse(dat.long$IMP<1, 'Original', 'Imputed')

#Plotting Joyous Affect Imputation
g1<-ggplot()+
  geom_freqpoly(data = dat.long, aes(x=JOY, group=IMP, color=Orig), stat = 'density')+
  theme(legend.title = element_blank())+
  ggtitle('Imputed Distributions of Momentary Joyous Affect')+
  ylab('Density')+
  xlab('Momentary Joyous Affect')
g1

png(paste0(study2.graphics, '/S2_Imputation_JoyousAffect.png'), res=300, 
    units='in', height=5, width = 8)
g1
dev.off()

#Plotting Calm Affect Imputation
g1<-ggplot()+
  geom_freqpoly(data = dat.long, aes(x=CALM, group=IMP, color=Orig), stat = 'density')+
  theme(legend.title = element_blank())+
  ggtitle('Imputed Distributions of Momentary Calm Affect')+
  ylab('Density')+
  xlab('Momentary Calm Affect')
g1

png(paste0(study2.graphics, '/S2_Imputation_CalmAffect.png'), res=300, 
    units='in', height=5, width = 8)
g1
dev.off()


#Plotting Anxious Affect Imputation
g1<-ggplot()+
  geom_freqpoly(data = dat.long, aes(x=ANX, group=IMP, color=Orig), stat = 'density')+
  theme(legend.title = element_blank())+
  ggtitle('Imputed Distributions of Momentary Anxious Affect')+
  ylab('Density')+
  xlab('Momentary Anxious Affect')
g1

png(paste0(study2.graphics, '/S2_Imputation_AnxiousAffect.png'), res=300, 
    units='in', height=5, width = 8)
g1
dev.off()

#Plotting Angry Affect Imputation
g1<-ggplot()+
  geom_freqpoly(data = dat.long, aes(x=ANG, group=IMP, color=Orig), stat = 'density')+
  theme(legend.title = element_blank())+
  ggtitle('Imputed Distributions of Momentary Angry Affect')+
  ylab('Density')+
  xlab('Momentary Angry Affect')
g1

png(paste0(study2.graphics, '/S2_Imputation_AngryAffect.png'), res=300, 
    units='in', height=5, width = 8)
g1
dev.off()

#Plotting Tired Affect Imputation
g1<-ggplot()+
  geom_freqpoly(data = dat.long, aes(x=TIRED, group=IMP, color=Orig), stat = 'density')+
  theme(legend.title = element_blank())+
  ggtitle('Imputed Distributions of Momentary Tired Affect')+
  ylab('Density')+
  xlab('Momentary Tired Affect')
g1

png(paste0(study2.graphics, '/S2_Imputation_TiredAffect.png'), res=300, 
    units='in', height=5, width = 8)
g1
dev.off()

#Plotting Depressed Affect Imputation
g1<-ggplot()+
  geom_freqpoly(data = dat.long, aes(x=DEP, group=IMP, color=Orig), stat = 'density')+
  theme(legend.title = element_blank())+
  ggtitle('Imputed Distributions of Momentary Depressed Affect')+
  ylab('Density')+
  xlab('Momentary Depressed Affect')
g1

png(paste0(study2.graphics, '/S2_Imputation_DepressedAffect.png'), res=300, 
    units='in', height=5, width = 8)
g1
dev.off()

#Plotting Positive Affect Imputation
g1<-ggplot()+
  geom_freqpoly(data = dat.long, aes(x=POS, group=IMP, color=Orig), stat = 'density')+
  theme(legend.title = element_blank())+
  ggtitle('Imputed Distributions of Momentary Positive Affect')+
  ylab('Density')+
  xlab('Momentary Positive Affect')
g1

png(paste0(study2.graphics, '/S2_Imputation_PositiveAffect.png'), res=300, 
    units='in', height=5, width = 8)
g1
dev.off()

#Plotting Negative Affect Imputation
g1<-ggplot()+
  geom_freqpoly(data = dat.long, aes(x=NEG, group=IMP, color=Orig), stat = 'density')+
  theme(legend.title = element_blank())+
  ggtitle('Imputed Distributions of Momentary Negative Affect')+
  ylab('Density')+
  xlab('Momentary Negative Affect')
g1

png(paste0(study2.graphics, '/S2_Imputation_NegativeAffect.png'), res=300, 
    units='in', height=5, width = 8)
g1
dev.off()

############################################################################################################
#summarizing imputations - checking for convergence
sink(paste0(study2.out, '/S2_Imputation_Summary.txt'))
summary(imp)
psych::describeBy(dat.long, group='IMP')
sink()    #warnings have to do with factor variables for the most part

dat.study2.list<-list()

for(i in 1:20){
  dat.study2.list[[i]]<-dat[[i]]
}

#Looping through and constraining values to be 1 <= x <= 5.  
#Note that 20 is the number of imputed data sets:
for(i in 1:20){
  dat.study2.list[[i]]$T.JOY<-ifelse(dat.study2.list[[i]]$JOY>5, 5, dat.study2.list[[i]]$JOY)
  dat.study2.list[[i]]$T.JOY<-ifelse(dat.study2.list[[i]]$T.JOY<1, 1, dat.study2.list[[i]]$T.JOY)
  dat.study2.list[[i]]$T.CALM<-ifelse(dat.study2.list[[i]]$CALM>5, 5, dat.study2.list[[i]]$CALM)
  dat.study2.list[[i]]$T.CALM<-ifelse(dat.study2.list[[i]]$T.CALM<1, 1, dat.study2.list[[i]]$T.CALM)
  dat.study2.list[[i]]$T.ANX<-ifelse(dat.study2.list[[i]]$ANX>5, 5, dat.study2.list[[i]]$ANX)
  dat.study2.list[[i]]$T.ANX<-ifelse(dat.study2.list[[i]]$T.ANX<1, 1, dat.study2.list[[i]]$T.ANX)
  dat.study2.list[[i]]$T.ANG<-ifelse(dat.study2.list[[i]]$ANG>5, 5, dat.study2.list[[i]]$ANG)
  dat.study2.list[[i]]$T.ANG<-ifelse(dat.study2.list[[i]]$T.ANG<1, 1, dat.study2.list[[i]]$T.ANG)
  dat.study2.list[[i]]$T.TIRED<-ifelse(dat.study2.list[[i]]$TIRED>5, 5, dat.study2.list[[i]]$TIRED)
  dat.study2.list[[i]]$T.TIRED<-ifelse(dat.study2.list[[i]]$T.TIRED<1, 1, dat.study2.list[[i]]$T.TIRED)
  dat.study2.list[[i]]$T.DEP<-ifelse(dat.study2.list[[i]]$DEP>5, 5, dat.study2.list[[i]]$DEP)
  dat.study2.list[[i]]$T.DEP<-ifelse(dat.study2.list[[i]]$T.DEP<1, 1, dat.study2.list[[i]]$T.DEP)
  dat.study2.list[[i]]$T.POS<-ifelse(dat.study2.list[[i]]$POS>5, 5, dat.study2.list[[i]]$POS)
  dat.study2.list[[i]]$T.POS<-ifelse(dat.study2.list[[i]]$T.POS<1, 1, dat.study2.list[[i]]$T.POS)
  dat.study2.list[[i]]$T.NEG<-ifelse(dat.study2.list[[i]]$NEG>5, 5, dat.study2.list[[i]]$NEG)
  dat.study2.list[[i]]$T.NEG<-ifelse(dat.study2.list[[i]]$T.NEG<1, 1, dat.study2.list[[i]]$T.NEG)
  dat.study2.list[[i]]$T.PosEvnt<-ifelse(dat.study2.list[[i]]$PosEvnt<=.5, 0, 1)
  dat.study2.list[[i]]$T.NegEvnt<-ifelse(dat.study2.list[[i]]$NegEvnt<=.5, 0, 1)
}

#---------------------------------------------------------------------------------------------------------
#NEGATIVE EVENT MODEL
#---------------------------------------------------------------------------------------------------------

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#Negative Event unconditional model (i.e., random intercepts only)
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
NegEvnt_ucm<-brms::brm_multiple(T.NegEvnt~1+(1|ID), 
                            data = dat.study2.list, 
                            family = 'bernoulli',
                            prior = c(set_prior('student_t(3,0,1)', class='sd'), 
                                      set_prior('normal(0,2)', class='Intercept')),
                            warmup = 3000, 
                            iter = 6000,
                            thin = 3, 
                            chains = 4,
                            control = list(adapt_delta=.90, 
                                           max_treedepth=10), 
                            save_model = paste0(stan.code, '/S2_Evnt_ucm.stan'))

#Examine Convergence in each imputed data set
round(NegEvnt_ucm$rhats[,1:3], 3)   #Can get random effects as well 
summary(NegEvnt_ucm)                #note that overall R-hat "looks" bad
#that is due to the fact that each model is working with different data
#need to look at the R-hats for the submodels (see code above)

#plot convergence
png(paste0(study2.graphics, '/S2_NegEvnt_ucm_traceplots.png'), 
    units='in', width = 8, height=8, res=300)
plot(NegEvnt_ucm, pars = c("b_Intercept", 'sigma', 'sd'))
dev.off()

#plot recovery of original distribution
png(paste0(study2.graphics, '/S2_NegEvnt_ucm_pp_dens.png'), 
    units='in', width = 8, height=8, res=300)
pp_check(NegEvnt_ucm, nsamples=100)
dev.off()

#Storing most important aspects of Summary
bayes.to.txt(model = NegEvnt_ucm, 
             out.folder = study2.model,
             file = paste0('S2_NegEvnt_ucm_Summary_', Sys.Date()),
             DF = dat.study2.list[[1]], 
             tot.pars = 3, 
             impute = TRUE)

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#Negative Event level 2 model 
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
NegEvnt_lv2<-brms::brm_multiple(T.NegEvnt~1+DN_comb+(1|ID), 
                                data = dat.study2.list, 
                                family = 'bernoulli',
                                prior = c(set_prior('student_t(3,0,1)', class='sd'), 
                                          set_prior('normal(0,2)', class='Intercept')),
                                warmup = 3000, 
                                iter = 6000,
                                thin = 3, 
                                chains = 4,
                                control = list(adapt_delta=.90, 
                                               max_treedepth=10), 
                                save_model = paste0(stan.code, '/S2_Evnt_lv2.stan'))


#Examine Convergence in each imputed data set
round(NegEvnt_lv2$rhats[,1:4], 3)   #Can get random effects as well 
summary(NegEvnt_lv2)                #note that overall R-hat "looks" bad
#that is due to the fact that each model is working with different data
#need to look at the R-hats for the submodels (see code above)

#plot convergence
png(paste0(study2.graphics, '/S2_NegEvnt_lv2_traceplots.png'), 
    units='in', width = 8, height=8, res=300)
plot(NegEvnt_lv2, pars = c("b_Intercept",'b_DN_comb', 'sigma', 'sd'))
dev.off()

#plot recovery of original distribution
png(paste0(study2.graphics, '/S2_NegEvnt_lv2_pp_dens.png'), 
    units='in', width = 8, height=8, res=300)
pp_check(NegEvnt_lv2, nsamples=100)
dev.off()

#Storing most important aspects of Summary
bayes.to.txt(model = NegEvnt_lv2, 
             out.folder = study2.model,
             file = paste0('S2_NegEvnt_lv2_Summary_', Sys.Date()),
             DF = dat.study2.list[[1]], 
             tot.pars = 4, 
             impute = TRUE)

#---------------------------------------------------------------------------------------------------------
#POSITIVE EVENT MODEL
#---------------------------------------------------------------------------------------------------------

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#Positive Event unconditional model (i.e., random intercepts only)
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
PosEvnt_ucm<-brms::brm_multiple(T.PosEvnt~1+(1|ID), 
                                data = dat.study2.list, 
                                family = 'bernoulli',
                                prior = c(set_prior('student_t(3,0,1)', class='sd'), 
                                          set_prior('normal(0,2)', class='Intercept')),
                                warmup = 3000, 
                                iter = 6000,
                                thin = 3, 
                                chains = 4,
                                control = list(adapt_delta=.90, 
                                               max_treedepth=10), 
                                save_model = paste0(stan.code, '/S2_Evnt_ucm.stan'))


#Examine Convergence in each imputed data set
round(PosEvnt_ucm$rhats[,1:3], 3)   #Can get random effects as well 
summary(PosEvnt_ucm)                #note that overall R-hat "looks" bad
#that is due to the fact that each model is working with different data
#need to look at the R-hats for the submodels (see code above)

#plot convergence
png(paste0(study2.graphics, '/S2_PosEvnt_ucm_traceplots.png'), 
    units='in', width = 8, height=8, res=300)
plot(PosEvnt_ucm, pars = c("b_Intercept", 'sigma', 'sd'))
dev.off()

#plot recovery of original distribution
png(paste0(study2.graphics, '/S2_PosEvnt_ucm_pp_dens.png'), 
    units='in', width = 8, height=8, res=300)
pp_check(PosEvnt_ucm, nsamples=100)
dev.off()

#Storing most important aspects of Summary
bayes.to.txt(model = PosEvnt_ucm, 
             out.folder = study2.model,
             file = paste0('S2_PosEvnt_ucm_Summary_', Sys.Date()),
             DF = dat.study2.list[[1]], 
             tot.pars = 3, 
             impute = TRUE)

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#Positive Event level 2 model 
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
PosEvnt_lv2<-brms::brm_multiple(T.PosEvnt~1+DN_comb+(1|ID), 
                                data = dat.study2.list, 
                                family = 'bernoulli',
                                prior = c(set_prior('student_t(3,0,1)', class='sd'), 
                                          set_prior('normal(0,2)', class='Intercept')),
                                warmup = 3000, 
                                iter = 6000,
                                thin = 3, 
                                chains = 4,
                                control = list(adapt_delta=.90, 
                                               max_treedepth=10), 
                                save_model = paste0(stan.code, '/S2_Evnt_lv2.stan'))


#Examine Convergence in each imputed data set
round(PosEvnt_lv2$rhats[,1:4], 3)   #Can get random effects as well 
summary(PosEvnt_lv2)                #note that overall R-hat "looks" bad
#that is due to the fact that each model is working with different data
#need to look at the R-hats for the submodels (see code above)

#plot convergence
png(paste0(study2.graphics, '/S2_PosEvnt_lv2_traceplots.png'), 
    units='in', width = 8, height=8, res=300)
plot(PosEvnt_lv2, pars = c("b_Intercept",'b_DN_comb', 'sigma', 'sd'))
dev.off()

#plot recovery of original distribution
png(paste0(study2.graphics, '/S2_PosEvnt_lv2_pp_dens.png'), 
    units='in', width = 8, height=8, res=300)
pp_check(PosEvnt_lv2, nsamples=100)
dev.off()

#Storing most important aspects of Summary
bayes.to.txt(model = PosEvnt_lv2, 
             out.folder = study2.model,
             file = paste0('S2_PosEvnt_lv2_Summary_', Sys.Date()),
             DF = dat.study2.list[[1]], 
             tot.pars = 4, 
             impute = TRUE)

#---------------------------------------------------------------------------------------------------------
#JOYOUS AFFECT MODEL
#---------------------------------------------------------------------------------------------------------

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#Joyous Affect unconditional model (i.e., random intercepts only)
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
JOY_ucm<-brms::brm_multiple(T.JOY~1+(1|ID), 
                                data = dat.study2.list, 
                                family = 'lognormal',
                                prior = c(set_prior('student_t(3,0,1)', class='sd'), 
                                          set_prior('normal(0,2)', class='Intercept')),
                                warmup = 3000, 
                                iter = 6000,
                                thin = 3, 
                                chains = 4,
                                control = list(adapt_delta=.90, 
                                               max_treedepth=10), 
                                save_model = paste0(stan.code, '/S2_Mood_ucm.stan'))


#Examine Convergence in each imputed data set
round(JOY_ucm$rhats[,1:3], 3)   #Can get random effects as well 
summary(JOY_ucm)                #note that overall R-hat "looks" bad
#that is due to the fact that each model is working with different data
#need to look at the R-hats for the submodels (see code above)

#plot convergence
png(paste0(study2.graphics, '/S2_JOY_ucm_traceplots.png'), 
    units='in', width = 8, height=8, res=300)
plot(JOY_ucm, pars = c("b_Intercept", 'sigma', 'sd'))
dev.off()

#plot recovery of original distribution
png(paste0(study2.graphics, '/S2_JOY_ucm_pp_dens.png'), 
    units='in', width = 8, height=8, res=300)
pp_check(JOY_ucm, nsamples=100)
dev.off()

#Storing most important aspects of Summary
bayes.to.txt(model = JOY_ucm, 
             out.folder = study2.model,
             file = paste0('S2_JOY_ucm_Summary_', Sys.Date()),
             DF = dat.study2.list[[1]], 
             tot.pars = 3, 
             impute = TRUE)

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#Joyous Affect level 1 model 
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
JOY_lv1<-brms::brm_multiple(T.JOY~1 + T.PosEvnt + T.NegEvnt + 
                              (1 + T.PosEvnt + T.NegEvnt|ID), 
                            data = dat.study2.list, 
                            family = 'lognormal',
                            prior = c(set_prior('student_t(3,0,1)', class='sd'), 
                                      set_prior('normal(0,2)', class='Intercept'), 
                                      set_prior('lkj(2)', class='cor')),
                            warmup = 3000, 
                            iter = 6000,
                            thin = 3, 
                            chains = 4,
                            control = list(adapt_delta=.90, 
                                           max_treedepth=10), 
                            save_model = paste0(stan.code, '/S2_Mood_lv1.stan'))


#Examine Convergence in each imputed data set
round(JOY_lv1$rhats[,1:7], 3)   #Can get random effects as well 
summary(JOY_lv1)                #note that overall R-hat "looks" bad
#that is due to the fact that each model is working with different data
#need to look at the R-hats for the submodels (see code above)

#plot convergence
png(paste0(study2.graphics, '/S2_JOY_lv1_traceplots.png'), 
    units='in', width = 8, height=8, res=300)
plot(JOY_lv1, pars = c("^b_", 'sigma', 'sd'))
dev.off()

#plot recovery of original distribution
png(paste0(study2.graphics, '/S2_JOY_lv1_pp_dens.png'), 
    units='in', width = 8, height=8, res=300)
pp_check(JOY_lv1, nsamples=100)
dev.off()

#Storing most important aspects of Summary
bayes.to.txt(model = JOY_lv1, 
             out.folder = study2.model,
             file = paste0('S2_JOY_lv1_Summary_', Sys.Date()),
             DF = dat.study2.list[[1]], 
             tot.pars = 7, 
             impute = TRUE)

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#Joyous Affect level 2 model 
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
JOY_lv2<-brms::brm_multiple(T.JOY~1+DN_comb+(1|ID), 
                            data = dat.study2.list, 
                            family = 'lognormal',
                            prior = c(set_prior('student_t(3,0,1)', class='sd'), 
                                      set_prior('normal(0,2)', class='Intercept')),
                            warmup = 3000, 
                            iter = 6000,
                            thin = 3, 
                            chains = 4,
                            control = list(adapt_delta=.90, 
                                           max_treedepth=10), 
                            save_model = paste0(stan.code, '/S2_Mood_cross.stan'))

#Examine Convergence in each imputed data set
round(JOY_lv2$rhats[,1:4], 3)   #Can get random effects as well 
summary(JOY_lv2)                #note that overall R-hat "looks" bad
#that is due to the fact that each model is working with different data
#need to look at the R-hats for the submodels (see code above)

#plot convergence
png(paste0(study2.graphics, '/S2_JOY_lv2_traceplots.png'), 
    units='in', width = 8, height=8, res=300)
plot(JOY_lv2, pars = c("b_Intercept",'b_DN_comb', 'sigma', 'sd'))
dev.off()

#plot recovery of original distribution
png(paste0(study2.graphics, '/S2_JOY_lv2_pp_dens.png'), 
    units='in', width = 8, height=8, res=300)
pp_check(JOY_lv2, nsamples=100)
dev.off()

#Storing most important aspects of Summary
bayes.to.txt(model = JOY_lv2, 
             out.folder = study2.model,
             file = paste0('S2_JOY_lv2_Summary_', Sys.Date()),
             DF = dat.study2.list[[1]], 
             tot.pars = 4, 
             impute = TRUE)

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#Joyous Affect cross-level
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
JOY_cross<-brms::brm_multiple(T.JOY~1+DN_comb*T.PosEvnt + DN_comb*T.NegEvnt + 
                                (1 + PosEvnt + NegEvnt|ID), 
                              data = dat.study2.list, 
                              family = 'lognormal',
                              prior = c(set_prior('student_t(3,0,1)', class='sd'), 
                                        set_prior('normal(0,2)', class='Intercept'), 
                                        set_prior('lkj(2)', class='cor')),
                              warmup = 3000, 
                              iter = 6000,
                              thin = 3, 
                              chains = 4,
                              control = list(adapt_delta=.90, 
                                             max_treedepth=10), 
                              save_model = paste0(stan.code, '/S2_Mood_lv2.stan'))

#Examine Convergence in each imputed data set
round(JOY_cross$rhats[,1:13], 3)   #Can get random effects as well 
summary(JOY_cross)                #note that overall R-hat "looks" bad
#that is due to the fact that each model is working with different data
#need to look at the R-hats for the submodels (see code above)

#plot convergence
png(paste0(study2.graphics, '/S2_JOY_cross_traceplots.png'), 
    units='in', width = 8, height=8, res=300)
plot(JOY_cross, pars = c("b_Intercept",'b_DN_comb', 'sigma', 'sd'))
dev.off()

#plot recovery of original distribution
png(paste0(study2.graphics, '/S2_JOY_cross_pp_dens.png'), 
    units='in', width = 8, height=8, res=300)
pp_check(JOY_cross, nsamples=100)
dev.off()

#Storing most important aspects of Summary
bayes.to.txt(model = JOY_cross, 
             out.folder = study2.model,
             file = paste0('S2_JOY_cross_Summary_', Sys.Date()),
             DF = dat.study2.list[[1]], 
             tot.pars = 13, 
             impute = TRUE)

#---------------------------------------------------------------------------------------------------------
#CALM AFFECT MODEL
#---------------------------------------------------------------------------------------------------------

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#Calm Affect unconditional model (i.e., random intercepts only)
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
CALM_ucm<-brms::brm_multiple(T.CALM~1+(1|ID), 
                            data = dat.study2.list, 
                            family = 'normal',
                            prior = c(set_prior('student_t(3,0,1)', class='sd'), 
                                      set_prior('normal(0,2)', class='Intercept')),
                            warmup = 3000, 
                            iter = 6000,
                            thin = 3, 
                            chains = 4,
                            control = list(adapt_delta=.90, 
                                           max_treedepth=10), 
                            save_model = paste0(stan.code, '/S2_Mood_ucm.stan'))


#Examine Convergence in each imputed data set
round(CALM_ucm$rhats[,1:3], 3)   #Can get random effects as well 
summary(CALM_ucm)                #note that overall R-hat "looks" bad
#that is due to the fact that each model is working with different data
#need to look at the R-hats for the submodels (see code above)

#plot convergence
png(paste0(study2.graphics, '/S2_CALM_ucm_traceplots.png'), 
    units='in', width = 8, height=8, res=300)
plot(CALM_ucm, pars = c("b_Intercept", 'sigma', 'sd'))
dev.off()

#plot recovery of original distribution
png(paste0(study2.graphics, '/S2_CALM_ucm_pp_dens.png'), 
    units='in', width = 8, height=8, res=300)
pp_check(CALM_ucm, nsamples=100)
dev.off()

#Storing most important aspects of Summary
bayes.to.txt(model = CALM_ucm, 
             out.folder = study2.model,
             file = paste0('S2_CALM_ucm_Summary_', Sys.Date()),
             DF = dat.study2.list[[1]], 
             tot.pars = 3, 
             impute = TRUE)

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#Calm Affect level 1 model 
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
CALM_lv1<-brms::brm_multiple(T.CALM~1 + T.PosEvnt + T.NegEvnt + 
                              (1 + T.PosEvnt + T.NegEvnt|ID), 
                            data = dat.study2.list, 
                            family = 'normal',
                            prior = c(set_prior('student_t(3,0,1)', class='sd'), 
                                      set_prior('normal(0,2)', class='Intercept'), 
                                      set_prior('lkj(2)', class='cor')),
                            warmup = 3000, 
                            iter = 6000,
                            thin = 3, 
                            chains = 4,
                            control = list(adapt_delta=.90, 
                                           max_treedepth=10), 
                            save_model = paste0(stan.code, '/S2_Mood_lv1.stan'))


#Examine Convergence in each imputed data set
round(CALM_lv1$rhats[,1:7], 3)   #Can get random effects as well 
summary(CALM_lv1)                #note that overall R-hat "looks" bad
#that is due to the fact that each model is working with different data
#need to look at the R-hats for the submodels (see code above)

#plot convergence
png(paste0(study2.graphics, '/S2_CALM_lv1_traceplots.png'), 
    units='in', width = 8, height=8, res=300)
plot(CALM_lv1, pars = c("^b_", 'sigma', 'sd'))
dev.off()

#plot recovery of original distribution
png(paste0(study2.graphics, '/S2_CALM_lv1_pp_dens.png'), 
    units='in', width = 8, height=8, res=300)
pp_check(CALM_lv1, nsamples=100)
dev.off()

#Storing most important aspects of Summary
bayes.to.txt(model = CALM_lv1, 
             out.folder = study2.model,
             file = paste0('S2_CALM_lv1_Summary_', Sys.Date()),
             DF = dat.study2.list[[1]], 
             tot.pars = 7, 
             impute = TRUE)

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#Calm Affect level 2 model 
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
CALM_lv2<-brms::brm_multiple(T.CALM~1+DN_comb+(1|ID), 
                             data = dat.study2.list, 
                             family = 'normal',
                             prior = c(set_prior('student_t(3,0,1)', class='sd'), 
                                       set_prior('normal(0,2)', class='Intercept')),
                             warmup = 3000, 
                             iter = 6000,
                             thin = 3, 
                             chains = 4,
                             control = list(adapt_delta=.90, 
                                            max_treedepth=10), 
                             save_model = paste0(stan.code, '/S2_Mood_cross.stan'))

#Examine Convergence in each imputed data set
round(CALM_lv2$rhats[,1:4], 3)   #Can get random effects as well 
summary(CALM_lv2)                #note that overall R-hat "looks" bad
#that is due to the fact that each model is working with different data
#need to look at the R-hats for the submodels (see code above)

#plot convergence
png(paste0(study2.graphics, '/S2_CALM_lv2_traceplots.png'), 
    units='in', width = 8, height=8, res=300)
plot(CALM_lv2, pars = c("b_Intercept",'b_DN_comb', 'sigma', 'sd'))
dev.off()

#plot recovery of original distribution
png(paste0(study2.graphics, '/S2_CALM_lv2_pp_dens.png'), 
    units='in', width = 8, height=8, res=300)
pp_check(CALM_lv2, nsamples=100)
dev.off()

#Storing most important aspects of Summary
bayes.to.txt(model = CALM_lv2, 
             out.folder = study2.model,
             file = paste0('S2_CALM_lv2_Summary_', Sys.Date()),
             DF = dat.study2.list[[1]], 
             tot.pars = 4, 
             impute = TRUE)

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#Calm Affect cross-level
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
CALM_cross<-brms::brm_multiple(T.CALM~1+DN_comb*T.PosEvnt + DN_comb*T.NegEvnt + 
                                 (1 + PosEvnt + NegEvnt|ID), 
                               data = dat.study2.list, 
                               family = 'normal',
                               prior = c(set_prior('student_t(3,0,1)', class='sd'), 
                                         set_prior('normal(0,2)', class='Intercept'), 
                                         set_prior('lkj(2)', class='cor')),
                               warmup = 3000, 
                               iter = 6000,
                               thin = 3, 
                               chains = 4,
                               control = list(adapt_delta=.90, 
                                              max_treedepth=10), 
                               save_model = paste0(stan.code, '/S2_Mood_lv2.stan'))

#Examine Convergence in each imputed data set
round(CALM_cross$rhats[,1:13], 3)   #Can get random effects as well 
summary(CALM_cross)                #note that overall R-hat "looks" bad
#that is due to the fact that each model is working with different data
#need to look at the R-hats for the submodels (see code above)

#plot convergence
png(paste0(study2.graphics, '/S2_CALM_cross_traceplots.png'), 
    units='in', width = 8, height=8, res=300)
plot(CALM_cross, pars = c("b_Intercept",'b_DN_comb', 'sigma', 'sd'))
dev.off()

#plot recovery of original distribution
png(paste0(study2.graphics, '/S2_CALM_cross_pp_dens.png'), 
    units='in', width = 8, height=8, res=300)
pp_check(CALM_cross, nsamples=100)
dev.off()

#Storing most important aspects of Summary
bayes.to.txt(model = CALM_cross, 
             out.folder = study2.model,
             file = paste0('S2_CALM_cross_Summary_', Sys.Date()),
             DF = dat.study2.list[[1]], 
             tot.pars = 13, 
             impute = TRUE)

#---------------------------------------------------------------------------------------------------------
#ANXIOUS AFFECT MODEL
#---------------------------------------------------------------------------------------------------------

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#Anxious Affect unconditional model (i.e., random intercepts only)
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
ANX_ucm<-brms::brm_multiple(T.ANX~1+(1|ID), 
                            data = dat.study2.list, 
                            family = 'lognormal',
                            prior = c(set_prior('student_t(3,0,1)', class='sd'), 
                                      set_prior('normal(0,2)', class='Intercept')),
                            warmup = 3000, 
                            iter = 6000,
                            thin = 3, 
                            chains = 4,
                            control = list(adapt_delta=.90, 
                                           max_treedepth=10), 
                            save_model = paste0(stan.code, '/S2_Mood_ucm.stan'))


#Examine Convergence in each imputed data set
round(ANX_ucm$rhats[,1:3], 3)   #Can get random effects as well 
summary(ANX_ucm)                #note that overall R-hat "looks" bad
#that is due to the fact that each model is working with different data
#need to look at the R-hats for the submodels (see code above)

#plot convergence
png(paste0(study2.graphics, '/S2_ANX_ucm_traceplots.png'), 
    units='in', width = 8, height=8, res=300)
plot(ANX_ucm, pars = c("b_Intercept", 'sigma', 'sd'))
dev.off()

#plot recovery of original distribution
png(paste0(study2.graphics, '/S2_ANX_ucm_pp_dens.png'), 
    units='in', width = 8, height=8, res=300)
pp_check(ANX_ucm, nsamples=100)
dev.off()

#Storing most important aspects of Summary
bayes.to.txt(model = ANX_ucm, 
             out.folder = study2.model,
             file = paste0('S2_ANX_ucm_Summary_', Sys.Date()),
             DF = dat.study2.list[[1]], 
             tot.pars = 3, 
             impute = TRUE)

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#Anxious Affect level 1 model 
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
ANX_lv1<-brms::brm_multiple(T.ANX~1 + T.PosEvnt + T.NegEvnt + 
                              (1 + T.PosEvnt + T.NegEvnt|ID), 
                            data = dat.study2.list, 
                            family = 'lognormal',
                            prior = c(set_prior('student_t(3,0,1)', class='sd'), 
                                      set_prior('normal(0,2)', class='Intercept'), 
                                      set_prior('lkj(2)', class='cor')),
                            warmup = 3000, 
                            iter = 6000,
                            thin = 3, 
                            chains = 4,
                            control = list(adapt_delta=.90, 
                                           max_treedepth=10), 
                            save_model = paste0(stan.code, '/S2_Mood_lv1.stan'))


#Examine Convergence in each imputed data set
round(ANX_lv1$rhats[,1:7], 3)   #Can get random effects as well 
summary(ANX_lv1)                #note that overall R-hat "looks" bad
#that is due to the fact that each model is working with different data
#need to look at the R-hats for the submodels (see code above)

#plot convergence
png(paste0(study2.graphics, '/S2_ANX_lv1_traceplots.png'), 
    units='in', width = 8, height=8, res=300)
plot(ANX_lv1, pars = c("^b_", 'sigma', 'sd'))
dev.off()

#plot recovery of original distribution
png(paste0(study2.graphics, '/S2_ANX_lv1_pp_dens.png'), 
    units='in', width = 8, height=8, res=300)
pp_check(ANX_lv1, nsamples=100)
dev.off()

#Storing most important aspects of Summary
bayes.to.txt(model = ANX_lv1, 
             out.folder = study2.model,
             file = paste0('S2_ANX_lv1_Summary_', Sys.Date()),
             DF = dat.study2.list[[1]], 
             tot.pars = 7, 
             impute = TRUE)

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#Anxious Affect level 2 model 
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
ANX_lv2<-brms::brm_multiple(T.ANX~1+DN_comb+(1|ID), 
                            data = dat.study2.list, 
                            family = 'lognormal',
                            prior = c(set_prior('student_t(3,0,1)', class='sd'), 
                                      set_prior('normal(0,2)', class='Intercept')),
                            warmup = 3000, 
                            iter = 6000,
                            thin = 3, 
                            chains = 4,
                            control = list(adapt_delta=.90, 
                                           max_treedepth=10), 
                            save_model = paste0(stan.code, '/S2_Mood_cross.stan'))

#Examine Convergence in each imputed data set
round(ANX_lv2$rhats[,1:4], 3)   #Can get random effects as well 
summary(ANX_lv2)                #note that overall R-hat "looks" bad
#that is due to the fact that each model is working with different data
#need to look at the R-hats for the submodels (see code above)

#plot convergence
png(paste0(study2.graphics, '/S2_ANX_lv2_traceplots.png'), 
    units='in', width = 8, height=8, res=300)
plot(ANX_lv2, pars = c("b_Intercept",'b_DN_comb', 'sigma', 'sd'))
dev.off()

#plot recovery of original distribution
png(paste0(study2.graphics, '/S2_ANX_lv2_pp_dens.png'), 
    units='in', width = 8, height=8, res=300)
pp_check(ANX_lv2, nsamples=100)
dev.off()

#Storing most important aspects of Summary
bayes.to.txt(model = ANX_lv2, 
             out.folder = study2.model,
             file = paste0('S2_ANX_lv2_Summary_', Sys.Date()),
             DF = dat.study2.list[[1]], 
             tot.pars = 4, 
             impute = TRUE)

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#Anxious Affect cross-level
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
ANX_cross<-brms::brm_multiple(T.ANX~1+DN_comb*T.PosEvnt + DN_comb*T.NegEvnt + 
                                (1 + PosEvnt + NegEvnt|ID), 
                              data = dat.study2.list, 
                              family = 'lognormal',
                              prior = c(set_prior('student_t(3,0,1)', class='sd'), 
                                        set_prior('normal(0,2)', class='Intercept'), 
                                        set_prior('lkj(2)', class='cor')),
                              warmup = 3000, 
                              iter = 6000,
                              thin = 3, 
                              chains = 4,
                              control = list(adapt_delta=.90, 
                                             max_treedepth=10), 
                              save_model = paste0(stan.code, '/S2_Mood_lv2.stan'))

#Examine Convergence in each imputed data set
round(ANX_cross$rhats[,1:13], 3)   #Can get random effects as well 
summary(ANX_cross)                #note that overall R-hat "looks" bad
#that is due to the fact that each model is working with different data
#need to look at the R-hats for the submodels (see code above)

#plot convergence
png(paste0(study2.graphics, '/S2_ANX_cross_traceplots.png'), 
    units='in', width = 8, height=8, res=300)
plot(ANX_cross, pars = c("b_Intercept",'b_DN_comb', 'sigma', 'sd'))
dev.off()

#plot recovery of original distribution
png(paste0(study2.graphics, '/S2_ANX_cross_pp_dens.png'), 
    units='in', width = 8, height=8, res=300)
pp_check(ANX_cross, nsamples=100)
dev.off()

#Storing most important aspects of Summary
bayes.to.txt(model = ANX_cross, 
             out.folder = study2.model,
             file = paste0('S2_ANX_cross_Summary_', Sys.Date()),
             DF = dat.study2.list[[1]], 
             tot.pars = 13, 
             impute = TRUE)

#---------------------------------------------------------------------------------------------------------
#ANGRY AFFECT MODEL
#---------------------------------------------------------------------------------------------------------

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#Angry Affect unconditional model (i.e., random intercepts only)
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
ANG_ucm<-brms::brm_multiple(T.ANG~1+(1|ID), 
                            data = dat.study2.list, 
                            family = 'lognormal',
                            prior = c(set_prior('student_t(3,0,1)', class='sd'), 
                                      set_prior('normal(0,2)', class='Intercept')),
                            warmup = 3000, 
                            iter = 6000,
                            thin = 3, 
                            chains = 4,
                            control = list(adapt_delta=.90, 
                                           max_treedepth=10), 
                            save_model = paste0(stan.code, '/S2_Mood_ucm.stan'))


#Examine Convergence in each imputed data set
round(ANG_ucm$rhats[,1:3], 3)   #Can get random effects as well 
summary(ANG_ucm)                #note that overall R-hat "looks" bad
#that is due to the fact that each model is working with different data
#need to look at the R-hats for the submodels (see code above)

#plot convergence
png(paste0(study2.graphics, '/S2_ANG_ucm_traceplots.png'), 
    units='in', width = 8, height=8, res=300)
plot(ANG_ucm, pars = c("b_Intercept", 'sigma', 'sd'))
dev.off()

#plot recovery of original distribution
png(paste0(study2.graphics, '/S2_ANG_ucm_pp_dens.png'), 
    units='in', width = 8, height=8, res=300)
pp_check(ANG_ucm, nsamples=100)
dev.off()

#Storing most important aspects of Summary
bayes.to.txt(model = ANG_ucm, 
             out.folder = study2.model,
             file = paste0('S2_ANG_ucm_Summary_', Sys.Date()),
             DF = dat.study2.list[[1]], 
             tot.pars = 3, 
             impute = TRUE)

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#Angry Affect level 1 model 
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
ANG_lv1<-brms::brm_multiple(T.ANG~1 + T.PosEvnt + T.NegEvnt + 
                              (1 + T.PosEvnt + T.NegEvnt|ID), 
                            data = dat.study2.list, 
                            family = 'lognormal',
                            prior = c(set_prior('student_t(3,0,1)', class='sd'), 
                                      set_prior('normal(0,2)', class='Intercept'), 
                                      set_prior('lkj(2)', class='cor')),
                            warmup = 3000, 
                            iter = 6000,
                            thin = 3, 
                            chains = 4,
                            control = list(adapt_delta=.90, 
                                           max_treedepth=10), 
                            save_model = paste0(stan.code, '/S2_Mood_lv1.stan'))


#Examine Convergence in each imputed data set
round(ANG_lv1$rhats[,1:7], 3)   #Can get random effects as well 
summary(ANG_lv1)                #note that overall R-hat "looks" bad
#that is due to the fact that each model is working with different data
#need to look at the R-hats for the submodels (see code above)

#plot convergence
png(paste0(study2.graphics, '/S2_ANG_lv1_traceplots.png'), 
    units='in', width = 8, height=8, res=300)
plot(ANG_lv1, pars = c("^b_", 'sigma', 'sd'))
dev.off()

#plot recovery of original distribution
png(paste0(study2.graphics, '/S2_ANG_lv1_pp_dens.png'), 
    units='in', width = 8, height=8, res=300)
pp_check(ANG_lv1, nsamples=100)
dev.off()

#Storing most important aspects of Summary
bayes.to.txt(model = ANG_lv1, 
             out.folder = study2.model,
             file = paste0('S2_ANG_lv1_Summary_', Sys.Date()),
             DF = dat.study2.list[[1]], 
             tot.pars = 7, 
             impute = TRUE)

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#Angry Affect level 2 model 
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
ANG_lv2<-brms::brm_multiple(T.ANG~1+DN_comb+(1|ID), 
                            data = dat.study2.list, 
                            family = 'lognormal',
                            prior = c(set_prior('student_t(3,0,1)', class='sd'), 
                                      set_prior('normal(0,2)', class='Intercept')),
                            warmup = 3000, 
                            iter = 6000,
                            thin = 3, 
                            chains = 4,
                            control = list(adapt_delta=.90, 
                                           max_treedepth=10), 
                            save_model = paste0(stan.code, '/S2_Mood_cross.stan'))

#Examine Convergence in each imputed data set
round(ANG_lv2$rhats[,1:4], 3)   #Can get random effects as well 
summary(ANG_lv2)                #note that overall R-hat "looks" bad
#that is due to the fact that each model is working with different data
#need to look at the R-hats for the submodels (see code above)

#plot convergence
png(paste0(study2.graphics, '/S2_ANG_lv2_traceplots.png'), 
    units='in', width = 8, height=8, res=300)
plot(ANG_lv2, pars = c("b_Intercept",'b_DN_comb', 'sigma', 'sd'))
dev.off()

#plot recovery of original distribution
png(paste0(study2.graphics, '/S2_ANG_lv2_pp_dens.png'), 
    units='in', width = 8, height=8, res=300)
pp_check(ANG_lv2, nsamples=100)
dev.off()

#Storing most important aspects of Summary
bayes.to.txt(model = ANG_lv2, 
             out.folder = study2.model,
             file = paste0('S2_ANG_lv2_Summary_', Sys.Date()),
             DF = dat.study2.list[[1]], 
             tot.pars = 4, 
             impute = TRUE)

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#Angry Affect cross-level
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
ANG_cross<-brms::brm_multiple(T.ANG~1+DN_comb*T.PosEvnt + DN_comb*T.NegEvnt + 
                                (1 + PosEvnt + NegEvnt|ID), 
                              data = dat.study2.list, 
                              family = 'lognormal',
                              prior = c(set_prior('student_t(3,0,1)', class='sd'), 
                                        set_prior('normal(0,2)', class='Intercept'), 
                                        set_prior('lkj(2)', class='cor')),
                              warmup = 3000, 
                              iter = 6000,
                              thin = 3, 
                              chains = 4,
                              control = list(adapt_delta=.90, 
                                             max_treedepth=10), 
                              save_model = paste0(stan.code, '/S2_Mood_lv2.stan'))

#Examine Convergence in each imputed data set
round(ANG_cross$rhats[,1:13], 3)   #Can get random effects as well 
summary(ANG_cross)                #note that overall R-hat "looks" bad
#that is due to the fact that each model is working with different data
#need to look at the R-hats for the submodels (see code above)

#plot convergence
png(paste0(study2.graphics, '/S2_ANG_cross_traceplots.png'), 
    units='in', width = 8, height=8, res=300)
plot(ANG_cross, pars = c("b_Intercept",'b_DN_comb', 'sigma', 'sd'))
dev.off()

#plot recovery of original distribution
png(paste0(study2.graphics, '/S2_ANG_cross_pp_dens.png'), 
    units='in', width = 8, height=8, res=300)
pp_check(ANG_cross, nsamples=100)
dev.off()

#Storing most important aspects of Summary
bayes.to.txt(model = ANG_cross, 
             out.folder = study2.model,
             file = paste0('S2_ANG_cross_Summary_', Sys.Date()),
             DF = dat.study2.list[[1]], 
             tot.pars = 13, 
             impute = TRUE)

#---------------------------------------------------------------------------------------------------------
#TIRED AFFECT MODEL
#---------------------------------------------------------------------------------------------------------

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#Tired Affect unconditional model (i.e., random intercepts only)
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
TIRED_ucm<-brms::brm_multiple(T.TIRED~1+(1|ID), 
                            data = dat.study2.list, 
                            family = 'lognormal',
                            prior = c(set_prior('student_t(3,0,1)', class='sd'), 
                                      set_prior('normal(0,2)', class='Intercept')),
                            warmup = 3000, 
                            iter = 6000,
                            thin = 3, 
                            chains = 4,
                            control = list(adapt_delta=.90, 
                                           max_treedepth=10), 
                            save_model = paste0(stan.code, '/S2_Mood_ucm.stan'))


#Examine Convergence in each imputed data set
round(TIRED_ucm$rhats[,1:3], 3)   #Can get random effects as well 
summary(TIRED_ucm)                #note that overall R-hat "looks" bad
#that is due to the fact that each model is working with different data
#need to look at the R-hats for the submodels (see code above)

#plot convergence
png(paste0(study2.graphics, '/S2_TIRED_ucm_traceplots.png'), 
    units='in', width = 8, height=8, res=300)
plot(TIRED_ucm, pars = c("b_Intercept", 'sigma', 'sd'))
dev.off()

#plot recovery of original distribution
png(paste0(study2.graphics, '/S2_TIRED_ucm_pp_dens.png'), 
    units='in', width = 8, height=8, res=300)
pp_check(TIRED_ucm, nsamples=100)
dev.off()

#Storing most important aspects of Summary
bayes.to.txt(model = TIRED_ucm, 
             out.folder = study2.model,
             file = paste0('S2_TIRED_ucm_Summary_', Sys.Date()),
             DF = dat.study2.list[[1]], 
             tot.pars = 3, 
             impute = TRUE)

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#Tired Affect level 1 model 
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
TIRED_lv1<-brms::brm_multiple(T.TIRED~1 + T.PosEvnt + T.NegEvnt + 
                              (1 + T.PosEvnt + T.NegEvnt|ID), 
                            data = dat.study2.list, 
                            family = 'lognormal',
                            prior = c(set_prior('student_t(3,0,1)', class='sd'), 
                                      set_prior('normal(0,2)', class='Intercept'), 
                                      set_prior('lkj(2)', class='cor')),
                            warmup = 3000, 
                            iter = 6000,
                            thin = 3, 
                            chains = 4,
                            control = list(adapt_delta=.90, 
                                           max_treedepth=10), 
                            save_model = paste0(stan.code, '/S2_Mood_lv1.stan'))


#Examine Convergence in each imputed data set
round(TIRED_lv1$rhats[,1:7], 3)   #Can get random effects as well 
summary(TIRED_lv1)                #note that overall R-hat "looks" bad
#that is due to the fact that each model is working with different data
#need to look at the R-hats for the submodels (see code above)

#plot convergence
png(paste0(study2.graphics, '/S2_TIRED_lv1_traceplots.png'), 
    units='in', width = 8, height=8, res=300)
plot(TIRED_lv1, pars = c("^b_", 'sigma', 'sd'))
dev.off()

#plot recovery of original distribution
png(paste0(study2.graphics, '/S2_TIRED_lv1_pp_dens.png'), 
    units='in', width = 8, height=8, res=300)
pp_check(TIRED_lv1, nsamples=100)
dev.off()

#Storing most important aspects of Summary
bayes.to.txt(model = TIRED_lv1, 
             out.folder = study2.model,
             file = paste0('S2_TIRED_lv1_Summary_', Sys.Date()),
             DF = dat.study2.list[[1]], 
             tot.pars = 7, 
             impute = TRUE)

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#Tired Affect level 2 model 
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
TIRED_lv2<-brms::brm_multiple(T.TIRED~1+DN_comb+(1|ID), 
                              data = dat.study2.list, 
                              family = 'lognormal',
                              prior = c(set_prior('student_t(3,0,1)', class='sd'), 
                                        set_prior('normal(0,2)', class='Intercept')),
                              warmup = 3000, 
                              iter = 6000,
                              thin = 3, 
                              chains = 4,
                              control = list(adapt_delta=.90, 
                                             max_treedepth=10), 
                              save_model = paste0(stan.code, '/S2_Mood_cross.stan'))

#Examine Convergence in each imputed data set
round(TIRED_lv2$rhats[,1:4], 3)   #Can get random effects as well 
summary(TIRED_lv2)                #note that overall R-hat "looks" bad
#that is due to the fact that each model is working with different data
#need to look at the R-hats for the submodels (see code above)

#plot convergence
png(paste0(study2.graphics, '/S2_TIRED_lv2_traceplots.png'), 
    units='in', width = 8, height=8, res=300)
plot(TIRED_lv2, pars = c("b_Intercept",'b_DN_comb', 'sigma', 'sd'))
dev.off()

#plot recovery of original distribution
png(paste0(study2.graphics, '/S2_TIRED_lv2_pp_dens.png'), 
    units='in', width = 8, height=8, res=300)
pp_check(TIRED_lv2, nsamples=100)
dev.off()

#Storing most important aspects of Summary
bayes.to.txt(model = TIRED_lv2, 
             out.folder = study2.model,
             file = paste0('S2_TIRED_lv2_Summary_', Sys.Date()),
             DF = dat.study2.list[[1]], 
             tot.pars = 4, 
             impute = TRUE)

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#Tired Affect cross-level
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
TIRED_cross<-brms::brm_multiple(T.TIRED~1+DN_comb*T.PosEvnt + DN_comb*T.NegEvnt + 
                                  (1 + PosEvnt + NegEvnt|ID), 
                                data = dat.study2.list, 
                                family = 'lognormal',
                                prior = c(set_prior('student_t(3,0,1)', class='sd'), 
                                          set_prior('normal(0,2)', class='Intercept'), 
                                          set_prior('lkj(2)', class='cor')),
                                warmup = 3000, 
                                iter = 6000,
                                thin = 3, 
                                chains = 4,
                                control = list(adapt_delta=.90, 
                                               max_treedepth=10), 
                                save_model = paste0(stan.code, '/S2_Mood_lv2.stan'))


#Examine Convergence in each imputed data set
round(TIRED_cross$rhats[,1:13], 3)   #Can get random effects as well 
summary(TIRED_cross)                #note that overall R-hat "looks" bad
#that is due to the fact that each model is working with different data
#need to look at the R-hats for the submodels (see code above)

#plot convergence
png(paste0(study2.graphics, '/S2_TIRED_cross_traceplots.png'), 
    units='in', width = 8, height=8, res=300)
plot(TIRED_cross, pars = c("b_Intercept",'b_DN_comb', 'sigma', 'sd'))
dev.off()

#plot recovery of original distribution
png(paste0(study2.graphics, '/S2_TIRED_cross_pp_dens.png'), 
    units='in', width = 8, height=8, res=300)
pp_check(TIRED_cross, nsamples=100)
dev.off()

#Storing most important aspects of Summary
bayes.to.txt(model = TIRED_cross, 
             out.folder = study2.model,
             file = paste0('S2_TIRED_cross_Summary_', Sys.Date()),
             DF = dat.study2.list[[1]], 
             tot.pars = 13, 
             impute = TRUE)

#---------------------------------------------------------------------------------------------------------
#DEPRESSED AFFECT MODEL
#---------------------------------------------------------------------------------------------------------

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#Depressed Affect unconditional model (i.e., random intercepts only)
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
DEP_ucm<-brms::brm_multiple(T.DEP~1+(1|ID), 
                            data = dat.study2.list, 
                            family = 'lognormal',
                            prior = c(set_prior('student_t(3,0,1)', class='sd'), 
                                      set_prior('normal(0,2)', class='Intercept')),
                            warmup = 3000, 
                            iter = 6000,
                            thin = 3, 
                            chains = 4,
                            control = list(adapt_delta=.90, 
                                           max_treedepth=10), 
                            save_model = paste0(stan.code, '/S2_Mood_ucm.stan'))


#Examine Convergence in each imputed data set
round(DEP_ucm$rhats[,1:3], 3)   #Can get random effects as well 
summary(DEP_ucm)                #note that overall R-hat "looks" bad
#that is due to the fact that each model is working with different data
#need to look at the R-hats for the submodels (see code above)

#plot convergence
png(paste0(study2.graphics, '/S2_DEP_ucm_traceplots.png'), 
    units='in', width = 8, height=8, res=300)
plot(DEP_ucm, pars = c("b_Intercept", 'sigma', 'sd'))
dev.off()

#plot recovery of original distribution
png(paste0(study2.graphics, '/S2_DEP_ucm_pp_dens.png'), 
    units='in', width = 8, height=8, res=300)
pp_check(DEP_ucm, nsamples=100)
dev.off()

#Storing most important aspects of Summary
bayes.to.txt(model = DEP_ucm, 
             out.folder = study2.model,
             file = paste0('S2_DEP_ucm_Summary_', Sys.Date()),
             DF = dat.study2.list[[1]], 
             tot.pars = 3, 
             impute = TRUE)

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#Depressed Affect level 1 model 
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
DEP_lv1<-brms::brm_multiple(T.DEP~1 + T.PosEvnt + T.NegEvnt + 
                              (1 + T.PosEvnt + T.NegEvnt|ID), 
                            data = dat.study2.list, 
                            family = 'lognormal',
                            prior = c(set_prior('student_t(3,0,1)', class='sd'), 
                                      set_prior('normal(0,2)', class='Intercept'), 
                                      set_prior('lkj(2)', class='cor')),
                            warmup = 3000, 
                            iter = 6000,
                            thin = 3, 
                            chains = 4,
                            control = list(adapt_delta=.90, 
                                           max_treedepth=10), 
                            save_model = paste0(stan.code, '/S2_Mood_lv1.stan'))


#Examine Convergence in each imputed data set
round(DEP_lv1$rhats[,1:7], 3)   #Can get random effects as well 
summary(DEP_lv1)                #note that overall R-hat "looks" bad
#that is due to the fact that each model is working with different data
#need to look at the R-hats for the submodels (see code above)

#plot convergence
png(paste0(study2.graphics, '/S2_DEP_lv1_traceplots.png'), 
    units='in', width = 8, height=8, res=300)
plot(DEP_lv1, pars = c("^b_", 'sigma', 'sd'))
dev.off()

#plot recovery of original distribution
png(paste0(study2.graphics, '/S2_DEP_lv1_pp_dens.png'), 
    units='in', width = 8, height=8, res=300)
pp_check(DEP_lv1, nsamples=100)
dev.off()

#Storing most important aspects of Summary
bayes.to.txt(model = DEP_lv1, 
             out.folder = study2.model,
             file = paste0('S2_DEP_lv1_Summary_', Sys.Date()),
             DF = dat.study2.list[[1]], 
             tot.pars = 7, 
             impute = TRUE)

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#Depressed Affect level 2 model 
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
DEP_lv2<-brms::brm_multiple(T.DEP~1+DN_comb+(1|ID), 
                            data = dat.study2.list, 
                            family = 'lognormal',
                            prior = c(set_prior('student_t(3,0,1)', class='sd'), 
                                      set_prior('normal(0,2)', class='Intercept')),
                            warmup = 3000, 
                            iter = 6000,
                            thin = 3, 
                            chains = 4,
                            control = list(adapt_delta=.90, 
                                           max_treedepth=10), 
                            save_model = paste0(stan.code, '/S2_Mood_cross.stan'))

#Examine Convergence in each imputed data set
round(DEP_lv2$rhats[,1:4], 3)   #Can get random effects as well 
summary(DEP_lv2)                #note that overall R-hat "looks" bad
#that is due to the fact that each model is working with different data
#need to look at the R-hats for the submodels (see code above)

#plot convergence
png(paste0(study2.graphics, '/S2_DEP_lv2_traceplots.png'), 
    units='in', width = 8, height=8, res=300)
plot(DEP_lv2, pars = c("b_Intercept",'b_DN_comb', 'sigma', 'sd'))
dev.off()

#plot recovery of original distribution
png(paste0(study2.graphics, '/S2_DEP_lv2_pp_dens.png'), 
    units='in', width = 8, height=8, res=300)
pp_check(DEP_lv2, nsamples=100)
dev.off()

#Storing most important aspects of Summary
bayes.to.txt(model = DEP_lv2, 
             out.folder = study2.model,
             file = paste0('S2_DEP_lv2_Summary_', Sys.Date()),
             DF = dat.study2.list[[1]], 
             tot.pars = 4, 
             impute = TRUE)

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#Depressed Affect cross-level
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
DEP_cross<-brms::brm_multiple(T.DEP~1+DN_comb*T.PosEvnt + DN_comb*T.NegEvnt + 
                                (1 + PosEvnt + NegEvnt|ID), 
                              data = dat.study2.list, 
                              family = 'lognormal',
                              prior = c(set_prior('student_t(3,0,1)', class='sd'), 
                                        set_prior('normal(0,2)', class='Intercept'), 
                                        set_prior('lkj(2)', class='cor')),
                              warmup = 3000, 
                              iter = 6000,
                              thin = 3, 
                              chains = 4,
                              control = list(adapt_delta=.90, 
                                             max_treedepth=10), 
                              save_model = paste0(stan.code, '/S2_Mood_lv2.stan'))
  
#Examine Convergence in each imputed data set
round(DEP_cross$rhats[,1:13], 3)   #Can get random effects as well 
summary(DEP_cross)                #note that overall R-hat "looks" bad
#that is due to the fact that each model is working with different data
#need to look at the R-hats for the submodels (see code above)

#plot convergence
png(paste0(study2.graphics, '/S2_DEP_cross_traceplots.png'), 
    units='in', width = 8, height=8, res=300)
plot(DEP_cross, pars = c("b_Intercept",'b_DN_comb', 'sigma', 'sd'))
dev.off()

#plot recovery of original distribution
png(paste0(study2.graphics, '/S2_DEP_cross_pp_dens.png'), 
    units='in', width = 8, height=8, res=300)
pp_check(DEP_cross, nsamples=100)
dev.off()

#Storing most important aspects of Summary
bayes.to.txt(model = DEP_cross, 
             out.folder = study2.model,
             file = paste0('S2_DEP_cross_Summary_', Sys.Date()),
             DF = dat.study2.list[[1]], 
             tot.pars = 13, 
             impute = TRUE)

#---------------------------------------------------------------------------------------------------------
#POSITIVE AFFECT MODEL
#---------------------------------------------------------------------------------------------------------

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#Positive Affect unconditional model (i.e., random intercepts only)
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
POS_ucm<-brms::brm_multiple(T.POS~1+(1|ID), 
                            data = dat.study2.list, 
                            family = 'normal',
                            prior = c(set_prior('student_t(3,0,1)', class='sd'), 
                                      set_prior('normal(0,2)', class='Intercept')),
                            warmup = 3000, 
                            iter = 6000,
                            thin = 3, 
                            chains = 4,
                            control = list(adapt_delta=.90, 
                                           max_treedepth=10), 
                            save_model = paste0(stan.code, '/S2_Mood_ucm.stan'))


#Examine Convergence in each imputed data set
round(POS_ucm$rhats[,1:3], 3)   #Can get random effects as well 
summary(POS_ucm)                #note that overall R-hat "looks" bad
#that is due to the fact that each model is working with different data
#need to look at the R-hats for the submodels (see code above)

#plot convergence
png(paste0(study2.graphics, '/S2_POS_ucm_traceplots.png'), 
    units='in', width = 8, height=8, res=300)
plot(POS_ucm, pars = c("b_Intercept", 'sigma', 'sd'))
dev.off()

#plot recovery of original distribution
png(paste0(study2.graphics, '/S2_POS_ucm_pp_dens.png'), 
    units='in', width = 8, height=8, res=300)
pp_check(POS_ucm, nsamples=100)
dev.off()

#Storing most important aspects of Summary
bayes.to.txt(model = POS_ucm, 
             out.folder = study2.model,
             file = paste0('S2_POS_ucm_Summary_', Sys.Date()),
             DF = dat.study2.list[[1]], 
             tot.pars = 3, 
             impute = TRUE)

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#Positive Affect level 1 model 
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
POS_lv1<-brms::brm_multiple(T.POS~1 + T.PosEvnt + T.NegEvnt + 
                              (1 + T.PosEvnt + T.NegEvnt|ID), 
                            data = dat.study2.list, 
                            family = 'normal',
                            prior = c(set_prior('student_t(3,0,1)', class='sd'), 
                                      set_prior('normal(0,2)', class='Intercept'), 
                                      set_prior('lkj(2)', class='cor')),
                            warmup = 3000, 
                            iter = 6000,
                            thin = 3, 
                            chains = 4,
                            control = list(adapt_delta=.90, 
                                           max_treedepth=10), 
                            save_model = paste0(stan.code, '/S2_Mood_lv1.stan'))


#Examine Convergence in each imputed data set
round(POS_lv1$rhats[,1:7], 3)   #Can get random effects as well 
summary(POS_lv1)                #note that overall R-hat "looks" bad
#that is due to the fact that each model is working with different data
#need to look at the R-hats for the submodels (see code above)

#plot convergence
png(paste0(study2.graphics, '/S2_POS_lv1_traceplots.png'), 
    units='in', width = 8, height=8, res=300)
plot(POS_lv1, pars = c("^b_", 'sigma', 'sd'))
dev.off()

#plot recovery of original distribution
png(paste0(study2.graphics, '/S2_POS_lv1_pp_dens.png'), 
    units='in', width = 8, height=8, res=300)
pp_check(POS_lv1, nsamples=100)
dev.off()

#Storing most important aspects of Summary
bayes.to.txt(model = POS_lv1, 
             out.folder = study2.model,
             file = paste0('S2_POS_lv1_Summary_', Sys.Date()),
             DF = dat.study2.list[[1]], 
             tot.pars = 7, 
             impute = TRUE)

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#Positive Affect level 2 model 
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
POS_lv2<-brms::brm_multiple(T.POS~1+DN_comb+(1|ID), 
                            data = dat.study2.list, 
                            family = 'normal',
                            prior = c(set_prior('student_t(3,0,1)', class='sd'), 
                                      set_prior('normal(0,2)', class='Intercept')),
                            warmup = 3000, 
                            iter = 6000,
                            thin = 3, 
                            chains = 4,
                            control = list(adapt_delta=.90, 
                                           max_treedepth=10), 
                            save_model = paste0(stan.code, '/S2_Mood_cross.stan'))

#Examine Convergence in each imputed data set
round(POS_lv2$rhats[,1:4], 3)   #Can get random effects as well 
summary(POS_lv2)                #note that overall R-hat "looks" bad
#that is due to the fact that each model is working with different data
#need to look at the R-hats for the submodels (see code above)

#plot convergence
png(paste0(study2.graphics, '/S2_POS_lv2_traceplots.png'), 
    units='in', width = 8, height=8, res=300)
plot(POS_lv2, pars = c("b_Intercept",'b_DN_comb', 'sigma', 'sd'))
dev.off()

#plot recovery of original distribution
png(paste0(study2.graphics, '/S2_POS_lv2_pp_dens.png'), 
    units='in', width = 8, height=8, res=300)
pp_check(POS_lv2, nsamples=100)
dev.off()

#Storing most important aspects of Summary
bayes.to.txt(model = POS_lv2, 
             out.folder = study2.model,
             file = paste0('S2_POS_lv2_Summary_', Sys.Date()),
             DF = dat.study2.list[[1]], 
             tot.pars = 4, 
             impute = TRUE)

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#Positive Affect cross-level
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
POS_cross<-brms::brm_multiple(T.POS~1+DN_comb*T.PosEvnt + DN_comb*T.NegEvnt + 
                                (1 + PosEvnt + NegEvnt|ID), 
                              data = dat.study2.list, 
                              family = 'normal',
                              prior = c(set_prior('student_t(3,0,1)', class='sd'), 
                                        set_prior('normal(0,2)', class='Intercept'), 
                                        set_prior('lkj(2)', class='cor')),
                              warmup = 3000, 
                              iter = 6000,
                              thin = 3, 
                              chains = 4,
                              control = list(adapt_delta=.90, 
                                             max_treedepth=10), 
                              save_model = paste0(stan.code, '/S2_Mood_lv2.stan'))
############################################################################################################
#Authors:     Barstead, M.G., DeYoung, K.D., Anderson, A. S., & Shackman, A. J.

#Title:       The moment-to-moment affective experience of dispositionally negative individuals 

#Contact:     barstead@umd.edu

#Contents:    Bayesian Regression Analyses - Study 2
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
study2.out<-paste0(wd, '/Study2')
study2.graphics<-paste0(study2.out, '/Graphics')
study2.model<-paste0(study2.out, '/Model_Summaries')
stan.code<-paste0(wd, '/Stan_code')

dat.lv1<-read.csv(paste0(data.folder, '/MB_PAX_T1_Baseline_EMA_FINAL MASTER_070318.csv'), 
                  header = TRUE, stringsAsFactors = FALSE)

dat.lv1$ID.int<-as.numeric(substr(dat.lv1$ID, start = 4, stop=6))
dat.lv1<-dat.lv1[dat.lv1$STUDY.YEAR=='Y1',]   #selecting only year 1

dat.lv1.use<-dat.lv1[dat.lv1$AvailableForAnalyses==1,]

joy.cols<-5:7
calm.cols<-8:10
anx.cols<-11:13
ang.cols<-14:15
tired.cols<-c(16,18)
dep.cols<-c(17,19)

dat.lv1.trunc<-data.frame(ID=dat.lv1.use$ID, 
                          JOY=rowMeans(dat.lv1.use[,joy.cols])+1,
                          CALM=rowMeans(dat.lv1.use[,calm.cols])+1,
                          ANX=rowMeans(dat.lv1.use[,anx.cols])+1,
                          ANG=rowMeans(dat.lv1.use[,ang.cols])+1,
                          TIRED=rowMeans(dat.lv1.use[,tired.cols])+1,
                          DEP=rowMeans(dat.lv1.use[,dep.cols])+1, 
                          POS=rowMeans(dat.lv1.use[,c(joy.cols, calm.cols)])+1, 
                          NEG=rowMeans(dat.lv1.use[,c(anx.cols, dep.cols)])+1,
                          PosEvnt = dat.lv1.use$EMA_12_PosEvent, 
                          NegEvnt = dat.lv1.use$EMA_13_NegEvent, 
                          PosEvnt.Amnt=dat.lv1.use$EMA_12a_PosAmt, 
                          NegEvnt.Amnt=dat.lv1.use$EMA_13a_NegAmt,
                          PosEvnt.Imp=dat.lv1.use$EMA_12b_PosImport, 
                          NegEvnt.Imp=dat.lv1.use$EMA_13b_NegImport,
                          DAY_Corrected = dat.lv1.use$DAY_Corrected, 
                          SIG_Corrected = dat.lv1.use$SIG_Corrected,
                          stringsAsFactors = FALSE)

psych::pairs.panels(dat.lv1.trunc[,2:9])

dat.lv1.trunc$MB_id<-paste0(dat.lv1.trunc$ID, dat.lv1.trunc$DAY_Corrected, dat.lv1.trunc$SIG_Corrected)

dat.lv1.temp<-data.frame(ID = rep(unique(dat.lv1.trunc$ID), each = 56), 
                         DAY = rep(rep(1:7, each=8), times=length(unique(dat.lv1.trunc$ID))),
                         SIG = rep(101:108), 
                         stringsAsFactors = F)

dat.lv1.temp$MB_id<-paste0(dat.lv1.temp$ID, dat.lv1.temp$DAY, dat.lv1.temp$SIG)

dat.lv1.trunc<-merge(dat.lv1.trunc, dat.lv1.temp, by='MB_id', all=TRUE)
#A little bit of cleanup - trimming down unncessary columns made in process
dat.lv1.trunc<-dat.lv1.trunc[,-17:-19]
colnames(dat.lv1.trunc)[2]<-'ID'

dat.lv1.trunc$ID<-ifelse(!is.na(dat.lv1.trunc$ID), dat.lv1.trunc$ID, substr(dat.lv1.trunc$MB_id, start=1, stop=6))
dat.lv2<-read.csv(paste0(data.folder, '/PAX_yr1_Level2_dataset.csv'), header=TRUE, 
                  stringsAsFactors = FALSE)

colnames(dat.lv2)[1]<-'ID'
dat.lv2.trunc<-dat.lv2[,c(1,3,72)]

dat.study2<-merge(dat.lv1.trunc, dat.lv2.trunc, by='ID')
table(dat.study2$ID)  #Everyone should be at 56 observations now - have added missingness... 

############################################################################################################
#Creating a Function to Finalize Imputed Model Output: 
bayes.to.txt<-function(model=NULL, 
                       out.folder=NULL, 
                       file=NULL, 
                       DF=NULL, 
                       tot.pars=3, 
                       impute=TRUE){
  require(parallel)
  require(benchmarkme)
  sink(paste0(out.folder, '/', file, '.txt'))
  cat('System Information:')
  cat('\n=====================================================================================')
  cat(paste0('\nProcessor:', '\t\t\t', benchmarkme::get_cpu()$model_name))
  cat(paste0('\nNumber of Threads:', '\t\t', detectCores(logical=T)))
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
#inspecting missing patterns in the data: 
md.pattern(dat.study2)  

#Will only be imputing for the Potential Dependent variables
#Create initial matrix for convenience
fml<- JOY + CALM + ANX + ANG + TIRED + DEP + POS + NEG + PosEvnt + NegEvnt ~ 
  1 + DN_comb + Sex + (1|ID)

imp<-panImpute(dat.study2, formula=fml, n.burn=10000, 
               n.iter = 5000, m=20, seed = 0320)

dat<-mitmlComplete(imp, force.list = TRUE)

#Inspecting imputation properties - especially interested in outcome measures  
dat.long<-data.frame()
for(i in 1:20){
  dat.temp<-dat[[i]]
  dat.temp$IMP<-rep(i, length(dat.temp[,1]))
  dat.long<-rbind(dat.long, dat.temp)
}

#Creating variables for plotting imputed vs. original distributions:
dat.study2$IMP<-rep(0, length(dat.study2[,1]))
dat.long<-rbind(dat.study2, dat.long)

dat.long$Orig<-ifelse(dat.long$IMP<1, 'Original', 'Imputed')

#Plotting Joyous Affect Imputation
g1<-ggplot()+
  geom_freqpoly(data = dat.long, aes(x=JOY, group=IMP, color=Orig), stat = 'density')+
  theme(legend.title = element_blank())+
  ggtitle('Imputed Distributions of Momentary Joyous Affect')+
  ylab('Density')+
  xlab('Momentary Joyous Affect')
g1

png(paste0(study2.graphics, '/S2_Imputation_JoyousAffect.png'), res=300, 
    units='in', height=5, width = 8)
g1
dev.off()

#Plotting Calm Affect Imputation
g1<-ggplot()+
  geom_freqpoly(data = dat.long, aes(x=CALM, group=IMP, color=Orig), stat = 'density')+
  theme(legend.title = element_blank())+
  ggtitle('Imputed Distributions of Momentary Calm Affect')+
  ylab('Density')+
  xlab('Momentary Calm Affect')
g1

png(paste0(study2.graphics, '/S2_Imputation_CalmAffect.png'), res=300, 
    units='in', height=5, width = 8)
g1
dev.off()


#Plotting Anxious Affect Imputation
g1<-ggplot()+
  geom_freqpoly(data = dat.long, aes(x=ANX, group=IMP, color=Orig), stat = 'density')+
  theme(legend.title = element_blank())+
  ggtitle('Imputed Distributions of Momentary Anxious Affect')+
  ylab('Density')+
  xlab('Momentary Anxious Affect')
g1

png(paste0(study2.graphics, '/S2_Imputation_AnxiousAffect.png'), res=300, 
    units='in', height=5, width = 8)
g1
dev.off()

#Plotting Angry Affect Imputation
g1<-ggplot()+
  geom_freqpoly(data = dat.long, aes(x=ANG, group=IMP, color=Orig), stat = 'density')+
  theme(legend.title = element_blank())+
  ggtitle('Imputed Distributions of Momentary Angry Affect')+
  ylab('Density')+
  xlab('Momentary Angry Affect')
g1

png(paste0(study2.graphics, '/S2_Imputation_AngryAffect.png'), res=300, 
    units='in', height=5, width = 8)
g1
dev.off()

#Plotting Tired Affect Imputation
g1<-ggplot()+
  geom_freqpoly(data = dat.long, aes(x=TIRED, group=IMP, color=Orig), stat = 'density')+
  theme(legend.title = element_blank())+
  ggtitle('Imputed Distributions of Momentary Tired Affect')+
  ylab('Density')+
  xlab('Momentary Tired Affect')
g1

png(paste0(study2.graphics, '/S2_Imputation_TiredAffect.png'), res=300, 
    units='in', height=5, width = 8)
g1
dev.off()

#Plotting Depressed Affect Imputation
g1<-ggplot()+
  geom_freqpoly(data = dat.long, aes(x=DEP, group=IMP, color=Orig), stat = 'density')+
  theme(legend.title = element_blank())+
  ggtitle('Imputed Distributions of Momentary Depressed Affect')+
  ylab('Density')+
  xlab('Momentary Depressed Affect')
g1

png(paste0(study2.graphics, '/S2_Imputation_DepressedAffect.png'), res=300, 
    units='in', height=5, width = 8)
g1
dev.off()

#Plotting Positive Affect Imputation
g1<-ggplot()+
  geom_freqpoly(data = dat.long, aes(x=POS, group=IMP, color=Orig), stat = 'density')+
  theme(legend.title = element_blank())+
  ggtitle('Imputed Distributions of Momentary Positive Affect')+
  ylab('Density')+
  xlab('Momentary Positive Affect')
g1

png(paste0(study2.graphics, '/S2_Imputation_PositiveAffect.png'), res=300, 
    units='in', height=5, width = 8)
g1
dev.off()

#Plotting Negative Affect Imputation
g1<-ggplot()+
  geom_freqpoly(data = dat.long, aes(x=NEG, group=IMP, color=Orig), stat = 'density')+
  theme(legend.title = element_blank())+
  ggtitle('Imputed Distributions of Momentary Negative Affect')+
  ylab('Density')+
  xlab('Momentary Negative Affect')
g1

png(paste0(study2.graphics, '/S2_Imputation_NegativeAffect.png'), res=300, 
    units='in', height=5, width = 8)
g1
dev.off()

############################################################################################################
#summarizing imputations - checking for convergence
sink(paste0(study2.out, '/S2_Imputation_Summary.txt'))
summary(imp)
psych::describeBy(dat.long, group='IMP')
sink()    #warnings have to do with factor variables for the most part

dat.study2.list<-list()

for(i in 1:20){
  dat.study2.list[[i]]<-dat[[i]]
}

#Looping through and constraining values to be 1 <= x <= 5.  
#Note that 20 is the number of imputed data sets:
for(i in 1:20){
  dat.study2.list[[i]]$T.JOY<-ifelse(dat.study2.list[[i]]$JOY>5, 5, dat.study2.list[[i]]$JOY)
  dat.study2.list[[i]]$T.JOY<-ifelse(dat.study2.list[[i]]$T.JOY<1, 1, dat.study2.list[[i]]$T.JOY)
  dat.study2.list[[i]]$T.CALM<-ifelse(dat.study2.list[[i]]$CALM>5, 5, dat.study2.list[[i]]$CALM)
  dat.study2.list[[i]]$T.CALM<-ifelse(dat.study2.list[[i]]$T.CALM<1, 1, dat.study2.list[[i]]$T.CALM)
  dat.study2.list[[i]]$T.ANX<-ifelse(dat.study2.list[[i]]$ANX>5, 5, dat.study2.list[[i]]$ANX)
  dat.study2.list[[i]]$T.ANX<-ifelse(dat.study2.list[[i]]$T.ANX<1, 1, dat.study2.list[[i]]$T.ANX)
  dat.study2.list[[i]]$T.ANG<-ifelse(dat.study2.list[[i]]$ANG>5, 5, dat.study2.list[[i]]$ANG)
  dat.study2.list[[i]]$T.ANG<-ifelse(dat.study2.list[[i]]$T.ANG<1, 1, dat.study2.list[[i]]$T.ANG)
  dat.study2.list[[i]]$T.TIRED<-ifelse(dat.study2.list[[i]]$TIRED>5, 5, dat.study2.list[[i]]$TIRED)
  dat.study2.list[[i]]$T.TIRED<-ifelse(dat.study2.list[[i]]$T.TIRED<1, 1, dat.study2.list[[i]]$T.TIRED)
  dat.study2.list[[i]]$T.DEP<-ifelse(dat.study2.list[[i]]$DEP>5, 5, dat.study2.list[[i]]$DEP)
  dat.study2.list[[i]]$T.DEP<-ifelse(dat.study2.list[[i]]$T.DEP<1, 1, dat.study2.list[[i]]$T.DEP)
  dat.study2.list[[i]]$T.POS<-ifelse(dat.study2.list[[i]]$POS>5, 5, dat.study2.list[[i]]$POS)
  dat.study2.list[[i]]$T.POS<-ifelse(dat.study2.list[[i]]$T.POS<1, 1, dat.study2.list[[i]]$T.POS)
  dat.study2.list[[i]]$T.NEG<-ifelse(dat.study2.list[[i]]$NEG>5, 5, dat.study2.list[[i]]$NEG)
  dat.study2.list[[i]]$T.NEG<-ifelse(dat.study2.list[[i]]$T.NEG<1, 1, dat.study2.list[[i]]$T.NEG)
  dat.study2.list[[i]]$T.PosEvnt<-ifelse(dat.study2.list[[i]]$PosEvnt<=.5, 0, 1)
  dat.study2.list[[i]]$T.NegEvnt<-ifelse(dat.study2.list[[i]]$NegEvnt<=.5, 0, 1)
}

#---------------------------------------------------------------------------------------------------------
#NEGATIVE EVENT MODEL
#---------------------------------------------------------------------------------------------------------

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#Negative Event unconditional model (i.e., random intercepts only)
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
NegEvnt_ucm<-brms::brm_multiple(T.NegEvnt~1+(1|ID), 
                                data = dat.study2.list, 
                                family = 'bernoulli',
                                prior = c(set_prior('student_t(3,0,1)', class='sd'), 
                                          set_prior('normal(0,2)', class='Intercept')),
                                warmup = 3000, 
                                iter = 6000,
                                thin = 3, 
                                chains = 4,
                                control = list(adapt_delta=.90, 
                                               max_treedepth=10), 
                                save_model = paste0(stan.code, '/S2_Evnt_ucm.stan'))

#Examine Convergence in each imputed data set
round(NegEvnt_ucm$rhats[,1:3], 3)   #Can get random effects as well 
summary(NegEvnt_ucm)                #note that overall R-hat "looks" bad
#that is due to the fact that each model is working with different data
#need to look at the R-hats for the submodels (see code above)

#plot convergence
png(paste0(study2.graphics, '/S2_NegEvnt_ucm_traceplots.png'), 
    units='in', width = 8, height=8, res=300)
plot(NegEvnt_ucm, pars = c("b_Intercept", 'sigma', 'sd'))
dev.off()

#plot recovery of original distribution
png(paste0(study2.graphics, '/S2_NegEvnt_ucm_pp_dens.png'), 
    units='in', width = 8, height=8, res=300)
pp_check(NegEvnt_ucm, nsamples=100)
dev.off()

#Storing most important aspects of Summary
bayes.to.txt(model = NegEvnt_ucm, 
             out.folder = study2.model,
             file = paste0('S2_NegEvnt_ucm_Summary_', Sys.Date()),
             DF = dat.study2.list[[1]], 
             tot.pars = 3, 
             impute = TRUE)

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#Negative Event level 2 model 
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
NegEvnt_lv2<-brms::brm_multiple(T.NegEvnt~1+DN_comb+(1|ID), 
                                data = dat.study2.list, 
                                family = 'bernoulli',
                                prior = c(set_prior('student_t(3,0,1)', class='sd'), 
                                          set_prior('normal(0,2)', class='Intercept')),
                                warmup = 3000, 
                                iter = 6000,
                                thin = 3, 
                                chains = 4,
                                control = list(adapt_delta=.90, 
                                               max_treedepth=10), 
                                save_model = paste0(stan.code, '/S2_Evnt_lv2.stan'))


#Examine Convergence in each imputed data set
round(NegEvnt_lv2$rhats[,1:4], 3)   #Can get random effects as well 
summary(NegEvnt_lv2)                #note that overall R-hat "looks" bad
#that is due to the fact that each model is working with different data
#need to look at the R-hats for the submodels (see code above)

#plot convergence
png(paste0(study2.graphics, '/S2_NegEvnt_lv2_traceplots.png'), 
    units='in', width = 8, height=8, res=300)
plot(NegEvnt_lv2, pars = c("b_Intercept",'b_DN_comb', 'sigma', 'sd'))
dev.off()

#plot recovery of original distribution
png(paste0(study2.graphics, '/S2_NegEvnt_lv2_pp_dens.png'), 
    units='in', width = 8, height=8, res=300)
pp_check(NegEvnt_lv2, nsamples=100)
dev.off()

#Storing most important aspects of Summary
bayes.to.txt(model = NegEvnt_lv2, 
             out.folder = study2.model,
             file = paste0('S2_NegEvnt_lv2_Summary_', Sys.Date()),
             DF = dat.study2.list[[1]], 
             tot.pars = 4, 
             impute = TRUE)

#---------------------------------------------------------------------------------------------------------
#POSITIVE EVENT MODEL
#---------------------------------------------------------------------------------------------------------

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#Positive Event unconditional model (i.e., random intercepts only)
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
PosEvnt_ucm<-brms::brm_multiple(T.PosEvnt~1+(1|ID), 
                                data = dat.study2.list, 
                                family = 'bernoulli',
                                prior = c(set_prior('student_t(3,0,1)', class='sd'), 
                                          set_prior('normal(0,2)', class='Intercept')),
                                warmup = 3000, 
                                iter = 6000,
                                thin = 3, 
                                chains = 4,
                                control = list(adapt_delta=.90, 
                                               max_treedepth=10), 
                                save_model = paste0(stan.code, '/S2_Evnt_ucm.stan'))


#Examine Convergence in each imputed data set
round(PosEvnt_ucm$rhats[,1:3], 3)   #Can get random effects as well 
summary(PosEvnt_ucm)                #note that overall R-hat "looks" bad
#that is due to the fact that each model is working with different data
#need to look at the R-hats for the submodels (see code above)

#plot convergence
png(paste0(study2.graphics, '/S2_PosEvnt_ucm_traceplots.png'), 
    units='in', width = 8, height=8, res=300)
plot(PosEvnt_ucm, pars = c("b_Intercept", 'sigma', 'sd'))
dev.off()

#plot recovery of original distribution
png(paste0(study2.graphics, '/S2_PosEvnt_ucm_pp_dens.png'), 
    units='in', width = 8, height=8, res=300)
pp_check(PosEvnt_ucm, nsamples=100)
dev.off()

#Storing most important aspects of Summary
bayes.to.txt(model = PosEvnt_ucm, 
             out.folder = study2.model,
             file = paste0('S2_PosEvnt_ucm_Summary_', Sys.Date()),
             DF = dat.study2.list[[1]], 
             tot.pars = 3, 
             impute = TRUE)

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#Positive Event level 2 model 
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
PosEvnt_lv2<-brms::brm_multiple(T.PosEvnt~1+DN_comb+(1|ID), 
                                data = dat.study2.list, 
                                family = 'bernoulli',
                                prior = c(set_prior('student_t(3,0,1)', class='sd'), 
                                          set_prior('normal(0,2)', class='Intercept')),
                                warmup = 3000, 
                                iter = 6000,
                                thin = 3, 
                                chains = 4,
                                control = list(adapt_delta=.90, 
                                               max_treedepth=10), 
                                save_model = paste0(stan.code, '/S2_Evnt_lv2.stan'))


#Examine Convergence in each imputed data set
round(PosEvnt_lv2$rhats[,1:4], 3)   #Can get random effects as well 
summary(PosEvnt_lv2)                #note that overall R-hat "looks" bad
#that is due to the fact that each model is working with different data
#need to look at the R-hats for the submodels (see code above)

#plot convergence
png(paste0(study2.graphics, '/S2_PosEvnt_lv2_traceplots.png'), 
    units='in', width = 8, height=8, res=300)
plot(PosEvnt_lv2, pars = c("b_Intercept",'b_DN_comb', 'sigma', 'sd'))
dev.off()

#plot recovery of original distribution
png(paste0(study2.graphics, '/S2_PosEvnt_lv2_pp_dens.png'), 
    units='in', width = 8, height=8, res=300)
pp_check(PosEvnt_lv2, nsamples=100)
dev.off()

#Storing most important aspects of Summary
bayes.to.txt(model = PosEvnt_lv2, 
             out.folder = study2.model,
             file = paste0('S2_PosEvnt_lv2_Summary_', Sys.Date()),
             DF = dat.study2.list[[1]], 
             tot.pars = 4, 
             impute = TRUE)

#---------------------------------------------------------------------------------------------------------
#JOYOUS AFFECT MODEL
#---------------------------------------------------------------------------------------------------------

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#Joyous Affect unconditional model (i.e., random intercepts only)
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
JOY_ucm<-brms::brm_multiple(T.JOY~1+(1|ID), 
                            data = dat.study2.list, 
                            family = 'lognormal',
                            prior = c(set_prior('student_t(3,0,1)', class='sd'), 
                                      set_prior('normal(0,2)', class='Intercept')),
                            warmup = 3000, 
                            iter = 6000,
                            thin = 3, 
                            chains = 4,
                            control = list(adapt_delta=.90, 
                                           max_treedepth=10), 
                            save_model = paste0(stan.code, '/S2_Mood_ucm.stan'))


#Examine Convergence in each imputed data set
round(JOY_ucm$rhats[,1:3], 3)   #Can get random effects as well 
summary(JOY_ucm)                #note that overall R-hat "looks" bad
#that is due to the fact that each model is working with different data
#need to look at the R-hats for the submodels (see code above)

#plot convergence
png(paste0(study2.graphics, '/S2_JOY_ucm_traceplots.png'), 
    units='in', width = 8, height=8, res=300)
plot(JOY_ucm, pars = c("b_Intercept", 'sigma', 'sd'))
dev.off()

#plot recovery of original distribution
png(paste0(study2.graphics, '/S2_JOY_ucm_pp_dens.png'), 
    units='in', width = 8, height=8, res=300)
pp_check(JOY_ucm, nsamples=100)
dev.off()

#Storing most important aspects of Summary
bayes.to.txt(model = JOY_ucm, 
             out.folder = study2.model,
             file = paste0('S2_JOY_ucm_Summary_', Sys.Date()),
             DF = dat.study2.list[[1]], 
             tot.pars = 3, 
             impute = TRUE)

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#Joyous Affect level 1 model 
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
JOY_lv1<-brms::brm_multiple(T.JOY~1 + T.PosEvnt + T.NegEvnt + 
                              (1 + T.PosEvnt + T.NegEvnt|ID), 
                            data = dat.study2.list, 
                            family = 'lognormal',
                            prior = c(set_prior('student_t(3,0,1)', class='sd'), 
                                      set_prior('normal(0,2)', class='Intercept'), 
                                      set_prior('lkj(2)', class='cor')),
                            warmup = 3000, 
                            iter = 6000,
                            thin = 3, 
                            chains = 4,
                            control = list(adapt_delta=.90, 
                                           max_treedepth=10), 
                            save_model = paste0(stan.code, '/S2_Mood_lv1.stan'))


#Examine Convergence in each imputed data set
round(JOY_lv1$rhats[,1:7], 3)   #Can get random effects as well 
summary(JOY_lv1)                #note that overall R-hat "looks" bad
#that is due to the fact that each model is working with different data
#need to look at the R-hats for the submodels (see code above)

#plot convergence
png(paste0(study2.graphics, '/S2_JOY_lv1_traceplots.png'), 
    units='in', width = 8, height=8, res=300)
plot(JOY_lv1, pars = c("^b_", 'sigma', 'sd'), N=7)
dev.off()

#plot recovery of original distribution
png(paste0(study2.graphics, '/S2_JOY_lv1_pp_dens.png'), 
    units='in', width = 8, height=8, res=300)
pp_check(JOY_lv1, nsamples=100)
dev.off()

#Storing most important aspects of Summary
bayes.to.txt(model = JOY_lv1, 
             out.folder = study2.model,
             file = paste0('S2_JOY_lv1_Summary_', Sys.Date()),
             DF = dat.study2.list[[1]], 
             tot.pars = 7, 
             impute = TRUE)

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#Joyous Affect level 2 model 
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
JOY_lv2<-brms::brm_multiple(T.JOY~1+DN_comb+(1|ID), 
                            data = dat.study2.list, 
                            family = 'lognormal',
                            prior = c(set_prior('student_t(3,0,1)', class='sd'), 
                                      set_prior('normal(0,2)', class='Intercept')),
                            warmup = 3000, 
                            iter = 6000,
                            thin = 3, 
                            chains = 4,
                            control = list(adapt_delta=.90, 
                                           max_treedepth=10), 
                            save_model = paste0(stan.code, '/S2_Mood_cross.stan'))

#Examine Convergence in each imputed data set
round(JOY_lv2$rhats[,1:4], 3)   #Can get random effects as well 
summary(JOY_lv2)                #note that overall R-hat "looks" bad
#that is due to the fact that each model is working with different data
#need to look at the R-hats for the submodels (see code above)

#plot convergence
png(paste0(study2.graphics, '/S2_JOY_lv2_traceplots.png'), 
    units='in', width = 8, height=8, res=300)
plot(JOY_lv2, pars = c("b_Intercept",'b_DN_comb', 'sigma', 'sd'))
dev.off()

#plot recovery of original distribution
png(paste0(study2.graphics, '/S2_JOY_lv2_pp_dens.png'), 
    units='in', width = 8, height=8, res=300)
pp_check(JOY_lv2, nsamples=100)
dev.off()

#Storing most important aspects of Summary
bayes.to.txt(model = JOY_lv2, 
             out.folder = study2.model,
             file = paste0('S2_JOY_lv2_Summary_', Sys.Date()),
             DF = dat.study2.list[[1]], 
             tot.pars = 4, 
             impute = TRUE)

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#Joyous Affect cross-level
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
JOY_cross<-brms::brm_multiple(T.JOY~1+DN_comb*T.PosEvnt + DN_comb*T.NegEvnt + 
                                (1 + PosEvnt + NegEvnt|ID), 
                              data = dat.study2.list, 
                              family = 'lognormal',
                              prior = c(set_prior('student_t(3,0,1)', class='sd'), 
                                        set_prior('normal(0,2)', class='Intercept'), 
                                        set_prior('lkj(2)', class='cor')),
                              warmup = 3000, 
                              iter = 6000,
                              thin = 3, 
                              chains = 4,
                              control = list(adapt_delta=.90, 
                                             max_treedepth=10), 
                              save_model = paste0(stan.code, '/S2_Mood_lv2.stan'))

#Examine Convergence in each imputed data set
round(JOY_cross$rhats[,1:13], 3)   #Can get random effects as well 
summary(JOY_cross)                #note that overall R-hat "looks" bad
#that is due to the fact that each model is working with different data
#need to look at the R-hats for the submodels (see code above)

#plot convergence
png(paste0(study2.graphics, '/S2_JOY_cross_traceplots.png'), 
    units='in', width = 8, height=8, res=300)
plot(JOY_cross, pars = c('^b_', 'sigma', 'sd'), N=10)
dev.off()

#plot recovery of original distribution
png(paste0(study2.graphics, '/S2_JOY_cross_pp_dens.png'), 
    units='in', width = 8, height=8, res=300)
pp_check(JOY_cross, nsamples=100)
dev.off()

#Storing most important aspects of Summary
bayes.to.txt(model = JOY_cross, 
             out.folder = study2.model,
             file = paste0('S2_JOY_cross_Summary_', Sys.Date()),
             DF = dat.study2.list[[1]], 
             tot.pars = 13, 
             impute = TRUE)

#---------------------------------------------------------------------------------------------------------
#CALM AFFECT MODEL
#---------------------------------------------------------------------------------------------------------

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#Calm Affect unconditional model (i.e., random intercepts only)
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
CALM_ucm<-brms::brm_multiple(T.CALM~1+(1|ID), 
                             data = dat.study2.list, 
                             family = 'normal',
                             prior = c(set_prior('student_t(3,0,1)', class='sd'), 
                                       set_prior('normal(0,2)', class='Intercept')),
                             warmup = 3000, 
                             iter = 6000,
                             thin = 3, 
                             chains = 4,
                             control = list(adapt_delta=.90, 
                                            max_treedepth=10), 
                             save_model = paste0(stan.code, '/S2_Mood_ucm.stan'))


#Examine Convergence in each imputed data set
round(CALM_ucm$rhats[,1:3], 3)   #Can get random effects as well 
summary(CALM_ucm)                #note that overall R-hat "looks" bad
#that is due to the fact that each model is working with different data
#need to look at the R-hats for the submodels (see code above)

#plot convergence
png(paste0(study2.graphics, '/S2_CALM_ucm_traceplots.png'), 
    units='in', width = 8, height=8, res=300)
plot(CALM_ucm, pars = c("b_Intercept", 'sigma', 'sd'))
dev.off()

#plot recovery of original distribution
png(paste0(study2.graphics, '/S2_CALM_ucm_pp_dens.png'), 
    units='in', width = 8, height=8, res=300)
pp_check(CALM_ucm, nsamples=100)
dev.off()

#Storing most important aspects of Summary
bayes.to.txt(model = CALM_ucm, 
             out.folder = study2.model,
             file = paste0('S2_CALM_ucm_Summary_', Sys.Date()),
             DF = dat.study2.list[[1]], 
             tot.pars = 3, 
             impute = TRUE)

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#Calm Affect level 1 model 
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
CALM_lv1<-brms::brm_multiple(T.CALM~1 + T.PosEvnt + T.NegEvnt + 
                               (1 + T.PosEvnt + T.NegEvnt|ID), 
                             data = dat.study2.list, 
                             family = 'normal',
                             prior = c(set_prior('student_t(3,0,1)', class='sd'), 
                                       set_prior('normal(0,2)', class='Intercept'), 
                                       set_prior('lkj(2)', class='cor')),
                             warmup = 3000, 
                             iter = 6000,
                             thin = 3, 
                             chains = 4,
                             control = list(adapt_delta=.90, 
                                            max_treedepth=10), 
                             save_model = paste0(stan.code, '/S2_Mood_lv1.stan'))


#Examine Convergence in each imputed data set
round(CALM_lv1$rhats[,1:7], 3)   #Can get random effects as well 
summary(CALM_lv1)                #note that overall R-hat "looks" bad
#that is due to the fact that each model is working with different data
#need to look at the R-hats for the submodels (see code above)

#plot convergence
png(paste0(study2.graphics, '/S2_CALM_lv1_traceplots.png'), 
    units='in', width = 8, height=8, res=300)
plot(CALM_lv1, pars = c("^b_", 'sigma', 'sd'), N=7)
dev.off()

#plot recovery of original distribution
png(paste0(study2.graphics, '/S2_CALM_lv1_pp_dens.png'), 
    units='in', width = 8, height=8, res=300)
pp_check(CALM_lv1, nsamples=100)
dev.off()

#Storing most important aspects of Summary
bayes.to.txt(model = CALM_lv1, 
             out.folder = study2.model,
             file = paste0('S2_CALM_lv1_Summary_', Sys.Date()),
             DF = dat.study2.list[[1]], 
             tot.pars = 7, 
             impute = TRUE)

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#Calm Affect level 2 model 
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
CALM_lv2<-brms::brm_multiple(T.CALM~1+DN_comb+(1|ID), 
                             data = dat.study2.list, 
                             family = 'normal',
                             prior = c(set_prior('student_t(3,0,1)', class='sd'), 
                                       set_prior('normal(0,2)', class='Intercept')),
                             warmup = 3000, 
                             iter = 6000,
                             thin = 3, 
                             chains = 4,
                             control = list(adapt_delta=.90, 
                                            max_treedepth=10), 
                             save_model = paste0(stan.code, '/S2_Mood_cross.stan'))

#Examine Convergence in each imputed data set
round(CALM_lv2$rhats[,1:4], 3)   #Can get random effects as well 
summary(CALM_lv2)                #note that overall R-hat "looks" bad
#that is due to the fact that each model is working with different data
#need to look at the R-hats for the submodels (see code above)

#plot convergence
png(paste0(study2.graphics, '/S2_CALM_lv2_traceplots.png'), 
    units='in', width = 8, height=8, res=300)
plot(CALM_lv2, pars = c("b_Intercept",'b_DN_comb', 'sigma', 'sd'))
dev.off()

#plot recovery of original distribution
png(paste0(study2.graphics, '/S2_CALM_lv2_pp_dens.png'), 
    units='in', width = 8, height=8, res=300)
pp_check(CALM_lv2, nsamples=100)
dev.off()

#Storing most important aspects of Summary
bayes.to.txt(model = CALM_lv2, 
             out.folder = study2.model,
             file = paste0('S2_CALM_lv2_Summary_', Sys.Date()),
             DF = dat.study2.list[[1]], 
             tot.pars = 4, 
             impute = TRUE)

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#Calm Affect cross-level
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
CALM_cross<-brms::brm_multiple(T.CALM~1+DN_comb*T.PosEvnt + DN_comb*T.NegEvnt + 
                                 (1 + PosEvnt + NegEvnt|ID), 
                               data = dat.study2.list, 
                               family = 'normal',
                               prior = c(set_prior('student_t(3,0,1)', class='sd'), 
                                         set_prior('normal(0,2)', class='Intercept'), 
                                         set_prior('lkj(2)', class='cor')),
                               warmup = 3000, 
                               iter = 6000,
                               thin = 3, 
                               chains = 4,
                               control = list(adapt_delta=.90, 
                                              max_treedepth=10), 
                               save_model = paste0(stan.code, '/S2_Mood_lv2.stan'))

#Examine Convergence in each imputed data set
round(CALM_cross$rhats[,1:13], 3)   #Can get random effects as well 
summary(CALM_cross)                #note that overall R-hat "looks" bad
#that is due to the fact that each model is working with different data
#need to look at the R-hats for the submodels (see code above)

#plot convergence
png(paste0(study2.graphics, '/S2_CALM_cross_traceplots.png'), 
    units='in', width = 8, height=8, res=300)
plot(CALM_cross, pars = c('^b_', 'sigma', 'sd'), N=10)
dev.off()

#plot recovery of original distribution
png(paste0(study2.graphics, '/S2_CALM_cross_pp_dens.png'), 
    units='in', width = 8, height=8, res=300)
pp_check(CALM_cross, nsamples=100)
dev.off()

#Storing most important aspects of Summary
bayes.to.txt(model = CALM_cross, 
             out.folder = study2.model,
             file = paste0('S2_CALM_cross_Summary_', Sys.Date()),
             DF = dat.study2.list[[1]], 
             tot.pars = 13, 
             impute = TRUE)

#---------------------------------------------------------------------------------------------------------
#ANXIOUS AFFECT MODEL
#---------------------------------------------------------------------------------------------------------

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#Anxious Affect unconditional model (i.e., random intercepts only)
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
ANX_ucm<-brms::brm_multiple(T.ANX~1+(1|ID), 
                            data = dat.study2.list, 
                            family = 'lognormal',
                            prior = c(set_prior('student_t(3,0,1)', class='sd'), 
                                      set_prior('normal(0,2)', class='Intercept')),
                            warmup = 3000, 
                            iter = 6000,
                            thin = 3, 
                            chains = 4,
                            control = list(adapt_delta=.90, 
                                           max_treedepth=10), 
                            save_model = paste0(stan.code, '/S2_Mood_ucm.stan'))


#Examine Convergence in each imputed data set
round(ANX_ucm$rhats[,1:3], 3)   #Can get random effects as well 
summary(ANX_ucm)                #note that overall R-hat "looks" bad
#that is due to the fact that each model is working with different data
#need to look at the R-hats for the submodels (see code above)

#plot convergence
png(paste0(study2.graphics, '/S2_ANX_ucm_traceplots.png'), 
    units='in', width = 8, height=8, res=300)
plot(ANX_ucm, pars = c("b_Intercept", 'sigma', 'sd'))
dev.off()

#plot recovery of original distribution
png(paste0(study2.graphics, '/S2_ANX_ucm_pp_dens.png'), 
    units='in', width = 8, height=8, res=300)
pp_check(ANX_ucm, nsamples=100)
dev.off()

#Storing most important aspects of Summary
bayes.to.txt(model = ANX_ucm, 
             out.folder = study2.model,
             file = paste0('S2_ANX_ucm_Summary_', Sys.Date()),
             DF = dat.study2.list[[1]], 
             tot.pars = 3, 
             impute = TRUE)

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#Anxious Affect level 1 model 
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
ANX_lv1<-brms::brm_multiple(T.ANX~1 + T.PosEvnt + T.NegEvnt + 
                              (1 + T.PosEvnt + T.NegEvnt|ID), 
                            data = dat.study2.list, 
                            family = 'lognormal',
                            prior = c(set_prior('student_t(3,0,1)', class='sd'), 
                                      set_prior('normal(0,2)', class='Intercept'), 
                                      set_prior('lkj(2)', class='cor')),
                            warmup = 3000, 
                            iter = 6000,
                            thin = 3, 
                            chains = 4,
                            control = list(adapt_delta=.90, 
                                           max_treedepth=10), 
                            save_model = paste0(stan.code, '/S2_Mood_lv1.stan'))


#Examine Convergence in each imputed data set
round(ANX_lv1$rhats[,1:7], 3)   #Can get random effects as well 
summary(ANX_lv1)                #note that overall R-hat "looks" bad
#that is due to the fact that each model is working with different data
#need to look at the R-hats for the submodels (see code above)

#plot convergence
png(paste0(study2.graphics, '/S2_ANX_lv1_traceplots.png'), 
    units='in', width = 8, height=8, res=300)
plot(ANX_lv1, pars = c("^b_", 'sigma', 'sd'), N=7)
dev.off()

#plot recovery of original distribution
png(paste0(study2.graphics, '/S2_ANX_lv1_pp_dens.png'), 
    units='in', width = 8, height=8, res=300)
pp_check(ANX_lv1, nsamples=100)
dev.off()

#Storing most important aspects of Summary
bayes.to.txt(model = ANX_lv1, 
             out.folder = study2.model,
             file = paste0('S2_ANX_lv1_Summary_', Sys.Date()),
             DF = dat.study2.list[[1]], 
             tot.pars = 7, 
             impute = TRUE)

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#Anxious Affect level 2 model 
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
ANX_lv2<-brms::brm_multiple(T.ANX~1+DN_comb+(1|ID), 
                            data = dat.study2.list, 
                            family = 'lognormal',
                            prior = c(set_prior('student_t(3,0,1)', class='sd'), 
                                      set_prior('normal(0,2)', class='Intercept')),
                            warmup = 3000, 
                            iter = 6000,
                            thin = 3, 
                            chains = 4,
                            control = list(adapt_delta=.90, 
                                           max_treedepth=10), 
                            save_model = paste0(stan.code, '/S2_Mood_cross.stan'))

#Examine Convergence in each imputed data set
round(ANX_lv2$rhats[,1:4], 3)   #Can get random effects as well 
summary(ANX_lv2)                #note that overall R-hat "looks" bad
#that is due to the fact that each model is working with different data
#need to look at the R-hats for the submodels (see code above)

#plot convergence
png(paste0(study2.graphics, '/S2_ANX_lv2_traceplots.png'), 
    units='in', width = 8, height=8, res=300)
plot(ANX_lv2, pars = c("b_Intercept",'b_DN_comb', 'sigma', 'sd'))
dev.off()

#plot recovery of original distribution
png(paste0(study2.graphics, '/S2_ANX_lv2_pp_dens.png'), 
    units='in', width = 8, height=8, res=300)
pp_check(ANX_lv2, nsamples=100)
dev.off()

#Storing most important aspects of Summary
bayes.to.txt(model = ANX_lv2, 
             out.folder = study2.model,
             file = paste0('S2_ANX_lv2_Summary_', Sys.Date()),
             DF = dat.study2.list[[1]], 
             tot.pars = 4, 
             impute = TRUE)

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#Anxious Affect cross-level
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
ANX_cross<-brms::brm_multiple(T.ANX~1+DN_comb*T.PosEvnt + DN_comb*T.NegEvnt + 
                                (1 + PosEvnt + NegEvnt|ID), 
                              data = dat.study2.list, 
                              family = 'lognormal',
                              prior = c(set_prior('student_t(3,0,1)', class='sd'), 
                                        set_prior('normal(0,2)', class='Intercept'), 
                                        set_prior('lkj(2)', class='cor')),
                              warmup = 3000, 
                              iter = 6000,
                              thin = 3, 
                              chains = 4,
                              control = list(adapt_delta=.90, 
                                             max_treedepth=10), 
                              save_model = paste0(stan.code, '/S2_Mood_lv2.stan'))

#Examine Convergence in each imputed data set
round(ANX_cross$rhats[,1:13], 3)   #Can get random effects as well 
summary(ANX_cross)                #note that overall R-hat "looks" bad
#that is due to the fact that each model is working with different data
#need to look at the R-hats for the submodels (see code above)

#plot convergence
png(paste0(study2.graphics, '/S2_ANX_cross_traceplots.png'), 
    units='in', width = 8, height=8, res=300)
plot(ANX_cross, pars = c('^b_', 'sigma', 'sd'), N=10)
dev.off()

#plot recovery of original distribution
png(paste0(study2.graphics, '/S2_ANX_cross_pp_dens.png'), 
    units='in', width = 8, height=8, res=300)
pp_check(ANX_cross, nsamples=100)
dev.off()

#Storing most important aspects of Summary
bayes.to.txt(model = ANX_cross, 
             out.folder = study2.model,
             file = paste0('S2_ANX_cross_Summary_', Sys.Date()),
             DF = dat.study2.list[[1]], 
             tot.pars = 13, 
             impute = TRUE)

#---------------------------------------------------------------------------------------------------------
#ANGRY AFFECT MODEL
#---------------------------------------------------------------------------------------------------------

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#Angry Affect unconditional model (i.e., random intercepts only)
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
ANG_ucm<-brms::brm_multiple(T.ANG~1+(1|ID), 
                            data = dat.study2.list, 
                            family = 'lognormal',
                            prior = c(set_prior('student_t(3,0,1)', class='sd'), 
                                      set_prior('normal(0,2)', class='Intercept')),
                            warmup = 3000, 
                            iter = 6000,
                            thin = 3, 
                            chains = 4,
                            control = list(adapt_delta=.90, 
                                           max_treedepth=10), 
                            save_model = paste0(stan.code, '/S2_Mood_ucm.stan'))


#Examine Convergence in each imputed data set
round(ANG_ucm$rhats[,1:3], 3)   #Can get random effects as well 
summary(ANG_ucm)                #note that overall R-hat "looks" bad
#that is due to the fact that each model is working with different data
#need to look at the R-hats for the submodels (see code above)

#plot convergence
png(paste0(study2.graphics, '/S2_ANG_ucm_traceplots.png'), 
    units='in', width = 8, height=8, res=300)
plot(ANG_ucm, pars = c("b_Intercept", 'sigma', 'sd'))
dev.off()

#plot recovery of original distribution
png(paste0(study2.graphics, '/S2_ANG_ucm_pp_dens.png'), 
    units='in', width = 8, height=8, res=300)
pp_check(ANG_ucm, nsamples=100)
dev.off()

#Storing most important aspects of Summary
bayes.to.txt(model = ANG_ucm, 
             out.folder = study2.model,
             file = paste0('S2_ANG_ucm_Summary_', Sys.Date()),
             DF = dat.study2.list[[1]], 
             tot.pars = 3, 
             impute = TRUE)

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#Angry Affect level 1 model 
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
ANG_lv1<-brms::brm_multiple(T.ANG~1 + T.PosEvnt + T.NegEvnt + 
                              (1 + T.PosEvnt + T.NegEvnt|ID), 
                            data = dat.study2.list, 
                            family = 'lognormal',
                            prior = c(set_prior('student_t(3,0,1)', class='sd'), 
                                      set_prior('normal(0,2)', class='Intercept'), 
                                      set_prior('lkj(2)', class='cor')),
                            warmup = 3000, 
                            iter = 6000,
                            thin = 3, 
                            chains = 4,
                            control = list(adapt_delta=.90, 
                                           max_treedepth=10), 
                            save_model = paste0(stan.code, '/S2_Mood_lv1.stan'))


#Examine Convergence in each imputed data set
round(ANG_lv1$rhats[,1:7], 3)   #Can get random effects as well 
summary(ANG_lv1)                #note that overall R-hat "looks" bad
#that is due to the fact that each model is working with different data
#need to look at the R-hats for the submodels (see code above)

#plot convergence
png(paste0(study2.graphics, '/S2_ANG_lv1_traceplots.png'), 
    units='in', width = 8, height=8, res=300)
plot(ANG_lv1, pars = c("^b_", 'sigma', 'sd'), N=7)
dev.off()

#plot recovery of original distribution
png(paste0(study2.graphics, '/S2_ANG_lv1_pp_dens.png'), 
    units='in', width = 8, height=8, res=300)
pp_check(ANG_lv1, nsamples=100)
dev.off()

#Storing most important aspects of Summary
bayes.to.txt(model = ANG_lv1, 
             out.folder = study2.model,
             file = paste0('S2_ANG_lv1_Summary_', Sys.Date()),
             DF = dat.study2.list[[1]], 
             tot.pars = 7, 
             impute = TRUE)

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#Angry Affect level 2 model 
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
ANG_lv2<-brms::brm_multiple(T.ANG~1+DN_comb+(1|ID), 
                            data = dat.study2.list, 
                            family = 'lognormal',
                            prior = c(set_prior('student_t(3,0,1)', class='sd'), 
                                      set_prior('normal(0,2)', class='Intercept')),
                            warmup = 3000, 
                            iter = 6000,
                            thin = 3, 
                            chains = 4,
                            control = list(adapt_delta=.90, 
                                           max_treedepth=10), 
                            save_model = paste0(stan.code, '/S2_Mood_cross.stan'))

#Examine Convergence in each imputed data set
round(ANG_lv2$rhats[,1:4], 3)   #Can get random effects as well 
summary(ANG_lv2)                #note that overall R-hat "looks" bad
#that is due to the fact that each model is working with different data
#need to look at the R-hats for the submodels (see code above)

#plot convergence
png(paste0(study2.graphics, '/S2_ANG_lv2_traceplots.png'), 
    units='in', width = 8, height=8, res=300)
plot(ANG_lv2, pars = c("b_Intercept",'b_DN_comb', 'sigma', 'sd'))
dev.off()

#plot recovery of original distribution
png(paste0(study2.graphics, '/S2_ANG_lv2_pp_dens.png'), 
    units='in', width = 8, height=8, res=300)
pp_check(ANG_lv2, nsamples=100)
dev.off()

#Storing most important aspects of Summary
bayes.to.txt(model = ANG_lv2, 
             out.folder = study2.model,
             file = paste0('S2_ANG_lv2_Summary_', Sys.Date()),
             DF = dat.study2.list[[1]], 
             tot.pars = 4, 
             impute = TRUE)

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#Angry Affect cross-level
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
ANG_cross<-brms::brm_multiple(T.ANG~1+DN_comb*T.PosEvnt + DN_comb*T.NegEvnt + 
                                (1 + PosEvnt + NegEvnt|ID), 
                              data = dat.study2.list, 
                              family = 'lognormal',
                              prior = c(set_prior('student_t(3,0,1)', class='sd'), 
                                        set_prior('normal(0,2)', class='Intercept'), 
                                        set_prior('lkj(2)', class='cor')),
                              warmup = 3000, 
                              iter = 6000,
                              thin = 3, 
                              chains = 4,
                              control = list(adapt_delta=.90, 
                                             max_treedepth=10), 
                              save_model = paste0(stan.code, '/S2_Mood_lv2.stan'))

#Examine Convergence in each imputed data set
round(ANG_cross$rhats[,1:13], 3)   #Can get random effects as well 
summary(ANG_cross)                #note that overall R-hat "looks" bad
#that is due to the fact that each model is working with different data
#need to look at the R-hats for the submodels (see code above)

#plot convergence
png(paste0(study2.graphics, '/S2_ANG_cross_traceplots.png'), 
    units='in', width = 8, height=8, res=300)
plot(ANG_cross, pars = c('^b_', 'sigma', 'sd'), N=10)
dev.off()

#plot recovery of original distribution
png(paste0(study2.graphics, '/S2_ANG_cross_pp_dens.png'), 
    units='in', width = 8, height=8, res=300)
pp_check(ANG_cross, nsamples=100)
dev.off()

#Storing most important aspects of Summary
bayes.to.txt(model = ANG_cross, 
             out.folder = study2.model,
             file = paste0('S2_ANG_cross_Summary_', Sys.Date()),
             DF = dat.study2.list[[1]], 
             tot.pars = 13, 
             impute = TRUE)

#---------------------------------------------------------------------------------------------------------
#TIRED AFFECT MODEL
#---------------------------------------------------------------------------------------------------------

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#Tired Affect unconditional model (i.e., random intercepts only)
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
TIRED_ucm<-brms::brm_multiple(T.TIRED~1+(1|ID), 
                              data = dat.study2.list, 
                              family = 'lognormal',
                              prior = c(set_prior('student_t(3,0,1)', class='sd'), 
                                        set_prior('normal(0,2)', class='Intercept')),
                              warmup = 3000, 
                              iter = 6000,
                              thin = 3, 
                              chains = 4,
                              control = list(adapt_delta=.90, 
                                             max_treedepth=10), 
                              save_model = paste0(stan.code, '/S2_Mood_ucm.stan'))


#Examine Convergence in each imputed data set
round(TIRED_ucm$rhats[,1:3], 3)   #Can get random effects as well 
summary(TIRED_ucm)                #note that overall R-hat "looks" bad
#that is due to the fact that each model is working with different data
#need to look at the R-hats for the submodels (see code above)

#plot convergence
png(paste0(study2.graphics, '/S2_TIRED_ucm_traceplots.png'), 
    units='in', width = 8, height=8, res=300)
plot(TIRED_ucm, pars = c("b_Intercept", 'sigma', 'sd'))
dev.off()

#plot recovery of original distribution
png(paste0(study2.graphics, '/S2_TIRED_ucm_pp_dens.png'), 
    units='in', width = 8, height=8, res=300)
pp_check(TIRED_ucm, nsamples=100)
dev.off()

#Storing most important aspects of Summary
bayes.to.txt(model = TIRED_ucm, 
             out.folder = study2.model,
             file = paste0('S2_TIRED_ucm_Summary_', Sys.Date()),
             DF = dat.study2.list[[1]], 
             tot.pars = 3, 
             impute = TRUE)

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#Tired Affect level 1 model 
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
TIRED_lv1<-brms::brm_multiple(T.TIRED~1 + T.PosEvnt + T.NegEvnt + 
                                (1 + T.PosEvnt + T.NegEvnt|ID), 
                              data = dat.study2.list, 
                              family = 'lognormal',
                              prior = c(set_prior('student_t(3,0,1)', class='sd'), 
                                        set_prior('normal(0,2)', class='Intercept'), 
                                        set_prior('lkj(2)', class='cor')),
                              warmup = 3000, 
                              iter = 6000,
                              thin = 3, 
                              chains = 4,
                              control = list(adapt_delta=.90, 
                                             max_treedepth=10), 
                              save_model = paste0(stan.code, '/S2_Mood_lv1.stan'))


#Examine Convergence in each imputed data set
round(TIRED_lv1$rhats[,1:7], 3)   #Can get random effects as well 
summary(TIRED_lv1)                #note that overall R-hat "looks" bad
#that is due to the fact that each model is working with different data
#need to look at the R-hats for the submodels (see code above)

#plot convergence
png(paste0(study2.graphics, '/S2_TIRED_lv1_traceplots.png'), 
    units='in', width = 8, height=8, res=300)
plot(TIRED_lv1, pars = c("^b_", 'sigma', 'sd'), N=7)
dev.off()

#plot recovery of original distribution
png(paste0(study2.graphics, '/S2_TIRED_lv1_pp_dens.png'), 
    units='in', width = 8, height=8, res=300)
pp_check(TIRED_lv1, nsamples=100)
dev.off()

#Storing most important aspects of Summary
bayes.to.txt(model = TIRED_lv1, 
             out.folder = study2.model,
             file = paste0('S2_TIRED_lv1_Summary_', Sys.Date()),
             DF = dat.study2.list[[1]], 
             tot.pars = 7, 
             impute = TRUE)

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#Tired Affect level 2 model 
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
TIRED_lv2<-brms::brm_multiple(T.TIRED~1+DN_comb+(1|ID), 
                              data = dat.study2.list, 
                              family = 'lognormal',
                              prior = c(set_prior('student_t(3,0,1)', class='sd'), 
                                        set_prior('normal(0,2)', class='Intercept')),
                              warmup = 3000, 
                              iter = 6000,
                              thin = 3, 
                              chains = 4,
                              control = list(adapt_delta=.90, 
                                             max_treedepth=10), 
                              save_model = paste0(stan.code, '/S2_Mood_cross.stan'))

#Examine Convergence in each imputed data set
round(TIRED_lv2$rhats[,1:4], 3)   #Can get random effects as well 
summary(TIRED_lv2)                #note that overall R-hat "looks" bad
#that is due to the fact that each model is working with different data
#need to look at the R-hats for the submodels (see code above)

#plot convergence
png(paste0(study2.graphics, '/S2_TIRED_lv2_traceplots.png'), 
    units='in', width = 8, height=8, res=300)
plot(TIRED_lv2, pars = c("b_Intercept",'b_DN_comb', 'sigma', 'sd'))
dev.off()

#plot recovery of original distribution
png(paste0(study2.graphics, '/S2_TIRED_lv2_pp_dens.png'), 
    units='in', width = 8, height=8, res=300)
pp_check(TIRED_lv2, nsamples=100)
dev.off()

#Storing most important aspects of Summary
bayes.to.txt(model = TIRED_lv2, 
             out.folder = study2.model,
             file = paste0('S2_TIRED_lv2_Summary_', Sys.Date()),
             DF = dat.study2.list[[1]], 
             tot.pars = 4, 
             impute = TRUE)

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#Tired Affect cross-level
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
TIRED_cross<-brms::brm_multiple(T.TIRED~1+DN_comb*T.PosEvnt + DN_comb*T.NegEvnt + 
                                  (1 + PosEvnt + NegEvnt|ID), 
                                data = dat.study2.list, 
                                family = 'lognormal',
                                prior = c(set_prior('student_t(3,0,1)', class='sd'), 
                                          set_prior('normal(0,2)', class='Intercept'), 
                                          set_prior('lkj(2)', class='cor')),
                                warmup = 3000, 
                                iter = 6000,
                                thin = 3, 
                                chains = 4,
                                control = list(adapt_delta=.90, 
                                               max_treedepth=10), 
                                save_model = paste0(stan.code, '/S2_Mood_lv2.stan'))


#Examine Convergence in each imputed data set
round(TIRED_cross$rhats[,1:13], 3)   #Can get random effects as well 
summary(TIRED_cross)                #note that overall R-hat "looks" bad
#that is due to the fact that each model is working with different data
#need to look at the R-hats for the submodels (see code above)

#plot convergence
png(paste0(study2.graphics, '/S2_TIRED_cross_traceplots.png'), 
    units='in', width = 8, height=8, res=300)
plot(TIRED_cross, pars = c('^b_', 'sigma', 'sd'), N=10)
dev.off()

#plot recovery of original distribution
png(paste0(study2.graphics, '/S2_TIRED_cross_pp_dens.png'), 
    units='in', width = 8, height=8, res=300)
pp_check(TIRED_cross, nsamples=100)
dev.off()

#Storing most important aspects of Summary
bayes.to.txt(model = TIRED_cross, 
             out.folder = study2.model,
             file = paste0('S2_TIRED_cross_Summary_', Sys.Date()),
             DF = dat.study2.list[[1]], 
             tot.pars = 13, 
             impute = TRUE)

#---------------------------------------------------------------------------------------------------------
#DEPRESSED AFFECT MODEL
#---------------------------------------------------------------------------------------------------------

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#Depressed Affect unconditional model (i.e., random intercepts only)
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
DEP_ucm<-brms::brm_multiple(T.DEP~1+(1|ID), 
                            data = dat.study2.list, 
                            family = 'lognormal',
                            prior = c(set_prior('student_t(3,0,1)', class='sd'), 
                                      set_prior('normal(0,2)', class='Intercept')),
                            warmup = 3000, 
                            iter = 6000,
                            thin = 3, 
                            chains = 4,
                            control = list(adapt_delta=.90, 
                                           max_treedepth=10), 
                            save_model = paste0(stan.code, '/S2_Mood_ucm.stan'))


#Examine Convergence in each imputed data set
round(DEP_ucm$rhats[,1:3], 3)   #Can get random effects as well 
summary(DEP_ucm)                #note that overall R-hat "looks" bad
#that is due to the fact that each model is working with different data
#need to look at the R-hats for the submodels (see code above)

#plot convergence
png(paste0(study2.graphics, '/S2_DEP_ucm_traceplots.png'), 
    units='in', width = 8, height=8, res=300)
plot(DEP_ucm, pars = c("b_Intercept", 'sigma', 'sd'))
dev.off()

#plot recovery of original distribution
png(paste0(study2.graphics, '/S2_DEP_ucm_pp_dens.png'), 
    units='in', width = 8, height=8, res=300)
pp_check(DEP_ucm, nsamples=100)
dev.off()

#Storing most important aspects of Summary
bayes.to.txt(model = DEP_ucm, 
             out.folder = study2.model,
             file = paste0('S2_DEP_ucm_Summary_', Sys.Date()),
             DF = dat.study2.list[[1]], 
             tot.pars = 3, 
             impute = TRUE)

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#Depressed Affect level 1 model 
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
DEP_lv1<-brms::brm_multiple(T.DEP~1 + T.PosEvnt + T.NegEvnt + 
                              (1 + T.PosEvnt + T.NegEvnt|ID), 
                            data = dat.study2.list, 
                            family = 'lognormal',
                            prior = c(set_prior('student_t(3,0,1)', class='sd'), 
                                      set_prior('normal(0,2)', class='Intercept'), 
                                      set_prior('lkj(2)', class='cor')),
                            warmup = 3000, 
                            iter = 6000,
                            thin = 3, 
                            chains = 4,
                            control = list(adapt_delta=.90, 
                                           max_treedepth=10), 
                            save_model = paste0(stan.code, '/S2_Mood_lv1.stan'))


#Examine Convergence in each imputed data set
round(DEP_lv1$rhats[,1:7], 3)   #Can get random effects as well 
summary(DEP_lv1)                #note that overall R-hat "looks" bad
#that is due to the fact that each model is working with different data
#need to look at the R-hats for the submodels (see code above)

#plot convergence
png(paste0(study2.graphics, '/S2_DEP_lv1_traceplots.png'), 
    units='in', width = 8, height=8, res=300)
plot(DEP_lv1, pars = c("^b_", 'sigma', 'sd'), N=7)
dev.off()

#plot recovery of original distribution
png(paste0(study2.graphics, '/S2_DEP_lv1_pp_dens.png'), 
    units='in', width = 8, height=8, res=300)
pp_check(DEP_lv1, nsamples=100)
dev.off()

#Storing most important aspects of Summary
bayes.to.txt(model = DEP_lv1, 
             out.folder = study2.model,
             file = paste0('S2_DEP_lv1_Summary_', Sys.Date()),
             DF = dat.study2.list[[1]], 
             tot.pars = 7, 
             impute = TRUE)

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#Depressed Affect level 2 model 
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
DEP_lv2<-brms::brm_multiple(T.DEP~1+DN_comb+(1|ID), 
                            data = dat.study2.list, 
                            family = 'lognormal',
                            prior = c(set_prior('student_t(3,0,1)', class='sd'), 
                                      set_prior('normal(0,2)', class='Intercept')),
                            warmup = 3000, 
                            iter = 6000,
                            thin = 3, 
                            chains = 4,
                            control = list(adapt_delta=.90, 
                                           max_treedepth=10), 
                            save_model = paste0(stan.code, '/S2_Mood_cross.stan'))

#Examine Convergence in each imputed data set
round(DEP_lv2$rhats[,1:4], 3)   #Can get random effects as well 
summary(DEP_lv2)                #note that overall R-hat "looks" bad
#that is due to the fact that each model is working with different data
#need to look at the R-hats for the submodels (see code above)

#plot convergence
png(paste0(study2.graphics, '/S2_DEP_lv2_traceplots.png'), 
    units='in', width = 8, height=8, res=300)
plot(DEP_lv2, pars = c("b_Intercept",'b_DN_comb', 'sigma', 'sd'))
dev.off()

#plot recovery of original distribution
png(paste0(study2.graphics, '/S2_DEP_lv2_pp_dens.png'), 
    units='in', width = 8, height=8, res=300)
pp_check(DEP_lv2, nsamples=100)
dev.off()

#Storing most important aspects of Summary
bayes.to.txt(model = DEP_lv2, 
             out.folder = study2.model,
             file = paste0('S2_DEP_lv2_Summary_', Sys.Date()),
             DF = dat.study2.list[[1]], 
             tot.pars = 4, 
             impute = TRUE)

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#Depressed Affect cross-level
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
DEP_cross<-brms::brm_multiple(T.DEP~1+DN_comb*T.PosEvnt + DN_comb*T.NegEvnt + 
                                (1 + PosEvnt + NegEvnt|ID), 
                              data = dat.study2.list, 
                              family = 'lognormal',
                              prior = c(set_prior('student_t(3,0,1)', class='sd'), 
                                        set_prior('normal(0,2)', class='Intercept'), 
                                        set_prior('lkj(2)', class='cor')),
                              warmup = 3000, 
                              iter = 6000,
                              thin = 3, 
                              chains = 4,
                              control = list(adapt_delta=.90, 
                                             max_treedepth=10), 
                              save_model = paste0(stan.code, '/S2_Mood_lv2.stan'))

#Examine Convergence in each imputed data set
round(DEP_cross$rhats[,1:13], 3)   #Can get random effects as well 
summary(DEP_cross)                #note that overall R-hat "looks" bad
#that is due to the fact that each model is working with different data
#need to look at the R-hats for the submodels (see code above)

#plot convergence
png(paste0(study2.graphics, '/S2_DEP_cross_traceplots.png'), 
    units='in', width = 8, height=8, res=300)
plot(DEP_cross, pars = c('^b_', 'sigma', 'sd'), N=10)
dev.off()

#plot recovery of original distribution
png(paste0(study2.graphics, '/S2_DEP_cross_pp_dens.png'), 
    units='in', width = 8, height=8, res=300)
pp_check(DEP_cross, nsamples=100)
dev.off()

#Storing most important aspects of Summary
bayes.to.txt(model = DEP_cross, 
             out.folder = study2.model,
             file = paste0('S2_DEP_cross_Summary_', Sys.Date()),
             DF = dat.study2.list[[1]], 
             tot.pars = 13, 
             impute = TRUE)

#---------------------------------------------------------------------------------------------------------
#POSITIVE AFFECT MODEL
#---------------------------------------------------------------------------------------------------------

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#Positive Affect unconditional model (i.e., random intercepts only)
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
POS_ucm<-brms::brm_multiple(T.POS~1+(1|ID), 
                            data = dat.study2.list, 
                            family = 'normal',
                            prior = c(set_prior('student_t(3,0,1)', class='sd'), 
                                      set_prior('normal(0,2)', class='Intercept')),
                            warmup = 3000, 
                            iter = 6000,
                            thin = 3, 
                            chains = 4,
                            control = list(adapt_delta=.90, 
                                           max_treedepth=10), 
                            save_model = paste0(stan.code, '/S2_Mood_ucm.stan'))


#Examine Convergence in each imputed data set
round(POS_ucm$rhats[,1:3], 3)   #Can get random effects as well 
summary(POS_ucm)                #note that overall R-hat "looks" bad
#that is due to the fact that each model is working with different data
#need to look at the R-hats for the submodels (see code above)

#plot convergence
png(paste0(study2.graphics, '/S2_POS_ucm_traceplots.png'), 
    units='in', width = 8, height=8, res=300)
plot(POS_ucm, pars = c("b_Intercept", 'sigma', 'sd'))
dev.off()

#plot recovery of original distribution
png(paste0(study2.graphics, '/S2_POS_ucm_pp_dens.png'), 
    units='in', width = 8, height=8, res=300)
pp_check(POS_ucm, nsamples=100)
dev.off()

#Storing most important aspects of Summary
bayes.to.txt(model = POS_ucm, 
             out.folder = study2.model,
             file = paste0('S2_POS_ucm_Summary_', Sys.Date()),
             DF = dat.study2.list[[1]], 
             tot.pars = 3, 
             impute = TRUE)

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#Positive Affect level 1 model 
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
POS_lv1<-brms::brm_multiple(T.POS~1 + T.PosEvnt + T.NegEvnt + 
                              (1 + T.PosEvnt + T.NegEvnt|ID), 
                            data = dat.study2.list, 
                            family = 'normal',
                            prior = c(set_prior('student_t(3,0,1)', class='sd'), 
                                      set_prior('normal(0,2)', class='Intercept'), 
                                      set_prior('lkj(2)', class='cor')),
                            warmup = 3000, 
                            iter = 6000,
                            thin = 3, 
                            chains = 4,
                            control = list(adapt_delta=.90, 
                                           max_treedepth=10), 
                            save_model = paste0(stan.code, '/S2_Mood_lv1.stan'))


#Examine Convergence in each imputed data set
round(POS_lv1$rhats[,1:7], 3)   #Can get random effects as well 
summary(POS_lv1)                #note that overall R-hat "looks" bad
#that is due to the fact that each model is working with different data
#need to look at the R-hats for the submodels (see code above)

#plot convergence
png(paste0(study2.graphics, '/S2_POS_lv1_traceplots.png'), 
    units='in', width = 8, height=8, res=300)
plot(POS_lv1, pars = c("^b_", 'sigma', 'sd'), N=7)
dev.off()

#plot recovery of original distribution
png(paste0(study2.graphics, '/S2_POS_lv1_pp_dens.png'), 
    units='in', width = 8, height=8, res=300)
pp_check(POS_lv1, nsamples=100)
dev.off()

#Storing most important aspects of Summary
bayes.to.txt(model = POS_lv1, 
             out.folder = study2.model,
             file = paste0('S2_POS_lv1_Summary_', Sys.Date()),
             DF = dat.study2.list[[1]], 
             tot.pars = 7, 
             impute = TRUE)

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#Positive Affect level 2 model 
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
POS_lv2<-brms::brm_multiple(T.POS~1+DN_comb+(1|ID), 
                            data = dat.study2.list, 
                            family = 'normal',
                            prior = c(set_prior('student_t(3,0,1)', class='sd'), 
                                      set_prior('normal(0,2)', class='Intercept')),
                            warmup = 3000, 
                            iter = 6000,
                            thin = 3, 
                            chains = 4,
                            control = list(adapt_delta=.90, 
                                           max_treedepth=10), 
                            save_model = paste0(stan.code, '/S2_Mood_cross.stan'))

#Examine Convergence in each imputed data set
round(POS_lv2$rhats[,1:4], 3)   #Can get random effects as well 
summary(POS_lv2)                #note that overall R-hat "looks" bad
#that is due to the fact that each model is working with different data
#need to look at the R-hats for the submodels (see code above)

#plot convergence
png(paste0(study2.graphics, '/S2_POS_lv2_traceplots.png'), 
    units='in', width = 8, height=8, res=300)
plot(POS_lv2, pars = c("b_Intercept",'b_DN_comb', 'sigma', 'sd'))
dev.off()

#plot recovery of original distribution
png(paste0(study2.graphics, '/S2_POS_lv2_pp_dens.png'), 
    units='in', width = 8, height=8, res=300)
pp_check(POS_lv2, nsamples=100)
dev.off()

#Storing most important aspects of Summary
bayes.to.txt(model = POS_lv2, 
             out.folder = study2.model,
             file = paste0('S2_POS_lv2_Summary_', Sys.Date()),
             DF = dat.study2.list[[1]], 
             tot.pars = 4, 
             impute = TRUE)

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#Positive Affect cross-level
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
POS_cross<-brms::brm_multiple(T.POS~1+DN_comb*T.PosEvnt + DN_comb*T.NegEvnt + 
                                (1 + PosEvnt + NegEvnt|ID), 
                              data = dat.study2.list, 
                              family = 'normal',
                              prior = c(set_prior('student_t(3,0,1)', class='sd'), 
                                        set_prior('normal(0,2)', class='Intercept'), 
                                        set_prior('lkj(2)', class='cor')),
                              warmup = 3000, 
                              iter = 6000,
                              thin = 3, 
                              chains = 4,
                              control = list(adapt_delta=.90, 
                                             max_treedepth=10), 
                              save_model = paste0(stan.code, '/S2_Mood_lv2.stan'))

#Examine Convergence in each imputed data set
round(POS_cross$rhats[,1:13], 3)   #Can get random effects as well 
summary(POS_cross)                #note that overall R-hat "looks" bad
#that is due to the fact that each model is working with different data
#need to look at the R-hats for the submodels (see code above)

#plot convergence
png(paste0(study2.graphics, '/S2_POS_cross_traceplots.png'), 
    units='in', width = 8, height=8, res=300)
plot(POS_cross, pars = c('^b_', 'sigma', 'sd'), N=10)
dev.off()

#plot recovery of original distribution
png(paste0(study2.graphics, '/S2_POS_cross_pp_dens.png'), 
    units='in', width = 8, height=8, res=300)
pp_check(POS_cross, nsamples=100)
dev.off()

#Storing most important aspects of Summary
bayes.to.txt(model = POS_cross, 
             out.folder = study2.model,
             file = paste0('S2_POS_cross_Summary_', Sys.Date()),
             DF = dat.study2.list[[1]], 
             tot.pars = 13, 
             impute = TRUE)

#---------------------------------------------------------------------------------------------------------
#NEGATIVE AFFECT MODEL
#---------------------------------------------------------------------------------------------------------

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#Negative Affect unconditional model (i.e., random intercepts only)
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
NEG_ucm<-brms::brm_multiple(T.NEG~1+(1|ID), 
                            data = dat.study2.list, 
                            family = 'lognormal',
                            prior = c(set_prior('student_t(3,0,1)', class='sd'), 
                                      set_prior('normal(0,2)', class='Intercept')),
                            warmup = 3000, 
                            iter = 6000,
                            thin = 3, 
                            chains = 4,
                            control = list(adapt_delta=.90, 
                                           max_treedepth=10), 
                            save_model = paste0(stan.code, '/S2_Mood_ucm.stan'))


#Examine Convergence in each imputed data set
round(NEG_ucm$rhats[,1:3], 3)   #Can get random effects as well 
summary(NEG_ucm)                #note that overall R-hat "looks" bad
#that is due to the fact that each model is working with different data
#need to look at the R-hats for the submodels (see code above)

#plot convergence
png(paste0(study2.graphics, '/S2_NEG_ucm_traceplots.png'), 
    units='in', width = 8, height=8, res=300)
plot(NEG_ucm, pars = c("b_Intercept", 'sigma', 'sd'))
dev.off()

#plot recovery of original distribution
png(paste0(study2.graphics, '/S2_NEG_ucm_pp_dens.png'), 
    units='in', width = 8, height=8, res=300)
pp_check(NEG_ucm, nsamples=100)
dev.off()

#Storing most important aspects of Summary
bayes.to.txt(model = NEG_ucm, 
             out.folder = study2.model,
             file = paste0('S2_NEG_ucm_Summary_', Sys.Date()),
             DF = dat.study2.list[[1]], 
             tot.pars = 3, 
             impute = TRUE)

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#Negative Affect level 1 model 
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
NEG_lv1<-brms::brm_multiple(T.NEG~1 + T.PosEvnt + T.NegEvnt + 
                              (1 + T.PosEvnt + T.NegEvnt|ID), 
                            data = dat.study2.list, 
                            family = 'lognormal',
                            prior = c(set_prior('student_t(3,0,1)', class='sd'), 
                                      set_prior('normal(0,2)', class='Intercept'), 
                                      set_prior('lkj(2)', class='cor')),
                            warmup = 3000, 
                            iter = 6000,
                            thin = 3, 
                            chains = 4,
                            control = list(adapt_delta=.90, 
                                           max_treedepth=10), 
                            save_model = paste0(stan.code, '/S2_Mood_lv1.stan'))


#Examine Convergence in each imputed data set
round(NEG_lv1$rhats[,1:7], 3)   #Can get random effects as well 
summary(NEG_lv1)                #note that overall R-hat "looks" bad
#that is due to the fact that each model is working with different data
#need to look at the R-hats for the submodels (see code above)

#plot convergence
png(paste0(study2.graphics, '/S2_NEG_lv1_traceplots.png'), 
    units='in', width = 8, height=8, res=300)
plot(NEG_lv1, pars = c("^b_", 'sigma', 'sd'), N=7)
dev.off()

#plot recovery of original distribution
png(paste0(study2.graphics, '/S2_NEG_lv1_pp_dens.png'), 
    units='in', width = 8, height=8, res=300)
pp_check(NEG_lv1, nsamples=100)
dev.off()

#Storing most important aspects of Summary
bayes.to.txt(model = NEG_lv1, 
             out.folder = study2.model,
             file = paste0('S2_NEG_lv1_Summary_', Sys.Date()),
             DF = dat.study2.list[[1]], 
             tot.pars = 7, 
             impute = TRUE)

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#Negative Affect level 2 model 
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
NEG_lv2<-brms::brm_multiple(T.NEG~1+DN_comb+(1|ID), 
                              data = dat.study2.list, 
                              family = 'lognormal',
                              prior = c(set_prior('student_t(3,0,1)', class='sd'), 
                                        set_prior('normal(0,2)', class='Intercept')),
                              warmup = 3000, 
                              iter = 6000,
                              thin = 3, 
                              chains = 4,
                              control = list(adapt_delta=.90, 
                                             max_treedepth=10), 
                              save_model = paste0(stan.code, '/S2_Mood_cross.stan'))

#Examine Convergence in each imputed data set
round(NEG_lv2$rhats[,1:4], 3)   #Can get random effects as well 
summary(NEG_lv2)                #note that overall R-hat "looks" bad
#that is due to the fact that each model is working with different data
#need to look at the R-hats for the submodels (see code above)

#plot convergence
png(paste0(study2.graphics, '/S2_NEG_lv2_traceplots.png'), 
    units='in', width = 8, height=8, res=300)
plot(NEG_lv2, pars = c("b_Intercept",'b_DN_comb', 'sigma', 'sd'))
dev.off()

#plot recovery of original distribution
png(paste0(study2.graphics, '/S2_NEG_lv2_pp_dens.png'), 
    units='in', width = 8, height=8, res=300)
pp_check(NEG_lv2, nsamples=100)
dev.off()

#Storing most important aspects of Summary
bayes.to.txt(model = NEG_lv2, 
             out.folder = study2.model,
             file = paste0('S2_NEG_lv2_Summary_', Sys.Date()),
             DF = dat.study2.list[[1]], 
             tot.pars = 4, 
             impute = TRUE)

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#Negative Affect cross-level
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
NEG_cross<-brms::brm_multiple(T.NEG~1+DN_comb*T.PosEvnt + DN_comb*T.NegEvnt + 
                              (1 + PosEvnt + NegEvnt|ID), 
                            data = dat.study2.list, 
                            family = 'lognormal',
                            prior = c(set_prior('student_t(3,0,1)', class='sd'), 
                                      set_prior('normal(0,2)', class='Intercept'), 
                                      set_prior('lkj(2)', class='cor')),
                            warmup = 3000, 
                            iter = 6000,
                            thin = 3, 
                            chains = 4,
                            control = list(adapt_delta=.90, 
                                           max_treedepth=10), 
                            save_model = paste0(stan.code, '/S2_Mood_lv2.stan'))



#Examine Convergence in each imputed data set
round(NEG_cross$rhats[,1:13], 3)   #Can get random effects as well 
summary(NEG_cross)                #note that overall R-hat "looks" bad
#that is due to the fact that each model is working with different data
#need to look at the R-hats for the submodels (see code above)

#plot convergence
png(paste0(study2.graphics, '/S2_NEG_cross_traceplots.png'), 
    units='in', width = 8, height=11, res=300)
plot(NEG_cross, pars = c("^b_",'sigma', 'sd'), N=10)
dev.off()

#plot recovery of original distribution
png(paste0(study2.graphics, '/S2_NEG_cross_pp_dens.png'), 
    units='in', width = 8, height=8, res=300)
pp_check(NEG_cross, nsamples=100)
dev.off()