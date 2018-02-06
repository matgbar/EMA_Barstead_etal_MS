#Barstead et al - EMA manuscript
X<-cbind(rep(1, length(DF.stan.lv1[,1])),
         DF.stan.lv1$PosOnly,
         DF.stan.lv1$NegOnly,
         DF.stan.lv1$Neg_Pos)
colnames(X)<-paste0('x', 1:4)

u1<-vector()
u2<-vector()
u3<-vector()

#note standardizing all values to have unit variance and mean of 0
#the rationale is to improve intreptability of final estimates
for(i in 1:length(IDs)){
  u1[i]<-DF.stan.lv1$DN_comb[DF.stan.lv1$id==IDs[i]][1]
  u2[i]<-DF.stan.lv1$Pos.prop[DF.stan.lv1$id==IDs[i]][1]
  u3[i]<-DF.stan.lv1$Neg.prop[DF.stan.lv1$id==IDs[i]][1]
}

U<-cbind(rep(1, length(IDs)),
         as.numeric(scale(u1)),
         as.numeric(scale(u2)),
         as.numeric(scale(u3))
)
DF.stan.lv1$ANX_r<-DF.stan.lv1$ANX+1
dat.complete<-list(N=length(DF.stan.lv1[,1]),
                   K=ncol(X),
                   J=length(IDs),
                   L=ncol(U),
                   x=X,
                   u=U,
                   ID=DF.stan.lv1$ID.bayes,
                   y=DF.stan.lv1$ANX_r
)

pars.monitor<-c('Sigma_beta',
                'sigma',
                'gamma',
                'beta',
                'y_pred',
                'High_sd_base',
                'Low_sd_base',
                'High_sd_pos',
                'Low_sd_pos',
                'High_sd_neg',
                'Low_sd_neg',
                'High_sd_both',
                'Low_sd_both'
)

fit.lv2 <- stan(file = paste0(stan.folder, 'STAN_lv2b.stan'),
                data = dat.complete, 
                warmup = 1000,
                iter = 2000, 
                chains = 4, 
                refresh=20,
                pars = pars.monitor,
                control = list(adapt_delta=.95)
)