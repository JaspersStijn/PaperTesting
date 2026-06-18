## Load packages

library(tidyverse)
library(gee)
library(pim)
library(matrixcalc)
#install.packages("geesmv")
library(geesmv)

## Functions

expit = function(x) return(exp(x)/(1+exp(x)))


## Data Diacerein on Blister Pain

dat = read.csv(file="Diacerein_blister_pain_pruritus.txt",sep="\t")


# Add treatment period that indicates whether the subjects got the treatment before (period = 1) or 
# after (period=2) the placebo. This might have an impact on the blister count after the wash-out time

dat = dat %>% mutate(Period= ifelse(((Group=="V")&(Time%in%c("t0","t2","t4","t7")))|((Group=="P")&!(Time%in%c("t0","t2","t4","t7"))),1,2))

# We select the patients that have 8 measurement times in total
use_ID = dat%>%filter(Completed.years == "Placebo + Verum" & Blister_count!= "n/a")%>%group_by(Id)%>%summarise(n=n())%>%filter(n==8)%>%dplyr::select(Id)


# Get the longitudinal measurements per ID, per treatment period. The measurements before washout are  
# referred to as Baseline, W2, W4 and FU. The measurements after washout are referred to as Baseline_CO,
# W2_CO, W4_CO and FU_CO. If the treatment period=1, the CO measurements are obtained under placebo, 
# if trt period = 2,  these CO measurements are obtained under treatment.

Period1_Treated = dat%>%dplyr::select(c(Id,Group,Period,Time,Blister_count))  %>% group_by(Id,Group) %>%
  pivot_wider(names_from=Time, values_from = Blister_count, names_prefix = "val")%>%
  filter(Id %in% t(use_ID)&Group=="V"&!is.na(valt0))%>%dplyr::select(Group,Period,Id,valt0,valt2,
                                                              valt4,valt7)%>%
  mutate(Baseline=valt0,W2=valt2,W4=valt4,FU=valt7)

Period1_Treated_CO = dat%>%dplyr::select(c(Id,Group,Period,Time,Blister_count))  %>% group_by(Id,Group) %>%
  pivot_wider(names_from=Time, values_from = Blister_count, names_prefix = "val")%>%
  filter(Id %in% t(use_ID)&Group=="P"&is.na(valt0))%>%dplyr::select(Group,Period,Id,valt8,valt10,
                                                              valt12,valt15)%>%
  mutate(Baseline_CO=valt8,W2_CO=valt10,W4_CO=valt12,FU_CO=valt15)


Period1_Placebo = dat%>%dplyr::select(c(Id,Group,Period,Time,Blister_count))  %>% group_by(Id,Group) %>%
  pivot_wider(names_from=Time, values_from = Blister_count, names_prefix = "val")%>%
  filter(Id %in% t(use_ID)&Group=="P"&!is.na(valt0))%>%dplyr::select(Group,Period,Id,valt0,valt2,
                                                              valt4,valt7)%>%
  mutate(Baseline=valt0,W2=valt2,W4=valt4,FU=valt7)


Period1_Placebo_CO = dat%>%dplyr::select(c(Id,Group,Period,Time,Blister_count))  %>% group_by(Id,Group) %>%
  pivot_wider(names_from=Time, values_from = Blister_count, names_prefix = "val")%>%
  filter(Id %in% t(use_ID)&Group=="V"&is.na(valt0))%>%dplyr::select(Group,Period,Id,valt8,valt10,
                                                             valt12,valt15)%>%
  mutate(Baseline_CO=valt8,W2_CO=valt10,W4_CO=valt12,FU_CO=valt15)


# Put them together in wide format

dat_wide = Period1_Treated%>%ungroup()%>%dplyr::select(Id,Period,Baseline, W2, W4, FU)%>%
  left_join(Period1_Treated_CO%>%ungroup()%>%dplyr::select(Id,Period,Baseline_CO, W2_CO, W4_CO, FU_CO))%>%
  bind_rows(Period1_Placebo%>%ungroup()%>%dplyr::select(Id,Period,Baseline, W2, W4, FU)%>%
              left_join(Period1_Placebo_CO%>%ungroup()%>%dplyr::select(Id,Period,Baseline_CO, W2_CO, W4_CO, FU_CO)))

dat_wide = dat_wide %>% mutate_if(is.character, as.numeric)

# Create the binary indicator which equals 1 if the percentage reduction in blister counts is
# more than 40%
# Also create a continuous variable (blister count relative vs. baseline)

dat_wide = dat_wide%>% mutate(W2_bin = (((Baseline-W2)/Baseline)>0.40)*1,
                              W4_bin = (((Baseline-W4)/Baseline)>0.40)*1,
                              FU_bin = (((Baseline-FU)/Baseline)>0.40)*1,
                              W2_CO_bin = ifelse(Baseline_CO==0,0,(((Baseline_CO-W2_CO)/Baseline_CO)>0.40)*1),
                              W4_CO_bin = ifelse(Baseline_CO==0,0,(((Baseline_CO-W4_CO)/Baseline_CO)>0.40)*1),
                              FU_CO_bin = ifelse(Baseline_CO==0,0,(((Baseline_CO-FU_CO)/Baseline_CO)>0.40)*1),
                              W2_rel =   W2/Baseline,
                              W4_rel =   W4/Baseline,
                              FU_rel =   FU/Baseline,
                              W2_CO_rel = ifelse(Baseline_CO==0,999,W2_CO/Baseline_CO),
                              W4_CO_rel = ifelse(Baseline_CO==0,999,W4_CO/Baseline_CO),
                              FU_CO_rel = ifelse(Baseline_CO==0,999,FU_CO/Baseline_CO)
                              )


## Patient 2006 received treatment in period 1 and
#  has 0 blister count after washout period, so at start of CO the count is 0 ==> drop?                                                                                                                                                                                                                            


## Compare measurements within patient

# 1. Pairwise comparison of timepoints (W2 treated with W2 Placebo) + Period covariate

# Pseudo-observations 0, 0.5 or 1 for loss, tie, win when comparing the placebo vs the treatment: 
# Y_bin_pl <= Y_bin_trt. Note that when Period = 1, the placebo is in CO vraiables

dat_wide = dat_wide%>%mutate(W2_PO = ifelse(Period==1, (W2_CO_bin < W2_bin)*1 + 0.5*(W2_CO_bin == W2_bin), (W2_bin < W2_CO_bin)*1 + 0.5*(W2_CO_bin == W2_bin)),
                             W4_PO = ifelse(Period==1, (W4_CO_bin < W4_bin)*1 + 0.5*(W4_CO_bin == W4_bin), (W4_bin < W4_CO_bin)*1 + 0.5*(W4_CO_bin == W4_bin)),
                             FU_PO = ifelse(Period==1, (FU_CO_bin < FU_bin)*1 + 0.5*(FU_CO_bin == FU_bin), (FU_bin < FU_CO_bin)*1 + 0.5*(FU_CO_bin == FU_bin)))

dat_long = dat_wide%>%dplyr::select(c(Id,Period,W2_PO,W4_PO,FU_PO)) %>%pivot_longer(!c(Id,Period),names_to="Time",values_to="PO")%>%mutate(Period=as.factor(Period))

# number of pseudo-observations: 13 subjects * 3 comparisons  = 39 POs
nrow(dat_long)

mod = glm(PO~1+Period,data=dat_long,family=binomial(link="logit"))
summary(mod)
expit(coef(mod)[1]) # A probability of 0.8095238 that at a specific timepoint, the binary indicator for the placebo outcome within a patient is smaller or equal than the outcome under treatment, when the treatment was given in period 1     
expit(coef(mod)[1]+coef(mod)[2]) # A probability of 0.5 that at a specific timepoint, the binary indicator for the placebo outcome within a patient is smaller or equal than the outcome under treatment, when the treatment was given in period 2 

mod_gee = gee(PO~1+Period,id=Id,data=dat_long,family=binomial(link="logit"))
summary(mod_gee)
expit(coef(mod_gee)[1])
expit(coef(mod_gee)[1]+coef(mod_gee)[2])


dat_long_PIM_1_1 = dat_wide%>%filter(Period==1)%>%dplyr::select(c(Id,Period,W2_CO_bin,W4_CO_bin,FU_CO_bin)) %>%pivot_longer(!c(Id,Period),names_to="Time",values_to="Y")
dat_long_PIM_1_2 = dat_wide%>%filter(Period==2)%>%dplyr::select(c(Id,Period,W2_bin,W4_bin,FU_bin)) %>%pivot_longer(!c(Id,Period),names_to="Time",values_to="Y")

dat_long_PIM_2_1 = dat_wide%>%filter(Period==1)%>%dplyr::select(c(Id,Period,W2_bin,W4_bin,FU_bin)) %>%pivot_longer(!c(Id,Period),names_to="Time",values_to="Y")
dat_long_PIM_2_2 = dat_wide%>%filter(Period==2)%>%dplyr::select(c(Id,Period,W2_CO_bin,W4_CO_bin,FU_CO_bin)) %>%pivot_longer(!c(Id,Period),names_to="Time",values_to="Y")%>%mutate(Period=Period+1)

PIM_data = dat_long_PIM_1_1%>%bind_rows(dat_long_PIM_1_2)%>%bind_rows(dat_long_PIM_2_1)%>%bind_rows(dat_long_PIM_2_2)

compare <- data.frame("Var1" = 1:(nrow(dat_long_PIM_1_1)+nrow(dat_long_PIM_1_2)),"Var2" = ((nrow(dat_long_PIM_1_1)+nrow(dat_long_PIM_1_2)))+( 1:(nrow(dat_long_PIM_2_1)+nrow(dat_long_PIM_2_2))))

mod_PIM = pim(Y~1+Period,compare = compare,data=PIM_data)
vcov(mod_PIM)
summary(mod_gee)
dat_long = dat_long%>%arrange(Id)
dat_long$id = rep(1:13,each=3)
gee_mbn <- GEE.var.mbn_new(PO~1+Period,id="id",data=dat_long,family=binomial,corstr="independence",d=2,r=1) ##Independence correlation structure;

mod_gee$robust.variance
gee_mbn$cov.beta.full


mod = gee_supressed(PO~1+Period,id=id,data=dat_long,family=binomial(link="logit"), 
              corstr = "independence")

mod$robust.variance
# 2. Pairwise comparison of timepoints (W2 treated with W2 Placebo) + Period and time difference as covariate

PIM_data = dat_long_PIM_1_1%>%bind_rows(dat_long_PIM_1_2)%>%bind_rows(dat_long_PIM_2_1)%>%bind_rows(dat_long_PIM_2_2)%>%mutate(Timing=rep(c(1,2,3),26))
compare <- data.frame("Var1" = rep(c(1:(nrow(dat_long_PIM_1_1)+nrow(dat_long_PIM_1_2))),rep(c(3,2,1),(nrow(dat_long_PIM_1_1)+nrow(dat_long_PIM_1_2))/3)),"Var2" = rep(seq((nrow(dat_long_PIM_1_1)+nrow(dat_long_PIM_1_2))+1,(nrow(dat_long_PIM_1_1)+nrow(dat_long_PIM_1_2))+(nrow(dat_long_PIM_2_1)+nrow(dat_long_PIM_2_2)),3),each=6)+rep(c(0,1,2,1,2,2),(nrow(dat_long_PIM_2_1)+nrow(dat_long_PIM_2_2))/3))
mod_PIM2 = pim(Y~1+Period+Timing,compare = compare,data=PIM_data)
summary(mod_PIM2)
vcov(mod_PIM2)

mod_PIM2b = pim(Y~1+Period+as.factor(Timing),compare = compare,data=PIM_data)
summary(mod_PIM2b)
vcov(mod_PIM2b)

PIM_data$Timing

expit(mod_PIM2b@fitted)


# Number of POs: 78 = 13 subjects with 6 comparisons (W2 vs. W2_CO, W4_CO, FU_CO ; W4 vs. W4_CO, FU_CO; FU vs. FU_CO)
length(response(mod_PIM2))

dat_GEE2 = data.frame(cbind(response(mod_PIM2),model.matrix(mod_PIM2)[,-1]))
dat_GEE2$Id = rep(1:(nrow(dat_GEE2)/6),each=6) ## cluster per subject
mod_gee2 = gee(V1~1+Period+Timing,id=Id,data=dat_GEE2,family=binomial(link="logit"))
summary(mod_gee2)
mod_gee2$robust.variance


gee_mbn <- GEE.var.mbn_new(V1~1+Period+Timing,id="Id",data=dat_GEE2,family=binomial,corstr="independence",d=2,r=1) ##Independence correlation structure;
gee_mbn$cov.beta           
gee_mbn$cov.beta.full           


model.matrix(mod_PIM2b)


dat_GEE2b = data.frame(cbind(response(mod_PIM2b),model.matrix(mod_PIM2b)[,-1]))
dat_GEE2b$Id = rep(1:(nrow(dat_GEE2b)/6),each=6) ## cluster per subject
mod_gee2b = gee(V1~1+Period+as.factor.Timing.2+as.factor.Timing.3,id=Id,data=dat_GEE2b,family=binomial(link="logit"))
summary(mod_gee2b)
mod_gee2b$robust.variance

summary(mod_PIM2b)

## Check whether GEE models give similar se for larger datasets
dat_GEE_test = rbind(dat_GEE2,dat_GEE2[sample(1:nrow(dat_GEE2),nrow(dat_GEE2),replace=TRUE),],dat_GEE2[sample(1:nrow(dat_GEE2),nrow(dat_GEE2),replace=TRUE),],dat_GEE2[sample(1:nrow(dat_GEE2),nrow(dat_GEE2),replace=TRUE),])
dat_GEE_test$Id = rep(1:(nrow(dat_GEE_test)/6),each=6) ## cluster per subject
mod_gee_test = gee(V1~1+Period+Timing,id=Id,data=dat_GEE_test,family=binomial(link="logit"))
sqrt(diag(mod_gee_test$robust.variance))
gee_mbn_test <- GEE.var.mbn_new(V1~1+Period+Timing,id="Id",data=dat_GEE_test,family=binomial,corstr="independence",d=2,r=1) ##Independence correlation structure;
sqrt(diag(gee_mbn_test$cov.beta.full))           



#3. Compare two treatments ignoring the crossover, summarising the binary indiator as sum

dat_trt_wide = dat_wide%>%mutate(O1 = W2_bin+W4_bin+FU_bin,
                  O2 = W2_CO_bin+W4_CO_bin+FU_CO_bin)%>%dplyr::select(Id,Period,O1,O2)

dat_long1_1 = dat_trt_wide%>%filter(Period==1)%>%dplyr::select(Id, Period, O1)%>%mutate(Treatment="TRT",Y=O1)%>%dplyr::select(Id, Period,Y,Treatment)
dat_long1_2 = dat_trt_wide%>%filter(Period==1)%>%dplyr::select(Id, Period, O2)%>%mutate(Treatment="CTRL",Y=O2)%>%dplyr::select(Id, Period,Y,Treatment)

dat_long2_1 = dat_trt_wide%>%filter(Period==2)%>%dplyr::select(Id, Period, O2)%>%mutate(Treatment="TRT",Y=O2)%>%dplyr::select(Id, Period,Y,Treatment)
dat_long2_2 = dat_trt_wide%>%filter(Period==2)%>%dplyr::select(Id, Period, O1)%>%mutate(Treatment="CTRL",Y=O1)%>%dplyr::select(Id, Period,Y,Treatment)


dat_long_trt = dat_long1_1 %>%bind_rows(dat_long1_2)%>%bind_rows(dat_long2_1) %>% bind_rows(dat_long2_2)

# We have two groups of 13 observations, so we use the small sample PIM correction based on GEE
# We control for period

fit = GEE_MH_fit(data=data.frame(dat_long_trt) ,
                 response="Y",
                 treatment = "Treatment",
                 control = "Period",
                 correction = "MBN")

# The PI is 
expit(fit$beta_est)

# In terms of the NTB: 2*PI - 1
2*expit(fit$beta_est) - 1 ## -0.5503787 cfr. also paper by Johan Table 3, unmatched univariate GPC for the dichotomized outcome: 0.55



id.nonfac <- which(dat_long_trt[,"Treatment"] == "TRT")
id.fac <- which(dat_long_trt[,"Treatment"] == "CTRL")

compare <- expand.grid(id.nonfac, id.fac)


mod = pim(Y~1+Period,data=dat_long_trt,compare=compare)
2*(expit(coef(mod)[1]+coef(mod)[2]))-1
2*(expit(coef(mod)[1]))-1
2*(expit(coef(mod)[1]-coef(mod)[2]))-1





## Compare continuous endpoints (blister count relative vs. baseline)

dat_long_PIM_1_1 = dat_wide%>%filter(Period==1)%>%dplyr::select(c(Id,Period,W2_CO_rel,W4_CO_rel,FU_CO_rel)) %>%pivot_longer(!c(Id,Period),names_to="Time",values_to="Y")%>%filter(Id!=2006)
dat_long_PIM_1_2 = dat_wide%>%filter(Period==2)%>%dplyr::select(c(Id,Period,W2_rel,W4_rel,FU_rel)) %>%pivot_longer(!c(Id,Period),names_to="Time",values_to="Y")%>%filter(Id!=2006)

dat_long_PIM_2_1 = dat_wide%>%filter(Period==1)%>%dplyr::select(c(Id,Period,W2_rel,W4_rel,FU_rel)) %>%pivot_longer(!c(Id,Period),names_to="Time",values_to="Y")%>%filter(Id!=2006)
dat_long_PIM_2_2 = dat_wide%>%filter(Period==2)%>%dplyr::select(c(Id,Period,W2_CO_rel,W4_CO_rel,FU_CO_rel)) %>%pivot_longer(!c(Id,Period),names_to="Time",values_to="Y")%>%mutate(Period=Period+1)%>%filter(Id!=2006)

PIM_data_rel = dat_long_PIM_1_1%>%bind_rows(dat_long_PIM_1_2)%>%bind_rows(dat_long_PIM_2_1)%>%bind_rows(dat_long_PIM_2_2)%>%mutate(Timing=rep(c(1,2,3),24))

compare <- data.frame("Var1" = rep(c(1:(nrow(dat_long_PIM_1_1)+nrow(dat_long_PIM_1_2))),rep(c(3,2,1),(nrow(dat_long_PIM_1_1)+nrow(dat_long_PIM_1_2))/3)),"Var2" = rep(seq((nrow(dat_long_PIM_1_1)+nrow(dat_long_PIM_1_2))+1,(nrow(dat_long_PIM_1_1)+nrow(dat_long_PIM_1_2))+(nrow(dat_long_PIM_2_1)+nrow(dat_long_PIM_2_2)),3),each=6)+rep(c(0,1,2,1,2,2),(nrow(dat_long_PIM_2_1)+nrow(dat_long_PIM_2_2))/3))
mod_PIM2 = pim(Y~1+Period+Timing,compare = compare,data=PIM_data_rel)
summary(mod_PIM2)
vcov(mod_PIM2)

mod_PIM2b = pim(Y~1+Period+as.factor(Timing),compare = compare,data=PIM_data_rel)
summary(mod_PIM2b)
vcov(mod_PIM2b)

expit(-1.6987)
model.matrix(mod_PIM2b)






dat_long_PIM_1_1 = dat_wide%>%filter(Period==1)%>%dplyr::select(c(Id,Period,Baseline_CO, W2_CO, W4_CO, FU_CO)) %>%pivot_longer(!c(Id,Period),names_to="Time",values_to="Y")%>%filter(Id!=2006)
dat_long_PIM_1_2 = dat_wide%>%filter(Period==2)%>%dplyr::select(c(Id,Period,Baseline ,   W2 ,   W4 ,   FU )) %>%pivot_longer(!c(Id,Period),names_to="Time",values_to="Y")%>%filter(Id!=2006)

dat_long_PIM_2_1 = dat_wide%>%filter(Period==1)%>%dplyr::select(c(Id,Period,Baseline ,   W2 ,   W4 ,   FU )) %>%pivot_longer(!c(Id,Period),names_to="Time",values_to="Y")%>%filter(Id!=2006)
dat_long_PIM_2_2 = dat_wide%>%filter(Period==2)%>%dplyr::select(c(Id,Period,Baseline_CO, W2_CO, W4_CO, FU_CO)) %>%pivot_longer(!c(Id,Period),names_to="Time",values_to="Y")%>%filter(Id!=2006)


ctrl = dat_long_PIM_1_1%>%bind_rows(dat_long_PIM_1_2)
ctrl$trt = "Control"
ctrl$Timing = rep(c(1,2,3,4),nrow(ctrl)/4)
trt = dat_long_PIM_2_1%>%bind_rows(dat_long_PIM_2_2)
trt$trt = "Treatment"
trt$Timing = rep(c(1,2,3,4),nrow(trt)/4)


check = ctrl%>%bind_rows(trt)

p=ggplot(check,aes(x=Timing,y=Y,linetype=as.factor(trt),colour=as.factor(Id)))+geom_line()
p + facet_wrap(vars(Period))

summar = check%>%group_by(trt,Timing,Period)%>%summarise(mean = mean(Y))
summar$Id = rep(c("Med_ctrl","Med_trt"),each=8)
p+ facet_wrap(vars(Period))+
  geom_line(data=summar, aes(x=Timing, y=mean,linetype=as.factor(trt)), colour="grey50")

ggplot(data=summar, aes(x=Timing, y=mean,linetype=as.factor(trt)), colour="grey50")+geom_line()+ facet_wrap(vars(Period))





dat_long_PIM_1_1 = dat_wide%>%filter(Period==1)%>%dplyr::select(c(Id,Period,W2_CO_bin,W4_CO_bin,FU_CO_bin)) %>%pivot_longer(!c(Id,Period),names_to="Time",values_to="Y")%>%filter(Id!=2006)
dat_long_PIM_1_2 = dat_wide%>%filter(Period==2)%>%dplyr::select(c(Id,Period,W2_bin,W4_bin,FU_bin)) %>%pivot_longer(!c(Id,Period),names_to="Time",values_to="Y")%>%filter(Id!=2006)

dat_long_PIM_2_1 = dat_wide%>%filter(Period==1)%>%dplyr::select(c(Id,Period,W2_bin,W4_bin,FU_bin)) %>%pivot_longer(!c(Id,Period),names_to="Time",values_to="Y")%>%filter(Id!=2006)
dat_long_PIM_2_2 = dat_wide%>%filter(Period==2)%>%dplyr::select(c(Id,Period,W2_CO_bin,W4_CO_bin,FU_CO_bin)) %>%pivot_longer(!c(Id,Period),names_to="Time",values_to="Y")%>%filter(Id!=2006)

ctrl = dat_long_PIM_1_1%>%bind_rows(dat_long_PIM_1_2)
ctrl$trt = "Control"
ctrl$Timing = rep(c(1,2,3),nrow(ctrl)/3)
trt = dat_long_PIM_2_1%>%bind_rows(dat_long_PIM_2_2)
trt$trt = "Treatment"
trt$Timing = rep(c(1,2,3),nrow(trt)/3)


check = ctrl%>%bind_rows(trt)

p=ggplot(check,aes(x=Timing,y=Y,linetype=as.factor(trt),colour=as.factor(Id)))+geom_line()
p + facet_wrap(vars(Period))

summar = check%>%group_by(trt,Timing,Period)%>%summarise(mean = mean(Y))
summar$Id = rep(c("Med_ctrl","Med_trt"),each=6)
p+ facet_wrap(vars(Period))+
  geom_line(data=summar, aes(x=Timing, y=mean,linetype=as.factor(trt)), colour="grey50")

ggplot(data=summar, aes(x=Timing, y=mean,linetype=as.factor(trt)), colour="grey50")+geom_line()+ facet_wrap(vars(Period))
  




dat_long_PIM_1_1 = dat_wide%>%filter(Period==1)%>%dplyr::select(c(Id,Period,W2_CO_rel,W4_CO_rel,FU_CO_rel)) %>%pivot_longer(!c(Id,Period),names_to="Time",values_to="Y")%>%filter(Id!=2006)
dat_long_PIM_1_2 = dat_wide%>%filter(Period==2)%>%dplyr::select(c(Id,Period,W2_rel,W4_rel,FU_rel)) %>%pivot_longer(!c(Id,Period),names_to="Time",values_to="Y")%>%filter(Id!=2006)

dat_long_PIM_2_1 = dat_wide%>%filter(Period==1)%>%dplyr::select(c(Id,Period,W2_rel,W4_rel,FU_rel)) %>%pivot_longer(!c(Id,Period),names_to="Time",values_to="Y")%>%filter(Id!=2006)
dat_long_PIM_2_2 = dat_wide%>%filter(Period==2)%>%dplyr::select(c(Id,Period,W2_CO_rel,W4_CO_rel,FU_CO_rel)) %>%pivot_longer(!c(Id,Period),names_to="Time",values_to="Y")%>%filter(Id!=2006)

ctrl = dat_long_PIM_1_1%>%bind_rows(dat_long_PIM_1_2)
ctrl$trt = "Control"
ctrl$Timing = rep(c(1,2,3),nrow(ctrl)/3)
trt = dat_long_PIM_2_1%>%bind_rows(dat_long_PIM_2_2)
trt$trt = "Treatment"
trt$Timing = rep(c(1,2,3),nrow(trt)/3)


check = ctrl%>%bind_rows(trt)

p=ggplot(check,aes(x=Timing,y=Y,linetype=as.factor(trt),colour=as.factor(Id)))+geom_line()
p + facet_wrap(vars(Period))

summar = check%>%group_by(trt,Timing,Period)%>%summarise(median = median(Y))
summar$Id = rep(c("Med_ctrl","Med_trt"),each=6)
p+ facet_wrap(vars(Period))+
  geom_line(data=summar, aes(x=Timing, y=median,linetype=as.factor(trt)), colour="grey50")

ggplot(data=summar, aes(x=Timing, y=median,linetype=as.factor(trt)), colour="grey50")+geom_line()+ facet_wrap(vars(Period))



g1 = rep(1:5,each=2)
g2 = 6:10

length(union(g1,g2))
(xsum <- rowsum(x, group))
crossprod(xsum )

t(xsum)%*%xsum



dat_plot = subset(dat,Time%in%paste0("t",c(0,2,4,7)))
ggplot(dat_plot,aes(x=Time,y=Blister_count,color=Group,group=Id))+geom_line()

dat_plot$group=dat_plot$Group
dat_plot=dat_plot%>%mutate(time = case_when(Time=="t0" ~1,Time=="t2" ~2,Time=="t4" ~3,Time=="t7" ~4))



dat_wide = Period1_Treated%>%ungroup()%>%dplyr::select(Group,Id,Period,Baseline, W2, W4, FU)%>%
  left_join(Period1_Treated_CO%>%ungroup()%>%dplyr::select(Id,Period,Baseline_CO, W2_CO, W4_CO, FU_CO))%>%
  bind_rows(Period1_Placebo%>%ungroup()%>%dplyr::select(Group,Id,Period,Baseline, W2, W4, FU)%>%
              left_join(Period1_Placebo_CO%>%ungroup()%>%dplyr::select(Id,Period,Baseline_CO, W2_CO, W4_CO, FU_CO)))


dat_wide[,2:ncol(dat_wide)] = dat_wide[,2:ncol(dat_wide)] %>% mutate_if(is.character, as.numeric)


dat_wide = dat_wide%>% mutate(W2_bin = (((Baseline-W2)/Baseline)>0.40)*1,
                              W4_bin = (((Baseline-W4)/Baseline)>0.40)*1,
                              FU_bin = (((Baseline-FU)/Baseline)>0.40)*1,
                              W2_red_rel =   (Baseline-W2)/Baseline,
                              W4_red_rel =   (Baseline-W4)/Baseline,
                              FU_red_rel =   (Baseline-FU)/Baseline,

)



dat_long = dat_wide%>%dplyr::select(c(Id,Period,W2_red_rel,W4_red_rel,FU_red_rel,Group)) %>%pivot_longer(!c(Id,Period,Group),names_to="Time",values_to="Red")%>%mutate(Period=as.factor(Period))
dat_long=dat_long%>%mutate(time = case_when(Time=="W2_red_rel" ~1,Time=="W4_red_rel" ~2,Time=="FU_red_rel" ~3))
dat_long$group = dat_long$Group

unique(dat_long$Id)
dat_long=dat_long%>%mutate(id= case_when(Id=="1001" ~1,Id=="1004" ~2,Id=="2002" ~3,
                                         Id=="2004" ~4,Id=="2005" ~5,Id=="2006" ~6,
                                         Id=="2008" ~7,Id=="1002" ~8,Id=="1005" ~9,
                                         Id=="2001" ~10,Id=="2003" ~11,Id=="2007" ~12,
                                         Id=="2009" ~13))

ggplot(dat_long,aes(x=time,y=Red,linetype=as.factor(group),group=as.factor(Id)))+geom_line()


id.fac1 <- which((dat_long[,"group"] == "P")&(dat_long[,"time"] == 1))
id.nonfac1 <- which((dat_long[,"group"] == "V")&(dat_long[,"time"] == 1))
compare1 <- expand.grid(id.fac1,id.nonfac1)

id.fac2 <- which((dat_long[,"group"] == "P")&(dat_long[,"time"] == 2))
id.nonfac2 <- which((dat_long[,"group"] == "V")&(dat_long[,"time"] == 2))
compare2 <- expand.grid(id.fac2,id.nonfac2)

id.fac3 <- which((dat_long[,"group"] == "P")&(dat_long[,"time"] == 3))
id.nonfac3 <- which((dat_long[,"group"] == "V")&(dat_long[,"time"] == 3))
compare3 <- expand.grid(id.fac3,id.nonfac3)


compare_between = rbind(compare1,compare2,compare3)
compare_between = compare_between[order(compare_between[,"Var2"]),]

start_timepoints <- c()
later_timepoints <- c()
k=1
for (i in 1:(nrow(dat_long) - 1)) {
  if(i%%3 != 0){
    start <- i           # Current time point
    later <- (i+1) : (3*k)  # Later time points
    # Store the results in vectors
    start_timepoints <- c(start_timepoints, rep(start, length(later)))
    later_timepoints <- c(later_timepoints, later)}
  if(i%%3 == 0){
    k=k+1}
}

within1 = start_timepoints
within2 = later_timepoints

compare_within = cbind(within1,within2)
colnames(compare_within) = c("Var1","Var2")
compare = rbind(compare_between,compare_within)
str(dat_long)
dat_long$group=ifelse(dat_long$group=="V",1,0)
library(pim)
assignInNamespace("sandwich.estimator", sandwich.estimator, ns = "pim")
mod2_orig = pim(Red~I((R(time) - L(time))*R(group)*L(group))+I((R(time) - L(time))*(1-R(group))*(1-L(group)))+I((R(group)-L(group))*(R(time)==1)*(L(time)==1))+I((R(group)-L(group))*(R(time)==2)*(L(time)==2))+I((R(group)-L(group))*(R(time)==3)*(L(time)==3)) ,data=dat_long,compare=compare,link="probit")
summary(mod2_orig)



individuals_var1=dat_long[compare$Var1,"id"]%>%unlist()
individuals_var2=dat_long[compare$Var2,"id"]%>%unlist()

assignInNamespace("sandwich.estimator", sandwich.estimator_clustered_mod1_mod2, ns = "pim")
mod2_adj = pim(Red~I((R(time) - L(time))*R(group)*L(group))+I((R(time) - L(time))*(1-R(group))*(1-L(group)))+I((R(group)-L(group))*(R(time)==1)*(L(time)==1))+I((R(group)-L(group))*(R(time)==2)*(L(time)==2))+I((R(group)-L(group))*(R(time)==3)*(L(time)==3)) ,data=dat_long,compare=compare,link="probit")
summary(mod2_adj)


pnorm(b[3:5])
V=vcov(mod2_adj)
b = coef(mod2_adj)

L_new <- rbind(
  c(0,0,1,-1,0),
  c(0,0,1,0,-1)
)

library(aod)

### Between-Within degrees of freedom
out = wald.test(V, b, Terms = NULL, L = L_new, H0 = NULL,  
                         df = 13-qr(model.matrix(mod2_adj))$rank, verbose = FALSE)

out$result$Ftest

qr(model.matrix(mod2_adj)[which(individuals_var1==8),])$rank


compare = rbind(compare_between,compare_within)

cbind(compare,individuals_var1,individuals_var2)

model.matrix(mod2_orig)


mod = lm(Red~time*group,data=dat_long)
summary(mod)
dat_sim = dat_long
beta1 = 0.22616
beta2= 0.27058
beta3 = 0.02710

dat_sim$Red_use = -0.13+beta1*dat_long$time+beta2*dat_long$group+beta3*dat_long$time*dat_long$group+rep(rnorm(13,0,0.05),each=nrow(dat_long)/13)+rnorm(nrow(dat_sim),0,0.3044)

dat_sim$time


dat_test = data.frame(time = rep(c(1:3),50),group = rep(c(1,0),each=75))
dat_test$Red_use = -0.13+(beta1+rep(rnorm(50,0,0.15),each=nrow(dat_test)/50))*dat_test$time+beta2*dat_test$group+beta3*dat_test$time*dat_test$group+rep(rnorm(50,0,0.15),each=nrow(dat_test)/50)+rnorm(nrow(dat_test),0,0.3044)
dat_test$id=rep(1:50,each=3)
id.fac1 <- which((dat_test[,"group"] == "0")&(dat_test[,"time"] == 1))
id.nonfac1 <- which((dat_test[,"group"] == "1")&(dat_test[,"time"] == 1))
compare1 <- expand.grid(id.fac1,id.nonfac1)

id.fac2 <- which((dat_test[,"group"] == "0")&(dat_test[,"time"] == 2))
id.nonfac2 <- which((dat_test[,"group"] == "1")&(dat_test[,"time"] == 2))
compare2 <- expand.grid(id.fac2,id.nonfac2)

id.fac3 <- which((dat_test[,"group"] == "0")&(dat_test[,"time"] == 3))
id.nonfac3 <- which((dat_test[,"group"] == "1")&(dat_test[,"time"] == 3))
compare3 <- expand.grid(id.fac3,id.nonfac3)


compare_between = rbind(compare1,compare2,compare3)
compare_between = compare_between[order(compare_between[,"Var2"]),]

start_timepoints <- c()
later_timepoints <- c()
k=1
for (i in 1:(nrow(dat_test) - 1)) {
  if(i%%3 != 0){
    start <- i           # Current time point
    later <- (i+1) : (3*k)  # Later time points
    # Store the results in vectors
    start_timepoints <- c(start_timepoints, rep(start, length(later)))
    later_timepoints <- c(later_timepoints, later)}
  if(i%%3 == 0){
    k=k+1}
}

within1 = start_timepoints
within2 = later_timepoints

compare_within = cbind(within1,within2)
colnames(compare_within) = c("Var1","Var2")
compare = rbind(compare_between,compare_within)

assignInNamespace("sandwich.estimator", sandwich.estimator, ns = "pim")
mod2_orig_sim = pim(Red_use~I((R(time) - L(time))*R(group)*L(group))+I((R(time) - L(time))*(1-R(group))*(1-L(group)))+I((R(group)-L(group))*(R(time)==1)*(L(time)==1))+I((R(group)-L(group))*(R(time)==2)*(L(time)==2))+I((R(group)-L(group))*(R(time)==3)*(L(time)==3)) ,data=dat_test,compare=compare,link="probit")
summary(mod2_orig_sim)

check = glm(response(mod2_orig_sim)~-1+model.matrix(mod2_orig_sim),family = binomial(link="probit"))
summary(check)

individuals_var1=dat_test[compare$Var1,"id"]%>%unlist()
individuals_var2=dat_test[compare$Var2,"id"]%>%unlist()

assignInNamespace("sandwich.estimator", sandwich.estimator_clustered_mod1_mod2, ns = "pim")
mod2_adj_sim = pim(Red_use~I((R(time) - L(time))*R(group)*L(group))+I((R(time) - L(time))*(1-R(group))*(1-L(group)))+I((R(group)-L(group))*(R(time)==1)*(L(time)==1))+I((R(group)-L(group))*(R(time)==2)*(L(time)==2))+I((R(group)-L(group))*(R(time)==3)*(L(time)==3)) ,data=dat_test,compare=compare,link="probit")
summary(mod2_adj_sim)
v_check = vcov(check)
rownames(v_check) = colnames(v_check)=NULL
v_pim = vcov(mod2_orig_sim)
rownames(v_pim) = colnames(v_pim)=NULL


sqrt(diag(vcov(mod2_adj_sim)))
sqrt(diag(vcov(mod2_orig_sim)))



tmp.new1 <- PIMM.new(formula = Red ~ time+group,  data = dat_long, id = dat_long$Id,
                     index.function = 1, var = TRUE)

tmp.new2 <- PIMM.new(formula = Red ~ time+group,  data = dat_long, id = dat_long$Id,
                     index.function = 2, var = TRUE)

tmp.new1$beta
tmp.new1$sigma2
tmp.new1$vcov
sqrt(diag(tmp.new1$vcov))





## Review paper Statistics in Biopharmaceutical Research (mail Johan)
library(geessbin)

m0 = c(20,15,9,11)
m1 = c(20,16,12,9)
sigma0 = matrix(c(12.8,7.4,3.6,7.1,7.4,43.2,22.7,23.8,3.6,22.7,21.8,18.8,7.1,23.8,18.8,22.4),nrow=4)
sigma1 = matrix(c(15.6,12.9,4.8,4.4,12.9,37.5,22.8,11.6,4.8,22.8,34.2,17.9,4.4,11.6,17.9,21.9),nrow=4)

dat0 = mvrnorm(50,m0,sigma0)
dat1 = mvrnorm(50,m1,sigma1)

dat_sim = data.frame(rbind(cbind(dat0,1:50,rep(0,50)),cbind(dat1,51:100,rep(1,50))))
dat_sim_long = dat_sim%>%pivot_longer(!c(X5,X6),names_to="Time",values_to="resp")%>%mutate(time=case_when(Time=="X1" ~1,
                                                                                                          Time=="X2" ~2,
                                                                                                          Time=="X3" ~3,
                                                                                                          Time=="X4" ~4),
                                                                                           group=X6)

ggplot(dat_sim_long,aes(x=time,y=resp,linetype=as.factor(group),colour=as.factor(X5)))+geom_line()

id.fac1 <- which((dat_sim_long[,"group"] == "0")&(dat_sim_long[,"time"] == 1))
id.nonfac1 <- which((dat_sim_long[,"group"] == "1")&(dat_sim_long[,"time"] == 1))
compare1 <- expand.grid(id.fac1,id.nonfac1)

id.fac2 <- which((dat_sim_long[,"group"] == "0")&(dat_sim_long[,"time"] == 2))
id.nonfac2 <- which((dat_sim_long[,"group"] == "1")&(dat_sim_long[,"time"] == 2))
compare2 <- expand.grid(id.fac2,id.nonfac2)

id.fac3 <- which((dat_sim_long[,"group"] == "0")&(dat_sim_long[,"time"] == 3))
id.nonfac3 <- which((dat_sim_long[,"group"] == "1")&(dat_sim_long[,"time"] == 3))
compare3 <- expand.grid(id.fac3,id.nonfac3)

id.fac4 <- which((dat_sim_long[,"group"] == "0")&(dat_sim_long[,"time"] == 4))
id.nonfac4 <- which((dat_sim_long[,"group"] == "1")&(dat_sim_long[,"time"] == 4))
compare4 <- expand.grid(id.fac4,id.nonfac4)


compare_between = rbind(compare1,compare2,compare3,compare4)
compare_between = compare_between[order(compare_between[,"Var2"]),]

start_timepoints <- c()
later_timepoints <- c()


for (i in unique(dat_sim_long$X5)) {
  subs_use = subset(dat_sim_long,X5==i)
    starting <- which(dat_sim_long$X5==i)[1]  
    ending <- which(dat_sim_long$X5==i)[length(which(dat_sim_long$X5==i))]# Current time point
    for(k in starting:(ending-1)){
      start = k
    later <- (k+1) : (ending)  # Later time points
    # Store the results in vectors
    start_timepoints <- c(start_timepoints, rep(start, length(later)))
    later_timepoints <- c(later_timepoints, later)}
    }



within1 = start_timepoints
within2 = later_timepoints

compare_within = cbind(within1,within2)
colnames(compare_within) = c("Var1","Var2")
compare = rbind(compare_between,compare_within)

assignInNamespace("sandwich.estimator", sandwich.estimator, ns = "pim")
mod2_orig_sim = pim(resp~I((R(time) - L(time))*R(group)*L(group))+I((R(time) - L(time))*(1-R(group))*(1-L(group)))+I((R(group)-L(group))*(R(time)==1)*(L(time)==1))+I((R(group)-L(group))*(R(time)==2)*(L(time)==2))+I((R(group)-L(group))*(R(time)==3)*(L(time)==3))+I((R(group)-L(group))*(R(time)==4)*(L(time)==4)) ,data=dat_sim_long,compare=compare,link="probit")
summary(mod2_orig_sim)

n_subjects = 10

parms=c()
se_orig = c()
se_adj = c()

se_gee = c()
parms_gee = c()

for(j in 1:1000){
  set.seed(j)
m0 = c(20,15,9,11)
m1 = c(20,16,12,9)
sigma0 = matrix(c(12.8,7.4,3.6,7.1,7.4,43.2,22.7,23.8,3.6,22.7,21.8,18.8,7.1,23.8,18.8,22.4),nrow=4)
sigma1 = matrix(c(15.6,12.9,4.8,4.4,12.9,37.5,22.8,11.6,4.8,22.8,34.2,17.9,4.4,11.6,17.9,21.9),nrow=4)

dat0 = mvrnorm(n_subjects,m0,sigma0)
dat1 = mvrnorm(n_subjects,m1,sigma1)

dat_sim = data.frame(rbind(cbind(dat0,1:n_subjects,rep(0,n_subjects)),cbind(dat1,(n_subjects+1):(2*n_subjects),rep(1,n_subjects))))
dat_sim_long = dat_sim%>%pivot_longer(!c(X5,X6),names_to="Time",values_to="resp")%>%mutate(time=case_when(Time=="X1" ~1,
                                                                                                          Time=="X2" ~2,
                                                                                                          Time=="X3" ~3,
                                                                                                          Time=="X4" ~4),
                                                                                           group=X6)

dat_sim_long$resp = dat_sim_long$resp
# missing = sample(1:100,30,replace = FALSE)
# print(sort(missing))
# missing1=missing[1:10]
# missing2=missing[11:20]
# missing3=missing[21:30]
# 
# dat_sim_long$resp[which((dat_sim_long$X5%in%missing1) & dat_sim_long$time>1)]=NA
# dat_sim_long$resp[which((dat_sim_long$X5%in%missing2) & dat_sim_long$time>2)]=NA
# dat_sim_long$resp[which((dat_sim_long$X5%in%missing3) & dat_sim_long$time>3)]=NA
# dat_sim_long = dat_sim_long%>%drop_na()

id.fac1 <- which((dat_sim_long[,"group"] == "0")&(dat_sim_long[,"time"] == 1))
id.nonfac1 <- which((dat_sim_long[,"group"] == "1")&(dat_sim_long[,"time"] == 1))
compare1 <- expand.grid(id.fac1,id.nonfac1)

id.fac2 <- which((dat_sim_long[,"group"] == "0")&(dat_sim_long[,"time"] == 2))
id.nonfac2 <- which((dat_sim_long[,"group"] == "1")&(dat_sim_long[,"time"] == 2))
compare2 <- expand.grid(id.fac2,id.nonfac2)

id.fac3 <- which((dat_sim_long[,"group"] == "0")&(dat_sim_long[,"time"] == 3))
id.nonfac3 <- which((dat_sim_long[,"group"] == "1")&(dat_sim_long[,"time"] == 3))
compare3 <- expand.grid(id.fac3,id.nonfac3)

id.fac4 <- which((dat_sim_long[,"group"] == "0")&(dat_sim_long[,"time"] == 4))
id.nonfac4 <- which((dat_sim_long[,"group"] == "1")&(dat_sim_long[,"time"] == 4))
compare4 <- expand.grid(id.fac4,id.nonfac4)


compare_between = rbind(compare1,compare2,compare3,compare4)
compare_between = compare_between[order(compare_between[,"Var2"]),]

start_timepoints <- c()
later_timepoints <- c()

for (i in unique(dat_sim_long$X5)) {
  subs_use = subset(dat_sim_long,X5==i)
  if(nrow(subs_use)>1){
  starting <- which(dat_sim_long$X5==i)[1]  
  ending <- which(dat_sim_long$X5==i)[length(which(dat_sim_long$X5==i))]# Current time point
  for(k in starting:(ending-1)){
    start = k
    later <- (k+1) : (ending)  # Later time points
    # Store the results in vectors
    start_timepoints <- c(start_timepoints, rep(start, length(later)))
    later_timepoints <- c(later_timepoints, later)}
}}


within1 = start_timepoints
within2 = later_timepoints

compare_within = cbind(within1,within2)
colnames(compare_within) = c("Var1","Var2")
compare = rbind(compare_between,compare_within)
nrow(compare_within)


assignInNamespace("sandwich.estimator", sandwich.estimator, ns = "pim")
mod2_orig_sim = pim(resp~I((R(time) - L(time))*R(group)*L(group))+I((R(time) - L(time))*(1-R(group))*(1-L(group)))+I((R(group)-L(group))*(R(time)==1)*(L(time)==1))+I((R(group)-L(group))*(R(time)==2)*(L(time)==2))+I((R(group)-L(group))*(R(time)==3)*(L(time)==3))+I((R(group)-L(group))*(R(time)==4)*(L(time)==4)) ,data=dat_sim_long,compare=compare,link="probit")

se_orig = c(se_orig,sqrt(diag(vcov(mod2_orig_sim)))[6])


individuals_var1=dat_sim_long[compare$Var1,"X5"]%>%unlist()
individuals_var2=dat_sim_long[compare$Var2,"X5"]%>%unlist()


assignInNamespace("sandwich.estimator", sandwich.estimator_clustered_mod1_mod2, ns = "pim")
mod2_adj_sim = pim(resp~I((R(time) - L(time))*R(group)*L(group))+I((R(time) - L(time))*(1-R(group))*(1-L(group)))+I((R(group)-L(group))*(R(time)==1)*(L(time)==1))+I((R(group)-L(group))*(R(time)==2)*(L(time)==2))+I((R(group)-L(group))*(R(time)==3)*(L(time)==3))+I((R(group)-L(group))*(R(time)==4)*(L(time)==4)) ,data=dat_sim_long,compare=compare,link="probit")
summary(mod2_adj_sim)
parms = c(parms,pnorm(-coef(mod2_adj_sim)[6]))
se_adj = c(se_adj,sqrt(diag(vcov(mod2_adj_sim)))[6])

## GEE

X = data.frame(model.matrix(mod2_adj_sim))
Y = response(mod2_adj_sim) 
X$id1 = individuals_var1
X$id2 = individuals_var2
X$id3 <- paste(X$id1, X$id2,sep="_")
dat_GEE = cbind(X,Y)

dat_GEE=dat_GEE[order(dat_GEE$id1),]
mod1 = geessbin(Y~.-1-id1-id2-id3,data=dat_GEE,id=id1, corstr = "independence",beta.method="PGEE",SE.method = "FW")
dat_GEE=dat_GEE[order(dat_GEE$id2),]
mod2 = geessbin(Y~.-1-id1-id2-id3,data=dat_GEE,id=id2, corstr = "independence",beta.method="PGEE",SE.method = "FW")
dat_GEE=dat_GEE[order(dat_GEE$id3),]
mod3 =geessbin(Y~.-1-id1-id2-id3,data=dat_GEE,id=id3, corstr = "independence",beta.method="PGEE",SE.method = "FW")

se_gee = c(se_gee ,sqrt(diag(mod1$covb+mod2$covb-mod3$covb)[6]))
parms_gee = c(parms_gee,coef(mod1)[6])

}

true = pnorm((11-9)/sqrt(21.9+22.4))

mean((parms-true)/true)*100
mean((plogis(-parms_gee)-true)/true)*100


lower = pnorm(qnorm(parms)-1.96*se_adj)
upper = pnorm(qnorm(parms)+1.96*se_adj)

low_ind = true>=lower;mean(1-low_ind)
up_ind = true<=upper;mean(1-up_ind)
cov_adj = mean(low_ind*up_ind)
cov_adj
lower = pnorm(qnorm(parms)-1.96*se_orig)
upper = pnorm(qnorm(parms)+1.96*se_orig)

low_ind = true>=lower;mean(1-low_ind)
up_ind = true<=upper;mean(1-up_ind)
cov_orig = mean(low_ind*up_ind)
cov_orig

mean((0.5<lower)|(0.5>upper))

max(se_adj-se_orig)



sqrt(var(parms))
mean(se_adj)
mean(se_orig)

sqrt(var(parms_gee))
mean(se_gee)




lower = plogis(-(parms_gee+1.96*se_gee))
upper = plogis(-(parms_gee-1.96*se_gee))

low_ind = true>=lower;mean(1-low_ind)
up_ind = true<=upper;mean(1-up_ind)
cov_gee = mean(low_ind*up_ind)
cov_gee


