## Data are simulated according to the following scenario:
####  - nr_of subjects in total, equally divided over group0 and group1
####  - each followed for nr_of_times timepoints
####  - model: Y = beta0+beta1*time+beta2*group+beta3*time*group+epsilon_ij+u_i
############## u_i is random intercept per subject
library(pim)
library(tidyverse)
nr_of_subjects=20
use_times = nr_of_times = 3
id = rep(1:nr_of_subjects,each=nr_of_times)
time = rep(1:nr_of_times,nr_of_subjects)
group = rep(rep(c(0,1),each=nr_of_subjects/2),each=nr_of_times)

beta0=0
beta1=sqrt(2)*0.5
beta2=sqrt(2)*0.5
beta3 = 2

design  = data.frame(id,time,group)

#set.seed(123)
design$Y = beta0+beta1*design$time+beta2*(design$group==1)+beta3*design$time*(design$group==1)+rnorm(nrow(design),0,1)+rep(rnorm(nr_of_subjects,0,1.5),each=nr_of_times)

## Below we fit several PIMs to model the effect over time (see overleaf for definition)


## Old sandwich estimator (only sums U over shared indices)
sandwich.estimator<-function(U, U.diff, 
                             g1, g2, 
                             shared.factor=1, switched.factor=1,
                             self.factor=1)
{

  #note: this estimator is greatly optimized, based on T Lumley's code!
  #print(U.diff)
  if(any(fullyequal <- g1==g2 ))
  {
    #the contribution of these terms to variance and covariance is zero anyway, because
    #I_{ii} = 0.5 always (ie: constant), so we drop them out!
    U<-U[!fullyequal,,drop=FALSE]
    g1<-g1[!fullyequal]
    g2<-g2[!fullyequal]
  }
  
  usum1.tmp <- rowsum(U,g1,reorder=FALSE)
  usum2.tmp <- rowsum(U,g2,reorder=FALSE)
  
  ng <- length(union(g1,g2))
  
  usum1  <- matrix(nrow = ng, ncol = ncol(usum1.tmp),0)
  usum2  <- matrix(nrow = ng, ncol = ncol(usum2.tmp),0)
  usum1[unique(g1),] <- usum1.tmp
  usum2[unique(g2),] <- usum2.tmp
  
  utushared <- crossprod(usum1) + crossprod(usum2) 
  utuswitched <- crossprod(usum1,usum2) + crossprod(usum2,usum1)
  uDiag<-crossprod(U) #Is counted twice as shared.factor, but needs to be counted as self.factor
  

  utu<-shared.factor*utushared  + 
    (switched.factor)*utuswitched +
    (self.factor-2*shared.factor)*(uDiag)
  
  #if the inverses occur (ij, ji), they are counted doubly as switched!!
  #However: they should be counted as - self
  
  mx<-ng+1
  
  uids<-g1*mx+g2
  invuids<-g2*mx+g1
  invrowperrow<-match(uids, invuids)
  rowsWithInv<-(!is.na(invrowperrow)) 
  if(any(rowsWithInv))
  {
    rowsWithInv<-which(rowsWithInv)
    tmp<-sapply(rowsWithInv, function(ri){
      res<-U[ri,] %*% t(U[invrowperrow[ri],]) 
      return(res)
    })
    if(is.null(dim(tmp)))
    {
      if(length(rowsWithInv)==1)
      {
        tmp<-matrix(tmp, ncol=1)
      }
      else
      {
        tmp<-matrix(tmp, nrow=1)
      }
    }
    
    uijji<-matrix(rowSums(tmp), ncol=ncol(usum1.tmp))
    
    utu<-utu-(2 *  switched.factor + self.factor)*uijji
  }
  
  
  
   individuals_var2 = g2
   individuals_var1 = g1
  #
  utu_test = matrix(0,nrow=ncol(U),ncol=ncol(U))
  phi = ((outer(individuals_var1, individuals_var1, FUN = "==")+outer(individuals_var2, individuals_var2, FUN = "==")+outer(individuals_var1, individuals_var2, FUN = "==")+outer(individuals_var2, individuals_var1, FUN = "=="))>0)*1
  for(j in 1:nrow(U)){
    for(k in 1:nrow(U)){
      utu_test = utu_test+ (phi[j,k]*U[j,]%*%t(U[k,]))
    }
   }
  rownames(utu)=rownames(utu_test)=colnames(utu_test) = colnames(utu) = c()
  print(utu)
  print(utu_test)
  
  #skip this last part if U.diff is identity... (note: NULL here will mean identity)
  if(is.null(U.diff)) return (utu)
  
  U.diff.inv<-solve(U.diff)
  
  return(U.diff.inv%*%utu%*%U.diff.inv)
}

# We will also adjust the sandwich estimate
sandwich.estimator_clustered_mod1_mod2<- function(U, U.diff, 
                                             g1, g2, 
                                             shared.factor = 1, switched.factor = 1,
                                             self.factor = 1) {
  # Note: this estimator is greatly optimized, based on T Lumley's code!

  # We transform the original g1 and g2 from the comparison of two groups to indicators corresponding to
  # the overarching subject ids. This leads to additional terms that are counted twice in utu (as compared to
  # only the diagonal elements in the original code). This is accounted for using the subtract_df below. All shared/switched/selffactors are fixed to 1
  
  g1 = individuals_var1
  g2 = individuals_var2
  
  use_n1 = length(unique(g1))
  tim = proc.time()[3]
  
  # if (any(fullyequal <- g1 == g2)) {
  #   U <- U[!fullyequal, , drop = FALSE]
  #   g1 <- g1[!fullyequal]
  #   g2 <- g2[!fullyequal]
  # }
  
  usum1.tmp <- rowsum(U, g1, reorder = FALSE)
  usum2.tmp <- rowsum(U, g2, reorder = FALSE)
  
  ng <- length(union(g1, g2))
  
  usum1 <- matrix(nrow = ng, ncol = ncol(usum1.tmp), 0)
  usum2 <- matrix(nrow = ng, ncol = ncol(usum2.tmp), 0)
  usum1[unique(g1), ] <- usum1.tmp
  usum2[unique(g2), ] <- usum2.tmp
  
  
  utushared <- crossprod(usum1) + crossprod(usum2) 
  utuswitched <- crossprod(usum1, usum2) + crossprod(usum2, usum1)
  
  
  # if (any(fullyequal <- g1 == g2)) {
  #   U <- U[!fullyequal, , drop = FALSE]
  #   g1 <- g1[!fullyequal]
  #   g2 <- g2[!fullyequal]
  # }
  
  subtract_df = data.frame(g1,g2)
  subtract_df = subtract_df %>%group_by(g1, g2)%>% mutate(ID = cur_group_id()) 
  if(length(unique(subtract_df$ID)) == length(g1))
  {uDiag <- crossprod(U)  
  utu <- shared.factor * utushared + 
    (switched.factor) * utuswitched +
    (self.factor - 2 * shared.factor) * (uDiag)
  }
  
  
  if(length(unique(subtract_df$ID)) < length(g1))
  { 
    usum3.tmp <- rowsum(U, subtract_df$ID, reorder = FALSE)
    subtract = crossprod(usum3.tmp[1:(use_n1/2*use_n1/2),]) 
    utu <- utushared - subtract 
    utu[ncol(U),ncol(U)] = utu[ncol(U),ncol(U)]+utuswitched[ncol(U),ncol(U)]
    #usum3.tmp <- rowsum(U[76:105,], subtract_df$ID[76:105], reorder = FALSE)
    subtract = crossprod(usum3.tmp[(use_n1/2*use_n1/2+1):(use_n1+use_n1/2*use_n1/2),])  
    utu <- utu - subtract
    
  }
  
  
  
  mx <- ng + 1
  
  # uids <- g1 * mx + g2
  # invuids <- g2 * mx + g1
  # invrowperrow <- match(uids, invuids)
  # rowsWithInv <- (!is.na(invrowperrow)) 
  # 
  # if (any(rowsWithInv)) {
  #   rowsWithInv <- which(rowsWithInv)
  #   tmp <- sapply(rowsWithInv, function(ri) {
  #     res <- U[ri, ] %*% t(U[invrowperrow[ri], ]) 
  #     return(res)
  #   })
  #   
  #   if (is.null(dim(tmp))) {
  #     if (length(rowsWithInv) == 1) {
  #       tmp <- matrix(tmp, ncol = 1)
  #     } else {
  #       tmp <- matrix(tmp, nrow = 1)
  #     }
  #   }
  # 
  #   uijji <- matrix(rowSums(tmp), ncol = ncol(usum1.tmp))
  #   print(uijji)
  #   utu <- utu - (2 * switched.factor + self.factor) * uijji
  # }
  
  # Skip this last part if U.diff is identity... (note: NULL here will mean identity)
  if (is.null(U.diff)) return(utu)
  
  U.diff.inv<-solve(U.diff)
  
  utu_test = matrix(0,nrow=ncol(U),ncol=ncol(U))

  phi = (((outer(individuals_var1, individuals_var1, FUN = "==")+outer(individuals_var2, individuals_var2, FUN = "==")+outer(individuals_var1, individuals_var2, FUN = "==")+outer(individuals_var2, individuals_var1, FUN = "==")))>0)*1

  for(j in 1:nrow(U)){
    for(k in 1:nrow(U)){
      utu_test = utu_test+ (phi[j,k]*(U[j,]%*%t(U[k,])))
    }
  }
  rownames(utu)=rownames(utu_test)=colnames(utu_test) = colnames(utu) = c()
  

  return(U.diff.inv%*%utu_test%*%U.diff.inv)
}

sandwich.estimator_clustered_mod3 <- function(U, U.diff, 
                                              g1, g2, 
                                              shared.factor = 1, switched.factor = 1,
                                              self.factor = 1) {
  # Note: this estimator is greatly optimized, based on T Lumley's code!
  #print(U.diff)
  # We transform the original g1 and g2 from the comparison of two groups to indicators corresponding to
  # the overarching subject ids. This leads to additional terms that are counted twice in utu (as compared to
  # only the diagonal elements in the original code). This is accounted for using the subtract_df below. All shared/switched/selffactors are fixed to 1
  
  g1 = individuals_var1
  g2 = individuals_var2
  
  use_n1 = length(unique(g1))
  
  tim = proc.time()[3]
  
  # if (any(fullyequal <- g1 == g2)) {
  #   U <- U[!fullyequal, , drop = FALSE]
  #   g1 <- g1[!fullyequal]
  #   g2 <- g2[!fullyequal]
  # }
  
  usum1.tmp <- rowsum(U, g1, reorder = FALSE)
  usum2.tmp <- rowsum(U, g2, reorder = FALSE)
  
  ng <- length(union(g1, g2))
  
  usum1 <- matrix(nrow = ng, ncol = ncol(usum1.tmp), 0)
  usum2 <- matrix(nrow = ng, ncol = ncol(usum2.tmp), 0)
  usum1[unique(g1), ] <- usum1.tmp
  usum2[unique(g2), ] <- usum2.tmp
  
  
  utushared <- crossprod(usum1) + crossprod(usum2) 
  utuswitched <- crossprod(usum1, usum2) + crossprod(usum2, usum1)
  
  
  # if (any(fullyequal <- g1 == g2)) {
  #   U <- U[!fullyequal, , drop = FALSE]
  #   g1 <- g1[!fullyequal]
  #   g2 <- g2[!fullyequal]
  # }
  
  subtract_df = data.frame(g1,g2)
  subtract_df = subtract_df %>%group_by(g1, g2)%>% mutate(ID = cur_group_id()) 
  if(length(unique(subtract_df$ID)) == length(g1))
  {uDiag <- crossprod(U)  
  utu <- shared.factor * utushared + 
    (switched.factor) * utuswitched +
    (self.factor - 2 * shared.factor) * (uDiag)
  }
  
  
  if(length(unique(subtract_df$ID)) < length(g1))
  { 
    
    usum3.tmp <- rowsum(U, subtract_df$ID, reorder = FALSE)
    subtract = crossprod(usum3.tmp[1:(use_n1/2*use_n1/2),]) 
    utu <- utushared - subtract 
    utu[ncol(U),ncol(U)] = utu[ncol(U),ncol(U)]+utuswitched[ncol(U),ncol(U)]
    #usum3.tmp <- rowsum(U[76:105,], subtract_df$ID[76:105], reorder = FALSE)
    subtract = crossprod(usum3.tmp[(use_n1/2*use_n1/2+1):(use_n1+use_n1/2*use_n1/2),])  
    utu <- utu - subtract
    
  }
  
  
  
  mx <- ng + 1
  
  # uids <- g1 * mx + g2
  # invuids <- g2 * mx + g1
  # invrowperrow <- match(uids, invuids)
  # rowsWithInv <- (!is.na(invrowperrow)) 
  # 
  # if (any(rowsWithInv)) {
  #   rowsWithInv <- which(rowsWithInv)
  #   tmp <- sapply(rowsWithInv, function(ri) {
  #     res <- U[ri, ] %*% t(U[invrowperrow[ri], ]) 
  #     return(res)
  #   })
  #   
  #   if (is.null(dim(tmp))) {
  #     if (length(rowsWithInv) == 1) {
  #       tmp <- matrix(tmp, ncol = 1)
  #     } else {
  #       tmp <- matrix(tmp, nrow = 1)
  #     }
  #   }
  # 
  #   uijji <- matrix(rowSums(tmp), ncol = ncol(usum1.tmp))
  #   print(uijji)
  #   utu <- utu - (2 * switched.factor + self.factor) * uijji
  # }
  
  # Skip this last part if U.diff is identity... (note: NULL here will mean identity)
  if (is.null(U.diff)) return(utu)
  
  U.diff.inv<-solve(U.diff)
  
  # utu_test = matrix(0,nrow=ncol(U),ncol=ncol(U))
  # 
  # phi = (((outer(individuals_var1, individuals_var1, FUN = "==")+outer(individuals_var2, individuals_var2, FUN = "==")+outer(individuals_var1, individuals_var2, FUN = "==")+outer(individuals_var2, individuals_var1, FUN = "==")))>0)*1
  # 
  # for(j in 1:nrow(U)){
  #   for(k in 1:nrow(U)){
  #     utu_test = utu_test+ (phi[j,k]*(U[j,]%*%t(U[k,])))
  #   }
  # }
  # colnames(utu_test) = colnames(utu) = c()
  # print(utu)
  # print(utu_test)
  
  return(U.diff.inv%*%utu%*%U.diff.inv)
}

## Model 1 and 2 have the same poset, defined via compare. This also translates in
# common patient indicators individuals_var1 and individuals_var2 (used to calculate new sandwich estimate)
subj_g = nr_of_subjects/2
individuals_var2 = c(rep((subj_g+1):nr_of_subjects ,each=subj_g*nr_of_times),rep(1:subj_g,each=nr_of_times-1),rep((subj_g+1):nr_of_subjects,each=nr_of_times-1),1:nr_of_subjects)
individuals_var1 = c(rep(1:subj_g,subj_g*nr_of_times),rep(1:subj_g,each=nr_of_times-1),rep((subj_g+1):nr_of_subjects,each=nr_of_times-1),1:nr_of_subjects)
# Poset definition

id.fac1 <- which((design[,"group"] == 0)&(design[,"time"] == 1))
id.nonfac1 <- which((design[,"group"] == 1)&(design[,"time"] == 1))
compare1 <- expand.grid(id.fac1,id.nonfac1)

id.fac2 <- which((design[,"group"] == 0)&(design[,"time"] == 2))
id.nonfac2 <- which((design[,"group"] == 1)&(design[,"time"] == 2))
compare2 <- expand.grid(id.fac2,id.nonfac2)

id.fac3 <- which((design[,"group"] == 0)&(design[,"time"] == 3))
id.nonfac3 <- which((design[,"group"] == 1)&(design[,"time"] == 3))
compare3 <- expand.grid(id.fac3,id.nonfac3)

compare_between = rbind(compare1,compare2,compare3)
compare_between = compare_between[order(compare_between[,"Var2"]),]

within1 = c((1:(nr_of_subjects*nr_of_times))[-seq(nr_of_times,nr_of_subjects*nr_of_times,nr_of_times)],seq(1,nr_of_subjects*nr_of_times-(nr_of_times-1),nr_of_times))
within2 = c((1:(nr_of_subjects*nr_of_times))[-seq(1,nr_of_subjects*nr_of_times-(nr_of_times-1),nr_of_times)],seq(nr_of_times,nr_of_subjects*nr_of_times,nr_of_times))

compare_within = cbind(within1,within2)
colnames(compare_within) = c("Var1","Var2")
compare = rbind(compare_between,compare_within)

assignInNamespace("sandwich.estimator", sandwich.estimator, ns = "pim")

mod1_orig = pim(Y~1+time+group,data=design,compare=compare,link="probit")
summary(mod1_orig)

assignInNamespace("sandwich.estimator", sandwich.estimator_clustered_mod1_mod2, ns = "pim")

mod1_adj = pim(Y~1+time+group,data=design,compare=compare,link="probit")
summary(mod1_adj)

assignInNamespace("sandwich.estimator", sandwich.estimator, ns = "pim")

mod2_orig = pim(Y~1+I((R(time) - L(time))*R(group)*L(group))+I((R(time) - L(time))*(1-R(group))*(1-L(group)))+I((R(group)-L(group))*(R(time)==1)*(L(time)==1))+I((R(group)-L(group))*(R(time)==2)*(L(time)==2))+I((R(group)-L(group))*(R(time)==3)*(L(time)==3)) ,data=design,compare=compare,link="probit")
summary(mod2_orig)

assignInNamespace("sandwich.estimator", sandwich.estimator_clustered_mod1_mod2, ns = "pim")

mod2_adj = pim(Y~1+I((R(time) - L(time))*R(group)*L(group))+I((R(time) - L(time))*(1-R(group))*(1-L(group)))+I((R(group)-L(group))*(R(time)==1)*(L(time)==1))+I((R(group)-L(group))*(R(time)==2)*(L(time)==2))+I((R(group)-L(group))*(R(time)==3)*(L(time)==3)) ,data=design,compare=compare,link="probit")
summary(mod2_adj)

## Model 3 has a different poset, but allows for an interaction

individuals_var2 = c(rep((subj_g+1):nr_of_subjects ,each=subj_g*nr_of_times))
individuals_var1 = c(rep(1:subj_g,subj_g*nr_of_times))
head(compare)
#poset is compare_between defined above

assignInNamespace("sandwich.estimator", sandwich.estimator, ns = "pim")

mod3_orig = pim(Y~I((R(group)-L(group))*((R(time)==1) & (L(time)==1)))+I((R(group)-L(group))*((R(time)==2) & (L(time)==2)))+I((R(group)-L(group))*((R(time)==3) & (L(time)==3))),data=design,compare=compare_between,link="probit")
summary(mod3_orig)

assignInNamespace("sandwich.estimator", sandwich.estimator_clustered_mod3, ns = "pim")
# mod3_OT = pim(Y~1+I((R(group)-L(group))*((R(time)==1) & (L(time)==1)))+I((R(group)-L(group))*((R(time)==2) & (L(time)==2)))+I((R(group)-L(group))*((R(time)==3) & (L(time)==3))),data=design,compare=compare_between,link="probit")
# summary(mod3_OT) # not able to fit

mod3_adj = pim(Y~I((R(group)-L(group))*((R(time)==1) & (L(time)==1)))+I((R(group)-L(group))*((R(time)==2) & (L(time)==2)))+I((R(group)-L(group))*((R(time)==3) & (L(time)==3))),data=design,compare=compare_between,link="probit")
summary(mod3_adj) # removed intercept



