path<-"C:/Users/DellPC/Dropbox/Abdoul/These/Article1/Code/MDSV/Code_These/Ess"
path<-"/home/maoudek/Rsim/Article1/MDSV/Returns"
path<-"C:/Users/DellPC/Dropbox/Abdoul/These/Article1/Code/MDSV/Code_These/MDSV_package/MDSV/test/These/Chapitre3"
path<-"C:/Users/DellPC/Dropbox/Abdoul/These/Article1/Code/MDSV/Code_These/Ess/forcastArticle"
setwd(path)
library(MDSV)
if(!require(Rsolnp)){install.packages("Rsolnp")}; library(Rsolnp)
if(!require(actuar)){install.packages("actuar")}; library(actuar)
if(!require(KScorrect)){install.packages("KScorrect")}; library(KScorrect)
if(!require(mhsmm)){install.packages("mhsmm")}; library(mhsmm)
if(!require(Rcpp)){install.packages("Rcpp")}; library(Rcpp)

index_set<-c("aex", "aord", "bell", "bsesn", "bvsp", "cac40", "dax", "dji", "ftmib", "ftse",
             "hsi", "ibex35", "kospi", "kse", "mxx", "n225", "nasdaq", "nifty50", "omxc20", "omxhpi",
             "omxspi", "oseax", "psi", "rut", "smsi", "sp500", "ssec", "ssmi", "sti", "stoxx50", "tsx")

index_set<-c("sp500", "nasdaq", "ftse", "n225")

start.date <- as.Date("2000-01-01") 
end.date   <- as.Date("2019-12-31")

ctrl <- list(TOL=1e-15, trace=0)

#-------------------------------------------------------------------------------
# Real-EGARCH : ESTIMATION
#-------------------------------------------------------------------------------

#### Real-EGARCH

para_names<-function(LEVIER){
  vars.names<-c("w","beta","tau1","tau2","gamma","xi","phi","sigma")
  if(LEVIER) vars.names<-c(vars.names,"tau1l","tau2l")
  return(vars.names)
}

# Quelques fonctions
workNat<-function(para_tilde,LEVIER){
  para<-para_tilde #c("w","beta","tau_1","tau_2","gamma","xi","phi",...)
  para[8]<-exp(para_tilde[8])    #sigma_{epsilon}
  return(para)
}

natWork<-function(para,LEVIER){
  para_tilde<-para #c("w","beta","tau_1","tau_2","gamma","xi","phi",...)
  para_tilde[8]<-log(para[8])    #sigma_{epsilon}
  return(para_tilde)
}

logLik<-function(ech,para_tilde,LEVIER,r=0){
  para<-workNat(para_tilde,LEVIER)
  n<-length(ech[,"r"])
  sigma<-((sqrt(1/(1-2*para[4])))*exp(para[1]+(para[5]^2)*para[8]/2-para[4]+(para[3]^2)/(2*(1-para[4]))))^(1/(1-para[2])) # variance inconditionnel
  mu_rv<-para[6] + para[7]*log(sigma)+para[3]*(ech[1,"r"]/sqrt(sigma))+para[4]*((ech[1,"r"]/sqrt(sigma))^2-1)
  aj<-dnorm(ech[1,"r"]/sqrt(sigma))*dlnorm(ech[1,"rv"],mu_rv,sqrt(para[8]))/sqrt(sigma)
  ajm<-dnorm(ech[1,"r"]/sqrt(sigma))/sqrt(sigma)
  lik<-log(aj)
  likm<-log(ajm)
  
  if(!LEVIER){
    para<-c(para,0,0)
  }
  
  for(i in 2:n){
    sigma<-exp(para[1]+para[2]*log(sigma)+para[5]*(log(ech[i-1,"rv"])-mu_rv)+para[9]*(ech[i-1,"r"]/sqrt(sigma))+para[10]*((ech[i-1,"r"]/sqrt(sigma))^2-1))
    mu_rv<-para[6] + para[7]*log(sigma)+para[3]*(ech[i,"r"]/sqrt(sigma))+para[4]*((ech[i,"r"]/sqrt(sigma))^2-1)
    aj<-dnorm(ech[i,"r"]/sqrt(sigma))*dlnorm(ech[i,"rv"],mu_rv,sqrt(para[8]))/sqrt(sigma)
    ajm<-dnorm(ech[i,"r"]/sqrt(sigma))/sqrt(sigma)
    lik<-lik+log(aj)
    likm <-likm+log(ajm)
  }
  sigma<-exp(para[1]+para[2]*log(sigma)+para[5]*(log(ech[n,"rv"])-mu_rv)+para[9]*(ech[n,"r"]/sqrt(sigma))+para[10]*((ech[n,"r"]/sqrt(sigma))^2-1))
  aj<-dnorm(r/sqrt(sigma))/sqrt(sigma)
  likp<-log(aj)
  
  attr(lik,"Pred_lik") <- likp
  attr(lik,"Marg_lik") <- likm
  return(-lik);
}

model<-expand.grid(index=index_set, start.date=as.Date("2000-01-01"), end.date=as.Date("2019-12-31"), length=0,
                   LEVIER=c(FALSE,TRUE))

vars<-c('loglikm','loglik', 'AIC', 'BIC', c("w","beta","tau1","tau2","gamma","xi","phi","sigma","tau1l","tau2l"),"conv","time")

model_add <- matrix(0, nrow=nrow(model), ncol=length(vars))
colnames(model_add) <- vars
model <- cbind(model, model_add)

filename <- paste("RealEGARCH","_all_index_",start.date,"_",end.date, sep="")

for(i in 1:nrow(model)){
  start_time <- Sys.time()
  index                 <- as.character(model[i,"index"])
  donne                 <- get(index)[rownames(get(index))>=start.date & rownames(get(index))<=end.date,]
  nam_                  <- rownames(donne)
  donne[,"r"]           <- donne[,"r"] - mean(donne[,"r"])
  donne                 <- cbind("rv"=donne[,"rv"], "r"=donne[,"r"])
  model[i,"start.date"] <- as.Date(nam_[1])
  model[i,"end.date"]   <- as.Date(nam_[length(nam_)])
  model[i,"length"]     <- length(nam_)
  dist                  <- as.character(model[i,"dist"])
  LEVIER                <- model[i,"LEVIER"]
  
  para<-c(-0.03,0.7,0.01,0.04,0.18,0.12,1.12,5)
  vars<-c("w","beta","tau1","tau2","gamma","xi","phi","sigma")
  if(LEVIER) {
    para<-c(para,-0.005,0.043)
    vars<-c(vars,"tau1l","tau2l")
  }
  para_tilde<-natWork(para,LEVIER)
  oldw <- getOption("warn")
  options(warn = -1)
  opt<-try(solnp(pars=para_tilde,fun=logLik,ech=donne,LEVIER=LEVIER,control=ctrl),silent=T)
  options(warn = oldw)
  
  if(is(opt,"try-error")) {
    model[i,'conv']  <- "No Convergence."
    print(paste("===",round(100*i/nrow(model),2) , "% : ERROR" ,"===="))
    next
  }
  
  params        <- workNat(opt$pars,LEVIER=LEVIER)
  names(params) <- para_names(LEVIER)
  model[i,colnames(model) %in% names(params)] <- round(params[vars],5)

  model[i,'conv']   <- "Convergence."
  model[i,'loglik'] <- -as.numeric(opt$values[length(opt$values)])
  model[i,"loglikm"]<- attr(logLik(donne,opt$pars,LEVIER),"Marg_lik")
  model[i,'AIC']    <- model[i,"loglik"] - length(params)
  model[i,'BIC']    <- model[i,"loglik"] - length(params)*log(length(nam_))/2
  model[i,"time"]   <- difftime(Sys.time(), start_time,units = "secs")[[1]]
  
  #write.csv(model, paste(filename,"csv",sep="."), row.names=FALSE)
  print(paste("===",round(100*i/nrow(model),2) , "%" ,"==== LOGLIK = ", model[i,"loglik"], "==="))
}

write.csv(model, "Estim_ALLRealEGARCH1.csv", row.names=FALSE)

#-------------------------------------------------------------------------------
# Real-EGARCH : PREVISION
#-------------------------------------------------------------------------------

g<-function(vector){
  Sortie<-matrix(NA,2,2)
  for(i in c(0,1)) for(j in c(0,1)) Sortie[i+1,j+1]<-sum((vector[-1]==j)*(vector[-length(vector)]==i))
  if(vector[1]==0) Sortie[1,1]<-Sortie[1,1]+1
  if(vector[1]==1) Sortie[1,2]<-Sortie[1,2]+1
  colnames(Sortie)<-c(0,1)
  rownames(Sortie)<-c(0,1)
  return(Sortie)
}

model_extern<-expand.grid(window = c(252, 504, 756, 1008, 1260, 1512), index=index_set, LEVIER=c(TRUE), n.ahead=c(100),
                          forecast.length=c(1512), refit.every=c(21), refit.window=c("recursive"),
                          calculate.VaR=c(TRUE), rseed=1050)

Loss.horizon    <- c(1,5,10,25,50,75,100)
n.ahead         <- 100
forecast.length <- 1512
refit.every     <- 21
refit.window    <- "recursive"
VaR.alpha       <- c(0.01,0.05,0.10)
rseed           <- 1050

vars<-c('start.date','end.date','pred_dens',paste('QLIKc_R',Loss.horizon), paste('RMSEc_R',Loss.horizon), 
        paste('MAEc_R',Loss.horizon), paste('QLIKm_R',Loss.horizon), paste('RMSEm_R',Loss.horizon),
        paste('MAEm_R',Loss.horizon), paste('QLIKc_RV',Loss.horizon), paste('RMSEc_RV',Loss.horizon), 
        paste('MAEc_RV',Loss.horizon), paste('QLIKm_RV',Loss.horizon), paste('RMSEm_RV',Loss.horizon),
        paste('MAEm_RV',Loss.horizon), paste('LR.uc_stat',100*(1-c(0.01,0.05,0.10))), 
        paste('LR.uc_pvalue',100*(1-c(0.01,0.05,0.10))), paste('LR.ind_stat',100*(1-c(0.01,0.05,0.10))), 
        paste('LR.ind_pvalue',100*(1-c(0.01,0.05,0.10))), paste('LR.cc_stat',100*(1-c(0.01,0.05,0.10))), 
        paste('LR.cc_pvalue',100*(1-c(0.01,0.05,0.10))),"time")

model_add <- matrix(0, nrow=nrow(model_extern), ncol=length(vars))
colnames(model_add) <- vars
model_extern <- cbind(model_extern, model_add)
R_var <- rep(0,length(Loss.horizon))

for(i in c(1,7,13,19)){
  start_time <- Sys.time()
  index                        <- as.character(model_extern[i,"index"])
  donne                        <- get(index)[rownames(get(index))>=start.date & rownames(get(index))<=end.date,]
  nam_                         <- rownames(donne)
  donne[,"r"]                  <- donne[,"r"] - mean(donne[,"r"])
  donne                        <- cbind("rv"=donne[,"rv"], "r"=donne[,"r"])
  model_extern[i:(i+5),"start.date"] <- as.character(as.Date(nam_[length(nam_)])-forecast.length)
  model_extern[i:(i+5),"end.date"]   <- as.character(nam_[length(nam_)])
  dist                         <- as.character(model_extern[i,"dist"])
  LEVIER                       <- model_extern[i,"LEVIER"]
  calculate.VaR                <- as.logical(model_extern[i,"calculate.VaR"])
  
  set.seed(rseed)
  
  filename <- paste("Forecast_RealEGARCH_LEVIER_",LEVIER,"_", index, "_length_",forecast.length,
                    "_", model_extern[i,"start.date"], "_", model_extern[i,"end.date"], sep="")
  
  model<-expand.grid(date = nam_[((length(nam_)-forecast.length+1):length(nam_))],rt=0,rvt=0,Levier=LEVIER)
  
  model$rt  <- donne[((length(nam_)-forecast.length+1):length(nam_)),"r"]
  model$rvt  <- donne[((length(nam_)-forecast.length+1):length(nam_)),"rv"]
  
  vars<-c(paste0("R_for",Loss.horizon),paste0("R_tru",Loss.horizon),
          paste0("R_for_m",Loss.horizon),paste0("R_tru_m",Loss.horizon),
          paste0("RV_for",Loss.horizon),paste0("RV_tru",Loss.horizon),
          paste0("RV_for_m",Loss.horizon),paste0("RV_tru_m",Loss.horizon),'pred_lik','loglik', 'AIC', 'BIC',
          c("w","beta","tau1","tau2","gamma","xi","phi","sigma","tau1l","tau2l"), paste0("VaR",100*(1-VaR.alpha)), 
          paste0("I",100*(1-VaR.alpha)))
  
  model_add <- matrix(0, nrow=nrow(model), ncol=length(vars))
  colnames(model_add) <- vars
  model <- cbind(model, model_add)
  
  para            <- c(-0.03,0.7,0.01,0.04,0.18,0.12,1.12,5)
  if(LEVIER) para <- c(para,-0.005,0.043)
  para_tilde      <- natWork(para,LEVIER)
  opt             <- NULL
  update_date     <- seq(0,forecast.length,by=refit.every)
  strt            <- 1
  
  for(t in 0:(forecast.length-1)){
    ech    <- donne[strt:(length(nam_)-forecast.length+t),]
    
    if(t %in% update_date){
      if((!is.null(opt)) & !(is(opt,"try-error"))) para_tilde<-opt$pars
      oldw <- getOption("warn")
      options(warn = -1)
      opt<-try(solnp(pars=para_tilde,fun=logLik,ech=ech,LEVIER=LEVIER,control=ctrl),silent=T)
      options(warn = oldw)
      
      if(refit.window == "moving") strt <- strt + refit.every
      
      para<-workNat(opt$pars,LEVIER=LEVIER)
      names(para)<-para_names(LEVIER)
    }

    model[t+1,colnames(model) %in% names(para)] <- round(para,5)
    l<-logLik(ech,opt$pars,LEVIER,r=donne[(length(nam_)-forecast.length+t+1),"r"])
    
    model[t+1,'loglik']   <- -as.numeric(opt$values[length(opt$values)])
    model[t+1,"loglikm"]  <- attr(l,"Marg_lik")
    model[t+1,"pred_lik"] <- attr(l,"Pred_lik")
    model[t+1,'AIC']      <- model[t+1,"loglik"] - length(para)
    model[t+1,'BIC']      <- model[t+1,"loglik"] - length(para)*log(nrow(ech))/2
    
    n <- 10000 #number of simulations
    H <- max(Loss.horizon)    #nb of years simulated
    
    w       <- para[1]
    beta    <- para[2]
    tau_1   <- para[3] 
    tau_2   <- para[4]
    gamma   <- para[5]
    xi      <- para[6]  
    phi     <- para[7]   
    sigma   <- para[8]
    if(LEVIER){
      tau1l <- para[9]    #tau1l
      tau2l <- para[10]    #tau2l  
    }else{
      tau1l <- 0    #tau1l
      tau2l <- 0    #tau2l  
    }
    
    z_t<-matrix(rnorm(n*H),nrow=n,ncol=H)
    
    e_t<-matrix(rnorm(n*H,0,sd=sqrt(sigma)),nrow=n,ncol=H)
    rv_t<-matrix(NA,n,H)
    r_t<-matrix(NA,n,H)
    
    h_t<-var(ech[,"r"])
    for(b in 1:nrow(ech)) h_t<-exp(w+beta*log(h_t)+tau1l*(ech[b,"r"]/sqrt(h_t))+tau2l*((ech[b,"r"]/sqrt(h_t))^2-1)+gamma*(log(ech[b,"rv"])-xi-phi*log(h_t)-tau_1*(ech[b,"r"]/sqrt(h_t))-tau_2*((ech[b,"r"]/sqrt(h_t))^2-1)))
    
    if(calculate.VaR) for(iter in 1:length(VaR.alpha)){
      model[t+1,paste0('VaR',100*(1-VaR.alpha[iter]))] <- qnorm(VaR.alpha[iter],0,sqrt(h_t))
    }
    
    h_t<-h_t*matrix(rep(1,n),n,1)
    if(!LEVIER){
      for(h in 1:(H)){
        r_t[,h]<-sqrt(h_t)*z_t[,h]
        rv_t[,h]<-exp(xi+phi*log(h_t)+tau_1*(z_t[,h])+tau_2*((z_t[,h])^2-1)+(e_t[,h]))
        h_t<-exp(w+beta*log(h_t)+gamma*(e_t[,h]))
      }
    }else{
      for(h in 1:(H)){
        rv_t[,h]<-exp(xi+phi*log(h_t)+tau_1*(z_t[,h])+tau_2*((z_t[,h])^2-1)+(e_t[,h]))
        r_t[,h]<-sqrt(h_t)*z_t[,h]
        h_t<-exp(w+beta*log(h_t)+gamma*(e_t[,h])+tau1l*(z_t[,h])+tau2l*((z_t[,h])^2-1))
      }
    }
    
    rt_sim<-colMeans(r_t^2)
    
    for(k in 1:length(Loss.horizon)) R_var[k] <- sum(rt_sim[1:Loss.horizon[k]])
    names(R_var) <- paste0("R_for",Loss.horizon)
    model[(t+1),colnames(model) %in% names(R_var)]<-R_var
    
    for(k in 1:length(Loss.horizon)) R_var[k] <- sum(rt_sim[Loss.horizon[k]])
    names(R_var) <- paste0("R_for_m",Loss.horizon)
    model[(t+1),colnames(model) %in% names(R_var)]<-R_var
    
    ech    <- donne[(length(nam_)-forecast.length+t+1):length(nam_),"r"]
    for(k in 1:length(Loss.horizon)) R_var[k] <- sum((ech[1:Loss.horizon[k]])^2)
    names(R_var) <- paste0("R_tru",Loss.horizon)
    model[(t+1),colnames(model) %in% names(R_var)]<-R_var
    
    for(k in 1:length(Loss.horizon)) R_var[k] <- sum((ech[Loss.horizon[k]])^2)
    names(R_var) <- paste0("R_tru_m",Loss.horizon)
    model[(t+1),colnames(model) %in% names(R_var)]<-R_var
    
    rvt_sim<-colMeans(rv_t)
    
    for(k in 1:length(Loss.horizon)) R_var[k] <- sum(rvt_sim[1:Loss.horizon[k]])
    names(R_var) <- paste0("RV_for",Loss.horizon)
    model[(t+1),colnames(model) %in% names(R_var)]<-R_var
    
    for(k in 1:length(Loss.horizon)) R_var[k] <- sum(rvt_sim[Loss.horizon[k]])
    names(R_var) <- paste0("RV_for_m",Loss.horizon)
    model[(t+1),colnames(model) %in% names(R_var)]<-R_var
    
    ech    <- donne[(length(nam_)-forecast.length+t+1):length(nam_),"rv"]
    for(k in 1:length(Loss.horizon)) R_var[k] <- sum(ech[1:Loss.horizon[k]])
    names(R_var) <- paste0("RV_tru",Loss.horizon)
    model[(t+1),colnames(model) %in% names(R_var)]<-R_var
    
    for(k in 1:length(Loss.horizon)) R_var[k] <- sum(ech[Loss.horizon[k]])
    names(R_var) <- paste0("RV_tru_m",Loss.horizon)
    model[(t+1),colnames(model) %in% names(R_var)]<-R_var
    
  }
  
  if(calculate.VaR)  for(iter in 1:length(VaR.alpha)){
      model[,paste0('I',100*(1-VaR.alpha[iter]))] <- (model[,'rt'] < model[,paste0('VaR',100*(1-VaR.alpha[iter]))])
    }
  
  for(j in 0:5){
    window = c(252, 504, 756, 1008, 1260, 1512)[j+1]
    
    mod <- model[(1512-window+1):1512,]
    model_extern[i+j,'pred_dens'] <- sum(mod[,"pred_lik"])
    
    for_R_var   <- mod[,grep('R_for', colnames(mod), fixed=TRUE)]
    for_R_var   <- for_R_var[,-grep('R_for_m', colnames(for_R_var), fixed=TRUE)]
    for_R_err   <- mod[,grep('R_tru', colnames(mod), fixed=TRUE)]
    for_R_err   <- for_R_err[,-grep('R_tru_m', colnames(for_R_err), fixed=TRUE)]
    model_extern[i+j,paste('QLIKc_R',Loss.horizon)] <- colMeans(log(for_R_var) + for_R_err/for_R_var, na.rm=TRUE)
    model_extern[i+j,paste('RMSEc_R',Loss.horizon)] <- sqrt(colMeans( (for_R_var - for_R_err)^2, na.rm=TRUE ))/Loss.horizon
    model_extern[i+j,paste('MAEc_R',Loss.horizon)] <- colMeans( abs(for_R_var - for_R_err), na.rm=TRUE )/Loss.horizon
    
    for_R_var   <- mod[,grep('R_for_m', colnames(mod), fixed=TRUE)]
    for_R_err   <- mod[,grep('R_tru_m', colnames(mod), fixed=TRUE)]
    model_extern[i+j,paste('QLIKm_R',Loss.horizon)] <- colMeans(log(for_R_var) + for_R_err/for_R_var, na.rm=TRUE)
    model_extern[i+j,paste('RMSEm_R',Loss.horizon)] <- sqrt(colMeans( (for_R_var - for_R_err)^2, na.rm=TRUE ))/Loss.horizon
    model_extern[i+j,paste('MAEm_R',Loss.horizon)]  <- colMeans( abs(for_R_var - for_R_err), na.rm=TRUE )/Loss.horizon
    
    for_R_var   <- mod[,grep('RV_for', colnames(mod), fixed=TRUE)]
    for_R_var   <- for_R_var[,-grep('RV_for_m', colnames(for_R_var), fixed=TRUE)]
    for_R_err   <- mod[,grep('RV_tru', colnames(mod), fixed=TRUE)]
    for_R_err   <- for_R_err[,-grep('RV_tru_m', colnames(for_R_err), fixed=TRUE)]
    model_extern[i+j,paste('QLIKc_RV',Loss.horizon)] <- colMeans(log(for_R_var) + for_R_err/for_R_var, na.rm=TRUE)
    model_extern[i+j,paste('RMSEc_RV',Loss.horizon)] <- sqrt(colMeans( (for_R_var - for_R_err)^2, na.rm=TRUE ))/Loss.horizon
    model_extern[i+j,paste('MAEc_RV',Loss.horizon)] <- colMeans( abs(for_R_var - for_R_err), na.rm=TRUE )/Loss.horizon
    
    for_R_var   <- mod[,grep('RV_for_m', colnames(mod), fixed=TRUE)]
    for_R_err   <- mod[,grep('RV_tru_m', colnames(mod), fixed=TRUE)]
    model_extern[i+j,paste('QLIKm_RV',Loss.horizon)] <- colMeans(log(for_R_var) + for_R_err/for_R_var, na.rm=TRUE)
    model_extern[i+j,paste('RMSEm_RV',Loss.horizon)] <- sqrt(colMeans( (for_R_var - for_R_err)^2, na.rm=TRUE ))/Loss.horizon
    model_extern[i+j,paste('MAEm_RV',Loss.horizon)]  <- colMeans( abs(for_R_var - for_R_err), na.rm=TRUE )/Loss.horizon
    
    LR.uc_vec <- LR.cc_vec <- LR.ind_vec <- NULL
    p.uc_vec  <- p.cc_vec  <- p.ind_vec  <- NULL
    
    for(iter in 1:length(VaR.alpha)){
      viol      <- sum(mod[,paste0('I',100*(1-VaR.alpha[iter]))])
      alpha_hat <- sum(mod[,paste0('I',100*(1-VaR.alpha[iter]))])/forecast.length
      LR.uc     <- 2*log(((alpha_hat^viol)*((1-alpha_hat)^(forecast.length-viol)))/((VaR.alpha[iter]^viol)*((1-VaR.alpha[iter])^(forecast.length-viol))))
      
      viol      <- g(mod[,paste0('I',100*(1-VaR.alpha[iter]))])
      pi        <- (viol[1,2]+viol[2,2])/sum(viol)
      pi0       <- (viol[1,2])/(viol[1,1]+viol[1,2])
      pi1       <- (viol[2,2])/(viol[2,2]+viol[2,1])
      LR.ind    <- - 2*log((((1-pi)^(viol[1,1]+viol[2,1]))*(pi^(viol[1,2]+viol[2,2])))/
                             ((((1-pi0)^viol[1,1])*(pi0^viol[1,2]))*(((1-pi1)^viol[2,1])*(pi1^viol[2,2]))))
      
      LR.uc_vec <- c(LR.uc_vec,LR.uc)
      LR.cc_vec <- c(LR.cc_vec,LR.uc+LR.ind)
      LR.ind_vec <- c(LR.ind_vec,LR.ind)
      p.uc_vec <- c(p.uc_vec,1-pchisq(LR.uc,1))
      p.cc_vec <- c(p.cc_vec,1-pchisq(LR.uc+LR.ind,2))
      p.ind_vec <- c(p.ind_vec,1-pchisq(LR.ind,1))
    }
    
    names(LR.uc_vec) <- names(LR.cc_vec) <- names(LR.ind_vec) <- paste(100*(1-VaR.alpha), "%")
    names(p.uc_vec)  <- names(p.cc_vec)  <- names(p.ind_vec) <- paste(100*(1-VaR.alpha), "%")
    
    model_extern[i+j,paste('LR.uc_stat',100*(1-c(0.01,0.05,0.10)))] <- LR.uc_vec
    model_extern[i+j,paste('LR.ind_stat',100*(1-c(0.01,0.05,0.10)))] <- LR.ind_vec
    model_extern[i+j,paste('LR.cc_stat',100*(1-c(0.01,0.05,0.10)))] <- LR.cc_vec
    
    model_extern[i+j,paste('LR.uc_pvalue',100*(1-c(0.01,0.05,0.10)))] <- p.uc_vec
    model_extern[i+j,paste('LR.ind_pvalue',100*(1-c(0.01,0.05,0.10)))] <- p.ind_vec
    model_extern[i+j,paste('LR.cc_pvalue',100*(1-c(0.01,0.05,0.10)))] <- p.cc_vec
    
    model_extern[i+j,"time"] <- difftime(Sys.time(), start_time,units = "secs")[[1]]
  }
  
  write.csv(model, paste(filename,"csv",sep="."), row.names=FALSE)
  print(paste("===",round(100*i/nrow(model_extern),2) , "%" ,"===="))
  write.csv(model_extern, "Forecast_RealEGARCH.csv", row.names=FALSE)
}

#write.csv(model_extern, "Forecast_RealEGARCH_756.csv", row.names=FALSE)


#-------------------------------------------------------------------------------
# MS-RV : ESTIMATION
#-------------------------------------------------------------------------------

para_init<-function(LEVIER,K){
  para<-c(0,0.80,0.12,-0.05,5.01)
  if(K==2){
    Vol<-c(0.02,3.5)
    P<-matrix(c(0.99,0.01,0.01,0.99),K,K)
  }
  if(K==3){
    Vol<-c(1,2,3)
    P<-matrix(c(0.9,0.05,0.05,0.05,0.9,0.05,0.05,0.05,0.9),K,K)
  }
  if(K==4){
    Vol<-c(0.02,0.05,0.5,1.5)
    P<-matrix(c(0.7,0.1,0.1,0.1,0.1,0.7,0.1,0.1,0.1,0.1,0.7,0.1,0.1,0.1,0.1,0.7),K,K)
  }
  if(LEVIER) para<-c(para,2.75,0.52)
  
  return(list(para=para,P=P,Vol=Vol))
}

model<-rbind(expand.grid(index=c('ftse','dax'), start.date=as.Date("2000-01-01"), end.date=as.Date("2019-12-31"), length=0,
                   K=c(3), LEVIER=c(TRUE)),
             expand.grid(index=c('psi','dji','nifty50','omxc20'), start.date=as.Date("2000-01-01"), end.date=as.Date("2019-12-31"), length=0,
                   K=c(4), LEVIER=c(TRUE)))

model<-rbind(expand.grid(index=c('ftse'), start.date=as.Date("2000-01-01"), end.date=as.Date("2019-12-31"), length=0,
                         K=c(4), LEVIER=c(TRUE)))

vars<-c('loglikm','loglik', 'AIC', 'BIC', paste0("param",1:23),"conv","time")

model_add <- matrix(0, nrow=nrow(model), ncol=length(vars))
colnames(model_add) <- vars
model <- cbind(model, model_add)

filename <- paste("MSRV","_all_index_",start.date,"_",end.date, sep="")
sourceCpp('C:/Users/DellPC/Dropbox/Abdoul/These/Article1/Code/MDSV/Code_These/MDSV_package/MDSV/test/These/Chapitre3/benchmarks/RealizedMS.cpp')

fn<-function(para_tilde,ech,K=K,LEVIER=LEVIER)
  logLik(ech=ech,para_tilde = para_tilde,K=K,LEVIER=LEVIER, Model_type = "logN")

for(i in 1:nrow(model)){
  start_time <- Sys.time()
  index                 <- as.character(model[i,"index"])
  donne                 <- get(index)[rownames(get(index))>=start.date & rownames(get(index))<=end.date,]
  nam_                  <- rownames(donne)
  donne[,"r"]           <- donne[,"r"] - mean(donne[,"r"])
  donne                 <- cbind("rv"=donne[,"rv"], "r"=donne[,"r"])
  model[i,"start.date"] <- as.Date(nam_[1])
  model[i,"end.date"]   <- as.Date(nam_[length(nam_)])
  model[i,"length"]     <- length(nam_)
  dist                  <- as.character(model[i,"dist"])
  LEVIER                <- model[i,"LEVIER"]
  K                     <- as.numeric(model[i,"K"])
  
  para<-para_init(LEVIER,K)
  para_tilde<-natWork(para=para$para,Vol=para$Vol,P=para$P,K=K,LEVIER,"logN")
  
  LB  = rep(-10, 5+K+K*(K-1)+2*LEVIER)#, -10, -10)
  UB  = rep(10, 5+K+K*(K-1)+2*LEVIER)#, 10, 10)
  
  oldw <- getOption("warn")
  options(warn = -1)
  opt<-try(gosolnp(pars=para_tilde,fun=fn,ech=donne,LEVIER=LEVIER,K=K,control=ctrl,
                   LB=LB,UB=UB,n.restarts=10,n.sim=5000,cluster=NULL),silent=T)
  options(warn = oldw)
  
  if(is(opt,"try-error")) {
    model[i,'conv']  <- "No Convergence."
    print(paste("===",round(100*i/nrow(model),2) , "% : ERROR" ,"===="))
    next
  }
  
  model[i,colnames(model) %in% paste0("param",1:(5+K+K*(K-1)+2*LEVIER))] <- round(opt$pars,5)
  l<-logLik2(ech=donne,para_tilde=opt$pars,LEVIER=LEVIER,K=K,Model_type="logN",t=nrow(donne),r=0)
  model[i,'conv']   <- "Convergence."
  model[i,'loglik'] <- -as.numeric(opt$values[length(opt$values)])
  model[i,"loglikm"]<- l$Marg_loglik
  model[i,'AIC']    <- model[i,"loglik"] - length(opt$pars)
  model[i,'BIC']    <- model[i,"loglik"] - length(opt$pars)*log(length(nam_))/2
  model[i,"time"]   <- difftime(Sys.time(), start_time,units = "secs")[[1]]
  
  #write.csv(model, paste(filename,"csv",sep="."), row.names=FALSE)
  print(paste("===",round(100*i/nrow(model),2) , "%" ,"==== LOGLIK = ", model[i,"loglik"], "==="))
}

write.csv(model, "Estim_ALLRealMSRV1.csv", row.names=FALSE)
X<-matrix(unlist(model[,11:33]),nrow=nrow(model))

#-------------------------------------------------------------------------------
# MS-RV : PREVISION
#-------------------------------------------------------------------------------

g<-function(vector){
  Sortie<-matrix(NA,2,2)
  for(i in c(0,1)) for(j in c(0,1)) Sortie[i+1,j+1]<-sum((vector[-1]==j)*(vector[-length(vector)]==i))
  if(vector[1]==0) Sortie[1,1]<-Sortie[1,1]+1
  if(vector[1]==1) Sortie[1,2]<-Sortie[1,2]+1
  colnames(Sortie)<-c(0,1)
  rownames(Sortie)<-c(0,1)
  return(Sortie)
}

pmist2n <- function(y,sigma,p){
  return(sum(p*pnorm(y, 0, sqrt(sigma))))
}#pmist2n

qmist2n <- function(q,sigma,p){
  # the min minmax below is computed to supply a range to the solver
  # the solution must be between the min and max
  # quantile of the mixed distributions
  minmax <- range(qnorm(q,0,sqrt(sigma)))
  uniroot(function(x) pmist2n(x,sigma=sigma,p=p)-q,
          interval = minmax,
          tol = 10^{-16})$root  
}#qmist2n

para_name<-function(params,K,LEVIER){
  vars.names<-c("shape")
  if(LEVIER) vars.names<-c(vars.names,'l','theta')
  
  names(params$para)<-vars.names
  names(params$Vol)<-paste('Vol',1:K)
  return(params)
}

ctrl <- list(TOL=1e-15, trace=0)

model_extern<-expand.grid(window = c(252, 504, 756, 1008, 1260, 1512), index=index_set, K=c(4), LEVIER=c(TRUE), n.ahead=c(100),
                          forecast.length=c(1512), refit.every=c(21), refit.window=c("recursive"),
                          calculate.VaR=c(TRUE), rseed=1050)

Loss.horizon    <- c(1,5,10,25,50,75,100)
n.ahead         <- 100
forecast.length <- 1512
refit.every     <- 21
refit.window    <- "recursive"
VaR.alpha       <- c(0.01,0.05,0.10)
rseed           <- 1050

vars<-c('start.date','end.date','pred_dens',paste('QLIKc_R',Loss.horizon), paste('RMSEc_R',Loss.horizon), 
        paste('MAEc_R',Loss.horizon), paste('QLIKm_R',Loss.horizon), paste('RMSEm_R',Loss.horizon),
        paste('MAEm_R',Loss.horizon), paste('QLIKc_RV',Loss.horizon), paste('RMSEc_RV',Loss.horizon), 
        paste('MAEc_RV',Loss.horizon), paste('QLIKm_RV',Loss.horizon), paste('RMSEm_RV',Loss.horizon),
        paste('MAEm_RV',Loss.horizon), paste('LR.uc_stat',100*(1-c(0.01,0.05,0.10))), 
        paste('LR.uc_pvalue',100*(1-c(0.01,0.05,0.10))), paste('LR.ind_stat',100*(1-c(0.01,0.05,0.10))), 
        paste('LR.ind_pvalue',100*(1-c(0.01,0.05,0.10))), paste('LR.cc_stat',100*(1-c(0.01,0.05,0.10))), 
        paste('LR.cc_pvalue',100*(1-c(0.01,0.05,0.10))),"time")

model_add <- matrix(0, nrow=nrow(model_extern), ncol=length(vars))
colnames(model_add) <- vars
model_extern <- cbind(model_extern, model_add)

R_var <- rep(0,length(Loss.horizon))

library(readr)

for(i in c(1,7,13)){
  start_time <- Sys.time()
  index                        <- as.character(model_extern[i,"index"])
  donne                        <- get(index)[rownames(get(index))>=start.date & rownames(get(index))<=end.date,]
  nam_                         <- rownames(donne)
  donne[,"r"]                  <- donne[,"r"] - mean(donne[,"r"])
  donne                        <- cbind("rv"=donne[,"rv"], "r"=donne[,"r"])
  model_extern[i:(i+5),"start.date"] <- as.character(as.Date(nam_[length(nam_)])-forecast.length)
  model_extern[i:(i+5),"end.date"]   <- as.character(nam_[length(nam_)])
  dist                         <- as.character(model_extern[i,"dist"])
  LEVIER                       <- model_extern[i,"LEVIER"]
  calculate.VaR                <- as.logical(model_extern[i,"calculate.VaR"])
  K                            <- as.numeric(model_extern[i,"K"])
  
  set.seed(rseed)
  
  filename <- paste("Forecast_RealMSRV_K_",K,"_LEVIER_",LEVIER,"_", index, "_length_",forecast.length,
                    "_", model_extern[i,"start.date"], "_", model_extern[i,"end.date"], sep="")
  
  model<-expand.grid(date = nam_[((length(nam_)-forecast.length+1):length(nam_))],rt=0,rvt=0,Levier=LEVIER)
  
  model$rt  <- donne[((length(nam_)-forecast.length+1):length(nam_)),"r"]
  model$rvt  <- donne[((length(nam_)-forecast.length+1):length(nam_)),"rv"]
  
  vars<-c(paste0("R_for",Loss.horizon),paste0("R_tru",Loss.horizon),
          paste0("R_for_m",Loss.horizon),paste0("R_tru_m",Loss.horizon),
          paste0("RV_for",Loss.horizon),paste0("RV_tru",Loss.horizon),
          paste0("RV_for_m",Loss.horizon),paste0("RV_tru_m",Loss.horizon),'pred_lik','loglik', 'AIC', 'BIC',
          paste0("param",1:23), paste0("VaR",100*(1-VaR.alpha)), 
          paste0("I",100*(1-VaR.alpha)))
  
  model_add <- matrix(0, nrow=nrow(model), ncol=length(vars))
  colnames(model_add) <- vars
  model <- cbind(model, model_add)
  X <- read_csv("MSRV_Estim.csv")
  
  #para<-para_init(LEVIER,K)
  #para_tilde<-natWork(para=para$para,Vol=para$Vol,P=para$P,K=K,LEVIER,"logN")
  # para_tilde <- X[i,]
  X <- unlist(X[(X$index=='stoxx50'),6:ncol(X)])
  # X <- c(-0.42403,	0.97708, -0.10733,	0.1151,	-1.51787,	-0.44547,
  #        -1.76274,	-1.4538,	0.30502,	-2.32556,	-0.70217,	3.53387,
  #        -9.41913,	-0.10517,	-10,	3.02247,	-10,	3.11789,	-9.62701,
  #        0.104503, 0.0225187, 3.343334, 0.238355)
  
  # X <- c(0.92239,1.66738,-0.90161,-2.0886,0.26459,-1.51787,	-0.44547,
  #        2.15614, -10,-10,-3.62345,3.15387,0.00935,-9.59061, 4.50915,
  #        6.46873,-3.63302,-2.94299,-10, 3.343334, 0.238355, 0.104503, 0.0225187)
  
  LB  = rep(-10, 5+K+K*(K-1)+2*LEVIER)#, -10, -10)
  UB  = rep(10, 5+K+K*(K-1)+2*LEVIER)#, 10, 10)
  opt             <- NULL
  update_date     <- seq(0,forecast.length,by=refit.every)
  strt            <- 1
  
  for(t in 0:(forecast.length-1)){
  # for(t in 441:(441+62)){
    ech    <- donne[strt:(length(nam_)-forecast.length+t),]
    
    if(t %in% update_date){
      
      if((!is.null(opt)) & !(is(opt,"try-error"))) para_tilde<-opt$pars
      oldw <- getOption("warn")
      options(warn = -1)
      # opt1 <- try(gosolnp(pars=para_tilde,fun=fn,ech=ech,LEVIER=LEVIER,K=K,control=ctrl,
      #                           LB=LB,UB=UB,n.restarts=10,n.sim=2000,cluster=NULL),silent=T)
      opt1 <- NULL
      opt2 <- try(solnp(pars=X,fun=fn,ech=ech,LEVIER=LEVIER,K=K,control=ctrl),silent=T)

      if((!is.null(opt1)) & !(is(opt1,"try-error"))){
        if((!is.null(opt2)) & !(is(opt2,"try-error"))){
          if((-as.numeric(opt1$values[length(opt1$values)])) < (-as.numeric(opt2$values[length(opt2$values)]))){
            opt<-opt2
          }else{
            opt<-opt1
          }
        }else{opt<-opt1}
      }else{
        if((!is.null(opt2)) & !(is(opt2,"try-error"))) opt<-opt2
      }
        
      options(warn = oldw)
      
      if(refit.window == "moving") strt <- strt + refit.every
      
    }
    
    
    model[t+1,colnames(model) %in% paste0("param",1:(5+K+K*(K-1)+2*LEVIER))] <- round(opt$pars,5)
    l<-logLik2(ech=ech,para_tilde=opt$pars,LEVIER=LEVIER,K=K,Model_type="logN",t=nrow(ech),
               r=donne[(length(nam_)-forecast.length+t+1),"r"])
  
    
    model[t+1,'loglik']   <- -as.numeric(opt$values[length(opt$values)])
    model[t+1,"loglikm"]  <- l$Marg_loglik
    model[t+1,"pred_lik"] <- l$Pred_loglik
    model[t+1,'AIC']      <- model[t+1,"loglik"] - length(opt$pars)
    model[t+1,'BIC']      <- model[t+1,"loglik"] - length(opt$pars)*log(nrow(ech))/2
    
    n <- 10000 #number of simulations
    H <- max(Loss.horizon)    #nb of years simulated
    
    para<-workNat(opt$pars,LEVIER,K,"logN")
    
    sig<- para$Vol
    pi_0 <- l$w_hat
    if(LEVIER){
      Levier<-levierVolatility(ech[(length(ech[,"r"])-500):(length(ech[,"r"])),"r"],Nl=70,para=para$para,Model_type="logN")$`levier`
      sig<-sig*Levier
    }
    
    if(calculate.VaR) for(iter in 1:length(VaR.alpha)){
      va <- try(qmist2n(VaR.alpha[iter], sigma=sig, p=pi_0),silent=T)
      if(class(va) =='try-error') {
        va <- try(qmixnorm(VaR.alpha[iter], rep(0,length(sig)), sqrt(sig), pi_0),silent=T)
        if(!(class(va) =='try-error')) model[t+1,paste0('VaR',100*(1-VaR.alpha[iter]))] <- va
      }else{
        model[t+1,paste0('VaR',100*(1-VaR.alpha[iter]))] <- va
      }
    }
    
    matP<-para$P
    sig<-para$Vol
    para<-para$para
    #if(para[7]>0.90)  para <- c(para[-c(7)], 0.8838735)
    
    if(LEVIER){
      MC_sim <- t(matrix(sim.mc(pi_0, matP, rep(H, n)),H,n,byrow=FALSE)) #simulation of Markov chain
      z_t<-matrix(rnorm(n*H),nrow=n,ncol=H)
      n_t<-matrix(rinvgamma(n*H,shape = para[5]+1,scale = para[5]),nrow=n,ncol=H)
      
      Levier<-rep(1,n)%*%t(levierVolatility(ech[(length(ech[,"r"])-500):(length(ech[,"r"])),"r"],Nl=70,para,"logN")$`Levier`)
      sim<-R_hat(H,ech[(length(ech[,"r"])-500):(length(ech[,"r"])),"r"],MC_sim,z_t,n_t=n_t,Levier,sig,para,K,"logN",Nl=70)
      rt2_sim<-sim$`rt2`
      rvt_sim<-sim$`rvt`
      rt2_sim <- rt2_sim[,(ncol(rt2_sim)-H+1):ncol(rt2_sim)]
      rvt_sim <- rvt_sim[,(ncol(rvt_sim)-H+1):ncol(rvt_sim)]
    }else {
      sim<-f_sim_logN(H,sig,pi_0,matP,varphi=para[2],xi=para[1],shape=para[5],delta1=para[3],delta2=para[4])
      rt2_sim<-sim$`rt2`
      rvt_sim<-sim$`rvt`
    }
    
    for(k in 1:length(Loss.horizon)) R_var[k] <- sum(rt2_sim[1:Loss.horizon[k]])
    names(R_var) <- paste0("R_for",Loss.horizon)
    model[(t+1),colnames(model) %in% names(R_var)]<-R_var
    
    for(k in 1:length(Loss.horizon)) R_var[k] <- sum(rt2_sim[Loss.horizon[k]])
    names(R_var) <- paste0("R_for_m",Loss.horizon)
    model[(t+1),colnames(model) %in% names(R_var)]<-R_var
    
    ech    <- donne[(length(nam_)-forecast.length+t+1):length(nam_),"r"]
    for(k in 1:length(Loss.horizon)) R_var[k] <- sum((ech[1:Loss.horizon[k]])^2)
    names(R_var) <- paste0("R_tru",Loss.horizon)
    model[(t+1),colnames(model) %in% names(R_var)]<-R_var
    
    for(k in 1:length(Loss.horizon)) R_var[k] <- sum((ech[Loss.horizon[k]])^2)
    names(R_var) <- paste0("R_tru_m",Loss.horizon)
    model[(t+1),colnames(model) %in% names(R_var)]<-R_var
    
    for(k in 1:length(Loss.horizon)) R_var[k] <- sum(rvt_sim[1:Loss.horizon[k]])
    names(R_var) <- paste0("RV_for",Loss.horizon)
    model[(t+1),colnames(model) %in% names(R_var)]<-R_var
    
    for(k in 1:length(Loss.horizon)) R_var[k] <- sum(rvt_sim[Loss.horizon[k]])
    names(R_var) <- paste0("RV_for_m",Loss.horizon)
    model[(t+1),colnames(model) %in% names(R_var)]<-R_var
    
    ech    <- donne[(length(nam_)-forecast.length+t+1):length(nam_),"rv"]
    for(k in 1:length(Loss.horizon)) R_var[k] <- sum(ech[1:Loss.horizon[k]])
    names(R_var) <- paste0("RV_tru",Loss.horizon)
    model[(t+1),colnames(model) %in% names(R_var)]<-R_var
    
    for(k in 1:length(Loss.horizon)) R_var[k] <- sum(ech[Loss.horizon[k]])
    names(R_var) <- paste0("RV_tru_m",Loss.horizon)
    model[(t+1),colnames(model) %in% names(R_var)]<-R_var
    
    print(paste0("========== i = ",i, " et pourcentage = ",100*(t+1)/756," % ==============="))
    
  }
  
  if(calculate.VaR)  for(iter in 1:length(VaR.alpha)){
    model[,paste0('I',100*(1-VaR.alpha[iter]))] <- (model[,'rt'] < model[,paste0('VaR',100*(1-VaR.alpha[iter]))])
  }
  
  for(j in 0:5){
    window = c(252, 504, 756, 1008, 1260, 1512)[j+1]
    
    mod <- model[(1512-window+1):1512,]
    model_extern[i+j,'pred_dens'] <- sum(mod[,"pred_lik"])
    
    for_R_var   <- mod[,grep('R_for', colnames(mod), fixed=TRUE)]
    for_R_var   <- for_R_var[,-grep('R_for_m', colnames(for_R_var), fixed=TRUE)]
    for_R_err   <- mod[,grep('R_tru', colnames(mod), fixed=TRUE)]
    for_R_err   <- for_R_err[,-grep('R_tru_m', colnames(for_R_err), fixed=TRUE)]
    model_extern[i+j,paste('QLIKc_R',Loss.horizon)] <- colMeans(log(for_R_var) + for_R_err/for_R_var, na.rm=TRUE)
    model_extern[i+j,paste('RMSEc_R',Loss.horizon)] <- sqrt(colMeans( (for_R_var - for_R_err)^2, na.rm=TRUE ))/Loss.horizon
    model_extern[i+j,paste('MAEc_R',Loss.horizon)] <- colMeans( abs(for_R_var - for_R_err), na.rm=TRUE )/Loss.horizon
    
    for_R_var   <- mod[,grep('R_for_m', colnames(mod), fixed=TRUE)]
    for_R_err   <- mod[,grep('R_tru_m', colnames(mod), fixed=TRUE)]
    model_extern[i+j,paste('QLIKm_R',Loss.horizon)] <- colMeans(log(for_R_var) + for_R_err/for_R_var, na.rm=TRUE)
    model_extern[i+j,paste('RMSEm_R',Loss.horizon)] <- sqrt(colMeans( (for_R_var - for_R_err)^2, na.rm=TRUE ))/Loss.horizon
    model_extern[i+j,paste('MAEm_R',Loss.horizon)]  <- colMeans( abs(for_R_var - for_R_err), na.rm=TRUE )/Loss.horizon
    
    for_R_var   <- mod[,grep('RV_for', colnames(mod), fixed=TRUE)]
    for_R_var   <- for_R_var[,-grep('RV_for_m', colnames(for_R_var), fixed=TRUE)]
    for_R_err   <- mod[,grep('RV_tru', colnames(mod), fixed=TRUE)]
    for_R_err   <- for_R_err[,-grep('RV_tru_m', colnames(for_R_err), fixed=TRUE)]
    model_extern[i+j,paste('QLIKc_RV',Loss.horizon)] <- colMeans(log(for_R_var) + for_R_err/for_R_var, na.rm=TRUE)
    model_extern[i+j,paste('RMSEc_RV',Loss.horizon)] <- sqrt(colMeans( (for_R_var - for_R_err)^2, na.rm=TRUE ))/Loss.horizon
    model_extern[i+j,paste('MAEc_RV',Loss.horizon)] <- colMeans( abs(for_R_var - for_R_err), na.rm=TRUE )/Loss.horizon
    
    for_R_var   <- mod[,grep('RV_for_m', colnames(mod), fixed=TRUE)]
    for_R_err   <- mod[,grep('RV_tru_m', colnames(mod), fixed=TRUE)]
    model_extern[i+j,paste('QLIKm_RV',Loss.horizon)] <- colMeans(log(for_R_var) + for_R_err/for_R_var, na.rm=TRUE)
    model_extern[i+j,paste('RMSEm_RV',Loss.horizon)] <- sqrt(colMeans( (for_R_var - for_R_err)^2, na.rm=TRUE ))/Loss.horizon
    model_extern[i+j,paste('MAEm_RV',Loss.horizon)]  <- colMeans( abs(for_R_var - for_R_err), na.rm=TRUE )/Loss.horizon
    
    LR.uc_vec <- LR.cc_vec <- LR.ind_vec <- NULL
    p.uc_vec  <- p.cc_vec  <- p.ind_vec  <- NULL
    
    for(iter in 1:length(VaR.alpha)){
      viol      <- sum(mod[,paste0('I',100*(1-VaR.alpha[iter]))])
      alpha_hat <- sum(mod[,paste0('I',100*(1-VaR.alpha[iter]))])/forecast.length
      LR.uc     <- 2*log(((alpha_hat^viol)*((1-alpha_hat)^(forecast.length-viol)))/((VaR.alpha[iter]^viol)*((1-VaR.alpha[iter])^(forecast.length-viol))))
      
      viol      <- g(mod[,paste0('I',100*(1-VaR.alpha[iter]))])
      pi        <- (viol[1,2]+viol[2,2])/sum(viol)
      pi0       <- (viol[1,2])/(viol[1,1]+viol[1,2])
      pi1       <- (viol[2,2])/(viol[2,2]+viol[2,1])
      LR.ind    <- - 2*log((((1-pi)^(viol[1,1]+viol[2,1]))*(pi^(viol[1,2]+viol[2,2])))/
                             ((((1-pi0)^viol[1,1])*(pi0^viol[1,2]))*(((1-pi1)^viol[2,1])*(pi1^viol[2,2]))))
      
      LR.uc_vec <- c(LR.uc_vec,LR.uc)
      LR.cc_vec <- c(LR.cc_vec,LR.uc+LR.ind)
      LR.ind_vec <- c(LR.ind_vec,LR.ind)
      p.uc_vec <- c(p.uc_vec,1-pchisq(LR.uc,1))
      p.cc_vec <- c(p.cc_vec,1-pchisq(LR.uc+LR.ind,2))
      p.ind_vec <- c(p.ind_vec,1-pchisq(LR.ind,1))
    }
    
    names(LR.uc_vec) <- names(LR.cc_vec) <- names(LR.ind_vec) <- paste(100*(1-VaR.alpha), "%")
    names(p.uc_vec)  <- names(p.cc_vec)  <- names(p.ind_vec) <- paste(100*(1-VaR.alpha), "%")
    
    model_extern[i+j,paste('LR.uc_stat',100*(1-c(0.01,0.05,0.10)))] <- LR.uc_vec
    model_extern[i+j,paste('LR.ind_stat',100*(1-c(0.01,0.05,0.10)))] <- LR.ind_vec
    model_extern[i+j,paste('LR.cc_stat',100*(1-c(0.01,0.05,0.10)))] <- LR.cc_vec
    
    model_extern[i+j,paste('LR.uc_pvalue',100*(1-c(0.01,0.05,0.10)))] <- p.uc_vec
    model_extern[i+j,paste('LR.ind_pvalue',100*(1-c(0.01,0.05,0.10)))] <- p.ind_vec
    model_extern[i+j,paste('LR.cc_pvalue',100*(1-c(0.01,0.05,0.10)))] <- p.cc_vec
    
    model_extern[i+j,"time"] <- difftime(Sys.time(), start_time,units = "secs")[[1]]
  }
  
  #model_extern[,c("RMSEc_RV 1","RMSEc_RV 5","RMSEc_RV 25","RMSEc_RV 100")]
  
  write.csv(model, paste(filename,"csv",sep="."), row.names=FALSE)
  print(paste("===",round(100*i/nrow(model_extern),2) , "%" ,"===="))
  write.csv(model_extern, "Forecast_RealMSRV.csv", row.names=FALSE)
}

write.csv(model_extern, "Forecast_RealMSRV_756.csv", row.names=FALSE)

