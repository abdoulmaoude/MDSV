path<-"C:/Users/DellPC/Dropbox/Abdoul/These/Article1/Code/MDSV/Code_These/Ess"
path<-"/home/maoudek/Rsim/Article1/MDSV/Forecast"
setwd(path)
library(MDSV)
if(!require(parallel)){install.packages("parallel")}; library(parallel)

index_set<-c("aex", "aord", "bell", "bsesn", "bvsp", "cac40", "dax", "dji", "ftmib", "ftse",
             "hsi", "ibex35", "kospi", "kse", "mxx", "n225", "nasdaq", "nifty50", "omxc20", "omxhpi",
             "omxspi", "oseax", "psi", "rut", "smsi", "sp500", "ssec", "ssmi", "sti", "stoxx50", "tsx")

index_set<-c("sp500", "nasdaq", "ftse", "n225")

start.date <- as.Date("2000-01-01") 
end.date   <- as.Date("2019-12-31")

#-------------------------------------------------------------------------------
# Mise en place
#-------------------------------------------------------------------------------

filenam <- paste("Forecast_MDSV_RandJ_Good_ALL_756_", start.date, "_",end.date, sep="")

# model_extern<-rbind(expand.grid(index=c("sp500"), K=c(2), N=c(3),
#                                 ModelType=c(0,2), LEVIER=c(FALSE),n.ahead=c(100),forecast.length=c(756),refit.every=c(63),
#                                 refit.window=c("recursive"),calculate.VaR=c(TRUE), rseed=1050))

model_extern<-rbind(expand.grid(index=index_set, K=c(2), N=c(10), LEVIER=c(FALSE,TRUE), ModelType=c(0,2),
                          n.ahead=c(100),forecast.length=c(756),refit.every=c(63),
                          refit.window=c("recursive"),calculate.VaR=c(TRUE), rseed=1050),
                    expand.grid(index=index_set, K=c(10), N=c(3), LEVIER=c(FALSE,TRUE), ModelType=c(0,2),
                          n.ahead=c(100),forecast.length=c(756),refit.every=c(63),
                          refit.window=c("recursive"),calculate.VaR=c(TRUE), rseed=1050),
                    expand.grid(index=index_set, K=c(3), N=c(6), LEVIER=c(FALSE,TRUE), ModelType=c(0,2),
                          n.ahead=c(100),forecast.length=c(756),refit.every=c(63),
                          refit.window=c("recursive"),calculate.VaR=c(TRUE), rseed=1050))

model_extern$nstate <- (model_extern$K)^(model_extern$N)

Loss.horizon <- c(1,5,10,25,50,75,100)

vars<-c('start.date','end.date','pred_dens',paste('QLIKc_R',Loss.horizon), paste('RMSEc_R',Loss.horizon), 
        paste('MAEc_R',Loss.horizon), paste('QLIKm_R',Loss.horizon), paste('RMSEm_R',Loss.horizon),
        paste('MAEm_R',Loss.horizon),paste('QLIKc_RV',Loss.horizon), paste('RMSEc_RV',Loss.horizon), 
        paste('MAEc_RV',Loss.horizon), paste('QLIKm_RV',Loss.horizon), paste('RMSEm_RV',Loss.horizon), 
        paste('MAEm_RV',Loss.horizon), paste('LR.uc_stat',100*(1-c(0.01,0.05,0.10))), 
        paste('LR.uc_pvalue',100*(1-c(0.01,0.05,0.10))), paste('LR.ind_stat',100*(1-c(0.01,0.05,0.10))), 
        paste('LR.ind_pvalue',100*(1-c(0.01,0.05,0.10))), paste('LR.cc_stat',100*(1-c(0.01,0.05,0.10))), 
        paste('LR.cc_pvalue',100*(1-c(0.01,0.05,0.10))),"time")

model_add <- matrix(0, nrow=nrow(model_extern), ncol=length(vars))
colnames(model_add) <- vars
model_extern <- cbind(model_extern, model_add)

model_extern<-model_extern[-c(4,7,23,26,39,40,11,26,41,42),]

ctrl <- list(TOL=1e-15, trace=0)

n.ahead         <- 100
forecast.length <- 756
refit.every     <- 63
refit.window    <- "recursive"
VaR.alpha       <- c(0.01,0.05,0.10)
rseed           <- 1050

for(i in 1:nrow(model_extern)){
  #8, 9, 18, 21 || (jaune) 2, 3, 7, 11, 12, 15, 25, 35, 37 
#for(i in c(1,237,473,709,945,1181)){
  start_time      <- Sys.time()
  index           <- as.character(model_extern[i,"index"])
  donne           <- get(index)[rownames(get(index))>=start.date & rownames(get(index))<=end.date,]
  nam_            <- rownames(donne)
  donne[,"r"]     <- donne[,"r"] - mean(donne[,"r"])
  ModelType       <- as.numeric(model_extern[i,"ModelType"])
  LEVIER          <- as.logical(model_extern[i,"LEVIER"])
  N               <- as.numeric(model_extern[i,"N"])
  K               <- as.numeric(model_extern[i,"K"])
  calculate.VaR   <- as.logical(model_extern[i,"calculate.VaR"])

  # tmp<-read.csv("All_estimations.csv")
  # if(ModelType==0) {
  #   modeltype <- "log-return"
  # }else{
  #   modeltype <- "Joint"
  # }
  # 
  # ij<-(tmp$index==index) & (tmp$ModelType == modeltype) & (tmp$LEVIER == LEVIER) & (tmp$N == N) & (tmp$K == K)
  # 
  # vars<-c("omega","a","b","sigma","v0")
  # if(ModelType==2) {
  #   vars <- c(vars,"xi","varphi","delta1","delta2","shape")
  # }
  # if(LEVIER) {
  #   vars <- c(vars,"l","theta")
  # }
  # 
  # start.pars<-unlist(tmp[ij,vars])
  
  start.pars <- c(0.25,0.87, 2.77, sqrt(var(donne[,"r"])),0.33)
  if(ModelType==1) start.pars <- c(start.pars,2.10)
  if(ModelType==2) start.pars <- c(start.pars,-1.5,	0.72,	-0.09,	0.04,	2.10)
  if(LEVIER)       start.pars <- c(start.pars,1.5,0.87568)
  
  cluster         <- makeCluster(6)
  out<-MDSVroll(N=N, K=K, data=donne, ModelType=ModelType, LEVIER=LEVIER, n.ahead = n.ahead, 
                forecast.length = forecast.length, refit.every = refit.every, refit.window = refit.window, 
                window.size=NULL,calculate.VaR = calculate.VaR, VaR.alpha = VaR.alpha, cluster = cluster, 
                rseed = rseed,start.pars=start.pars)
  parallel::stopCluster(cluster)
  
  out<-print(out, VaR.test=calculate.VaR, Loss.horizon = c(1,5,10,25,50,75,100), Loss.window = 756)
  
  model_extern[i,"start.date"] <- as.character(as.Date(nam_[length(nam_)])-forecast.length)
  model_extern[i,"end.date"]   <- as.character(nam_[length(nam_)])
  
  Model<-cbind(out$Loss[,1:7],predict_loglik=out$estimates[,"predict_loglik"],out$Loss[,8:ncol(out$Loss)],
               out$estimates[,grep('VaR', colnames(out$estimates), fixed=TRUE)],
               out$estimates[,grep('I', colnames(out$estimates), fixed=TRUE)])
  
  if(ModelType==0) model_type <- "R"
  if(ModelType==1) model_type <- "RV"
  if(ModelType==2) model_type <- "J"
  
  filename <- paste("Forecast_MDSV_",model_type,"_N_",N,"_K_",K,"_LEVIER_",LEVIER,"_length_",forecast.length,
                    "_", index,"_", model_extern[i,"start.date"], "_", model_extern[i,"end.date"], sep="")

  write.csv(Model, paste0(filename,".csv"), row.names=FALSE)
  
  if(!(ModelType==1)){
    model_extern[i,'pred_dens'] <- sum(out$estimates[,"predict_loglik"])
    model_extern[i,paste('QLIKc_R',Loss.horizon)] <- out$Loss_functions_R["QLIK",]
    model_extern[i,paste('RMSEc_R',Loss.horizon)] <- out$Loss_functions_R["RMSE",]
    model_extern[i,paste('MAEc_R',Loss.horizon)] <- out$Loss_functions_R["MAE",]
    
    model_extern[i,paste('QLIKm_R',Loss.horizon)] <- out$Loss_functions_marg_R["QLIK",]
    model_extern[i,paste('RMSEm_R',Loss.horizon)] <- out$Loss_functions_marg_R["RMSE",]
    model_extern[i,paste('MAEm_R',Loss.horizon)] <- out$Loss_functions_marg_R["MAE",]
    
    model_extern[i,paste('LR.uc_stat',100*(1-c(0.01,0.05,0.10)))] <- out$LR.uc
    model_extern[i,paste('LR.ind_stat',100*(1-c(0.01,0.05,0.10)))] <- out$LR.ind
    model_extern[i,paste('LR.cc_stat',100*(1-c(0.01,0.05,0.10)))] <- out$LR.cc
    
    model_extern[i,paste('LR.uc_pvalue',100*(1-c(0.01,0.05,0.10)))] <- out$p.uc
    model_extern[i,paste('LR.ind_pvalue',100*(1-c(0.01,0.05,0.10)))] <- out$p.ind
    model_extern[i,paste('LR.cc_pvalue',100*(1-c(0.01,0.05,0.10)))] <- out$p.cc
  }
  
  if(!(ModelType==0)){
    model_extern[i,paste('QLIKc_RV',Loss.horizon)] <- out$Loss_functions_RV["QLIK",]
    model_extern[i,paste('RMSEc_RV',Loss.horizon)] <- out$Loss_functions_RV["RMSE",]
    model_extern[i,paste('MAEc_RV',Loss.horizon)] <- out$Loss_functions_RV["MAE",]
    
    model_extern[i,paste('QLIKm_RV',Loss.horizon)] <- out$Loss_functions_marg_RV["QLIK",]
    model_extern[i,paste('RMSEm_RV',Loss.horizon)] <- out$Loss_functions_marg_RV["RMSE",]
    model_extern[i,paste('MAEm_RV',Loss.horizon)] <- out$Loss_functions_marg_RV["MAE",]
  }
  
  model_extern[i,"time"] <- difftime(Sys.time(), start_time,units = "secs")[[1]]
  
  if(file.exists(paste0(filenam,".csv"))){
    tmp <- read.csv(paste0(filenam,".csv"))
    tmp[i,13:ncol(tmp)] <- model_extern[i,13:ncol(model_extern)]
    write.csv(tmp, paste0(filenam,".csv"), row.names=FALSE)
  }else{
    write.csv(model_extern, paste0(filenam,".csv"), row.names=FALSE)
  }
  
  # print(paste("===",round(100*i/nrow(model_extern),2) , "%" ,"===="))
}

