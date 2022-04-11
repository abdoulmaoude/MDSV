path<-"C:/Users/DellPC/Dropbox/Abdoul/These/Article1/Code/MDSV/Code_These/Ess"
setwd(path)

library(MDSV)
if(!require(timeDate)){install.packages("timeDate")}; library(timeDate)

index_set<-c("aex", "aord", "bell", "bsesn", "bvsp", "cac40", "dax", "dji", "ftmib", "ftse",
             "hsi", "ibex35", "kospi", "kse", "mxx", "n225", "nasdaq", "nifty50", "omxc20", "omxhpi",
             "omxspi", "oseax", "psi", "rut", "smsi", "sp500", "ssec", "ssmi", "sti", "stoxx50", "tsx")

start.date <- as.Date("2000-01-01") 
end.date   <- as.Date("2019-12-31")

#-------------------------------------------------------------------------------
# Descriptive statitics
#-------------------------------------------------------------------------------

des_stats<-NULL
for(index in index_set){
  donne<-get(index)[rownames(get(index))>=start.date & rownames(get(index))<=end.date,]
  ech<-donne[,"r"]
  ech<-ech-mean(ech)
  q<-quantile(ech)
  #descriptive statistics
  des_stat <-    c("N" = length(ech),
                   "Min." = q[[1]],
                   "Q1" = q[[2]],
                   "Mean" = mean(ech),
                   "Median" = q[[3]],
                   "StDev" = sd(ech),
                   "Skewness" = timeDate::skewness(ech, method="moment"),
                   "Kurtosis" = timeDate::kurtosis(ech, method="moment"),
                   "Q3" = q[[4]],
                   "Max." = q[[5]]
  )
  des_stats<-rbind(des_stats,des_stat)
  ech<-donne[,"rv"]
  q<-quantile(ech)
  #descriptive statistics
  des_stat <-    c("N" = length(ech),
                   "Min." = q[[1]],
                   "Q1" = q[[2]],
                   "Mean" = mean(ech),
                   "Median" = q[[3]],
                   "StDev" = sd(ech),
                   "Skewness" = timeDate::skewness(ech, method="moment"),
                   "Kurtosis" = timeDate::kurtosis(ech, method="moment"),
                   "Q3" = q[[4]],
                   "Max." = q[[5]]
  )
  des_stats<-rbind(des_stats,des_stat)
}

rownames(des_stats)<-rep(index_set,1,each=2)
write.csv(round(des_stats,3), "stat.csv")


View(round(des_stats,2))

#-------------------------------------------------------------------------------
# Log return and realized variance graphs
#-------------------------------------------------------------------------------
pdf( file ="RV.pdf", width =6.5 , height =4)
index_set<-c("S&P 500","NASDAQ 100", "FTSE 100", "NIKKEI 225")
par(mfrow=c(2,2),mar=c(2.5, 2.5, 2, 2))
for(index in c("sp500","nasdaq","ftse","n225")){
  donne<-get(index)[rownames(get(index))>=start.date & rownames(get(index))<=end.date,]
  nam_<-rownames(donne)
  ech<-donne[,"rv"]
  
  # Traitement des dates
  t_ <- as.POSIXlt(c(as.character(start.date),nam_), format ="%Y -%m -%d")
  n.days <- rep (365,length (t_))
  n.days [((1900+ t_$ year) %% 4)==0] <- 366
  t_ <- 1900 + t_$ year + (t_$ yday +1)/n.days
  
  x.start <- t_[1] # date de dÃ©but pour les figures
  i.start <- which.min(abs(t_-x.start ))
  x.end <- t_[length (t_)] # date de fin pour les figures
  i.end <- which.min(abs(t_-x.end))
  t_<- t_[(i.start+1):i.end]
  donne <- donne [i.start:(i.end-1),]
  
  index<-index_set[which(index == c("sp500","nasdaq","ftse","n225"))]
  
  t.at <- seq(from = ceiling(x.start ),to=floor(x.end),by =2)
  xlim <- range(t_)
  cex.axis <- 0.85
  cex.mtext <- 0.85
  
  plot(t_,donne[,"rv"],type="l",xlab="", ylab="", xaxt ="n",main="", xlim =xlim , xaxs ="i")
  axis (side =1, cex.axis=cex.axis, padj=-0.6,at=t.at)
  mtext (index, side =3, line=0.5, font=2, las=0, cex=cex.mtext)
  
}
dev.off ()

#-------------------------------------------------------------------------------
# Nevolution, Nprevision, Kevolution, Kprevision
#-------------------------------------------------------------------------------

filename  <- "All_estimations.csv"
index     <- "sp500"
LEVIER    <- FALSE
ModelType <- "Joint"
start.date <- as.Date("2000-01-01") 
end.date   <- as.Date("2019-12-31")

# Nevolution
out  <- read.csv(filename)
out  <- out[((out$index==index) & (out$LEVIER==LEVIER) & (out$ModelType=="Joint")),]
out  <- out[order(out$N),]

pdf( file ="Nevolution.pdf", width =6.5 , height =4)
t.at <- seq(from = 0,to=11,by =1)
xlim <- range(0:11)
cex.axis <- 0.85
cex.mtext <- 0.85

plot(x=out$N[out$K==2],y=out$loglik[out$K==2],type='l',xlab="", ylab="", xaxt ="n",main="",xlim=xlim , xaxs ="i",lty=1)
lines(x=out$N[out$K==3],y=out$loglik[out$K==3],lty=2)
lines(x=out$N[out$K==4],y=out$loglik[out$K==4],lty=3)
lines(x=out$N[out$K==5],y=out$loglik[out$K==5],lty=4)
lines(x=out$N[out$K==6],y=out$loglik[out$K==6],lty=5)
legend("bottomright",legend=c("K = 2","K = 3","K = 4","K = 5","K = 6"),lty=1:5)
axis (side =1, cex.axis=cex.axis, padj=-0.6,at=t.at)
dev.off ()


# Kevolution
out  <- read.csv(filename)
out  <- out[((out$index==index) & (out$LEVIER==LEVIER) & (out$N==1) & (out$ModelType=="Joint")),]
out  <- out[order(out$K),]

pdf( file ="Kevolution.pdf", width =6.5 , height =4)
t.at <- seq(from = 0,to=50,by =5)
xlim <- range(0:50)
cex.axis <- 0.85
cex.mtext <- 0.85

plot(x=out$K,y=out$loglik,type='l',xlab="", ylab="", xaxt ="n",main="",xlim=xlim , xaxs ="i")
axis (side =1, cex.axis=cex.axis, padj=-0.6,at=t.at)
dev.off ()


#Nprevision

filename  <- "Forecast_MDSV_ALL_756_2_2000-01-01_2019-12-31.csv"
out  <- read.csv(filename)
out  <- out[((out$index==index) & (out$LEVIER==LEVIER) & (out$ModelType==2)),]

pdf( file ="Nprevision.pdf", width =6.5 , height =4)
par(mfrow=c(2,2),mar=c(2.5, 2.5, 2, 2))
t.at <- seq(from = 1,to=11,by =1)
xlim <- range(1:11)
cex.axis <- 0.85
cex.mtext <- 0.85

plot(x=out$N,y=out$pred_dens,type='l',xlab="", ylab="", xaxt ="n",main="", xlim=xlim , xaxs ="i")
axis (side =1, cex.axis=cex.axis, padj=-0.6,at=t.at)
mtext ("Densité prédictive cumulative", side =3, line=0.5, font=2, las=0, cex=cex.mtext)

plot(x=out$N,y=out$QLIKc_RV.1,type='l',xlab="", ylab="", xaxt ="n",main="", xlim=xlim , xaxs ="i")
axis (side =1, cex.axis=cex.axis, padj=-0.6,at=t.at)
mtext ("QLIK à 1 pas de temps", side =3, line=0.5, font=2, las=0, cex=cex.mtext)

plot(x=out$N,y=out$RMSEc_RV.1,type='l',xlab="", ylab="", xaxt ="n",main="", xlim=xlim , xaxs ="i")
axis (side =1, cex.axis=cex.axis, padj=-0.6,at=t.at)
mtext ("RMSFE à 1 pas de temps", side =3, line=0.5, font=2, las=0, cex=cex.mtext)

plot(x=out$N,y=out$LR.cc_pvalue.95,type='l',xlab="", ylab="", xaxt ="n",main="", xlim=xlim , xaxs ="i")
axis (side =1, cex.axis=cex.axis, padj=-0.6,at=t.at)
mtext ("Valeur-p du test de couverture conditionnelle", side =3, line=0.5, font=2, las=0, cex=cex.mtext)

dev.off ()


#Kprevision

filename  <- "Forecast_MDSV_Selection_length_756_2000-01-01_2019-12-31.csv"
out  <- read.csv(filename)
out  <- out[((out$index==index) & (out$LEVIER==LEVIER) & (out$N==1) & (out$ModelType==2)),]
out  <- out[out$K<100,]
out  <- out[order(out$K),]

pdf( file ="Kprevision.pdf", width =6.5 , height =4)
par(mfrow=c(2,2),mar=c(2.5, 2.5, 2, 2))
t.at <- seq(from = 0,to=50,by =5)
xlim <- range(0:50)
cex.axis <- 0.85
cex.mtext <- 0.85

plot(x=out$K,y=out$pred_dens,type='l',xlab="", ylab="", xaxt ="n",main="", xlim=xlim , xaxs ="i")
axis (side =1, cex.axis=cex.axis, padj=-0.6,at=t.at)
mtext ("Densité prédictive cumulative", side =3, line=0.5, font=2, las=0, cex=cex.mtext)

plot(x=out$K,y=out$QLIKc_RV.1,type='l',xlab="", ylab="", xaxt ="n",main="", xlim=xlim , xaxs ="i")
axis (side =1, cex.axis=cex.axis, padj=-0.6,at=t.at)
mtext ("QLIK à 1 pas de temps", side =3, line=0.5, font=2, las=0, cex=cex.mtext)

plot(x=out$K,y=out$RMSEc_RV.1,type='l',xlab="", ylab="", xaxt ="n",main="", xlim=xlim , xaxs ="i")
axis (side =1, cex.axis=cex.axis, padj=-0.6,at=t.at)
mtext ("RMSFE à 1 pas de temps", side =3, line=0.5, font=2, las=0, cex=cex.mtext)

plot(x=out$K,y=out$LR.cc_pvalue.95,type='l',xlab="", ylab="", xaxt ="n",main="", xlim=xlim , xaxs ="i")
axis (side =1, cex.axis=cex.axis, padj=-0.6,at=t.at)
mtext ("Valeur-p du test de couverture conditionnelle", side =3, line=0.5, font=2, las=0, cex=cex.mtext)

dev.off ()

#-------------------------------------------------------------------------------
# Valeur A risque
#-------------------------------------------------------------------------------

path<-"C:/Users/DellPC/Dropbox/Abdoul/These/Article1/Code/MDSV/Code_These/Ess/Forecasts/JvsR"

setwd(path)

files.all<-Sys.glob("Forecast_*")

index<-"sp500"
LEVIER<-TRUE
N<-10

# files.all<-files.all[grepl(paste0("N_",N),files.all, fixed = TRUE)]
files.all<-files.all[grepl(index,files.all, fixed = TRUE)]
files.all<-files.all[grepl(LEVIER,files.all, fixed = TRUE)]
files.all<-files.all[grepl("J",files.all, fixed = TRUE)]

start.date <- as.Date("2000-01-01") 
pdf( file ="VaRSP500.pdf", width =6.5 , height =6)
index_set<-c("S&P 500","NASDAQ 100", "FTSE 100", "NIKKEI 225")
par(mfrow=c(3,1),mar=c(2.5, 2.5, 2, 2))
par(mfrow=c(3,1))
#for(index in c("sp500","nasdaq","ftse","n225")){
for(filename in files.all){
  # filename<-files.all[grepl(index,files.all, fixed = TRUE)]
  out <- read.csv(filename)
  ind <- 1:(nrow(out))
  
  model <- gsub("Forecast_","",filename)
  model <- as.character(gsub("\\_.*", "", model))
  temp  <- paste0("Forecast_",model)
  
  form  <- gsub(paste0(temp,"_"),"",filename)
  form  <- as.character(gsub("\\_.*", "", form))
  temp  <- paste0(temp,"_",form)
  
  N     <- gsub(paste0(temp,"_N_"),"",filename)
  N     <- as.numeric(gsub("\\_.*", "", N))
  K     <- gsub(paste0(temp,"_N_",N,"_K_"),"",filename)
  K     <- as.numeric(gsub("\\_.*", "", K))
  temp  <- paste0(temp,"_N_",N,"_K_",K)
  
  LEVIER <- gsub(paste0(temp,"_LEVIER_"),"",filename)
  LEVIER <- as.logical(gsub("\\_.*", "", LEVIER))
  temp<-paste0(temp,"_LEVIER_",LEVIER,"_length_756")
  
  index <- gsub(paste0(temp,"_"),"",filename)
  index <- gsub("\\_.*", "", index)
  
  out2<-paste0("Forecast_MDSV_R_N_",N,"_K_",K,"_LEVIER_",LEVIER,"_length_756_",index,"*.csv")
  filename<-Sys.glob(out2)
  out2<-read.csv(filename)
  
  t_ <- as.POSIXlt(c(as.character(start.date),out$date), format ="%Y -%m -%d")
  n.days <- rep (365,length (t_))
  n.days [((1900+ t_$ year) %% 4)==0] <- 366
  t_ <- 1900 + t_$ year + (t_$ yday +1)/n.days
  
  x.start <- t_[1] # date de d??but pour les figures
  i.start <- which.min(abs(t_-x.start ))
  x.end <- t_[length (t_)] # date de fin pour les figures
  i.end <- which.min(abs(t_-x.end))
  t_<- t_[(i.start+1):i.end]
  
  t.at <- seq(from = ceiling(x.start ),to=floor(x.end),by =1)
  xlim <- range(t_)
  cex.axis <- 0.85
  cex.mtext <- 0.85
  index<-index_set[which(index == c("sp500","nasdaq","ftse","n225"))]
  
  r_t<-out$rt
  sum(r_t==out2$rt)
  VaR95_1<-out$VaR95
  VaR95_2<-out2$VaR95
  x<-out$date
  
  plot(t_,r_t,type = 'l',xlab="", ylab="", xaxt ="n",main="",xlim=xlim , xaxs ="i")
  lines(t_,VaR95_1,col="red")
  lines(t_,VaR95_2,col="blue")
  axis (side =1, cex.axis=cex.axis, padj=-0.6,at=t.at)
  mtext (paste0(index, " : N = ", N, " et K = ",K), side =3, line=0.5, font=2, las=0, cex=cex.mtext)
  legend("topleft", legend=c("log-rendements", paste0("VaR95 jointe : ", sum(VaR95_1>r_t), " violations"), 
                             paste0("VaR95 univariée : ", sum(VaR95_2>r_t), " violations")),
         col=c("black","red", "blue"), lty=1, cex=0.8)

}
dev.off ()


#-------------------------------------------------------------------------------
# Stationnary distribution and New Impact Curve
#-------------------------------------------------------------------------------

library(MDSV)
library(ggplot2)
library(ggformula)
theme_set(theme_classic())

pdf( file ="Stationnary_NIC.pdf", width =6.5 , height =4)
# par(mfrow=c(1,1))
par(mfrow=c(1,2),mar=c(2.5, 2.5, 2, 2))

path<-"C:/Users/DellPC/Dropbox/Abdoul/These/Article1/Code/MDSV/Code_These/Ess"
setwd(path)
filename  <- "All_estimations.csv"
index     <- "sp500"
LEVIER    <- TRUE
ModelType <- "Joint"
start.date <- as.Date("2000-01-01") 
end.date   <- as.Date("2019-12-31")
donne<-get(index)[rownames(get(index))>=start.date & rownames(get(index))<=end.date,]
donne[,"r"]<-donne[,"r"] - mean(donne[,"r"])
n<-nrow(donne)

out <- read.csv(filename)
vars<-c("omega","a","b","sigma","v0","xi","varphi","delta1","delta2","shape","l","theta")

i_sample <- which((out$ModelType == ModelType) & (out$LEVIER == LEVIER) & (out$index == index) &
                    (out$N %in% c(10,3,6)) & (out$K %in% c(2,10,3)) & (out$nstate > 700))

plotDat<-NULL
for(i in i_sample){
  K<-as.numeric(out[i,"K"])
  N<-as.numeric(out[i,"N"])
  params<-as.numeric(out[i,vars])
  names(params)<-vars

  X<-paste0(vars," = ",as.character(round(params,2)))

  xlab<-c(paste(X[1:(length(X)/2+1)],collapse = ", "),paste(X[(length(X)/2+1):length(X)],collapse = ", "))

  sig<-sqrt(MDSV:::volatilityVector(para=params, K=K, N=N))

  prob<-MDSV:::probapi(params["omega"], K=K, N=N)

  temp<-aggregate(x=prob,by=list(round(sig,4)),FUN="sum")
  temp<-cbind(temp,paste0("N = ",N,", K = ",K))
  plotDat<-rbind(plotDat,temp)
}
names(plotDat)<-c("x","y","Model")
plotDat<-plotDat[(plotDat$y>0.0001),]

t.at <- seq(from = 0, to=max(ceiling(max(plotDat$x))),by =0.3)
xlim <- range(plotDat$x)
cex.axis <- 0.85
cex.mtext <- 0.85

Model_range <- names(table(plotDat$Model))
Model_range <- Model_range[-1]
y<-plotDat$y[plotDat$Model ==  "N = 10, K = 2" ]
x<-plotDat$x[plotDat$Model ==  "N = 10, K = 2" ]

plot(x=x,y=y,type='l',xlab="", ylab="", xaxt ="n",main="",xlim=xlim , xaxs ="i")
for(j in 1:2){
  Model <- Model_range[j]
  y<-plotDat$y[plotDat$Model ==  Model ]
  x<-plotDat$x[plotDat$Model ==  Model ]
  lines(x=x,y=y,lty=j+1,type='l')
}
axis (side =1, cex.axis=cex.axis, padj=-0.6,at=t.at)
mtext ("Distribution stationnaire de la volatilité", side =3, line=0.5, font=2, las=0, cex=cex.mtext)
legend('topright',legend=c("MDSV : N = 10, K = 2", "MDSV : N = 3, K = 10", "MDSV : N = 6, K = 3"),
col=c("black"), lty=1:3, cex=0.8)

proba_lissees_real_EGARCH<-function(ech,para){
  n<-length(ech[,"r"])
  sigma<-numeric(n)
  sigma[1]<-((sqrt(1/(1-2*para[4])))*exp(para[1]+(para[5]^2)*para[8]/2-para[4]+(para[3]^2)/(2*(1-para[4]))))^(1/(1-para[2])) # variance inconditionnel
  mu_rv<-para[6] + para[7]*log(sigma[1])+para[3]*(ech[1,"r"]/sqrt(sigma[1]))+para[4]*((ech[1,"r"]/sqrt(sigma[1]))^2-1)
  aj<-dnorm(ech[1,"r"],0,sqrt(sigma[1]))*dlnorm(ech[1,"rv"],mu_rv,sqrt(para[8]))
  lik<-log(aj)
  
  ajm<-dnorm(ech[1,"r"],0,sqrt(sigma[1]))
  likm<-log(ajm)
  
  for(i in 2:n){
    sigma[i]<-exp(para[1]+para[2]*log(sigma[i-1])+para[5]*(log(ech[i-1,"rv"])-mu_rv)+para[9]*(ech[i-1,"r"]/sqrt(sigma[i-1]))+para[10]*((ech[i-1,"r"]/sqrt(sigma[i-1]))^2-1))
    mu_rv<-para[6] + para[7]*log(sigma[i])+para[3]*(ech[i,"r"]/sqrt(sigma[i]))+para[4]*((ech[i,"r"]/sqrt(sigma[i]))^2-1)
    aj<-dnorm(ech[i,"r"],0,sqrt(sigma[i]))*dlnorm(ech[i,"rv"],mu_rv,sqrt(para[8]))
    lik<-c(lik+log(aj))
    
    ajm<-dnorm(ech[i,"r"],0,sqrt(sigma[i]))
    likm<-c(likm+log(ajm))
  }
  attr(lik,"Marginal")<-likm
  attr(lik,"sig")<-sigma
  return(-lik);
}

r_t <- donne[1:(n-1),"r"]
ind <- order(r_t)

plotDat<-NULL
for(i in i_sample){
  K<-as.numeric(out[i,"K"])
  N<-as.numeric(out[i,"N"])
  
  params<-as.numeric(out[i,vars])
  
  proba_lis<-MDSVfilter(N=N,K=K, data=donne, para = params, ModelType = 2, LEVIER = TRUE)
  
  proba_lis<-proba_lis$smoothed_proba
    
  s<-numeric(n)
  for(i in 1:n) s[i]<-which.max(proba_lis[,i])
  
  Levier <- MDSV:::levierVolatility(ech=donne[,'r'],Nl=70,para=params, Model_type = 2)$`Levier`
  sig    <- MDSV:::volatilityVector(para=params, K=K, N=N)
  
  V_t<-sig[s]*Levier
  y<-V_t[-1]
  # y<-V_t[ind]
  
  X<-data.frame(x=r_t,y=y,Model=paste0("N = ",N,", K = ",K))
  
  plotDat<-rbind(plotDat,X)
}

names(plotDat)<-c("x","y","Model")
t.at <- seq(from = ceiling(min(r_t)), to=max(floor(max(r_t))),by =3)
xlim <- range(r_t)

Model_range <- names(table(plotDat$Model))
Model_range <- Model_range[-4]

Model <- Model_range[1]
y<-plotDat$y[(plotDat$Model == Model)]
lo <- loess(y~r_t)
plot(x=r_t[ind],y=predict(lo)[ind],type='l',xlab="", ylab="", xaxt ="n",main="",xlim=xlim , xaxs ="i",lty=1)

for(j in 2:3){
  Model <- Model_range[j]
  y<-plotDat$y[(plotDat$Model == Model)]
  lo <- loess(y~r_t)
  lines(x=r_t[ind],y=predict(lo)[ind],lty=j)
}
axis (side =1, cex.axis=cex.axis, padj=-0.6,at=t.at)
mtext ("Courbe d'impacts des nouvelles", side =3, line=0.5, font=2, las=0, cex=cex.mtext)
# legend('topleft',legend=c("MDSV : N = 10, K = 2", "MDSV : N = 3, K = 10", "MDSV : N = 6, K = 3", "Real EGARCH"),
#        col=c("black"), lty=1:4, cex=0.8)

filename<-"Estim_ALLRealEGARCH.csv"
out<-read.csv(filename)
vars<-c("w","beta","tau1",'tau2', "gamma",	"xi",	"phi",	"sigma",	"tau1l"	,"tau2l")

i<- which((out$index==index) & (out$LEVIER==LEVIER))
params<-as.numeric(out[i,vars])
names(params)<-vars

sigma<-attr(proba_lissees_real_EGARCH(donne,params),"sig")

plotDat<-rbind(plotDat,data.frame(x=r_t,y=sigma[-1],Model="RealEGARCH"))
ind2<-plotDat$Model=="RealEGARCH"

y<-sigma[-1]
lo <- loess(y~r_t)
lines(x=r_t[ind],y=predict(lo)[ind],lty=4)

dev.off ()


# Other graphs of Leverage
#-------------------------

path<-"C:/Users/DellPC/Dropbox/Abdoul/These/Article1/Code/MDSV/Code_These/Ess"
setwd(path)
filename  <- "All_estimations.csv"
index     <- "sp500"
LEVIER    <- TRUE
ModelType <- "Joint"
start.date <- as.Date("2000-01-01") 
end.date   <- as.Date("2019-12-31")
donne<-get(index)[rownames(get(index))>=start.date & rownames(get(index))<=end.date,]
n<-nrow(donne)
nam_<-rownames(donne)[1:n]
donne[,"r"]<-donne[,"r"] - mean(donne[,"r"])
out <- read.csv(filename)
vars<-c("omega","a","b","sigma","v0","xi","varphi","delta1","delta2","shape","l","theta")

i <- which((out$ModelType == ModelType) & (out$LEVIER == LEVIER) & (out$index == index) & 
                    (out$N %in% c(3)) & (out$K %in% c(10)))

K<-as.numeric(out[i,"K"])
N<-as.numeric(out[i,"N"])

params<-as.numeric(out[i,vars])
names(params)<-vars

Levier1<-MDSV:::levierVolatility(ech=donne[,"r"],Nl=70,para=params,Model_type = 2)$`Levier`
Levier2<-params["l"]*(params["theta"])^(0:69)

par(mfrow=c(1,1))
pdf( file ="LeverageComponent.pdf", width =6.5 , height =6)

t_ <- as.POSIXlt(c(as.character(start.date),nam_), format ="%Y -%m -%d")
n.days <- rep (365,length (t_))
n.days [((1900+ t_$ year) %% 4)==0] <- 366
t_ <- 1900 + t_$ year + (t_$ yday +1)/n.days

x.start <- t_[1] # date de d??but pour les figures
i.start <- which.min(abs(t_-x.start ))
x.end <- t_[length (t_)] # date de fin pour les figures
i.end <- which.min(abs(t_-x.end))
t_<- t_[(i.start+1):i.end]

t.at <- seq(from = ceiling(x.start ),to=floor(x.end),by =1)
xlim <- range(t_)

plot(t_,Levier1,type = 'l',xlab="", ylab="", xaxt ="n",main="",xlim=xlim , xaxs ="i")
axis (side =1, cex.axis=cex.axis, padj=-0.6,at=t.at)
# mtext ("Effet levier journalier", side =3, line=0.5, font=2, las=0, cex=cex.mtext)
dev.off ()


pdf( file ="LeverageParameters.pdf", width =6.5 , height =6)

t.at <- seq(from = 0,to=69,by = 20)
xlim <- range(0:69)

plot(0:69,Levier2,type = 'l',xlab="", ylab="", xaxt ="n",main="",xlim=xlim , xaxs ="i")
axis (side =1, cex.axis=cex.axis, padj=-0.6,at=t.at)
# mtext ("Param??tres d'impact", side =3, line=0.5, font=2, las=0, cex=cex.mtext)

dev.off ()


