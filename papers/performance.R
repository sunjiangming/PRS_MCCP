#!/usr/bin/env Rscript

library(rms)
library(caret)
library(MLmetrics)

args = commandArgs(trailingOnly=TRUE)
wdir=ags[1]
gname="MCCP"
ifile=args[2]
n_min=49
score="score"
disease=args[3]
ofile=paste(disease,"_",sep="")

returnConMatrix <- function(p0,p1,y,alpha) {
 tp=0
 fp=0
 tn=0
 fn=0
 uncertains=0
 emptys=0
 tp=sum(p1>alpha & p0<=alpha & y==1)
 fp=sum(p1>alpha & p0<=alpha & y==0)
 tn=sum(p0>alpha & p1<=alpha & y==0)
 fn=sum(p0>alpha & p1<=alpha & y==1)
 uncertains=sum(p1>alpha & p0>alpha)
 emptys=sum(p1<=alpha & p0<=alpha)

 err0=sum(p0<=alpha & y==0)
 err1=sum(p1<=alpha & y==1)

 #validity
 val0=err0/max(sum(y==0),1.0)
 val1=err1/max(sum(y==1),1.0)
 val=(err0+err1)/length(y)
 #efficiency
 eff=(tp+fp+tn+fn)/length(y) 
 output=c(tp,fp,tn,fn,uncertains,emptys,1.0-val,1.0-val0,1.0-val1,eff)
 return(output)
}

setwd(wdir)
df=read.table(ifile,head=T,sep="\t")

len=100
val=rep(NA,len)
val0=val
val1=val
eff=val
auc=val
mccp_ppv=val
mccp_npv=val

cm=returnConMatrix(df$p0,df$p1,df$Y_test_label,0.05)
tp=cm[1]
fp=cm[2]
tn=cm[3]
fn=cm[4]
pp=tp+fp
pn=tn+fn

x_end=0.2
x_start=max(min(df$p0),min(df$p1))*1.05
ranges=seq(x_start,x_end,(x_end-x_start)/(len-1))

for( i in 1:length(ranges)) {
 err=ranges[i]
 cm=returnConMatrix(df$p0,df$p1,df$Y_test_label,err)
 tp=cm[1]
 fp=cm[2]
 tn=cm[3]
 fn=cm[4]
 uncertain=cm[5]
 empty=cm[6]
 val[i]=cm[7]
 val0[i]=cm[8]
 val1[i]=cm[9]
 eff[i]=cm[10]
 
 singleton= (df$p0>err & df$p1<=err) | (df$p1>err & df$p0<=err)

 if( eff[i]>0 & length(df$Y_test_label[singleton])>n_min) {

  mod_full <- lrm(Y_test_label ~ score + Batch + Age + Genetic_Sex + PC1 + PC2 + PC3 + PC4 + PC5 + PC6,x=TRUE, y=TRUE,data=df[singleton,])
  mod_covar <- lrm(Y_test_label ~ Batch + Age Genetic_Sex + PC1 + PC2 + PC3 + PC4+ PC5 + PC6,data=df[singleton,])

  auc[i]= mod_full$stats[6]
 
  mccp_ppv[i]=tp/max(tp+fp,1)
  mccp_npv[i]=tn/max(tn+fn,1)
 }
}


ppv=rep(NA,len)
npv=rep(NA,len)
auc_top=ppv

for(i in 1:len){
 prs=df[,score]
 alpha1=quantile(prs,1-eff[i]/2) # top #% PRS
 alpha2=quantile(prs,eff[i]/2) # bottom #%PRS
 
 sdf=df[prs>=alpha1 | prs<alpha2,]
 
 if( !is.na(mccp_ppv[i]) & length(unique( df$Y_test_label[prs>=alpha1 | prs<alpha2]))>1 ) {
  mod_full <- lrm(Y_test_label ~ score + Batch + Age + Genetic_Sex + PC1 + PC2 + PC3 + PC4 + PC5 + PC6,x=TRUE, y=TRUE,data=sdf )
  auc_top[i]= mod_full$stats[6]
  }
 
 idx=is.na(auc)
 auc_top[idx]=NA
 
 ##PPV and NPV
 
 sdf$y=sdf$Y_test_label
 sdf$y[sdf$y==0]="neg"
 sdf$y[sdf$y==1]="pos"
 sdf$y=factor(sdf$y)
 
 set.seed(2020)
 folds <- createFolds(sdf$y, k=5, list = FALSE)
 true_label=NA
 pred_label=NA
 
 model_weights <- ifelse(sdf$y == "neg",
                         (1/table(sdf$y)[1]) * 0.5,
                         (1/table(sdf$y)[2]) * 0.5)
 # train the model on training set
 for( j in 1:5) {
   model <- train(y ~ score + Batch + Age + PC1 + PC2 + PC3 + PC4 + PC5 + PC6,
                  data = sdf[folds!=j,],
                  trControl = trainControl(method = "none",summaryFunction = twoClassSummary,classProbs=TRUE),
                  weights = model_weights[folds!=j],
                  metric="ROC",
                  method = "glm",
                  family=binomial()
   )
   pred_label=c(pred_label,predict(model, newdata = sdf[folds==j,]))
   true_label=c(true_label,sdf$y[folds==j])
 }
 
 true_label=true_label[-1]
 pred_label=pred_label[-1]
 
 lcm=confusionMatrix(factor(pred_label),factor(true_label))
 
 ppv[i]=lcm$byClass[4] #0: control,1: case
 npv[i]=lcm$byClass[3] #0: control,1: case
}


d1=rbind(data.frame(alpha=eff,ppv=mccp_ppv,g=rep(gname,len),col=ranges,disease=rep(disease,len)),data.frame(alpha=eff,ppv=ppv,g=rep("Empirical",len),col=rep(NA,len),disease=rep(disease,len) ))
d2=rbind(data.frame(alpha=eff,npv=mccp_npv,g=rep(gname,len),col=ranges,disease=rep(disease,len)),data.frame(alpha=eff,npv=npv,g=rep("Empirical",len ),col=rep(NA,len),disease=rep(disease,len)) )
d3=rbind(data.frame(alpha=eff,auc=auc,g=rep(gname,len),col=ranges,disease=rep(disease,len)),data.frame(alpha=eff,auc=auc_top,g=rep("Empirical",len),col=rep(NA,len),disease=rep(disease,len)) )

save(d1, d2,d3, file = paste(disease,".RData",sep=""))
