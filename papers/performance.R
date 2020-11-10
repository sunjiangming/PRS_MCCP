#!/usr/bin/env Rscript

library("pROC")
library("rms")

library(viridis)
library(ggplot2)


args = commandArgs(trailingOnly=TRUE)
wdir=args[1]
ifile=args[2]
score=args[3]
covar=eval(parse(text=args[4]))
output=args[5]

disease="pseudoPhenotypes"
gname="MCCP"
n_min=9

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

mccp_auc=val
mccp_auc_lower=val
mccp_auc_upper=val

mccp_ppv=val
mccp_ppv_lower=val
mccp_ppv_upper=val

mccp_npv=val
mccp_npv_lower=val
mccp_npv_upper=val


x_end=0.2
x_start=max(min(df$p0),min(df$p1))*1.05
ranges=seq(x_start,x_end,(x_end-x_start)/(100-1))


for( i in 1:len) {
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

  formula1<-as.formula(paste("Y_test_label", paste(c(score,colnames(covar)), collapse=" + "), sep=" ~ "))
  mod_full <- lrm(formula1,x=TRUE, y=TRUE,data=df[singleton,])
   
  y_fitted <- predict(mod_full, type="fitted")

  roc1<-roc(df$Y_test_label[singleton], y_fitted, ci=TRUE, plot = FALSE)
  
  mccp_auc[i]=roc1$auc[1]
  mccp_auc_lower[i]=roc1$ci[1]
  mccp_auc_upper[i]=roc1$ci[3]
  
  ci_s=ci.coords(roc1, "best", ret=c("ppv", "npv"),conf.level=0.95, progress="none",best.policy="random")
  
  mccp_ppv[i]=ci_s$ppv[2]
  mccp_npv[i]=ci_s$npv[2]
  
  mccp_ppv_lower[i]=ci_s$ppv[1]
  mccp_ppv_upper[i]=ci_s$ppv[3]
  
  mccp_npv_lower[i]=ci_s$npv[1]
  mccp_npv_upper[i]=ci_s$npv[3]
  
 }
}


ppv=rep(NA,len)
npv=rep(NA,len)

auc=ppv
auc_lower=ppv
auc_upper=ppv

ppv_lower=ppv
ppv_upper=ppv

npv_lower=ppv
npv_upper=ppv

for(i in 1:len){
 prs=df[,score]
 
 #top and bottom x%
 alpha1=quantile(prs,1-eff[i]/2)
 alpha2=quantile(prs,eff[i]/2)
 
 if( length(unique( df$Y_test_label[prs>=alpha1 | prs<alpha2]))>1 ) {

  formula1<- as.formula(paste("Y_test_label", paste(c(score,colnames(covar)), collapse=" + "), sep=" ~ "))
  mod_full <- lrm(formula1,x=TRUE, y=TRUE,data=df[prs>=alpha1 | prs<alpha2,])
  
  y_fitted <- predict(mod_full, type="fitted")
  roc2<-roc(df$Y_test_label[prs>=alpha1 | prs<alpha2], y_fitted, ci=TRUE, plot = FALSE)
  
  auc[i]=roc2$auc[1]
  auc_lower[i]=roc2$ci[1]
  auc_upper[i]=roc2$ci[3]
  
  ci_s2=ci.coords(roc2, "best", ret=c("ppv", "npv"),conf.level=0.95, progress="none",best.policy="random")
  
  ppv[i]=ci_s2$ppv[2]
  npv[i]=ci_s2$npv[2]
 
  ppv_lower[i]=ci_s2$ppv[1]
  ppv_upper[i]=ci_s2$ppv[3]
  
  npv_lower[i]=ci_s2$npv[1]
  npv_upper[i]=ci_s2$npv[3]
  
  }
}


PPVs=rbind(data.frame(alpha=eff,value=mccp_ppv,value_lower=mccp_ppv_lower, value_upper=mccp_ppv_upper, g=rep(gname,len),col=ranges,disease=rep(disease,len)),data.frame(alpha=eff,value=ppv,value_lower=ppv_lower, value_upper=ppv_upper,g=rep("Empirical",len),col=rep(NA,len),disease=rep(disease,len) ))
NPVs=rbind(data.frame(alpha=eff,value=mccp_npv,value_lower=mccp_npv_lower, value_upper=mccp_npv_upper, g=rep(gname,len),col=ranges,disease=rep(disease,len)),data.frame(alpha=eff,value=npv,value_lower=npv_lower, value_upper=npv_upper,g=rep("Empirical",len),col=rep(NA,len),disease=rep(disease,len) ))
AUCs=rbind(data.frame(alpha=eff,value=mccp_auc,value_lower=mccp_auc_lower, value_upper=mccp_auc_upper, g=rep(gname,len),col=ranges,disease=rep(disease,len)),data.frame(alpha=eff,value=auc,value_lower=auc_lower, value_upper=auc_upper,g=rep("Empirical",len),col=rep(NA,len),disease=rep(disease,len) ))

save(PPVs, NPVs, AUCs, file = paste(output,".RData",sep=""))


#plot

p1<-ggplot(PPVs,aes(alpha,value,group=g)) +
  geom_line(aes(linetype=g,color=g)) +
  geom_ribbon(aes(ymin=value_lower,ymax=value_upper,fill=g),alpha=0.1)+
  geom_point(aes(color=g,size=col)) +
  labs(x="Sample coverage",y="PPV") +
  scale_fill_viridis(discrete=TRUE,end = 0.5,direction=-1, name=" ") +
  scale_color_viridis(discrete=TRUE,end = 0.5,direction=-1,name=" ")+
  theme_minimal()+
  scale_linetype(name=NULL)+
  guides(linetype=FALSE)+
  scale_size_continuous(name="Expected error",range=c(0,1))+
  theme(legend.position="bottom",legend.box="vertical", 
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 10),
        legend.key.width=unit(1,"line"),
        legend.key.height=unit(1,"line"),
        legend.spacing.y = unit(-0.1, "cm"),
        strip.text.x = element_text(size = 8)
  )

p2<-ggplot(NPVs,aes(alpha,value,group=g)) +
  geom_line(aes(linetype=g,color=g)) +
  geom_ribbon(aes(ymin=value_lower,ymax=value_upper,fill=g),alpha=0.1)+
  geom_point(aes(color=g,size=col)) +
  labs(x="Sample coverage",y="NPV") +
  scale_color_discrete(name=" ")+
  scale_fill_viridis(discrete=TRUE,end = 0.5,direction=-1, name=" ") +
  scale_color_viridis(discrete=TRUE,end = 0.5,direction=-1,name=" ")+
  theme_minimal()+
  scale_linetype(name=" ")+
  guides(linetype=FALSE)+
  scale_size_continuous(name="Expected error",range=c(0,1))+
  theme(legend.position="bottom",legend.box="vertical", 
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 10),
        legend.key.width=unit(1,"line"),
        legend.key.height=unit(1,"line"),
        legend.spacing.y = unit(-0.1, "cm"),
        strip.text.x = element_text(size = 8)
  )

p3<-ggplot(AUCs,aes(alpha,value,group=g)) +
  geom_line(aes(linetype=g,color=g)) +
  geom_ribbon(aes(ymin=value_lower,ymax=value_upper,fill=g),alpha=0.1)+
  geom_point(aes(color=g,size=col)) +
  labs(x="Sample coverage",y="AUC") +
  scale_color_discrete(name=" ")+
  scale_fill_viridis(discrete=TRUE,end = 0.5,direction=-1, name=" ") +
  scale_color_viridis(discrete=TRUE,end = 0.5,direction=-1, name=" ")+
  theme_minimal()+
  scale_linetype(name=" ")+
  guides(linetype=FALSE)+
  scale_size_continuous(name="Expected error",range=c(0,1))+
  theme(legend.position="bottom",legend.box="vertical", 
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 10),
        legend.key.width=unit(1,"line"),
        legend.key.height=unit(1,"line"),
        legend.spacing.y = unit(-0.1, "cm")
  )
  
ggsave(paste0(output,"_PPVs.pdf"),p1)
ggsave(paste0(output,"_NPVs.pdf"),p2)
ggsave(paste0(output,"_AUCs.pdf"),p3)
