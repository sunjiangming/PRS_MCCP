#!/usr/bin/env Rscript

pkgs <- list("glmnet", "doParallel", "foreach", "data.table","caret","impute","Zelig")
lapply(pkgs, require, character.only = T)

args = commandArgs(trailingOnly=TRUE)
wdir=args[1]
ifile=args[2]
y_col=eval(parse(text=args[3]))
prs_col=eval(parse(text=args[4]))
covar_col=eval(parse(text=args[5]))
ofile=args[6]

computePValues <- function(parameters,X_properTrain,y_properTrain,X_calibration,y_calibration,X_test,y_test){

if(parameters$multiPRS!=0){
	registerDoParallel(cores = 20)
	a <- seq(0, 1, 0.05)
	search <- foreach(i = a, .combine = rbind) %dopar% {
		cv <- cv.glmnet(X_properTrain, y_properTrain, family = "binomial", nfold = 5, type.measure = "auc", paralle = TRUE, alpha = i)
		data.frame(cvm = cv$cvm[cv$lambda == cv$lambda.1se], lambda.1se = cv$lambda.1se, alpha = i)
	}
	cv2 <- search[search$cvm == max(search$cvm), ]
	md2 <- glmnet(X_properTrain, y_properTrain, family = "binomial", lambda = cv2$lambda.1se, alpha = cv2$alpha)
	y_calibration_score = predict.lognet(md2,X_calibration,type="link")

        y_test_score = predict.lognet(md2,X_test,type="link")
        y_test_prob = predict.lognet(md2,X_test,type="response")
	
} else {
	df_glm = data.frame(y_properTrain,X_properTrain)
	fit <- glm(y_properTrain ~., family=binomial(link='logit'),data=df_glm)
	if(dim(df_glm)[2]==2){
		y_calibration_score = predict(fit,data.frame(X_properTrain=X_calibration),type="link")
		y_test_score = predict(fit,data.frame(X_properTrain=X_test),type="link")
		y_test_prob = predict(fit,data.frame(X_properTrain=X_test),type="response")
	} else {
                y_calibration_score = predict(fit,data.frame(X_calibration),type="link")
                y_test_score = predict(fit,data.frame(X_test),type="link")
                y_test_prob = predict(fit,data.frame(X_test),type="response")
	}
}

score_zero = y_calibration_score[y_calibration==0]
score_one = -1.0*y_calibration_score[y_calibration==1]

p0 = sapply(y_test_score,function(x) 1.0-1.0*sum(score_zero <= x)/(length(score_zero)+1) )
p1 = sapply(y_test_score,function(x) 1.0-1.0*sum(score_one <= -1.0*x)/(length(score_one)+1) )
p = list("p0"=p0,"p1"=p1)
return(p)
}

predictCCP <- function(parameters,X_train,y_train,X_test,y_test){
	set.seed(2018)
	idx_inner = createFolds(y_train, k = 5, list = F, returnTrain = F)
	for(i in 1:5) {
		X_properTrain = data.matrix(X_train)[idx_inner!=i,]
		y_properTrain = y_train[idx_inner!=i]
		X_calibration = data.matrix(X_train)[idx_inner==i,]
		y_calibration = y_train[idx_inner==i]
	
		pt <- computePValues(parameters,X_properTrain,y_properTrain,X_calibration,y_calibration,X_test,y_test)
		
		if(i==1) {
			p0 <- pt$p0
			p1 <- pt$p1
		} else {
			p0=cbind(p0,pt$p0)
			p1=cbind(p1,pt$p1)
		}
	}
	p = list("p0"=apply(p0,1,mean),"p1"=apply(p1,1,mean))
	return(p)
}

predictLR <- function(X_train,y_train,X_test,y_test){
  set.seed(2018)
  
  df_glm = data.frame(y_train,X_train)
  
  df_glm_weights=rep(1,length(y_train))
  # assume 1 is minority, 0 is majority
  df_glm_weights[y_train==1] = round( ( (2 * sum(y_train == 0)) /length(y_train) )*length(y_train) / (2 * sum(y_train == 1)) )
  
  fit <- zelig(factor(y_train) ~ ., model = "logit", weights = df_glm_weights,
                  data = df_glm, cite = FALSE)
  #
  #fit <- glm(factor(y_train) ~. , weights=df_glm_weights, family = binomial(link='logit'), data=df_glm)
  if(dim(df_glm)[2]==2){
    y_test_prob = predict(fit,data.frame(X_train=X_test),type="response")
  } else {
    y_test_prob = predict(fit,data.frame(X_test),type="response")
  }
  
  names(y_test_prob)="logit_prob"
  p = y_test_prob
  return(p)
}

parameters <- list(
	# General Parameters
	seed  = 2018,
	multiPRS   = 0 
)

setwd(wdir)

df<-fread(ifile,head=T)
df=as.data.frame(df)

c_idx=c(prs_col,covar_col)

if(length(prs_col)==1) {
	NA_PRS=is.na(df[,prs_col])
} else {
	NA_PRS=is.na(df[,prs_col[1]])
	parameters$multiPRS=1
}

label=df[!NA_PRS,y_col]
IID=df[!NA_PRS,1]
df=df[!NA_PRS,c_idx]

if(sum(is.na(df))>0){
	imputed2 <- impute.knn(as.matrix(df),k=10,rng.seed=2018)
	df=imputed2$data
}

#split into 5 folds
set.seed(2018)
idx = createFolds(label, k = 5, list = F, returnTrain = F)

system(paste("if [ -e ",ofile," ]; then rm ",ofile,"; fi",sep=""))

for(i in 1:5) {
	# split train and test
	test_idx = idx == i
	X_train = data.matrix(df)[!test_idx,]
	y_train = label[!test_idx]
		
	X_test = data.matrix(df)[test_idx,]
	y_test = label[test_idx]
	test_IID=IID[test_idx]
	
	logit_p = predictLR(X_train, y_train, X_test, y_test)
	
	p=predictCCP(parameters, X_train, y_train, X_test, y_test)
	
	if( length(c_idx)>1 ) {
	  output=data.frame(IID=test_IID,Y_test_label=y_test, p0=p$p0, p1=p$p1, logit_prob=logit_p, score=as.numeric(X_test[,1]))
    	}
    	if( length(c_idx)==1 ) {
	  output=data.frame(IID=test_IID,Y_test_label=y_test, p0=p$p0, p1=p$p1, logit_prob=logit_p, score=as.numeric(X_test))
 	}

	if(i==1) { write.table(output,ofile,quote=F,sep="\t",row.names=F)
	} else write.table(output,ofile,append=T, quote=F,sep="\t",row.names=F,col.names=F)
}
