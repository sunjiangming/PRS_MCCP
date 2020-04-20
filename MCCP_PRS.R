#!/usr/bin/env Rscript

pkgs <- list("glmnet", "doParallel", "foreach", "data.table","caret","impute")
lapply(pkgs, require, character.only = T)

args = commandArgs(trailingOnly=TRUE)
wdir=args[1]
tr_file=args[2]
te_file=args[3]
y_col=eval(parse(text=args[4]))
prs_col=eval(parse(text=args[5]))
covar_col=eval(parse(text=args[6]))
ofile=args[7]

computePValues <- function(parameters,X_properTrain,y_properTrain,X_calibration,y_calibration,X_test,y_test){

if(parameters$multiPRS!=0){
	registerDoParallel(cores = parameters$cores)
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
	set.seed(parameters$seed)
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

parameters <- list(
	seed = 2018
	multiPRS   = 0 
	cores = 1
)

setwd(wdir)
set.seed(parameters$seed)

df<-fread(tr_file,head=T)
df=as.data.frame(df)
		    
df2<-fread(te_file,head=T)
df2=as.data.frame(df2)		 

c_idx=c(prs_col,covar_col)
if(length(prs_col)==1) {
	NA_PRS=is.na(df[,prs_col])
} else {
	NA_PRS=is.na(df[,prs_col[1]])
	parameters$multiPRS=1
}

label=df[!NA_PRS,y_col]
pid=df[!NA_PRS,1]
df=df[!NA_PRS,c_idx]

## impute missing values by 10-NN
if(sum(is.na(df))>0){
	imputed2 <- impute.knn(as.matrix(df),k=10,rng.seed=parameters$seed)
	df=imputed2$data
}

## predicting
X_train = data.matrix(df)
y_train = label
		    
p=predictCCP(parameters, X_train, y_train, X_test, y_test)

if( length(c_idx)>1 ) {
	output=data.frame(pid=test_pid,Y_test_label=y_test, p0=p$p0, p1=p$p1, score=as.numeric(X_test[,1]))
}
if( length(c_idx)==1 ) {
	output=data.frame(pid=test_pid,Y_test_label=y_test, p0=p$p0, p1=p$p1, score=as.numeric(X_test))
}
		    
write.table(output,ofile,append=T, quote=F,sep="\t",row.names=F,col.names=F)
		    ##K fold

idx = createFolds(label, k = 5, list = F, returnTrain = F)
for(i in 1:5) {
	#split into training and test set
	test_idx = idx == i
	X_train = data.matrix(df)[!test_idx,]
	y_train = label[!test_idx]
		
	X_test = data.matrix(df)[test_idx,]
	y_test = label[test_idx]
	test_pid=pid[test_idx]
	
	p=predictCCP(parameters, X_train, y_train, X_test, y_test)
	
	if( length(c_idx)>1 ) {
	  output=data.frame(pid=test_pid,Y_test_label=y_test, p0=p$p0, p1=p$p1, score=as.numeric(X_test[,1]))
    	}
    	if( length(c_idx)==1 ) {
	  output=data.frame(pid=test_pid,Y_test_label=y_test, p0=p$p0, p1=p$p1, score=as.numeric(X_test))
 	}

	if(i==1) { write.table(output,ofile,quote=F,sep="\t",row.names=F)
	} else write.table(output,ofile,append=T, quote=F,sep="\t",row.names=F,col.names=F)
}
