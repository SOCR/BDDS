---
title: "Predictive Big Data Analytics (PBDA) for Parkinson's Disease"
author: "SOCR Team"
date: "September 18, 2016"
output: html_document
---
  
  
  This is an R Markdown document. Markdown is a simple formatting syntax for authoring HTML, PDF, and MS Word documents. For more details on using R Markdown see <http://rmarkdown.rstudio.com>.

###install required packages #######
```{r}
install.packages("geepack",repos="http://cran.us.r-project.org")
install.packages(c("lme4","mi","MASS"),repos="http://cran.us.r-project.org")
install.packages("unbalanced",repos="http://cran.us.r-project.org") 
install.packages("crossval",repos="http://cran.us.r-project.org") 
install.packages("ada",repos="http://cran.us.r-project.org")
install.packages("e1071",repos="http://cran.us.r-project.org")
install.packages("rpart",repos="http://cran.us.r-project.org")
install.packages("class",repos="http://cran.us.r-project.org")
install.packages(c("RWeka","rJava"),repos="http://cran.us.r-project.org")
install.packages("lme4",repos="http://cran.us.r-project.org") 
install.packages("geepack",repos="http://cran.us.r-project.org") 
install.packages("mi",repos="http://cran.us.r-project.org") 
install.packages("unbalanced",repos="http://cran.us.r-project.org") 
options(warn=-1) 

```


```{r}
####data manipulation #######
setwd("~/documents/PPMI") #set working director
clinical<-read.csv("ppmi_clinical.csv",header=TRUE) #read raw data
colnames(clinical)[2]<-"FID_IID" 
ppmi_new<-read.csv("PPMI_GSA_clinical.csv",header=TRUE)
merged_dataset<-merge(ppmi_new,clinical,by="FID_IID",all.x=TRUE) #merge clinical data with brain image
write.csv(merged_dataset,"ppmi_clinical_complete.csv")
ppmi<-merged_dataset[,c(1:399,820:823,833,832,816)] #Select "tier 1" UPDRSs
head(ppmi)
#longitudinal analysis
ppmi_longit<-read.csv("ppmi_top_wide_UPDRS.csv",header=TRUE)
clinical<-read.csv("ppmi_clinical.csv",header=TRUE)
merged_dataset<-merge(ppmi_longit,clinical,by="FID_IID",all.x=TRUE)
ppmi_long<-merged_dataset[,c(1:342,351:354,364,365,347)]
head(ppmi_long)
```






####Longtitudinal analysis with different time points of UPDRS measurements ####
```{r}
# use complete dataset with the result
DataMissOm<-na.omit(subset(ppmi_long,select=c(time_visit,PD,FID_IID,COGSTATE,COGDECLN,FNCDTCOG,COGDXCL,EDUCYRS,L_superior_parietal_gyrus_ComputeArea
                                              ,L_superior_parietal_gyrus_Volume,R_superior_parietal_gyrus_ComputeArea,R_superior_parietal_gyrus_Volume,L_putamen_ComputeArea,L_putamen_Volume,R_putamen_Volume 
                                              ,R_putamen_ShapeIndex,L_caudate_ComputeArea,L_caudate_Volume,R_caudate_ComputeArea,R_caudate_Volume,chr12_rs34637584_GT 
                                              ,chr17_rs11868035_GT,chr17_rs11012_GT,chr17_rs393152_GT,chr17_rs12185268_GT,chr17_rs199533_GT,Sex,Weight,Age,UPDRS_part_I,UPDRS_part_II,UPDRS_part_III)))
#normalize continuous variables
for(i in c(8:20,28:32)){
  DataMissOm[,i]<-scale(DataMissOm[,i])
}

library(geepack)
DataMissOm$COGSTATE<-as.numeric(DataMissOm$COGSTATE)
DataMissOm$COGDECLN<-as.numeric(DataMissOm$COGDECLN)
DataMissOm$FNCDTCOG<-as.numeric(DataMissOm$FNCDTCOG)
DataMissOm$COGDXCL<-as.numeric(DataMissOm$COGDXCL)
model_gee_UPDRS<-geeglm(PD~L_superior_parietal_gyrus_ComputeArea
                        +L_superior_parietal_gyrus_Volume+R_superior_parietal_gyrus_ComputeArea+R_superior_parietal_gyrus_Volume+L_putamen_ComputeArea+L_putamen_Volume+R_putamen_Volume 
                        +R_putamen_ShapeIndex+L_caudate_ComputeArea+L_caudate_Volume+R_caudate_ComputeArea+R_caudate_Volume+chr12_rs34637584_GT 
                        +chr17_rs11868035_GT+chr17_rs11012_GT+chr17_rs393152_GT+chr17_rs12185268_GT+chr17_rs199533_GT+Sex+Weight+Age
                        +UPDRS_part_I+UPDRS_part_II+UPDRS_part_III+FID_IID+COGSTATE+COGDECLN+FNCDTCOG+COGDXCL+EDUCYRS,
                        data=DataMissOm,waves=time_visit,family="binomial",id=FID_IID,corstr="exchangeable",scale.fix=TRUE)
summary(model_gee_UPDRS)

library(lme4)
model_lmm_UPDRS<-glmer(PD~L_superior_parietal_gyrus_ComputeArea
                       +L_superior_parietal_gyrus_Volume+R_superior_parietal_gyrus_ComputeArea+R_superior_parietal_gyrus_Volume+L_putamen_ComputeArea+L_putamen_Volume+R_putamen_Volume 
                       +R_putamen_ShapeIndex+L_caudate_ComputeArea+L_caudate_Volume+R_caudate_ComputeArea+R_caudate_Volume+chr12_rs34637584_GT 
                       +chr17_rs11868035_GT+chr17_rs11012_GT+chr17_rs393152_GT+chr17_rs12185268_GT+chr17_rs199533_GT+Sex+Weight+Age+UPDRS_part_I+UPDRS_part_II+UPDRS_part_III+FID_IID+COGSTATE+COGDECLN+FNCDTCOG+COGDXCL+EDUCYRS
                       +(1|FID_IID),control = glmerControl(optCtrl=list(maxfun=20000)),
                       data=DataMissOm,family="binomial")
summary(model_lmm_UPDRS)

length(unique(DataMissOm$FID_IID)) #406 patients
nrow(DataMissOm) #2562 observations

```

####Longtitudinal analysis with different time point of brain imaging measures ####
```{r}
ppmi_base<-read.csv("ppmi_ROI_baseline.csv",header=TRUE)
ppmi_base<-ppmi_base[,c(2,24:26,29:33)] #may add 27 28
ppmi_new<-read.csv("PPMI_GSA_clinical.csv",header=TRUE)
ppmi_imaging<-ppmi_new[,1:335]

merged_dataset<-merge(ppmi_imaging,ppmi_base,by="FID_IID",all.x=TRUE)
DataMissOm<-na.omit(subset(merged_dataset,select=c(VisitID,ResearchGroup,FID_IID,COGSTATE,COGDECLN,FNCDTCOG,COGDXCL,EDUCYRS,L_superior_parietal_gyrus_ComputeArea
                                                   ,L_superior_parietal_gyrus_Volume,R_superior_parietal_gyrus_ComputeArea,R_superior_parietal_gyrus_Volume,L_putamen_ComputeArea,L_putamen_Volume,R_putamen_Volume 
                                                   ,R_putamen_ShapeIndex,L_caudate_ComputeArea,L_caudate_Volume,R_caudate_ComputeArea,R_caudate_Volume,chr12_rs34637584_GT 
                                                   ,chr17_rs11868035_GT,chr17_rs11012_GT,chr17_rs393152_GT,chr17_rs12185268_GT,chr17_rs199533_GT,Sex,Weight,Age,UPDRS_Part_I_Summary_Score_Baseline,UPDRS_Part_II_Patient_Questionnaire_Summary_Score_Baseline,UPDRS_Part_III_Summary_Score_Baseline)))
#normalize continuous variables
for(i in c(9:20,28:32)){
  DataMissOm[,i]<-scale(DataMissOm[,i])
}
library(geepack)
DataMissOm$COGSTATE<-as.numeric(DataMissOm$COGSTATE)
DataMissOm$COGDECLN<-as.numeric(DataMissOm$COGDECLN)
DataMissOm$FNCDTCOG<-as.numeric(DataMissOm$FNCDTCOG)
DataMissOm$COGDXCL<-as.numeric(DataMissOm$COGDXCL)
DataMissOm$PD<-ifelse(DataMissOm$ResearchGroup=="Control",0,1)
model_gee_imaging<-geeglm(PD~L_superior_parietal_gyrus_ComputeArea
                          +L_superior_parietal_gyrus_Volume+R_superior_parietal_gyrus_ComputeArea+R_superior_parietal_gyrus_Volume+L_putamen_ComputeArea+L_putamen_Volume+R_putamen_Volume 
                          +R_putamen_ShapeIndex+L_caudate_ComputeArea+L_caudate_Volume+R_caudate_ComputeArea+R_caudate_Volume+chr12_rs34637584_GT 
                          +chr17_rs11868035_GT+chr17_rs11012_GT+chr17_rs393152_GT+chr17_rs12185268_GT+chr17_rs199533_GT+Sex+Weight+Age
                          +UPDRS_Part_I_Summary_Score_Baseline+UPDRS_Part_II_Patient_Questionnaire_Summary_Score_Baseline+UPDRS_Part_III_Summary_Score_Baseline+FID_IID+COGSTATE+COGDECLN+FNCDTCOG+COGDXCL+EDUCYRS,
                          data=DataMissOm,waves=VisitID,family="binomial",id=FID_IID,corstr="exchangeable",scale.fix=TRUE)
summary(model_gee_imaging)

library(lme4)
#model_lmm_imaging<-glmer(PD~L_superior_parietal_gyrus_ComputeArea
 #                        +L_superior_parietal_gyrus_Volume+R_superior_parietal_gyrus_ComputeArea+R_superior_parietal_gyrus_Volume+L_putamen_ComputeArea+L_putamen_Volume+R_putamen_Volume 
  #                       +R_putamen_ShapeIndex+L_caudate_ComputeArea+L_caudate_Volume+R_caudate_ComputeArea+R_caudate_Volume+chr12_rs34637584_GT 
   #                      +chr17_rs11868035_GT+chr17_rs11012_GT+chr17_rs393152_GT+chr17_rs12185268_GT+chr17_rs199533_GT+Sex+Weight+Age+UPDRS_Part_I_Summary_Score_Baseline+UPDRS_Part_II_Patient_Questionnaire_Summary_Score_Baseline+UPDRS_Part_III_Summary_Score_Baseline+FID_IID+COGSTATE+COGDECLN+FNCDTCOG+COGDXCL+EDUCYRS
    #                     +(1|FID_IID),#control = glmerControl(optCtrl=list(maxfun=20000)),
     #                    data=DataMissOm,family="binomial")
#summary(model_lmm_imaging)
length(unique(DataMissOm$FID_IID))
nrow(DataMissOm) 


```

#### baseline analysis ####
```{r}
ppmi_1<-ppmi[ppmi$VisitID==1,]
ppmi_1_complete<-ppmi_1[,c(1,
                           3:282, # imaging index
                           283:284,287, #demographic index
                           288,291,294,297,300,303, #genotype
                           336,348,360,384,392, #UPDRS+asessment non motor may delete 372
                           400:404, #clinical index
                           405:406)]
ppmi_1_complete$ind<-is.na(ppmi_1_complete$chr17_rs199533_GT)
ppmi_1_complete<-ppmi_1_complete[ppmi_1_complete$ind==FALSE,-303]
ppmi_MI<-ppmi_1_complete[,c(291:300)]

# Multiple Imputation
library(mi)
set.seed(1414)
mdf<-missing_data.frame(ppmi_MI)
summary(mdf)
imputation<-mi(mdf)
complete<-complete(imputation) 
c4<-complete[[4]] #extract complete dataset from it from imputed datasets. 
# in this case, we randomly chose 4th. 
c4<-c4[,1:10]
ppmi_1<-cbind(ppmi_1_complete[,1:290],c4,ppmi_1_complete[,301:302]) #complete dataset 

ppmi_1_ROI<-ppmi_1[,c(1,23,24 ,28,29,33,34,38,39,143,144,148,149,282:302)] #complete dataset of regions of interest

write.csv(ppmi_1,"ppmi_complete_baseline.csv")
write.csv(ppmi_1_ROI,"ppmi_ROI_baseline.csv")

#####explore the distribution of baseline data ####
#categorical
col_n<-c(14,17:22,28:31)
#ppmi_1_ROI$APPRDX<-as.numeric(ppmi_1_ROI$APPRDX)

# frequency table
for (i in col_n){
  print(colnames(ppmi_1_ROI[i]))
  m<-table(ppmi_1_ROI[,i],ppmi_1_ROI[,33])
  print(m)
  print(fisher.test(m)$p.value)
}
for (i in col_n){
  print(colnames(ppmi_1_ROI[i]))
  m<-table(ppmi_1_ROI[,i],ppmi_1_ROI[,34])
  print(m)
  print(fisher.test(m)$p.value)
}

#continuous
#normality test
col_n<-c(2:13,15,16,23:27,32)
aa<-matrix(NA,20,1)
k<-1
for(i in col_n){
  aa[k,1]<-shapiro.test(ppmi_1_ROI[,i])$p.value
  k<-k+1
}
k<-1
#for(i in col_n){
# x<-ppmi_1_ROI[ppmi_1_ROI$RECRUITMENT_CAT=="HC",i] 
#  y<-ppmi_1_ROI[ppmi_1_ROI$RECRUITMENT_CAT=="PD",i] 
#  aa[k,2]<-t.test(x,y)$p.value
#  k<-k+1
#}

aa<-as.data.frame(aa)
rownames(aa)<-colnames(ppmi_1_ROI)[col_n]
colnames(aa)<-c("normality-test (>0.1 normal)")
write.csv(aa,"continuous_descriptive.csv")
```



###model-based analysis####
```{r}
#logistic regression
for(i in c(2:13,15,16,23:27,32)){
  ppmi_1_ROI[,i]<-scale(ppmi_1_ROI[,i])
}
ppmi_1_ROI<-as.data.frame(ppmi_1_ROI)
levels(ppmi_1_ROI$RECRUITMENT_CAT)<-c(0,1)
ppmi_1_ROI$RECRUITMENT_CAT<-as.numeric(ppmi_1_ROI$RECRUITMENT_CAT)
ppmi_1_ROI$RECRUITMENT_CAT<-ppmi_1_ROI$RECRUITMENT_CAT-1
logistic_ROI<-glm(RECRUITMENT_CAT~.,data=ppmi_1_ROI[,-c(1,34)],maxit=50,family="binomial")
library(MASS)
stepp<-stepAIC(logistic_ROI,trace=FALSE) #stepwise logistic regression
summary(logistic_ROI)
summary(stepp)


#logistic_complete<-glm(RECRUITMENT_CAT~.,data=ppmi_1[,-c(1,302)],maxit=50,family="binomial")
#stepp<-stepAIC(logistic_complete,trace=FALSE) #stepwise model selection of ROI data
#summary(stepp)


#ppmi_1 complete
ppmi_1_PD1<-ppmi_1[ppmi_1$APPRDX!="psychiatric",]
ppmi_1_PD1<-ppmi_1[ppmi_1$APPRDX!="SWEDD",]
ppmi_1_PD1<-ppmi_1_PD1[,-c(1,302)]
ppmi_1_PD2<-ppmi_1[ppmi_1$APPRDX!="psychiatric",]
ppmi_1_PD2<-ppmi_1[,-c(1,302)]

ppmi_1_PD1$PD1<-ifelse(ppmi_1_PD1$RECRUITMENT_CAT=="HC",1,0)
ppmi_1_PD1<-ppmi_1_PD1[,-300]
ppmi_1_PD1$PD1<-as.factor(ppmi_1_PD1$PD1)
```


####Balance cases####
```{r}
#PD1=PD PD2=PD+SWEDD
set.seed(623)

library(unbalanced)
n<-ncol(ppmi_1_PD1)
output<-ppmi_1_PD1$PD1
input<-ppmi_1_PD1[ ,-n]

#balance the dataset
data<-ubBalance(X= input, Y=output,type="ubSMOTE", verbose=TRUE)
balancedData<-cbind(data$X,data$Y)
ppmi_1_PD1_balanced<-balancedData
colnames(ppmi_1_PD1_balanced)<-colnames(ppmi_1_PD1)


ppmi_1_PD2$PD2<-ifelse(ppmi_1_PD2$RECRUITMENT_CAT=="HC",1,0)
ppmi_1_PD2<-ppmi_1_PD2[,-300]
ppmi_1_PD2$PD2<-as.factor(ppmi_1_PD2$PD2)
#balance cases
set.seed(623)
require(unbalanced)
n<-ncol(ppmi_1_PD2)
output<-ppmi_1_PD2$PD2
input<-ppmi_1_PD2[ ,-n]

#balance the dataset
data<-ubBalance(X= input, Y=output,type="ubSMOTE", verbose=TRUE)
balancedData<-cbind(data$X,data$Y)
ppmi_1_PD2_balanced<-balancedData
colnames(ppmi_1_PD2_balanced)<-colnames(ppmi_1_PD2)

#numerize the 4 datasets
for(i in 295:298){
  ppmi_1_PD1[,i]<-as.numeric(ppmi_1_PD1[,i])
}
for(i in 295:298){
  ppmi_1_PD2[,i]<-as.numeric(ppmi_1_PD2[,i])
}
for(i in 295:298){
  ppmi_1_PD1_balanced[,i]<-as.numeric(ppmi_1_PD1_balanced[,i])
}
for(i in 295:298){
  ppmi_1_PD2_balanced[,i]<-as.numeric(ppmi_1_PD2_balanced[,i])
}

write.csv(ppmi_1_PD1,"PD1_un.csv")
write.csv(ppmi_1_PD2,"PD2_un.csv")
write.csv(ppmi_1_PD1_balanced,"PD1_b.csv")
write.csv(ppmi_1_PD2_balanced,"PD2_b.csv")
```


```{r}
### Prediction functions ####
# cv 5 functions
library(crossval)
#adaboost
library(ada)
library(rpart)
my.ada <- function (train.x, train.y, test.x, test.y, negative, formula){
  ada.fit <- ada(train.x, train.y)
  predict.y <- predict(ada.fit, test.x)
  #count TP,NP..
  out <- confusionMatrix(test.y, predict.y, negative = negative)
  return (out)
}


# svm
library(e1071)
my.svm <- function (train.x, train.y, test.x, test.y, negative, formula){
  svm.fit <- svm(train.x, train.y)
  predict.y <- predict(svm.fit, test.x)
  #count TP,NP..
  out <- confusionMatrix(test.y, predict.y, negative = negative)
  return (out)
}
# naive bayesian
my.nb <- function (train.x, train.y, test.x, test.y, negative, formula){
  nb.fit <- naiveBayes(train.x, train.y)
  predict.y <- predict(nb.fit, test.x)
  #count TP,NP..
  out <- confusionMatrix(test.y, predict.y, negative = negative)
  return (out)
}
# decision trees

library(rpart)
my.detree <- function (train.x, train.y, test.x, test.y, negative, formula){
  train.matrix<-cbind(train.x,train.y)
  detree.fit <- rpart(train.y~.,method="class",data=train.matrix)
  predict.y <- predict(detree.fit, test.x)
  out <- confusionMatrix(test.y, predict.y, negative = negative)
  return (out)
}
##k-Nearest Neighbor ClassiÃ¯Â¬cation

library("class")
my.knn <- function (train.x, train.y, test.x, test.y, negative, formula){
  train.matrix<-cbind(train.x,train.y)
  test.matrix<-cbind(test.x,test.y)
  predict.y <- knn(train=train.matrix,test=test.matrix,cl=train.y,k=2)
  out <- confusionMatrix(test.y, predict.y, negative = negative)
  return (out)
}
## Kmeans
my.kmeans <- function (train.x, train.y, test.x, test.y, negative, formula){
  train.matrix<-cbind(train.x,train.y)
  test.matrix<-cbind(test.x,test.y)
  kmeans.fit <- kmeans(train.x,2)
  predict.y <- kmeans.fit$cluster-1
  out <- confusionMatrix(test.y,predict.y, negative = negative)
  return (out)
}
```



### function dealt with complete dataset ####
```{r}
classfi<-function(a=my.ada){
  a1<-matrix(NA,4,4)
  a2<-matrix(NA,4,6)
  #PD1 unbalanced
  X <- ppmi_1_PD1[, 1:299]
  Y <- ppmi_1_PD1[,300] 
  form <- "train.y~." 
  neg <- 0
  set.seed(714)
  cv.out <- crossval(a, X, Y, K = 5, B = 1, negative = neg, formula = form)
  out <- diagnosticErrors(cv.out$stat) 
  a1[1,]<-cv.out$stat
  a2[1,]<-out
  
  
  
  #PD2 unbalanced
  X <- ppmi_1_PD2[, 1:299]
  Y <- ppmi_1_PD2[,300] 
  form <- "train.y~." 
  neg <- 0
  set.seed(714)
  cv.out <- crossval(a, X, Y, K = 5, B = 1, negative = neg, formula = form)
  out <- diagnosticErrors(cv.out$stat) 
  
  a1[2,]<-cv.out$stat
  a2[2,]<-out
  
  #PD1 balanced
  X <- ppmi_1_PD1_balanced[, 1:299]
  Y <- ppmi_1_PD1_balanced[,300] 
  form <- "train.y~." 
  neg <- 0
  set.seed(714)
  cv.out <- crossval(a, X, Y, K = 5, B = 1, negative = neg, formula = form)
  out <- diagnosticErrors(cv.out$stat) 
  
  a1[3,]<-cv.out$stat
  a2[3,]<-out
  
  #PD2 balanced
  X <- ppmi_1_PD2_balanced[, 1:299]
  Y <- ppmi_1_PD2_balanced[,300] 
  form <- "train.y~." 
  neg <- 0
  set.seed(714)
  cv.out <- crossval(a, X, Y, K = 5, B = 1, negative = neg, formula = form)
  out <- diagnosticErrors(cv.out$stat) 
  a1[4,]<-cv.out$stat
  a2[4,]<-out
  classification<-as.data.frame(cbind(a1,a2))
  classification$Group<-c("PD1_un","PD2_un","PD1_b","PD2_b")
  write.csv(classification,"classification.csv")
  print(classification)
}
```


### Classification result
```{r}
# adaboost
classfi()
# detree
classfi(a=my.detree)
#naive bayes
classfi(a=my.nb)
#svm
classfi(a=my.svm)
# knn
classfi(a=my.knn)
# kmeans
classfi(a=my.kmeans)
```


####function dealt with dataset without UPDRS (and clinical features)
```{r}
classfi1<-function(a=my.ada){
  a1<-matrix(NA,4,4)
  a2<-matrix(NA,4,6)
  #PD1 unbalanced
  X <- ppmi_1_PD1[, 1:289]
  Y <- ppmi_1_PD1[,300] 
  form <- "train.y~." 
  neg <- 0
  set.seed(714)
  cv.out <- crossval(a, X, Y, K = 5, B = 1, negative = neg, formula = form)
  out <- diagnosticErrors(cv.out$stat) 
  a1[1,]<-cv.out$stat
  a2[1,]<-out
  
  
  
  #PD2 unbalanced
  X <- ppmi_1_PD2[, 1:289]
  Y <- ppmi_1_PD2[,300] 
  form <- "train.y~." 
  neg <- 0
  set.seed(714)
  cv.out <- crossval(a, X, Y, K = 5, B = 1, negative = neg, formula = form)
  out <- diagnosticErrors(cv.out$stat) 
  
  a1[2,]<-cv.out$stat
  a2[2,]<-out
  
  #PD1 balanced
  X <- ppmi_1_PD1_balanced[, 1:289]
  Y <- ppmi_1_PD1_balanced[,300] 
  form <- "train.y~." 
  neg <- 0
  set.seed(714)
  cv.out <- crossval(a, X, Y, K = 5, B = 1, negative = neg, formula = form)
  out <- diagnosticErrors(cv.out$stat) 
  
  a1[3,]<-cv.out$stat
  a2[3,]<-out
  
  #PD2 balanced
  X <- ppmi_1_PD2_balanced[, 1:289]
  Y <- ppmi_1_PD2_balanced[,300] 
  form <- "train.y~." 
  neg <- 0
  set.seed(714)
  cv.out <- crossval(a, X, Y, K = 5, B = 1, negative = neg, formula = form)
  out <- diagnosticErrors(cv.out$stat) 
  a1[4,]<-cv.out$stat
  a2[4,]<-out
  classification<-cbind(a1,a2)
  write.csv(classification,"classification2.csv")
  print(classification)
}
#adaboost
classfi1(a=my.ada)
#detree
classfi1(a=my.detree)
#naive bayes
classfi1(a=my.nb)
# svm
classfi1(a=my.svm)
# knn
classfi1(a=my.knn)
# kmeans
classfi1(a=my.kmeans)
```


###adaboost scores and results ####
```{r}
# Varplot, table S2 and table S3
require(ada)
#unbalanced PD1
set.seed(729) 
c1<-ada(PD1~.,data=ppmi_1_PD1)
```

```{r, echo=FALSE}
varplot(c1)
```

```{r}
n_c<-ncol(ppmi_1_PD1)-1
PD1_score<-varplot(c1,type="scores",max.var.show=n_c)
top30<-varplot(c1,type="scores")
data_top30<-cbind(names(top30),"PD1","unBalanced")
write.csv(data_top30,"top30.csv")

#unbalanced PD2
set.seed(729) 
c1<-ada(PD2~.,data=ppmi_1_PD2)
```

```{r, echo=FALSE}
varplot(c1)
```

```{r}
n_c<-ncol(ppmi_1_PD2)-1
PD2_score<-varplot(c1,type="scores",max.var.show=n_c)
top30<-varplot(c1,type="scores")
data_top30<-cbind(names(top30),"PD2","unBalanced")
write.csv(data_top30,"top30.csv")

#balanced PD1
set.seed(729) 
c1<-ada(PD1~.,data=ppmi_1_PD1_balanced)
```

```{r, echo=FALSE}
varplot(c1)
```

```{r}
n_c<-ncol(ppmi_1_PD1_balanced)-1
PD1_score_b<-varplot(c1,type="scores",max.var.show=n_c)
top30<-varplot(c1,type="scores")
data_top30<-cbind(names(top30),"PD1","Balanced")
write.csv(data_top30,"top30.csv")


#balanced PD2
set.seed(729) 
c1<-ada(PD2~.,data=ppmi_1_PD2_balanced)
```

```{r, echo=FALSE}
varplot(c1)
```

```{r}
n_c<-ncol(ppmi_1_PD2_balanced)-1
PD2_score_b<-varplot(c1,type="scores",max.var.show=n_c)
top30<-varplot(c1,type="scores")
data_top30<-cbind(names(top30),"PD2","Balanced")
write.csv(data_top30,"top30.csv")


write.csv(PD1_score,"trial1.csv")
write.csv(PD2_score,"trial2.csv")
write.csv(PD1_score_b,"trial3.csv")
write.csv(PD2_score_b,"trial4.csv")
########
```



###without UPDRS ####
```{r}
#unbalanced PD1
set.seed(729) 
c1<-ada(PD1~.,data=ppmi_1_PD1[,-c(290:294)])
```

```{r, echo=FALSE}
varplot(c1)
```

```{r}
n_c<-ncol(ppmi_1_PD1[,-c(290:294)])-1
PD1_score_noUPDRS<-varplot(c1,type="scores",max.var.show=n_c)
top30<-varplot(c1,type="scores")
data_top30<-cbind(names(top30),"PD1","unBalanced")
write.csv(data_top30,"top30.csv")

#unbalanced PD2
set.seed(729) 
c1<-ada(PD2~.,data=ppmi_1_PD2[,-c(290:294)])
```

```{r, echo=FALSE}
varplot(c1)
```

```{r}
n_c<-ncol(ppmi_1_PD2[,-c(290:294)])-1
PD2_score_noUPDRS<-varplot(c1,type="scores",max.var.show=n_c)
top30<-varplot(c1,type="scores")
data_top30<-cbind(names(top30),"PD2","unBalanced")
write.csv(data_top30,"top30.csv")


#balanced PD1
set.seed(729) 
c1<-ada(PD1~.,data=ppmi_1_PD1_balanced[,-c(290:294)])
```

```{r, echo=FALSE}
varplot(c1)
```

```{r}
n_c<-ncol(ppmi_1_PD1_balanced[,-c(290:294)])-1
PD1_score_b_noUPDRS<-varplot(c1,type="scores",max.var.show=n_c)
top30<-varplot(c1,type="scores")
data_top30<-cbind(names(top30),"PD1","Balanced")
write.csv(data_top30,"top30.csv")

#balanced PD2
set.seed(729) 
c1<-ada(PD2~.,data=ppmi_1_PD2_balanced[,-c(290:294)])
```

```{r, echo=FALSE}
varplot(c1)
```

```{r}
n_c<-ncol(ppmi_1_PD2_balanced[,-c(290:294)])-1
PD2_score_b_noUPDRS<-varplot(c1,type="scores",max.var.show=n_c)
top30<-varplot(c1,type="scores")
data_top30<-cbind(names(top30),"PD2","Balanced")
write.csv(data_top30,"top30.csv")

write.csv(PD1_score_noUPDRS,"trial5.csv")
write.csv(PD2_score_noUPDRS,"trial6.csv")
write.csv(PD1_score_b_noUPDRS,"trial7.csv")
write.csv(PD2_score_b_noUPDRS,"trial8.csv")
```

Note that the `echo = FALSE` parameter was added to the code chunk to prevent printing of the R code that generated the plot.
