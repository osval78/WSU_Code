#install pls package (if not already installed)
#install.packages("pls")
#load pls package
rm(list = ls())
library(pls)
library(SKM)
library(BGLR)

#make this example reproducible
set.seed(1)
source("PLS_Individual.R")

data <- load("Data_Wheat_Carter.RData", verbose = TRUE)
ls()
names(data)
BLUEs=Data_Wheat_Carter$BLUEs_Final
head(BLUEs)
Geno=Data_Wheat_Carter$Geno
dim(Geno)
unique(BLUEs$Year)
unique(BLUEs$Loc)
boxplot( GY~Year,
        data=BLUEs,
        main="Different boxplots for each year",
        xlab="Year",
        ylab="Grain Yield",
        col="orange",
        border="brown"
)
Data_names=c("2019","2020","2021")
Data_names
Summary_All_EC_Line=data.frame()
Summary_All_EC_Image=data.frame()
Summary_All_EC_Joint=data.frame()
Summary_All_EC_Line_trn=data.frame()
Summary_All_EC_Image_trn=data.frame()
Summary_All_EC_Joint_trn=data.frame()

All_summary=data.frame()
All_summary_trn=data.frame()

for(w in 1:3){
 # w=1
  name <-Data_names[w]
  
  Pheno=BLUEs[BLUEs$Year==name,]
  Pheno=droplevels(Pheno)
  rownames(Pheno)=1:nrow(Pheno)
  Gids=sort(unique(Pheno$Hybrid_Name))
  length(Gids)
  Geno1=Geno[Gids,Gids]
  dim(Geno1)
  
  Pheno = Pheno[order(Pheno$Env, Pheno$Hybrid_Name),] #ordenar pheno por environment
  rownames(Pheno) =1:nrow(Pheno)
  
  # Data preparation
  Line <- model.matrix(~0 + as.factor(Hybrid_Name), data = Pheno)
  
  X_Env <- model.matrix(~0 + Env, data = Pheno)
  
  #####Environmenal covariates

  EC=Pheno[,7:17]
  head(EC)
  
  X_Image=data.matrix(EC)
  rownames(X_Image)=1:nrow(X_Image)
  str(X_Image)
  
  X_Env=X_Env
  K.E=X_Env%*%t(X_Env)/ncol(X_Env)
  Geno1=Geno1
  diag(Geno1)=diag(Geno1)+0.001
  L_G <- t(chol(Geno1))
  X_Line<- Line%*%L_G
  typeof(X_Line)
  Geno11=as.matrix(Geno1)
  K.L=Line%*%Geno11%*%t(Line)
  K.LE=K.L*K.E
  X_LE=t(chol(K.LE))
  X1=cbind(X_Env,X_Line,X_Image,X_LE)
  X_EC_Line <- model.matrix(~0+ X_Line:X_Image, data = Pheno)
  X2=cbind(X1,X_EC_Line)
  #head(Pheno)
  Trait_names=colnames(Pheno)[4:6]
  #Trait_names

  All_summary_Lines=data.frame()
  All_summary_Image=data.frame()
  All_summary_Joint=data.frame()
  All_summary_Lines_trn=data.frame()
  All_summary_Image_trn=data.frame()
  All_summary_Joint_trn=data.frame()
  # 7-fold Partition
  set.seed(2022)
  folds <- cv_kfold(records_number = nrow(Pheno),
                    k = 7)
for (t in 1:3){
#  t=1
  Trait_t=Trait_names[t]
  y<- Pheno[, 3+t]
  pos_Na=which(is.na(y)==TRUE)
  pos_Na
  if (length(pos_Na)>0) {
    y[ pos_Na]=median(y)
  } else {
    y=y
  }
  
  Envs=unique(Pheno$Env)
  Predictions_Joint<- data.frame()
  Predictions_Image<- data.frame()
  Predictions_Line<- data.frame()
  Predictions_Image_trn<- data.frame()
  Predictions_Line_trn<- data.frame()
  Predictions_Joint_trn<- data.frame()
 for (e in 1:length(Envs)){
    cat("*** Fold:", e, " ***\n")
#    e=1
   fold <- folds[[e]]
   
   
   pos_test_e=fold$testing
    ####Response variable
    y_trn <- y[-pos_test_e]
    y_tst<-y[pos_test_e]
    y2=y
    
    ######input of lines
    X_trn_Line <-X1[-pos_test_e, ]
    X_tst_Line <-X1[ pos_test_e, ]
   
    PLS_res_Line=PLS_Individual(X_trn=X_trn_Line,y_trn=y_trn,X_tst=X_tst_Line)
    LV_trn_Line=PLS_res_Line$LV_trn
    LV_tst_Line=PLS_res_Line$LV_tst
    Pred_Line=PLS_res_Line$Pred
    Pred_Line_trn=PLS_res_Line$Pred_trn
  #plot(y_tst,Pred_Line)
  cbind(y_tst,Pred_Line)
  cor(y_tst,Pred_Line)
  ###Saving predictions of lines
  Predictions_Line_e=data.frame(Line=Pheno$Hybrid_Name[pos_test_e], Fold=e, Env=Pheno$Env[pos_test_e], Observed=y_tst,Predicted=Pred_Line)
  Predictions_Line=rbind(Predictions_Line,Predictions_Line_e)
  

  ######input of Images
  X_trn_Image <-X2[-pos_test_e, ]
  X_tst_Image <-X2[ pos_test_e, ]
  
  PLS_res_Image=PLS_Individual(X_trn=X_trn_Image,y_trn=y_trn,X_tst=X_tst_Image)
  LV_trn_Image=PLS_res_Image$LV_trn
  LV_tst_Image=PLS_res_Image$LV_tst
  Pred_Image=PLS_res_Image$Pred
  Pred_Image_trn=PLS_res_Image$Pred_trn
  #plot(y_tst,Pred_Image)
  cbind(y_tst,Pred_Image)
  cor(y_tst,Pred_Image)
  
  ###Saving predictions of Images
  Predictions_Image_e=data.frame(Line=Pheno$Hybrid_Name[pos_test_e], Fold=e, Env=Pheno$Env[pos_test_e], Observed=y_tst,Predicted=Pred_Image)
  Predictions_Image=rbind(Predictions_Image,Predictions_Image_e)
  
  }
  
  Pred_Env_Summary_Line=gs_summaries(Predictions_Line)$env
  All_summary_Lines=rbind(All_summary_Lines,data.frame(cbind(Modality="Line",EC=w,Trait_t=Trait_t,Pred_Env_Summary_Line)))


  Pred_Env_Summary_Image=gs_summaries(Predictions_Image)$env
  All_summary_Image=rbind(All_summary_Image,data.frame(cbind(Modality="Image",EC=w, Trait_t=Trait_t,Pred_Env_Summary_Image)))
  

}
  
  Summary_All_EC_Line=rbind(Summary_All_EC_Line,All_summary_Lines)
  Summary_All_EC_Image=rbind(Summary_All_EC_Image,All_summary_Image)

  }
All_summary=rbind(Summary_All_EC_Line,Summary_All_EC_Image)
write.csv(All_summary,file="Carter_PLS_Model2_Predictor_E+G+EC_7FCV_tst.csv")


 