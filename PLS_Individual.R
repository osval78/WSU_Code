PLS_Individual=function(X_trn,y_trn,X_tst){
  
  ######Fit the model with the optimal number of principal components
  Model_Ind<-plsr(y_trn~X_trn, scale=F, validation="CV")
  
  n_principal_components<- selectNcomp(Model_Ind, method = "onesigma")
  
  ######Fit the model with the optimal number of principal components
  Model_Ind_2<-plsr(y_trn~X_trn, scale=F, validation="none")
  X_tst2=scale(X_tst, center=colMeans(X_tst)-colMeans(X_trn), scale=F)
  X_training=scale(X_trn, center=colMeans(X_trn), scale=F)
  if(n_principal_components<2){
    n_principal_components=1
  }else{
    n_principal_components=n_principal_components
  }
    W_loading=Model_Ind_2$loading.weights[,1:n_principal_components]
  Y_predicted <-c(predict(Model_Ind_2,X_tst2,n_principal_components))
  Y_predicted_trn <-c(predict(Model_Ind_2,X_trn,2,n_principal_components))
  
  #####Matrix containing the regression coefficients
  P_loadings=Model_Ind_2$loadings[,1:n_principal_components]
  Pt_loadings=t(P_loadings)
  ####DtW inverse
  PtW_loadings_inv=solve(Pt_loadings%*%W_loading)
  ###R is equal to W multiplied by DtW_inv
  R=W_loading%*%PtW_loadings_inv
  b=Model_Ind_2$coefficients[,,n_principal_components]
  T_training=X_training%*%R
  X_tst3=scale(X_tst, center=colMeans(X_tst), scale=F)
  T_testing=X_tst3%*%R  
  #####Computing B= beta coefficients with 5 components
  Y_scale=scale(y_trn, center=T, scale=F)
  #Y_scale=as.matrix(train[,1])
  T_manual_2=T_training
  TtT_inv=solve(t(T_manual_2)%*%T_manual_2)
  tTY=t(T_manual_2)%*%Y_scale
  C=TtT_inv%*%tTY
  return(list(LV_trn=T_training,LV_tst=T_testing,R=R,b=b,beta=C,W=W_loading,P=P_loadings,Pred =Y_predicted,Pred_trn=Y_predicted_trn))
}

