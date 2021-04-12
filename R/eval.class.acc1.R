eval.class.acc1 <-
function(){
  
  
  #Working here: classification accuracy training and test set
  if(num_levels==2){
    if(split.train.test==TRUE){
      
      res<-get_classification.accuracy(kfold=kfold,featuretable=data_m_fc_withfeats[goodip,],classlabels=Y,classifier="logit",testfeaturetable=data_m_fc_withfeats_test[goodip,],testclasslabels=Ytest,errortype="BAR",kernelname=svm_kernel,svm.cost=svm.cost,svm.gamma=svm.gamma)
      
      
    }else{
      #get_roc(dataA=good_feats,classlabels=Y,classifier="svm",kname=svm_kernel,rocfeatlist=seq(2,10,1),rocfeatincrement=TRUE,testset=X,testclasslabels=Y,mainlabel="Using training set based on SVM")
      
      #get_roc(dataA=good_feats,classlabels=Y,classifier="svm",kname=svm_kernel,rocfeatlist=c(dim(good_feats)[2]),rocfeatincrement=FALSE,testset=good_feats_test,testclasslabels=Ytest,mainlabel="Test set")
      
      res<-get_classification.accuracy(kfold=kfold,featuretable=good_feats,classlabels=Y,classifier="logit",testfeaturetable=good_feats_test,testclasslabels=Ytest,errortype="BAR",kernelname=svm_kernel,svm.cost=svm.cost,svm.gamma=svm.gamma)
      
      
    }
  }
  
  test_acc_mat<-{}
  class_levels_vec<-levels(as.factor(Y))
  
  
  learnmatrix <- matrix(seq(1,nrow(X)), nrow = 1)
  fiveCV10iter<-new("learningsets", learnmatrix = learnmatrix, method = "none",ntrain = ncol(learnmatrix), iter = nrow(learnmatrix))
  X<-rbind(X,Xtest)
  Y<-c(Y,Ytest)
  classifier_names<-c("PLSLDA","PLSRF","SCDA","SVM","RF","NNet","pLR","PLSLR","pLRlasso","pLRelasticnet")
  
  result_cma_list<-new("list")
  
  
  ###savelist=ls(),file="debug3.Rda")
  
  confusion_matrix_list<-new("list")
  
  if(length(good_feats_index)>2)
  {
    
    if(tune_classifiers==TRUE){
      s1<-ldply(tune_plslda@tuneres,rbind)
      
      s2<-apply(s1,2,median)
      
      
      t1<-new("list")
      confusion_matrix_list<-new("list")
      t1[[1]]<-s2
      tune_plslda1<-tune_plslda
      tune_plslda1@tuneres<-t1
      
      set.seed(seedvalue)
      class_res <- classification(X = X, y = Y, learningsets = fiveCV10iter, classifier = pls_ldaCMA, tuneres = tune_plslda1,trace=FALSE)
      
      b1<-best(tune_plslda1)
      
      
      
      learnmatrix<-as.numeric(learnmatrix)
      
      set.seed(seedvalue)
      class_res2<-pls_ldaCMA(X = X, y = Y, learnind = learnmatrix, comp = median(unlist(b1)))
      
      
    }else{
      
      set.seed(seedvalue)
      class_res <- classification(X = X, y = Y, learningsets = fiveCV10iter, classifier = pls_ldaCMA, trace=FALSE)
      learnmatrix<-as.numeric(learnmatrix)
      
      set.seed(seedvalue)
      class_res2<-pls_ldaCMA(X = X, y = Y, learnind = learnmatrix)
      
      
    }
    ##saveclass_res2,file="pls_ldaCMA.Rda")
    result_cma_list[[1]]<-class_res2
    
    ###saveconfusion_matrix_1,file="confusion_matrix_plslda.Rda")
    confusion_matrix_list[[1]]<-table(class_res2@y,class_res2@yhat)
    
    if(length(class_levels_vec)==2){
      testeval_res_auc<-evaluation(class_res,measure = "auc")
      ###savetesteval_res_auc,file="testeval_res_auc.Rda")
      
      test_auc<-100*(mean(testeval_res_auc@score))
      test_acc_mat<-c(test_acc_mat,test_auc)
      
    }else{
      
      test_acc<-evaluation(class_res)
      
      test_acc<-100*(1-mean(testeval_res_auc@score))
      test_acc_mat<-c(test_acc_mat,test_acc) #cbind(best_kfold_acc,permkfold_acc,test_acc)
      
      
    }
    
    
    
    #if(best_classifier_name==classifier_names[2]){
    
    if(tune_classifiers==TRUE){
      s1<-ldply(tune_plsrf@tuneres,rbind)
      s2<-apply(s1,2,median)
      t1<-new("list")
      t1[[1]]<-s2
      
      tune_plsrf1<-tune_plsrf
      tune_plsrf1@tuneres<-t1
      set.seed(seedvalue)
      class_res <- classification(X = X, y = Y, learningsets = fiveCV10iter, classifier = pls_rfCMA, tuneres = tune_plsrf1,trace=FALSE)
      
      b1<-best(tune_plsrf1)
      # ###saveb1,file="b1_pls_rfCMA.Rda")
      
      learnmatrix<-as.numeric(learnmatrix)
      
      set.seed(seedvalue)
      class_res2<-pls_rfCMA(X = X, y = Y, learnind = learnmatrix, comp = median(unlist(b1)))
    }else{
      set.seed(seedvalue)
      class_res <- classification(X = X, y = Y, learningsets = fiveCV10iter, classifier = pls_rfCMA, trace=FALSE)
      learnmatrix<-as.numeric(learnmatrix)
      
      set.seed(seedvalue)
      class_res2<-pls_rfCMA(X = X, y = Y, learnind = learnmatrix)
      
      
    }
    
    ##saveclass_res2,file="pls_rfCMA.Rda")
    result_cma_list[[2]]<-class_res2
    #confusion_matrix_1<-table(class_res2@y,class_res2@yhat)
    ###saveconfusion_matrix_1,file="confusion_matrix_rfCMA.Rda")
    
    confusion_matrix_list[[2]]<-table(class_res2@y,class_res2@yhat)
    
    if(length(class_levels_vec)==2){
      
      testeval_res_auc<-evaluation(class_res,measure = "auc")
      ###savetesteval_res_auc,file="testeval_res_auc.Rda")
      
      test_auc<-100*(mean(testeval_res_auc@score))
      
      test_acc_mat<-c(test_acc_mat,test_auc)
      
      
      
    }else{
      
      test_acc<-evaluation(class_res)
      
      test_acc<-100*(1-mean(test_acc@score))
      test_acc_mat<-c(test_acc_mat,test_acc) #cbind(best_kfold_acc,permkfold_acc,test_acc)
      
      
    }
    
    
    if(tune_classifiers==TRUE){
      
      s1<-ldply(tune_scda@tuneres,rbind)
      s2<-apply(s1,2,median)
      
      t1<-new("list")
      t1[[1]]<-s2
      
      tune_scda1<-tune_scda
      tune_scda1@tuneres<-t1
      
      #tune_scda<-new("tuningresult",tuneres=t1)
      
      set.seed(seedvalue)
      class_res <- classification(X = X, y = Y, learningsets = fiveCV10iter, classifier = scdaCMA, tuneres = tune_scda1,trace=FALSE)
      #scdaCMA(X, y, f, learnind, delta = 0.5, models=FALSE,...)
      
      b1<-best(tune_scda1)
      
      
      learnmatrix<-as.numeric(learnmatrix)
      set.seed(seedvalue)
      class_res2<-scdaCMA(X = X, y = Y, learnind = learnmatrix, delta = median(unlist(b1)))
      
    }else{
      
      if(is.na(class_scda)==FALSE){
        set.seed(seedvalue)
        class_res <- classification(X = X, y = Y, learningsets = fiveCV10iter, classifier = scdaCMA,trace=FALSE)
        set.seed(seedvalue)
        class_res2<-scdaCMA(X = X, y = Y, learnind = learnmatrix)
        
      }
    }
    
    if(is.na(class_scda)==FALSE){
      
      ##saveclass_res2,file="scdaCMA.Rda")
      result_cma_list[[3]]<-class_res2
      #confusion_matrix_1<-table(class_res2@y,class_res2@yhat)
      
      ###saveconfusion_matrix_1,file="confusion_matrix_scdaCMA.Rda")
      confusion_matrix_list[[3]]<-table(class_res2@y,class_res2@yhat)
      
      if(length(class_levels_vec)==2){
        
        testeval_res_auc<-evaluation(class_res,measure = "auc")
        ###savetesteval_res_auc,file="testeval_res_auc.Rda")
        
        
        test_auc<-100*(mean(testeval_res_auc@score))
        
        test_acc_mat<-c(test_acc_mat,test_auc)
        
        
        
      }else{
        
        test_acc<-evaluation(class_res)
        
        test_acc<-100*(1-mean(test_acc@score))
        test_acc_mat<-c(test_acc_mat,test_acc) #cbind(best_kfold_acc,permkfold_acc,test_acc)
        
        
      }
    }else{
      result_cma_list[[3]]<-NA
      confusion_matrix_list[[3]]<-NA
      test_acc_mat<-c(test_acc_mat,0)
      
    }
    
    
    
    #if(best_classifier_name==classifier_names[4])
    {
      
      if(tune_classifiers==TRUE){
        s1<-ldply(tune_svm@tuneres,rbind)
        
        s2<-apply(s1,2,median)
        
        t1<-new("list")
        t1[[1]]<-s2
        
        
        tune_svm1<-tune_svm
        tune_svm1@tuneres<-t1
        
        b1<-best(tune_svm1)
        
        
        #class_res <- classification(X = X, y = Y, learningsets = fiveCV10iter, classifier = svmCMA, tuneres=tune_svm1, kernel ="radial",trace=FALSE,probability=TRUE)
        
        
        #           ##save(X,Y,learnmatrix,seedvalue,fiveCV10iter,tune_svm1,b1,svmCMA,file="svmdebug.Rda")
        set.seed(seedvalue)
        class_res <- classification(X = X, y = Y, learningsets = fiveCV10iter, classifier = svmCMA, tuneres=tune_svm1, kernel ="radial",trace=FALSE,probability=TRUE)
        
        
        cost_1=b1[[1]]$cost
        gamma_1=b1[[1]]$gamma
        
        learnmatrix<-as.numeric(learnmatrix)
        
        set.seed(seedvalue)
        #class_res2<-svmCMA(X = X, y = Y, learnind = learnmatrix, cost = cost_1,gamma=gamma_1,probability=TRUE)
        class_res2<-svmCMA(X = X, y = Y, learnind = learnmatrix,probability=TRUE)
      }else{
        
        set.seed(seedvalue)
        class_res <- classification(X = X, y = Y, learningsets = fiveCV10iter, classifier = svmCMA, kernel ="radial",trace=FALSE,probability=TRUE)
        set.seed(seedvalue)
        class_res2<-svmCMA(X = X, y = Y, learnind = learnmatrix,probability=TRUE)
        
      }
      
      ##saveclass_res2,file="svmCMA.Rda")
      result_cma_list[[4]]<-class_res2
      confusion_matrix_1<-table(class_res2@y,class_res2@yhat)
      
      ###saveconfusion_matrix_1,file="confusion_matrix_svmCMA.Rda")
      confusion_matrix_list[[4]]<-table(class_res2@y,class_res2@yhat)
      
      if(length(class_levels_vec)==2){
        
        testeval_res_auc<-evaluation(class_res,measure = "auc")
        ####savetesteval_res_auc,file="testeval_res_auc.Rda")
        
        test_auc<-100*(mean(testeval_res_auc@score))
        
        test_acc_mat<-c(test_acc_mat,test_auc)
        
        
        
      }else{
        
        test_acc<-evaluation(class_res)
        
        test_acc<-100*(1-mean(test_acc@score))
        test_acc_mat<-c(test_acc_mat,test_acc) #cbind(best_kfold_acc,permkfold_acc,test_acc)
        
        
      }
      
    }
    
    #if(best_classifier_name==classifier_names[5])
    {
      
      if(FALSE){
        s1<-ldply(tune_rf@tuneres,rbind)
        s2<-apply(s1,2,median)
        
        t1<-new("list")
        t1[[1]]<-s2
        
        tune_rf1<-tune_rf
        tune_rf1@tuneres<-t1
      }
      
      set.seed(seedvalue)
      class_res <- classification(X =X, y = Y, learningsets = fiveCV10iter, classifier = rfCMA,trace=FALSE) #,tuneres=tune_rf1
      
      learnmatrix<-as.numeric(learnmatrix)
      
      
      set.seed(seedvalue)
      class_res2 <- rfCMA(X =X, y = Y, learnind=learnmatrix, varimp = FALSE) #mtry=tune_rf1$mtry,nodesize=tune_rf1$nodesize
      
      ##saveclass_res2,file="rfCMA.Rda")
      result_cma_list[[5]]<-class_res2
      
      confusion_matrix_1<-table(class_res2@y,class_res2@yhat)
      
      ###saveconfusion_matrix_1,file="confusion_matrix_rfCMA.Rda")
      confusion_matrix_list[[5]]<-table(class_res2@y,class_res2@yhat)
      
      if(length(class_levels_vec)==2){
        
        testeval_res_auc<-evaluation(class_res,measure = "auc")
        ####savetesteval_res_auc,file="testeval_res_auc.Rda")
        
        test_auc<-100*(mean(testeval_res_auc@score))
        
        test_acc_mat<-c(test_acc_mat,test_auc)
        
        
        
      }else{
        
        test_acc<-evaluation(class_res)
        
        test_acc<-100*(1-mean(test_acc@score))
        test_acc_mat<-c(test_acc_mat,test_acc) #cbind(best_kfold_acc,permkfold_acc,test_acc)
        
        
      }
      
    }
    
    #if(best_classifier_name==classifier_names[6])
    {
      
      #if(FALSE){
      
      if(tune_classifiers==TRUE){
        s1<-ldply(tune_svm@tuneres,rbind)
        
        s2<-apply(s1,2,median)
        
        t1<-new("list")
        t1[[1]]<-s2
        
        
        tune_svm1<-tune_svm
        tune_svm1@tuneres<-t1
        
        b1<-best(tune_svm1)
        s1<-ldply(tune_nnet@tuneres,rbind)
        s2<-apply(s1,2,median)
        
        t1<-new("list")
        t1[[1]]<-s2
        
        tune_nnet1<-tune_nnet
        tune_nnet1@tuneres<-t1
        
        
        b1<-best(tune_nnet1)
        #}
        
        set.seed(seedvalue)
        class_res<- classification(X =X, y = Y, learningsets = fiveCV10iter, classifier = nnetCMA,trace=FALSE,tuneres=tune_nnet1)
        #size=3,decay=0.1) #,
        testeval_res_auc<-evaluation(class_res,measure = "auc")
        
        set.seed(seedvalue)
        #nnetresult <- nnetCMA(X=golubX, y=golubY, learnind=learnind, size = 3, decay = 0.01)
        class_res2 <- nnetCMA(X =X, y = Y, learnind=as.numeric(learnmatrix),sigma=median(unlist(best(tune_nnet))))
        #size=b1[[1]]$size,decay=b1[[1]]$decay)
        #size=3,decay=0.1)
      }else{
        
        
        set.seed(seedvalue)
        class_res<- classification(X =X, y = Y, learningsets = fiveCV10iter, classifier = nnetCMA,trace=FALSE)
        
        testeval_res_auc<-evaluation(class_res,measure = "auc")
        
        set.seed(seedvalue)
        #nnetresult <- nnetCMA(X=golubX, y=golubY, learnind=learnind, size = 3, decay = 0.01)
        class_res2 <- nnetCMA(X =X, y = Y, learnind=as.numeric(learnmatrix))
        
        
        
      }
      ##saveclass_res2,file="nnetCMA.Rda")
      result_cma_list[[6]]<-class_res2
      confusion_matrix_1<-table(class_res2@y,class_res2@yhat)
      
      ###saveconfusion_matrix_1,file="confusion_matrix_nnetCMA.Rda")
      confusion_matrix_list[[6]]<-table(class_res2@y,class_res2@yhat)
      
      if(length(class_levels_vec)==2){
        
        testeval_res_auc<-evaluation(class_res,measure = "auc")
        ####savetesteval_res_auc,file="testeval_res_auc.Rda")
        
        test_auc<-100*(mean(testeval_res_auc@score))
        
        test_acc_mat<-c(test_acc_mat,test_auc)
        
        
      }else{
        
        test_acc<-evaluation(class_res)
        
        test_acc<-100*(1-mean(test_acc@score))
        test_acc_mat<-c(test_acc_mat,test_acc) #cbind(best_kfold_acc,permkfold_acc,test_acc)
        
        
      }
      
    }
    
    #if(best_classifier_name==classifier_names[7])
    {
      
      if(FALSE){
        s1<-ldply(tune_plr@tuneres,rbind)
        s2<-apply(s1,2,median)
        
        t1<-new("list")
        t1[[1]]<-s2
        
        tune_plr1<-tune_plr
        tune_plr1@tuneres<-t1
        
        
        b1<-best(tune_plr1)
      }
      
      set.seed(seedvalue)
      class_res<- classification(X =X, y = Y, learningsets = fiveCV10iter, classifier = plrCMA,trace=FALSE) #,tuneres=tune_plr
      
      set.seed(seedvalue)
      class_res2 <- plrCMA(X =X, y = Y, learnind=learnmatrix)
      
      ##saveclass_res2,file="plrCMA.Rda")
      result_cma_list[[7]]<-class_res2
      
      confusion_matrix_1<-table(class_res2@y,class_res2@yhat)
      
      ###saveconfusion_matrix_1,file="confusion_matrix_plrCMA.Rda")
      confusion_matrix_list[[7]]<-table(class_res2@y,class_res2@yhat)
      
      if(length(class_levels_vec)==2){
        
        testeval_res_auc<-evaluation(class_res,measure = "auc")
        ####savetesteval_res_auc,file="testeval_res_auc.Rda")
        
        test_auc<-100*(mean(testeval_res_auc@score))
        
        test_acc_mat<-c(test_acc_mat,test_auc)
        
        
        
      }else{
        
        test_acc<-evaluation(class_res)
        
        test_acc<-100*(1-mean(test_acc@score))
        test_acc_mat<-c(test_acc_mat,test_acc) #cbind(best_kfold_acc,permkfold_acc,test_acc)
        
        
      }
      
    }
    
    
    if(length(class_levels_vec)==2){
      
      
      
      acc_mat<-cbind(emat1,test_acc_mat) #cbind(best_kfold_acc,permkfold_acc,test_acc)
      
      
      
      if(is.na(test_data_check)==TRUE){
        
        colnames(acc_mat)<-c("Training kfold CV Accuracy (AUC)","Training permuted kfold CV accuracy","Score", "Training accuracy(AUC)")
        
        print("Training set evaluation using selected features and the best classifier based on AUC measure")
      }else{
        colnames(acc_mat)<-c("Training kfold CV Accuracy (AUC)","Training permuted kfold CV accuracy","Score", "Test accuracy(AUC)")
        
        print("Test set evaluation using selected features and the best classifier based on AUC measure")
      }
      
      print(acc_mat[,c(1,4)])
      ##save(acc_mat,file="acc_mat.Rda")
      #rownames(acc_mat)<-best_classifier_name[1]
      
      write.table(acc_mat,file="Classification_evaluation_results_AUC.txt",sep="\t")
      
      
    }else{
      
      acc_mat<-cbind(emat1,test_acc_mat) #cbind(best_kfold_acc,permkfold_acc,test_acc)
      
      
      if(is.na(test_data_check)==TRUE){
        
        colnames(acc_mat)<-c("Training kfold CV Accuracy (Misclassification)","Training permuted kfold CV Accuracy (Misclassification)","Score", "Training accuracy (Misclassification)")
        print("Training set evaluation using selected features and the best classifier based on misclassification rate measure")
        
      }else{
        colnames(acc_mat)<-c("Training kfold CV Accuracy (Misclassification)","Training permuted kfold CV Accuracy (Misclassification)","Score", "Test accuracy (Misclassification)")
        print("Test set evaluation using selected features and the best classifier based on misclassification rate measure")
      }
      
      
      #rownames(acc_mat)<-best_classifier_name[1]
      
      write.table(acc_mat,file="Classification_evaluation_results_misclassification.txt",sep="\t")
      
    }
    
    
    
    acc_mat<-acc_mat[,-c(3)]
    
    
    mainlab<-paste("Performance evaluation using classifiers and selected features",sep="")
    
    if(output.device.type!="pdf"){
      
      temp_filename_1<-"Figures/Barplot_classification_accuracy.png"
      
      png(temp_filename_1,width=plots.width,height=plots.height,res=plots.res,type=plots.type,units="in")
    }
    
    
    acc_mat1<-t(acc_mat)
    #xaxt="n",
    #
    w <- 0.1
    par(omd=c(0, 1-w, 0, 1),cex.main=0.7)
    if(length(class_levels_vec)==2){
      barplot(acc_mat1[c(1,3),],main=mainlab,ylab="Classification accuracy(%)",ylim=c(50,100),xpd=FALSE,beside=TRUE,col=barplot.col.opt,cex.axis=0.7,cex.names=0.7)
    }else{
      
      barplot(acc_mat1[c(1,3),],main=mainlab,ylab="Classification accuracy(%)",ylim=c(50,100),xpd=FALSE,beside=TRUE,col=barplot.col.opt,cex.axis=0.7,cex.names=0.7)
      
    }
    if(FALSE){
      if(length(class_levels_vec)==2){
        axis(side=1,at=seq(1,4),labels=c("kfold CV","Permuted kfold CV","Test set","AUC"),cex.axis=cex.plots)
      }else{
        axis(side=1,at=seq(1,3),labels=c("kfold CV","Permuted kfold CV","Test set"),cex.axis=cex.plots)
      }
    }
    
    if(is.na(test_data_check)==TRUE){
      le1<-(legend(par('usr')[2], par('usr')[4], bty='n', xpd=NA,c("Training (kfold)","Training (overall)"), pch = rep(19,2), pt.cex = 0.6, title = "",cex=0.7,col=barplot.col.opt))
      
    }else{
      le1<-(legend(par('usr')[2], par('usr')[4], bty='n', xpd=NA,c("Training (kfold)","Test"), pch = rep(19,2), pt.cex = 0.6, title = "",cex=0.7,col=barplot.col.opt))
    }
    
    if(output.device.type!="pdf"){
      
      try(dev.off(),silent=TRUE)
    }
    
    
    try(dev.off(),silent=TRUE)
    
  }else{
    
    test_acc_mat<-c(0,0)
    
    if(tune_classifiers==TRUE){
      
      s1<-ldply(tune_scda@tuneres,rbind)
      s2<-apply(s1,2,median)
      
      t1<-new("list")
      t1[[1]]<-s2
      
      tune_scda1<-tune_scda
      tune_scda1@tuneres<-t1
      
      #tune_scda<-new("tuningresult",tuneres=t1)
      
      set.seed(seedvalue)
      class_res <- classification(X = X, y = Y, learningsets = fiveCV10iter, classifier = scdaCMA, tuneres = tune_scda1,trace=FALSE)
      #scdaCMA(X, y, f, learnind, delta = 0.5, models=FALSE,...)
      
      b1<-best(tune_scda1)
      
      
      learnmatrix<-as.numeric(learnmatrix)
      set.seed(seedvalue)
      class_res2<-scdaCMA(X = X, y = Y, learnind = learnmatrix, delta = median(unlist(b1)))
      
    }else{
      
      if(is.na(class_scda)==FALSE){
        set.seed(seedvalue)
        class_res <- classification(X = X, y = Y, learningsets = fiveCV10iter, classifier = scdaCMA,trace=FALSE)
        set.seed(seedvalue)
        class_res2<-scdaCMA(X = X, y = Y, learnind = learnmatrix)
        
      }
    }
    
    if(is.na(class_scda)==FALSE){
      
      ##saveclass_res2,file="scdaCMA.Rda")
      result_cma_list[[3]]<-class_res2
      #confusion_matrix_1<-table(class_res2@y,class_res2@yhat)
      
      ###saveconfusion_matrix_1,file="confusion_matrix_scdaCMA.Rda")
      confusion_matrix_list[[3]]<-table(class_res2@y,class_res2@yhat)
      
      if(length(class_levels_vec)==2){
        
        testeval_res_auc<-evaluation(class_res,measure = "auc")
        ###savetesteval_res_auc,file="testeval_res_auc.Rda")
        
        
        test_auc<-100*(mean(testeval_res_auc@score))
        
        test_acc_mat<-c(test_acc_mat,test_auc)
        
        
        
      }else{
        
        test_acc<-evaluation(class_res)
        
        test_acc<-100*(1-mean(test_acc@score))
        test_acc_mat<-c(test_acc_mat,test_acc) #cbind(best_kfold_acc,permkfold_acc,test_acc)
        
        
      }
    }else{
      result_cma_list[[3]]<-NA
      confusion_matrix_list[[3]]<-NA
      test_acc_mat<-c(test_acc_mat,0)
      
    }
    
    
    
    #if(best_classifier_name==classifier_names[4])
    {
      
      if(tune_classifiers==TRUE){
        s1<-ldply(tune_svm@tuneres,rbind)
        
        s2<-apply(s1,2,median)
        
        t1<-new("list")
        t1[[1]]<-s2
        
        
        tune_svm1<-tune_svm
        tune_svm1@tuneres<-t1
        
        b1<-best(tune_svm1)
        
        
        #class_res <- classification(X = X, y = Y, learningsets = fiveCV10iter, classifier = svmCMA, tuneres=tune_svm1, kernel ="radial",trace=FALSE,probability=TRUE)
        
        
        #                ##save(X,Y,learnmatrix,seedvalue,fiveCV10iter,tune_svm1,b1,svmCMA,file="svmdebug.Rda")
        set.seed(seedvalue)
        class_res <- classification(X = X, y = Y, learningsets = fiveCV10iter, classifier = svmCMA, tuneres=tune_svm1, kernel ="radial",trace=FALSE,probability=TRUE)
        
        
        cost_1=b1[[1]]$cost
        gamma_1=b1[[1]]$gamma
        
        learnmatrix<-as.numeric(learnmatrix)
        
        
        set.seed(seedvalue)
        #class_res2<-svmCMA(X = X, y = Y, learnind = learnmatrix, cost = cost_1,gamma=gamma_1,probability=TRUE)
        class_res2<-svmCMA(X = X, y = Y, learnind = learnmatrix,probability=TRUE)
        
      }else{
        
        set.seed(seedvalue)
        class_res <- classification(X = X, y = Y, learningsets = fiveCV10iter, classifier = svmCMA, kernel ="radial",trace=FALSE,probability=TRUE)
        set.seed(seedvalue)
        class_res2<-svmCMA(X = X, y = Y, learnind = learnmatrix,probability=TRUE)
        
      }
      
      ##saveclass_res2,file="svmCMA.Rda")
      result_cma_list[[4]]<-class_res2
      confusion_matrix_1<-table(class_res2@y,class_res2@yhat)
      
      ###saveconfusion_matrix_1,file="confusion_matrix_svmCMA.Rda")
      confusion_matrix_list[[4]]<-table(class_res2@y,class_res2@yhat)
      
      if(length(class_levels_vec)==2){
        
        testeval_res_auc<-evaluation(class_res,measure = "auc")
        ####savetesteval_res_auc,file="testeval_res_auc.Rda")
        
        test_auc<-100*(mean(testeval_res_auc@score))
        
        test_acc_mat<-c(test_acc_mat,test_auc)
        
        
        
      }else{
        
        test_acc<-evaluation(class_res)
        
        test_acc<-100*(1-mean(test_acc@score))
        test_acc_mat<-c(test_acc_mat,test_acc) #cbind(best_kfold_acc,permkfold_acc,test_acc)
        
        
      }
      
    }
    
    #if(best_classifier_name==classifier_names[5])
    {
      
      if(FALSE){
        s1<-ldply(tune_rf@tuneres,rbind)
        s2<-apply(s1,2,median)
        
        t1<-new("list")
        t1[[1]]<-s2
        
        tune_rf1<-tune_rf
        tune_rf1@tuneres<-t1
      }
      
      set.seed(seedvalue)
      class_res <- classification(X =X, y = Y, learningsets = fiveCV10iter, classifier = rfCMA,trace=FALSE) #,tuneres=tune_rf1
      
      learnmatrix<-as.numeric(learnmatrix)
      
      
      set.seed(seedvalue)
      class_res2 <- rfCMA(X =X, y = Y, learnind=learnmatrix, varimp = FALSE) #mtry=tune_rf1$mtry,nodesize=tune_rf1$nodesize
      
      ##saveclass_res2,file="rfCMA.Rda")
      result_cma_list[[5]]<-class_res2
      
      confusion_matrix_1<-table(class_res2@y,class_res2@yhat)
      
      ###saveconfusion_matrix_1,file="confusion_matrix_rfCMA.Rda")
      confusion_matrix_list[[5]]<-table(class_res2@y,class_res2@yhat)
      
      if(length(class_levels_vec)==2){
        
        testeval_res_auc<-evaluation(class_res,measure = "auc")
        ####savetesteval_res_auc,file="testeval_res_auc.Rda")
        
        test_auc<-100*(mean(testeval_res_auc@score))
        
        test_acc_mat<-c(test_acc_mat,test_auc)
        
        
        
      }else{
        
        test_acc<-evaluation(class_res)
        
        test_acc<-100*(1-mean(test_acc@score))
        test_acc_mat<-c(test_acc_mat,test_acc) #cbind(best_kfold_acc,permkfold_acc,test_acc)
        
        
      }
      
    }
    
    #if(best_classifier_name==classifier_names[6])
    {
      
      #if(FALSE){
      if(tune_classifiers==TRUE){
        s1<-ldply(tune_nnet@tuneres,rbind)
        s2<-apply(s1,2,median)
        
        t1<-new("list")
        t1[[1]]<-s2
        
        tune_nnet1<-tune_nnet
        tune_nnet1@tuneres<-t1
        
        
        b1<-best(tune_nnet1)
        
        b2=unlist(unlist(b1))
        
        
        # }
        
        set.seed(seedvalue)
        class_res<- classification(X =X, y = Y, learningsets = fiveCV10iter, classifier = nnetCMA,trace=FALSE,tuneres=tune_nnet) #size=3,decay=0.1) #,tuneres=tune_nnet1
        testeval_res_auc<-evaluation(class_res,measure = "auc")
        
        set.seed(seedvalue)
        #nnetresult <- nnetCMA(X=golubX, y=golubY, learnind=learnind, size = 3, decay = 0.01)
        class_res2 <- nnetCMA(X =X, y = Y, learnind=as.numeric(learnmatrix),size=median(b2[seq(1,length(b2),2)]),decay=median(b2[seq(2,length(b2),2)]))
        #sigma=median(unlist(best(tune_nnet)))) #) #size=3,decay=0.1)
      }else{
        
        set.seed(seedvalue)
        class_res<- classification(X =X, y = Y, learningsets = fiveCV10iter, classifier = nnetCMA,trace=FALSE) #size=3,decay=0.1) #,tuneres=tune_nnet1
        testeval_res_auc<-evaluation(class_res,measure = "auc")
        
        set.seed(seedvalue)
        #nnetresult <- nnetCMA(X=golubX, y=golubY, learnind=learnind, size = 3, decay = 0.01)
        class_res2 <- nnetCMA(X =X, y = Y, learnind=as.numeric(learnmatrix))
        
      }
      ##saveclass_res2,file="nnetCMA.Rda")
      result_cma_list[[6]]<-class_res2
      confusion_matrix_1<-table(class_res2@y,class_res2@yhat)
      
      ###saveconfusion_matrix_1,file="confusion_matrix_nnetCMA.Rda")
      confusion_matrix_list[[6]]<-table(class_res2@y,class_res2@yhat)
      
      if(length(class_levels_vec)==2){
        
        testeval_res_auc<-evaluation(class_res,measure = "auc")
        ####savetesteval_res_auc,file="testeval_res_auc.Rda")
        
        test_auc<-100*(mean(testeval_res_auc@score))
        
        test_acc_mat<-c(test_acc_mat,test_auc)
        
        
      }else{
        
        test_acc<-evaluation(class_res)
        
        test_acc<-100*(1-mean(test_acc@score))
        test_acc_mat<-c(test_acc_mat,test_acc) #cbind(best_kfold_acc,permkfold_acc,test_acc)
        
        
      }
      
      
    }
    
    
    
    
    
    
    
    
    
    
    #if(best_classifier_name==classifier_names[7])
    {
      
      if(FALSE){
        s1<-ldply(tune_plr@tuneres,rbind)
        s2<-apply(s1,2,median)
        
        t1<-new("list")
        t1[[1]]<-s2
        
        tune_plr1<-tune_plr
        tune_plr1@tuneres<-t1
        
        
        b1<-best(tune_plr1)
      }
      
      set.seed(seedvalue)
      class_res<- classification(X =X, y = Y, learningsets = fiveCV10iter, classifier = plrCMA,trace=FALSE) #,tuneres=tune_plr
      
      set.seed(seedvalue)
      class_res2 <- plrCMA(X =X, y = Y, learnind=learnmatrix)
      
      ##saveclass_res2,file="plrCMA.Rda")
      result_cma_list[[7]]<-class_res2
      
      confusion_matrix_1<-table(class_res2@y,class_res2@yhat)
      
      ###saveconfusion_matrix_1,file="confusion_matrix_plrCMA.Rda")
      confusion_matrix_list[[7]]<-table(class_res2@y,class_res2@yhat)
      
      if(length(class_levels_vec)==2){
        
        testeval_res_auc<-evaluation(class_res,measure = "auc")
        ####savetesteval_res_auc,file="testeval_res_auc.Rda")
        
        test_auc<-100*(mean(testeval_res_auc@score))
        
        test_acc_mat<-c(test_acc_mat,test_auc)
        
        
        
      }else{
        
        test_acc<-evaluation(class_res)
        
        test_acc<-100*(1-mean(test_acc@score))
        test_acc_mat<-c(test_acc_mat,test_acc) #cbind(best_kfold_acc,permkfold_acc,test_acc)
        
        
      }
      
    }
    
    
    if(length(class_levels_vec)==2){
      
      
      
      acc_mat<-cbind(emat1,test_acc_mat) #cbind(best_kfold_acc,permkfold_acc,test_acc)
      
      
      if(is.na(test_data_check)==TRUE){
        colnames(acc_mat)<-c("Training kfold CV Accuracy (AUC)","Training permuted kfold CV accuracy","Score", "Training accuracy(AUC)")
        
        print("Training set evaluation using selected features and the best classifier based on AUC measure")
      }else{
        colnames(acc_mat)<-c("Training kfold CV Accuracy (AUC)","Training permuted kfold CV accuracy","Score", "Test accuracy(AUC)")
        
        print("Test set evaluation using selected features and the best classifier based on AUC measure")
        
      }
      
      print(acc_mat[,c(1,4)])
      
      #rownames(acc_mat)<-best_classifier_name[1]
      
      write.table(acc_mat,file="Classification_evaluation_results_AUC.txt",sep="\t")
      
      
    }else{
      
      acc_mat<-cbind(emat1,test_acc_mat) #cbind(best_kfold_acc,permkfold_acc,test_acc)
      
      if(is.na(test_data_check)==TRUE){
        colnames(acc_mat)<-c("Training kfold CV Accuracy (Misclassification)","Training permuted kfold CV Accuracy (Misclassification)","Score", "Training accuracy (Misclassification)")
        
        print("Training set evaluation using selected features and the best classifier based on misclassification rate measure")
      }else{
        
        colnames(acc_mat)<-c("Training kfold CV Accuracy (Misclassification)","Training permuted kfold CV Accuracy (Misclassification)","Score", "Test accuracy (Misclassification)")
        
        print("Test set evaluation using selected features and the best classifier based on misclassification rate measure")
        
        
      }
      
#      print(acc_mat[,c(1,4)])
      
      #rownames(acc_mat)<-best_classifier_name[1]
      
      write.table(acc_mat,file="Classification_evaluation_results_misclassification.txt",sep="\t")
      
    }
    
    
    ##save(acc_mat,file="acc_mat.Rda")
    
    acc_mat<-acc_mat[,-c(3)]
    
    
    mainlab<-paste("Performance evaluation using classifiers and selected features",sep="")
    
    if(output.device.type!="pdf"){
      
      temp_filename_1<-"Figures/Barplot_classification_accuracy.png"
      
      png(temp_filename_1,width=plots.width,height=plots.height,res=plots.res,type=plots.type,units="in")
    }
    
    
    acc_mat1<-t(acc_mat)
    #xaxt="n",
    #
    w <- 0.1
    par(omd=c(0, 1-w, 0, 1),cex.main=0.7)
    if(length(class_levels_vec)==2){
      barplot(acc_mat1[c(1,3),],main=mainlab,ylab="Classification accuracy(%)",ylim=c(50,100),xpd=FALSE,beside=TRUE,col=barplot.col.opt,cex.axis=0.7,cex.names=0.7)
    }else{
      
      barplot(acc_mat1[c(1,3),],main=mainlab,ylab="Classification accuracy(%)",ylim=c(50,100),xpd=FALSE,beside=TRUE,col=barplot.col.opt,cex.axis=0.7,cex.names=0.7)
      
    }
    if(FALSE){
      if(length(class_levels_vec)==2){
        axis(side=1,at=seq(1,4),labels=c("kfold CV","Permuted kfold CV","Test set","AUC"),cex.axis=cex.plots)
      }else{
        axis(side=1,at=seq(1,3),labels=c("kfold CV","Permuted kfold CV","Test set"),cex.axis=cex.plots)
      }
    }
    #col = col_vec[1:length(t1)],
    le1<-(legend(par('usr')[2], par('usr')[4], bty='n', xpd=NA,c("Training (kfold)","Test"), pch = rep(19,2), pt.cex = 0.6, title = "",cex=0.7,col=barplot.col.opt))
    
    
    if(output.device.type!="pdf"){
      
      try(dev.off(),silent=TRUE)
    }
    
    
    try(dev.off(),silent=TRUE)
    print("number of features too small.")
  }
  
}
