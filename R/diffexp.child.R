diffexp.child <-
function(Xmat,Ymat,feature_table_file,parentoutput_dir,class_labels_file,num_replicates,feat.filt.thresh,summarize.replicates,summary.method,
                        summary.na.replacement,missing.val,rep.max.missing.thresh,
                        all.missing.thresh,group.missing.thresh,input.intensity.scale,
                        log2transform,medcenter,znormtransform,quantile_norm,lowess_norm,madscaling,TIC_norm,rangescaling,mstus,paretoscaling,sva_norm,eigenms_norm,vsn_norm,
                        normalization.method,rsd.filt.list,
                        pairedanalysis,featselmethod,fdrthresh,fdrmethod,cor.method,networktype,network.label.cex,abs.cor.thresh,cor.fdrthresh,kfold,pred.eval.method,feat_weight,globalcor,
                        target.metab.file,target.mzmatch.diff,target.rtmatch.diff,max.cor.num, samplermindex,pcacenter,pcascale,
                        numtrees,analysismode,net_node_colors,net_legend,svm_kernel,heatmap.col.opt,manhattanplot.col.opt,boxplot.col.opt,barplot.col.opt,sample.col.opt,lineplot.col.opt,scatterplot.col.opt,hca_type,alphacol,pls_vip_thresh,num_nodes,max_varsel,
                        pls_ncomp,pca.stage2.eval,scoreplot_legend,pca.global.eval,rocfeatlist,rocfeatincrement,rocclassifier,foldchangethresh,wgcnarsdthresh,WGCNAmodules,optselect,max_comp_sel,saveRda,legendlocation,degree_rank_method,
                        pca.cex.val,pca.ellipse,ellipse.conf.level,pls.permut.count,svm.acc.tolerance,limmadecideTests,pls.vip.selection,globalclustering,plots.res,plots.width,plots.height,plots.type,output.device.type,pvalue.thresh,individualsampleplot.col.opt,
                        pamr.threshold.select.max,mars.gcv.thresh,error.bar,cex.plots,modeltype,barplot.xaxis,lineplot.lty.option,match_class_dist,timeseries.lineplots,alphabetical.order,kegg_species_code,database,reference_set,target.data.annot,
                        add.pvalues=TRUE,add.jitter=TRUE,fcs.permutation.type,fcs.method,
                        fcs.min.hits,names_with_mz_time,ylab_text,xlab_text,boxplot.type,
                        degree.centrality.method,log2.transform.constant,balance.classes,
                        balance.classes.sizefactor,balance.classes.method,balance.classes.seed,
                        cv.perm.count=100,multiple.figures.perpanel=TRUE,labRow.value = TRUE, labCol.value = TRUE,
                        alpha.col=1,similarity.matrix,outlier.method,removeRda=TRUE,color.palette=c("journal"),
                        plot_DiNa_graph=FALSE,limma.contrasts.type=c("contr.sum","contr.treatment"),hca.cex.legend=0.7,differential.network.analysis.method,
                        plot.boxplots.raw=FALSE,vcovHC.type,ggplot.type1,facet.nrow,facet.ncol,...)
{
  
  
  
  #############
  options(warn=-1)
  
  lme.modeltype=modeltype
  remove_firstrun=FALSE #TRUE or FALSE
  run_number=1
  minmaxtransform=FALSE
  pca.CV=TRUE
  max_rf_var=5000
  alphacol=alpha.col
  
  hca.labRow.value=labRow.value
  hca.labCol.value=labCol.value
  logistic_reg=FALSE
  poisson_reg=FALSE
  goodfeats_allfields={}
  mwan_fdr={}
  targetedan_fdr={}
  data_m_fc_withfeats={}
  classlabels_orig={}
  robust.estimate=FALSE
  #alphabetical.order=FALSE
  
  analysistype="oneway"
  
  plot.ylab_text=ylab_text
  
  limmarobust=FALSE
  
  featselmethod<-unique(featselmethod)
  
  parentfeatselmethod=featselmethod
  
  if(featselmethod=="limmarobust"){
    
    featselmethod="limma"
    limmarobust=TRUE
  }else{
    
    if(featselmethod=="limma1wayrepeatrobust"){
      
      featselmethod="limma1wayrepeat"
      limmarobust=TRUE
    }else{
      if(featselmethod=="limma2wayrepeatrobust"){
        
        featselmethod="limma2wayrepeat"
        limmarobust=TRUE
      }else{
        
        if(featselmethod=="limma2wayrobust"){
          
          featselmethod="limma2way"
          limmarobust=TRUE
        }else{
          if(featselmethod=="limma1wayrobust"){
            
            featselmethod="limma1way"
            limmarobust=TRUE
          }
          
        }
        
      }
      
    }
    
  }
  
  
  if(normalization.method=="log2quantilenorm" || normalization.method=="log2quantnorm"){
    print("Performing log2 transformation and quantile normalization")
    log2transform=TRUE
    quantile_norm=TRUE
    
  }else{
    if(normalization.method=="log2transform"){
      print("Performing log2 transformation")
      log2transform=TRUE
    }else{
      if(normalization.method=="znormtransform"){
        print("Performing autoscaling")
        znormtransform=TRUE
        
      }else{
        if(normalization.method=="quantile_norm"){
          suppressMessages(library(limma))
          print("Performing quantile normalization")
          quantile_norm=TRUE
        }else{
          if(normalization.method=="lowess_norm"){
            suppressMessages(library(limma))
            print("Performing Cyclic Lowess normalization")
            lowess_norm=TRUE
          }else{
            
            if(normalization.method=="rangescaling"){
              print("Performing Range scaling")
              rangescaling=TRUE
            }else{
              if(normalization.method=="paretoscaling"){
                print("Performing Pareto scaling")
                paretoscaling=TRUE
              }else{
                
                if(normalization.method=="mstus"){
                  
                  print("Performing MS Total Useful Signal (MSTUS) normalization")
                  mstus=TRUE
                }else{
                  
                  if(normalization.method=="sva_norm"){
                    suppressMessages(library(sva))
                    print("Performing Surrogate Variable Analysis (SVA) normalization")
                    sva_norm=TRUE
                    log2transform=TRUE
                  }else{
                    if(normalization.method=="eigenms_norm"){
                      print("Performing EigenMS normalization")
                      eigenms_norm=TRUE
                      if(input.intensity.scale=="raw"){
                        log2transform=TRUE
                      }
                      
                    }else{
                      if(normalization.method=="vsn_norm"){
                        suppressMessages(library(limma))
                        print("Performing variance stabilizing normalization")
                        vsn_norm=TRUE
                        
                      }
                    }
                  }
                  
                }
                
                
              }
              
            }
            
          }
          
        }
        
        
      }
    }
    
    
  }
  
  if(input.intensity.scale=="log2"){
    
    log2transform=FALSE
  }
  
  rfconditional=FALSE
  
  print("############################")
  print("Feature selection method:")
  
  print(featselmethod)
  print("############################")
  if(featselmethod=="rf" | featselmethod=="RF"){
    
    suppressMessages(library(randomForest))
    suppressMessages(library(Boruta))
    
    featselmethod="RF"
    
    rfconditional=FALSE
  }else{
    
    if(featselmethod=="rfconditional" | featselmethod=="RFconditional" | featselmethod=="RFcond" | featselmethod=="rfcond"){
      
      
      suppressMessages(library(party))
      
      featselmethod="RF"
      
      rfconditional=TRUE
    }
  }
  
  
  if(featselmethod=="rf"){
    
    featselmethod="RF"
  }else{
    if(featselmethod=="mars"){
      
      
      suppressMessages(library(earth))
      featselmethod="MARS"
    }
  }
  
  if(featselmethod=="lmregrobust"){
    
    suppressMessages(library(sandwich))
    robust.estimate=TRUE
    featselmethod="lmreg"
  }else{
    
    if(featselmethod=="logitregrobust"){
      robust.estimate=TRUE
      suppressMessages(library(sandwich))
      featselmethod="logitreg"
    }else{
      
      if(featselmethod=="poissonregrobust"){
        robust.estimate=TRUE
        suppressMessages(library(sandwich))
        featselmethod="poissonreg"
      }
    }
  }
  
  if(featselmethod=="plsrepeat"){
    
    featselmethod="pls"
    pairedanalysis=TRUE
    
  }else{
    if(featselmethod=="splsrepeat"){
      featselmethod="spls"
      pairedanalysis=TRUE
    }else{
      if(featselmethod=="o1plsrepeat"){
        featselmethod="o1pls"
        pairedanalysis=TRUE
      }else{
        if(featselmethod=="o1splsrepeat"){
          featselmethod="o1spls"
          pairedanalysis=TRUE
        }
        
      }
      
    }
    
  }
  
  
  log2.fold.change.thresh_list<-rsd.filt.list
  if(featselmethod=="limma" | featselmethod=="limma2way" | featselmethod=="limma2wayrepeat" | featselmethod=="limma1wayrepeat"){
    
    if(analysismode=="regression"){
      
      stop("Invalid analysis mode. Please set analysismode=\"classification\".")
    }else{
      suppressMessages(library(limma))
      print("##############Level 1: Using LIMMA function to find differentially expressed metabolites###########")
    }
  }else{
    if(featselmethod=="RF"){
      
      print("##############Level 1: Using random forest function to find discriminatory metabolites###########")
      
      
    }else{
      if(featselmethod=="RFcond"){
        suppressMessages(library(party))
        print("##############Level 1: Using conditional random forest function to find discriminatory metabolites###########")
        #stop("Please use \"limma\", \"RF\", or \"MARS\".")
        
      }else{
        if(featselmethod=="MARS"){
          suppressMessages(library(earth))
          
          print("##############Level 1: Using MARS to find discriminatory metabolites###########")
          #log2.fold.change.thresh_list<-c(0)
        }else{
          
          if(featselmethod=="lmreg" | featselmethod=="logitreg" | featselmethod=="poissonreg" | featselmethod=="lm1wayanova" | featselmethod=="lm2wayanova" | featselmethod=="lm1wayanovarepeat" | featselmethod=="lm2wayanovarepeat" | featselmethod=="rfesvm" | 
             featselmethod=="wilcox" | featselmethod=="ttest" | featselmethod=="pamr" | featselmethod=="ttestrepeat" | featselmethod=="wilcoxrepeat" | featselmethod=="lmregrepeat"){
            print("##########Level 1: Finding discriminatory metabolites ###########")
            
            if(featselmethod=="logitreg"){
              
              featselmethod="lmreg"
              logistic_reg=TRUE
              poisson_reg=FALSE
            }else{
              
              if(featselmethod=="poissonreg"){
                poisson_reg=TRUE
                featselmethod="lmreg"
                logistic_reg=FALSE
              }else{
                logistic_reg=FALSE
                poisson_reg=FALSE
                
                if(featselmethod=="rfesvm"){
                  
                  suppressMessages(library(e1071))
                }else{
                  if(featselmethod=="pamr"){
                    
                    suppressMessages(library(pamr))
                  }else{
                    
                    if(featselmethod=="lm2wayanovarepeat" | featselmethod=="lm1wayanovarepeat"){
                      
                      suppressMessages(library(nlme))
                      suppressMessages(library(lsmeans))
                    }
                    
                  }
                  
                }
                
              }
              
            }
          }else{
            
            if(featselmethod=="pls" | featselmethod=="o1pls" | featselmethod=="o2pls" | featselmethod=="spls" | featselmethod=="spls1wayrepeat" | featselmethod=="spls2wayrepeat" | featselmethod=="pls2way" | featselmethod=="spls2way" | featselmethod=="o1spls" | featselmethod=="o2spls"){
              
              suppressMessages(library(mixOmics))
              suppressMessages(library(pls))
              suppressMessages(library(plsgenomics))
              
              print("##########Level 1: Finding discriminatory metabolites ###########")
              
            }else{
              
              
              stop("Invalid featselmethod specified.")
            }
            
          }
          
          #stop("Invalid featselmethod specified. Please use \"limma\", \"RF\", or \"MARS\".")
        }
        
      }
      
    }
    
  }
  ####################################################################################
  
  
  dir.create(parentoutput_dir,showWarnings=FALSE)
  parentoutput_dir1<-paste(parentoutput_dir,"/Stage1/",sep="")
  
  dir.create(parentoutput_dir1,showWarnings=FALSE)
  
  setwd(parentoutput_dir1)	
  if(is.na(Xmat[1])==TRUE){
    X<-read.table(feature_table_file,sep="\t",header=TRUE,stringsAsFactors=FALSE,check.names=FALSE)
    
  
    cnames<-colnames(X)
    
    cnames<- gsub(cnames,pattern="[\\s]*",replacement="",perl=TRUE)
    cnames<- gsub(cnames,pattern="[(|)|\\[|\\]]",replacement="",perl=TRUE)
    
    cnames<-gsub(cnames,pattern="\\||-|;|,|\\.",replacement="_",perl=TRUE)
    
    
    colnames(X)<-cnames
    
    cnames<-tolower(cnames)
    
    check_names<-grep(cnames,pattern="^name$")
    
    #if the Name column exists
    if(length(check_names)>0){
      
      if(check_names==1){
        
        check_names1<-grep(cnames,pattern="^mz$")
        check_names2<-grep(cnames,pattern="^time$")
        
        if(length(check_names1)<1 & length(check_names2)<1){
          mz<-seq(1.00001,nrow(X)+1,1)
          time<-seq(1.01,nrow(X)+1,1.00)
          check_ind<-gregexpr(cnames,pattern="^name$")
          check_ind<-which(check_ind>0)
          X<-as.data.frame(X)	
          Name<-as.character(X[,check_ind])
          if(length(which(duplicated(Name)==TRUE))>0){
            stop("Duplicate variable names are not allowed.")
          }
          X<-cbind(mz,time,X[,-check_ind])
          names_with_mz_time=cbind(Name,mz,time)		
          
          names_with_mz_time<-as.data.frame(names_with_mz_time)
          X<-as.data.frame(X)
          write.table(names_with_mz_time,file="Name_mz_time_mapping.txt",sep="\t",row.names=FALSE)
          
        }else{
          
          if(length(check_names1)>0 & length(check_names2)>0){
            
            check_ind<-gregexpr(cnames,pattern="^name$")
            check_ind<-which(check_ind>0)
            Name<-as.character(X[,check_ind])
            X<-X[,-check_ind]
            
            names_with_mz_time=cbind(Name,X$mz,X$time)
            colnames(names_with_mz_time)<-c("Name","mz","time")
            names_with_mz_time<-as.data.frame(names_with_mz_time)
            X<-as.data.frame(X)
            write.table(names_with_mz_time,file="Name_mz_time_mapping.txt",sep="\t",row.names=FALSE)
          }
        }		
        
      }
    }else{
      
      #mz time format
      check_names1<-grep(cnames[1],pattern="^mz$")
      check_names2<-grep(cnames[2],pattern="^time$")
      if(length(check_names1)<1 || length(check_names2)<1){
        stop("Invalid feature table format. First two columns should be mz and time. Please check example files.")			
      }
      
      X[,1]<-round(X[,1],5)
      X[,2]<-round(X[,2],2)
      
      mz_time<-paste(round(X[,1],5),"_",round(X[,2],2),sep="")
      if(length(which(duplicated(mz_time)==TRUE))>0){
        
        stop("Duplicate variable names are not allowed.")
      }
      Name<-mz_time
      names_with_mz_time=cbind(Name,X$mz,X$time)
      colnames(names_with_mz_time)<-c("Name","mz","time")
      names_with_mz_time<-as.data.frame(names_with_mz_time)
      X<-as.data.frame(X)
      write.table(names_with_mz_time,file="Name_mz_time_mapping.txt",sep="\t",row.names=FALSE)
      
      
    } 
    
    X[,1]<-round(X[,1],5)
    X[,2]<-round(X[,2],2)
    
    Xmat<-t(X[,-c(1:2)])
    
    rownames(Xmat)<-colnames(X[,-c(1:2)])
    
    Xmat<-as.data.frame(Xmat)
    
    colnames(Xmat)<-names_with_mz_time$Name
    
  }else{
    X<-Xmat
    
    
    cnames<-colnames(X)
    
   
    cnames<- gsub(cnames,pattern="[\\s]*",replacement="",perl=TRUE)
    cnames<- gsub(cnames,pattern="[(|)|\\[|\\]]",replacement="",perl=TRUE)
    
    cnames<-gsub(cnames,pattern="\\||-|;|,|\\.",replacement="_",perl=TRUE)
    
    colnames(X)<-cnames
    
    
    cnames<-tolower(cnames)
    
    check_names<-grep(cnames,pattern="^name$")
    
    
    if(length(check_names)>0){
      
      if(check_names==1){
        
        check_names1<-grep(cnames,pattern="^mz$")
        check_names2<-grep(cnames,pattern="^time$")
        
        
        if(length(check_names1)<1 & length(check_names2)<1){
          mz<-seq(1.00001,nrow(X)+1,1)
          time<-seq(1.01,nrow(X)+1,1.00)
          check_ind<-gregexpr(cnames,pattern="^name$")
          check_ind<-which(check_ind>0)
          X<-as.data.frame(X)
          
          
          Name<-as.character(X[,check_ind])
          
          
          X<-cbind(mz,time,X[,-check_ind])
          names_with_mz_time=cbind(Name,mz,time)
          
          names_with_mz_time<-as.data.frame(names_with_mz_time)
          X<-as.data.frame(X)
          
          # print(getwd())
          
          
          write.table(names_with_mz_time,file="Name_mz_time_mapping.txt",sep="\t",row.names=FALSE)
          
        }else{
          
          if(length(check_names1)>0 & length(check_names2)>0){
            
            check_ind<-gregexpr(cnames,pattern="^name$")
            check_ind<-which(check_ind>0)
            Name<-as.character(X[,check_ind])
            X<-X[,-check_ind]
            names_with_mz_time=cbind(Name,X$mz,X$time)
            colnames(names_with_mz_time)<-c("Name","mz","time")
            names_with_mz_time<-as.data.frame(names_with_mz_time)
            X<-as.data.frame(X)
            write.table(names_with_mz_time,file="Name_mz_time_mapping.txt",sep="\t",row.names=FALSE)
          }
        }
        
      }
    }else{
      
      
      
      
      check_names1<-grep(cnames[1],pattern="^mz$")
      check_names2<-grep(cnames[2],pattern="^time$")
      if(length(check_names1)<1 || length(check_names2)<1){
        stop("Invalid feature table format. First two columns should be mz and time. Please check example files.")
      }
      
      X[,1]<-round(X[,1],5)
      X[,2]<-round(X[,2],3)
      
      mz_time<-paste(round(X[,1],5),"_",round(X[,2],3),sep="")
      if(length(which(duplicated(mz_time)==TRUE))>0){
        
        stop("Duplicate variable names are not allowed.")
      }
      Name<-mz_time
      names_with_mz_time=cbind(Name,X$mz,X$time)
      colnames(names_with_mz_time)<-c("Name","mz","time")
      names_with_mz_time<-as.data.frame(names_with_mz_time)
      X<-as.data.frame(X)
      write.table(names_with_mz_time,file="Name_mz_time_mapping.txt",sep="\t",row.names=FALSE)
      
    }
    
    
    
    Xmat<-t(X[,-c(1:2)])
    
    rownames(Xmat)<-colnames(X[,-c(1:2)])
    Xmat<-as.data.frame(Xmat)
    colnames(Xmat)<-names_with_mz_time$Name
  }
  
  
  
  ####saveXmat,file="Xmat.Rda")
  
  if(analysismode=="regression")
  {
    
    #log2.fold.change.thresh_list<-c(0)
    
    print("Performing regression analysis")
    if(is.na(Ymat[1])==TRUE){
      classlabels<-read.table(class_labels_file,sep="\t",header=TRUE)
    
      
      Ymat<-classlabels			
    }else{
      classlabels<-Ymat
      
    }
    
    
    
    classlabels[,1]<- gsub(classlabels[,1],pattern="[\\s]*",replacement="",perl=TRUE)
    classlabels[,1]<- gsub(classlabels[,1],pattern="[(|)|\\[|\\]]",replacement="",perl=TRUE)
    
    classlabels[,1]<-gsub(classlabels[,1],pattern="\\||-|;|,|\\.",replacement="_",perl=TRUE)
    
    #classlabels[,1]<-gsub(classlabels[,1],pattern=" |-",replacement=".")
   # Ymat[,1]<-gsub(Ymat[,1],pattern=" |-",replacement=".")
    
    Ymat<-classlabels
    
    classlabels_orig<-classlabels
    classlabels_sub<-classlabels
    class_labels_levels<-c("A")
    
    if(featselmethod=="lmregrepeat" || featselmethod=="splsrepeat" || featselmethod=="plsrepeat" || featselmethod=="spls" || featselmethod=="pls" || featselmethod=="o1pls" || featselmethod=="o1splsrepeat"){
      if(pairedanalysis==TRUE){
        colnames(classlabels)<-c("SampleID","SubjectNum",paste("Response",sep=""))
        
        
        #Xmat<-chocolate[,1]
        Xmat_temp<-Xmat #t(Xmat)
        Xmat_temp<-cbind(classlabels,Xmat_temp)
        
        #Xmat_temp<-Xmat_temp[order(Xmat_temp[,3],Xmat_temp[,2]),]
        
        cnames<-colnames(Xmat_temp)
        
        factor_lastcol<-grep("^Response", cnames)
        
        classlabels<-Xmat_temp[,c(1:factor_lastcol[length(factor_lastcol)])]
        
        subject_inf<-classlabels[,2]
        classlabels<-classlabels[,-c(2)]
        
        Xmat<-Xmat_temp[,-c(1:factor_lastcol[length(factor_lastcol)])]
      }
      
    }
    
    
    classlabels<-as.data.frame(classlabels)
    classlabels_response_mat<-classlabels[,-c(1)]
    classlabels_response_mat<-as.data.frame(classlabels_response_mat)
    Ymat<-classlabels
    Ymat<-as.data.frame(Ymat)
    
    
    
    
    rnames_xmat<-as.character(rownames(Xmat))
    rnames_ymat<-as.character(Ymat[,1])
    
    
    
    if(length(which(duplicated(rnames_ymat)==TRUE))>0){
      
      stop("Duplicate sample IDs are not allowed. Please represent replicates by _1,_2,_3.")
    }
    
    check_ylabel<-regexpr(rnames_ymat[1],pattern="^[0-9]*",perl=TRUE)
    check_xlabel<-regexpr(rnames_xmat[1],pattern="^X[0-9]*",perl=TRUE)
    
    
    if(length(check_ylabel)>0 && length(check_xlabel)>0){
      if(attr(check_ylabel,"match.length")>0 && attr(check_xlabel,"match.length")>0){
        
        rnames_ymat<-paste("X",rnames_ymat,sep="")
      }
    }
    
    match_names<-match(rnames_xmat,rnames_ymat)
    
    bad_colnames<-length(which(is.na(match_names)==TRUE))
    
   # save(rnames_xmat,rnames_ymat,Xmat,Ymat,file="debugnames.Rda")
    #   print("Check here2")
    #if(is.na()==TRUE){
    
    bool_names_match_check<-all(rnames_xmat==rnames_ymat)
    
    if(bad_colnames>0 | bool_names_match_check==FALSE){
      print("Sample names do not match between feature table and class labels files.\n Please try replacing any \"-\" with \".\" in sample names.")
      print("Sample names in feature table")
      print(head(rnames_xmat))
      print("Sample names in classlabels file")
      
      print(head(rnames_ymat))
      stop("Sample names do not match between feature table and class labels files.\n Please try replacing any \"-\" with \".\" in sample names. Please try again.")
    }
    
    Xmat<-t(Xmat)
    Xmat<-cbind(X[,c(1:2)],Xmat)
    Xmat<-as.data.frame(Xmat)
    
    rownames(Xmat)<-names_with_mz_time$Name
    
    if(is.na(all(diff(match(rnames_xmat,rnames_ymat))))==FALSE){
      if(all(diff(match(rnames_xmat,rnames_ymat)) > 0)==TRUE){
        
        setwd("../")
        
        
        #data preprocess regression
        data_matrix<-data_preprocess(Xmat=Xmat,Ymat=Ymat,feature_table_file=feature_table_file,parentoutput_dir=parentoutput_dir,class_labels_file=NA,num_replicates=num_replicates,feat.filt.thresh=NA,summarize.replicates=summarize.replicates,summary.method=summary.method,
                                     all.missing.thresh=all.missing.thresh,group.missing.thresh=NA,
                                     log2transform=log2transform,medcenter=medcenter,znormtransform=znormtransform,,quantile_norm=quantile_norm,lowess_norm=lowess_norm,
                                     rangescaling=rangescaling,paretoscaling=paretoscaling,mstus=mstus,sva_norm=sva_norm,eigenms_norm=eigenms_norm,
                                     vsn_norm=vsn_norm,madscaling=madscaling,missing.val=0,samplermindex=NA, rep.max.missing.thresh=rep.max.missing.thresh,
                                     summary.na.replacement=summary.na.replacement,featselmethod=featselmethod,TIC_norm=TIC_norm,normalization.method=normalization.method,
                                     input.intensity.scale=input.intensity.scale,log2.transform.constant=log2.transform.constant,alphabetical.order=alphabetical.order)
        
      }
    }else{
      
      #print(diff(match(rnames_xmat,rnames_ymat)))
      stop("Orders of feature table and classlabels do not match")
    }
    
    
  }else{
    if(analysismode=="classification")
    {
      
      analysistype="oneway"
      
      classlabels_sub<-NA
      
      if(featselmethod=="limma2way" | featselmethod=="lm2wayanova" | featselmethod=="spls2way"){
        analysistype="twoway"
      }else{
        
        if(featselmethod=="limma2wayrepeat" | featselmethod=="lm2wayanovarepeat" | featselmethod=="spls2wayrepeat"){
          analysistype="twowayrepeat"
          pairedanalysis=TRUE
        }else{
          
          if(featselmethod=="limma1wayrepeat" | featselmethod=="lm1wayanovarepeat" | featselmethod=="spls1wayrepeat" |  featselmethod=="lmregrepeat"){
            analysistype="onewayrepeat"
            pairedanalysis=TRUE
          }
          
        }
        
      }
      
      

      
      
      if(is.na(Ymat)==TRUE){
        classlabels<-read.table(class_labels_file,sep="\t",header=TRUE)
        
 
        Ymat<-classlabels
        
      }else{
        classlabels<-Ymat
        
      }
      
      classlabels[,1]<- gsub(classlabels[,1],pattern="[\\s]*",replacement="",perl=TRUE)
      classlabels[,1]<- gsub(classlabels[,1],pattern="[(|)|\\[|\\]]",replacement="",perl=TRUE)
      
      classlabels[,1]<-gsub(classlabels[,1],pattern="\\||-|;|,|\\.",replacement="_",perl=TRUE)
      
      #classlabels[,1]<-gsub(classlabels[,1],pattern=" |-",replacement=".")
      # Ymat[,1]<-gsub(Ymat[,1],pattern=" |-",replacement=".")
      
      Ymat<-classlabels
      
     # classlabels[,1]<-gsub(classlabels[,1],pattern=" |-",replacement=".")
      Ymat[,1]<-gsub(Ymat[,1],pattern=" |-",replacement=".")
      
      print(paste("Number of samples in class labels file:",dim(Ymat)[1],sep=""))
      print(paste("Number of samples in feature table:",dim(Xmat)[1],sep=""))
      
      if(dim(Ymat)[1]!=(dim(Xmat)[1]))
      {
        
        stop("Number of samples are different in feature table and class labels file.")
      }
      
      
      if(fdrmethod=="none"){
        
        fdrthresh=pvalue.thresh
      }
      
      if(featselmethod=="limma" | featselmethod=="limma2way" | featselmethod=="limma2wayrepeat" | featselmethod=="limma1way" | featselmethod=="limma1wayrepeat" | 
         featselmethod=="MARS" | featselmethod=="RF" | featselmethod=="pls" | featselmethod=="o1pls" | featselmethod=="o2pls" | featselmethod=="lmreg" | featselmethod=="logitreg" | 
         featselmethod=="spls" | featselmethod=="pls1wayrepeat" | featselmethod=="spls1wayrepeat" | featselmethod=="pls2wayrepeat" | 
         featselmethod=="spls2wayrepeat" | featselmethod=="pls2way" | featselmethod=="spls2way" | featselmethod=="o1spls" | 
         featselmethod=="o2spls" | featselmethod=="lm1wayanova" | featselmethod=="lm2wayanova" | featselmethod=="lm1wayanovarepeat" | 
         featselmethod=="lm2wayanovarepeat" | featselmethod=="rfesvm" | featselmethod=="wilcox" | featselmethod=="ttest" | 
         featselmethod=="pamr" | featselmethod=="ttestrepeat" | featselmethod=="poissonreg" | featselmethod=="wilcoxrepeat" | featselmethod=="lmregrepeat")
      {
        #analysismode="classification"
        
        save(classlabels,file="thisclasslabels.Rda")
        #if(is.na(Ymat)==TRUE)
        {
          #classlabels<-read.table(class_labels_file,sep="\t",header=TRUE)
          
          if(analysismode=="classification"){
            if(featselmethod=="lmreg" | featselmethod=="logitreg" | featselmethod=="poissonreg")
            {
              if(alphabetical.order==FALSE){
                classlabels[,2] <- factor(classlabels[,2], levels=unique(classlabels[,2]))
              }
              
              levels_classA<-levels(factor(classlabels[,2]))
              for(l1 in levels_classA){
                g1<-grep(x=l1,pattern="[0-9]")
                
                if(length(g1)>0){
                  #stop("Class labels or factor levels should not have any numbers.")
                }
              }
              
            }else{
              
              if(featselmethod=="lmregrepeat"){
                if(alphabetical.order==FALSE){
                  classlabels[,3] <- factor(classlabels[,3], levels=unique(classlabels[,3]))
                }
                
                levels_classA<-levels(factor(classlabels[,3]))
                for(l1 in levels_classA){
                  g1<-grep(x=l1,pattern="[0-9]")
                  
                  if(length(g1)>0){
                    #stop("Class labels or factor levels should not have any numbers.")
                  }
                }
                  
              }else{
              for(c1 in 2:dim(classlabels)[2]){
                
                if(alphabetical.order==FALSE){
                  classlabels[,c1] <- factor(classlabels[,c1], levels=unique(classlabels[,c1]))
                }
                levels_classA<-levels(factor(classlabels[,c1]))
                for(l1 in levels_classA){
                  g1<-grep(x=l1,pattern="[0-9]")
                  
                  if(length(g1)>0){
                    #stop("Class labels or factor levels should not have any numbers.")
                  }
                  
                }
                
              }
            }
            }
          }
          
          
          classlabels_orig<-classlabels
          
          
          if(featselmethod=="limma1way"){
            
            featselmethod="limma"
          }
          
          
          
          
          # | featselmethod=="limma1wayrepeat"
          if(featselmethod=="limma" | featselmethod=="limma1way" | featselmethod=="MARS" | featselmethod=="RF" | featselmethod=="pls" | featselmethod=="o1pls" | featselmethod=="o2pls" | featselmethod=="lmreg" | 
             featselmethod=="logitreg" | featselmethod=="spls" | featselmethod=="o1spls" | featselmethod=="o2spls" | featselmethod=="rfesvm" | featselmethod=="pamr" | 
             featselmethod=="poissonreg" | featselmethod=="ttest" | featselmethod=="wilcox" | featselmethod=="lm1wayanova")
          {
            
            
            if(featselmethod=="lmreg" | featselmethod=="logitreg" | featselmethod=="poissonreg")
            {
              factor_inf<-classlabels[,-c(1)]
              factor_inf<-as.data.frame(factor_inf)
              #print(factor_inf)
              
              classlabels_orig<-colnames(classlabels[,-c(1)])
              colnames(classlabels)<-c("SampleID",paste("Factor",seq(1,dim(factor_inf)[2]),sep=""))
              
              Xmat_temp<-Xmat #t(Xmat)
              
              #print(Xmat_temp[1:2,1:3])
              Xmat_temp<-cbind(classlabels,Xmat_temp)
              #print("here")				
              
              if(alphabetical.order==TRUE){
              Xmat_temp<-Xmat_temp[order(Xmat_temp[,2]),]
              }else{
                if(analysismode=="classification"){
                    Xmat_temp[,2] <- factor(Xmat_temp[,2], levels=unique(Xmat_temp[,2]))
                }
              }
              
              cnames<-colnames(Xmat_temp)
              
              factor_lastcol<-grep("^Factor", cnames)
              
              classlabels<-Xmat_temp[,c(1:factor_lastcol[length(factor_lastcol)])]
              
              levels_classA<-levels(factor(classlabels[,2]))
              
             
              print(paste("Factor 1 levels: ",paste(levels_classA,collapse=","),sep=""))
              
            
              classlabels_class<-as.factor(classlabels[,2])
              classtable1<-table(classlabels[,2])
              
              classlabels_xyplots<-classlabels
              #classlabels_orig<-classlabels
              
              # classlabels_orig<-classlabels_orig[seq(1,dim(classlabels)[1],num_replicates),]
              classlabels<-cbind(as.data.frame(classlabels[,1]),as.data.frame(classlabels_class))
              classlabels_xyplots<-classlabels
              
              rownames(Xmat_temp)<-as.character(Xmat_temp[,1])
              Xmat<-Xmat_temp[,-c(1:factor_lastcol[length(factor_lastcol)])]
              
              classlabels_response_mat<-classlabels[,-c(1)]
              
              classlabels<-as.data.frame(classlabels)
              
              #keeps the class order as in the input file
              if(alphabetical.order==FALSE){
                classlabels[,2] <- factor(classlabels[,2], levels=unique(classlabels[,2]))
              }
              classlabels_response_mat<-classlabels[,-c(1)]
              
              classlabels_response_mat<-as.data.frame(classlabels_response_mat)
              
              colnames(classlabels_response_mat)<-as.character(classlabels_orig)
              
              
              
              
              Ymat<-classlabels
              
              classlabels_orig<-classlabels
              
            }else
            {
              
              
              
              if(dim(classlabels)[2]>2){
                if(pairedanalysis==FALSE){	
                  #print("Invalid classlabels file format. Correct format: \nColumnA: SampleID\nColumnB: Class")
                  print("Using the first column as sample ID and second column as Class. Ignoring additional columns.")
                  classlabels<-classlabels[,c(1:2)]
                }
              }
              
              if(analysismode=="classification")
              {
                factor_inf<-classlabels[,-c(1)]
                factor_inf<-as.data.frame(factor_inf)
                
                colnames(classlabels)<-c("SampleID",paste("Factor",seq(1,dim(factor_inf)[2]),sep=""))
                
                Xmat_temp<-Xmat #t(Xmat)
                
                Xmat_temp<-cbind(classlabels,Xmat_temp)
                
                #  	##save(Xmat_temp,file="Xmat_temp.Rda")
                rownames(Xmat_temp)<-as.character(Xmat_temp[,1])						
                
                   
                if(alphabetical.order==TRUE){
                  
                  Xmat_temp<-Xmat_temp[order(Xmat_temp[,2]),]
                  
              
                }else{
                  Xmat_temp[,2] <- factor(Xmat_temp[,2], levels=unique(Xmat_temp[,2]))
                  
                 
                }
                
                cnames<-colnames(Xmat_temp)
                
                factor_lastcol<-grep("^Factor", cnames)
                
                classlabels<-Xmat_temp[,c(1:factor_lastcol[length(factor_lastcol)])]
                
                
                Xmat<-Xmat_temp[,-c(1:factor_lastcol[length(factor_lastcol)])]
                
                
                levels_classA<-levels(factor(classlabels[,2]))
                
               
                print(paste("Factor 1 levels: ",paste(levels_classA,collapse=","),sep=""))
                
               
                
                classlabels_class<-as.factor(classlabels[,2])
                
                classtable1<-table(classlabels[,2])
                
                classlabels_xyplots<-classlabels
                #classlabels_orig<-classlabels
                
                # classlabels_orig<-classlabels_orig[seq(1,dim(classlabels)[1],num_replicates),]
                classlabels<-cbind(as.data.frame(classlabels[,1]),as.data.frame(classlabels_class))
                
                #rownames(Xmat)<-rownames(Xmat_temp)
                classlabels_xyplots<-classlabels
                
                classlabels_sub<-classlabels[,-c(1)]
               
                if(alphabetical.order==FALSE){
                  classlabels[,2] <- factor(classlabels[,2], levels=unique(classlabels[,2]))
                  
                  
                  if(dim(classlabels)[2]>2){
                    #classlabels[,3] <- factor(classlabels[,3], levels=unique(classlabels[,3]))
                    stop("Invalid classlabels format.")
                  }
                }
              }
              
              classlabels_response_mat<-classlabels[,-c(1)]
              
              
              classlabels<-as.data.frame(classlabels)
              classlabels_response_mat<-classlabels[,-c(1)]
              classlabels_response_mat<-as.data.frame(classlabels_response_mat)
              
              #classlabels[,1]<-as.factor(classlabels[,1])
              Ymat<-classlabels
              
              classlabels_orig<-classlabels
              
            }
            
            #print("here 2")
            
          }
          
          if(featselmethod=="limma1wayrepeat"){
            factor_inf<-classlabels[,-c(1:2)]
            factor_inf<-as.data.frame(factor_inf)
            
            # print("here")
            colnames(classlabels)<-c("SampleID","SubjectNum",paste("Factor",seq(1,length(factor_inf)),sep=""))
            
            
            #Xmat<-chocolate[,1]
            Xmat_temp<-Xmat #t(Xmat)
            Xmat_temp<-cbind(classlabels,Xmat_temp)
            
            if(alphabetical.order==TRUE){
                Xmat_temp<-Xmat_temp[order(Xmat_temp[,3],Xmat_temp[,2]),]
            }else{
              
              Xmat_temp[,3] <- factor(Xmat_temp[,3], levels=unique(Xmat_temp[,3]))
              
              
            }
            
            cnames<-colnames(Xmat_temp)
            
            factor_lastcol<-grep("^Factor", cnames)
            
            classlabels<-Xmat_temp[,c(1:factor_lastcol[length(factor_lastcol)])]
            
            if(alphabetical.order==FALSE){
             
              classlabels[,3] <- factor(classlabels[,3], levels=unique(classlabels[,3]))
              
              
            }
            
            subject_inf<-classlabels[,2]
            classlabels_sub<-classlabels[,-c(1)]
            subject_inf<-subject_inf[seq(1,dim(classlabels)[1],num_replicates)]
            classlabels<-classlabels[,-c(2)]
            
            levels_classA<-levels(factor(classlabels[,2]))
            
           
            print(paste("Factor 1 levels: ",paste(levels_classA,collapse=","),sep=""))
            
         
            classlabels_class<-as.factor(classlabels[,2])
            classtable1<-table(classlabels[,2])
            
            classlabels_xyplots<-classlabels
            #classlabels_orig<-classlabels
            
            # classlabels_orig<-classlabels_orig[seq(1,dim(classlabels)[1],num_replicates),]
            classlabels<-cbind(as.data.frame(classlabels[,1]),as.data.frame(classlabels_class))
            
            classlabels_xyplots<-classlabels
            
            Xmat<-Xmat_temp[,-c(1:factor_lastcol[length(factor_lastcol)])]
            
            classlabels_response_mat<-classlabels[,-c(1)]
            
            classlabels<-as.data.frame(classlabels)
            classlabels_response_mat<-classlabels[,-c(1)]
            classlabels_response_mat<-as.data.frame(classlabels_response_mat)
            Ymat<-classlabels
            
            
            
            if(featselmethod=="limma1wayrepeat"){	
              featselmethod="limma"
              pairedanalysis = TRUE
            }else{
              
              if(featselmethod=="spls1wayrepeat"){
                featselmethod="spls"
                pairedanalysis = TRUE
              }else{
                if(featselmethod=="pls1wayrepeat"){
                  featselmethod="pls"
                  pairedanalysis = TRUE
                }
              }
            }
            pairedanalysis = TRUE
            
          }
          
          
          
          if(featselmethod=="limma2way"){
            
            factor_inf<-classlabels[,-c(1)]
            factor_inf<-as.data.frame(factor_inf)
            
            colnames(classlabels)<-c("SampleID",paste("Factor",seq(1,dim(factor_inf)[2]),sep=""))
            
            Xmat_temp<-Xmat #t(Xmat)
            
            ####saveXmat,file="Xmat.Rda")
            
            ####saveclasslabels,file="Xmat_classlabels.Rda")
            
            if(dim(classlabels)[2]>2){
              
              
             # save(Xmat_temp,classlabels,file="Xmat_temp_limma.Rda")
              
              
              Xmat_temp<-cbind(classlabels,Xmat_temp)
              
              
              
              
              #		print(Xmat_temp[1:10,1:10])
              
             if(alphabetical.order==TRUE){
                  Xmat_temp<-Xmat_temp[order(Xmat_temp[,2],Xmat_temp[,3]),]
              }else{
                Xmat_temp[,2] <- factor(Xmat_temp[,2], levels=unique(Xmat_temp[,2]))
                Xmat_temp[,3] <- factor(Xmat_temp[,3], levels=unique(Xmat_temp[,3]))
                
              }
              #		print(Xmat_temp[1:10,1:10])
              cnames<-colnames(Xmat_temp)
              
              factor_lastcol<-grep("^Factor", cnames)
              
              classlabels<-Xmat_temp[,c(1:factor_lastcol[length(factor_lastcol)])]
              
              Xmat<-Xmat_temp[,-c(1:factor_lastcol[length(factor_lastcol)])]
              classlabels_sub<-classlabels[,-c(1)]
              
              classlabels_response_mat<-classlabels[,-c(1)]
              classlabels<-as.data.frame(classlabels)
              
              classlabels_response_mat<-as.data.frame(classlabels_response_mat)
              
              if(alphabetical.order==FALSE){
                classlabels[,2] <- factor(classlabels[,2], levels=unique(classlabels[,2]))
                classlabels[,3] <- factor(classlabels[,3], levels=unique(classlabels[,3]))
              }
              levels_classA<-levels(factor(classlabels[,2]))
              
              levels_classB<-levels(factor(classlabels[,3]))
              print(paste("Factor 1 levels: ",paste(levels_classA,collapse=","),sep=""))
              
              print(paste("Factor 2 levels: ",paste(levels_classB,collapse=","),sep=""))
              
              classlabels_class<-as.factor(classlabels[,2]):as.factor(classlabels[,3])
              classtable1<-table(classlabels[,2],classlabels[,3])
              
              classlabels_xyplots<-classlabels
              #classlabels_orig<-classlabels
              
              # classlabels_orig<-classlabels_orig[seq(1,dim(classlabels)[1],num_replicates),]
              classlabels<-cbind(as.data.frame(classlabels[,1]),as.data.frame(classlabels_class))
              Ymat<-classlabels
              
              
              #classlabels_response_mat<-classlabels[,-c(1)]
              classlabels<-as.data.frame(classlabels)
              #classlabels_response_mat<-classlabels[,-c(1)]
              #classlabels_response_mat<-as.data.frame(classlabels_response_mat)
              Ymat<-classlabels
              
              
              #classlabels_orig<-classlabels
              
              
              
            }
            else{
              stop("Only one factor specificied in the class labels file.")			
            }
            
          }
          
          
          
          if(featselmethod=="limma2wayrepeat"){
            factor_inf<-classlabels[,-c(1:2)]
            factor_inf<-as.data.frame(factor_inf)
            
            
            
            colnames(classlabels)<-c("SampleID","SubjectNum",paste("Factor",seq(1,dim(factor_inf)[2]),sep=""))
            
            Xmat_temp<-Xmat
            if(dim(classlabels)[2]>2)
            {
              
              levels_classA<-levels(factor(classlabels[,3]))
              
              if(length(levels_classA)>2){
                
                #stop("Factor 1 can only have two levels/categories. Factor 2 can have upto 6 levels. \nPlease rearrange the factors in your classlabels file.")
                #	classtemp<-classlabels[,3]
                #	classlabels[,3]<-classlabels[,4]
                #	classlabels[,4]<-classtemp
              }
              
              levels_classA<-levels(factor(classlabels[,3]))
              
              if(length(levels_classA)>2){
                #stop("Only one of the factors can have more than 2 levels/categories. \nPlease rearrange the factors in your classlabels file or use lm2wayanovarepeat.")
                #stop("Please select lm2wayanova or lm2wayanovarepeat option for greater than 2x2 designs.")
                stop("Factor 1 can only have two levels/categories. Factor 2 can have upto 6 levels. \nPlease rearrange the factors in your classlabels file. Or use lm2wayanova option.")
              }
              
              levels_classB<-levels(factor(classlabels[,4]))
              if(length(levels_classB)>7){
                #stop("Only one of the factors can have more than 2 levels/categories. \nPlease rearrange the factors in your classlabels file or use lm2wayanova.")
                
                stop("Please select lm2wayanovarepeat option for greater than 2x7 designs.")
              }							
              
              
              Xmat_temp<-cbind(classlabels,Xmat_temp)
              
              
              
              if(alphabetical.order==TRUE){
                #Xmat_temp<-Xmat_temp[order(Xmat_temp[,2],Xmat_temp[,3]),]
                
                Xmat_temp<-Xmat_temp[order(Xmat_temp[,3],Xmat_temp[,4],Xmat_temp[,2]),]
                
              }else{
                Xmat_temp[,4] <- factor(Xmat_temp[,4], levels=unique(Xmat_temp[,4]))
                Xmat_temp[,3] <- factor(Xmat_temp[,3], levels=unique(Xmat_temp[,3]))
                
              }
              
            
              cnames<-colnames(Xmat_temp)
              
              factor_lastcol<-grep("^Factor", cnames)
              
              classlabels<-Xmat_temp[,c(1:factor_lastcol[length(factor_lastcol)])]
              
              classlabels_sub<-classlabels[,-c(1)]
              subject_inf<-classlabels[,2]
              classlabels<-classlabels[,-c(2)]
              
              classlabels_response_mat<-classlabels[,-c(1)]
              classlabels<-as.data.frame(classlabels)
              
              classlabels_response_mat<-as.data.frame(classlabels_response_mat)
              
              classlabels_xyplots<-classlabels
              subject_inf<-subject_inf[seq(1,dim(classlabels)[1],num_replicates)]
              
              #write.table(classlabels,file="organized_classlabelsA1.txt",sep="\t",row.names=FALSE)
              Xmat<-Xmat_temp[,-c(1:factor_lastcol[length(factor_lastcol)])]
              
              
              #write.table(Xmat_temp,file="organized_featuretableA1.txt",sep="\t",row.names=TRUE)
              
              if(alphabetical.order==FALSE){
                classlabels[,2] <- factor(classlabels[,2], levels=unique(classlabels[,2]))
                classlabels[,3] <- factor(classlabels[,3], levels=unique(classlabels[,3]))
                
              }
              
              levels_classA<-levels(factor(classlabels[,2]))
              
              levels_classB<-levels(factor(classlabels[,3]))
              print(paste("Factor 1 levels: ",paste(levels_classA,collapse=","),sep=""))
              
              print(paste("Factor 2 levels: ",paste(levels_classB,collapse=","),sep=""))
              classlabels_class<-as.factor(classlabels[,2]):as.factor(classlabels[,3])
              classtable1<-table(classlabels[,2],classlabels[,3])
              
              #classlabels_orig<-classlabels
              #classlabels<-cbind(as.character(classlabels[,1]),as.character(classlabels_class))
              
              classlabels<-cbind(as.data.frame(classlabels[,1]),as.data.frame(classlabels_class))
              
              Ymat<-classlabels
              
              
              print("Class labels file limma2wayrep:")
              print(head(classlabels))
              #rownames(Xmat)<-as.character(classlabels[,1])
              
              #write.table(classlabels,file="organized_classlabels.txt",sep="\t",row.names=FALSE)
              
              Xmat1<-cbind(classlabels,Xmat)
              #write.table(Xmat1,file="organized_featuretable.txt",sep="\t",row.names=TRUE)
              
              featselmethod="limma2way"
              pairedanalysis = TRUE
              
            }
            else{
              stop("Only one factor specificied in the class labels file.")			
            }
          }
          
          
          
        }
        
        classlabels<-as.data.frame(classlabels)
        
        
        
        
        
        
        
        
        if(featselmethod=="lm2wayanova" | featselmethod=="pls2way" | featselmethod=="spls2way"){
          
          analysismode="classification"
          
          #classlabels<-read.table(class_labels_file,sep="\t",header=TRUE)
          
          if(is.na(Ymat)==TRUE){
            classlabels<-read.table(class_labels_file,sep="\t",header=TRUE)
            Ymat<-classlabels
            
          }else{
            classlabels<-Ymat
            
          }
          
          
          
          #cnames[2]<-"Factor1"
          
          cnames<-colnames(classlabels)
          
          factor_inf<-classlabels[,-c(1)]
          factor_inf<-as.data.frame(factor_inf)
          
          colnames(classlabels)<-c("SampleID",paste("Factor",seq(1,dim(factor_inf)[2]),sep=""))
          
          analysismode="classification"
          
          Xmat_temp<-Xmat #t(Xmat)
          
        #  save(Xmat_temp,classlabels,file="Xmat_temp_lm2way.Rda")
          
          
          
          Xmat_temp<-cbind(classlabels,Xmat_temp)
          
          
          rnames_xmat<-rownames(Xmat)
          rnames_ymat<-as.character(Ymat[,1])
          
          
          # ###saveXmat_temp,file="Xmat_temp.Rda")
          
          
            if(featselmethod=="lm2wayanova" | featselmethod=="pls2way" | featselmethod=="spls2way"){
              
              if(alphabetical.order==TRUE){
                Xmat_temp<-Xmat_temp[order(Xmat_temp[,2],Xmat_temp[,3]),]
              }
              
            }
          
          cnames<-colnames(Xmat_temp)
          
          factor_lastcol<-grep("^Factor", cnames)
          
      #   save(Xmat_temp,classlabels,factor_lastcol,file="debudsort.Rda")
          
          if(alphabetical.order==FALSE){
            
            Xmat_temp[,2] <- factor(Xmat_temp[,2], levels=unique(Xmat_temp[,2]))
            Xmat_temp[,3] <- factor(Xmat_temp[,3], levels=unique(Xmat_temp[,3]))
            
            classlabels<-Xmat_temp[,c(1:factor_lastcol[length(factor_lastcol)])]
            
            classlabels[,2] <- factor(classlabels[,2], levels=unique(classlabels[,2]))
            classlabels[,3] <- factor(classlabels[,3], levels=unique(classlabels[,3]))
            
          }else{
            
            classlabels<-Xmat_temp[,c(1:factor_lastcol[length(factor_lastcol)])]
            
          }
          
          levels_classA<-levels(factor(classlabels[,2]))
          levels_classB<-levels(factor(classlabels[,3]))
          
          print(paste("Factor 1 levels: ",paste(levels_classA,collapse=","),sep=""))
          print(paste("Factor 2 levels: ",paste(levels_classB,collapse=","),sep=""))
          
          classlabels_sub<-classlabels[,-c(1)]
          
          classlabels_response_mat<-classlabels[,-c(1)]
          
          
          Ymat<-classlabels
          
          classlabels_orig<-classlabels
          
          #Xmat<-Xmat_temp[,-c(1:factor_lastcol[length(factor_lastcol)])]
          ###save(Xmat,file="Xmat2.Rda")
          if(featselmethod=="lm2wayanova" | featselmethod=="pls2way" | featselmethod=="spls2way"){
            
            
            classlabels_class<-as.factor(classlabels[,2]):as.factor(classlabels[,3])
            classtable1<-table(classlabels[,2],classlabels[,3])
            
            classlabels_xyplots<-classlabels
            #classlabels_orig<-classlabels
            
            # classlabels_orig<-classlabels_orig[seq(1,dim(classlabels)[1],num_replicates),]
            classlabels<-cbind(as.data.frame(classlabels[,1]),as.data.frame(classlabels_class))
            Ymat<-classlabels
            if(featselmethod=="pls2way"){
              featselmethod="pls"
            }else{
              
              if(featselmethod=="spls2way"){
                featselmethod="spls"
              }
            }
            
          }
          
          
          # write.table(classlabels,file="organized_classlabelsB.txt",sep="\t",row.names=FALSE)
          Xmat<-Xmat_temp[,-c(1:factor_lastcol[length(factor_lastcol)])]
          
          
          #write.table(Xmat_temp,file="organized_featuretableA.txt",sep="\t",row.names=TRUE)
          #write.table(classlabels,file="organized_classlabelsA.txt",sep="\t",row.names=FALSE)
          
          
        }
        
        if(featselmethod=="lm1wayanovarepeat" | featselmethod=="lm2wayanovarepeat" | featselmethod=="pls1wayrepeat" | featselmethod=="spls1wayrepeat" | featselmethod=="pls2wayrepeat" | 
           featselmethod=="spls2wayrepeat" | featselmethod=="ttestrepeat" | featselmethod=="wilcoxrepeat" | featselmethod=="lmregrepeat"){
          
    
          #analysismode="classification"
          pairedanalysis=TRUE
          #							classlabels<-read.table(class_labels_file,sep="\t",header=TRUE)
          
          
          if(is.na(Ymat)==TRUE){
            classlabels<-read.table(class_labels_file,sep="\t",header=TRUE)
            Ymat<-classlabels
            
          }else{
            classlabels<-Ymat
            
          }
          
          
          cnames<-colnames(classlabels)
          
          factor_inf<-classlabels[,-c(1:2)]
          factor_inf<-as.data.frame(factor_inf)
          
          colnames(classlabels)<-c("SampleID","SubjectNum",paste("Factor",seq(1,dim(factor_inf)[2]),sep=""))
          
          classlabels_orig<-classlabels
          #Xmat<-chocolate[,1]
          Xmat_temp<-Xmat #t(Xmat)
          Xmat_temp<-cbind(classlabels,Xmat_temp)
          
          pairedanalysis=TRUE
          if(featselmethod=="lm1wayanovarepeat" | featselmethod=="pls1wayrepeat" | featselmethod=="spls1wayrepeat" | featselmethod=="ttestrepeat" | featselmethod=="wilcoxrepeat" | featselmethod=="lmregrepeat"){
            
            if(alphabetical.order==TRUE){
            Xmat_temp<-Xmat_temp[order(Xmat_temp[,3],Xmat_temp[,2]),]
            }else{
              Xmat_temp[,3] <- factor(Xmat_temp[,3], levels=unique(Xmat_temp[,3]))
              
            }
            
            
           
            
            
            cnames<-colnames(Xmat_temp)
            
            factor_lastcol<-grep("^Factor", cnames)
            
            classlabels<-Xmat_temp[,c(1:factor_lastcol[length(factor_lastcol)])]
            
            subject_inf<-classlabels[,2]
            
            subject_inf<-subject_inf[seq(1,dim(classlabels)[1],num_replicates)]
            
            classlabels_response_mat<-classlabels[,-c(1:2)]
            
            
            
            # classlabels_orig<-classlabels
            classlabels_sub<-classlabels[,-c(1)]
            
            if(alphabetical.order==FALSE){
              classlabels[,3] <- factor(classlabels[,3], levels=unique(classlabels[,3]))
            }
            
            levels_classA<-levels(factor(classlabels[,3]))
            print(paste("Factor 1 levels: ",paste(levels_classA,collapse=","),sep=""))
            
            
            classlabels<-classlabels[,-c(2)]
            
            if(alphabetical.order==FALSE){
              classlabels[,2] <- factor(classlabels[,2], levels=unique(classlabels[,2]))
            }
            
            classlabels_class<-classlabels[,2]
            
            classtable1<-table(classlabels[,2])
            
            
            #classlabels<-cbind(as.character(classlabels[,1]),as.character(classlabels_class))
            
            classlabels<-cbind(as.data.frame(classlabels[,1]),as.data.frame(classlabels_class))
            
            Ymat<-classlabels
            
            classlabels_xyplots<-classlabels
            
            
            # classlabels<-classlabels[seq(1,dim(classlabels)[1],num_replicates),]
            Ymat<-classlabels
            
            Xmat<-Xmat_temp[,-c(1:factor_lastcol[length(factor_lastcol)])]
            
            
            
            # write.table(Xmat_temp,file="organized_featuretableA.txt",sep="\t",row.names=FALSE)
            
            ####saveYmat,file="Ymat.Rda")
            #                                       ###saveXmat,file="Xmat.Rda")
            
            if(featselmethod=="spls1wayrepeat"){
              featselmethod="spls"
              
            }else{
              if(featselmethod=="pls1wayrepeat"){
                featselmethod="pls"
              }														
            }
            
            
            if(featselmethod=="wilcoxrepeat"){
              
              featselmethod=="wilcox"
              pairedanalysis=TRUE
            }
            
            if(featselmethod=="ttestrepeat"){
              
              featselmethod=="ttest"
              pairedanalysis=TRUE
            }
          }
          if(featselmethod=="lm2wayanovarepeat" | featselmethod=="pls2wayrepeat" | featselmethod=="spls2wayrepeat"){
            
            if(alphabetical.order==TRUE){
            Xmat_temp<-Xmat_temp[order(Xmat_temp[,3],Xmat_temp[,4],Xmat_temp[,2]),]
            }else{
              
              Xmat_temp[,3] <- factor(Xmat_temp[,3], levels=unique(Xmat_temp[,3]))
              Xmat_temp[,4] <- factor(Xmat_temp[,4], levels=unique(Xmat_temp[,4]))
            }
            
            
            
            
            cnames<-colnames(Xmat_temp)
            
            factor_lastcol<-grep("^Factor", cnames)
            
            classlabels<-Xmat_temp[,c(1:factor_lastcol[length(factor_lastcol)])]
            classlabels_sub<-classlabels[,-c(1)]
            
            subject_inf<-classlabels[,2]
            
            subject_inf<-subject_inf[seq(1,dim(classlabels)[1],num_replicates)]
            classlabels_response_mat<-classlabels[,-c(1:2)]
            
            Ymat<-classlabels
            
            
            
            classlabels_xyplots<-classlabels[,-c(2)]
            
            if(alphabetical.order==FALSE){
              classlabels[,4] <- factor(classlabels[,4], levels=unique(classlabels[,4]))
              classlabels[,3] <- factor(classlabels[,3], levels=unique(classlabels[,3]))
            }
            
            levels_classA<-levels(factor(classlabels[,3]))
            print(paste("Factor 1 levels: ",paste(levels_classA,collapse=","),sep=""))
            levels_classB<-levels(factor(classlabels[,4]))
            
            print(paste("Factor 2 levels: ",paste(levels_classB,collapse=","),sep=""))
            Ymat<-classlabels
            
            
            #print(head(classlabels))
            
            
            classlabels<-classlabels[,-c(2)]
            
            classlabels_class<-paste(classlabels[,2],":",classlabels[,3],sep="")
            
            classtable1<-table(classlabels[,2],classlabels[,3])
            
            
            #classlabels<-cbind(as.character(classlabels[,1]),as.character(classlabels_class))
            
            classlabels<-cbind(as.data.frame(classlabels[,1]),as.data.frame(classlabels_class))
            
            Ymat<-classlabels
            
            
            # write.table(classlabels,file="organized_classlabelsA1.txt",sep="\t",row.names=FALSE)
            Xmat<-Xmat_temp[,-c(1:factor_lastcol[length(factor_lastcol)])]
            
            
            
            
            #write.table(Xmat_temp,file="organized_featuretableA.txt",sep="\t",row.names=FALSE)
            #write.table(Xmat,file="organized_featuretableB1.txt",sep="\t",row.names=FALSE)
            pairedanalysis=TRUE
            if(featselmethod=="spls2wayrepeat"){
              featselmethod="spls"
              
            }
          }
      
          
          
        }
        
      }
      
      
      rownames(Xmat)<-as.character(Xmat_temp[,1])
      
     # save(Xmat,Xmat_temp,file="Xmat1.Rda")
      #save(Ymat,file="Ymat1.Rda")
      rnames_xmat<-rownames(Xmat)
      rnames_ymat<-as.character(Ymat[,1])
      
      
      if(length(which(duplicated(rnames_ymat)==TRUE))>0){
        
        stop("Duplicate sample IDs are not allowed. Please represent replicates by _1,_2,_3.")
      }
      
      check_ylabel<-regexpr(rnames_ymat[1],pattern="^[0-9]*",perl=TRUE)
      check_xlabel<-regexpr(rnames_xmat[1],pattern="^X[0-9]*",perl=TRUE)
      
      if(length(check_ylabel)>0 && length(check_xlabel)>0){
        if(attr(check_ylabel,"match.length")>0 && attr(check_xlabel,"match.length")>0){
          
          rnames_ymat<-paste("X",rnames_ymat,sep="") #gsub(rnames_ymat,pattern="\\.[0-9]*",replacement="")
          
          
        }
      }
      
      
      
      Xmat<-t(Xmat)
      
      
      
      colnames(Xmat)<-as.character(Ymat[,1])
      
      Xmat<-cbind(X[,c(1:2)],Xmat)
      
      Xmat<-as.data.frame(Xmat)
      Ymat<-as.data.frame(Ymat)
      
      
      match_names<-match(rnames_xmat,rnames_ymat)
      
      bad_colnames<-length(which(is.na(match_names)==TRUE))
      
      #print(match_names) 
      #if(is.na()==TRUE){
     #save(rnames_xmat,rnames_ymat,Xmat,Ymat,file="debugnames.Rda")
      
      bool_names_match_check<-all(rnames_xmat==rnames_ymat)
      
      if(bad_colnames>0 | bool_names_match_check==FALSE){
        
     # if(bad_colnames>0){
        print("Sample names do not match between feature table and class labels files.\n Please try replacing any \"-\" with \".\" in sample names.")
        print("Sample names in feature table")
        print(head(rnames_xmat))
        print("Sample names in classlabels file")
        
        print(head(rnames_ymat))
        stop("Sample names do not match between feature table and class labels files.\n Please try replacing any \"-\" with \".\" in sample names. Please try again.")
      }
      
      if(is.na(all(diff(match(rnames_xmat,rnames_ymat))))==FALSE){
        if(all(diff(match(rnames_xmat,rnames_ymat)) > 0)==TRUE){
          
          setwd("../")
          
          
          
         #save(Xmat,Ymat,names_with_mz_time,feature_table_file,parentoutput_dir,class_labels_file,num_replicates,feat.filt.thresh,summarize.replicates,
          #    summary.method,all.missing.thresh,group.missing.thresh,missing.val,samplermindex,rep.max.missing.thresh,summary.na.replacement,featselmethod,pairedanalysis,input.intensity.scale,file="data_preprocess_in.Rda")
          ######
          
         rownames(Xmat)<-names_with_mz_time$Name
         
          #data preprocess classification
          data_matrix<-data_preprocess(Xmat=Xmat,Ymat=Ymat,feature_table_file=feature_table_file,parentoutput_dir=parentoutput_dir,class_labels_file=NA,num_replicates=num_replicates,feat.filt.thresh=NA,summarize.replicates=summarize.replicates,summary.method=summary.method,
                                       all.missing.thresh=all.missing.thresh,group.missing.thresh=group.missing.thresh,
                                       log2transform=log2transform,medcenter=medcenter,znormtransform=znormtransform,,quantile_norm=quantile_norm,lowess_norm=lowess_norm,rangescaling=rangescaling,paretoscaling=paretoscaling,
                                       mstus=mstus,sva_norm=sva_norm,eigenms_norm=eigenms_norm,vsn_norm=vsn_norm,madscaling=madscaling,missing.val=missing.val, rep.max.missing.thresh=rep.max.missing.thresh,
                                       summary.na.replacement=summary.na.replacement,featselmethod=featselmethod,TIC_norm=TIC_norm,normalization.method=normalization.method,
                                       input.intensity.scale=input.intensity.scale,log2.transform.constant=log2.transform.constant,alphabetical.order=alphabetical.order)
          
        #  save(data_matrix,names_with_mz_time,file="data_preprocess_out.Rda")
          
        }else{
          
          
          stop("Orders of feature table and classlabels do not match")
        }
        
        
        
      }else{
        
        #print(diff(match(rnames_xmat,rnames_ymat)))
        stop("Orders of feature table and classlabels do not match")
      }
      
      
      
      if(FALSE){
        data_matrix<-data_preprocess(Xmat,Ymat,
                                     feature_table_file,
                                     parentoutput_dir="C:/Users/kuppal2/Documents/Projects/EGCG_pos//xmsPANDA_preprocess3/",
                                     class_labels_file=NA,num_replicates=1,feat.filt.thresh=NA,summarize.replicates=TRUE,
                                     summary.method="mean", all.missing.thresh=0.5,group.missing.thresh=0.5,
                                     log2transform =FALSE, medcenter=FALSE, znormtransform = FALSE, 
                                     quantile_norm = FALSE, lowess_norm = FALSE, madscaling = FALSE, 
                                     missing.val=0, samplermindex=NA,rep.max.missing.thresh=0.5,summary.na.replacement="zeros")
        
      }
    }else{
      stop("Invalid value for analysismode parameter. Please use regression or classification.")
    }
    
    
  }
  
  if(is.na(names_with_mz_time)==TRUE){
    names_with_mz_time=data_matrix$names_with_mz_time
  }
  #  #save(data_matrix,file="data_matrix.Rda")
  data_matrix_beforescaling<-data_matrix$data_matrix_prescaling
  
  data_matrix_beforescaling<-as.data.frame( data_matrix_beforescaling)
  data_matrix<-data_matrix$data_matrix_afternorm_scaling
  
  
  
  
  
  
  #classlabels<-as.data.frame(classlabels)
  if(dim(classlabels)[2]<2){
    
    stop("The class labels/response matrix should have two columns: SampleID, Class/Response. Please see the example.")
  }
  
  
  
  data_m<-data_matrix[,-c(1:2)]
  classlabels<-classlabels[seq(1,dim(classlabels)[1],num_replicates),]
  
  # #save(classlabels,data_matrix,classlabels_orig,Ymat,file="Stage1/datarose.Rda")
  
  classlabels_raw_boxplots<-classlabels
  

  
  if(dim(classlabels)[2]==2){
    if(length(levels(as.factor(classlabels[,2])))==2){
      if(balance.classes==TRUE){
        
        table_classes<-table(classlabels[,2])
        
        
        suppressWarnings(library(ROSE))
        Ytrain<-classlabels[,2]
        data1=cbind(Ytrain,t(data_matrix[,-c(1:2)]))
        
        ##save(data1,classlabels,data_matrix,file="Stage1/data1.Rda")
        
        #   data_matrix_presim<-data_matrix
        
        data1<-as.data.frame(data1)
        
        colnames(data1)<-c("Ytrain",paste("var",seq(1,ncol(data1)-1),sep=""))
        
        data1$Ytrain<-classlabels[,2]
        
        if(table_classes[1]==table_classes[2])
        {
          
          
          set.seed(balance.classes.seed)
          
          data1[,-c(1)]<-apply(data1[,-c(1)],2,as.numeric)
          new_sample<-aggregate(x=data1[,-c(1)],by=list(as.factor(data1$Ytrain)),mean)
          colnames(new_sample)<-colnames(data1)
          data1<-rbind(data1,new_sample[1,])
          set.seed(balance.classes.seed)
          
          # #save(data1,classlabels,file="Stage1/dataB.Rda")
          
          newData <- ROSE((Ytrain) ~ ., data1, seed = balance.classes.seed,N=nrow(data1)*balance.classes.sizefactor)$data
          
          # newData <- SMOTE(Ytrain ~ ., data=data1, perc.over = 100)
          #*balance.classes.sizefactor,perc.under=200*(balance.classes.sizefactor/(balance.classes.sizefactor/0.5)))
          
        }else{
          if(balance.classes.method=="ROSE"){
            set.seed(balance.classes.seed)
            data1[,-c(1)]<-apply(data1[,-c(1)],2,as.numeric)
            
            
            newData <- ROSE((Ytrain) ~ ., data1, seed = balance.classes.seed,N=nrow(data1)*balance.classes.sizefactor)$data
          }else{
            
            set.seed(balance.classes.seed)
            newData <- SMOTE(Ytrain ~ ., data=data1, perc.over = 100)
            #*balance.classes.sizefactor,perc.under=200*(balance.classes.sizefactor/(balance.classes.sizefactor/0.5)))
            
          }
        }
        newData<-na.omit(newData)
        Xtrain<-newData[,-c(1)]
        Xtrain<-as.matrix(Xtrain)
        Ytrain<-newData[,c(1)]
        
        Ytrain_mat<-cbind((rownames(Xtrain)),(Ytrain))
        Ytrain_mat<-as.data.frame(Ytrain_mat)
        print("new data")
        print(dim(Xtrain))
        print(dim(Ytrain_mat))
        print(table(newData$Ytrain))
        
        data_m<-t(Xtrain)
        data_matrix<-cbind(data_matrix[,c(1:2)],data_m)
        classlabels<-cbind(paste("S",seq(1,nrow(newData)),sep=""),Ytrain)
        classlabels<-as.data.frame(classlabels)
        print(dim(classlabels))
        classlabels_orig<-classlabels
        classlabels_sub<-classlabels[,-c(1)]
        Ymat<-classlabels
        
        ##save(newData,file="Stage1/newData.Rda")
        
      }
    }
  }
  
  
  classlabelsA<-classlabels
  Xmat<-data_matrix
  
  
  #if(dim(classlabels_orig)==TRUE){
  
  
  
  
  
  classlabels_orig<-classlabels_orig[seq(1,dim(classlabels_orig)[1],num_replicates),]
  
  classlabels_response_mat<-as.data.frame(classlabels_response_mat)
  
  classlabels_response_mat<-classlabels_response_mat[seq(1,dim(classlabels_response_mat)[1],num_replicates),]
  
  
  
  class_labels_levels_main<-c("S")
  Ymat<-classlabels
  
  
  rnames1<-as.character(Ymat[,1])
  rnames2<-as.character(classlabels_orig[,1])
  
  sorted_index<-{}
  for(i in 1:length(rnames1)){
    
    
    sorted_index<-c(sorted_index,grep(x=rnames2,pattern=paste("^",rnames1[i],"$",sep="")))
    
  }
  classlabels_orig<-classlabels_orig[sorted_index,]
  
  #write.table(classlabels_response_mat,file="original_classlabelsB.txt",sep="\t",row.names=TRUE)
  classlabelsA<-classlabels
  
  
  if(length(which(duplicated(classlabels)==TRUE))>0){
    rownames(classlabels)<-paste("S",seq(1,dim(classlabels)[1]),sep="")
  }else{
    rownames(classlabels)<-as.character(classlabels[,1])
  }#as.character(classlabels[,1])
  #print(classlabels)
  #print(classlabels[1:10,])
  
  #           ###saveclasslabels,file="classlabels.Rda")
  #          ###saveclasslabels_orig,file="classlabels_orig.Rda")
  #         ###saveclasslabels_response_mat,file="classlabels_response_mat.Rda")
  
  if(pairedanalysis==TRUE){
    
    ###savesubject_inf,file="subjectinf.Rda")
  }
  
  if(analysismode=="classification")
  {	
    
    
    class_labels_levels<-levels(as.factor(classlabels[,2]))
    
   # print("Using the following class labels")			
    #print(class_labels_levels)
    
    class_labels_levels_main<-class_labels_levels
    
    class_labels_levels<-unique(class_labels_levels)
    
    
    bad_rows<-which(class_labels_levels=="")
    if(length(bad_rows)>0){
      class_labels_levels<-class_labels_levels[-bad_rows]
    }
    ordered_labels={}
    num_samps_group<-new("list")
    num_samps_group[[1]]<-0
    groupwiseindex<-new("list")
    groupwiseindex[[1]]<-0
    
    for(c in 1:length(class_labels_levels))
    {
      
      classlabels_index<-which(classlabels[,2]==class_labels_levels[c])
      ordered_labels<-c(ordered_labels,as.character(classlabels[classlabels_index,2]))
      num_samps_group[[c]]<-length(classlabels_index)
      groupwiseindex[[c]]<-classlabels_index
    }
    
    Ymatorig<-classlabels
    
    #debugclasslabels
    #save(classlabels,class_labels_levels,num_samps_group,Ymatorig,data_matrix,data_m_fc_withfeats,data_m,file="classlabels_1.Rda")
    ####saveclass_labels_levels,file="class_labels_levels.Rda")
    
   # print("HERE1")

    classlabels_dataframe<-classlabels
      
      class_label_alphabets<-class_labels_levels
      classlabels<-{}
      
      if(length(class_labels_levels)==2){
        #num_samps_group[[1]]=length(which(ordered_labels==class_labels_levels[1]))
        #num_samps_group[[2]]=length(which(ordered_labels==class_labels_levels[2]))
        class_label_A<-class_labels_levels[[1]]
        class_label_B<-class_labels_levels[[2]]
        #classlabels<-c(rep("ClassA",num_samps_group[[1]]),rep("ClassB",num_samps_group[[2]]))
        classlabels<-c(rep(class_label_A,num_samps_group[[1]]),rep(class_label_B,num_samps_group[[2]]))
      }else{
        if(length(class_labels_levels)==3){
          
          class_label_A<-class_labels_levels[[1]]
          class_label_B<-class_labels_levels[[2]]
          class_label_C<-class_labels_levels[[3]]
          classlabels<-c(rep(class_label_A,num_samps_group[[1]]),rep(class_label_B,num_samps_group[[2]]),rep(class_label_C,num_samps_group[[3]]))
        }else{
          
          for(c in 1:length(class_labels_levels)){
            
            num_samps_group_cur=length(which(Ymatorig[,2]==class_labels_levels[c]))
            
            classlabels<-c(classlabels,rep(paste(class_labels_levels[c],sep=""),num_samps_group_cur))
            #,rep("ClassB",num_samps_group[[2]]),rep("ClassC",num_samps_group[[3]]))
            
          }
          
        }
      }
      #   print("Class mapping:")
      # print(cbind(class_labels_levels,classlabels))
      
      
      
      classlabels<-classlabels_dataframe[,2]
      
      
  
    classlabels_2=classlabels
    
    
    #save(classlabels_2,class_labels_levels,Ymatorig,data_matrix,data_m_fc_withfeats,data_m,file="classlabels_2.Rda")
    
    
    ####################################################################################
    #print(head(data_m))
    
    snames<-colnames(data_m)
    
    Ymat<-as.data.frame(classlabels)
    m1<-match(snames,Ymat[,1])
    #Ymat<-Ymat[m1,]
    
    data_temp<-data_matrix_beforescaling[,-c(1:2)]
    
    
   rnames<-paste("mzid_",seq(1,nrow(data_matrix)),sep="")
   rownames(data_m)=rnames
    
    mzid_mzrt<-data_matrix[,c(1:2)]
    colnames(mzid_mzrt)<-c("mz","time")
    rownames(mzid_mzrt)=rnames
    write.table(mzid_mzrt, file="Stage1/mzid_mzrt.txt",sep="\t",row.names=TRUE)
    
    
    cl<-makeCluster(num_nodes)
    
    
    
    mean_overall<-apply(data_temp,1,do_mean)
    
    #clusterExport(cl,"do_mean")
    #mean_overall<-parApply(cl,data_temp,1,do_mean)
    
    #stopCluster(cl)
    
    #mean_overall<-unlist(mean_overall)
   # print("mean overall")
    #print(summary(mean_overall))
    bad_feat<-which(mean_overall==0)
    
    if(length(bad_feat)>0){
      
      data_matrix_beforescaling<-data_matrix_beforescaling[-bad_feat,]
      data_m<-data_m[-bad_feat,]
      data_matrix<-data_matrix[-bad_feat,]
      
    } 
    
    
    #Step 5) RSD/CV calculation
    
    
  }else{
    
    classlabels<-(classlabels[,-c(1)])
  }
  
  #	print("######classlabels#########")
  #print(classlabels)
  
  
  class_labels_levels_new<-levels(classlabels)
  
  if(analysismode=="classification"){
    test_classlabels<-cbind(class_labels_levels_main,class_labels_levels_new)
  }
  
  if(featselmethod=="ttest" | featselmethod=="wilcox"){
    
    if(length(class_labels_levels)>2){
      
      print("#######################")
      print(paste("Warning: More than two classes detected. Invalid feature selection option. Skipping the feature selection for option ",featselmethod,sep=""))
      
      print("#######################")
      
      return("More than two classes detected. Invalid feature selection option.")
      
    }
    
  }
  
  
  #print("here 2")
  ######################################################################################
  
  #Step 6) Log2 mean fold change criteria from 0 to 1 with step of 0.1
  feat_eval<-{}
  feat_sigfdrthresh<-{}
  feat_sigfdrthresh_cv<-{}
  feat_sigfdrthresh_permut<-{}
  
  permut_acc<-{}
  feat_sigfdrthresh<-rep(0,length(log2.fold.change.thresh_list))
  feat_sigfdrthresh_cv<-rep(NA,length(log2.fold.change.thresh_list))
  
  feat_sigfdrthresh_permut<-rep(NA,length(log2.fold.change.thresh_list))
  res_score_vec<-rep(0,length(log2.fold.change.thresh_list))
  #feat_eval<-seq(0,1,0.1)
  
  if(analysismode=="classification"){
    
    best_cv_res<-(-1)*10^30
  }else{
    best_cv_res<-(1)*10^30
    
  }
  
  
  best_feats<-{}
  
  goodfeats<-{}
  mwan_fdr<-{}
  targetedan_fdr<-{}
  best_limma_res<-{}
  best_acc<-{}
  termA<-{}
  
  fheader="transformed_log2fc_threshold_"
  
  
  X<-t(data_m)
  
  X<-replace(as.matrix(X),which(is.na(X)==TRUE),0)
  
  
  
  
  
  # rm(pcaMethods)
  #try(detach("package:pcaMethods",unload=TRUE),silent=TRUE)
  library(mixOmics)
  
  if(featselmethod=="lmreg" || featselmethod=="lmregrobust" || featselmethod=="logitreg" || featselmethod=="logitregrobust"){
    
    if(length(class_labels_levels)>2){
      
      stop(paste(featselmethod, " feature selection option is only available for 2 class comparisons."),sep="")
      
      
      
    }
  }
  
  
  
  if(sample.col.opt=="default"){
    
    col_vec<-c("#CC0000","#AAC000","blue","mediumpurple4","mediumpurple1","blueviolet","cornflowerblue","cyan4","skyblue",
               "darkgreen", "seagreen1", "green","yellow","orange","pink", "coral1", "palevioletred2",
               "red","saddlebrown","brown","brown3","white","darkgray","aliceblue",
               "aquamarine","aquamarine3","bisque","burlywood1","lavender","khaki3","black")
    
  }else{ 
    if(sample.col.opt=="topo"){
      #col_vec<-topo.colors(256) #length(class_labels_levels)) 
      
      #col_vec<-col_vec[seq(1,length(col_vec),)]
      
      col_vec <- topo.colors(length(class_labels_levels), alpha=alphacol)
    }else{
      if(sample.col.opt=="heat"){
        #col_vec<-heat.colors(256) #length(class_labels_levels))
        
        col_vec <- heat.colors(length(class_labels_levels), alpha=alphacol)
      }else{
        if(sample.col.opt=="rainbow"){
          #col_vec<-heat.colors(256) #length(class_labels_levels))
          col_vec<-rainbow(length(class_labels_levels), start = 0, end = alphacol)
          
          #col_vec <- heat.colors(length(class_labels_levels), alpha=alphacol)
        }else{
          
          if(sample.col.opt=="terrain"){
            #col_vec<-heat.colors(256) #length(class_labels_levels))
            #col_vec<-rainbow(length(class_labels_levels), start = 0, end = alphacol)
            
            col_vec <- cm.colors(length(class_labels_levels), alpha=alphacol)
          }else{
            
            if(sample.col.opt=="colorblind"){
              #col_vec <-c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")
              # col_vec <- c("#0072B2", "#E69F00", "#009E73", "gold1", "#56B4E9", "#D55E00", "#CC79A7","black")
              
              if(length(class_labels_levels)<9){
                
                col_vec <- c("#0072B2", "#E69F00", "#009E73", "#56B4E9", "#D55E00", "#CC79A7", "#E64B35FF", "grey57")
                
              }else{
                
                #col_vec<-colorRampPalette(brewer.pal(10, "RdBu"))(length(class_labels_levels))
                col_vec<-c("#0072B2", "#E69F00", "#009E73", "#56B4E9", "#D55E00", "#CC79A7","#E64B35B2", "#4DBBD5B2","#00A087B2","#3C5488B2","#F39B7FB2","#8491B4B2","#91D1C2B2","#DC0000B2","#7E6148B2",
                           "#374E55B2","#DF8F44B2","#00A1D5B2","#B24745B2","#79AF97B2","#6A6599B2","#80796BB2","#0073C2B2","#EFC000B2", "#868686B2","#CD534CB2","#7AA6DCB2","#003C67B2","grey57")
                
              }
              
              
            }else{
              
              check_brewer<-grep(pattern="brewer",x=sample.col.opt)
              
              if(length(check_brewer)>0){
                
                sample.col.opt_temp=gsub(x=sample.col.opt,pattern="brewer.",replacement="")
                col_vec <- colorRampPalette(brewer.pal(10, sample.col.opt_temp))(length(class_labels_levels))
                
              }else{
                
                if(sample.col.opt=="journal"){
                  
                  col_vec<-c("#0072B2", "#E69F00", "#009E73", "#56B4E9", "#D55E00", "#CC79A7","#E64B35FF","#3C5488FF","#F39B7FFF",
                             "#8491B4FF","#91D1C2FF","#DC0000FF","#B09C85FF","#5F559BFF",
                             "#808180FF","#20854EFF","#FFDC91FF","#B24745FF",
                             
                             "#374E55FF","#8F7700FF","#5050FFFF","#6BD76BFF",
                             "#E64B3519","#4DBBD519","#631879E5","grey75")
                  if(length(class_labels_levels)<8){
                    col_vec<-c("#0072B2", "#E69F00", "#009E73", "#56B4E9", "#D55E00", "#CC79A7","grey75")
                    
                    #col_vec2<-brewer.pal(n = 8, name = "Dark2")
                    
                  }else{
                    if(length(class_labels_levels)<=28){
                      # col_vec<-c("#0072B2", "#E69F00", "#009E73", "#56B4E9", "#D55E00", "#CC79A7", "grey75","#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D", "#666666","#1B9E77", "#7570B3", "#E7298A", "#A6761D", "#666666", "#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D", "#666666")
                      
                      col_vec<-c("#0072B2", "#E69F00", "#009E73", "#56B4E9", "#D55E00", "#CC79A7","#E64B35FF","#3C5488FF","#F39B7FFF",
                                 "#8491B4FF","#91D1C2FF","#DC0000FF","#B09C85FF","#5F559BFF",
                                 "#808180FF","#20854EFF","#FFDC91FF","#B24745FF",
                                 
                                 "#374E55FF","#8F7700FF","#5050FFFF","#6BD76BFF", "#8BD76BFF",
                                 "#E64B3519","#9DBBD0FF","#631879E5","#666666","grey75")
                      
                    }else{
                      
                      
                      
                      
                      colfunc <-colorRampPalette(c("#0072B2", "#E69F00", "#009E73", "#56B4E9", "#D55E00", "#CC79A7","grey75"));col_vec<-colfunc(length(class_labels_levels))
                      
                      col_vec<-col_vec[sample(col_vec)]
                      
                      
                    }
                  }
                  
                }else{
                  #colfunc <-colorRampPalette(sample.col.opt);col_vec<-colfunc(length(class_labels_levels))
                  if(length(sample.col.opt)==1){
                    col_vec <-rep(sample.col.opt,length(class_labels_levels))
                  }else{
                    
                    colfunc <-colorRampPalette(sample.col.opt);col_vec<-colfunc(length(class_labels_levels))
                    
                    col_vec<-col_vec[sample(col_vec)]
                    
                  }
                }
                
              }
              
            }
          }
          
          
        }
        
      }
      
    }	
  }
  #pca_col_vec<-col_vec
  
  pca_col_vec<-c("mediumpurple4","mediumpurple1","blueviolet","darkblue","blue","cornflowerblue","cyan4","skyblue",
                 "darkgreen", "seagreen1", "green","yellow","orange","pink", "coral1", "palevioletred2",
                 "red","saddlebrown","brown","brown3","white","darkgray","aliceblue",
                 "aquamarine","aquamarine3","bisque","burlywood1","lavender","khaki3","black")
  
  
  if(is.na(individualsampleplot.col.opt)==TRUE){
    
    individualsampleplot.col.opt=col_vec
  }
  
  
  
  
  #cl<-makeCluster(num_nodes)
  #feat_sds<-parApply(cl,data_m,1,sd)
  feat_sds<-apply(data_m,1,function(x){sd(x,na.rm=TRUE)})
  
  #stopCluster(cl)
  
  bad_sd_ind<-c(which(feat_sds==0),which(is.na(feat_sds)==TRUE))
  
  bad_sd_ind<-unique(bad_sd_ind)
  
  if(length(bad_sd_ind)>0){
    
    data_matrix<-data_matrix[-c(bad_sd_ind),]
    
    data_m<-data_m[-c(bad_sd_ind),]
    
    data_matrix_beforescaling<-data_matrix_beforescaling[-c(bad_sd_ind),]
  }
  
  data_temp<-data_matrix_beforescaling[,-c(1:2)]
  
  
  
  
  
  #cl<-makeCluster(num_nodes)
  
  #clusterExport(cl,"do_rsd")
  #feat_rsds<-parApply(cl,data_temp,1,do_rsd)
  feat_rsds<-apply(data_temp,1,do_rsd)
  #stopCluster(cl)
  #  #save(feat_rsds,data_temp,data_matrix_beforescaling,data_m,file="rsds.Rda")
  sum_rsd<-summary(feat_rsds,na.rm=TRUE)
  max_rsd<-max(feat_rsds,na.rm=TRUE)
  max_rsd<-round(max_rsd,2)
  
  
  print("Summary of RSD across all features:")
  print(sum_rsd)
  
  if(log2.fold.change.thresh_list[length(log2.fold.change.thresh_list)]>max_rsd){
    stop(paste("The maximum relative standard deviation threshold in rsd.filt.list should be below ",max_rsd,sep=""))
  }
  
  classlabels_parent<-classlabels
  classlabels_sub_parent<-classlabels_sub
  classlabels_orig_parent<-classlabels_orig
  
  #write.table(classlabels_orig,file="classlabels.txt",sep="\t",row.names=FALSE)
  classlabels_response_mat_parent<-classlabels_response_mat
  
  
  
  
  
  parent_data_m<-round(data_m,5)
  
  
  res_score<-0
  #best_cv_res<-0
  best_feats<-{}
  
  best_acc<-0
  best_limma_res<-{}
  best_logfc_ind<-1
  
  
  
  output_dir1<-paste(parentoutput_dir,"/Stage2/",sep="")
  dir.create(output_dir1,showWarnings=FALSE)
  
  setwd(output_dir1)
  
  classlabels_sub_parent<-classlabels_sub
  classlabels_orig_parent<-classlabels_orig
  
  #write.table(classlabels_orig,file="classlabels.txt",sep="\t",row.names=FALSE)
  classlabels_response_mat_parent<-classlabels_response_mat
  
  rocfeatlist<-rocfeatlist+1
  
  if(pairedanalysis==TRUE){
    #print(subject_inf)
    write.table(subject_inf,file="subject_inf.txt",sep="\t")
    paireddesign=subject_inf
  }else{
    paireddesign=NA
    
    
  }
  #write.table(classlabels_orig,file="classlabels_orig.txt",sep="\t")
  #write.table(classlabels,file="classlabels.txt",sep="\t")
  #write.table(classlabels_response_mat,file="classlabels_response_mat.txt",sep="\t")
  
  if(is.na(max_varsel)==TRUE){
    
    max_varsel=dim(data_m)[1]
  }
  
  for(lf in 1:length(log2.fold.change.thresh_list))
  {
    
    allmetabs_res<-{}
    classlabels_response_mat<-classlabels_response_mat_parent
    classlabels_sub<-classlabels_sub_parent
    classlabels_orig<-classlabels_orig_parent
    
    setwd(parentoutput_dir)
    log2.fold.change.thresh=log2.fold.change.thresh_list[lf]
    
    
    output_dir1<-paste(parentoutput_dir,"/Stage2/",sep="")
    dir.create(output_dir1,showWarnings=FALSE)
    
    setwd(output_dir1)
    
    
    
    if(logistic_reg==TRUE){
      
      if(robust.estimate==FALSE){ output_dir<-paste(output_dir1,"logitreg","signalthresh",group.missing.thresh,"RSD",log2.fold.change.thresh,"/",sep="")
      }else{
        
        if(robust.estimate==TRUE){ output_dir<-paste(output_dir1,"logitregrobust","signalthresh",group.missing.thresh,"RSD",log2.fold.change.thresh,"/",sep="")
        }
      }
    }else{
      
      
      if(poisson_reg==TRUE){
        
        if(robust.estimate==FALSE){ output_dir<-paste(output_dir1,"poissonreg","signalthresh",group.missing.thresh,"RSD",log2.fold.change.thresh,"/",sep="")
        }else{
          if(robust.estimate==TRUE){
            output_dir<-paste(output_dir1,"poissonregrobust","signalthresh",group.missing.thresh,"RSD",log2.fold.change.thresh,"/",sep="")
            
          }
          
        }
      }else{
        
        if(featselmethod=="lmreg"){
          
          if(robust.estimate==TRUE){
            output_dir<-paste(output_dir1,"lmregrobust","signalthresh",group.missing.thresh,"RSD",log2.fold.change.thresh,"/",sep="")
            
          }else{
            
            output_dir<-paste(output_dir1,"lmreg","signalthresh",group.missing.thresh,"RSD",log2.fold.change.thresh,"/",sep="")
            
            
          }
          
        }else{
          
          
          
          output_dir<-paste(output_dir1,parentfeatselmethod,"signalthresh",group.missing.thresh,"RSD",log2.fold.change.thresh,"/",sep="")
        }
        
      }
    }
    dir.create(output_dir,showWarnings=FALSE)
    
    setwd(output_dir)
    
    dir.create("Figures",showWarnings = FALSE)
    
    dir.create("Tables",showWarnings = FALSE)
    
    
    
    data_m<-parent_data_m
    
    #print("dim of data_m")
    #print(dim(data_m))
    
    pdf_fname<-paste("Figures/Results_RSD",log2.fold.change.thresh,".pdf",sep="")
    
    #zip_fname<-paste("Results_RSD",log2.fold.change.thresh,".zip",sep="")
    
    if(output.device.type=="pdf"){
      pdf(pdf_fname,width=10,height=10)
    }
    
    if(analysismode=="classification" | analysismode=="regression"){
      
      print(paste("Performing RSD filtering using ",log2.fold.change.thresh, " as threshold",sep=""))
      if(log2.fold.change.thresh>=0){
        
        if(log2.fold.change.thresh==0){
          log2.fold.change.thresh=0.001
        }
        
        
        #good_metabs<-which(abs(mean_groups)>log2.fold.change.thresh)
        abs_feat_rsds<-abs(feat_rsds)
        
        good_metabs<-which(abs_feat_rsds>log2.fold.change.thresh)
        
        #print("length of good_metabs")
        #print(good_metabs)
        
      }else{
        good_metabs<-seq(1,dim(data_m)[1])
      }
      
      if(length(good_metabs)>0){
        
        data_m_fc<-data_m[good_metabs,]
        
        data_m_fc_withfeats<-data_matrix[good_metabs,c(1:2)]
        
        data_matrix_beforescaling_rsd<-data_matrix_beforescaling[good_metabs,]
        
      }else{
        #data_m_fc<-data_m
        #data_m_fc_withfeats<-data_matrix[,c(1:2)]
        
        stop(paste("Please decrease the maximum relative standard deviation (rsd.filt.thresh) threshold to ",max_rsd,sep=""))
        
      }
    }else{
      
      data_m_fc<-data_m
      data_m_fc_withfeats<-data_matrix[,c(1:2)]
    }
   # save(data_m_fc_withfeats,data_m_fc,data_m,data_matrix,file="datadebug.Rda")
    
    
    ylab_text_raw<-ylab_text
    
    if(log2transform==TRUE || input.intensity.scale=="log2"){
      
      if(znormtransform==TRUE){
        ylab_text_2="scale normalized"
      }else{
        if(quantile_norm==TRUE){
          
          ylab_text_2="quantile normalized"
        }else{
          ylab_text_2=""
        }
      }
      ylab_text=paste("log2 ",ylab_text," ",ylab_text_2,sep="")
    }else{
      if(znormtransform==TRUE){
        ylab_text_2="scale normalized"
      }else{
        if(quantile_norm==TRUE){
          
          ylab_text_2="quantile normalized"
        }else{
          ylab_text_2=""
        }
      }
      ylab_text=paste("Raw ",ylab_text," ",ylab_text_2,sep="") #paste("Raw intensity ",ylab_text_2,sep="")
    }
    #ylab_text=paste("Abundance",sep="")
    
    
    if(is.na(names_with_mz_time)==FALSE){
      data_m_fc_with_names<-merge(names_with_mz_time,data_m_fc_withfeats,by=c("mz","time"))
      data_m_fc_with_names<-data_m_fc_with_names[match(data_m_fc_withfeats$mz,data_m_fc_with_names$mz),]
      #save(names_with_mz_time,goodfeats,goodfeats_with_names,file="goodfeats_with_names.Rda")
      
      # goodfeats_name<-goodfeats_with_names$Name
      #}
    }
    
  #  save(data_m_fc_withfeats,data_matrix,data_m,data_m_fc,data_m_fc_with_names,names_with_mz_time,file="debugnames.Rda")
    
    
    
    if(dim(data_m_fc)[2]>50){
      
      if(output.device.type!="pdf"){
        
        temp_filename_1<-"Figures/SampleIntensityDistribution.png"
        
        png(temp_filename_1,width=plots.width,height=plots.height,res=plots.res,type=plots.type,units="in")
        
      }
      
      size_num<-min(100,dim(data_m_fc)[2])
      
      par(mfrow=c(1,1),family="sans",cex=cex.plots)
      samp_index<-sample(x=1:dim(data_m_fc)[2],size=size_num)
      #  try(boxplot(data_m_fc[,samp_index],main="Intensity distribution across samples after preprocessing",xlab="Samples",ylab=ylab_text,col=boxplot.col.opt),silent=TRUE)
      
      #samp_dist_col<-get_boxplot_colors(boxplot.col.opt,class_labels_levels=c(1))
      
      boxplot(data_m_fc[,samp_index],main="Intensity distribution across samples after preprocessing",xlab="Samples",ylab=ylab_text,col="white")
      
      if(output.device.type!="pdf"){
        
        try(dev.off(),silent=TRUE)
      }
      
    }else{
      
      if(output.device.type!="pdf"){
        
        temp_filename_1<-"Figures/SampleIntensityDistribution.png"
        
        png(temp_filename_1,width=plots.width,height=plots.height,res=plots.res,type=plots.type,units="in")
      }
      par(mfrow=c(1,1),family="sans",cex=cex.plots)
      try(boxplot(data_m_fc,main="Intensity distribution across samples after preprocessing",xlab="Samples",ylab=ylab_text,col="white"),silent=TRUE)
      if(output.device.type!="pdf"){
        
        try(dev.off(),silent=TRUE)
      }
      
    }
    
    
    if(output.device.type!="pdf"){
      
      temp_filename_1<-paste("Figures/OutlierDetection",outlier.method,".png",sep="")
      
      png(temp_filename_1,width=plots.width,height=plots.height,res=plots.res,type=plots.type,units="in")
    }
    par(mfrow=c(1,1),family="sans",cex=cex.plots)
    
    ##save(data_matrix,file="dm1.Rda")
    outlier_detect(data_matrix=data_matrix,ncomp=2,column.rm.index=c(1,2),outlier.method=outlier.method)
    
    print("done outlier")
    
    if(output.device.type!="pdf"){
      
      try(dev.off(),silent=TRUE)
    }
    
    
    
    
    data_m_fc_withfeats<-cbind(data_m_fc_withfeats,data_m_fc)
    
    allmetabs_res_withnames<-{}
    
    feat_eval[lf]<-0
    res_score_vec[lf]<-0
    #feat_sigfdrthresh_cv[lf]<-0
    
    filename<-paste(fheader,log2.fold.change.thresh,".txt",sep="")
    #write.table(data_m_fc_withfeats, file=filename,sep="\t",row.names=FALSE)
    
    
    if(length(data_m_fc)>=dim(parent_data_m)[2])
    {
      
      
      if(dim(data_m_fc)[1]>0){
        
        
        
        feat_eval[lf]<-dim(data_m_fc)[1]
        
        # col_vec<-c("#CC0000","#AAC000","blue","mediumpurple4","mediumpurple1","blueviolet","darkblue","blue","cornflowerblue","cyan4","skyblue",
        #"darkgreen", "seagreen1", "green","yellow","orange","pink", "coral1", "palevioletred2",
        #"red","saddlebrown","brown","brown3","white","darkgray","aliceblue",
        #"aquamarine","aquamarine3","bisque","burlywood1","lavender","khaki3","black")
        
        if(analysismode=="classification")
        {
          
          sampleclass<-{}
          patientcolors<-{}
          #
          classlabels<-as.data.frame(classlabels)
          #print(classlabels)
          
          f<-factor(classlabels[,1])
          
          for(c in 1:length(class_labels_levels)){
            
            num_samps_group_cur=length(which(ordered_labels==class_labels_levels[c]))
            
            #classlabels<-c(classlabels,rep(paste("Class",class_label_alphabets,sep=""),num_samps_group_cur))
            #,rep("ClassB",num_samps_group[[2]]),rep("ClassC",num_samps_group[[3]]))
            sampleclass<-c(sampleclass,rep(paste("Class",class_label_alphabets[c],sep=""),num_samps_group_cur))
            #sampleclass<-classlabels[,1] #c(sampleclass,rep(paste("Class",class_labels_levels[c],sep=""),num_samps_group_cur))
            
            patientcolors <-c(patientcolors,rep(col_vec[c],num_samps_group_cur))
          }
          
          
          
          
          library(pcaMethods)
          
          #p1<-pcaMethods::pca(data_m_fc,method="rnipals",center=TRUE,scale="uv",cv="q2",nPcs=3)
          
          tempX<-t(data_m_fc)
          
          
          
          #  p1<-pcaMethods::pca(tempX,method="rnipals",center=TRUE,scale="uv",cv="q2",nPcs=10)
          
          
          
          if(output.device.type!="pdf"){
            
            temp_filename_2<-"Figures/PCAdiagnostics_allfeats.png"
            
            # png(temp_filename_2,width=plots.width,height=plots.height,res=plots.res,type=plots.type,units="in")
          }
          
          
          
          if(output.device.type!="pdf"){
            
            # dev.off()
          }
          
          try(detach("package:pcaMethods",unload=TRUE),silent=TRUE)
          
          
          
          if(dim(classlabels)[2]>2){
            classgroup<-paste(classlabels[,1],":",classlabels[,2],sep="") #classlabels[,1]:classlabels[,2]
          }else{
            
            classgroup<-classlabels
          }
          
          classlabels_orig<-classlabels_orig_parent
          
          if(pairedanalysis==TRUE){
            
            #classlabels_orig<-classlabels_orig[,-c(2)]
            
            
          }else{
            
            if(featselmethod=="lmreg" || featselmethod=="logitreg" ||  featselmethod=="poissonreg"){
              classlabels_orig<-classlabels_orig[,c(1:2)]
              classlabels_orig<-as.data.frame(classlabels_orig)
            }
          }
          
          if(analysismode=="classification"){
            if(dim(classlabels_orig)[2]==2){
              if(alphabetical.order==FALSE){
                classlabels_orig[,2] <- factor(classlabels_orig[,2], levels=unique(classlabels_orig[,2]))
              }
            }
            if(dim(classlabels_orig)[2]==3){
              
              if(pairedanalysis==TRUE){
                if(alphabetical.order==FALSE){
                  classlabels_orig[,3] <- factor(classlabels_orig[,3], levels=unique(classlabels_orig[,3]))
                }
              }else{
                if(alphabetical.order==FALSE){
                  classlabels_orig[,2] <- factor(classlabels_orig[,2], levels=unique(classlabels_orig[,2]))
                  classlabels_orig[,3] <- factor(classlabels_orig[,3], levels=unique(classlabels_orig[,3]))
                }
                
              }
              
            }else{
              
              if(dim(classlabels_orig)[2]==4){
                
                if(pairedanalysis==TRUE){
                  if(alphabetical.order==FALSE){
                    classlabels_orig[,3] <- factor(classlabels_orig[,3], levels=unique(classlabels_orig[,3]))
                    
                    classlabels_orig[,4] <- factor(classlabels_orig[,4], levels=unique(classlabels_orig[,4]))
                  }
                }
                
              }
              
            }
            
          }
          
          
          if(output.device.type!="pdf"){
            
            temp_filename_1<-"Figures/PCAplots_allfeats.pdf"
            
            #png(temp_filename_1,width=plots.width,height=plots.height,res=plots.res,type=plots.type,units="in")
            
            pdf(temp_filename_1,width=plots.width,height=plots.height)
          }
          
          
          
          plot(0:10, type = "n", xaxt="n", yaxt="n", bty="n", xlab = "", ylab = "")
          
          
          text(5, 8, "PCA using all features left after pre-processing")
          text(5, 7, "The figures include: ")
          text(5, 6, "a. pairwise PC score plots ")
          text(5, 5, "b. scores for individual samples on each PC")
          text(5, 4, "c. Lineplots using PC scores for data with repeated measurements")
          
          ###savelist=ls(),file="pcaplotsall.Rda")
          
          
        #  save(data_m_fc_withfeats,classlabels_orig,sample.col.opt,col_vec,pairedanalysis,pca.cex.val,legendlocation,pca.ellipse,ellipse.conf.level,paireddesign,
         #      lineplot.col.opt,lineplot.lty.option,timeseries.lineplots,pcacenter,pcascale,alphabetical.order,
          #     analysistype,lme.modeltype,file="pcaplotsall.Rda")
          
          
          if(length(which(duplicated(data_m_fc_with_names$Name)==TRUE))>0){
            
            print("Duplicate features detected")
            print("Removing duplicate entries for the following features:")
           # print(data_m_fc_with_names$Name[which(duplicated(data_m_fc_with_names$Name)==TRUE)])
            
            data_m_fc_withfeats<-data_m_fc_withfeats[-which(duplicated(data_m_fc_with_names$Name)==TRUE),]
            
            data_m_fc<-data_m_fc[-which(duplicated(data_m_fc_with_names$Name)==TRUE),]
            
            data_matrix<-data_matrix[-which(duplicated(data_m_fc_with_names$Name)==TRUE),]
            
            data_m<-data_m[-which(duplicated(data_m_fc_with_names$Name)==TRUE),]
            
            data_m_fc_with_names<-data_m_fc_with_names[-which(duplicated(data_m_fc_with_names$Name)==TRUE),]
            
            #parent_data_m<-parent_data_m[-which(duplicated(data_m_fc_with_names$Name)==TRUE),]
            
            
           
          }
          
          rownames(data_m_fc_withfeats)<-data_m_fc_with_names$Name
          
          classlabels_orig_pca<-classlabels_orig
          c1=try(get_pcascoredistplots(X=data_m_fc_withfeats,Y=classlabels_orig,feature_table_file=NA,parentoutput_dir=getwd(),class_labels_file=NA,sample.col.opt=sample.col.opt,
                                       plots.width=2000,plots.height=2000,plots.res=300, alphacol=0.3,col_vec=col_vec,pairedanalysis=pairedanalysis,pca.cex.val=pca.cex.val,legendlocation=legendlocation,
                                       pca.ellipse=pca.ellipse,ellipse.conf.level=ellipse.conf.level,
                                       filename="all",paireddesign=paireddesign,lineplot.col.opt=lineplot.col.opt,
                                       lineplot.lty.option=lineplot.lty.option,timeseries.lineplots=timeseries.lineplots,
                                       pcacenter=pcacenter,pcascale=pcascale,alphabetical.order=alphabetical.order,
                                       study.design=analysistype,lme.modeltype=lme.modeltype),silent=TRUE)
          
          
          
          
          if(output.device.type!="pdf"){
            
            try(dev.off(),silent=TRUE)
          }
          
          
          classlabels_orig<-classlabels_orig_parent
        }else{
          #regression
          tempgroup<-rep("A",dim(data_m_fc)[2]) #cbind(classlabels_orig[,1],
          col_vec1<-rep("black",dim(data_m_fc)[2])
          class_labels_levels_main1<-c("A")
          
          analysistype="regression"
          
          
          
          #   get_pca(X=data_m_fc,samplelabels=tempgroup,legendlocation=legendlocation,filename="all",
          #          ncomp=3,pcacenter=pcacenter,pcascale=pcascale,legendcex=0.5,outloc=getwd(),col_vec=col_vec1,
          #         sample.col.opt=sample.col.opt,alphacol=0.3,class_levels=NA,pca.cex.val=pca.cex.val,pca.ellipse=FALSE,
          #        paireddesign=paireddesign,alphabetical.order=alphabetical.order,pairedanalysis=pairedanalysis,classlabels_orig=classlabels_orig,analysistype=analysistype) #,silent=TRUE)
          if(output.device.type!="pdf"){
            
            temp_filename_1<-"Figures/PCAplots_allfeats.pdf"
            
            #png(temp_filename_1,width=plots.width,height=plots.height,res=plots.res,type=plots.type,units="in")
            
            pdf(temp_filename_1,width=plots.width,height=plots.height)
          }
          
          
          
          plot(0:10, type = "n", xaxt="n", yaxt="n", bty="n", xlab = "", ylab = "")
          
          
          text(5, 8, "PCA using all features left after pre-processing")
          text(5, 7, "The figures include: ")
          text(5, 6, "a. pairwise PC score plots ")
          text(5, 5, "b. scores for individual samples on each PC")
          text(5, 4, "c. Lineplots using PC scores for data with repeated measurements")
          
          ###savelist=ls(),file="pcaplotsall.Rda")
          
          if(length(which(duplicated(data_m_fc_with_names$Name)==TRUE))>0){
            
            data_m_fc_withfeats<-data_m_fc_withfeats[-which(duplicated(data_m_fc_with_names$Name)==TRUE),]
            
            data_m_fc<-data_m_fc[-which(duplicated(data_m_fc_with_names$Name)==TRUE),]
            
            data_matrix<-data_matrix[-which(duplicated(data_m_fc_with_names$Name)==TRUE),]
            
            data_m<-data_m[-which(duplicated(data_m_fc_with_names$Name)==TRUE),]
            
            data_m_fc_with_names<-data_m_fc_with_names[-which(duplicated(data_m_fc_with_names$Name)==TRUE),]
            
           # parent_data_m<-parent_data_m[-which(duplicated(data_m_fc_with_names$Name)==TRUE),]
            
            
            print("Duplicate features detected")
            print("Removing duplicate entries for the following features:")
            print(data_m_fc_with_names$Name[which(duplicated(data_m_fc_with_names$Name)==TRUE)])
            print(nrow(data_m_fc_withfeats))
            print(nrow(data_m_fc))
            print(nrow(data_m))
            print(nrow(data_matrix))
          }
          
          
          ###save(data_m_fc_withfeats,classlabels_orig,sample.col.opt,col_vec,pairedanalysis,pca.cex.val,legendlocation,pca.ellipse,ellipse.conf.level,paireddesign,lineplot.col.opt,lineplot.lty.option,timeseries.lineplots,pcacenter,pcascale,file="pcaplotsall.Rda")
          rownames(data_m_fc_withfeats)<-data_m_fc_with_names$Name
          c1=try(get_pcascoredistplots(X=data_m_fc_withfeats,Y=classlabels_orig,feature_table_file=NA,parentoutput_dir=getwd(),class_labels_file=NA,
                                       sample.col.opt=sample.col.opt,
                                       plots.width=2000,plots.height=2000,plots.res=300, alphacol=0.3,col_vec=col_vec,pairedanalysis=pairedanalysis,pca.cex.val=pca.cex.val,legendlocation=legendlocation,
                                       pca.ellipse=pca.ellipse,ellipse.conf.level=ellipse.conf.level,filename="all",
                                       paireddesign=paireddesign,lineplot.col.opt=lineplot.col.opt,lineplot.lty.option=lineplot.lty.option,
                                       timeseries.lineplots=timeseries.lineplots,pcacenter=pcacenter,pcascale=pcascale,alphabetical.order=alphabetical.order,
                                       study.design=analysistype,lme.modeltype=lme.modeltype),silent=TRUE)
          
          if(output.device.type!="pdf"){
            
            try(dev.off(),silent=TRUE)
          }
          
        }
        
        
        if(featselmethod=="pamr"){
          
          #print("HERE")
          #savedata_m_fc,classlabels,file="pamdebug.Rda")
          
          if(is.na(fdrthresh)==FALSE){
            if(fdrthresh>0.5){
              
              pamrthresh=pvalue.thresh
              
            }else{
              pamrthresh=fdrthresh
              
            }
          }else{
            
            pamrthresh=pvalue.thresh
            
          }
          pamr.res<-do_pamr(X=data_m_fc,Y=classlabels,fdrthresh=pamrthresh,nperms=1000,pamr.threshold.select.max=pamr.threshold.select.max,kfold=kfold)
          
          ###save(pamr.res,file="pamr.res.Rda")
          goodip<-pamr.res$feature.list
          
          if(length(goodip)<1){
            goodip=NA
          }
          pamr.threshold_value<-pamr.res$threshold_value
          
          feature_rowindex<-seq(1,nrow(data_m_fc))
          
          discore<-rep(0,nrow(data_m_fc))
          
          discore_all<-pamr.res$max.discore.allfeats
          
          if(is.na(goodip)==FALSE){
            discore[goodip]<-pamr.res$max.discore.sigfeats
            
            sel.diffdrthresh<-feature_rowindex%in%goodip
            
            max_absolute_standardized_centroids_thresh0<-pamr.res$max.discore.allfeats[goodip]
            
            selected_id_withmztime<-cbind(data_m_fc_withfeats[goodip,c(1:2)],pamr.res$pam_toplist,max_absolute_standardized_centroids_thresh0)
            ###savepamr.res,file="pamr.res.Rda")
            write.csv(selected_id_withmztime,file="dscores.selectedfeats.csv",row.names=FALSE)
            
            rank_vec<-rank(-discore_all)
            
            max_absolute_standardized_centroids_thresh0<-pamr.res$max.discore.allfeats
            
            data_limma_fdrall_withfeats<-cbind(max_absolute_standardized_centroids_thresh0,data_m_fc_withfeats)
            write.table(data_limma_fdrall_withfeats,file="Tables/pamr_ranked_feature_table.txt",sep="\t",row.names=FALSE)
            
            
          }else{
            goodip<-{}
            sel.diffdrthresh<-rep(FALSE,length(feature_rowindex))
          }
          
          rank_vec<-rank(-discore_all)
          
          pamr_ythresh<-pamr.res$max.discore.all.thresh-0.00000001
          
          
          
        }
        
        
        if(featselmethod=="rfesvm"){
          
          svm_classlabels<-classlabels[,1]
          
          if(analysismode=="classification"){
            svm_classlabels<-as.data.frame(svm_classlabels)
          }
          
          
          # ##save(data_m_fc,svm_classlabels,svm_kernel,file="svmdebug.Rda")
          if(length(class_labels_levels)<3){
            rfesvmres = diffexpsvmrfe(x=t(data_m_fc),y=svm_classlabels,svmkernel=svm_kernel)
            
            
            featureRankedList=rfesvmres$featureRankedList
            featureWeights=rfesvmres$featureWeights 
            #best_subset<-featureRankedList$best_subset
            
          }else{
            
            rfesvmres = diffexpsvmrfemulticlass(x=t(data_m_fc),y=svm_classlabels,svmkernel=svm_kernel)
            featureRankedList=rfesvmres$featureRankedList
            featureWeights=rfesvmres$featureWeights  
            
          }
          
          #  ##save(rfesvmres,file="rfesvmres.Rda")
          rank_vec<-seq(1,dim(data_m_fc_withfeats)[1])
          goodip<-featureRankedList[1:max_varsel]
          
          #dtemp1<-data_m_fc_withfeats[goodip,]
          
          sel.diffdrthresh<-rank_vec%in%goodip
          
          
          rank_vec<-sort(featureRankedList,index.return=TRUE)$ix
          
          weight_vec<-featureWeights #[rank_vec]
          
          data_limma_fdrall_withfeats<-cbind(featureWeights,data_m_fc_withfeats)
          
          
        }
        
        
        f1={}
        corfit={}
        
        if(featselmethod=="limma" | featselmethod=="limma1way")
        {
        #  save(classlabels,classlabels_orig,classlabels_dataframe,classlabels_response_mat,file="cldebug.Rda")
          
          classlabels_temp1<-classlabels
          
          classlabels<-classlabels_dataframe #classlabels_orig
          colnames(classlabels)<-c("SampleID","Factor1")
          if(alphabetical.order==FALSE){
            classlabels$Factor1<-factor(classlabels$Factor1,levels=unique(classlabels$Factor1))
            
            Factor1<-factor(classlabels$Factor1,levels=unique(classlabels$Factor1))
            
            
          }else{
            
            Factor1<-factor(classlabels$Factor1)
            
            
          }
          
          
          
          
          if(limma.contrasts.type=="contr.sum"){
            contrasts_factor1<-contr.sum(length(levels(factor(Factor1))))
            
            rownames(contrasts_factor1)<-levels(factor(Factor1))
            
            cnames_contr_factor1<-apply(contrasts_factor1,2,function(x){paste(names(x[which(abs(x)==1)]),collapse = "-")})
            
            
            
          }else{
            
            contrasts_factor1<-contr.treatment(length(levels(factor(Factor1))))
            
            rownames(contrasts_factor1)<-levels(factor(Factor1))
            
            cnames_contr_factor1<-apply(contrasts_factor1,2,function(x){paste(names(x[1]),names(x[which(abs(x)==1)]),sep = "-")})
            
            
            
          }
          
          
          
          
          colnames(contrasts_factor1)<-cnames_contr_factor1
          
          contrasts(Factor1) <- contrasts_factor1
          
          
          design <- model.matrix(~Factor1)
          
          
          classlabels<-classlabels_temp1
          
         # design <- model.matrix(~ -1+f)
          #colnames(design) <- levels(f)
          
          
          options(digit=3)
          #parameterNames<-colnames(design)
          design_mat_names=colnames(design)
          
          design_mat_names<-design_mat_names[-c(1)]
          # limma paired analysis
          if(pairedanalysis==TRUE){
            
            f1<-{}
            for(c in 1:length(class_labels_levels)){
              
              f1<-c(f1,seq(1,num_samps_group[[c]]))
              
            }
            
            print("Paired samples order")
            
            f1<-subject_inf
            print(subject_inf)
            print("Design matrix")
            print(design)
            
            ####savelist=ls(),file="limma.Rda")
            
            
            ##save(subject_inf,file="subject_inf.Rda")
            
            corfit<-duplicateCorrelation(data_m_fc,design=design,block=subject_inf,ndups=1)
            
            if(limmarobust==TRUE)
            {
              fit<-lmFit(data_m_fc,design,block=f1,cor=corfit$consensus,method="robust")
              
            }else{
              fit<-lmFit(data_m_fc,design,block=f1,cor=corfit$consensus)
            }
          }else{
            
            #not paired analysis
            if(limmarobust==TRUE)
            {
              fit <- lmFit(data_m_fc,design,method="robust")
              
            }else{
              fit <- lmFit(data_m_fc,design)
            }
            #fit<-treat(fit,lfc=lf)
            
          }
          
          cont.matrix=attributes(design)$contrasts
          #print(data_m_fc[1:3,])
          #fit2  <- contrasts.fit(fit, cont.matrix)
          
          #remove the intercept coefficient
          fit<-fit[,-1]
          fit2 <- eBayes(fit)
          
         
          
         # save(fit2,fit,data_m_fc,design,f1,corfit,classlabels,Factor1,cnames_contr_factor1,file="limma.eBayes.fit.Rda")
          # Various ways of summarising or plotting the results
          #topTable(fit,coef=2)
          
          
          #write.table(t1,file="topTable_limma.txt",sep="\t")
          
          
          if(dim(design)[2]>2){
            pvalues<-fit2$F.p.value
            p.value<-fit2$F.p.value
            
          }else{
            pvalues<-fit2$p.value
            p.value<-fit2$p.value
          }
          
          if(fdrmethod=="BH"){
            fdr_adjust_pvalue<-p.adjust(pvalues,method="BH")
          }else{
            if(fdrmethod=="ST"){
              
              fdr_adjust_pvalue<-try(qvalue(pvalues),silent=TRUE)
              
              if(is(fdr_adjust_pvalue,"try-error")){
                
                fdr_adjust_pvalue<-qvalue(pvalues,lambda=max(pvalues,na.rm=TRUE))
              }
              
              fdr_adjust_pvalue<-fdr_adjust_pvalue$qvalues
              
              
            }else{
              if(fdrmethod=="Strimmer"){
                pdf("fdrtool.pdf")
                fdr_adjust_pvalue<-suppressWarnings(fdrtool(as.vector(pvalues),statistic="pvalue",verbose=FALSE))
                fdr_adjust_pvalue<-fdr_adjust_pvalue$qval
                try(dev.off(),silent=TRUE)
              }else{
                if(fdrmethod=="none"){
                  fdr_adjust_pvalue<-pvalues
                  #fdr_adjust_pvalue<-p.adjust(pvalues,method="none")
                }else{
                  if(fdrmethod=="BY"){
                    fdr_adjust_pvalue<-p.adjust(pvalues,method="BY")
                  }else{
                    if(fdrmethod=="bonferroni"){
                      fdr_adjust_pvalue<-p.adjust(pvalues,method="bonferroni")
                    }
                  }
                  
                  
                }
              }
            }
            
          }
          
          if(dim(design)[2]<3){
            
            if(fdrmethod=="none"){
              filename<-paste("Tables/",parentfeatselmethod,"_pvalall_withfeats.txt",sep="")
            }else{
              filename<-paste("Tables/",parentfeatselmethod,"_fdrall_withfeats.txt",sep="")
            }
            cnames_tab<-colnames(data_m_fc_withfeats)
            cnames_tab<-c("P.value","adjusted.P.value",cnames_tab)
            
            
            data_limma_fdrall_withfeats<-cbind(p.value,fdr_adjust_pvalue,data_m_fc_withfeats)
            
            colnames(data_limma_fdrall_withfeats)<-as.character(cnames_tab)
            
            pvalues<-p.value
            
            #data_limma_fdrall_withfeats<-data_limma_fdrall_withfeats[order(fdr_adjust_fpvalue),]
           # write.table(data_limma_fdrall_withfeats,file=filename,sep="\t",row.names=FALSE)
            
            
            final.pvalues<-pvalues
            sel.diffdrthresh<-fdr_adjust_pvalue<fdrthresh & final.pvalues<pvalue.thresh
            
            
            
            goodip<-which(sel.diffdrthresh==TRUE)
            d4<-as.data.frame(data_limma_fdrall_withfeats)
            logp<-(-1)*log((d4[,1]+(10^-20)),10)
            
            #tiff("pval_dist.tiff",compression="lzw")
            #hist(d4[,1],xlab="p",main="Distribution of p-values")
            #dev.off()
            
            
          }else{
            
            
            adjusted.P.value<-fdr_adjust_pvalue
            if(limmadecideTests==TRUE){
              results2<-decideTests(fit2,method="nestedF",adjust.method="BH",p.value=fdrthresh)
              #tiff("comparison_contrast_overlap.tiff",width=plots.width,height=plots.height,res=plots.res, compression="lzw")
              #if(length(class_labels_levels)<4){
               if(ncol(results2)<5){ 
                if(output.device.type!="pdf"){
                  
                  temp_filename_5<-"Figures/LIMMA_venn_diagram.png"
                  
                  png(temp_filename_5,width=plots.width,height=plots.height,res=plots.res,type=plots.type,units="in")
                }
                
                
                vennDiagram(results2,cex=0.8)
                
                if(output.device.type!="pdf"){
                  
                  try(dev.off(),silent=TRUE)
                }
                
              }
            }else{
              #dev.off()
              
              results2<-fit2$p.value[,-c(1)]
            }
            cnames_tab<-colnames(data_m_fc_withfeats)
            cnames_tab2<-colnames(results2)
            
            cnames_tab<-c("P.value","adjusted.P.value",cnames_tab2,cnames_tab)
            
            data_limma_fdrall_withfeats<-cbind(p.value,adjusted.P.value,results2,data_m_fc_withfeats)
            data_limma_fdrall_withfeats<-as.data.frame(data_limma_fdrall_withfeats)
            
            if(limmarobust==FALSE){
              filename<-"Tables/limma_posthoc1wayanova_results.txt"
            }else{
              filename<-"Tables/limmarobust_posthoc1wayanova_results.txt"
              
            }
            colnames(data_limma_fdrall_withfeats)<-as.character(cnames_tab)
            
            
            if(length(check_names)>0){
              
              data_limma_fdrall_withfeats<-cbind(p.value,adjusted.P.value,results2,data_m_fc_with_names,data_m_fc_withfeats[,-c(1:2)])
              data_limma_fdrall_withfeats<-as.data.frame(data_limma_fdrall_withfeats)
              #data_limma_fdrall_withfeats<-cbind(p.value,adjusted.p.value,results2,data_m_fc_with_names,data_m_fc_withfeats[,-c(1:2)])
              
              rem_col_ind1<-grep(colnames(data_limma_fdrall_withfeats),pattern=c("mz"))
              
              rem_col_ind2<-grep(colnames(data_limma_fdrall_withfeats),pattern=c("time"))
              
              rem_col_ind<-c(rem_col_ind1,rem_col_ind2)
              
            }else{
              rem_col_ind<-{}
            }
            
            if(length(rem_col_ind)>0){
              write.table(data_limma_fdrall_withfeats[,-c(rem_col_ind)], file=filename,sep="\t",row.names=FALSE)
            }else{
              
              write.table(data_limma_fdrall_withfeats, file=filename,sep="\t",row.names=FALSE)
            }
            
            #data_limma_fdrall_withfeats<-data_limma_fdrall_withfeats[order(fdr_adjust_fpvalue),]
            write.table(data_limma_fdrall_withfeats, file=filename,sep="\t",row.names=FALSE)
            
            data_limma_fdrall_withfeats<-cbind(p.value,adjusted.P.value,data_m_fc_withfeats)
            
            
            if(fdrmethod=="none"){
              filename<-paste("limma_posthoc1wayanova_pval",fdrthresh,"_results.txt",sep="")
            }else{
              filename<-paste("limma_posthoc1wayanova_fdr",fdrthresh,"_results.txt",sep="")
            }
            if(length(which(data_limma_fdrall_withfeats$adjusted.P.value<fdrthresh & data_limma_fdrall_withfeats$p.value<pvalue.thresh))>0){
              data_limma_sig_withfeats<-data_limma_fdrall_withfeats[data_limma_fdrall_withfeats$adjusted.P.value<fdrthresh & data_limma_fdrall_withfeats$p.value<pvalue.thresh,]
              #write.table(data_limma_sig_withfeats, file=filename,sep="\t",row.names=FALSE)
            }
            
            # data_limma_fdrall_withfeats<-cbind(p.value,adjusted.p.value,data_m_fc_withfeats)
            
            data_limma_fdrall_withfeats<-cbind(pvalues,fdr_adjust_pvalue,data_m_fc_withfeats)
            final.pvalues<-pvalues
            
            cnames_tab<-colnames(data_m_fc_withfeats)
            cnames_tab<-c("P.value","adjusted.P.value",cnames_tab)
            colnames(data_limma_fdrall_withfeats)<-as.character(cnames_tab)
          }
          
          #pvalues<-data_limma_fdrall_withfeats$p.value
          
          #final.pvalues<-pvalues
          
          # print("checking here")
          sel.diffdrthresh<-fdr_adjust_pvalue<fdrthresh & final.pvalues<pvalue.thresh
          
          goodip<-which(sel.diffdrthresh==TRUE)
          d4<-as.data.frame(data_limma_fdrall_withfeats)
          logp<-(-1)*log((d4[,1]+(10^-20)),10)
          
          #tiff("pval_dist.tiff",compression="lzw")
          #hist(d4[,1],xlab="p",main="Distribution of p-values")
          #dev.off()
          
          
        }
        
        
        if(featselmethod=="limma2way")
        {
          #design <- cbind(Grp1vs2=c(rep(1,num_samps_group[[1]]),rep(0,num_samps_group[[2]])),Grp2vs1=c(rep(0,num_samps_group[[1]]),rep(1,num_samps_group[[2]])))
         # print("here")
          
         # save(f,sampleclass,data_m_fc,classlabels,classlabels_orig,file="limma2way.Rda")
          classlabels_temp<-classlabels
          
          
          colnames(classlabels_orig)<-c("SampleID","Factor1","Factor2")
          classlabels<- classlabels_orig #classlabels_dataframe #
          colnames(classlabels)<-c("SampleID","Factor1","Factor2")
          
          #design <- model.matrix(~ -1+f)
          
         
          
          #classlabels<-read.table("/Users/karanuppal/Documents/Emory/JonesLab/Projects/DifferentialExpression/xmsPaNDA/examples_and_manual/Example_feature_table_and_classlabels/classlabels_two_way_anova.txt",sep="\t",header=TRUE)
          #classlabels<-classlabels[order(classlabels$Factor2,decreasing = T),]
          
         if(alphabetical.order==FALSE){
          classlabels$Factor1<-factor(classlabels$Factor1,levels=unique(classlabels$Factor1))
          classlabels$Factor2<-factor(classlabels$Factor2,levels=unique(classlabels$Factor2))
          Factor1<-factor(classlabels$Factor1,levels=unique(classlabels$Factor1))
          Factor2<-factor(classlabels$Factor2,levels=unique(classlabels$Factor2))
          
         }else{
           
           Factor1<-factor(classlabels$Factor1)
           
           Factor2<-factor(classlabels$Factor2)
           
         }
         
          
         
          #this will create sum to zero parametrization. Coefficient Comparison Interpretation
          #contrasts(Strain) <- contr.sum(2)
          #contrasts(Treatment) <- contr.sum(2)
          #design <- model.matrix(~Strain*Treatment)
          #Intercept (WT.U+WT.S+Mu.U+Mu.S)/4; Grand mean
          #Strain1 (WT.U+WT.S-Mu.U-Mu.S)/4; strain main effect
          #Treatment1 (WT.U-WT.S+Mu.U-Mu.S)/4; treatment main effect
          #Strain1:Treatment1 (WT.U-WT.S-Mu.U+Mu.S)/4; Interaction
          
         
         
          
          if(limma.contrasts.type=="contr.sum"){
              contrasts_factor1<-contr.sum(length(levels(factor(Factor1))))
              contrasts_factor2<-contr.sum(length(levels(factor(Factor2))))
              
              rownames(contrasts_factor1)<-levels(factor(Factor1))
              rownames(contrasts_factor2)<-levels(factor(Factor2))
              
              cnames_contr_factor1<-apply(contrasts_factor1,2,function(x){paste(names(x[which(abs(x)==1)]),collapse = "-")})
              
              cnames_contr_factor2<-apply(contrasts_factor2,2,function(x){paste(names(x[which(abs(x)==1)]),collapse = "-")})
              
              
          }else{
            
            contrasts_factor1<-contr.treatment(length(levels(factor(Factor1))))
            contrasts_factor2<-contr.treatment(length(levels(factor(Factor2))))
            
            rownames(contrasts_factor1)<-levels(factor(Factor1))
            rownames(contrasts_factor2)<-levels(factor(Factor2))
            
            cnames_contr_factor1<-apply(contrasts_factor1,2,function(x){paste(names(x[1]),names(x[which(abs(x)==1)]),sep = "-")})
            
            cnames_contr_factor2<-apply(contrasts_factor2,2,function(x){paste(names(x[1]),names(x[which(abs(x)==1)]),sep= "-")})
            
            
          }
          
         
          
         
          colnames(contrasts_factor1)<-cnames_contr_factor1
          colnames(contrasts_factor2)<-cnames_contr_factor2
          
          contrasts(Factor1) <- contrasts_factor1
          contrasts(Factor2) <- contrasts_factor2
        
          design <- model.matrix(~Factor1*Factor2)
          
         # fit<-lmFit(data_m_fc,design=design)
          
        
          #2. this will create contrasts with respect to the reference group (first level in each factor)
        if(FALSE){
            contrasts(Factor1) <- contr.treatment(length(levels(factor(Factor1))))
          contrasts(Factor2) <- contr.treatment(length(levels(factor(Factor2))))
          design.trt <- model.matrix(~Factor1*Factor2)
          
          fit.trt<-lmFit(data_m_fc,design=design.trt)
          
          s1=apply(fit.trt$coefficients,2,function(x){
            length(which(is.na(x))==TRUE)/length(x)
          })
          
        }
          
          #3. this will create design matrix with all factors
          call<-lapply(classlabels[,c(2:3)],contrasts,contrasts=FALSE)
          design.all<-model.matrix(~Factor1*Factor2,data=classlabels,contrasts.arg=call)
          
          #grand mean: mean of means (mean of each level)
          #mean_per_level<-lapply(2:ncol(design.all),function(x){mean(data_m_fc[1,which(design.all[,x]==1)])})
          #mean_per_level<-unlist(mean_per_level)
          #names(mean_per_level)<-colnames(design.all[,-1])
          #grand_mean<-mean(mean_per_level,na.rm=TRUE)
          
          #grand_mean<-with(d,tapply(data_m_fc[1,],list(Factor1,Factor2),mean))
          
          
          
         # colnames(design)<-gsub(colnames(design),pattern="Factor1",replacement="")
          #colnames(design)<-gsub(colnames(design),pattern="Factor2",replacement="")
          
          
       #   save(design,f,sampleclass,data_m_fc,classlabels,classlabels_orig,file="limma2way.Rda")
          
          
          classlabels<-classlabels_temp
          
          # print(data_m_fc[1:4,])
          #colnames(design) <- levels(f)
          #colnames(design)<-levels(factor(sampleclass))
          
          
          options(digit=3)
          parameterNames<-colnames(design)
          
          print("Design matrix")
          print(design)
          
          
          
          if(pairedanalysis==TRUE)
          {
            
            
            print("Paired design is")
            #print(f1)
            print(subject_inf)
            
            f1<-subject_inf
            
            #print(data_m_fc[1:10,1:10])
            
            save(design,subject_inf,file="limmadesign.Rda")
          }
          
          if(dim(design)[2]>=1){
            
            #cont.matrix <- makeContrasts(Grp1vs2="ClassA-ClassB",Grp1vs3="ClassC-ClassD",Grp2vs3=("ClassA-ClassB")-("ClassC-ClassD"),levels=design)
            #cont.matrix <- makeContrasts(Grp1vs2=ClassA-ClassB,Grp1vs3=ClassC-ClassD,Grp2vs3=(ClassA-ClassB)-(ClassC-ClassD),Grp3vs4=ClassA-ClassC,Group2vs4=ClassB-ClassD,levels=design)
            
            #cont.matrix <- makeContrasts(Factor1=(ClassA+ClassB)-(ClassC+ClassD),Factor2=(ClassA+ClassC)-(ClassB+ClassD),Factor1x2=(ClassA-ClassB)-(ClassC-ClassD),levels=design)
            
            design.pairs <-
              function(levels) {
                n <- length(levels)
                design <- matrix(0,n,choose(n,2))
                rownames(design) <- levels
                colnames(design) <- 1:choose(n,2)
                k <- 0
                for (i in 1:(n-1))
                  for (j in (i+1):n) {
                    k <- k+1
                    design[i,k] <- 1
                    design[j,k] <- -1
                    colnames(design)[k] <- paste(levels[i],"-",levels[j],sep="")
                  }
                design
              }
            
            
            #levels_1<-levels(factor(classlabels[,2]))
            #levels_2<-levels(factor(classlabels[,3]))
            
            #design2<-design.pairs(c(as.character(levels_1),as.character(levels_2)))
            
          
            #cont.matrix<-makeContrasts(contrasts=colnames(design2),levels=c(as.character(levels_1),as.character(levels_2)))
            
            if(pairedanalysis==TRUE){
              
              #class_table_facts<-table(classlabels)
              
              #f1<-c(seq(1,num_samps_group[[1]]),seq(1,num_samps_group[[2]]),seq(1,num_samps_group[[1]]),seq(1,num_samps_group[[2]]))
              
              
              corfit<-duplicateCorrelation(data_m_fc,design=design,block=subject_inf,ndups=1)
              
              #print(f1)
              
              if(limmarobust==TRUE)
              {
                fit<-lmFit(data_m_fc,design,block=f1,cor=corfit$consensus,method="robust")
              }else
              {   
                
                fit<-lmFit(data_m_fc,design,block=f1,cor=corfit$consensus)
              }
              
              s1=apply(fit$coefficients,2,function(x){
                length(which(is.na(x))==TRUE)/length(x)
              })
              
              if(length(which(s1==1))>0){
                design<-design[,-which(s1==1)]
                #fit <- lmFit(data_m_fc,design)
                
                if(limmarobust==TRUE)
                {
                  fit<-lmFit(data_m_fc,design,block=f1,cor=corfit$consensus,method="robust")
                }else{
                  fit<-lmFit(data_m_fc,design,block=f1,cor=corfit$consensus)
              
                }  
              }
              
            }
            else{
              
              
             # fit <- lmFit(data_m_fc,design)
              
              if(limmarobust==TRUE)
              {
                fit<-lmFit(data_m_fc,design,method="robust")
              }else{
                fit <- lmFit(data_m_fc,design)
              }
                
                s1=apply(fit$coefficients,2,function(x){
                     length(which(is.na(x))==TRUE)/length(x)
                 })
              
                if(length(which(s1==1))>0){
                    design<-design[,-which(s1==1)]
                    
                    if(limmarobust==TRUE)
                    {
                      fit<-lmFit(data_m_fc,design,method="robust")
                      
                      
                    }else{
                      
                        fit<-lmFit(data_m_fc,design)
                    }
                        
                }
                
              
              
      
            }
            
            
          
            }
          
          
          print("Limma classlabels")
          
          print(head(classlabels))
          print(head(classlabels_orig))
          
          fit<-fit[,-1]
          fit2=eBayes(fit)
          
          
          
          results <- topTableF(fit2, n=Inf)
         # decideresults<-decideTests(fit2)
          
          
          # Ordinary fit
          
        #  save(fit2,fit,results,file="limma.eBayes.fit.Rda")
          
          
          #fit2  <- contrasts.fit(fit, cont.matrix)
          
          #fit2 <- eBayes(fit2)
          #as.data.frame(fit2[1:10,])
          
          # Various ways of summarising or plotting the results
          #topTable(fit2,coef=2)
          
          #    ##save(fit2,file="fit2.Rda")
          
          if(dim(design)[2]>2){
            pvalues<-fit2$F.p.value
            p.value<-fit2$F.p.value
            
            
            
          }else{
            pvalues<-fit2$p.value
            p.value<-fit2$p.value
          }
          
          if(fdrmethod=="BH"){
            fdr_adjust_pvalue<-p.adjust(pvalues,method="BH")
          }else{
            if(fdrmethod=="ST"){
              #fdr_adjust_pvalue<-qvalue(pvalues)
              #fdr_adjust_pvalue<-fdr_adjust_pvalue$qvalues
              
              fdr_adjust_pvalue<-try(qvalue(pvalues),silent=TRUE)
              
              if(is(fdr_adjust_pvalue,"try-error")){
                
                fdr_adjust_pvalue<-qvalue(pvalues,lambda=max(pvalues,na.rm=TRUE))
              }
              
              fdr_adjust_pvalue<-fdr_adjust_pvalue$qvalues
              
              
            }else{
              if(fdrmethod=="Strimmer"){
                pdf("fdrtool.pdf")
                #par_rows=1
                #par(mfrow=c(par_rows,1))
                fdr_adjust_pvalue<-suppressWarnings(fdrtool(as.vector(pvalues),statistic="pvalue",verbose=FALSE))
                fdr_adjust_pvalue<-fdr_adjust_pvalue$qval
                try(dev.off(),silent=TRUE)
              }else{
                if(fdrmethod=="none"){
                  #	fdr_adjust_pvalue<-pvalues
                  fdr_adjust_pvalue<-p.adjust(pvalues,method="none")	
                }else{
                  if(fdrmethod=="BY"){
                    fdr_adjust_pvalue<-p.adjust(pvalues,method="BY")
                  }else{
                    if(fdrmethod=="bonferroni"){
                      fdr_adjust_pvalue<-p.adjust(pvalues,method="bonferroni")
                    }
                  }
                }
              }
            }
            
          }
          
          #print("Doing this:")
          
          adjusted.p.value<-fdr_adjust_pvalue
          
          data_limma_fdrall_withfeats<-cbind(p.value,adjusted.p.value,data_m_fc_withfeats)
          
          if(limmadecideTests==TRUE){
            results2<-decideTests(fit2,adjust.method="BH",method="nestedF",p.value=fdrthresh) #
            #tiff("comparison_contrast_overlap.tiff",width=plots.width,height=plots.height,res=plots.res, compression="lzw")
            
           # save(results2,file="results2.Rda")
            
            
            cnames_tab<-colnames(data_m_fc_withfeats)
            cnames_tab2<-colnames(results2)
            cnames_tab<-c("P.value","adjusted.P.value",cnames_tab2,cnames_tab)
            
            data_limma_fdrall_withfeats<-cbind(p.value,adjusted.p.value,results2,data_m_fc_withfeats)
            if(limmarobust==FALSE){
              filename<-"Tables/limma_2wayposthoc_decideresults.txt"
            }else{
              
              filename<-"Tables/limmarobust_2wayposthoc_decideresults.txt"
            }
            colnames(data_limma_fdrall_withfeats)<-as.character(cnames_tab)
            
            #data_limma_fdrall_withfeats<-data_limma_fdrall_withfeats[order(fdr_adjust_fpvalue),]
         #   write.table(data_limma_fdrall_withfeats, file=filename,sep="\t",row.names=FALSE)
            
            
            #if(length(class_labels_levels)<5){
             if(ncol(results2)<6){ 
              
              if(output.device.type!="pdf"){
                
                temp_filename_1<-"Figures/LIMMA_venn_diagram.png"
                
                png(temp_filename_1,width=plots.width,height=plots.height,res=plots.res,type=plots.type,units="in")
              }
              
              vennDiagram(results2,cex=0.8)
              
              if(output.device.type!="pdf"){
                
                try(dev.off(),silent=TRUE)
              }
              
              
             }
            
            
            
          }
          else{
            #dev.off()
            
            results2<-fit2$p.value[,-c(1)]
          }
          
            cnames_tab<-colnames(data_m_fc_withfeats)
            cnames_tab2<-colnames(results2)
            cnames_tab<-c("P.value","adjusted.P.value",cnames_tab2,cnames_tab)
            
            
            #save(data_m_fc_withfeats,names)
            data_limma_fdrall_withfeats<-cbind(p.value,adjusted.p.value,results2,data_m_fc_withfeats)
            if(limmarobust==FALSE){
              filename<-"Tables/limma_2wayposthoc_pvalues.txt"
            }else{
              
              filename<-"Tables/limmarobust_2wayposthoc_pvalues.txt"
            }
            colnames(data_limma_fdrall_withfeats)<-as.character(cnames_tab)
            
            
            if(length(check_names)>0){
            
              
              data_limma_fdrall_withfeats<-cbind(p.value,adjusted.p.value,results2,data_m_fc_with_names,data_m_fc_withfeats[,-c(1:2)])
              
              rem_col_ind1<-grep(colnames(data_limma_fdrall_withfeats),pattern=c("mz"))
              
              rem_col_ind2<-grep(colnames(data_limma_fdrall_withfeats),pattern=c("time"))
              
              rem_col_ind<-c(rem_col_ind1,rem_col_ind2)
              
            }else{
              rem_col_ind<-{}
            }
            
            if(length(rem_col_ind)>0){
              write.table(data_limma_fdrall_withfeats[,-c(rem_col_ind)], file=filename,sep="\t",row.names=FALSE)
            }else{
              
              write.table(data_limma_fdrall_withfeats, file=filename,sep="\t",row.names=FALSE)
            }
            
            
            #data_limma_fdrall_withfeats<-data_limma_fdrall_withfeats[order(fdr_adjust_fpvalue),]
           # write.table(data_limma_fdrall_withfeats, file=filename,sep="\t",row.names=FALSE)
            
            
          
          #tiff("comparison_contrast_overlap.tiff",width=plots.width,height=plots.height,res=plots.res, compression="lzw")
          
          #dev.off()
          
          #results2<-fit2$p.value
          
         
          
          
          
          classlabels_orig<-as.data.frame(classlabels_orig)
          
          
          
        
          
          data_limma_fdrall_withfeats<-cbind(p.value,adjusted.p.value,data_m_fc_withfeats)
          
          # data_limma_fdrall_withfeats<-cbind(pvalues,fdr_adjust_pvalue,data_m_fc_withfeats)
          
          cnames_tab<-colnames(data_m_fc_withfeats)
          cnames_tab<-c("P.value","adjusted.P.value",cnames_tab)
          colnames(data_limma_fdrall_withfeats)<-as.character(cnames_tab)
          #write.table(data_limma_fdrall_withfeats,file="Limma_posthoc2wayanova_results.txt",sep="\t",row.names=FALSE)
          #print("checking here")
          pvalues<-p.value
          
          final.pvalues<-pvalues
          
          sel.diffdrthresh<-fdr_adjust_pvalue<fdrthresh & final.pvalues<pvalue.thresh
          
          goodip<-which(sel.diffdrthresh==TRUE)
          d4<-as.data.frame(data_limma_fdrall_withfeats)
          logp<-(-1)*log((d4[,1]+(10^-20)),10)
          
          
          
          #results2<-decideTests(fit2,method="nestedF",adjust.method=fdrmethod,p.value=fdrthresh)
          
          
          
        }
        
        
        
        if(featselmethod=="RF")
        {
          
          maxint<-apply(data_m_fc,1,max)
          
          
          data_m_fc_withfeats<-as.data.frame(data_m_fc_withfeats)
          
          data_m_fc<-as.data.frame(data_m_fc)
          #write.table(classlabels,file="classlabels_rf.txt",sep="\t",row.names=FALSE)
          
          save(data_m_fc,classlabels,numtrees,analysismode,file="rfdebug.Rda")
          
          
          
          
          if(rfconditional==TRUE){
            
            print("Performing random forest analysis using the cforest")
            
            #rfcondres1<-do_rf_conditional(X=data_m_fc,rf_classlabels,ntrees=numtrees,analysismode) #,silent=TRUE)
            filename<-"RFconditional_VIM_allfeats.txt"
          }else{
            
            
            
            
            #varimp_res2<-do_rf(X=data_m_fc,classlabels=rf_classlabels,ntrees=numtrees,analysismode)
            
            if(analysismode=="classification"){
              rf_classlabels<-classlabels[,1]
              print("Performing random forest analysis using the randomForest and Boruta functions")
              varimp_res2<-do_rf_boruta(X=data_m_fc,classlabels=rf_classlabels) #,ntrees=numtrees,analysismode)
              filename<-"RF_VIM_Boruta_allfeats.txt"
            }else{
              rf_classlabels<-classlabels
              print("Performing random forest analysis using the randomForest function")
              varimp_res2<-do_rf(X=data_m_fc,classlabels=rf_classlabels,ntrees=numtrees,analysismode)
            }
            
          }
          
          varimp_res3<-cbind(data_m_fc_withfeats[,c(1:2)],varimp_res2)
          
          filename<-paste("Tables/",filename,sep="")
          write.table(varimp_res3, file=filename,sep="\t",row.names=TRUE)
          
          
          
          goodip<-which(varimp_res2>0)
          
          if(length(goodip)<1){
            print("No features were selected using the selection criteria.")
          }
          var_names<-paste(sprintf("%.3f",data_m_fc_withfeats[,1]),sprintf("%.1f",data_m_fc_withfeats[,2]),sep="_")
          
          names(varimp_res2)<-as.character(var_names)
          sel.diffdrthresh<-varimp_res2>0
          
          if(length(which(sel.diffdrthresh==TRUE))<1){
            print("No features were selected using the selection criteria")
          }
          
          
          
          num_var_rf<-length(which(sel.diffdrthresh==TRUE))
          
          if(num_var_rf>10){
            
            num_var_rf=10
          }
          sorted_varimp_res<-varimp_res2[order(varimp_res2,decreasing=TRUE)[1:(num_var_rf)]]
          
          sorted_varimp_res<-sort(sorted_varimp_res)
          
          barplot_text=paste("Variable Importance measures (top ",length(sorted_varimp_res)," shown)\n",sep="")
          
          if(output.device.type!="pdf"){
            
            temp_filename_1<-"Figures/RF_selectfeats_VIMbarplot.png"
            
            png(temp_filename_1,width=plots.width,height=plots.height,res=plots.res,type=plots.type,units="in")
          }
          
          #    ##save(varimp_res2,data_m_fc,rf_classlabels,sorted_varimp_res,file="test_rf.Rda")
          
          barplot(sorted_varimp_res, xlab="Selected features", main=barplot_text,cex.axis=0.5,cex.names=0.4, ylab="VIM",range(pretty(c(0,sorted_varimp_res))),space=0.1)
          
          if(output.device.type!="pdf"){
            
            try(dev.off(),silent=TRUE)
          }
          
          
          
          
          rank_num<-rank(-varimp_res2)
          
          data_limma_fdrall_withfeats<-cbind(varimp_res2,rank_num,data_m_fc_withfeats)
          
          cnames_tab<-colnames(data_m_fc_withfeats)
          cnames_tab<-c("VIM","Rank",cnames_tab)
          
          goodip<-which(sel.diffdrthresh==TRUE)
          
          feat_sigfdrthresh[lf]<-length(which(sel.diffdrthresh==TRUE))
          
          colnames(data_limma_fdrall_withfeats)<-as.character(cnames_tab)
          
          #data_limma_fdrall_withfeats<-data_limma_fdrall_withfeats[order(fdr_adjust_fpvalue),]
          #write.table(data_limma_fdrall_withfeats, file=filename,sep="\t",row.names=FALSE)
          
          
        }
        
        if(featselmethod=="MARS"){
          
          mars_classlabels<-classlabels[,1]
          marsres1<-do_mars(X=data_m_fc,mars_classlabels, analysismode,kfold)
          
          varimp_marsres1<-marsres1$mars_varimp
          
          mars_mznames<-rownames(varimp_marsres1)
          
          
          all_names<-paste("mz",seq(1,dim(data_m_fc)[1]),sep="")
          
          com1<-match(all_names,mars_mznames)
          
          
          filename<-"MARS_variable_importance.txt"
          
          
          if(is.na(max_varsel)==FALSE){
            
            if(max_varsel>dim(data_m_fc)[1]){
              max_varsel=dim(data_m_fc)[1]
            }
            varimp_res2<-varimp_marsres1[,4]
            
            #sort by VIM; and keep the top max_varsel scores
            sorted_varimp_res<-varimp_res2[order(varimp_res2,decreasing=TRUE)[1:(max_varsel)]]
            
            #get the minimum VIM from the top max_varsel scores
            min_thresh<-min(sorted_varimp_res[which(sorted_varimp_res>=mars.gcv.thresh)],na.rm=TRUE)
            
            
            row_num_vec<-seq(1,length(varimp_res2))
            
            #only use the top max_varsel scores
            #goodip<-order(varimp_res2,decreasing=TRUE)[1:(max_varsel)]
            #sel.diffdrthresh<-row_num_vec%in%goodip
            
            #use a threshold of mars.gcv.thresh
            sel.diffdrthresh<-varimp_marsres1[,4]>=min_thresh
            
            goodip<-which(sel.diffdrthresh==TRUE)
            
          }else{
            
            #use a threshold of mars.gcv.thresh
            sel.diffdrthresh<-varimp_marsres1[,4]>=mars.gcv.thresh
            
            goodip<-which(sel.diffdrthresh==TRUE)
          }
          
          
          num_var_rf<-length(which(sel.diffdrthresh==TRUE))
          
          if(num_var_rf>10){
            
            num_var_rf=10
          }
          sorted_varimp_res<-varimp_res2[order(varimp_res2,decreasing=TRUE)[1:(num_var_rf)]]
          
          sorted_varimp_res<-sort(sorted_varimp_res)
          
          barplot_text=paste("Generalized cross validation (top ",length(sorted_varimp_res)," shown)\n",sep="")
          
          if(output.device.type!="pdf"){
            
            temp_filename_1<-"Figures/MARS_selectfeats_GCVbarplot.png"
            
            png(temp_filename_1,width=plots.width,height=plots.height,res=plots.res,type=plots.type,units="in")
          }
          
          
          barplot(sorted_varimp_res, xlab="Selected features", main=barplot_text,cex.axis=0.5,cex.names=0.4, ylab="GCV",range(pretty(c(0,sorted_varimp_res))),space=0.1)
          
          if(output.device.type!="pdf"){
            
            try(dev.off(),silent=TRUE)
          }
          
          
          data_limma_fdrall_withfeats<-cbind(varimp_marsres1[,c(4,6)],data_m_fc_withfeats)
          
          cnames_tab<-colnames(data_m_fc_withfeats)
          cnames_tab<-c("GCV importance","RSS importance",cnames_tab)
          feat_sigfdrthresh[lf]<-length(which(sel.diffdrthresh==TRUE))
          
          colnames(data_limma_fdrall_withfeats)<-as.character(cnames_tab)
          
          
          goodip<-which(sel.diffdrthresh==TRUE)
          
          
          
        }
        
        if(featselmethod=="pls" | featselmethod=="o1pls" | featselmethod=="o2pls" | featselmethod=="spls" | featselmethod=="o1spls" | featselmethod=="o2spls")
        {
          
          
          
          classlabels<-as.data.frame(classlabels)
          
          
          if(is.na(max_comp_sel)==TRUE){
            max_comp_sel=pls_ncomp
          }
          
          rand_pls_sel<-{} #new("list")
          if(featselmethod=="spls" | featselmethod=="o1spls" | featselmethod=="o2spls"){
            
            
            if(featselmethod=="o1spls"){
              
              featselmethod="o1pls"
              
            }else{
              
              if(featselmethod=="o2spls"){
                featselmethod="o2pls"
              }
              
            }
            
            if(pairedanalysis==TRUE){
              
              classlabels_temp<-cbind(classlabels_sub[,2],classlabels)
              
              
              set.seed(999)
              
              plsres1<-do_plsda(X=data_m_fc,Y=classlabels_sub,oscmode=featselmethod,numcomp=pls_ncomp,kfold=kfold,evalmethod=pred.eval.method,keepX=max_varsel,sparseselect=TRUE,
                                analysismode,sample.col.opt=sample.col.opt,sample.col.vec=col_vec,scoreplot_legend=scoreplot_legend,pairedanalysis=pairedanalysis,
                                optselect=optselect,class_labels_levels_main=class_labels_levels_main,legendlocation=legendlocation,output.device.type=output.device.type,
                                plots.res=plots.res,plots.width=plots.width,plots.height=plots.height,plots.type=plots.type,pls.ellipse=pca.ellipse)
              
              if (is(plsres1, "try-error")){
                print(paste("sPLS could not be performed at RSD threshold: ",log2.fold.change.thresh,sep=""))
                #break;
              }
              
              opt_comp<-plsres1$opt_comp
              #for(randindex in 1:100)
              save(plsres1,file="plsres1.Rda")
              
              if(is.na(pls.permut.count)==FALSE){
                set.seed(999)
                seedvec<-runif(pls.permut.count,10,10*pls.permut.count)
                
                
                
                if(pls.permut.count>0){
                  
                  cl <- parallel::makeCluster(getOption("cl.cores", num_nodes))
                  clusterEvalQ(cl,library(plsgenomics))
                  clusterEvalQ(cl,library(dplyr))
                  
                  clusterEvalQ(cl,library(plyr))
                  clusterExport(cl,"pls.lda.cv",envir = .GlobalEnv)
                  clusterExport(cl,"plsda_cv",envir = .GlobalEnv)
                  #clusterExport(cl,"%>%",envir = .GlobalEnv) #%>%
                  clusterExport(cl,"do_plsda_rand",envir = .GlobalEnv)
                  clusterEvalQ(cl,library(mixOmics))
                  clusterEvalQ(cl,library(pls))
                  
                  
                  rand_pls_sel<-parLapply(cl,1:pls.permut.count,function(x)
                  {
                    
                    set.seed(seedvec[x])
                    
                    
                    plsresrand<-do_plsda_rand(X=data_m_fc,Y=classlabels_sub[sample(x=seq(1,dim(classlabels_sub)[1]),size=dim(classlabels_sub)[1]),],oscmode=featselmethod,numcomp=opt_comp,kfold=kfold,evalmethod=pred.eval.method,keepX=max_varsel,sparseselect=TRUE,analysismode,sample.col.vec=col_vec,scoreplot_legend=scoreplot_legend,pairedanalysis=pairedanalysis,optselect=FALSE,class_labels_levels_main=class_labels_levels_main,legendlocation=legendlocation,plotindiv=FALSE) #,silent=TRUE)
                    
                    #rand_pls_sel<-cbind(rand_pls_sel,plsresrand$vip_res[,1])
                    if (is(plsresrand, "try-error")){
                      
                      return(rep(0,dim(data_m_fc)[1]))
                    }else{
                      return(plsresrand$vip_res[,1])
                    }
                  })
                  
                  stopCluster(cl)
                }
              }
              
              
              
            }else{	
              #plsres1<-try(do_plsda(X=data_m_fc,Y=classlabels,oscmode=featselmethod,numcomp=pls_ncomp,kfold=kfold,evalmethod=pred.eval.method,keepX=max_varsel,sparseselect=TRUE,analysismode,sample.col.vec=col_vec,scoreplot_legend=scoreplot_legend,pairedanalysis=pairedanalysis,optselect=optselect,class_labels_levels_main=class_labels_levels_main,legendlocation=legendlocation,pls.vip.selection=pls.vip.selection),silent=TRUE)
              
              #  ##save(data_m_fc,classlabels,pls_ncomp,kfold,file="pls1.Rda")
              set.seed(999)
              plsres1<-do_plsda(X=data_m_fc,Y=classlabels,oscmode=featselmethod,numcomp=pls_ncomp,kfold=kfold,evalmethod=pred.eval.method,
                                keepX=max_varsel,sparseselect=TRUE,analysismode,sample.col.opt=sample.col.opt,sample.col.vec=col_vec,
                                scoreplot_legend=scoreplot_legend,pairedanalysis=pairedanalysis,optselect=optselect,
                                class_labels_levels_main=class_labels_levels_main,legendlocation=legendlocation,
                                pls.vip.selection=pls.vip.selection,output.device.type=output.device.type,
                                plots.res=plots.res,plots.width=plots.width,plots.height=plots.height,plots.type=plots.type,pls.ellipse=pca.ellipse)
              
              opt_comp<-plsres1$opt_comp
              
              if (is(plsres1, "try-error")){
                print(paste("sPLS could not be performed at RSD threshold: ",log2.fold.change.thresh,sep=""))
                break;
              }
              #for(randindex in 1:100)
              if(is.na(pls.permut.count)==FALSE){
                
                set.seed(999)
                seedvec<-runif(pls.permut.count,10,10*pls.permut.count)
                
                
                
                if(pls.permut.count>0){
                  
                  cl <- parallel::makeCluster(getOption("cl.cores", num_nodes))
                  clusterEvalQ(cl,library(plsgenomics))
                  clusterEvalQ(cl,library(dplyr))
                  
                  clusterEvalQ(cl,library(plyr))
                  clusterExport(cl,"pls.lda.cv",envir = .GlobalEnv)
                  clusterExport(cl,"plsda_cv",envir = .GlobalEnv)
                  #clusterExport(cl,"%>%",envir = .GlobalEnv) #%>%
                  clusterExport(cl,"do_plsda_rand",envir = .GlobalEnv)
                  clusterEvalQ(cl,library(mixOmics))
                  clusterEvalQ(cl,library(pls))
                  
                  
                  rand_pls_sel<-parLapply(cl,1:pls.permut.count,function(x)
                  {
                    
                    set.seed(seedvec[x])
                    
                    
                    plsresrand<-do_plsda_rand(X=data_m_fc,Y=classlabels[sample(x=seq(1,dim(classlabels)[1]),size=dim(classlabels)[1]),],oscmode=featselmethod,numcomp=opt_comp,kfold=kfold,evalmethod=pred.eval.method,keepX=max_varsel,sparseselect=TRUE,analysismode,sample.col.vec=col_vec,scoreplot_legend=scoreplot_legend,pairedanalysis=pairedanalysis,optselect=FALSE,class_labels_levels_main=class_labels_levels_main,legendlocation=legendlocation,plotindiv=FALSE)
                    
                    #rand_pls_sel<-cbind(rand_pls_sel,plsresrand$vip_res[,1])
                    #return(plsresrand$vip_res[,1])		
                    if (is(plsresrand, "try-error")){
                      
                      return(rep(0,dim(data_m_fc)[1]))
                    }else{
                      return(plsresrand$vip_res[,1])
                    }
                  })
                  
                  stopCluster(cl)
                  
                }
              }
              
            }
            pls_vip_thresh<-0
            
            if (is(plsres1, "try-error")){
              print(paste("sPLS could not be performed at RSD threshold: ",log2.fold.change.thresh,sep=""))
              break;
            }else{	
              opt_comp<-plsres1$opt_comp
            }
            
          }else{
            #PLS
            if(pairedanalysis==TRUE){
              classlabels_temp<-cbind(classlabels_sub[,2],classlabels)
              plsres1<-do_plsda(X=data_m_fc,Y=classlabels_temp,oscmode=featselmethod,numcomp=pls_ncomp,kfold=kfold,evalmethod=pred.eval.method,
                                keepX=max_varsel,sparseselect=FALSE,analysismode=analysismode,vip.thresh=pls_vip_thresh,sample.col.opt=sample.col.opt,
                                sample.col.vec=col_vec,scoreplot_legend=scoreplot_legend,pairedanalysis=pairedanalysis,optselect=optselect,
                                class_labels_levels_main=class_labels_levels_main,legendlocation=legendlocation,pls.vip.selection=pls.vip.selection,
                                output.device.type=output.device.type,plots.res=plots.res,plots.width=plots.width,plots.height=plots.height,plots.type=plots.type,pls.ellipse=pca.ellipse)
              
              if (is(plsres1, "try-error")){
                print(paste("PLS could not be performed at RSD threshold: ",log2.fold.change.thresh,sep=""))
                break;
              }else{
                opt_comp<-plsres1$opt_comp
              }
              
            }else{
              
              plsres1<-do_plsda(X=data_m_fc,Y=classlabels,oscmode=featselmethod,numcomp=pls_ncomp,kfold=kfold,evalmethod=pred.eval.method,keepX=max_varsel,
                                sparseselect=FALSE,analysismode=analysismode,vip.thresh=pls_vip_thresh,sample.col.opt=sample.col.opt,
                                sample.col.vec=col_vec,scoreplot_legend=scoreplot_legend,pairedanalysis=pairedanalysis,optselect=optselect,
                                class_labels_levels_main=class_labels_levels_main,legendlocation=legendlocation,pls.vip.selection=pls.vip.selection,
                                output.device.type=output.device.type,plots.res=plots.res,plots.width=plots.width,plots.height=plots.height,
                                plots.type=plots.type,pls.ellipse=pca.ellipse)
              
              if (is(plsres1, "try-error")){
                print(paste("PLS could not be performed at RSD threshold: ",log2.fold.change.thresh,sep=""))
                break;
              }else{
                opt_comp<-plsres1$opt_comp
              }
              #for(randindex in 1:100){
              if(is.na(pls.permut.count)==FALSE){
                
                
                set.seed(999)
                seedvec<-runif(pls.permut.count,10,10*pls.permut.count)
                
                
                
                if(pls.permut.count>0){
                  
                  cl <- parallel::makeCluster(getOption("cl.cores", num_nodes))
                  clusterEvalQ(cl,library(plsgenomics))
                  clusterEvalQ(cl,library(dplyr))
                  
                  clusterEvalQ(cl,library(plyr))
                  clusterExport(cl,"pls.lda.cv",envir = .GlobalEnv)
                  clusterExport(cl,"plsda_cv",envir = .GlobalEnv)
                  #clusterExport(cl,"%>%",envir = .GlobalEnv) #%>%
                  clusterExport(cl,"do_plsda_rand",envir = .GlobalEnv)
                  clusterEvalQ(cl,library(mixOmics))
                  clusterEvalQ(cl,library(pls))
                  
                  #here
                  rand_pls_sel<-parLapply(cl,1:pls.permut.count,function(x)
                  {
                    
                    set.seed(seedvec[x])
                    #t1fname<-paste("ranpls",x,".Rda",sep="")
                    ####savelist=ls(),file=t1fname)
                    print(paste("PLSDA permutation number: ",x,sep=""))
                    
                    plsresrand<-do_plsda_rand(X=data_m_fc,Y=classlabels[sample(x=seq(1,dim(classlabels)[1]),size=dim(classlabels)[1]),],oscmode=featselmethod,numcomp=opt_comp,kfold=kfold,evalmethod=pred.eval.method,keepX=max_varsel,sparseselect=FALSE,analysismode,sample.col.vec=col_vec,scoreplot_legend=scoreplot_legend,pairedanalysis=pairedanalysis,optselect=FALSE,class_labels_levels_main=class_labels_levels_main,legendlocation=legendlocation,plotindiv=FALSE) #,silent=TRUE)
                    
                    if (is(plsresrand, "try-error")){
                      
                      
                      return(1)
                    }else{
                      return(plsresrand$vip_res[,1])
                      
                    }
                    
                    
                  })
                  
                  stopCluster(cl)
                }
                ####saverand_pls_sel,file="rand_pls_sel1.Rda")
                
              }
              
            }
            opt_comp<-plsres1$opt_comp
          }
          
          if(length(plsres1$bad_variables)>0){
            
            data_m_fc_withfeats<-data_m_fc_withfeats[-c(plsres1$bad_variables),]
            data_m_fc<-data_m_fc[-c(plsres1$bad_variables),]
          }
          
          
          
          if(is.na(pls.permut.count)==FALSE){
            
            if(pls.permut.count>0){
              
              ###saverand_pls_sel,file="rand_pls_sel.Rda")
              
              #rand_pls_sel<-ldply(rand_pls_sel,rbind)  #unlist(rand_pls_sel)
              rand_pls_sel<-as.data.frame(rand_pls_sel)
              rand_pls_sel<-t(rand_pls_sel)
              rand_pls_sel<-as.data.frame(rand_pls_sel)
              
              if(featselmethod=="spls"){
                
                rand_pls_sel[rand_pls_sel!=0]<-1
              }else{
                
                rand_pls_sel[rand_pls_sel<pls_vip_thresh]<-0
                rand_pls_sel[rand_pls_sel>=pls_vip_thresh]<-1
                
              }
              
              ####saverand_pls_sel,file="rand_pls_sel2.Rda")
              rand_pls_sel_prob<-apply(rand_pls_sel,2,sum)/pls.permut.count
              #rand_pls_sel_fdr<-p.adjust(rand_pls_sel_prob,method=fdrmethod)
              pvalues<-rand_pls_sel_prob
              if(fdrmethod=="BH"){
                fdr_adjust_pvalue<-p.adjust(pvalues,method="BH")
              }else{
                if(fdrmethod=="ST"){
                  #fdr_adjust_pvalue<-qvalue(pvalues)
                  #fdr_adjust_pvalue<-fdr_adjust_pvalue$qvalues
                  
                  fdr_adjust_pvalue<-try(qvalue(pvalues),silent=TRUE)
                  
                  if(is(fdr_adjust_pvalue,"try-error")){
                    
                    fdr_adjust_pvalue<-qvalue(pvalues,lambda=max(pvalues,na.rm=TRUE))
                  }
                  
                  fdr_adjust_pvalue<-fdr_adjust_pvalue$qvalues
                  
                  
                }else{
                  if(fdrmethod=="Strimmer"){
                    pdf("fdrtool.pdf")
                    #par_rows=1
                    #par(mfrow=c(par_rows,1))
                    fdr_adjust_pvalue<-suppressWarnings(fdrtool(as.vector(pvalues),statistic="pvalue",verbose=FALSE))
                    fdr_adjust_pvalue<-fdr_adjust_pvalue$qval
                    try(dev.off(),silent=TRUE)
                  }else{
                    if(fdrmethod=="none"){
                      fdr_adjust_pvalue<-pvalues
                      #fdr_adjust_pvalue<-p.adjust(pvalues,method="none")
                    }else{
                      if(fdrmethod=="BY"){
                        fdr_adjust_pvalue<-p.adjust(pvalues,method="BY")
                      }else{
                        if(fdrmethod=="bonferroni"){
                          fdr_adjust_pvalue<-p.adjust(pvalues,method="bonferroni")
                        }
                      }
                      
                    }
                  }
                }
                
              }
              
              rand_pls_sel_fdr<-fdr_adjust_pvalue
              
              
              vip_res<-cbind(data_m_fc_withfeats[,c(1:2)],plsres1$vip_res,rand_pls_sel_prob,rand_pls_sel_fdr)
              
            }else{
              vip_res<-cbind(data_m_fc_withfeats[,c(1:2)],plsres1$vip_res)
              rand_pls_sel_fdr<-rep(0,dim(data_m_fc_withfeats[,c(1:2)])[1])
              rand_pls_sel_prob<-rep(0,dim(data_m_fc_withfeats[,c(1:2)])[1])
            }
          }else{
            vip_res<-cbind(data_m_fc_withfeats[,c(1:2)],plsres1$vip_res)
            rand_pls_sel_fdr<-rep(0,dim(data_m_fc_withfeats[,c(1:2)])[1])
            rand_pls_sel_prob<-rep(0,dim(data_m_fc_withfeats[,c(1:2)])[1])
          }
          
          write.table(vip_res,file="Tables/vip_res.txt",sep="\t",row.names=FALSE)
          
          
          #				write.table(r2_q2_valid_res,file="pls_r2_q2_res.txt",sep="\t",row.names=TRUE)
          
          
          varimp_plsres1<-plsres1$selected_variables
          
          opt_comp<-plsres1$opt_comp				
          if(max_comp_sel>opt_comp){
            
            max_comp_sel<-opt_comp
          }
          
          #	print("opt comp is")
          #print(opt_comp)
          if(featselmethod=="spls"){
            
            cnames_tab<-colnames(data_m_fc_withfeats)
            cnames_tab<-c("Loading (absolute)","Rank",cnames_tab)
            
            #
            if(opt_comp>1){
              
              #abs
              vip_res1<-abs(plsres1$vip_res)
              
              if(max_comp_sel>1){
                vip_res1<-apply(vip_res1[,c(1:max_comp_sel)],1,max)
                
              }else{
                
                vip_res1<-vip_res1[,c(1)]
              }		
            }else{
              
              vip_res1<-abs(plsres1$vip_res)
            }	
            
            pls_vip<-vip_res1 #(plsres1$vip_res)
            
            
            if(is.na(pls.permut.count)==FALSE){
            #based on loadings for sPLS
            sel.diffdrthresh<-pls_vip!=0 & rand_pls_sel_fdr<fdrthresh & rand_pls_sel_prob<pvalue.thresh
            }else{
              
              print("DOING SPLS #here999")
              sel.diffdrthresh<-pls_vip!=0
            }
            
            goodip<-which(sel.diffdrthresh==TRUE)
            
            save(goodip,pls_vip,rand_pls_sel_fdr,rand_pls_sel_prob,sel.diffdrthresh,file="splsdebug1.Rda")
            
            
            
            
          }else{
            
            cnames_tab<-colnames(data_m_fc_withfeats)
            cnames_tab<-c("VIP","Rank",cnames_tab)
            
            
            if(max_comp_sel>opt_comp){
              
              max_comp_sel<-opt_comp
            }
            
            
            #pls_vip<-plsres1$vip_res[,c(1:max_comp_sel)]
            
            
            if(opt_comp>1){
              vip_res1<-(plsres1$vip_res)
              if(max_comp_sel>1){
                
                if(pls.vip.selection=="mean"){
                  vip_res1<-apply(vip_res1[,c(1:max_comp_sel)],1,mean)
                  
                }else{
                  
                  vip_res1<-apply(vip_res1[,c(1:max_comp_sel)],1,max)
                }
              }else{
                
                vip_res1<-vip_res1[,c(1)]
              }
            }else{
              
              vip_res1<-plsres1$vip_res
            }
            
            
            
            #vip_res1<-plsres1$vip_res
            pls_vip<-vip_res1
            
            
            
            
            #pls
            sel.diffdrthresh<-pls_vip>=pls_vip_thresh & rand_pls_sel_fdr<fdrthresh & rand_pls_sel_prob<pvalue.thresh
            
            goodip<-which(sel.diffdrthresh==TRUE)
            
          }
          
          rank_vec<-order(pls_vip,decreasing=TRUE)
          rank_vec2<-seq(1,length(rank_vec))
          
          ranked_vec<-pls_vip[rank_vec]
          rank_num<-match(pls_vip,ranked_vec)
          
          
          
          data_limma_fdrall_withfeats<-cbind(pls_vip,rank_num,data_m_fc_withfeats)
          
          
          feat_sigfdrthresh[lf]<-length(which(sel.diffdrthresh==TRUE)) #length(plsres1$selected_variables) #length(which(sel.diffdrthresh==TRUE))
          
          filename<-paste("Tables/",parentfeatselmethod,"_variable_importance.txt",sep="")
          
          colnames(data_limma_fdrall_withfeats)<-as.character(cnames_tab)
          
          #data_limma_fdrall_withfeats<-data_limma_fdrall_withfeats[order(fdr_adjust_fpvalue),]
          write.table(data_limma_fdrall_withfeats, file=filename,sep="\t",row.names=FALSE)
          
          
          
          
        }
        
        #stop("Please choose limma, RF, RFcond, or MARS for featselmethod.")
        if(featselmethod=="lmreg" | featselmethod=="lm1wayanova" | featselmethod=="lm2wayanova" | featselmethod=="lm1wayanovarepeat" | featselmethod=="lm2wayanovarepeat"| featselmethod=="logitreg" | featselmethod=="wilcox" | featselmethod=="ttest" | 
           featselmethod=="ttestrepeat" |  featselmethod=="poissonreg" | featselmethod=="wilcoxrepeat" | featselmethod=="lmregrepeat")
        {
          pvalues<-{}
          
          classlabels_response_mat<-as.data.frame(classlabels_response_mat)
          
          if(featselmethod=="ttestrepeat"){
            featselmethod="ttest"
            pairedanalysis=TRUE
          }
          
          if(featselmethod=="wilcoxrepeat"){
            
            featselmethod="wilcox"
            pairedanalysis=TRUE
          }
          
          if(featselmethod=="lm1wayanova")
          {
            
            print("Performing one-way ANOVA analysis")
            
            #print(dim(data_m_fc))
            #print(dim(classlabels_response_mat))
            #print(dim(classlabels))
            
            #data_mat_anova<-cbind(t(data_m_fc),classlabels_response_mat)
            
            
            numcores<-round(detectCores()*0.6)
            
            cl <- parallel::makeCluster(getOption("cl.cores", num_nodes))
            
            clusterExport(cl,"diffexponewayanova",envir = .GlobalEnv)
            
            clusterExport(cl,"anova",envir = .GlobalEnv)
            
            
            clusterExport(cl,"TukeyHSD",envir = .GlobalEnv)
            
            clusterExport(cl,"aov",envir = .GlobalEnv)
            
            
            
            res1<-parApply(cl,data_m_fc,1,function(x,classlabels_response_mat){
              xvec<-x	
              
              
              data_mat_anova<-cbind(xvec,classlabels_response_mat)
              
              data_mat_anova<-as.data.frame(data_mat_anova)
              cnames<-colnames(data_mat_anova)
              
              cnames[1]<-"Response"
              
              colnames(data_mat_anova)<-c("Response","Factor1")
              
              
              
              data_mat_anova$Factor1<-as.factor(data_mat_anova$Factor1)
              
              anova_res<-diffexponewayanova(dataA=data_mat_anova)
              
              
              
              return(anova_res)
            },classlabels_response_mat)
            
            stopCluster(cl)
            main_pval_mat<-{}
            
            posthoc_pval_mat<-{}
            pvalues<-{}
            
            
            #print(head(res1))
            
            for(i in 1:length(res1)){
              
              
              main_pval_mat<-rbind(main_pval_mat,res1[[i]]$mainpvalues)
              pvalues<-c(pvalues,res1[[i]]$mainpvalues[1])
              posthoc_pval_mat<-rbind(posthoc_pval_mat,res1[[i]]$posthocfactor1)
              
              
            }
            pvalues<-unlist(pvalues)
            
            #print(summary(pvalues))
            
            if(fdrmethod=="BH"){
              fdr_adjust_pvalue<-p.adjust(pvalues,method="BH")
            }else{
              if(fdrmethod=="ST"){
                #fdr_adjust_pvalue<-qvalue(pvalues)
                #fdr_adjust_pvalue<-fdr_adjust_pvalue$qvalues
                
                fdr_adjust_pvalue<-try(qvalue(pvalues),silent=TRUE)
                
                if(is(fdr_adjust_pvalue,"try-error")){
                  
                  fdr_adjust_pvalue<-qvalue(pvalues,lambda=max(pvalues,na.rm=TRUE))
                }
                
                fdr_adjust_pvalue<-fdr_adjust_pvalue$qvalues
                
                
              }else{
                if(fdrmethod=="Strimmer"){
                  pdf("fdrtool.pdf")
                  #par_rows=1
                  #par(mfrow=c(par_rows,1))
                  fdr_adjust_pvalue<-suppressWarnings(fdrtool(as.vector(pvalues),statistic="pvalue",verbose=FALSE))
                  fdr_adjust_pvalue<-fdr_adjust_pvalue$qval
                  try(dev.off(),silent=TRUE)
                }else{
                  if(fdrmethod=="none"){
                    #fdr_adjust_pvalue<-pvalues
                    fdr_adjust_pvalue<-p.adjust(pvalues,method="none")
                  }else{
                    if(fdrmethod=="BY"){
                      fdr_adjust_pvalue<-p.adjust(pvalues,method="BY")
                    }else{
                      if(fdrmethod=="bonferroni"){
                        fdr_adjust_pvalue<-p.adjust(pvalues,method="bonferroni")
                      }
                    }
                  }
                }
              }
              
            }
            
            if(fdrmethod=="none"){
              filename<-"lm1wayanova_pvalall_posthoc.txt"
            }else{
              filename<-"lm1wayanova_fdrall_posthoc.txt"
            }
            cnames_tab<-colnames(data_m_fc_withfeats)
            
            
            posthoc_names<-colnames(posthoc_pval_mat)
            if(length(posthoc_names)<1){
              
              posthoc_names<-c("Factor1vs2")
            }																
            
            cnames_tab<-c("P.value","adjusted.P.value",posthoc_names,cnames_tab)
            
            #cnames_tab<-c("P.value","adjusted.P.value","posthoc.pvalue",cnames_tab)
            
            pvalues<-as.data.frame(pvalues)
            
            #pvalues<-t(pvalues)
            
            pvalues<-as.data.frame(pvalues)
            
            final.pvalues<-pvalues
            #final.pvalues<-pvalues
            
            data_limma_fdrall_withfeats<-cbind(pvalues,fdr_adjust_pvalue,posthoc_pval_mat,data_m_fc_withfeats)
            
            colnames(data_limma_fdrall_withfeats)<-as.character(cnames_tab)
            
            #gohere
            
            if(length(check_names)>0){
              
             # data_limma_fdrall_withfeats<-cbind(pvalues1,fdr_adjust_pvalue1,pvalues2,fdr_adjust_pvalue2,pvalues3,fdr_adjust_pvalue3,posthoc_pval_mat,data_m_fc_with_names,data_m_fc_withfeats[,-c(1:2)])
              
              #colnames(data_limma_fdrall_withfeats)<-as.character(cnames_tab)
              
              data_limma_fdrall_withfeats<-cbind(pvalues,fdr_adjust_pvalue,posthoc_pval_mat,data_m_fc_with_names,data_m_fc_withfeats[,-c(1:2)])
              
              colnames(data_limma_fdrall_withfeats)<-as.character(cnames_tab)
              
              
              data_limma_fdrall_withfeats<-as.data.frame(data_limma_fdrall_withfeats)
              #data_limma_fdrall_withfeats<-cbind(p.value,adjusted.p.value,results2,data_m_fc_with_names,data_m_fc_withfeats[,-c(1:2)])
              
              rem_col_ind1<-grep(colnames(data_limma_fdrall_withfeats),pattern=c("mz"))
              
              rem_col_ind2<-grep(colnames(data_limma_fdrall_withfeats),pattern=c("time"))
              
              rem_col_ind<-c(rem_col_ind1,rem_col_ind2)
              
            }else{
              rem_col_ind<-{}
            }
            
            #data_limma_fdrall_withfeats<-data_limma_fdrall_withfeats[order(fdr_adjust_fpvalue),]
            filename<-paste("Tables/",filename,sep="")
            
            
            
            if(length(rem_col_ind)>0){
              #write.table(data_limma_fdrall_withfeats[,-c(rem_col_ind)], file="Tables/twowayanova_with_posthoc_comparisons.txt",sep="\t",row.names=FALSE)
              
              write.table(data_limma_fdrall_withfeats[,-c(rem_col_ind)], file=filename,sep="\t",row.names=FALSE)
            }else{
              
              
              #write.table(data_limma_fdrall_withfeats,file="Tables/twowayanova_with_posthoc_comparisons.txt",sep="\t",row.names=FALSE)
              
              write.table(data_limma_fdrall_withfeats, file=filename,sep="\t",row.names=FALSE)
              
            }
            
           
            
            data_limma_fdrall_withfeats<-cbind(pvalues,fdr_adjust_pvalue,data_m_fc_withfeats)
            
          }
          
          
          
          if(featselmethod=="ttest" && pairedanalysis==TRUE)
          {
            
            print("Performing paired t-test analysis")
            
            #print(dim(data_m_fc))
            #print(dim(classlabels_response_mat))
            #print(dim(classlabels))
            
            #data_mat_anova<-cbind(t(data_m_fc),classlabels_response_mat)
            
            numcores<-round(detectCores()*0.5)
            
            cl <- parallel::makeCluster(getOption("cl.cores", num_nodes))
            
            clusterExport(cl,"t.test",envir = .GlobalEnv)
            
            
            res1<-parApply(cl,data_m_fc,1,function(x,classlabels_response_mat){
              
              xvec<-x
              
              
              data_mat_anova<-cbind(xvec,classlabels_response_mat)
              
              data_mat_anova<-as.data.frame(data_mat_anova)
              cnames<-colnames(data_mat_anova)
              
              cnames[1]<-"Response"
              
              colnames(data_mat_anova)<-c("Response","Factor1")
              
              #print(data_mat_anova)
              
              data_mat_anova$Factor1<-as.factor(data_mat_anova$Factor1)
              
              #anova_res<-diffexponewayanova(dataA=data_mat_anova)
              
              x1<-data_mat_anova$Response[which(data_mat_anova$Factor1==class_labels_levels[1])]
              
              y1<-data_mat_anova$Response[which(data_mat_anova$Factor1==class_labels_levels[2])]
              
              w1<-t.test(x=x1,y=y1,alternative="two.sided",paired=TRUE)
              
              return(w1$p.value)
            },classlabels_response_mat)
            
            stopCluster(cl)
            
            main_pval_mat<-{}
            
            posthoc_pval_mat<-{}
            pvalues<-{}
            
            
            
            pvalues<-unlist(res1)
            
            
            #print(summary(pvalues))
            
            if(fdrmethod=="BH"){
              fdr_adjust_pvalue<-p.adjust(pvalues,method="BH")
            }else{
              if(fdrmethod=="ST"){
                #fdr_adjust_pvalue<-qvalue(pvalues)
                #fdr_adjust_pvalue<-fdr_adjust_pvalue$qvalues
                
                fdr_adjust_pvalue<-try(qvalue(pvalues),silent=TRUE)
                
                if(is(fdr_adjust_pvalue,"try-error")){
                  
                  fdr_adjust_pvalue<-qvalue(pvalues,lambda=max(pvalues,na.rm=TRUE))
                }
                
                fdr_adjust_pvalue<-fdr_adjust_pvalue$qvalues
                
                
              }else{
                if(fdrmethod=="Strimmer"){
                  pdf("fdrtool.pdf")
                  #par_rows=1
                  #par(mfrow=c(par_rows,1))
                  fdr_adjust_pvalue<-suppressWarnings(fdrtool(as.vector(pvalues),statistic="pvalue",verbose=FALSE))
                  fdr_adjust_pvalue<-fdr_adjust_pvalue$qval
                  try(dev.off(),silent=TRUE)
                }else{
                  if(fdrmethod=="none"){
                    #fdr_adjust_pvalue<-pvalues
                    fdr_adjust_pvalue<-p.adjust(pvalues,method="none")
                  }else{
                    if(fdrmethod=="BY"){
                      fdr_adjust_pvalue<-p.adjust(pvalues,method="BY")
                    }else{
                      if(fdrmethod=="bonferroni"){
                        fdr_adjust_pvalue<-p.adjust(pvalues,method="bonferroni")
                      }
                    }
                  }
                }
              }
              
            }
            
            if(fdrmethod=="none"){
              filename<-"pairedttest_pvalall_withfeats.txt"
            }else{
              filename<-"pairedttest_fdrall_withfeats.txt"
            }
            cnames_tab<-colnames(data_m_fc_withfeats)
            
            
            posthoc_names<-colnames(posthoc_pval_mat)
            if(length(posthoc_names)<1){
              
              posthoc_names<-c("Factor1vs2")
            }
            
            cnames_tab<-c("P.value","adjusted.P.value",cnames_tab)
            
            #cnames_tab<-c("P.value","adjusted.P.value","posthoc.pvalue",cnames_tab)
            
            pvalues<-as.data.frame(pvalues)
            
            #pvalues<-t(pvalues)
            # print(dim(pvalues))
            #print(dim(data_m_fc_withfeats))
            
            final.pvalues<-pvalues
            sel.diffdrthresh<-fdr_adjust_pvalue<fdrthresh & final.pvalues<pvalue.thresh
            
            
            data_limma_fdrall_withfeats<-cbind(pvalues,fdr_adjust_pvalue,data_m_fc_withfeats)
            
            colnames(data_limma_fdrall_withfeats)<-as.character(cnames_tab)
            
            #data_limma_fdrall_withfeats<-data_limma_fdrall_withfeats[order(fdr_adjust_fpvalue),]
            #  write.table(data_limma_fdrall_withfeats, file=filename,sep="\t",row.names=FALSE)
            
            
            data_limma_fdrall_withfeats<-cbind(pvalues,fdr_adjust_pvalue,data_m_fc_withfeats)
            
          }
          
          if(featselmethod=="ttest" && pairedanalysis==FALSE)
          {
            
            print("Performing t-test analysis")
            
            #print(dim(data_m_fc))
            #print(dim(classlabels_response_mat))
            #print(dim(classlabels))
            
            #data_mat_anova<-cbind(t(data_m_fc),classlabels_response_mat)
            
            numcores<-round(detectCores()*0.5)
            
            cl <- parallel::makeCluster(getOption("cl.cores", num_nodes))
            
            clusterExport(cl,"t.test",envir = .GlobalEnv)
            
            
            res1<-parApply(cl,data_m_fc,1,function(x,classlabels_response_mat){
              
              xvec<-x
              
              
              data_mat_anova<-cbind(xvec,classlabels_response_mat)
              
              data_mat_anova<-as.data.frame(data_mat_anova)
              cnames<-colnames(data_mat_anova)
              
              cnames[1]<-"Response"
              
              colnames(data_mat_anova)<-c("Response","Factor1")
              
              #print(data_mat_anova)
              
              data_mat_anova$Factor1<-as.factor(data_mat_anova$Factor1)
              
              #anova_res<-diffexponewayanova(dataA=data_mat_anova)
              
              x1<-data_mat_anova$Response[which(data_mat_anova$Factor1==class_labels_levels[1])]
              
              y1<-data_mat_anova$Response[which(data_mat_anova$Factor1==class_labels_levels[2])]
              
              w1<-t.test(x=x1,y=y1,alternative="two.sided")
              
              return(w1$p.value)
            },classlabels_response_mat)
            
            stopCluster(cl)
            
            main_pval_mat<-{}
            
            posthoc_pval_mat<-{}
            pvalues<-{}
            
            
            
            pvalues<-unlist(res1)
            
            
            #print(summary(pvalues))
            
            if(fdrmethod=="BH"){
              fdr_adjust_pvalue<-p.adjust(pvalues,method="BH")
            }else{
              if(fdrmethod=="ST"){
                #fdr_adjust_pvalue<-qvalue(pvalues)
                #fdr_adjust_pvalue<-fdr_adjust_pvalue$qvalues
                
                fdr_adjust_pvalue<-try(qvalue(pvalues),silent=TRUE)
                
                if(is(fdr_adjust_pvalue,"try-error")){
                  
                  fdr_adjust_pvalue<-qvalue(pvalues,lambda=max(pvalues,na.rm=TRUE))
                }
                
                fdr_adjust_pvalue<-fdr_adjust_pvalue$qvalues
                
                
              }else{
                if(fdrmethod=="Strimmer"){
                  pdf("fdrtool.pdf")
                  #par_rows=1
                  #par(mfrow=c(par_rows,1))
                  fdr_adjust_pvalue<-suppressWarnings(fdrtool(as.vector(pvalues),statistic="pvalue",verbose=FALSE))
                  fdr_adjust_pvalue<-fdr_adjust_pvalue$qval
                  try(dev.off(),silent=TRUE)
                }else{
                  if(fdrmethod=="none"){
                    #fdr_adjust_pvalue<-pvalues
                    fdr_adjust_pvalue<-p.adjust(pvalues,method="none")
                  }else{
                    if(fdrmethod=="BY"){
                      fdr_adjust_pvalue<-p.adjust(pvalues,method="BY")
                    }else{
                      if(fdrmethod=="bonferroni"){
                        fdr_adjust_pvalue<-p.adjust(pvalues,method="bonferroni")
                      }
                    }
                  }
                }
              }
              
            }
            
            if(fdrmethod=="none"){
              filename<-"ttest_pvalall_withfeats.txt"
            }else{
              filename<-"ttest_fdrall_withfeats.txt"
            }
            cnames_tab<-colnames(data_m_fc_withfeats)
            
            
            posthoc_names<-colnames(posthoc_pval_mat)
            if(length(posthoc_names)<1){
              
              posthoc_names<-c("Factor1vs2")
            }
            
            cnames_tab<-c("P.value","adjusted.P.value",cnames_tab)
            
            #cnames_tab<-c("P.value","adjusted.P.value","posthoc.pvalue",cnames_tab)
            
            pvalues<-as.data.frame(pvalues)
            
            #pvalues<-t(pvalues)
            # print(dim(pvalues))
            #print(dim(data_m_fc_withfeats))
            
            final.pvalues<-pvalues
            sel.diffdrthresh<-fdr_adjust_pvalue<fdrthresh & final.pvalues<pvalue.thresh
            
            
            data_limma_fdrall_withfeats<-cbind(pvalues,fdr_adjust_pvalue,data_m_fc_withfeats)
            
            colnames(data_limma_fdrall_withfeats)<-as.character(cnames_tab)
            
            #data_limma_fdrall_withfeats<-data_limma_fdrall_withfeats[order(fdr_adjust_fpvalue),]
            #  write.table(data_limma_fdrall_withfeats, file=filename,sep="\t",row.names=FALSE)
            
            
            data_limma_fdrall_withfeats<-cbind(pvalues,fdr_adjust_pvalue,data_m_fc_withfeats)
            
          }
          
          
          if(featselmethod=="wilcox")
          {
            
            print("Performing Wilcox rank-sum analysis")
            
            #print(dim(data_m_fc))
            #print(dim(classlabels_response_mat))
            #print(dim(classlabels))
            
            #data_mat_anova<-cbind(t(data_m_fc),classlabels_response_mat)
            
            numcores<-round(detectCores()*0.5)
            
            cl <- parallel::makeCluster(getOption("cl.cores", num_nodes))
            
            clusterExport(cl,"wilcox.test",envir = .GlobalEnv)
            
            
            res1<-parApply(cl,data_m_fc,1,function(x,classlabels_response_mat){
              
              xvec<-x
              
              
              data_mat_anova<-cbind(xvec,classlabels_response_mat)
              
              data_mat_anova<-as.data.frame(data_mat_anova)
              cnames<-colnames(data_mat_anova)
              
              cnames[1]<-"Response"
              
              colnames(data_mat_anova)<-c("Response","Factor1")
              
              #print(data_mat_anova)
              
              data_mat_anova$Factor1<-as.factor(data_mat_anova$Factor1)
              
              #anova_res<-diffexponewayanova(dataA=data_mat_anova)
              
              x1<-data_mat_anova$Response[which(data_mat_anova$Factor1==class_labels_levels[1])]
              
              y1<-data_mat_anova$Response[which(data_mat_anova$Factor1==class_labels_levels[2])]
              
              w1<-wilcox.test(x=x1,y=y1,alternative="two.sided")
              
              return(w1$p.value)
            },classlabels_response_mat)
            
            stopCluster(cl)
            
            main_pval_mat<-{}
            
            posthoc_pval_mat<-{}
            pvalues<-{}
            
            
            
            pvalues<-unlist(res1)
            
            
            #print(summary(pvalues))
            
            if(fdrmethod=="BH"){
              fdr_adjust_pvalue<-p.adjust(pvalues,method="BH")
            }else{
              if(fdrmethod=="ST"){
                #fdr_adjust_pvalue<-qvalue(pvalues)
                #fdr_adjust_pvalue<-fdr_adjust_pvalue$qvalues
                
                fdr_adjust_pvalue<-try(qvalue(pvalues),silent=TRUE)
                
                if(is(fdr_adjust_pvalue,"try-error")){
                  
                  fdr_adjust_pvalue<-qvalue(pvalues,lambda=max(pvalues,na.rm=TRUE))
                }
                
                fdr_adjust_pvalue<-fdr_adjust_pvalue$qvalues
                
                
              }else{
                if(fdrmethod=="Strimmer"){
                  pdf("fdrtool.pdf")
                  #par_rows=1
                  #par(mfrow=c(par_rows,1))
                  fdr_adjust_pvalue<-suppressWarnings(fdrtool(as.vector(pvalues),statistic="pvalue",verbose=FALSE))
                  fdr_adjust_pvalue<-fdr_adjust_pvalue$qval
                  try(dev.off(),silent=TRUE)
                }else{
                  if(fdrmethod=="none"){
                    #fdr_adjust_pvalue<-pvalues
                    fdr_adjust_pvalue<-p.adjust(pvalues,method="none")
                  }else{
                    if(fdrmethod=="BY"){
                      fdr_adjust_pvalue<-p.adjust(pvalues,method="BY")
                    }else{
                      if(fdrmethod=="bonferroni"){
                        fdr_adjust_pvalue<-p.adjust(pvalues,method="bonferroni")
                      }
                    }
                  }
                }
              }
              
            }
            
            if(fdrmethod=="none"){
              filename<-"wilcox_pvalall_withfeats.txt"
            }else{
              filename<-"wilcox_fdrall_withfeats.txt"
            }
            cnames_tab<-colnames(data_m_fc_withfeats)
            
            
            posthoc_names<-colnames(posthoc_pval_mat)
            if(length(posthoc_names)<1){
              
              posthoc_names<-c("Factor1vs2")
            }																
            
            cnames_tab<-c("P.value","adjusted.P.value",cnames_tab)
            
            #cnames_tab<-c("P.value","adjusted.P.value","posthoc.pvalue",cnames_tab)
            
            pvalues<-as.data.frame(pvalues)
            
            
            final.pvalues<-pvalues
            sel.diffdrthresh<-fdr_adjust_pvalue<fdrthresh & final.pvalues<pvalue.thresh
            
            
            data_limma_fdrall_withfeats<-cbind(pvalues,fdr_adjust_pvalue,data_m_fc_withfeats)
            
            colnames(data_limma_fdrall_withfeats)<-as.character(cnames_tab)
            
            #data_limma_fdrall_withfeats<-data_limma_fdrall_withfeats[order(fdr_adjust_fpvalue),]
            #  write.table(data_limma_fdrall_withfeats, file=filename,sep="\t",row.names=FALSE)
            
            
            data_limma_fdrall_withfeats<-cbind(pvalues,fdr_adjust_pvalue,data_m_fc_withfeats)
            
          }
          
          #lmreg:feature selections
          if(featselmethod=="lmreg")
          {
            
            if(logistic_reg==TRUE){
              
              
              if(length(levels(classlabels_response_mat[,1]))>2){
                
                print("More than 2 classes found. Skipping logistic regression analysis.")
                next;
              }
              
              print("Performing logistic regression analysis:")
              
              classlabels_response_mat[,1]<-as.numeric((classlabels_response_mat[,1]))-1
              
              fileheader="logitreg"
              
              
            }else{
              
              if(poisson_reg==TRUE){
                
                
                print("Performing poisson regression analysis")
                fileheader="poissonreg"
                classlabels_response_mat[,1]<-as.numeric((classlabels_response_mat[,1]))
                
              }else{
                print("Performing linear regression analysis:")
                fileheader="lmreg"
              }
            }
            
            
            
            numcores<-num_nodes #round(detectCores()*0.5)
            
            cl <- parallel::makeCluster(getOption("cl.cores", num_nodes))
            
            clusterExport(cl,"diffexplmreg",envir = .GlobalEnv)
            clusterExport(cl,"lm",envir = .GlobalEnv)
            clusterExport(cl,"glm",envir = .GlobalEnv)
            clusterExport(cl,"summary",envir = .GlobalEnv)
            clusterExport(cl,"anova",envir = .GlobalEnv)
            clusterEvalQ(cl,library(sandwich))
            
            
            #data_mat_anova<-cbind(t(data_m_fc),classlabels_response_mat)
            res1<-parApply(cl,data_m_fc,1,function(x,classlabels_response_mat,logistic_reg,poisson_reg,robust.estimate,vcovHC.type){
              
              xvec<-x
              
              
              
              data_mat_anova<-cbind(xvec,classlabels_response_mat)
              
              cnames<-colnames(data_mat_anova)
              cnames[1]<-"Response"
              
              colnames(data_mat_anova)<-cnames
              
              #lmreg feature selection
              anova_res<-diffexplmreg(dataA=data_mat_anova,logistic_reg,poisson_reg,robust.estimate,vcovHC.type)
              
              return(anova_res)
            },classlabels_response_mat,logistic_reg,poisson_reg,robust.estimate,vcovHC.type)
            
            stopCluster(cl)
            main_pval_mat<-{}
            
            posthoc_pval_mat<-{}
            pvalues<-{}
            
            #save(res1,file="res1.Rda")
            
            all_inf_mat<-{}
            
            for(i in 1:length(res1)){
              
              main_pval_mat<-rbind(main_pval_mat,res1[[i]]$mainpvalues)
              pvalues<-c(pvalues,res1[[i]]$mainpvalues[1])
              #posthoc_pval_mat<-rbind(posthoc_pval_mat,res1[[i]]$posthocfactor1)
              
              cur_pvals<-t(res1[[i]]$mainpvalues)
              cur_est<-t(res1[[i]]$estimates)
              cur_stderr<-t(res1[[i]]$stderr)
              cur_tstat<-t(res1[[i]]$statistic)
              
              #cur_pvals<-as.data.frame(cur_pvals)
              
              cur_res<-cbind(cur_pvals,cur_est,cur_stderr,cur_tstat)
              
              all_inf_mat<-rbind(all_inf_mat,cur_res)
              
              
            }
            
            
            
            cnames_1<-c(paste("P.value_",colnames(cur_pvals),sep=""),paste("Estimate_",colnames(cur_pvals),sep=""),paste("StdError_var_",colnames(cur_pvals),sep=""),paste("t-statistic_",colnames(cur_pvals),sep=""))
            
            
            
            #	print("here after lm reg")
            
            #print(summary(pvalues))
            
            if(fdrmethod=="BH"){
              fdr_adjust_pvalue<-p.adjust(pvalues,method="BH")
            }else{
              if(fdrmethod=="ST"){
                #fdr_adjust_pvalue<-qvalue(pvalues)
                #fdr_adjust_pvalue<-fdr_adjust_pvalue$qvalues
                
                fdr_adjust_pvalue<-try(qvalue(pvalues),silent=TRUE)
                
                if(is(fdr_adjust_pvalue,"try-error")){
                  
                  fdr_adjust_pvalue<-qvalue(pvalues,lambda=max(pvalues,na.rm=TRUE))
                }
                
                fdr_adjust_pvalue<-fdr_adjust_pvalue$qvalues
                
                
              }else{
                if(fdrmethod=="Strimmer"){
                  pdf("fdrtool.pdf")
                  #par_rows=1
                  #par(mfrow=c(par_rows,1))
                  fdr_adjust_pvalue<-suppressWarnings(fdrtool(as.vector(pvalues),statistic="pvalue",verbose=FALSE))
                  fdr_adjust_pvalue<-fdr_adjust_pvalue$qval
                  try(dev.off(),silent=TRUE)
                }else{
                  if(fdrmethod=="none"){
                    #fdr_adjust_pvalue<-pvalues
                    fdr_adjust_pvalue<-p.adjust(pvalues,method="none")
                  }else{
                    if(fdrmethod=="BY"){
                      fdr_adjust_pvalue<-p.adjust(pvalues,method="BY")
                    }else{
                      if(fdrmethod=="bonferroni"){
                        fdr_adjust_pvalue<-p.adjust(pvalues,method="bonferroni")
                      }
                    }
                  }
                }
              }
              
              
              
            }
            
            
            if(fdrmethod=="none"){
              filename<-paste(fileheader,"_pvalall_withfeats.txt",sep="")
              
            }else{
              filename<-paste(fileheader,"_fdrall_withfeats.txt",sep="")
            }
            cnames_tab<-colnames(data_m_fc_withfeats)
            cnames_tab<-c("P.value","adjusted.P.value",cnames_tab)
            
            pvalues<-as.data.frame(pvalues)
            
            
            
            
            final.pvalues<-pvalues
            sel.diffdrthresh<-fdr_adjust_pvalue<fdrthresh & final.pvalues<pvalue.thresh
            
            
            data_limma_fdrall_withfeats<-cbind(pvalues,fdr_adjust_pvalue,data_m_fc_withfeats)
            
            colnames(data_limma_fdrall_withfeats)<-as.character(cnames_tab)
            
            #data_limma_fdrall_withfeats<-data_limma_fdrall_withfeats[order(fdr_adjust_fpvalue),]
            #write.table(data_limma_fdrall_withfeats, file=filename,sep="\t",row.names=FALSE)
            
            
            if(analysismode=="regression"){
              
              filename<-paste(fileheader,"_results_allfeatures.txt",sep="")
              filename<-paste("Tables/",filename,sep="")
              
              #                  write.table(data_limma_fdrall_withfeats, file=filename,sep="\t",row.names=FALSE)
            }
            
            filename<-paste(fileheader,"_pval_coef_stderr.txt",sep="")
            
            
            
            data_allinf_withfeats<-cbind(all_inf_mat,data_m_fc_withfeats)
            filename<-paste("Tables/",filename,sep="")
            #     write.table(data_allinf_withfeats, file=filename,sep="\t",row.names=FALSE)
            
            
            cnames_tab<-colnames(data_m_fc_withfeats)
            
            
            cnames_tab<-c(cnames_1,cnames_tab)
            
            class_column_names<-colnames(classlabels_response_mat)
            
            
            colnames(data_allinf_withfeats)<-as.character(cnames_tab)
            
            ###save(data_allinf_withfeats,cnames_tab,cnames_1,file="data_allinf_withfeats.Rda")
            
            
            pval_columns<-grep(colnames(data_allinf_withfeats),pattern="P.value")
            
            fdr_adjusted_pvalue<-get_fdr_adjusted_pvalue(data_matrix=data_allinf_withfeats,fdrmethod=fdrmethod)
            
            #   data_allinf_withfeats1<-cbind(data_allinf_withfeats[,pval_columns],fdr_adjusted_pvalue,data_allinf_withfeats[,-c(pval_columns)])
            
            cnames_tab1<-c(cnames_tab[pval_columns],colnames(fdr_adjusted_pvalue),cnames_tab[-pval_columns])
            pval_columns<-grep(colnames(data_allinf_withfeats),pattern="P.value")
            
            fdr_adjusted_pvalue<-get_fdr_adjusted_pvalue(data_matrix=data_allinf_withfeats,fdrmethod=fdrmethod)
            
            data_allinf_withfeats<-cbind(data_allinf_withfeats[,pval_columns],fdr_adjusted_pvalue,data_allinf_withfeats[,-c(pval_columns)])
            
            cnames_tab1<-c(cnames_tab[pval_columns],colnames(fdr_adjusted_pvalue),cnames_tab[-pval_columns])
            
            
            filename<-paste(fileheader,"_pval_coef_stderr.txt",sep="")
            filename<-paste("Tables/",filename,sep="")
            colnames(data_allinf_withfeats)<-cnames_tab1
            
            ###save(data_allinf_withfeats,file="d2.Rda")
            write.table(data_allinf_withfeats, file=filename,sep="\t",row.names=FALSE)
            
            
            
            
            
            
          }
          
          
          
          if(featselmethod=="lm2wayanova")
          {
            
            print("Performing two-way ANOVA analysis with Tukey post hoc comparisons")
            
            #print(dim(data_m_fc))
            #			print(dim(classlabels_response_mat))
            
           
            numcores<-num_nodes #round(detectCores()*0.5)
            
            cl <- parallel::makeCluster(getOption("cl.cores", num_nodes))
            
            
            clusterExport(cl,"diffexplmtwowayanova",envir = .GlobalEnv)
            
            clusterExport(cl,"TukeyHSD",envir = .GlobalEnv)
            clusterExport(cl,"plotTukeyHSD1",envir = .GlobalEnv)
            
            clusterExport(cl,"aov",envir = .GlobalEnv)
            clusterExport(cl,"anova",envir = .GlobalEnv)
            clusterEvalQ(cl,library(ggpubr))
            clusterEvalQ(cl,library(ggplot2))
           # clusterEvalQ(cl,library(cowplot))
            
            
            #res1<-apply(data_m_fc,1,function(x){
            
            res1<-parRapply(cl,data_m_fc,function(x,classlabels_response_mat){
              
              xvec<-x
              
              colnames(classlabels_response_mat)<-paste("Factor",seq(1,dim(classlabels_response_mat)[2]),sep="")
              
              data_mat_anova<-cbind(xvec,classlabels_response_mat)
              #print("2way anova")
              #	print(data_mat_anova[1:2,])
              cnames<-colnames(data_mat_anova)
              cnames[1]<-"Response"
              
              colnames(data_mat_anova)<-cnames
              
              ####savedata_mat_anova,file="data_mat_anova.Rda")
              
              #diffexplmtwowayanova
              anova_res<-diffexplmtwowayanova(dataA=data_mat_anova)
              
              
              return(anova_res)
            },classlabels_response_mat)
            
            
            stopCluster(cl)
            #	print("done")
            
            ####saveres1,file="res1.Rda")
            
            main_pval_mat<-{}
            posthoc_pval_mat<-{}
            pvalues1<-{}
            pvalues2<-{}
            pvalues3<-{}
            
            save(res1,file="tukeyhsd_plots.Rda")
            
            for(i in 1:length(res1)){
              
              #print(i)
              #print(res1[[i]]$mainpvalues)
              #print(res1[[i]]$posthoc)
              main_pval_mat<-rbind(main_pval_mat,res1[[i]]$mainpvalues)
              pvalues1<-c(pvalues1,as.vector(res1[[i]]$mainpvalues[1]))
              pvalues2<-c(pvalues2,as.vector(res1[[i]]$mainpvalues[2]))
              pvalues3<-c(pvalues3,as.vector(res1[[i]]$mainpvalues[3]))
              posthoc_pval_mat<-rbind(posthoc_pval_mat,res1[[i]]$posthoc)
              
              
            }
            twoanova_res<-cbind(data_m_fc_withfeats[,c(1:2)],main_pval_mat,posthoc_pval_mat)
            
            #write.table(twoanova_res,file="Tables/twoanova_with_posthoc_pvalues.txt",sep="\t",row.names=FALSE)
            pvalues1<-main_pval_mat[,1]
            pvalues2<-main_pval_mat[,2]
            pvalues3<-main_pval_mat[,3]
            
            if(fdrmethod=="none"){
              fdr_adjust_pvalue1<-p.adjust(pvalues1,method="none")
              fdr_adjust_pvalue2<-p.adjust(pvalues2,method="none")
              fdr_adjust_pvalue3<-p.adjust(pvalues3,method="none")
            }
            
            if(fdrmethod=="BH"){
              fdr_adjust_pvalue1<-p.adjust(pvalues1,method="BH")
              fdr_adjust_pvalue2<-p.adjust(pvalues2,method="BH")
              fdr_adjust_pvalue3<-p.adjust(pvalues3,method="BH")
            }else{
              if(fdrmethod=="ST"){
                
                
                fdr_adjust_pvalue1<-try(qvalue(pvalues1),silent=TRUE)
                
                if(is(fdr_adjust_pvalue1,"try-error")){
                  
                  fdr_adjust_pvalue1<-qvalue(pvalues1,lambda=max(pvalues1,na.rm=TRUE))
                }
                
                fdr_adjust_pvalue1<-fdr_adjust_pvalue1$qvalues
                
                fdr_adjust_pvalue2<-try(qvalue(pvalues2),silent=TRUE)
                
                if(is(fdr_adjust_pvalue2,"try-error")){
                  
                  fdr_adjust_pvalue2<-qvalue(pvalues2,lambda=max(pvalues2,na.rm=TRUE))
                }
                
                fdr_adjust_pvalue2<-fdr_adjust_pvalue2$qvalues
                
                fdr_adjust_pvalue3<-try(qvalue(pvalues3),silent=TRUE)
                
                if(is(fdr_adjust_pvalue3,"try-error")){
                  
                  fdr_adjust_pvalue3<-qvalue(pvalues3,lambda=max(pvalues3,na.rm=TRUE))
                }
                
                fdr_adjust_pvalue3<-fdr_adjust_pvalue3$qvalues
                
                
              }else{
                if(fdrmethod=="Strimmer"){
                  pdf("fdrtool.pdf")
                  #par_rows=1
                  #par(mfrow=c(par_rows,1))
                  fdr_adjust_pvalue1<-fdrtool(as.vector(pvalues1),statistic="pvalue",verbose=FALSE)
                  fdr_adjust_pvalue1<-fdr_adjust_pvalue1$qval
                  
                  fdr_adjust_pvalue2<-fdrtool(as.vector(pvalues2),statistic="pvalue",verbose=FALSE)
                  fdr_adjust_pvalue2<-fdr_adjust_pvalue2$qval
                  
                  fdr_adjust_pvalue3<-fdrtool(as.vector(pvalues3),statistic="pvalue",verbose=FALSE)
                  fdr_adjust_pvalue3<-fdr_adjust_pvalue3$qval
                  try(dev.off(),silent=TRUE)
                }else{
                  if(fdrmethod=="none"){
                    fdr_adjust_pvalue1<-p.adjust(pvalues1,method="none")
                    fdr_adjust_pvalue2<-p.adjust(pvalues2,method="none")
                    fdr_adjust_pvalue3<-p.adjust(pvalues3,method="none")
                    
                    
                  }else{
                    if(fdrmethod=="BY"){
                      fdr_adjust_pvalue1<-p.adjust(pvalues1,method="BY")
                      fdr_adjust_pvalue2<-p.adjust(pvalues2,method="BY")
                      fdr_adjust_pvalue3<-p.adjust(pvalues3,method="BY")
                      
                      
                    }else{
                      if(fdrmethod=="bonferroni"){
                        # fdr_adjust_pvalue<-p.adjust(pvalues,method="bonferroni")
                        fdr_adjust_pvalue1<-p.adjust(pvalues1,method="bonferroni")
                        fdr_adjust_pvalue2<-p.adjust(pvalues2,method="bonferroni")
                        fdr_adjust_pvalue3<-p.adjust(pvalues3,method="bonferroni")
                      }
                    }
                    
                  }
                }
              }
              
              
              
            }
            
            if(fdrmethod=="none"){
              filename<-paste(featselmethod,"_pvalall_withfeats.txt",sep="")
            }else{
              filename<-paste(featselmethod,"_fdrall_withfeats.txt",sep="")
            }
            cnames_tab<-colnames(data_m_fc_withfeats)
            
            posthoc_names<-colnames(posthoc_pval_mat)
            
            cnames_tab<-c("Factor1.P.value","Factor1.adjusted.P.value","Factor2.P.value","Factor2.adjusted.P.value","Interact.P.value","Interact.adjusted.P.value",posthoc_names,cnames_tab)
            
            if(FALSE)
            {
              pvalues1<-as.data.frame(pvalues1)
              pvalues1<-t(pvalues1)
              fdr_adjust_pvalue1<-as.data.frame(fdr_adjust_pvalue1)
              pvalues2<-as.data.frame(pvalues2)
              pvalues2<-t(pvalues2)
              fdr_adjust_pvalue2<-as.data.frame(fdr_adjust_pvalue2)
              pvalues3<-as.data.frame(pvalues3)
              pvalues3<-t(pvalues3)
              fdr_adjust_pvalue3<-as.data.frame(fdr_adjust_pvalue3)
              posthoc_pval_mat<-as.data.frame(posthoc_pval_mat)
            }
            
            # ###savedata_m_fc_withfeats,file="data_m_fc_withfeats.Rda")
            data_limma_fdrall_withfeats<-cbind(pvalues1,fdr_adjust_pvalue1,pvalues2,fdr_adjust_pvalue2,pvalues3,fdr_adjust_pvalue3,posthoc_pval_mat,data_m_fc_withfeats)
            
            fdr_adjust_pvalue<-cbind(fdr_adjust_pvalue1,fdr_adjust_pvalue2,fdr_adjust_pvalue3)
            fdr_adjust_pvalue<-apply(fdr_adjust_pvalue,1,function(x){min(x,na.rm=TRUE)})
            
            
            colnames(data_limma_fdrall_withfeats)<-as.character(cnames_tab)
            
            #data_limma_fdrall_withfeats<-data_limma_fdrall_withfeats[order(fdr_adjust_fpvalue),]
            
            if(length(check_names)>0){
              
              data_limma_fdrall_withfeats<-cbind(pvalues1,fdr_adjust_pvalue1,pvalues2,fdr_adjust_pvalue2,pvalues3,fdr_adjust_pvalue3,posthoc_pval_mat,data_m_fc_with_names,data_m_fc_withfeats[,-c(1:2)])
              
              colnames(data_limma_fdrall_withfeats)<-as.character(cnames_tab)
              
              data_limma_fdrall_withfeats<-as.data.frame(data_limma_fdrall_withfeats)
              #data_limma_fdrall_withfeats<-cbind(p.value,adjusted.p.value,results2,data_m_fc_with_names,data_m_fc_withfeats[,-c(1:2)])
              
              rem_col_ind1<-grep(colnames(data_limma_fdrall_withfeats),pattern=c("mz"))
              
              rem_col_ind2<-grep(colnames(data_limma_fdrall_withfeats),pattern=c("time"))
              
              rem_col_ind<-c(rem_col_ind1,rem_col_ind2)
              
            }else{
              rem_col_ind<-{}
            }
            
            if(length(rem_col_ind)>0){
              write.table(data_limma_fdrall_withfeats[,-c(rem_col_ind)], file="Tables/twowayanova_with_posthoc_comparisons.txt",sep="\t",row.names=FALSE)
            }else{
              
              
              write.table(data_limma_fdrall_withfeats,file="Tables/twowayanova_with_posthoc_comparisons.txt",sep="\t",row.names=FALSE)
              
            }
            
            
            filename<-paste("Tables/",filename,sep="")
            
            #write.table(data_limma_fdrall_withfeats,file="Tables/twowayanova_with_posthoc_comparisons.txt",sep="\t",row.names=FALSE)
            
            #write.table(data_limma_fdrall_withfeats, file=filename,sep="\t",row.names=FALSE)
            
            
            fdr_matrix<-cbind(fdr_adjust_pvalue1,fdr_adjust_pvalue2,fdr_adjust_pvalue3)
            
            fdr_matrix<-as.data.frame(fdr_matrix)
            
            fdr_adjust_pvalue_all<-apply(fdr_matrix,1,function(x){return(min(x,na.rm=TRUE)[1])})
            
            pvalues<-cbind(pvalues1,pvalues2,pvalues3)
            pvalues<-apply(pvalues,1,function(x){min(x,na.rm=TRUE)[1]})
            
            #pvalues1<-t(pvalues1)
            
            #print("here")
            pvalues1<-as.data.frame(pvalues1)
            pvalues1<-t(pvalues1)
            #print(dim(pvalues1))
            
            #pvalues2<-t(pvalues2)
            pvalues2<-as.data.frame(pvalues2)
            pvalues2<-t(pvalues2)
            
            #pvalues3<-t(pvalues3)
            pvalues3<-as.data.frame(pvalues3)
            pvalues3<-t(pvalues3)
            
          
            
            final.pvalues<-pvalues
            
            
            
            sel.diffdrthresh<-fdr_adjust_pvalue_all<fdrthresh & final.pvalues<pvalue.thresh
            
            if(length(which(fdr_adjust_pvalue1<fdrthresh))>0){
              
              
              
              
              X1=data_m_fc_withfeats[which(fdr_adjust_pvalue1<fdrthresh),]
              Y1=cbind(classlabels_orig[,1],as.character(classlabels_response_mat[,1]))
              Y1<-as.data.frame(Y1)
              
             
              
              if(output.device.type!="pdf"){
                
                temp_filename_1<-"Figures/HCA_Factor1selectedfeats.png"
                
                png(temp_filename_1,width=plots.width,height=plots.height,res=plots.res,type=plots.type,units="in")
              }
              
              hca_f1<-get_hca(feature_table_file=NA,parentoutput_dir=output_dir,class_labels_file=NA,X=X1,Y=Y1,heatmap.col.opt=heatmap.col.opt,cor.method=cor.method,is.data.znorm=FALSE,
                              analysismode="classification",
                      sample.col.opt=sample.col.opt,plots.width=2000,plots.height=2000,plots.res=300,
                      alphacol=0.3, hca_type=hca_type,newdevice=FALSE,input.type="intensity",mainlab="Factor1",
                      alphabetical.order=alphabetical.order,study.design=analysistype,labRow.value = labRow.value, labCol.value = labCol.value,similarity.matrix=similarity.matrix,cexLegend=hca.cex.legend)
              
              if(output.device.type!="pdf"){
                
                try(dev.off(),silent=TRUE)
              }
            }else{
              
              print("No significant features for Factor 1.")
            }
            
            if(length(which(fdr_adjust_pvalue2<fdrthresh))>0){
              X2=data_m_fc_withfeats[which(fdr_adjust_pvalue2<fdrthresh),]
              
              
              Y2=cbind(classlabels_orig[,1],as.character(classlabels_response_mat[,2]))
              Y2<-as.data.frame(Y2)
              
             
              
              if(output.device.type!="pdf"){
                
                temp_filename_1<-"Figures/HCA_Factor2selectedfeats.png"
                
                png(temp_filename_1,width=plots.width,height=plots.height,res=plots.res,type=plots.type,units="in")
              }
              
              
             hca_f2<-get_hca(feature_table_file=NA,parentoutput_dir=output_dir,class_labels_file=NA,X=X2,Y=Y2,heatmap.col.opt=heatmap.col.opt,cor.method=cor.method,is.data.znorm=FALSE,analysismode="classification",
                      sample.col.opt=sample.col.opt,plots.width=2000,plots.height=2000,plots.res=300, alphacol=0.3,
                      hca_type=hca_type,newdevice=FALSE,input.type="intensity",mainlab="Factor2",
                      alphabetical.order=alphabetical.order,study.design=analysistype,labRow.value = labRow.value, labCol.value = labCol.value,similarity.matrix=similarity.matrix,cexLegend=hca.cex.legend)
              
              
              if(output.device.type!="pdf"){
                
                try(dev.off(),silent=TRUE)
              }
              
            }else{
              
              print("No significant features for Factor 2.")
            }
            class_interact<-paste(classlabels_response_mat[,1],":",classlabels_response_mat[,2],sep="")
            #classlabels_response_mat[,1]:classlabels_response_mat[,2]
            
            if(length(which(fdr_adjust_pvalue3<fdrthresh))>0){
              X3=data_m_fc_withfeats[which(fdr_adjust_pvalue3<fdrthresh),]
              Y3=cbind(classlabels_orig[,1],class_interact)
              Y3<-as.data.frame(Y3)
              
              
              if(output.device.type!="pdf"){
                
                temp_filename_1<-"Figures/HCA_Factor1xFactor2selectedfeats.png"
                
                png(temp_filename_1,width=plots.width,height=plots.height,res=plots.res,type=plots.type,units="in")
              }
              
              hca_f3<-get_hca(feature_table_file=NA,parentoutput_dir=output_dir,class_labels_file=NA,X=X3,Y=Y3,heatmap.col.opt=heatmap.col.opt,cor.method=cor.method,is.data.znorm=FALSE,analysismode="classification",
                      sample.col.opt=sample.col.opt,plots.width=2000,plots.height=2000,plots.res=300, alphacol=0.3,
                      hca_type=hca_type,newdevice=FALSE,input.type="intensity",mainlab="Factor1 x Factor2",
                      alphabetical.order=alphabetical.order,study.design=analysistype,labRow.value = labRow.value, labCol.value = labCol.value,similarity.matrix=similarity.matrix,cexLegend=hca.cex.legend)
              
              
              if(output.device.type!="pdf"){
                
                try(dev.off(),silent=TRUE)
              }
              
              
            }else{
              
              print("No significant features for the interaction.")
            }
            
            
            
            data_limma_fdrall_withfeats<-cbind(final.pvalues,fdr_adjust_pvalue,data_m_fc_withfeats)
            
            cnames_tab<-colnames(data_m_fc_withfeats)
            cnames_tab<-c("P.value.Min(Factor1,Factor2,Interaction)","adjusted.P.value.Min(Factor1,Factor2,Interaction)",cnames_tab)
            colnames(data_limma_fdrall_withfeats)<-as.character(cnames_tab)
            
            
            
            #filename2<-"test2.txt"
            #data_limma_fdrsig_withfeats<-data_limma_fdrall_withfeats[sel.diffdrthresh==TRUE,]
            #write.table(data_limma_fdrsig_withfeats, file=filename2,sep="\t",row.names=FALSE)
            
            fdr_adjust_pvalue<-fdr_adjust_pvalue_all
            
          }
          
          if(featselmethod=="lm1wayanovarepeat"| featselmethod=="lmregrepeat"){
            
            
            
            
            
           # save(data_m_fc,classlabels_response_mat,subject_inf,modeltype,file="1waydebug.Rda")
            
            #clusterExport(cl,"classlabels_response_mat",envir = .GlobalEnv)
            #clusterExport(cl,"subject_inf",envir = .GlobalEnv)
            
            #res1<-apply(data_m_fc,1,function(x){
            
            if(featselmethod=="lm1wayanovarepeat"){
              
              print("Performing one-way ANOVA with repeated measurements analysis using nlme::lme()")
              
              numcores<-num_nodes #round(detectCores()*0.5)
              
              cl <- parallel::makeCluster(getOption("cl.cores", num_nodes))
              
              clusterExport(cl,"diffexplmonewayanovarepeat",envir = .GlobalEnv)
              clusterEvalQ(cl,library(nlme))
              clusterEvalQ(cl,library(multcomp))
              clusterEvalQ(cl,library(lsmeans))
              clusterExport(cl,"lme",envir = .GlobalEnv)
              clusterExport(cl,"interaction",envir = .GlobalEnv)
              clusterExport(cl,"anova",envir = .GlobalEnv)
              res1<-parApply(cl,data_m_fc,1,function(x,classlabels_response_mat,subject_inf,modeltype){
                
                #res1<-apply(data_m_fc,1,function(x){
                
                xvec<-x
                
                colnames(classlabels_response_mat)<-paste("Factor",seq(1,dim(classlabels_response_mat)[2]),sep="")
                
                data_mat_anova<-cbind(xvec,classlabels_response_mat)
                
                cnames<-colnames(data_mat_anova)
                cnames[1]<-"Response"
                
                colnames(data_mat_anova)<-cnames
                
                anova_res<-diffexplmonewayanovarepeat(dataA=data_mat_anova,subject_inf=subject_inf,modeltype=modeltype)
                
                return(anova_res)
                
                
              },classlabels_response_mat,subject_inf,modeltype)
              
              main_pval_mat<-{}
              
              posthoc_pval_mat<-{}
              pvalues<-{}
             
              
              bad_lm1feats<-{}
              
              ###saveres1,file="res1.Rda")
              for(i in 1:length(res1)){
                
                if(is.na(res1[[i]]$mainpvalues)==FALSE){
                  main_pval_mat<-rbind(main_pval_mat,res1[[i]]$mainpvalues)
                  pvalues<-c(pvalues,res1[[i]]$mainpvalues[1])
                  posthoc_pval_mat<-rbind(posthoc_pval_mat,res1[[i]]$posthoc)
                }else{
                  
                  bad_lm1feats<-c(bad_lm1feats,i)
                  
                  
                }
              }
              
              if(length(bad_lm1feats)>0){
                
                data_m_fc_withfeats<-data_m_fc_withfeats[-c(bad_lm1feats),]
                
                data_m_fc<-data_m_fc[-c(bad_lm1feats),]
              }
              #twoanovarepeat_res<-cbind(data_m_fc_withfeats[,c(1:2)],main_pval_mat,posthoc_pval_mat)
              
              #write.table(twoanovarepeat_res,file="Tables/lm2wayanovarepeat_with_posthoc_pvalues.txt",sep="\t",row.names=FALSE)
              
              
              pvalues1<-main_pval_mat[,1]
             
              
              onewayanova_res<-cbind(data_m_fc_withfeats[,c(1:2)],main_pval_mat,posthoc_pval_mat)
              
              # write.table(twoanova_res,file="twoanova_with_posthoc_pvalues.txt",sep="\t",row.names=FALSE)
              
              if(fdrmethod=="none"){
                fdr_adjust_pvalue1<-p.adjust(pvalues1,method="none")
               
              }
              
              if(fdrmethod=="BH"){
                fdr_adjust_pvalue1<-p.adjust(pvalues1,method="BH")
               
              }else{
                if(fdrmethod=="ST"){
                  #print(head(pvalues1))
                  #print(head(pvalues2))
                  #print(head(pvalues3))
                  #print(summary(pvalues1))
                  #print(summary(pvalues2))
                  #print(summary(pvalues3))
                  fdr_adjust_pvalue1<-try(qvalue(pvalues1),silent=TRUE)
                  
                  if(is(fdr_adjust_pvalue1,"try-error")){
                    
                    fdr_adjust_pvalue1<-qvalue(pvalues1,lambda=max(pvalues1,na.rm=TRUE))
                  }
                  
                  
                  
                  
                  fdr_adjust_pvalue1<-fdr_adjust_pvalue1$qvalues
                  
                 
                }else{
                  if(fdrmethod=="Strimmer"){
                    pdf("fdrtool.pdf")
                    #par_rows=1
                    #par(mfrow=c(par_rows,1))
                    fdr_adjust_pvalue1<-fdrtool(as.vector(pvalues1),statistic="pvalue",verbose=FALSE)
                    fdr_adjust_pvalue1<-fdr_adjust_pvalue1$qval
                    
                    
                    try(dev.off(),silent=TRUE)
                  }else{
                    if(fdrmethod=="none"){
                      fdr_adjust_pvalue1<-p.adjust(pvalues1,method="none")
                     
                      
                      
                    }else{
                      if(fdrmethod=="BY"){
                        fdr_adjust_pvalue1<-p.adjust(pvalues1,method="BY")
                     
                        
                      }else{
                        if(fdrmethod=="bonferroni"){
                          # fdr_adjust_pvalue<-p.adjust(pvalues,method="bonferroni")
                          fdr_adjust_pvalue1<-p.adjust(pvalues1,method="bonferroni")
                          
                        }
                      }
                    }
                  }
                }
                
                
                
              }
              
              if(fdrmethod=="none"){
                
                filename<-paste("Tables/",featselmethod,"_pvalall_withfeats.txt",sep="")
              }else{
                filename<-paste("Tables/",featselmethod,"_fdrall_withfeats.txt",sep="")
              }
              cnames_tab<-colnames(data_m_fc_withfeats)
              
              posthoc_names<-colnames(posthoc_pval_mat)
              
              #
              cnames_tab<-c("Factor1.P.value","Factor1.adjusted.P.value",posthoc_names,cnames_tab)
              
              
              data_limma_fdrall_withfeats<-cbind(pvalues1,fdr_adjust_pvalue1,posthoc_pval_mat,data_m_fc_withfeats)
              
              
              
              colnames(data_limma_fdrall_withfeats)<-as.character(cnames_tab)
              
              #gohere
              
              if(length(check_names)>0){
                
                
                data_limma_fdrall_withfeats<-cbind(pvalues1,fdr_adjust_pvalue1,posthoc_pval_mat,data_m_fc_with_names,data_m_fc_withfeats[,-c(1:2)])
                 
                colnames(data_limma_fdrall_withfeats)<-as.character(cnames_tab)
                
                data_limma_fdrall_withfeats<-as.data.frame(data_limma_fdrall_withfeats)
                #data_limma_fdrall_withfeats<-cbind(p.value,adjusted.p.value,results2,data_m_fc_with_names,data_m_fc_withfeats[,-c(1:2)])
                
                rem_col_ind1<-grep(colnames(data_limma_fdrall_withfeats),pattern=c("mz"))
                
                rem_col_ind2<-grep(colnames(data_limma_fdrall_withfeats),pattern=c("time"))
                
                rem_col_ind<-c(rem_col_ind1,rem_col_ind2)
                
              }else{
                rem_col_ind<-{}
              }
              
              if(length(rem_col_ind)>0){
                
                write.table(data_limma_fdrall_withfeats[,-c(rem_col_ind)],file="Tables/onewayanovarepeat_with_posthoc_comparisons.txt",sep="\t",row.names=FALSE)
                
              }else{
                
                write.table(data_limma_fdrall_withfeats,file="Tables/onewayanovarepeat_with_posthoc_comparisons.txt",sep="\t",row.names=FALSE)
                
                
              }
              
              #data_limma_fdrall_withfeats<-data_limma_fdrall_withfeats[order(fdr_adjust_fpvalue),]
              
              filename<-paste("Tables/",filename,sep="")
              
              
              
              
             
            
              fdr_adjust_pvalue<-fdr_adjust_pvalue1
              
              final.pvalues<-pvalues1
              
              sel.diffdrthresh<-fdr_adjust_pvalue1<fdrthresh & final.pvalues<pvalue.thresh
              
              
              
            }else{
              
              print("Performing linear regression with repeated measurements analysis using nlme::lme()")
              numcores<-num_nodes #round(detectCores()*0.5)
              
              cl <- parallel::makeCluster(getOption("cl.cores", num_nodes))
              
              clusterExport(cl,"diffexplmregrepeat",envir = .GlobalEnv)
              clusterEvalQ(cl,library(nlme))
              clusterEvalQ(cl,library(multcomp))
              clusterEvalQ(cl,library(lsmeans))
              clusterExport(cl,"lme",envir = .GlobalEnv)
              clusterExport(cl,"interaction",envir = .GlobalEnv)
              clusterExport(cl,"anova",envir = .GlobalEnv)
              res1<-parApply(cl,data_m_fc,1,function(x,classlabels_response_mat,subject_inf,modeltype){
                
                #res1<-apply(data_m_fc,1,function(x){
                
                xvec<-x
                
                colnames(classlabels_response_mat)<-paste("Factor",seq(1,dim(classlabels_response_mat)[2]),sep="")
                data_mat_anova<-cbind(xvec,classlabels_response_mat)
                
                cnames<-colnames(data_mat_anova)
                cnames[1]<-"Response"
                
                colnames(data_mat_anova)<-cnames
                
            #    save(data_mat_anova,subject_inf,modeltype,file="lmregdebug.Rda")
                
                if(ncol(data_mat_anova)>2){
                  
                  covar.matrix=classlabels_response_mat[,-c(1)]
                  
                }else{
                  covar.matrix=NA 
                }
                
                anova_res<-diffexplmregrepeat(dataA=data_mat_anova,subject_inf=subject_inf,modeltype=modeltype,covar.matrix = covar.matrix)
                
                return(anova_res)
                
                
              },classlabels_response_mat,subject_inf,modeltype)
              
            
            stopCluster(cl)
            
            main_pval_mat<-{}
            pvalues<-{}
            
           # save(res1,file="lmres1.Rda")
            
            posthoc_pval_mat<-{}
            
            bad_lm1feats<-{}
            
            res2<-t(res1)
            
            res2<-as.data.frame(res2)
            colnames(res2)<-c("pvalue","coefficient","std.error","t.value")
            pvalues<-res2$pvalue
            
            pvalues<-unlist(pvalues)
            if(fdrmethod=="BH"){
              fdr_adjust_pvalue<-p.adjust(pvalues,method="BH")
            }else{
              if(fdrmethod=="ST"){
                
                fdr_adjust_pvalue<-try(qvalue(pvalues),silent=TRUE)
                
                if(is(fdr_adjust_pvalue,"try-error")){
                  
                  fdr_adjust_pvalue<-qvalue(pvalues,lambda=max(pvalues,na.rm=TRUE))
                }
                
                fdr_adjust_pvalue<-fdr_adjust_pvalue$qvalues
                
                
              }else{
                if(fdrmethod=="Strimmer"){
                  pdf("fdrtool.pdf")
                  #par_rows=1
                  #par(mfrow=c(par_rows,1))
                  fdr_adjust_pvalue<-suppressWarnings(fdrtool(as.vector(pvalues),statistic="pvalue",verbose=FALSE))
                  fdr_adjust_pvalue<-fdr_adjust_pvalue$qval
                  try(dev.off(),silent=TRUE)
                }else{
                  if(fdrmethod=="none"){
                    fdr_adjust_pvalue<-pvalues
                    
                  }else{
                    if(fdrmethod=="BY"){
                      fdr_adjust_pvalue<-p.adjust(pvalues,method="BY")
                    }else{
                      if(fdrmethod=="bonferroni"){
                        fdr_adjust_pvalue<-p.adjust(pvalues,method="bonferroni")
                        
                      }
                    }
                  }
                }
              }
              
              
              
            }
            if(fdrmethod=="none"){
              filename<-paste(featselmethod,"_pvalall_withfeats.txt",sep="")
            }else{
              filename<-paste(featselmethod,"_fdrall_withfeats.txt",sep="")
            }
            cnames_tab<-colnames(data_m_fc_withfeats)
            
            
            
            cnames_tab<-c("P.value","adjusted.P.value",c("coefficient","std.error","t.value"),cnames_tab)
            
            
            pvalues<-as.data.frame(pvalues)
            
            final.pvalues<-pvalues
            sel.diffdrthresh<-fdr_adjust_pvalue<fdrthresh & final.pvalues<pvalue.thresh
            #pvalues<-t(pvalues)
            #print(dim(pvalues))
            #print(dim(data_m_fc_withfeats))
            
            if(length(bad_lm1feats)>0){
              data_m_fc_withfeats<-data_m_fc_withfeats[-c(bad_lm1feats),]
              
              data_m_fc<-data_m_fc[-c(bad_lm1feats),]
            }
            
            data_limma_fdrall_withfeats<-cbind(pvalues,fdr_adjust_pvalue,res2[,-c(1)],data_m_fc_withfeats)
            
            colnames(data_limma_fdrall_withfeats)<-as.character(cnames_tab)
            
            
            
            #data_limma_fdrall_withfeats<-data_limma_fdrall_withfeats[order(fdr_adjust_fpvalue),]
            
            filename<-paste("Tables/",filename,sep="")
          #  write.table(data_limma_fdrall_withfeats, file=filename,sep="\t",row.names=FALSE)
            
            data_limma_fdrall_withfeats<-cbind(pvalues,fdr_adjust_pvalue,data_m_fc_withfeats)
          }
            
          }
          
          if(featselmethod=="lm2wayanovarepeat"){
            
            print("Performing two-way ANOVA with repeated measurements analysis using nlme::lme()")
            
            
            numcores<-num_nodes #round(detectCores()*0.5)
            
            cl <- parallel::makeCluster(getOption("cl.cores", num_nodes))
            
            clusterExport(cl,"diffexplmtwowayanovarepeat",envir = .GlobalEnv)
            clusterEvalQ(cl,library(nlme))
            clusterEvalQ(cl,library(multcomp))
            clusterEvalQ(cl,library(lsmeans))
            clusterExport(cl,"lme",envir = .GlobalEnv)
            clusterExport(cl,"interaction",envir = .GlobalEnv)
            clusterExport(cl,"anova",envir = .GlobalEnv)
            
            #clusterExport(cl,"classlabels_response_mat",envir = .GlobalEnv)
            #clusterExport(cl,"subject_inf",envir = .GlobalEnv)
            
            #res1<-apply(data_m_fc,1,function(x){
            
            #	print(dim(data_m_fc))
            #	print(dim(classlabels_response_mat))
            
            res1<-parApply(cl,data_m_fc,1,function(x,classlabels_response_mat,subject_inf,modeltype){
              
              #  res1<-apply(data_m_fc,1,function(x){
              
              # ###saveclasslabels_response_mat,file="classlabels_response_mat.Rda")
              
              #       ###savesubject_inf,file="subject_inf.Rda")
              
              xvec<-x
              
              ####savexvec,file="xvec.Rda")
              colnames(classlabels_response_mat)<-paste("Factor",seq(1,dim(classlabels_response_mat)[2]),sep="")
              
              data_mat_anova<-cbind(xvec,classlabels_response_mat)
              
              cnames<-colnames(data_mat_anova)
              cnames[1]<-"Response"
              
              colnames(data_mat_anova)<-cnames
              
              #print(subject_inf)
              #print(dim(data_mat_anova))
              
              
              subject_inf<-as.data.frame(subject_inf)
              #print(dim(subject_inf))
              
              anova_res<-diffexplmtwowayanovarepeat(dataA=data_mat_anova,subject_inf=subject_inf[,1],modeltype=modeltype)
              
              return(anova_res)
            },classlabels_response_mat,subject_inf,modeltype)
            
            
            main_pval_mat<-{}
            
            stopCluster(cl)
            posthoc_pval_mat<-{}
            
            #print(head(res1))
            #	print("here")
            
            pvalues<-{}
            
            
            bad_lm1feats<-{}
            
            ###saveres1,file="res1.Rda")
            for(i in 1:length(res1)){
              
              if(is.na(res1[[i]]$mainpvalues)==FALSE){
                main_pval_mat<-rbind(main_pval_mat,res1[[i]]$mainpvalues)
                pvalues<-c(pvalues,res1[[i]]$mainpvalues[1])
                posthoc_pval_mat<-rbind(posthoc_pval_mat,res1[[i]]$posthoc)
              }else{
                
                bad_lm1feats<-c(bad_lm1feats,i)
                
                
              }
            }
            
            if(length(bad_lm1feats)>0){
              
              data_m_fc_withfeats<-data_m_fc_withfeats[-c(bad_lm1feats),]
              
              data_m_fc<-data_m_fc[-c(bad_lm1feats),]
            }
            twoanovarepeat_res<-cbind(data_m_fc_withfeats[,c(1:2)],main_pval_mat,posthoc_pval_mat)
            
            #write.table(twoanovarepeat_res,file="Tables/lm2wayanovarepeat_with_posthoc_pvalues.txt",sep="\t",row.names=FALSE)
            
            
            pvalues1<-main_pval_mat[,1]
            pvalues2<-main_pval_mat[,2]
            pvalues3<-main_pval_mat[,3]
            
            twoanova_res<-cbind(data_m_fc_withfeats[,c(1:2)],main_pval_mat,posthoc_pval_mat)
            
            # write.table(twoanova_res,file="twoanova_with_posthoc_pvalues.txt",sep="\t",row.names=FALSE)
            
            if(fdrmethod=="none"){
              fdr_adjust_pvalue1<-p.adjust(pvalues1,method="none")
              fdr_adjust_pvalue2<-p.adjust(pvalues2,method="none")
              fdr_adjust_pvalue3<-p.adjust(pvalues3,method="none")
            }
            
            if(fdrmethod=="BH"){
              fdr_adjust_pvalue1<-p.adjust(pvalues1,method="BH")
              fdr_adjust_pvalue2<-p.adjust(pvalues2,method="BH")
              fdr_adjust_pvalue3<-p.adjust(pvalues3,method="BH")
            }else{
              if(fdrmethod=="ST"){
                #print(head(pvalues1))
                #print(head(pvalues2))
                #print(head(pvalues3))
                #print(summary(pvalues1))
                #print(summary(pvalues2))
                #print(summary(pvalues3))
                fdr_adjust_pvalue1<-try(qvalue(pvalues1),silent=TRUE)
                fdr_adjust_pvalue2<-try(qvalue(pvalues2),silent=TRUE)
                fdr_adjust_pvalue3<-try(qvalue(pvalues3),silent=TRUE)
                
                if(is(fdr_adjust_pvalue1,"try-error")){
                  
                  fdr_adjust_pvalue1<-qvalue(pvalues1,lambda=max(pvalues1,na.rm=TRUE))
                }
                
                if(is(fdr_adjust_pvalue2,"try-error")){
                  fdr_adjust_pvalue2<-qvalue(pvalues2,lambda=max(pvalues2,na.rm=TRUE))
                }
                
                if(is(fdr_adjust_pvalue3,"try-error")){
                  fdr_adjust_pvalue3<-qvalue(pvalues3,lambda=max(pvalues3,na.rm=TRUE))
                }
                
                
                fdr_adjust_pvalue1<-fdr_adjust_pvalue1$qvalues
                
                fdr_adjust_pvalue2<-fdr_adjust_pvalue2$qvalues
                
                fdr_adjust_pvalue3<-fdr_adjust_pvalue3$qvalues
              }else{
                if(fdrmethod=="Strimmer"){
                  pdf("fdrtool.pdf")
                  #par_rows=1
                  #par(mfrow=c(par_rows,1))
                  fdr_adjust_pvalue1<-fdrtool(as.vector(pvalues1),statistic="pvalue",verbose=FALSE)
                  fdr_adjust_pvalue1<-fdr_adjust_pvalue1$qval
                  
                  fdr_adjust_pvalue2<-fdrtool(as.vector(pvalues2),statistic="pvalue",verbose=FALSE)
                  fdr_adjust_pvalue2<-fdr_adjust_pvalue2$qval
                  
                  fdr_adjust_pvalue3<-fdrtool(as.vector(pvalues3),statistic="pvalue",verbose=FALSE)
                  fdr_adjust_pvalue3<-fdr_adjust_pvalue3$qval
                  try(dev.off(),silent=TRUE)
                }else{
                  if(fdrmethod=="none"){
                    fdr_adjust_pvalue1<-p.adjust(pvalues1,method="none")
                    fdr_adjust_pvalue2<-p.adjust(pvalues2,method="none")
                    fdr_adjust_pvalue3<-p.adjust(pvalues3,method="none")
                    
                    
                  }else{
                    if(fdrmethod=="BY"){
                      fdr_adjust_pvalue1<-p.adjust(pvalues1,method="BY")
                      fdr_adjust_pvalue2<-p.adjust(pvalues2,method="BY")
                      fdr_adjust_pvalue3<-p.adjust(pvalues3,method="BY")
                      
                      
                    }else{
                      if(fdrmethod=="bonferroni"){
                        # fdr_adjust_pvalue<-p.adjust(pvalues,method="bonferroni")
                        fdr_adjust_pvalue1<-p.adjust(pvalues1,method="bonferroni")
                        fdr_adjust_pvalue2<-p.adjust(pvalues2,method="bonferroni")
                        fdr_adjust_pvalue3<-p.adjust(pvalues3,method="bonferroni")
                      }
                    }
                  }
                }
              }
              
              
              
            }
            
            if(fdrmethod=="none"){
              
              filename<-paste("Tables/",featselmethod,"_pvalall_withfeats.txt",sep="")
            }else{
              filename<-paste("Tables/",featselmethod,"_fdrall_withfeats.txt",sep="")
            }
            cnames_tab<-colnames(data_m_fc_withfeats)
            
            posthoc_names<-colnames(posthoc_pval_mat)
            
            #
            cnames_tab<-c("Factor1.P.value","Factor1.adjusted.P.value","Factor2.P.value","Factor2.adjusted.P.value","Interact.P.value","Interact.adjusted.P.value",posthoc_names,cnames_tab)
            
            
            data_limma_fdrall_withfeats<-cbind(pvalues1,fdr_adjust_pvalue1,pvalues2,fdr_adjust_pvalue2,pvalues3,fdr_adjust_pvalue3,posthoc_pval_mat,data_m_fc_withfeats)
            
            
            
            colnames(data_limma_fdrall_withfeats)<-as.character(cnames_tab)
            
            
            if(length(check_names)>0){
              
              data_limma_fdrall_withfeats<-cbind(pvalues1,fdr_adjust_pvalue1,pvalues2,fdr_adjust_pvalue2,pvalues3,fdr_adjust_pvalue3,posthoc_pval_mat,data_m_fc_with_names,data_m_fc_withfeats[,-c(1:2)])
              
              colnames(data_limma_fdrall_withfeats)<-as.character(cnames_tab)
            
              data_limma_fdrall_withfeats<-as.data.frame(data_limma_fdrall_withfeats)
              #data_limma_fdrall_withfeats<-cbind(p.value,adjusted.p.value,results2,data_m_fc_with_names,data_m_fc_withfeats[,-c(1:2)])
              
              rem_col_ind1<-grep(colnames(data_limma_fdrall_withfeats),pattern=c("mz"))
              
              rem_col_ind2<-grep(colnames(data_limma_fdrall_withfeats),pattern=c("time"))
              
              rem_col_ind<-c(rem_col_ind1,rem_col_ind2)
              
            }else{
              rem_col_ind<-{}
            }
            
            if(length(rem_col_ind)>0){
              write.table(data_limma_fdrall_withfeats[,-c(rem_col_ind)], file="Tables/twowayanovarepeat_with_posthoc_comparisons.txt",sep="\t",row.names=FALSE)
            }else{
              
              
              #write.table(data_limma_fdrall_withfeats,file="Tables/twowayanova_with_posthoc_comparisons.txt",sep="\t",row.names=FALSE)
              
              write.table(data_limma_fdrall_withfeats,file="Tables/twowayanovarepeat_with_posthoc_comparisons.txt",sep="\t",row.names=FALSE)
              
            }
            
            #filename<-paste("Tables/",filename,sep="")
            
            
            #data_limma_fdrall_withfeats<-data_limma_fdrall_withfeats[order(fdr_adjust_fpvalue),]
            
            filename<-paste("Tables/",filename,sep="")
            
            
            
            fdr_matrix<-cbind(fdr_adjust_pvalue1,fdr_adjust_pvalue2,fdr_adjust_pvalue3)
            
            fdr_matrix<-as.data.frame(fdr_matrix)
            
            fdr_adjust_pvalue_all<-apply(fdr_matrix,1,function(x){return(min(x,na.rm=TRUE))})
            
            pvalues_all<-cbind(pvalues1,pvalues2,pvalues3)
            pvalue_matrix<-as.data.frame(pvalues_all)
            
            pvalue_all<-apply(pvalue_matrix,1,function(x){return(min(x,na.rm=TRUE)[1])})
            
            
            #pvalues1<-t(pvalues1)
            
            #print("here")
            #pvalues1<-as.data.frame(pvalues1)
            #pvalues1<-t(pvalues1)
            #print(dim(pvalues1))
            
            #pvalues2<-t(pvalues2)
            #pvalues2<-as.data.frame(pvalues2)
            #pvalues2<-t(pvalues2)
            
            #pvalues3<-t(pvalues3)
            #pvalues3<-as.data.frame(pvalues3)
            #pvalues3<-t(pvalues3)
            
            #pvalues<-t(pvalues)
            #print(dim(pvalues1))
            #print(dim(pvalues2))
            #print(dim(pvalues3))
            #print(dim(data_m_fc_withfeats))
            
            pvalues<-pvalue_all
            
            final.pvalues<-pvalues
            sel.diffdrthresh<-fdr_adjust_pvalue_all<fdrthresh & final.pvalues<pvalue.thresh
            
            
            if(length(which(fdr_adjust_pvalue1<fdrthresh))>0){
              X1=data_m_fc_withfeats[which(fdr_adjust_pvalue1<fdrthresh),]
              Y1=cbind(classlabels_orig[,1],as.character(classlabels_response_mat[,1]))
              Y1<-as.data.frame(Y1)
              
              ###saveclasslabels_orig,file="classlabels_orig.Rda")
              ###saveclasslabels_response_mat,file="classlabels_response_mat.Rda")
              
              print("Performing HCA using features selected for Factor1")
              if(output.device.type!="pdf"){
                
                temp_filename_1<-"Figures/HCA_Factor1selectedfeats.png"
                
                png(temp_filename_1,width=plots.width,height=plots.height,res=plots.res,type=plots.type,units="in")
              }
              
              
              hca_f1<-get_hca(feature_table_file=NA,parentoutput_dir=output_dir,class_labels_file=NA,X=X1,Y=Y1,heatmap.col.opt=heatmap.col.opt,cor.method=cor.method,is.data.znorm=FALSE,analysismode="classification",
                      sample.col.opt=sample.col.opt,plots.width=2000,plots.height=2000,plots.res=300, 
                      alphacol=0.3, hca_type=hca_type,newdevice=FALSE,input.type="intensity",mainlab="Factor 1",
                      alphabetical.order=alphabetical.order,study.design="oneway",labRow.value = labRow.value, labCol.value = labCol.value,similarity.matrix=similarity.matrix,cexLegend=hca.cex.legend)
              
              if(output.device.type!="pdf"){
                
                try(dev.off(),silent=TRUE)
              }
            }else{
              print("No significant features for Factor 1.")
              
            }
            
            if(length(which(fdr_adjust_pvalue2<fdrthresh))>0){
              
              X2=data_m_fc_withfeats[which(fdr_adjust_pvalue2<fdrthresh),]
              Y2=cbind(classlabels_orig[,1],as.character(classlabels_response_mat[,2]))
              Y2<-as.data.frame(Y2)
              
              
              if(output.device.type!="pdf"){
                
                temp_filename_1<-"Figures/HCA_Factor2selectedfeats.png"
                
                png(temp_filename_1,width=plots.width,height=plots.height,res=plots.res,type=plots.type,units="in")
              }
              
              print("Performing HCA using features selected for Factor2")
              hca_f2<-get_hca(feature_table_file=NA,parentoutput_dir=output_dir,class_labels_file=NA,X=X2,Y=Y2,heatmap.col.opt=heatmap.col.opt,cor.method=cor.method,is.data.znorm=FALSE,analysismode="classification",
                      sample.col.opt=sample.col.opt,plots.width=2000,plots.height=2000,plots.res=300, 
                      alphacol=0.3, hca_type=hca_type,newdevice=FALSE,input.type="intensity",mainlab="Factor 2",
                      alphabetical.order=alphabetical.order,study.design="oneway",labRow.value = labRow.value, labCol.value = labCol.value,similarity.matrix=similarity.matrix,cexLegend=hca.cex.legend)
              
              if(output.device.type!="pdf"){
                
                try(dev.off(),silent=TRUE)
              }
              
            }else{
              
              print("No significant features for Factor 2.")
            }
            
            class_interact<-paste(classlabels_response_mat[,1],":",classlabels_response_mat[,2],sep="") #classlabels_response_mat[,1]:classlabels_response_mat[,2]
            
            if(length(which(fdr_adjust_pvalue3<fdrthresh))>0){
              
              X3=data_m_fc_withfeats[which(fdr_adjust_pvalue3<fdrthresh),]
              Y3=cbind(classlabels_orig[,1],class_interact)
              Y3<-as.data.frame(Y3)
              
              
              if(output.device.type!="pdf"){
                
                temp_filename_1<-"Figures/HCA_Factor1xFactor2selectedfeats.png"
                
                png(temp_filename_1,width=plots.width,height=plots.height,res=plots.res,type=plots.type,units="in")
              }
              
              print("Performing HCA using features selected for Factor1x2")
              hca_f3<-get_hca(feature_table_file=NA,parentoutput_dir=output_dir,class_labels_file=NA,X=X3,Y=Y3,heatmap.col.opt=heatmap.col.opt,cor.method=cor.method,is.data.znorm=FALSE,analysismode="classification",
                      sample.col.opt=sample.col.opt,plots.width=2000,plots.height=2000,plots.res=300,
                      alphacol=0.3, hca_type=hca_type,newdevice=FALSE,input.type="intensity",mainlab="Factor 1 x Factor 2",
                      alphabetical.order=alphabetical.order,study.design="oneway",labRow.value = labRow.value, labCol.value = labCol.value,similarity.matrix=similarity.matrix,cexLegend=hca.cex.legend)
              
              
              if(output.device.type!="pdf"){
                
                try(dev.off(),silent=TRUE)
              }
            }else{
              print("No significant features for Factor 1x2 interaction.")
            }
            
            
            #data_limma_fdrall_withfeats<-cbind(pvalues,fdr_adjust_pvalue,posthoc_pval_mat,data_m_fc_withfeats)
            
            
            #
            data_limma_fdrall_withfeats<-cbind(pvalues1,fdr_adjust_pvalue1,pvalues2,fdr_adjust_pvalue2,pvalues3,fdr_adjust_pvalue3,posthoc_pval_mat,data_m_fc_withfeats)
            
            fdr_adjust_pvalue<-cbind(fdr_adjust_pvalue1,fdr_adjust_pvalue2,fdr_adjust_pvalue3)
            fdr_adjust_pvalue<-apply(fdr_adjust_pvalue,1,function(x){min(x,na.rm=TRUE)})
            
            colnames(data_limma_fdrall_withfeats)<-as.character(cnames_tab)
            
            #data_limma_fdrall_withfeats<-data_limma_fdrall_withfeats[order(fdr_adjust_fpvalue),]
            #write.table(data_limma_fdrall_withfeats, file=filename,sep="\t",row.names=FALSE)
            
            data_limma_fdrall_withfeats<-cbind(final.pvalues,fdr_adjust_pvalue,data_m_fc_withfeats)
            
            cnames_tab<-colnames(data_m_fc_withfeats)
            cnames_tab<-c("P.value.Min(Factor1,Factor2,Interaction)","adjusted.P.value.Min(Factor1,Factor2,Interaction)",cnames_tab)
            colnames(data_limma_fdrall_withfeats)<-as.character(cnames_tab)
            
            #filename2<-"test2.txt"
            #data_limma_fdrsig_withfeats<-data_limma_fdrall_withfeats[sel.diffdrthresh==TRUE,]
            #write.table(data_limma_fdrsig_withfeats, file=filename2,sep="\t",row.names=FALSE)
            
            fdr_adjust_pvalue<-fdr_adjust_pvalue_all
            
          }
          
          
          
          
          
          
          
        }
        
        
        
        
        
        #end of feature selection methods
        
        
        
        
        
        
        
        
        if(featselmethod=="lmreg" | featselmethod=="lm1wayanova" | featselmethod=="lm2wayanova" | featselmethod=="lm1wayanovarepeat" | featselmethod=="lm2wayanovarepeat"
           | featselmethod=="limma" | featselmethod=="limma2way" | featselmethod=="logitreg" | featselmethod=="limma2wayrepeat" | featselmethod=="wilcox" | featselmethod=="ttest" |  featselmethod=="poissonreg" | featselmethod=="lmregrepeat")
        {
          
          
          
          sel.diffdrthresh<-fdr_adjust_pvalue<fdrthresh & final.pvalues<pvalue.thresh
          
          
          
          goodip<-which(sel.diffdrthresh==TRUE)
          
          
          classlabels<-as.data.frame(classlabels)
          
          
          
          
        #  if(featselmethod=="limma2way"){
         #   vennDiagram(results2,cex=0.8)
        #  }
          
          #print(summary(fdr_adjust_pvalue))
          #pheadrint(summary(final.pvalues))
          
        }
        
        pred_acc<-0 #("NA")
        
        
        
        #print("here")
        feat_sigfdrthresh[lf]<-length(goodip) #which(sel.diffdrthresh==TRUE))
        if(kfold>dim(data_m_fc)[2]){
          kfold=dim(data_m_fc)[2]
        }
        if(analysismode=="classification"){
          
          
          
          #print("classification")
          
          if(length(goodip)>0 & dim(data_m_fc)[2]>=kfold){
            
            #save(classlabels,classlabels_orig, data_m_fc,file="debug2.rda")
                 
                 
            if(alphabetical.order==FALSE){
              Targetvar <- factor(classlabels[,1], levels=unique(classlabels[,1]))
            }else{
              
              Targetvar<-factor(classlabels[,1])
            }
            
            
            dataA<-cbind(Targetvar,t(data_m_fc))
            
            dataA<-as.data.frame(dataA)
            
           
            
            dataA$Targetvar<-factor(Targetvar)
            
            
            #df.summary <- dataA %>% group_by(Targetvar) %>%  summarize_all(funs(mean))
            
            # df.summary <- dataA %>% group_by(Targetvar) %>%  summarize_all(funs(mean))
            
            dataA[,-c(1)]<-apply(dataA[,-c(1)],2,function(x){as.numeric(as.character(x))})
            
            if(alphabetical.order==FALSE){
              dataA$Targetvar <- factor(dataA$Targetvar, levels=unique(dataA$Targetvar))
            }
            df.summary <-aggregate(x=dataA,by=list(as.factor(dataA$Targetvar)),function(x){mean(x,na.rm=TRUE)})
            
            #save(dataA,file="errordataA.Rda")
            df.summary.sd <-aggregate(x=dataA[,-c(1)],by=list(as.factor(dataA$Targetvar)),function(x){sd(x,na.rm=TRUE)})
            
            df2<-as.data.frame(df.summary[,-c(1:2)])
            
            group_means<-t(df.summary)
            
           # save(classlabels,classlabels_orig, classlabels_class,Targetvar,dataA,data_m_fc,df.summary,df2,group_means,file="debugfoldchange.Rda")
            
            colnames(group_means)<-paste("mean",levels(as.factor(dataA$Targetvar)),sep="") #paste("Group",seq(1,length(unique(dataA$Targetvar))),sep="")
            
            group_means<-cbind(data_m_fc_withfeats[,c(1:2)],group_means[-c(1:2),])
            
            group_sd<-t(df.summary.sd)
            
            colnames(group_sd)<-paste("std.dev",levels(as.factor(dataA$Targetvar)),sep="") #paste("Group",seq(1,length(unique(dataA$Targetvar))),sep="")
            
            group_sd<-cbind(data_m_fc_withfeats[,c(1:2)],group_sd[-c(1),])
            
            # write.table(group_means,file="group_means.txt",sep="\t",row.names=FALSE)
            
            
            #  ###savedf2,file="df2.Rda")
            #    ###savedataA,file="dataA.Rda")
            #   ###saveTargetvar,file="Targetvar.Rda")
            
            if(log2transform==TRUE || input.intensity.scale=="log2"){
              
              
              cl <- parallel::makeCluster(getOption("cl.cores", num_nodes))
              foldchangeres<-parApply(cl,df2,2,function(x){
                
                
                res<-lapply(1:length(x),function(i){
                  return((x[i]-x[-i]))
                  
                })
                res<-unlist(res)
                
                tempres<-abs(res)
                res_ind<-which(tempres==max(tempres,na.rm=TRUE))
                return(res[res_ind[1]])
                
              })
              
              stopCluster(cl)
              
              
              print("Using log2 fold change threshold of")
              print(foldchangethresh)
              
              
            }else{
              
              #raw intensities
              if(znormtransform==FALSE)
              {
                #   foldchangeres<-apply(log2(df2+1),2,function(x){res<-{};for(i in 1:length(x)){res<-c(res,(x[i]-x[-i]));};tempres<-abs(res);res_ind<-which(tempres==max(tempres,na.rm=TRUE));return(res[res_ind[1]]);})
                
                if(FALSE){
                  foldchangeres<-apply(log2(df2+log2.transform.constant),2,dist)
                  
                  if(length(nrow(foldchangeres))>0){
                    
                    foldchangeres<-apply(foldchangeres,2,function(x)
                    {
                      
                      max_ind<-which(x==max(abs(x)))[1];
                      return(x[max_ind])
                      
                    }
                    )
                    
                  }
                  
                }
                
                cl <- parallel::makeCluster(getOption("cl.cores", num_nodes))
                foldchangeres<-parApply(cl,log2(df2+0.0000001),2,function(x){
                  
                  
                  res<-lapply(1:length(x),function(i){
                    return((x[i]-x[-i]))
                    
                  })
                  res<-unlist(res)
                  
                  tempres<-abs(res)
                  res_ind<-which(tempres==max(tempres,na.rm=TRUE))
                  return(res[res_ind[1]])
                  
                })
                
                stopCluster(cl)
                
                
                foldchangethresh=foldchangethresh
                print("Using raw fold change threshold of")
                print(foldchangethresh)
                
              }else{
                
                #  foldchangeres<-apply(df2,2,function(x){res<-{};for(i in 1:length(x)){res<-c(res,(x[i]-(x[-i])));};tempres<-abs(res);res_ind<-which(tempres==max(tempres,na.rm=TRUE));return(res[res_ind[1]]);})
                
                if(FALSE){
                  foldchangeres<-apply(df2,2,dist)
                  
                  if(length(nrow(foldchangeres))>0){
                    
                    foldchangeres<-apply(foldchangeres,2,function(x)
                    {
                      
                      max_ind<-which(x==max(abs(x)))[1];
                      return(x[max_ind])
                      
                    }
                    )
                    
                  }
                  
                }
                
                cl <- parallel::makeCluster(getOption("cl.cores", num_nodes))
                foldchangeres<-parApply(cl,df2,2,function(x){
                  
                  
                  res<-lapply(1:length(x),function(i){
                    return((x[i]-x[-i]))
                    
                  })
                  res<-unlist(res)
                  
                  tempres<-abs(res)
                  res_ind<-which(tempres==max(tempres,na.rm=TRUE))
                  return(res[res_ind[1]])
                  
                })
                
                stopCluster(cl)
                
                #print(summary(foldchangeres))
                
                
                #foldchangethresh=2^foldchangethresh
                print("Using Z-score change threshold of")
                print(foldchangethresh)
                
              }
              
              
            }
            
            
            
            
            if(length(class_labels_levels)==2){
              
              
              zvec=foldchangeres
            }else{
              
              zvec=NA
              
              if(featselmethod=="lmreg" && analysismode=="regression"){
                
                cnames_matrix<-colnames(data_limma_fdrall_withfeats)
                cnames_colindex<-grep("Estimate_",cnames_matrix)
                
                
                
                zvec<-data_limma_fdrall_withfeats[,c(cnames_colindex[1])]
              }
              
              
            }
            
            
            maxfoldchange<-foldchangeres
            
            goodipfoldchange<-which(abs(maxfoldchange)>foldchangethresh)
            
            #if(FALSE)
            {
              if(input.intensity.scale=="raw" && log2transform==FALSE && znormtransform==FALSE){
                
                foldchangeres<-2^((foldchangeres))
                
                
                
                
                
              }
            }
            
            maxfoldchange1<-foldchangeres
            
            roundUpNice <- function(x, nice=c(1,2,4,5,6,8,10))
            {
              if(length(x) != 1) stop("'x' must be of length 1")
              10^floor(log10(x)) * nice[[which(x <= 10^floor(log10(x)) * nice)[[1]]]]
            }
            
            d4<-as.data.frame(data_limma_fdrall_withfeats)
            
            max_mz_val<-roundUpNice(max(d4$mz)[1])
            max_time_val<-roundUpNice(max(d4$time)[1])
            
            x1increment=round_any(max_mz_val/10,10,f=floor)
            
            
            x2increment=round_any(max_time_val/10,10,f=floor)
            
            if(x2increment<1){
              x2increment=0.5
              
            }
            
            if(x1increment<1){
              x1increment=0.5
              
            }
            
            
            if(featselmethod=="lmreg" | featselmethod=="lm1wayanova" | featselmethod=="lm2wayanova" | featselmethod=="lm1wayanovarepeat" | featselmethod=="lm2wayanovarepeat"
               | featselmethod=="limma" | featselmethod=="limma2way" | featselmethod=="logitreg" | featselmethod=="limma2wayrepeat" | featselmethod=="wilcox" | featselmethod=="ttest" |  featselmethod=="poissonreg" | featselmethod=="lmregrepeat")
            {
              
              
              
              # print("Plotting manhattan plots")
              
              sel.diffdrthresh<-fdr_adjust_pvalue<fdrthresh & final.pvalues<pvalue.thresh
              
              
              goodip<-which(sel.diffdrthresh==TRUE)
              
              
              classlabels<-as.data.frame(classlabels)
              
              
              
              logp<-(-1)*log((d4[,1]+(10^-20)),10)
              
              
              if(fdrmethod=="none"){
                ythresh<-(-1)*log10(pvalue.thresh)
              }else{
                
                ythresh<-min(logp[goodip],na.rm=TRUE)
              }
              maintext1="Type 1 manhattan plot (-logp vs mz) \n m/z features above the dashed horizontal line meet the selection criteria"
              maintext2="Type 2 manhattan plot (-logp vs time) \n m/z features above the dashed horizontal line meet the selection criteria"
              
              
              
              if(is.na(zvec[1])==FALSE){
                maintext1=paste(maintext1,"\n",manhattanplot.col.opt[2],": lower in class ",class_labels_levels_main[1]," & ",manhattanplot.col.opt[1],": higher in class ",class_labels_levels_main[1],sep="")
                maintext2=paste(maintext2,"\n",manhattanplot.col.opt[2],": lower in class ",class_labels_levels_main[1]," & ",manhattanplot.col.opt[1],": higher in class ",class_labels_levels_main[1],sep="")
              }
              
              yvec_val=logp
              ylabel="(-)log10p"
              yincrement=1
              y2thresh=(-1)*log10(pvalue.thresh)
              
              
           # save(list=c("d4","logp","yvec_val","ythresh","zvec","x1increment","yincrement","maintext1","x2increment","maintext2","ylabel","y2thresh"),file="manhattanplot_objects.Rda")
              
              
              if(output.device.type!="pdf"){
                
                temp_filename_1<-"Figures/ManhattanPlot_Type1.png"
                
                png(temp_filename_1,width=plots.width,height=plots.height,res=plots.res,type=plots.type,units="in")
              }
              
              # get_manhattanplots(xvec=d4$mz,yvec=logp,ythresh=ythresh,up_or_down=zvec,xlab="mass-to-charge (m/z)",ylab=ylabel,xincrement=x1increment,yincrement=yincrement,maintext=maintext1,col_seq=c("black"),y2thresh=y2thresh,colorvec=manhattanplot.col.opt)
              
              ####savelist=ls(),file="m1.Rda")
              try(get_manhattanplots(xvec=d4$mz,yvec=logp,ythresh=ythresh,up_or_down=zvec,xlab="mass-to-charge (m/z)",ylab=ylabel,
                                     xincrement=x1increment,yincrement=yincrement,maintext=maintext1,col_seq=c("black"),y2thresh=y2thresh,colorvec=manhattanplot.col.opt),silent=TRUE)
              
              
              if(output.device.type!="pdf"){
                
                try(dev.off(),silent=TRUE)
              }
              
              if(output.device.type!="pdf"){
                
                temp_filename_1<-"Figures/ManhattanPlot_Type2.png"
                
                png(temp_filename_1,width=plots.width,height=plots.height,res=plots.res,type=plots.type,units="in")
              }
              
              try(get_manhattanplots(xvec=d4$time,yvec=logp,ythresh=ythresh,up_or_down=zvec,xlab="Retention time",ylab="-log10p",xincrement=x2increment,yincrement=1,maintext=maintext2,col_seq=c("black"),y2thresh=y2thresh,colorvec=manhattanplot.col.opt),silent=TRUE)
              
              if(output.device.type!="pdf"){
                
                try(dev.off(),silent=TRUE)
              }
              
              if(length(class_labels_levels)==2){
                if(output.device.type!="pdf"){
                  
                  temp_filename_1<-"Figures/VolcanoPlot.png"
                  
                  png(temp_filename_1,width=plots.width,height=plots.height,res=plots.res,type=plots.type,units="in")
                }
                
                maintext1="Volcano plot (-logp vs log2(fold change)) \n colored m/z features meet the selection criteria"
                if(is.na(zvec[1])==FALSE){
                  maintext1=paste(maintext1,"\n",manhattanplot.col.opt[2],": lower in class ",class_labels_levels_main[1]," & ",manhattanplot.col.opt[1],": higher in class ",class_labels_levels_main[1],sep="")
                  maintext2=paste(maintext2,"\n",manhattanplot.col.opt[2],": lower in class ",class_labels_levels_main[1]," & ",manhattanplot.col.opt[1],": higher in class ",class_labels_levels_main[1],sep="")
                }
                
                
                ##save(maxfoldchange,logp,zvec,ythresh,y2thresh,foldchangethresh,manhattanplot.col.opt,d4,file="debugvolcano.Rda")
                try(get_volcanoplots(xvec=maxfoldchange,yvec=logp,up_or_down=zvec,ythresh=ythresh,y2thresh=y2thresh,xthresh=foldchangethresh,maintext=maintext1,ylab="-log10(p-value)",xlab="log2(fold change)",colorvec=manhattanplot.col.opt),silent=TRUE)
                
                if(output.device.type!="pdf"){
                  
                  try(dev.off(),silent=TRUE)
                }
              }
              
            }else{
              
              if(featselmethod=="pls" | featselmethod=="o1pls"){
                
                
                # print("Time 2")
                #print(Sys.time())
                
                maintext1="Type 1 manhattan plot (VIP vs mz) \n m/z features above the dashed horizontal line meet the selection criteria"
                maintext2="Type 2 manhattan plot (VIP vs time) \n m/z features above the dashed horizontal line meet the selection criteria"
                
                if(is.na(zvec[1])==FALSE){
                  maintext1=paste(maintext1,"\n",manhattanplot.col.opt[2],": lower in class ",class_labels_levels_main[1]," & ",manhattanplot.col.opt[1],": higher in class ",class_labels_levels_main[1],sep="")
                  maintext2=paste(maintext2,"\n",manhattanplot.col.opt[2],": lower in class ",class_labels_levels_main[1]," & ",manhattanplot.col.opt[1],": higher in class ",class_labels_levels_main[1],sep="")
                }
                
                
                yvec_val<-data_limma_fdrall_withfeats[,1]
                ythresh=pls_vip_thresh
                vip_res<-as.data.frame(vip_res)
                bad.feature.index={}
                if(is.na(pls.permut.count)==FALSE){
                  #yvec_val[which(vip_res$rand_pls_sel_prob>=pvalue.thresh | vip_res$rand_pls_sel_fdr>=fdrthresh)]<-0 #(ythresh)*0.5
                  bad.feature.index=which(vip_res$rand_pls_sel_prob>=pvalue.thresh | vip_res$rand_pls_sel_fdr>=fdrthresh)
                }
                ylabel="VIP"
                yincrement=0.5
                y2thresh=NA
                
                # ###savelist=ls(),file="manhattandebug.Rda")
                
                if(output.device.type!="pdf"){
                  
                  temp_filename_1<-"Figures/ManhattanPlot_Type1.png"
                  png(temp_filename_1,width=plots.width,height=plots.height,res=plots.res,type=plots.type,units="in")
                  
                }
                
                try(get_manhattanplots(xvec=d4$mz,yvec=yvec_val,ythresh=pls_vip_thresh,up_or_down=zvec,xlab="mass-to-charge (m/z)",ylab="VIP",xincrement=x1increment,yincrement=0.5,maintext=maintext1,col_seq=c("black"),colorvec=manhattanplot.col.opt,bad.feature.index=bad.feature.index),silent=TRUE)
                
                if(output.device.type!="pdf"){
                  
                  try(dev.off(),silent=TRUE)
                }
                
                if(output.device.type!="pdf"){
                  
                  temp_filename_1<-"Figures/ManhattanPlot_Type2.png"
                  png(temp_filename_1,width=plots.width,height=plots.height,res=plots.res,type=plots.type,units="in")
                }
                
                try(get_manhattanplots(xvec=d4$time,yvec=yvec_val,ythresh=pls_vip_thresh,up_or_down=zvec,xlab="Retention time",ylab="VIP",xincrement=x2increment,yincrement=0.5,maintext=maintext2,col_seq=c("black"),colorvec=manhattanplot.col.opt,bad.feature.index=bad.feature.index),silent=TRUE)
                
                if(output.device.type!="pdf"){
                  
                  try(dev.off(),silent=TRUE)
                }
                
                if(length(class_labels_levels)==2){
                  
                  if(output.device.type!="pdf"){
                    temp_filename_1<-"Figures/VolcanoPlot_VIP_vs_foldchange.png"
                    
                    png(temp_filename_1,width=plots.width,height=plots.height,res=plots.res,type=plots.type,units="in")
                  }
                  
                  maintext1="Volcano plot (VIP vs log2(fold change)) \n colored m/z features meet the selection criteria"
                  
                  maintext1=paste(maintext1,"\n",manhattanplot.col.opt[2],": lower in class ",class_labels_levels_main[1]," & ",manhattanplot.col.opt[1],": higher in class ",class_labels_levels_main[1],sep="")
                  
                  maintext2=paste(maintext2,"\n",manhattanplot.col.opt[2],": lower in class ",class_labels_levels_main[1]," & ",manhattanplot.col.opt[1],": higher in class ",class_labels_levels_main[1],sep="")
                  
                  #  ###savelist=ls(),file="volcanodebug.Rda")
                  try(get_volcanoplots(xvec=maxfoldchange,yvec=yvec_val,up_or_down=maxfoldchange,ythresh=ythresh,xthresh=foldchangethresh,maintext=maintext1,ylab="VIP",xlab="log2(fold change)",bad.feature.index=bad.feature.index,colorvec=manhattanplot.col.opt),silent=TRUE)
                  
                  if(output.device.type!="pdf"){
                    
                    try(dev.off(),silent=TRUE)
                  }
                }
                
                
              }else{
                if(featselmethod=="spls" | featselmethod=="o1spls"){
                  
                  
                  maintext1="Type 1 manhattan plot (|loading| vs mz) \n m/z features with non-zero loadings meet the selection criteria"
                  maintext2="Type 2 manhattan plot (|loading| vs time) \n m/z features with non-zero loadings meet the selection criteria"
                  if(is.na(zvec[1])==FALSE){
                    maintext1=paste(maintext1,"\n",manhattanplot.col.opt[2],": lower in class ",class_labels_levels_main[1]," & ",manhattanplot.col.opt[1],": higher in class ",class_labels_levels_main[1],sep="")
                    maintext2=paste(maintext2,"\n",manhattanplot.col.opt[2],": lower in class ",class_labels_levels_main[1]," & ",manhattanplot.col.opt[1],": higher in class ",class_labels_levels_main[1],sep="")
                  }					
                  yvec_val<-data_limma_fdrall_withfeats[,1]
                  vip_res<-as.data.frame(vip_res)
                  
                  bad.feature.index={}
                  if(is.na(pls.permut.count)==FALSE){
                    # yvec_val[which(vip_res$rand_pls_sel_prob>=pvalue.thresh | vip_res$rand_pls_sel_fdr>=fdrthresh)]<-0
                    bad.feature.index=which(vip_res$rand_pls_sel_prob>=pvalue.thresh | vip_res$rand_pls_sel_fdr>=fdrthresh)
                  }
                  ythresh=0
                  ylabel="Loading (absolute)"
                  yincrement=0.1
                  y2thresh=NA
                  ####savelist=c("d4","yvec_val","ythresh","zvec","x1increment","yincrement","maintext1","x2increment","maintext2","ylabel","y2thresh"),file="manhattanplot_objects.Rda")
                  
                  if(output.device.type!="pdf"){
                    
                    temp_filename_1<-"Figures/ManhattanPlot_Type1.png"
                    
                    png(temp_filename_1,width=plots.width,height=plots.height,res=plots.res,type=plots.type,units="in")
                  }
                  
                  try(get_manhattanplots(xvec=d4$mz,yvec=yvec_val,ythresh=0,up_or_down=zvec,xlab="mass-to-charge (m/z)",ylab="Loading (absolute)",xincrement=x1increment,yincrement=0.1,maintext=maintext1,col_seq=c("black"),colorvec=manhattanplot.col.opt,bad.feature.index=bad.feature.index),silent=TRUE)
                  
                  if(output.device.type!="pdf"){
                    
                    try(dev.off(),silent=TRUE)
                  }
                  
                  if(output.device.type!="pdf"){
                    
                    temp_filename_1<-"Figures/ManhattanPlot_Type2.png"
                    
                    png(temp_filename_1,width=plots.width,height=plots.height,res=plots.res,type=plots.type,units="in")
                  }
                  
                  
                  try(get_manhattanplots(xvec=d4$time,yvec=yvec_val,ythresh=0,up_or_down=zvec,xlab="Retention time",ylab="Loading (absolute)",xincrement=x2increment,yincrement=0.1,maintext=maintext2,col_seq=c("black"),colorvec=manhattanplot.col.opt,bad.feature.index=bad.feature.index),silent=TRUE)
                  
                  
                  if(output.device.type!="pdf"){
                    
                    try(dev.off(),silent=TRUE)
                  }
                  
                  #volcanoplot
                  if(length(class_labels_levels)==2){
                    if(output.device.type!="pdf"){
                      
                      temp_filename_1<-"Figures/VolcanoPlot_Loading_vs_foldchange.png"
                      
                      png(temp_filename_1,width=plots.width,height=plots.height,res=plots.res,type=plots.type,units="in")
                    }
                    
                    maintext1="Volcano plot (absolute) Loading vs log2(fold change)) \n colored m/z features meet the selection criteria"
                    
                    maintext1=paste(maintext1,"\n",manhattanplot.col.opt[2],": lower in class ",class_labels_levels_main[1]," & ",manhattanplot.col.opt[1],": higher in class ",class_labels_levels_main[1],sep="")
                    maintext2=paste(maintext2,"\n",manhattanplot.col.opt[2],": lower in class ",class_labels_levels_main[1]," & ",manhattanplot.col.opt[1],": higher in class ",class_labels_levels_main[1],sep="")
                    
                    
                    
                    try(get_volcanoplots(xvec=maxfoldchange,yvec=yvec_val,up_or_down=maxfoldchange,ythresh=ythresh,xthresh=foldchangethresh,maintext=maintext1,ylab="(absolute) Loading",xlab="log2(fold change)",yincrement=0.1,bad.feature.index=bad.feature.index,colorvec=manhattanplot.col.opt),silent=TRUE)
                    
                    if(output.device.type!="pdf"){
                      
                      try(dev.off(),silent=TRUE)
                    }
                  }
                  
                  
                }else{
                  
                  if(featselmethod=="pamr"){
                    
                    
                    
                    
                    maintext1="Type 1 manhattan plot (max |standardized centroids (d-statistic)| vs mz) \n m/z features with above the horizontal line meet the selection criteria"
                    maintext2="Type 2 manhattan plot (max |standardized centroids (d-statistic)| vs time) \n m/z features with above the horizontal line meet the selection criteria"
                    if(is.na(zvec[1])==FALSE){
                      maintext1=paste(maintext1,"\n",manhattanplot.col.opt[2],": lower in class ",class_labels_levels_main[1]," & ",manhattanplot.col.opt[1],": higher in class ",class_labels_levels_main[1],sep="")
                      maintext2=paste(maintext2,"\n",manhattanplot.col.opt[2],": lower in class ",class_labels_levels_main[1]," & ",manhattanplot.col.opt[1],": higher in class ",class_labels_levels_main[1],sep="")
                    }
                    
                    
                    yvec_val<-data_limma_fdrall_withfeats[,1]
                    
                    ##error point
                    #vip_res<-as.data.frame(vip_res)
                    
                    discore<-as.data.frame(discore)
                    
                    bad.feature.index={}
                    if(is.na(pls.permut.count)==FALSE){
                      # yvec_val[which(vip_res$rand_pls_sel_prob>=pvalue.thresh | vip_res$rand_pls_sel_fdr>=fdrthresh)]<-0
                      #   bad.feature.index=which(vip_res$rand_pls_sel_prob>=pvalue.thresh | vip_res$rand_pls_sel_fdr>=fdrthresh)
                    }
                    ythresh=pamr_ythresh
                    ylabel="d-statistic (absolute)"
                    yincrement=0.1
                    y2thresh=NA
                    ####savelist=c("d4","yvec_val","ythresh","zvec","x1increment","yincrement","maintext1","x2increment","maintext2","ylabel","y2thresh"),file="manhattanplot_objects.Rda")
                    
                    if(output.device.type!="pdf"){
                      
                      temp_filename_1<-"Figures/ManhattanPlot_Type1.png"
                      
                      png(temp_filename_1,width=plots.width,height=plots.height,res=plots.res,type=plots.type,units="in")
                    }
                    
                    try(get_manhattanplots(xvec=d4$mz,yvec=yvec_val,ythresh=pamr_ythresh,up_or_down=zvec,xlab="mass-to-charge (m/z)",ylab="d-statistic (absolute) at threshold=0",xincrement=x1increment,yincrement=0.1,maintext=maintext1,col_seq=c("black"),colorvec=manhattanplot.col.opt,bad.feature.index=NA),silent=TRUE)
                    
                    if(output.device.type!="pdf"){
                      
                      try(dev.off(),silent=TRUE)
                    }
                    
                    if(output.device.type!="pdf"){
                      
                      temp_filename_1<-"Figures/ManhattanPlot_Type2.png"
                      
                      png(temp_filename_1,width=plots.width,height=plots.height,res=plots.res,type=plots.type,units="in")
                    }
                    
                    
                    try(get_manhattanplots(xvec=d4$time,yvec=yvec_val,ythresh=pamr_ythresh,up_or_down=zvec,xlab="Retention time",ylab="d-statistic (absolute) at threshold=0",xincrement=x2increment,yincrement=0.1,maintext=maintext2,col_seq=c("black"),colorvec=manhattanplot.col.opt,bad.feature.index=NA),silent=TRUE)
                    
                    
                    if(output.device.type!="pdf"){
                      
                      try(dev.off(),silent=TRUE)
                    }
                    
                    #volcanoplot
                    if(length(class_labels_levels)==2){
                      if(output.device.type!="pdf"){
                        
                        temp_filename_1<-"Figures/VolcanoPlot_Dstatistic_vs_foldchange.png"
                        
                        png(temp_filename_1,width=plots.width,height=plots.height,res=plots.res,type=plots.type,units="in")
                      }
                      
                      maintext1="Volcano plot (absolute) max standardized centroid (d-statistic) vs log2(fold change)) \n colored m/z features meet the selection criteria"
                      
                      maintext1=paste(maintext1,"\n",manhattanplot.col.opt[2],": lower in class ",class_labels_levels_main[1]," & ",manhattanplot.col.opt[1],": higher in class ",class_labels_levels_main[1],sep="")
                      maintext2=paste(maintext2,"\n",manhattanplot.col.opt[2],": lower in class ",class_labels_levels_main[1]," & ",manhattanplot.col.opt[1],": higher in class ",class_labels_levels_main[1],sep="")
                      
                      
                      
                      try(get_volcanoplots(xvec=maxfoldchange,yvec=yvec_val,up_or_down=maxfoldchange,ythresh=pamr_ythresh,xthresh=foldchangethresh,maintext=maintext1,ylab="(absolute) d-statistic at threshold=0",xlab="log2(fold change)",yincrement=0.1,bad.feature.index=NA,colorvec=manhattanplot.col.opt),silent=TRUE)
                      
                      if(output.device.type!="pdf"){
                        
                        try(dev.off(),silent=TRUE)
                      }
                    }
                    
                    
                  }
                  
                }
              }
              
            }
            
            
            
            
            goodip<-intersect(goodip,goodipfoldchange)
            
            
            dataA<-cbind(maxfoldchange,data_m_fc_withfeats)
            #write.table(dataA,file="foldchange.txt",sep="\t",row.names=FALSE)
            
            goodfeats_allfields<-{}
            
            
            if(length(goodip)>0){
              feat_sigfdrthresh[lf]<-length(goodip)
              subdata<-t(data_m_fc[goodip,])
              
              #save(parent_data_m,file="parent_data_m.Rda")
              
              data_minval<-min(parent_data_m[,-c(1:2)],na.rm=TRUE)*0.5
              
              #svm_model<-svm_cv(v=kfold,x=subdata,y=classlabels,kname=svm_kernel,errortype=pred.eval.method,conflevel=95)
              
              
              exp_fp<-1
              
              best_feats<-goodip
              
             
              
            }else{
              
              print("No features meet the fold change criteria.")
              
             
              
              
            }
            
          }else{
            
            if(dim(data_m_fc)[2]<kfold){
              print("Number of samples is too small to calculate cross-validation accuracy.")
            }
            
          }
          #feat_sigfdrthresh_cv<-c(feat_sigfdrthresh_cv,pred_acc)
          
          
          
          print("########################################")
          print(paste("Relative standard deviation (RSD) threshold: ", log2.fold.change.thresh," %",sep=""))
          #print(paste("FDR threshold: ", fdrthresh,sep=""))
          print(paste("Number of features left after RSD filtering: ", dim(data_m_fc)[1],sep=""))
          print(paste("Number of selected features: ", length(goodip),sep=""))
          
          if(length(goodip)<1){
            try(dev.off(),silent=TRUE)
            
            next
            
          }
          
          
          
          
         
          # save(data_m_fc_withfeats,data_matrix,data_m,goodip,names_with_mz_time,file="gdebug.Rda")
          
          
          print("######################################")
          suppressMessages(library(cluster))
          
          t1<-table(classlabels)
          
          
          if(is.na(names_with_mz_time)==FALSE){
            data_m_fc_withfeats_A1<-merge(names_with_mz_time,data_m_fc_withfeats,by=c("mz","time"))
            
            rownames(data_m_fc_withfeats)<-as.character(data_m_fc_withfeats_A1$Name)
            
            
          }else{
            
            rownames(data_m_fc_withfeats)<-as.character(paste(data_m_fc_withfeats[,1],data_m_fc_withfeats[,2],sep="_"))
          }
          
          
          #patientcolors <- unlist(lapply(sampleclass, color.map))
          if(length(goodip)>2){
            
            goodfeats<-as.data.frame(data_m_fc_withfeats[goodip,]) #[sel.diffdrthresh==TRUE,])
            
            goodfeats<-unique(goodfeats)
            
            
            rnames_goodfeats<-rownames(goodfeats) #as.character(paste(goodfeats[,1],goodfeats[,2],sep="_"))
            
            if(length(which(duplicated(rnames_goodfeats)==TRUE))>0){
              
              print("WARNING: Duplicated features found. Removing duplicate entries.")
              goodfeats<-goodfeats[-which(duplicated(rnames_goodfeats)==TRUE),]
              rnames_goodfeats<-rnames_goodfeats[-which(duplicated(rnames_goodfeats)==TRUE)]
            }
            
            
            #rownames(goodfeats)<-as.character(paste(goodfeats[,1],goodfeats[,2],sep="_"))
            
            
            
            data_m<-as.matrix(goodfeats[,-c(1:2)])
            
            
            rownames(data_m)<-rownames(goodfeats) #as.character(paste(goodfeats[,1],goodfeats[,2],sep="_"))
            
          
            
            
            data_m<-unique(data_m)			
            
            X<-t(data_m)
            
       
              {
              
              heatmap_file<-paste("heatmap_",featselmethod,".tiff",sep="")
              
              heatmap_mainlabel="" #2-way HCA using all significant features"
              
            if(FALSE)
              {
               # print("this step")
              #  save(hc,file="hc.Rda")
               # save(hr,file="hr.Rda")
                #save(distc,file="distc.Rda")
                #save(distr,file="distr.Rda")
              #  save(data_m,heatmap.col.opt,hca_type,classlabels,classlabels_orig,outloc,goodfeats,data_m_fc_withfeats,goodip,names_with_mz_time,plots.height,plots.width,plots.res,file="hcadata_m.Rda")
                #save(classlabels,file="classlabels.Rda")
              }
              
             
             # pdf("Testhca.pdf")
              
              #try(
               # 
              #dev.off()
              
              
              if(is.na(names_with_mz_time)==FALSE){
              
                goodfeats_with_names<-merge(names_with_mz_time,goodfeats,by=c("mz","time"))
                
                goodfeats_with_names<-goodfeats_with_names[match(paste(goodfeats$mz,"_",goodfeats$time,sep=""),paste(goodfeats_with_names$mz,"_",goodfeats_with_names$time,sep="")),]
                
            #    save(names_with_mz_time,goodfeats,goodfeats_with_names,file="goodfeats_with_names.Rda")
                
                goodfeats_name<-goodfeats_with_names$Name
                
                rownames(goodfeats)<-goodfeats_name
                
              }else{
                
                #print(head(names_with_mz_time))
               # print(head(goodfeats))
                #goodfeats_name<-NA
              }
              
              
          
              if(output.device.type!="pdf"){
                
              #  print(getwd())
                
        #  save(data_m,heatmap.col.opt,hca_type,classlabels,classlabels_orig,output_dir,goodfeats,names_with_mz_time,data_m_fc_withfeats,goodip,goodfeats_name,names_with_mz_time,
         #          plots.height,plots.width,plots.res,alphabetical.order,analysistype,labRow.value, labCol.value,hca.cex.legend,file="hcadata_mD.Rda")
                
              
                temp_filename_1<-"Figures/HCA_All_selectedfeats.png"
                png(temp_filename_1,width=plots.width,height=plots.height,res=plots.res,type="cairo",units="in")
                
                
                #Generate HCA for selected features
                hca_res<-get_hca(feature_table_file=NA,parentoutput_dir=output_dir,class_labels_file=NA,X=goodfeats,Y=classlabels_orig,heatmap.col.opt=heatmap.col.opt,
                                 cor.method=cor.method,is.data.znorm=FALSE,analysismode="classification",
                                 sample.col.opt=sample.col.opt,plots.width=2000,plots.height=2000,plots.res=300, alphacol=0.3, hca_type=hca_type,newdevice=FALSE,
                                 input.type="intensity",mainlab="",alphabetical.order=alphabetical.order,study.design=analysistype,
                                 labRow.value = labRow.value, labCol.value = labCol.value,similarity.matrix=similarity.matrix,cexLegend=hca.cex.legend)
                
                dev.off()
                
              }else{
                
                #Generate HCA for selected features
                hca_res<-get_hca(feature_table_file=NA,parentoutput_dir=output_dir,class_labels_file=NA,X=goodfeats,Y=classlabels_orig,heatmap.col.opt=heatmap.col.opt,cor.method=cor.method,is.data.znorm=FALSE,analysismode="classification",
                                 sample.col.opt=sample.col.opt,plots.width=2000,plots.height=2000,plots.res=300, alphacol=0.3, hca_type=hca_type,newdevice=FALSE,
                                 input.type="intensity",mainlab="",alphabetical.order=alphabetical.order,study.design=analysistype,
                                 labRow.value = labRow.value, labCol.value = labCol.value,similarity.matrix=similarity.matrix,cexLegend=hca.cex.legend)
                
                  
                 # get_hca(parentoutput_dir=getwd(),X=goodfeats,Y=classlabels_orig,heatmap.col.opt=heatmap.col.opt,cor.method="spearman",is.data.znorm=FALSE,analysismode="classification",
                       # sample.col.opt="rainbow",plots.width=2000,plots.height=2000,plots.res=300, alphacol=0.3, hca_type=hca_type,newdevice=FALSE) #,silent=TRUE)
                
              }
              
             
              
         
              
             
            }
            
            
            print("Done with HCA.")
            
          }
         
          
          
        }
        
        else
        {
          
          
          #print("regression")
          
          print("########################################")
          print(paste("RSD threshold: ", log2.fold.change.thresh,sep=""))
          #print(paste("FDR threshold: ", fdrthresh,sep=""))
          print(paste("Number of metabolites left after RSD filtering: ", dim(data_m_fc)[1],sep=""))
          print(paste("Number of sig metabolites: ", length(goodip),sep=""))
          
          
          
          if(featselmethod=="lmreg"){
            #d4<-read.table(paste(parentoutput_dir,"/Stage2/lmreg_pval_coef_stderr.txt",sep=""),sep="\t",header=TRUE,quote = "")
            
            d4<-read.table("Tables/lmreg_pval_coef_stderr.txt",sep="\t",header=TRUE)
            
          }
          
          
          
          if(length(goodip)>=1){
            
            
            subdata<-t(data_m_fc[goodip,])
            
            
            if(length(class_labels_levels)==2){
              
              
              #zvec=foldchangeres
            }else{
              
              zvec=NA
              
              if(featselmethod=="lmreg" && analysismode=="regression"){
                
                cnames_matrix<-colnames(d4)
                
                
                
                cnames_colindex<-grep("Estimate_",cnames_matrix)
                
                
                zvec<-d4[,c(cnames_colindex[1])]
                #zvec<-d4$Estimate_var1
                
                #if(length(zvec)<1){
                #   zvec<-d4$X.Estimate_var1.
                
                #}
              }
              
              
            }
            
            
            
            
            
            
            roundUpNice <- function(x, nice=c(1,2,4,5,6,8,10)) {
              if(length(x) != 1) stop("'x' must be of length 1")
              10^floor(log10(x)) * nice[[which(x <= 10^floor(log10(x)) * nice)[[1]]]]
            }
            
            d4<-as.data.frame(data_limma_fdrall_withfeats)
            # d4<-as.data.frame(d1)
            
            
           # save(d4,file="mtype1.rda")
            x1increment=round_any(max(d4$mz)/10,10,f=floor)
            x2increment=round_any(max(d4$time)/10,10,f=floor)
            
            #manplots
            
            if(featselmethod=="lmreg" | featselmethod=="lm1wayanova" | featselmethod=="lm2wayanova" | featselmethod=="lm1wayanovarepeat" | featselmethod=="lm2wayanovarepeat"
               | featselmethod=="limma" | featselmethod=="limma2way" | featselmethod=="logitreg" | featselmethod=="limma2wayrepeat" | featselmethod=="wilcox" | featselmethod=="ttest" | featselmethod=="poissonreg" | featselmethod=="lmregrepeat")
            {
              
              
              
              #print("Plotting manhattan plots")
              
              sel.diffdrthresh<-fdr_adjust_pvalue<fdrthresh & final.pvalues<pvalue.thresh
              
              
              
              goodip<-which(sel.diffdrthresh==TRUE)
              
              
              
              
              classlabels<-as.data.frame(classlabels)
              
              
              
              logp<-(-1)*log((d4[,1]+(10^-20)),10)
              ythresh<-min(logp[goodip],na.rm=TRUE)
              
              maintext1="Type 1 manhattan plot (-logp vs mz) \n m/z features above the dashed horizontal line meet the selection criteria"
              maintext2="Type 2 manhattan plot (-logp vs time) \n m/z features above the dashed horizontal line meet the selection criteria"
              
              
              # print("here1 A")
              #print(zvec)
              
              if(is.na(zvec[1])==FALSE){
                maintext1=paste(maintext1,"\n",manhattanplot.col.opt[2],": negative association "," & ",manhattanplot.col.opt[1],": positive association ",sep="")
                maintext2=paste(maintext2,"\n",manhattanplot.col.opt[2],": negative association "," & ",manhattanplot.col.opt[1],": positive association ",sep="")
              }
              
              
              if(output.device.type!="pdf"){
                
                temp_filename_1<-"Figures/ManhattanPlot_Type1.png"
                
                png(temp_filename_1,width=plots.width,height=plots.height,res=plots.res,type=plots.type,units="in")
              }
              
             
              try(get_manhattanplots(xvec=d4$mz,yvec=logp,ythresh=ythresh,up_or_down=zvec,xlab="mass-to-charge (m/z)",ylab="-logP",xincrement=x1increment,yincrement=1,
                                     maintext=maintext1,col_seq=c("black"),y2thresh=(-1)*log10(pvalue.thresh),colorvec=manhattanplot.col.opt),silent=TRUE)
              
              
              
              if(output.device.type!="pdf"){
                
                try(dev.off(),silent=TRUE)
              }
              
              
              if(output.device.type!="pdf"){
                
                temp_filename_1<-"Figures/ManhattanPlot_Type2.png"
                
                png(temp_filename_1,width=plots.width,height=plots.height,res=plots.res,type=plots.type,units="in")
              }
              
              
              
              try(get_manhattanplots(xvec=d4$time,yvec=logp,ythresh=ythresh,up_or_down=zvec,xlab="Retention time",ylab="-logP",xincrement=x2increment,yincrement=1,
                                     maintext=maintext2,col_seq=c("black"),y2thresh=(-1)*log10(pvalue.thresh),colorvec=manhattanplot.col.opt),silent=TRUE)
              
              
              if(output.device.type!="pdf"){
                
                try(dev.off(),silent=TRUE)
              }
              
              
              #print("Plotting manhattan plots")
              #get_manhattanplots(xvec=d4$mz,yvec=logp,ythresh=ythresh,up_or_down=zvec,xlab="mass-to-charge (m/z)",ylab="-logP",xincrement=x1increment,yincrement=1,maintext=maintext1)
              #get_manhattanplots(xvec=d4$time,yvec=logp,ythresh=ythresh,up_or_down=zvec,xlab="Retention time",ylab="-logP",xincrement=x2increment,yincrement=1,maintext=maintext2)
              
            }else{
              
              if(featselmethod=="pls" | featselmethod=="o1pls"){
                
                maintext1="Type 1 manhattan plot (VIP vs mz) \n m/z features above the dashed horizontal line meet the selection criteria"
                maintext2="Type 2 manhattan plot (VIP vs time) \n m/z features above the dashed horizontal line meet the selection criteria"
                
                if(is.na(zvec[1])==FALSE){
                  maintext1=paste(maintext1,"\n",manhattanplot.col.opt[2],": negative association "," & ",manhattanplot.col.opt[1],": positive association ",sep="")
                  maintext2=paste(maintext2,"\n",manhattanplot.col.opt[2],": negative association "," & ",manhattanplot.col.opt[1],": positive association ",sep="")
                }
                
                
                if(output.device.type!="pdf"){
                  
                  temp_filename_1<-"Figures/ManhattanPlot_Type1.png"
                  
                  png(temp_filename_1,width=plots.width,height=plots.height,res=plots.res,type=plots.type,units="in")
                }
                
                
                
                
                try(get_manhattanplots(xvec=d4$mz,yvec=data_limma_fdrall_withfeats[,1],ythresh=pls_vip_thresh,up_or_down=zvec,xlab="mass-to-charge (m/z)",ylab="VIP",xincrement=x1increment,yincrement=0.5,maintext=maintext1,col_seq=c("black"),colorvec=manhattanplot.col.opt),silent=TRUE)
                
                
                
                if(output.device.type!="pdf"){
                  
                  try(dev.off(),silent=TRUE)
                }
                
                
                if(output.device.type!="pdf"){
                  
                  temp_filename_1<-"Figures/ManhattanPlot_Type2.png"
                  
                  png(temp_filename_1,width=plots.width,height=plots.height,res=plots.res,type=plots.type,units="in")
                }
                
                
                
                try(get_manhattanplots(xvec=d4$time,yvec=data_limma_fdrall_withfeats[,1],ythresh=pls_vip_thresh,up_or_down=zvec,xlab="Retention time",ylab="VIP",xincrement=x2increment,yincrement=0.5,maintext=maintext2,col_seq=c("black"),colorvec=manhattanplot.col.opt),silent=TRUE)
                
                
                
                if(output.device.type!="pdf"){
                  
                  try(dev.off(),silent=TRUE)
                }
                
              }
              else{
                if(featselmethod=="spls" | featselmethod=="o1spls"){
                  
                  
                  maintext1="Type 1 manhattan plot (|loading| vs mz) \n m/z features with non-zero loadings meet the selection criteria"
                  maintext2="Type 2 manhattan plot (|loading| vs time) \n m/z features with non-zero loadings meet the selection criteria"
                  if(is.na(zvec[1])==FALSE){
                    maintext1=paste(maintext1,"\n",manhattanplot.col.opt[2],": negative association "," & ",manhattanplot.col.opt[1],": positive association ",sep="")
                    maintext2=paste(maintext2,"\n",manhattanplot.col.opt[2],": negative association "," & ",manhattanplot.col.opt[1],": positive association ",sep="")
                  }
                  
                  
                  if(output.device.type!="pdf"){
                    
                    temp_filename_1<-"Figures/ManhattanPlot_Type1.png"
                    
                    png(temp_filename_1,width=plots.width,height=plots.height,res=plots.res,type=plots.type,units="in")
                  }
                  
                  
                  try(get_manhattanplots(xvec=d4$mz,yvec=data_limma_fdrall_withfeats[,1],ythresh=0,up_or_down=zvec,xlab="mass-to-charge (m/z)",ylab="Loading",xincrement=x1increment,yincrement=0.1,maintext=maintext1,col_seq=c("black"),colorvec=manhattanplot.col.opt),silent=TRUE)
                  
                  
                  if(output.device.type!="pdf"){
                    
                    try(dev.off(),silent=TRUE)
                  }
                  
                  
                  if(output.device.type!="pdf"){
                    
                    temp_filename_1<-"Figures/ManhattanPlot_Type2.png"
                    
                    png(temp_filename_1,width=plots.width,height=plots.height,res=plots.res,type=plots.type,units="in")
                  }
                  
                  
                  
                  try(get_manhattanplots(xvec=d4$time,yvec=data_limma_fdrall_withfeats[,1],ythresh=0,up_or_down=zvec,xlab="Retention time",ylab="Loading",xincrement=x2increment,yincrement=0.1,maintext=maintext2,col_seq=c("black"),colorvec=manhattanplot.col.opt),silent=TRUE)
                  
                  
                  if(output.device.type!="pdf"){
                    
                    try(dev.off(),silent=TRUE)
                  }
                  
                }
              }
              
            }
            
            
            data_minval<-min(parent_data_m[,-c(1:2)],na.rm=TRUE)*0.5
            #subdata<-apply(subdata,2,function(x){naind<-which(is.na(x)==TRUE);if(length(naind)>0){ x[naind]<-median(x,na.rm=TRUE)};return(x)})
            subdata<-apply(subdata,2,function(x){naind<-which(is.na(x)==TRUE);if(length(naind)>0){ x[naind]<-data_minval};return(x)})
            
            
            #print(head(subdata))
            #print(dim(subdata))
            #print(dim(classlabels))
            
            #print(dim(classlabels))
            
            classlabels_response_mat<-as.data.frame(classlabels_response_mat)
            
            if(length(classlabels)>dim(parent_data_m)[2]){
              #classlabels<-as.data.frame(classlabels[,1])
              classlabels_response_mat<-as.data.frame(classlabels_response_mat[,1])
            }
            
            
            svm_model_reg<-try(svm(x=subdata,y=(classlabels_response_mat[,1]),type="eps",cross=kfold),silent=TRUE)
            
            if(is(svm_model_reg,"try-error")){
              print("SVM could not be performed. Skipping to the next step.")
              termA<-(-1)
              pred_acc<-termA
            }else{
              termA<-svm_model_reg$tot.MSE
              
              pred_acc<-termA
              print(paste(kfold,"-fold mean squared error: ", pred_acc,sep=""))
              
            }
            
            
            
            print("######################################")
          }else{
            print("Number of selected variables is too small to perform CV.")
            
          }
          
          #print("termA is ")
          #print(termA)
          
          # print("dim of goodfeats")
          goodfeats<-as.data.frame(data_m_fc_withfeats[sel.diffdrthresh==TRUE,])
          
          goodip<-which(sel.diffdrthresh==TRUE)
          
          #print(length(goodip))
          
          res_score<-termA
          
          #if(res_score<best_cv_res){
          
          if(length(which(sel.diffdrthresh==TRUE))>0){
            if(res_score<best_cv_res){
              
              
              best_logfc_ind<-lf
              
              best_feats<-goodip
              best_cv_res<-res_score
              best_acc<-pred_acc
              best_limma_res<-data_limma_fdrall_withfeats[sel.diffdrthresh==TRUE,]
              
            }
          }else{
            res_score<-(9999999)
          }
          res_score_vec[lf]<-res_score
   
          goodfeats<-unique(goodfeats)
          
          if(is.na(names_with_mz_time)==FALSE){
            goodfeats_with_names<-merge(names_with_mz_time,goodfeats,by=c("mz","time"))
            goodfeats_with_names<-goodfeats_with_names[match(goodfeats$mz,goodfeats_with_names$mz),]
            
            #save(names_with_mz_time,goodfeats,goodfeats_with_names,file="goodfeats_with_names.Rda")
            
            goodfeats_name<-goodfeats_with_names$Name
            #}
          }else{
            
            
            goodfeats_name<-as.character(paste(goodfeats[,1],goodfeats[,2],sep="_"))
          }
          
          
          if(length(which(sel.diffdrthresh==TRUE))>2){
            
            
            
            
            ##save(goodfeats,file="goodfeats.Rda")
            #rownames(goodfeats)<-as.character(goodfeats[,1])
            rownames(goodfeats)<-goodfeats_name #as.character(paste(goodfeats[,1],goodfeats[,2],sep="_"))
            
            
            data_m<-as.matrix(goodfeats[,-c(1:2)])
            
            rownames(data_m)<-rownames(goodfeats) #as.character(paste(goodfeats[,1],goodfeats[,2],sep="_"))
            
            
            X<-t(data_m)
            
            pca_comp<-min(dim(X)[1],dim(X)[2])
            
            t1<-seq(1,dim(data_m)[2])
            
            col <-col_vec[1:length(t1)]
            
            
            
            hr <- try(hclust(as.dist(1-cor(t(data_m),method=cor.method,use="pairwise.complete.obs"))),silent=TRUE) #metabolites
            hc <- try(hclust(as.dist(1-cor(data_m,method=cor.method,use="pairwise.complete.obs"))),silent=TRUE) #samples
            
            
            if(heatmap.col.opt=="RdBu"){
              
              heatmap.col.opt="redblue"
            }
            
            heatmap_cols <- colorRampPalette(brewer.pal(10, "RdBu"))(256)
            heatmap_cols<-rev(heatmap_cols)
            
            if(heatmap.col.opt=="topo"){
              heatmap_cols<-topo.colors(256)
              heatmap_cols<-rev(heatmap_cols)
            }else{
              if(heatmap.col.opt=="heat"){
                heatmap_cols<-heat.colors(256)
                heatmap_cols<-rev(heatmap_cols)
              }else{
                
                if(heatmap.col.opt=="yellowblue"){
                  
                  heatmap_cols<-colorRampPalette(c("yellow","blue"))(256) #colorRampPalette(c("yellow","white","blue"))(256)
                  #heatmap_cols<-blue2yellow(256) #colorRampPalette(c("yellow","blue"))(256)
                  heatmap_cols<-rev(heatmap_cols)
                }else{
                  
                  if(heatmap.col.opt=="redblue"){
                    
                    heatmap_cols <- colorRampPalette(brewer.pal(10, "RdBu"))(256)
                    heatmap_cols<-rev(heatmap_cols)
                  }else{
                    
                    #my_palette <- colorRampPalette(c("red", "yellow", "green"))(n = 299)
                    if(heatmap.col.opt=="redyellowgreen"){
                      
                      heatmap_cols <- colorRampPalette(c("red", "yellow", "green"))(n = 299)
                      heatmap_cols<-rev(heatmap_cols)
                    }else{
                      if(heatmap.col.opt=="yellowwhiteblue"){
                        
                        heatmap_cols<-colorRampPalette(c("yellow2","white","blue"))(256) #colorRampPalette(c("yellow","white","blue"))(256)
                        heatmap_cols<-rev(heatmap_cols)
                      }else{
                        
                        if(heatmap.col.opt=="redwhiteblue"){
                          
                          heatmap_cols<-colorRampPalette(c("red","white","blue"))(256) #colorRampPalette(c("yellow","white","blue"))(256)
                          heatmap_cols<-rev(heatmap_cols)
                        }else{
                          
                          
                          
                          heatmap_cols <- colorRampPalette(brewer.pal(10, heatmap.col.opt))(256)
                          heatmap_cols<-rev(heatmap_cols)
                          
                        }
                        
                      }
                      
                    }
                    
                  }
                  
                }
              }
              
            }
            
            
            
            if(is(hr,"try-error") || is(hc,"try-error")){
              
              print("Hierarchical clustering can not be performed. ")
            }else{
              
              
              mycl_samples <- cutree(hc, h=max(hc$height)/2)
              t1<-table(mycl_samples)
              col_clust<-topo.colors(length(t1))
              patientcolors=rep(col_clust,t1) #mycl_samples[col_clust]
              heatmap_file<-paste("heatmap_",featselmethod,"_imp_features.tiff",sep="")
              
              #tiff(heatmap_file,width=plots.width,height=plots.height,res=plots.res, compression="lzw")
              
              
              if(output.device.type!="pdf"){
                
                temp_filename_1<-"Figures/HCA_all_selectedfeats.png"
                
                png(temp_filename_1,width=plots.width,height=plots.height,res=plots.res,type=plots.type,units="in")
              }
              
              
              
              
              if(znormtransform==FALSE){
                h73<-heatmap.2(data_m, Rowv=as.dendrogram(hr), Colv=as.dendrogram(hc),  col=heatmap_cols, scale="row",key=TRUE, symkey=FALSE, 
                               density.info="none", trace="none", cexRow=1, cexCol=1,xlab="",ylab="", main="Using all selected features",labRow = hca.labRow.value, labCol = hca.labCol.value)
              }else{
                h73<-heatmap.2(data_m, Rowv=as.dendrogram(hr), Colv=as.dendrogram(hc),  col=heatmap_cols, scale="none",key=TRUE, 
                               symkey=FALSE, density.info="none", trace="none", cexRow=1, cexCol=1,xlab="",ylab="", main="Using all selected features",labRow = hca.labRow.value, labCol = hca.labCol.value)
              }
              
              
              if(output.device.type!="pdf"){
                
                try(dev.off(),silent=TRUE)
              }
              
              
              
              mycl_samples <- cutree(hc, h=max(hc$height)/2)
              mycl_metabs <- cutree(hr, h=max(hr$height)/2)
              
              ord_data<-cbind(mycl_metabs[rev(h73$rowInd)],goodfeats[rev(h73$rowInd),c(1:2)],data_m[rev(h73$rowInd),h73$colInd])
              
              cnames1<-colnames(ord_data)
              cnames1[1]<-"mz_cluster_label"
              colnames(ord_data)<-cnames1
              fname1<-paste("Tables/Clustering_based_sorted_intensity_data.txt",sep="")
              write.table(ord_data,file=fname1,sep="\t",row.names=FALSE)
              
              fname2<-paste("Tables/Sample_clusterlabels.txt",sep="")
              
              sample_clust_num<-mycl_samples[h73$colInd]
              
              
              
              classlabels<-as.data.frame(classlabels)
              
              
              temp1<-classlabels[h73$colInd,]
              
              temp3<-cbind(temp1,sample_clust_num)
              
              rnames1<-rownames(temp3)
              temp4<-cbind(rnames1,temp3)
              temp4<-as.data.frame(temp4)
              
              
              if(analysismode=="regression"){
                
                
                #names(temp3[,1)<-as.character(temp4[,1])
                
                
                
                temp3<-temp4[,-c(1)]
                temp3<-as.data.frame(temp3)
                temp3<-apply(temp3,2,as.numeric)
                
                
                temp_vec<-as.vector(temp3[,1])
                
                
                
                names(temp_vec)<-as.character(temp4[,1])
                
                
                if(output.device.type!="pdf"){
                  
                  temp_filename_1<-"Figures/Barplot_dependent_variable_ordered_by_HCA.png"
                  
                  png(temp_filename_1,width=plots.width,height=plots.height,res=plots.res,type=plots.type,units="in")
                }
                
                
                
                #tiff("Barplot_sample_cluster_ymat.tiff", width=plots.width,height=plots.height,res=plots.res, compression="lzw")
                barplot(temp_vec,col="brown",ylab="Y",cex.axis=0.5,cex.names=0.5,main="Dependent variable levels in samples; \n ordered based on hierarchical clustering")
                #dev.off()
                
                
                if(output.device.type!="pdf"){
                  
                  try(dev.off(),silent=TRUE)
                }
                
                
                
              }
              
              #	print(head(temp_vec))
              #temp4<-temp4[,-c(2)]
              write.table(temp4,file=fname2,sep="\t",row.names=FALSE)
              
              
              
              fname3<-paste("Metabolite_clusterlabels.txt",sep="")
              
              mycl_metabs_ord<-mycl_metabs[rev(h73$rowInd)]
              
            }
          }
          
          
          
          
          
        }
        
        classlabels_orig<-classlabels_orig_parent
        if(pairedanalysis==TRUE){
          
          classlabels_orig<-classlabels_orig[,-c(2)]
          
          
        }else{
          
          if(featselmethod=="lmreg" || featselmethod=="logitreg" || featselmethod=="poissonreg"){
            classlabels_orig<-classlabels_orig[,c(1:2)]
            classlabels_orig<-as.data.frame(classlabels_orig)
          }
        }
        
        
        node_names=rownames(data_m_fc_withfeats)
        #save(data_limma_fdrall_withfeats,goodip,data_m_fc_withfeats,data_matrix,names_with_mz_time,file="data_limma_fdrall_withfeats1.Rda")
       
        classlabels_orig_wgcna<-classlabels_orig
        
        if(analysismode=="classification"){
          
          classlabels_temp<-classlabels_orig_wgcna #cbind(classlabels_sub[,1],classlabels)
          
          sigfeats=data_m_fc_withfeats[goodip,c(1:2)]
         
          # save(data_m_fc_withfeats,classlabels_temp,sigfeats,goodip,num_nodes,abs.cor.thresh,cor.fdrthresh,alphabetical.order,
           #    plot_DiNa_graph,degree.centrality.method,node_names,networktype,file="debugdiffrank_eval.Rda")
          
          if(degree_rank_method=="diffrank"){
           # degree_eval_res<-try(diffrank_eval(X=data_m_fc_withfeats,Y=classlabels_temp,sigfeats=data_m_fc_withfeats[goodip,c(1:2)],sigfeatsind=goodip,
            #                                   num_nodes=num_nodes,abs.cor.thresh=abs.cor.thresh,cor.fdrthresh=cor.fdrthresh,alphabetical.order=alphabetical.order),silent=TRUE)
        
            degree_eval_res<-diffrank_eval(X=data_m_fc_withfeats,Y=classlabels_temp,sigfeats=sigfeats,sigfeatsind=goodip,
                                               num_nodes=num_nodes,abs.cor.thresh=abs.cor.thresh,cor.fdrthresh=cor.fdrthresh,alphabetical.order=alphabetical.order,
                                           node_names=node_names,plot_graph_bool=plot_DiNa_graph,
                                           degree.centrality.method = degree.centrality.method,networktype=networktype) #,silent=TRUE)
            
            
            
          }else{
            
            degree_eval_res<-{}
          }
          
        }
        
       
        
        sample_names_vec<-colnames(data_m_fc_withfeats[,-c(1:2)])
        
      #  save(degree_eval_res,file="DiNa_results.Rda")
        
       # save(data_limma_fdrall_withfeats,goodip,sample_names_vec,data_m_fc_withfeats,data_matrix,names_with_mz_time,file="data_limma_fdrall_withfeats.Rda")
        
        if(analysismode=="classification")
        {
          
          degree_rank<-rep(1,dim(data_m_fc_withfeats)[1])
          
          if(is(degree_eval_res,"try-error")){
            
            degree_rank<-rep(1,dim(data_m_fc_withfeats)[1])
            
            
          }else{
            if(degree_rank_method=="diffrank"){
              
              diff_degree_measure<-degree_eval_res$all
              
              degree_rank<-diff_degree_measure$DiffRank #rank((1)*diff_degree_measure)
            }
            
            
          }
          
        #  save(degree_rank,file="degree_rank.Rda")
          
          
          if(featselmethod=="lmreg" | featselmethod=="limma" | featselmethod=="limma2way" | featselmethod=="limma1way" | featselmethod=="lmreg" | featselmethod=="logitreg" | featselmethod=="limma1wayrepeat" | featselmethod=="limma2wayrepeat" | featselmethod=="lm1wayanova" | featselmethod=="lm2wayanova" | featselmethod=="lm1wayanovarepeat" | featselmethod=="lm2wayanovarepeat" | featselmethod=="wilcox" | featselmethod=="ttest" | featselmethod=="poissonreg" | featselmethod=="lmregrepeat")
          {
            diffexp_rank<-rank(data_limma_fdrall_withfeats[,2]) #order(data_limma_fdrall_withfeats[,2],decreasing=FALSE)
            
            type.statistic="pvalue"
            x1=Sys.time()
            stat_val<-(-1)*log10(data_limma_fdrall_withfeats[,2])
            
            if(output.device.type!="pdf"){
              pdf("Figures/pvalue.distribution.pdf",width=10,height=8)
            }
            par(mfrow=c(1,2))
            kstest_res<-ks.test(data_limma_fdrall_withfeats[,2],"punif",0,1)
            kstest_res<-round(kstest_res$p.value,3)
            
            hist(as.numeric(data_limma_fdrall_withfeats[,2]),main=paste("Distribution of pvalues\n","Kolmogorov-Smirnov test for uniform distribution, p=",kstest_res,sep=""),cex.main=0.75,xlab="pvalues")
            #if(FALSE)
            {
              simpleQQPlot = function (observedPValues,mainlab) {
                plot(-log10(1:length(observedPValues)/length(observedPValues)),
                     -log10(sort(observedPValues)),main=mainlab,xlab=paste("Expected -log10pvalue",sep=""),ylab=paste("Observed -logpvalue",sep=""),cex.main=0.75)
                abline(0, 1, col = "brown")
              }
              
              inflation <- function(pvalue) {
                chisq <- qchisq(1 - pvalue, 1)
                lambda <- median(chisq) / qchisq(0.5, 1)
                lambda
              }
              inflation_res<-round(inflation(data_limma_fdrall_withfeats[,2]),2)
              
              simpleQQPlot(data_limma_fdrall_withfeats[,2],mainlab=paste("QQplot pvalues","\np-value inflation factor: ",inflation_res," (no inflation: close to 1; bias: greater than 1)",sep=""))
              
              x2=Sys.time()
              print(x2-x1)
            }
            
            
            if(output.device.type!="pdf"){
              dev.off()
            }
            
            par(mfrow=c(1,1))
            
            
            
          }else{
            
            
            if(featselmethod=="rfesvm"){
              
              
              diffexp_rank<-rank((1)*abs(data_limma_fdrall_withfeats[,2]))
              #diffexp_rank<-rank_vec
              #data_limma_fdrall_withfeats<-cbind(rank_vec,data_limma_fdrall_withfeats)
              
              
            }else{
              
              
              if(featselmethod=="pamr"){
                
                diffexp_rank<-rank_vec
                #data_limma_fdrall_withfeats<-cbind(rank_vec,data_limma_fdrall_withfeats[,-c(1)])
                
                
              }else{
                
                if(featselmethod=="MARS"){
                  
                  diffexp_rank<-rank((-1)*data_limma_fdrall_withfeats[,2])
                  
                }else{
                  diffexp_rank<-rank((1)*data_limma_fdrall_withfeats[,2])
                }
              }
            }
          }
          
          if(input.intensity.scale=="raw" && log2transform==FALSE){
            
            
            
            fold.change.log2<-maxfoldchange
            data_limma_fdrall_withfeats_2<-cbind(fold.change.log2,degree_rank,diffexp_rank,data_limma_fdrall_withfeats)
            
            
          }else{
            
            if(input.intensity.scale=="log2" || log2transform==TRUE){
              
              fold.change.log2<-maxfoldchange
              data_limma_fdrall_withfeats_2<-cbind(fold.change.log2,degree_rank,diffexp_rank,data_limma_fdrall_withfeats)
            }
            
            
          }
          
       #   save(data_limma_fdrall_withfeats_2,file="data_limma_fdrall_withfeats_2.Rda")
          
         
          
          
          allmetabs_res<-data_limma_fdrall_withfeats_2
          
    
          
       
            
          if(analysismode=="classification"){
            
            if(logistic_reg==TRUE){
              
              fname4<-paste("logitreg","results_allfeatures.txt",sep="")
            }else{
              
              
              if(poisson_reg==TRUE){
                fname4<-paste("poissonreg","results_allfeatures.txt",sep="")
              }else{
                fname4<-paste(parentfeatselmethod,"results_allfeatures.txt",sep="")
              }
            }
            
            fname4<-paste("Tables/",fname4,sep="")
            
            
            
            if(is.na(names_with_mz_time)==FALSE){
              
              
              group_means1<-merge(group_means,group_sd,by=c("mz","time"))
              
              allmetabs_res_temp<-merge(group_means1,allmetabs_res,by=c("mz","time"))
              
              
              allmetabs_res_withnames<-merge(names_with_mz_time,allmetabs_res_temp,by=c("mz","time"))
              
            #  allmetabs_res_withnames<-merge(diff_degree_measure[,c("mz","time","DiffRank")],allmetabs_res_withnames,by=c("mz","time"))
              #allmetabs_res_withnames<-cbind(degree_rank,diffexp_rank,allmetabs_res_withnames)
             # allmetabs_res_withnames<-allmetabs_res_withnames[,c("DiffRank")]
              
             # save(allmetabs_res_withnames,file="allmetabs_res_withnames.Rda")
              
             # allmetabs_res_withnames<-allmetabs_res_withnames[order(allmetabs_res_withnames$mz,allmetabs_res_withnames$time),]
              allmetabs_res_withnames<-allmetabs_res_withnames[order(as.numeric(as.character(allmetabs_res_withnames$mz)),as.numeric(as.character(allmetabs_res_withnames$time))),]
              
             
              if(length(check_names)>0){
              rem_col_ind1<-grep(colnames(allmetabs_res_withnames),pattern=c("mz"))
              
              rem_col_ind2<-grep(colnames(allmetabs_res_withnames),pattern=c("time"))
              
              rem_col_ind<-c(rem_col_ind1,rem_col_ind2)
              }else{
                rem_col_ind<-{}
                
              }
              if(length(rem_col_ind)>0){
              write.table(allmetabs_res_withnames[,-c(rem_col_ind)], file=fname4,sep="\t",row.names=FALSE)
              }else{
                
                write.table(allmetabs_res_withnames, file=fname4,sep="\t",row.names=FALSE)
              }
              #rm(data_allinf_withfeats_withnames)
              #}
            }else{
              
              group_means1<-merge(group_means,group_sd,by=c("mz","time"))
              
              allmetabs_res_temp<-merge(group_means1,allmetabs_res,by=c("mz","time"))
              
              #allmetabs_res_temp<-merge(group_means,allmetabs_res,by=c("mz","time"))
              
             # allmetabs_res_temp<-cbind(degree_rank,diffexp_rank,allmetabs_res_temp)
              
              Name<-paste(allmetabs_res_temp$mz,allmetabs_res_temp$time,sep="_")
              allmetabs_res_withnames<-cbind(Name,allmetabs_res_temp)
              allmetabs_res_withnames<-as.data.frame(allmetabs_res_withnames)
             # allmetabs_res_withnames<-allmetabs_res_withnames[order(allmetabs_res_withnames$mz,allmetabs_res_withnames$time),]
              allmetabs_res_withnames<-allmetabs_res_withnames[order(as.numeric(as.character(allmetabs_res_withnames$mz)),as.numeric(as.character(allmetabs_res_withnames$time))),]
              
              write.table(allmetabs_res_withnames,file=fname4,sep="\t",row.names=FALSE)
            }
            
            rm(allmetabs_res_temp)
          }else{
            
            
          }
          
          #rm(allmetabs_res)
          
          if(length(goodip)>=1){
            
           # data_limma_fdrall_withfeats_2<-data_limma_fdrall_withfeats_2[goodip,]
            
            #data_limma_fdrall_withfeats_2<-as.data.frame(data_limma_fdrall_withfeats_2)
           # save(allmetabs_res_withnames,goodip,file="allmetabs_res_withnames.Rda")
            allmetabs_res_withnames<-allmetabs_res_withnames[order(as.numeric(as.character(allmetabs_res_withnames$mz)),as.numeric(as.character(allmetabs_res_withnames$time))),]
            goodfeats<-as.data.frame(allmetabs_res_withnames[goodip,]) #data_limma_fdrall_withfeats_2)
            
           
          #  write.table(allmetabs_res_withnames,file=fname4,sep="\t",row.names=FALSE)
          
            if(logistic_reg==TRUE){
              
              fname4<-paste("logitreg","results_selectedfeatures.txt",sep="")
              
            }else{
              
              if(poisson_reg==TRUE){
                
                fname4<-paste("poissonreg","results_selectedfeatures.txt",sep="")
                
              }else{
                fname4<-paste(featselmethod,"results_selectedfeatures.txt",sep="")
              }
            }
            
            #fname4<-paste("Tables/",fname4,sep="")
            
            write.table(goodfeats,file=fname4,sep="\t",row.names=FALSE)
            
            if(length(rocfeatlist)>length(goodip)){
              
              rocfeatlist<-seq(1,(length(goodip)))
              numselect<-length(goodip)
              rocfeatlist<-rocfeatlist+1
            }else{
              
              numselect<-length(rocfeatlist)
            }
            
            
          }
          
        }else{
          
          #analysismode=="regression"
          if(featselmethod=="lmreg" | featselmethod=="limma" | featselmethod=="limma2way" | featselmethod=="limma1way" | featselmethod=="lmreg" | featselmethod=="logitreg" | featselmethod=="limma1wayrepeat" | featselmethod=="limma2wayrepeat" | featselmethod=="lm1wayanova" | featselmethod=="lm2wayanova" | featselmethod=="lm1wayanovarepeat" | featselmethod=="lm2wayanovarepeat" | featselmethod=="wilcox" | featselmethod=="ttest" | featselmethod=="poissonreg" | featselmethod=="lmregrepeat")
          {
            diffexp_rank<-rank(data_limma_fdrall_withfeats[,1]) #order(data_limma_fdrall_withfeats[,2],decreasing=FALSE)
          }else{
            
            
            if(featselmethod=="rfesvm"){
              
              
              diffexp_rank<-rank_vec
              data_limma_fdrall_withfeats<-data_limma_fdrall_withfeats
              
              
            }else{
              
              
              if(featselmethod=="pamr"){
                
                diffexp_rank<-rank_vec
              #  data_limma_fdrall_withfeats<-cbind(rank_vec,data_limma_fdrall_withfeats)
                
                
              }else{
                
                if(featselmethod=="MARS"){
                  
                  diffexp_rank<-rank((-1)*data_limma_fdrall_withfeats[,1])
                  
                }else{
                  diffexp_rank<-rank((1)*data_limma_fdrall_withfeats[,2])
                }
              }
            }
            
          }
          
        
          #save(goodfeats,diffexp_rank,data_limma_fdrall_withfeats,file="t3.Rda")
          
          data_limma_fdrall_withfeats_2<-cbind(diffexp_rank,data_limma_fdrall_withfeats)
          
          
          # fname4<-paste(featselmethod,"_sigfeats.txt",sep="")
          
          if(logistic_reg==TRUE){
            
            fname4<-paste("logitreg","results_allfeatures.txt",sep="")
            
            
          }else{
            
            
            if(poisson_reg==TRUE){
              fname4<-paste("poissonreg","results_allfeatures.txt",sep="")
            }else{
              fname4<-paste(parentfeatselmethod,"results_allfeatures.txt",sep="")
            }
          }
          
          fname4<-paste("Tables/",fname4,sep="")
          
          allmetabs_res<-data_limma_fdrall_withfeats_2
         
          
          if(is.na(names_with_mz_time)==FALSE){
            
            
            allmetabs_res_withnames<-merge(names_with_mz_time,data_limma_fdrall_withfeats_2,by=c("mz","time"))
            
            # allmetabs_res_withnames<-cbind(degree_rank,diffexp_rank,allmetabs_res_withnames)
            allmetabs_res_withnames<-allmetabs_res_withnames[order(as.numeric(as.character(allmetabs_res_withnames$mz)),as.numeric(as.character(allmetabs_res_withnames$time))),]
            
           # allmetabs_res_withnames<-allmetabs_res_withnames[order(allmetabs_res_withnames$mz,allmetabs_res_withnames$time),]
            
            #write.table(allmetabs_res_withnames[,-c("mz","time")], file=fname4,sep="\t",row.names=FALSE)
           # save(allmetabs_res_withnames,file="allmetabs_res_withnames.Rda")
            #rem_col_ind<-grep(colnames(allmetabs_res_withnames),pattern=c("mz","time"))
            
            if(length(check_names)>0){
            rem_col_ind1<-grep(colnames(allmetabs_res_withnames),pattern=c("mz"))
            
            rem_col_ind2<-grep(colnames(allmetabs_res_withnames),pattern=c("time"))
            
            rem_col_ind<-c(rem_col_ind1,rem_col_ind2)
            }else{
              rem_col_ind<-{}
            }
            
            if(length(rem_col_ind)>0){
              write.table(allmetabs_res_withnames[,-c(rem_col_ind)], file=fname4,sep="\t",row.names=FALSE)
            }else{
              
              write.table(allmetabs_res_withnames, file=fname4,sep="\t",row.names=FALSE)
            }
            
            
           # rm(data_allinf_withfeats_withnames)
            
          }else{
            
            # allmetabs_res_temp<-cbind(degree_rank,diffexp_rank,allmetabs_res)
            
            allmetabs_res_withnames<-allmetabs_res
            write.table(allmetabs_res,file=fname4,sep="\t",row.names=FALSE)
          }
          goodfeats<-allmetabs_res_withnames[goodip,] #data_limma_fdrall_withfeats_2[goodip,] #[sel.diffdrthresh==TRUE,]
         # save(allmetabs_res_withnames,goodip,file="allmetabs_res_withnames.Rda")
          goodfeats<-as.data.frame(allmetabs_res_withnames[goodip,]) #data_limma_fdrall_withfeats_2)
          
          if(logistic_reg==TRUE){
            
            fname4<-paste("logitreg","results_selectedfeatures.txt",sep="")
            
          }else{
            
            if(poisson_reg==TRUE){
              
              fname4<-paste("poissonreg","results_selectedfeatures.txt",sep="")
              
            }else{
              fname4<-paste(featselmethod,"results_selectedfeatures.txt",sep="")
            }
          }
          
         # fname4<-paste("Tables/",fname4,sep="")
          
          write.table(goodfeats,file=fname4,sep="\t",row.names=FALSE)
          
          
          
          fname4<-paste("Tables/",parentfeatselmethod,"results_allfeatures.txt",sep="")
          
          #allmetabs_res<-goodfeats #data_limma_fdrall_withfeats_2
         
        }
        
        
      }
      
      
  #   save(goodfeats,file="goodfeats455.Rda")
      
      if(length(goodip)>1){
        goodfeats_by_DICErank<-{}
        
        if(analysismode=="classification"){
          if(featselmethod=="lmreg" | featselmethod=="limma" | featselmethod=="limma2way" | featselmethod=="limma1way" | featselmethod=="lmreg" | featselmethod=="logitreg" | featselmethod=="limma1wayrepeat" | featselmethod=="limma2wayrepeat" | featselmethod=="lm1wayanova" | featselmethod=="lm2wayanova" | featselmethod=="lm1wayanovarepeat" | featselmethod=="lm2wayanovarepeat" | featselmethod=="wilcox" | featselmethod=="ttest" | featselmethod=="poissonreg")
          {
            goodfeats<-goodfeats[order(goodfeats$diffexp_rank,decreasing=FALSE),]
            
            if(length(goodip)>1){
              # goodfeats_by_DICErank<-data_limma_fdrall_withfeats_2[r1$top.list,]
            }
          }else{
            
            goodfeats<-goodfeats[order(goodfeats$diffexp_rank,decreasing=FALSE),]
            
            if(length(goodip)>1){
              #goodfeats_by_DICErank<-data_limma_fdrall_withfeats_2[r1$top.list,]
              
            }
          }
          
          cnamesd1<-colnames(goodfeats)
          time_ind<-which(cnamesd1=="time")
          
          
          mz_ind<-which(cnamesd1=="mz")
          
          
          goodfeats_name<-goodfeats$Name
          
          
          goodfeats_temp<-cbind(goodfeats[,mz_ind],goodfeats[,time_ind],goodfeats[,which(colnames(goodfeats)%in%sample_names_vec)]) #goodfeats[,-c(1:time_ind)])
          
         # save(goodfeats_temp,file="goodfeats_temp.Rda")
          
          cnames_temp<-colnames(goodfeats_temp)
          cnames_temp<-c("mz","time",cnames_temp[-c(1:2)])
          colnames(goodfeats_temp)<-cnames_temp
          
          goodfeats<-goodfeats_temp
          
          
        }else{
          if(analysismode=="regression"){
            
           # save(goodfeats,file="goodfeats455.Rda")
            
            try(dev.off(),silent=TRUE)
            
            if(featselmethod=="lmreg" | featselmethod=="pls" | featselmethod=="spls" |      featselmethod=="o1pls" | featselmethod=="RF" | featselmethod=="MARS"){
              
              ####savegoodfeats,file="goodfeats.Rda")
              goodfeats<-goodfeats[order(goodfeats$diffexp_rank,decreasing=FALSE),]
              
            }else{
              
              
              #goodfeats<-goodfeats[order(goodfeats[,1],decreasing=TRUE),]
              
            }
            
            goodfeats<-as.data.frame(goodfeats)
            
            cnamesd1<-colnames(goodfeats)
            time_ind<-which(cnamesd1=="time")
            
            
            mz_ind<-which(cnamesd1=="mz")
            
            
            goodfeats_name<-goodfeats$Name
            
            
            goodfeats_temp<-cbind(goodfeats[,mz_ind],goodfeats[,time_ind],goodfeats[,which(colnames(goodfeats)%in%sample_names_vec)]) #goodfeats[,-c(1:time_ind)])
            
            #save(goodfeats_temp,goodfeats,goodfeats_name,file="goodfeats_temp.Rda")
            cnames_temp<-colnames(goodfeats_temp)
            cnames_temp<-c("mz","time",cnames_temp[-c(1:2)])
            colnames(goodfeats_temp)<-cnames_temp
            
            goodfeats<-goodfeats_temp
            
            rm(goodfeats_temp)
            print("PCA selected features")
            #      #save(goodfeats,goodfeats_temp,mz_ind,time_ind,classlabels_orig,analysistype,alphabetical.order,col_vec,file="pca1.Rda")
            
            num_sig_feats<-nrow(goodfeats)
            if(num_sig_feats>=3){
              if(output.device.type!="pdf"){
                
                temp_filename_1<-"Figures/PCAplots_selectedfeats.pdf"
                
                #png(temp_filename_1,width=plots.width,height=plots.height,res=plots.res,type=plots.type,units="in")
                #pdf(temp_filename_1)
                pdf(temp_filename_1,width=plots.width,height=plots.height)
              }
              
              plot(0:10, type = "n", xaxt="n", yaxt="n", bty="n", xlab = "", ylab = "")
              
              
              text(5, 8, "PCA using selected features after feature selection")
              text(5, 7, "The figures include: ")
              text(5, 6, "a. pairwise PC score plots ")
              text(5, 5, "b. scores for individual samples on each PC")
              text(5, 4, "c. Lineplots using PC scores for data with repeated measurements")
              
              
              par(mfrow=c(1,1),family="sans",cex=cex.plots)
              
              rownames(goodfeats)<-goodfeats$Name
              
              get_pcascoredistplots(X=goodfeats,Y=classlabels_orig,feature_table_file=NA,parentoutput_dir=getwd(),
                                    class_labels_file=NA,sample.col.opt=sample.col.opt,plots.width=2000,plots.height=2000,
                                    plots.res=300, alphacol=0.3,col_vec=col_vec,pairedanalysis=pairedanalysis,
                                    pca.cex.val=pca.cex.val,legendlocation=legendlocation,pca.ellipse=pca.ellipse,
                                    ellipse.conf.level=ellipse.conf.level,filename="selected",paireddesign=paireddesign,
                                    lineplot.col.opt=lineplot.col.opt,lineplot.lty.option=lineplot.lty.option,
                                    timeseries.lineplots=timeseries.lineplots,pcacenter=pcacenter,pcascale=pcascale,
                                    alphabetical.order=alphabetical.order,study.design=analysistype,lme.modeltype=modeltype) #,silent=TRUE)
              
              if(output.device.type!="pdf"){
                
                try(dev.off(),silent=TRUE)
              }
            }
          }
          
        }
        
       
        
        
      }
      
      
      class_label_A<-class_labels_levels[1]
      class_label_B<-class_labels_levels[2]
      
      goodfeats_allfields<-{}
      
      
      
      if(length(which(sel.diffdrthresh==TRUE))>1){
        
        goodfeats<-as.data.frame(goodfeats)
        mzvec<-goodfeats$mz
        timevec<-goodfeats$time
        
        if(length(mzvec)>4){
          max_per_row<-3
          
          
          par_rows<-ceiling(9/max_per_row)
          
        }else{
          max_per_row<-length(mzvec)
          par_rows<-1
        }
        
        
        goodfeats<-as.data.frame(goodfeats)
        
        cnamesd1<-colnames(goodfeats)
        time_ind<-which(cnamesd1=="time")
        
        goodfeats_allfields<-as.data.frame(goodfeats)
        
        
        file_ind<-1
        
        
        mz_ind<-which(cnamesd1=="mz")
        
       
        goodfeats_temp<-cbind(goodfeats[,mz_ind],goodfeats[,time_ind],goodfeats[,-c(1:time_ind)])
       
        cnames_temp<-colnames(goodfeats_temp)
        cnames_temp[1]<-"mz"
        cnames_temp[2]<-"time"
        colnames(goodfeats_temp)<-cnames_temp
        
        
        
        
        
        #if(length(class_labels_levels)<10)
        if(analysismode=="classification" && nrow(goodfeats)>=1 && length(goodip)>=1)
        {
          
          
          
          if(length(class_labels_levels)==2){
            
            #print("Generating ROC curve using top features on training set")
            
            
            #try(get_roc(dataA=goodfeats_temp,classlabels=classlabels,classifier=rocclassifier,kname="radial",rocfeatlist=rocfeatlist,rocfeatincrement=rocfeatincrement,mainlabel="Training set ROC curve using top features"),silent=TRUE)
            
          }
          
          
          subdata=t(goodfeats[,-c(1:time_ind)])
          
        # save(kfold,subdata,goodfeats,classlabels,svm_kernel,pred.eval.method,match_class_dist,file="svmdebug.Rda")
          
          svm_model<-try(svm_cv(v=kfold,x=subdata,y=classlabels,kname=svm_kernel,errortype=pred.eval.method,conflevel=95,match_class_dist=match_class_dist),silent=TRUE)
          
          #svm_model<-try(svm_cv(v=kfold,x=subdata,y=classlabels,kname=svm_kernel,errortype=pred.eval.method,conflevel=95,match_class_dist=match_class_dist),silent=TRUE)
          #svm_model<-try(svm(x=subdata,y=(classlabels),type="nu-classification",cross=kfold,kernel=svm_kernel),silent=TRUE)
          
          #svm_model<-try(svm_cv(v=kfold,x=subdata,y=classlabels,kname=svm_kernel,errortype=pred.eval.method,conflevel=95,match_class_dist=match_class_dist),silent=TRUE)
          
          classlabels<-as.data.frame(classlabels)
          
          if(is(svm_model,"try-error")){
            print("SVM could not be performed. Please try lowering the kfold or set kfold=total number of samples for Leave-one-out CV. Skipping to the next step.")
            termA<-(-1)
            pred_acc<-termA
            permut_acc<-(-1)
          }else{
            
            pred_acc<-svm_model$avg_acc
            
            print("Accuracy is:")
            print(pred_acc)
            
            print("Calculating permuted CV accuracy")
            
            permut_acc<-{}
            #permut_acc<-lapply(1:100,function(j){
            numcores<-num_nodes #round(detectCores()*0.5)
            cl <- parallel::makeCluster(getOption("cl.cores", num_nodes))
            clusterEvalQ(cl,library(e1071))
            clusterEvalQ(cl,library(pROC))
            clusterEvalQ(cl,library(ROCR))
            clusterEvalQ(cl,library(CMA))
            clusterExport(cl,"svm_cv",envir = .GlobalEnv)
            permut_acc<-parLapply(cl,1:cv.perm.count,function(p1){
              
              rand_order<-sample(1:dim(classlabels)[1],size=dim(classlabels)[1])
              classlabels_permut<-classlabels[rand_order,]
              classlabels_permut<-as.data.frame(classlabels_permut)
              svm_permut_res<-try(svm_cv(v=kfold,x=subdata,y=classlabels_permut,kname=svm_kernel,errortype=pred.eval.method,conflevel=95,match_class_dist=match_class_dist),silent=TRUE)
              
              #svm_permut_res<-try(svm(x=subdata,y=(classlabels_permut),type="nu-classification",cross=kfold,kernel=svm_kernel),silent=TRUE)
              #svm_permut_res<-svm_cv(v=kfold,x=subdata,y=classlabels_permut,kname=svm_kernel,errortype=pred.eval.method,conflevel=95,match_class_dist=match_class_dist)
              
              
              if(is(svm_permut_res,"try-error")){
                
                cur_perm_acc<-NA
              }else{
                cur_perm_acc<-svm_permut_res$avg_acc #tot.accuracy #
              }
              return(cur_perm_acc)
            })
            
            stopCluster(cl)
            
            permut_acc<-unlist(permut_acc)
            permut_acc<-mean(permut_acc,na.rm=TRUE)
            permut_acc<-round(permut_acc,2)
            
            print("mean Permuted accuracy is:")
            print(permut_acc)
            
          }
          
          
          
          
          termA<-100*pred_acc
          
          if(featselmethod=="limma" | featselmethod=="limma2way" | featselmethod=="limma2wayrepeat" | featselmethod=="lmreg" | featselmethod=="logitreg"
             | featselmethod=="lm2wayanova" | featselmethod=="lm1wayanova" | featselmethod=="lm1wayanovarepeat" | featselmethod=="lm2wayanovarepeat" | featselmethod=="wilcox" | featselmethod=="ttest" |  featselmethod=="poissonreg" | featselmethod=="lmregrepeat")
          {
            if(fdrmethod=="none"){
              exp_fp<-(dim(data_m_fc)[1]*fdrthresh)+1
            }else{
              exp_fp<-(feat_sigfdrthresh[lf]*fdrthresh)+1
            }
          }
          
          
          
          termB<-(dim(parent_data_m)[1]*dim(parent_data_m)[1])/(dim(data_m_fc)[1]*dim(data_m_fc)[1]*100)
          
          
          res_score<-(100*(termA-permut_acc))-(feat_weight*termB*exp_fp)
          res_score<-round(res_score,2)
          
          
          
          if(lf==0)
          {
            best_logfc_ind<-lf
            
            best_feats<-goodip
            best_cv_res<-res_score
            best_acc<-pred_acc
            best_limma_res<-data_limma_fdrall_withfeats[goodip,] #[sel.diffdrthresh==TRUE,]
            
            
            
          }else{
            
            if(res_score>best_cv_res){
              
              best_logfc_ind<-lf
              
              best_feats<-goodip
              best_cv_res<-res_score
              best_acc<-pred_acc
              best_limma_res<-data_limma_fdrall_withfeats[goodip,] #[sel.diffdrthresh==TRUE,]
            }
            
          }
          
          res_score_vec[lf]<-res_score
          if(pred.eval.method=="CV"){
            feat_sigfdrthresh_cv[lf]<-pred_acc
            feat_sigfdrthresh_permut[lf]<-permut_acc
            print(paste(kfold,"-fold CV accuracy: ", pred_acc,sep=""))
            print(paste("Permuted ",kfold,"-fold CV accuracy: ", permut_acc,sep=""))
            
          }else{
            if(pred.eval.method=="AUC"){
              feat_sigfdrthresh_cv[lf]<-pred_acc
              feat_sigfdrthresh_permut[lf]<-permut_acc
              print(paste("ROC area under the curve (AUC) is : ", pred_acc,sep=""))
              print(paste("Permuted ROC area under the curve (AUC) is : ", permut_acc,sep=""))
            }else{
              if(pred.eval.method=="BER"){
                feat_sigfdrthresh_cv[lf]<-pred_acc
                feat_sigfdrthresh_permut[lf]<-permut_acc
                print(paste("Balanced accuracy rate is : ", pred_acc,sep=""))
                print(paste("Permuted balanced accuracy rate is : ", permut_acc,sep=""))
              }
            }
          }
          
          #print("ROC done")
          best_subset<-{}
          best_acc<-0
          
          xvec<-{}
          yvec<-{}
          #for(i in 2:max_varsel)
          
          if(nrow(goodfeats_temp)>30){
            
            max_cv_varsel<-30
          }else{
            max_cv_varsel<-nrow(goodfeats_temp)
          }
          
          cv_yvec<-lapply(2:max_cv_varsel,function(i)
          {
            
            subdata<-t(goodfeats_temp[1:i,-c(1:2)])
            svm_model<-try(svm_cv(v=kfold,x=subdata,y=classlabels,kname=svm_kernel,errortype=pred.eval.method,conflevel=95,match_class_dist=match_class_dist),silent=TRUE)
            
            #svm_model<-svm_cv(v=kfold,x=subdata,y=classlabels,kname=svm_kernel,errortype=pred.eval.method,conflevel=95,match_class_dist=match_class_dist)
            
            
            if(is(svm_model,"try-error")){
              
              res1<-NA
              
            }else{
              
              
              
              res1<-svm_model$avg_acc
              
              
            }
            
            return(res1)
          })
          
          xvec<-seq(2:max_cv_varsel)
          
          yvec<-unlist(cv_yvec)
          
         
          
          
          if(pred.eval.method=="CV"){
            ylab_text=paste(pred.eval.method," accuracy (%)",sep="")
            
          }else{
            if(pred.eval.method=="BER"){
              ylab_text=paste("Balanced accuracy"," (%)",sep="")
            }else{
              
              ylab_text=paste("AUC"," (%)",sep="")
            }
          }
          
          
          if(length(yvec)>0){
            
            
            if(output.device.type!="pdf"){
              
              temp_filename_1<-"Figures/kfoldCV_forward_selection.png"
              
              png(temp_filename_1,width=plots.width,height=plots.height,res=plots.res,type=plots.type,units="in")
            }else{
              temp_filename_1<-"Figures/kfoldCV_forward_selection.pdf"
              
              pdf(temp_filename_1)
              
            }
            
            
            
            try(plot(x=xvec,y=yvec,main="k-fold CV classification accuracy based on forward selection of top features",xlab="Feature index",ylab=ylab_text,type="b",col="#0072B2",cex.main=0.7),silent=TRUE)
            
            
            if(output.device.type!="pdf"){
              
              try(dev.off(),silent=TRUE)
            }else{
              try(dev.off(),silent=TRUE)
              
            }
            
            cv_mat<-cbind(xvec,yvec)
            colnames(cv_mat)<-c("Feature Index",ylab_text)
            
            write.table(cv_mat,file="Tables/kfold_cv_mat.txt",sep="\t")
          }
          
          if(pairedanalysis==TRUE)
          {
            
            if(featselmethod=="pls" | featselmethod=="spls"){
              classlabels_sub<-classlabels_sub[,-c(1)]
              classlabels_temp<-cbind(classlabels_sub)
            }else{
              classlabels_sub<-classlabels_sub[,-c(1)]
              classlabels_temp<-cbind(classlabels_sub)
              
            }
            
          }else{
            classlabels_temp<-cbind(classlabels_sub,classlabels)
          }
          
          
          num_sig_feats<-nrow(goodfeats)
          
          if(num_sig_feats<3){
            
            pca.stage2.eval=FALSE
            
          }
          
          if(pca.stage2.eval==TRUE)
          {
            
            
            pca_comp<-min(10,dim(X)[2])
            
            #dev.off()
            
            # print("plotting")
            #pdf("sig_features_evaluation.pdf", height=2000,width=2000)
            library(pcaMethods)
            
            p1<-pcaMethods::pca(X,method="rnipals",center=TRUE,scale="uv",cv="q2",nPcs=pca_comp)
            
            if(output.device.type!="pdf"){
              
              temp_filename_1<-"Figures/PCAdiagnostics_selectedfeats.png"
              
              png(temp_filename_1,width=plots.width,height=plots.height,res=plots.res,type=plots.type,units="in")
            }
            
            
            
            p2<-plot(p1,col=c("darkgrey","grey"),main="PCA diagnostics after variable selection")
            
            print(p2)
            if(output.device.type!="pdf"){
              
              try(dev.off(),silent=TRUE)
            }
            #dev.off()
            
          }
          
          classlabels_orig<-classlabels_orig_parent
          
          
          if(pairedanalysis==TRUE){
            
            classlabels_orig<-classlabels_orig[,-c(2)]
            
            
          }else{
            
            if(featselmethod=="lmreg" || featselmethod=="logitreg" || featselmethod=="poissonreg"){
              classlabels_orig<-classlabels_orig[,c(1:2)]
              classlabels_orig<-as.data.frame(classlabels_orig)
            }
          }
          
          
          
          classlabels_orig_wgcna<-classlabels_orig
          
          
          goodfeats_temp<-cbind(goodfeats[,mz_ind],goodfeats[,time_ind],goodfeats[,-c(1:time_ind)])
          cnames_temp<-colnames(goodfeats_temp)
          cnames_temp<-c("mz","time",cnames_temp[-c(1:2)])
          colnames(goodfeats_temp)<-cnames_temp
          
          goodfeats_temp_with_names<-merge(names_with_mz_time,goodfeats_temp,by=c("mz","time"))
          
          
          goodfeats_temp_with_names<-goodfeats_temp_with_names[match(paste(goodfeats_temp$mz,"_",goodfeats_temp$time,sep=""),paste(goodfeats_temp_with_names$mz,"_",goodfeats_temp_with_names$time,sep="")),]
          
            
         # save(goodfeats,goodfeats_temp,names_with_mz_time,goodfeats_temp_with_names,file="goodfeats_pca.Rda")
          
          rownames(goodfeats_temp)<-goodfeats_temp_with_names$Name
          
          
          if(num_sig_feats>=3){
            if(output.device.type!="pdf"){
              
              temp_filename_1<-"Figures/PCAplots_selectedfeats.pdf"
              
              #png(temp_filename_1,width=plots.width,height=plots.height,res=plots.res,type=plots.type,units="in")
              #pdf(temp_filename_1)
              pdf(temp_filename_1,width=plots.width,height=plots.height)
            }
            
            plot(0:10, type = "n", xaxt="n", yaxt="n", bty="n", xlab = "", ylab = "")
            
            
            text(5, 8, "PCA using selected features after feature selection")
            text(5, 7, "The figures include: ")
            text(5, 6, "a. pairwise PC score plots ")
            text(5, 5, "b. scores for individual samples on each PC")
            text(5, 4, "c. Lineplots using PC scores for data with repeated measurements")
            
            
            par(mfrow=c(1,1),family="sans",cex=cex.plots)
            
            
            
            get_pcascoredistplots(X=goodfeats_temp,Y=classlabels_orig_pca,
                                  feature_table_file=NA,parentoutput_dir=getwd(),class_labels_file=NA,
                                  sample.col.opt=sample.col.opt,plots.width=2000,plots.height=2000,plots.res=300, alphacol=0.3,col_vec=col_vec,pairedanalysis=pairedanalysis,pca.cex.val=pca.cex.val,legendlocation=legendlocation,pca.ellipse=pca.ellipse,ellipse.conf.level=ellipse.conf.level,filename="selected",paireddesign=paireddesign,lineplot.col.opt=lineplot.col.opt,lineplot.lty.option=lineplot.lty.option,timeseries.lineplots=timeseries.lineplots,pcacenter=pcacenter,pcascale=pcascale,alphabetical.order=alphabetical.order,study.design=analysistype,lme.modeltype=modeltype) #,silent=TRUE)
            
            if(output.device.type!="pdf"){
              
              try(dev.off(),silent=TRUE)
            }
          }
          
          
          ####savelist=ls(),file="timeseries.Rda")
          
          #if(FALSE)
          {
            #if(FALSE)
            {
              if(log2transform==TRUE || input.intensity.scale=="log2"){
                
                if(znormtransform==TRUE){
                  ylab_text_2="scale normalized"
                }else{
                  if(quantile_norm==TRUE){
                    
                    ylab_text_2="quantile normalized"
                  }else{
                    if(eigenms_norm==TRUE){
                      
                      ylab_text_2="EigenMS normalized"
                    }else{
                      
                      if(sva_norm==TRUE){
                        
                        ylab_text_2="SVA normalized"
                      }else{
                        ylab_text_2=""
                      }
                    }
                  }
                }
                ylab_text=paste("log2 intensity ",ylab_text_2,sep="")
              }else{
                if(znormtransform==TRUE){
                  ylab_text_2="scale normalized"
                }else{
                  if(quantile_norm==TRUE){
                    
                    ylab_text_2="quantile normalized"
                  }else{
                    #ylab_text_2=""
                    if(medcenter==TRUE){
                      
                      ylab_text_2="median centered"
                    }else{
                      
                      if(lowess_norm==TRUE){
                        ylab_text_2="LOWESS normalized"
                        
                      }else{
                        if(rangescaling==TRUE){
                          ylab_text_2="range scaling normalized"
                          
                        }else{
                          
                          if(paretoscaling==TRUE){
                            ylab_text_2="pareto scaling normalized"
                            
                          }else{
                            
                            if(mstus==TRUE){
                              ylab_text_2="MSTUS normalized"
                              
                            }else{
                              
                              if(vsn_norm==TRUE){
                                ylab_text_2="VSN normalized"
                                
                              }else{
                                ylab_text_2=""
                                
                              }
                              
                            }
                          }
                          
                        }
                        
                      }
                      
                    }
                    
                  }
                }
                ylab_text=paste("Intensity ",ylab_text_2,sep="")
              }
              
            }
            
            #ylab_text_2=""
            #ylab_text=paste("Abundance",ylab_text_2,sep="")
            
            par(mfrow=c(1,1),family="sans",cex=cex.plots)
            
            
            
            if(pairedanalysis==TRUE || timeseries.lineplots==TRUE)
            {
              
              if(output.device.type!="pdf"){
                
                temp_filename_1<-"Figures/Lineplots_selectedfeats.pdf"
                
                #png(temp_filename_1,width=plots.width,height=plots.height,res=plots.res,type=plots.type,units="in")
                #pdf(temp_filename_1)
                pdf(temp_filename_1,width=plots.width,height=plots.height)
               # par(mfrow=c(1,1))
                
                par(mfrow=c(1,1),family="sans",cex=cex.plots)
              }
              
              
              #plot(0:10, type = "n", xaxt="n", yaxt="n", bty="n", xlab = "", ylab = "")
              
              
              #text(5, 8, "Lineplots using selected features")
            #  text(5, 7, "The error bars represent the 95% \nconfidence interval in each group (or timepoint)")
              
              
         #    save(goodfeats_temp,classlabels_orig,lineplot.col.opt,col_vec,pairedanalysis,
          #         pca.cex.val,pca.ellipse,ellipse.conf.level,legendlocation,ylab_text,error.bar,
           #       cex.plots,lineplot.lty.option,timeseries.lineplots,analysistype,goodfeats_name,alphabetical.order,
            #      multiple.figures.perpanel,plot.ylab_text,plots.height,plots.width,file="debuga_lineplots.Rda")
              
              
              #try(
                var_sum_list<-get_lineplots(X=goodfeats_temp,Y=classlabels_orig,feature_table_file=NA,
                                parentoutput_dir=getwd(),class_labels_file=NA,
                                lineplot.col.opt=lineplot.col.opt,alphacol=alphacol,col_vec=col_vec,
                                pairedanalysis=pairedanalysis,point.cex.val=pca.cex.val,
                                legendlocation=legendlocation,pca.ellipse=pca.ellipse,
                                ellipse.conf.level=ellipse.conf.level,filename="selected",
                                ylabel=plot.ylab_text,error.bar=error.bar,cex.plots=cex.plots,
                                lineplot.lty.option=lineplot.lty.option,timeseries.lineplots=timeseries.lineplots,
                                name=goodfeats_name,study.design=analysistype,
                              alphabetical.order=alphabetical.order,multiple.figures.perpanel=multiple.figures.perpanel,
                              plot.height = plots.height,plot.width=plots.width)
                #,silent=TRUE)  #,silent=TRUE)
              #save(var_sum_list,file="var_sum_list.Rda")
              var_sum_mat<-{}
             # for(i in 1:length(var_sum_list))
              #{
               #   var_sum_mat<-rbind(var_sum_mat,var_sum_list[[i]]$df_write_temp)
                
              #}
               # var_sum_mat<-ldply(var_sum_list,rbind)
                                   
               # write.table(var_sum_mat,file="Tables/data_summary.txt",sep="\t",row.names=FALSE)
                
              if(output.device.type!="pdf"){
                
                try(dev.off(),silent=TRUE)
              }
              
            }
            
            
            
            
          }
          
          
     #     save(goodfeats_temp,classlabels_orig,lineplot.col.opt,alphacol,col_vec,pairedanalysis,pca.cex.val,legendlocation,pca.ellipse,ellipse.conf.level,plot.ylab_text,error.bar,cex.plots,
     #          lineplot.lty.option,timeseries.lineplots,goodfeats_name,analysistype,alphabetical.order,multiple.figures.perpanel,plots.height,plots.width,file="var_sum.Rda")
        
            var_sum_list<-get_data_summary(X=goodfeats_temp,Y=classlabels_orig,feature_table_file=NA,
                                        parentoutput_dir=getwd(),class_labels_file=NA,
                                        lineplot.col.opt=lineplot.col.opt,alphacol=alphacol,col_vec=col_vec,
                                        pairedanalysis=pairedanalysis,point.cex.val=pca.cex.val,
                                        legendlocation=legendlocation,pca.ellipse=pca.ellipse,
                                        ellipse.conf.level=ellipse.conf.level,filename="selected",
                                        ylabel=plot.ylab_text,error.bar=error.bar,cex.plots=cex.plots,
                                        lineplot.lty.option=lineplot.lty.option,timeseries.lineplots=timeseries.lineplots,
                                        name=goodfeats_name,study.design=analysistype,
                                        alphabetical.order=alphabetical.order,multiple.figures.perpanel=multiple.figures.perpanel,plot.height = plots.height,plot.width=plots.width)
          
          
         # save(var_sum_list,file="var_sum_list2.Rda")
          
          #save(var_sum_list,file="var_sum_list.Rda")
         # var_sum_mat<-{}
          #for(i in 1:length(var_sum_list))
          #{
           # var_sum_mat<-rbind(var_sum_mat,var_sum_list[[i]]$df_write_temp)
            
          #}
          # var_sum_mat<-ldply(var_sum_list,rbind)
          
         # write.table(var_sum_mat,file="Tables/data_summary2.txt",sep="\t",row.names=FALSE)
          
          if(nrow(goodfeats)<1){
            
            print(paste("No features selected for ",featselmethod,sep=""))
          }
          #else
          {
            
            
            #write.table(goodfeats_temp,file="Tables/boxplots_file.normalized.txt",sep="\t",row.names=FALSE)
            
            goodfeats<-goodfeats[,-c(1:time_ind)]
            
            
            goodfeats_raw<-data_matrix_beforescaling_rsd[goodip,]
            #write.table(goodfeats_raw,file="Tables/boxplots_file.raw.txt",sep="\t",row.names=FALSE)
            
            goodfeats_raw<-goodfeats_raw[match(paste(goodfeats_temp$mz,"_",goodfeats_temp$time,sep=""),paste(goodfeats_raw$mz,"_",goodfeats_raw$time,sep="")),] 
            
            goodfeats_name<-as.character(goodfeats_name)
       #    save(goodfeats_name,goodfeats_temp,classlabels_orig,output_dir,boxplot.col.opt,cex.plots,ylab_text,file="boxplotdebug.Rda")
            
          #  if(FALSE)
            {
              if(output.device.type!="pdf"){
                
                try(dev.off(),silent=TRUE)
                
                temp_filename_1<-"Figures/Pairwise.correlation.plots.pdf"
                
               # pdf(temp_filename_1)
                pdf(temp_filename_1,width=plots.width,height=plots.height)
              }
              
              par(mfrow=c(1,1),family="sans",cex=cex.plots,cex.main=0.7)
              
             # cor1<-WGCNA::cor(t(goodfeats_temp[,-c(1:2)]))
              rownames(goodfeats_temp)<-goodfeats_name
              
              #Pairwise correlations between selected features
              cor1<-WGCNA::cor(t(goodfeats_temp[,-c(1:2)]),nThreads=num_nodes,method=cor.method,use = 'p')
              
              corpval1=apply(cor1,2,function(x){corPvalueStudent(x,n=ncol(goodfeats_temp[,-c(1:2)]))})
              
             # save(cor1,goodfeats_temp,corpval1,goodfeats_name,file="cor1.Rda")
              
              fdr_adjust_pvalue<-try(suppressWarnings(fdrtool(as.vector(cor1[upper.tri(cor1)]),statistic="correlation",verbose=FALSE,plot=FALSE)),silent=TRUE)
              
              cor1[(abs(cor1)<abs.cor.thresh)]<-0
              newnet <- cor1
              newnet[upper.tri(newnet)][fdr_adjust_pvalue$qval > cor.fdrthresh] <- 0
            #  newnet[upper.tri(newnet)][as.vector(corpval1[upper.tri(corpval1)]) > pvalue.thresh] <- 0
              
              newnet[lower.tri(newnet)] <- t(newnet)[lower.tri(newnet)]
              newnet <- as.matrix(newnet)
              
              corqval1=newnet
              diag(corqval1)<-0
              upperTriangle<-upper.tri(cor1, diag=F)
              lowerTriangle<-lower.tri(cor1, diag=F)
              corqval1[upperTriangle]<-fdr_adjust_pvalue$qval
              corqval1[lowerTriangle]<-corqval1[upperTriangle]
              
              cor1=newnet
              rm(newnet)
              
           #   rownames(cor1)<-paste(goodfeats_temp[,c(1)],goodfeats_temp[,c(2)],sep="_")
            #  colnames(cor1)<-rownames(cor1)
              
              #dendrogram="none",
              h1<-heatmap.2(cor1,col=rev(brewer.pal(11,"RdBu")),Rowv=TRUE,Colv=TRUE,scale="none",key=TRUE, symkey=FALSE, density.info="none", trace="none",main="Pairwise correlations between selected features",cexRow = 0.5,cexCol = 0.5,cex.main=0.7)
              
              upperTriangle<-upper.tri(cor1, diag=F) #turn into a upper triangle
              cor1.upperTriangle<-cor1 #take a copy of the original cor-mat
              cor1.upperTriangle[!upperTriangle]<-NA#set everything not in upper triangle o NA
              correlations_melted<-na.omit(melt(cor1.upperTriangle, value.name ="correlationCoef")) #use melt to reshape the matrix into triplets, na.omit to get rid of the NA rows
              colnames(correlations_melted)<-c("from", "to", "weight")
              
              correlations_melted<-as.data.frame(correlations_melted)
              
              correlations_melted$from<-paste("X",correlations_melted$from,sep="")
              correlations_melted$to<-paste("Y",correlations_melted$to,sep="")
              
              write.table(correlations_melted,file="Tables/pairwise.correlations.selectedfeatures.linkmatrix.txt",sep="\t",row.names=FALSE)
              
              if(ncol(cor1)>1000){
              netres<-plot_graph(correlations_melted,filename="sigfeats_top1000pairwisecor",interactive=FALSE,maxnodesperclass=1000,label.cex=network.label.cex,mtext.val="Top 1000 pairwise correlations between selected features")
              }
              netres<-plot_graph(correlations_melted,filename="sigfeats_pairwisecorrelations",interactive=FALSE,maxnodesperclass=NA,label.cex=network.label.cex,mtext.val="Pairwise correlations between selected features")
              
              if(output.device.type!="pdf"){
                
                try(dev.off(),silent=TRUE)
              }
              
            }
            
            
            
            
            
            
            
            
            if(output.device.type!="pdf"){
              
              try(dev.off(),silent=TRUE)
              
              temp_filename_1<-"Figures/Boxplots.selectedfeats.normalized.pdf"
              
              if(boxplot.type=="simple"){
                pdf(temp_filename_1,height=plots.height,width=plots.width)
              }
            }
            
            goodfeats_name<-as.character(goodfeats_name)
        #   save(goodfeats_name,goodfeats_temp,classlabels_orig,output_dir,boxplot.col.opt,cex.plots,ylab_text,plot.ylab_text,
         #        analysistype,boxplot.type,alphabetical.order,goodfeats_name,add.pvalues,add.jitter,file="boxplotdebug.Rda")
            
            par(mfrow=c(1,1),family="sans",cex=cex.plots)
            
            print("Generating boxplots")
           # plot(0:10, type = "n", xaxt="n", yaxt="n", bty="n", xlab = "", ylab = "")
            
            
          #  text(5, 8, "Boxplots of selected features using the\n normalized intensities/abundance levels",cex=1.5,font=2)
            
            #plot.ylab_text1=paste("(Normalized) ",ylab_text,sep="")
            
            #classlabels_paired<-cbind(as.character(classlabels[,1]),as.character(subject_inf),as.character(classlabels[,2]))
           
            #classlabels_paired<-as.data.frame(classlabels_paired)
            
            
            if(normalization.method!="none"){
              
              if(pairedanalysis==TRUE){
                
               #classlabels_paired<-cbind(classlabels[,1],subject_inf,classlabels[,2])
                
                res<-get_boxplots(X=goodfeats_temp,Y=classlabels_orig,parentoutput_dir=output_dir,boxplot.col.opt=boxplot.col.opt,
                                  newdevice=FALSE,cex.plots=cex.plots,ylabel=plot.ylab_text,name=goodfeats_name,add.pvalues=add.pvalues,add.jitter=add.jitter,
                                  alphabetical.order=alphabetical.order,boxplot.type=boxplot.type,study.design=gsub(analysistype,pattern="repeat",replacement=""),
                                  multiple.figures.perpanel=multiple.figures.perpanel,numnodes=num_nodes,
                                  plot.height = plots.height,plot.width=plots.width,
                                  filename="Figures/Boxplots.selectedfeats.normalized",alphacol = alpha.col,ggplot.type1=ggplot.type1,facet.nrow=facet.nrow)
                
              }else{
              
              res<-get_boxplots(X=goodfeats_temp,Y=classlabels_orig,parentoutput_dir=output_dir,boxplot.col.opt=boxplot.col.opt,
                          newdevice=FALSE,cex.plots=cex.plots,ylabel=plot.ylab_text,name=goodfeats_name,add.pvalues=add.pvalues,add.jitter=add.jitter,
                           alphabetical.order=alphabetical.order,boxplot.type=boxplot.type,study.design=analysistype,
                           multiple.figures.perpanel=multiple.figures.perpanel,numnodes=num_nodes,
                          plot.height = plots.height,plot.width=plots.width,
                           filename="Figures/Boxplots.selectedfeats.normalized",alphacol = alpha.col,ggplot.type1=ggplot.type1,facet.nrow=facet.nrow)
              }
              
            }else{
              
              plot.boxplots.raw=TRUE
              goodfeats_raw=goodfeats_temp
            }
            
            
            
            if(output.device.type!="pdf"){
              
              try(dev.off(),silent=TRUE)
            }
            
            if(output.device.type!="pdf"){
              
              temp_filename_1<-"Figures/Boxplots.selectedfeats.raw.pdf"
              if(boxplot.type=="simple"){
                  pdf(temp_filename_1,height=plots.height,width=plots.width)
              }
            }
            
        if(plot.boxplots.raw==TRUE){
            
      #    save(goodfeats_raw,goodfeats_temp,classlabels_raw_boxplots,classlabels_orig,
       #          output_dir,boxplot.col.opt,cex.plots,ylab_text,boxplot.type,ylab_text_raw,
        #        analysistype,multiple.figures.perpanel,alphabetical.order,goodfeats_name,plots.height,plots.width,file="boxplotrawdebug.Rda")
            
            par(mfrow=c(1,1),family="sans",cex=cex.plots)
                        par(mfrow=c(1,1),family="sans",cex=cex.plots)
            
            #get_boxplots(X=goodfeats_raw,Y=classlabels_raw_boxplots,parentoutput_dir=output_dir,boxplot.col.opt=boxplot.col.opt,alphacol=0.3,newdevice=FALSE,cex.plots=cex.plots,ylabel=" Intensity",name=goodfeats_name,add.pvalues=add.pvalues,
                     #    add.jitter=add.jitter,alphabetical.order=alphabetical.order,boxplot.type=boxplot.type,study.design=analysistype)
            plot.ylab_text1=paste("",ylab_text,sep="")
            
            if(pairedanalysis==TRUE){
              
              #classlabels_paired<-cbind(classlabels[,1],subject_inf,classlabels[,2])
              
              get_boxplots(X=goodfeats_raw,Y=classlabels_orig,parentoutput_dir=output_dir,boxplot.col.opt=boxplot.col.opt,
                           newdevice=FALSE,cex.plots=cex.plots,ylabel=ylab_text_raw,name=goodfeats_name,add.pvalues=add.pvalues,add.jitter=add.jitter,
                           alphabetical.order=alphabetical.order,boxplot.type=boxplot.type,
                           study.design=gsub(analysistype,pattern="repeat",replacement=""),multiple.figures.perpanel=multiple.figures.perpanel,numnodes=num_nodes,
                           plot.height = plots.height,plot.width=plots.width,
                           filename="Figures/Boxplots.selectedfeats.raw",alphacol = alpha.col,ggplot.type1=ggplot.type1,facet.nrow=facet.nrow)
              
            
            }else{
            get_boxplots(X=goodfeats_raw,Y=classlabels_orig,parentoutput_dir=output_dir,boxplot.col.opt=boxplot.col.opt,
                         newdevice=FALSE,cex.plots=cex.plots,ylabel=ylab_text_raw,name=goodfeats_name,add.pvalues=add.pvalues,add.jitter=add.jitter,
                         alphabetical.order=alphabetical.order,boxplot.type=boxplot.type,
                         study.design=analysistype,multiple.figures.perpanel=multiple.figures.perpanel,numnodes=num_nodes,plot.height = plots.height,plot.width=plots.width,
                         filename="Figures/Boxplots.selectedfeats.raw",alphacol = alpha.col,ggplot.type1=ggplot.type1,facet.nrow=facet.nrow)
            
            }
            
            #try(dev.off(),silent=TRUE)
            
            if(output.device.type!="pdf"){
              
              try(dev.off(),silent=TRUE)
            }
            
        }  
            if(FALSE)
              {
              if(output.device.type!="pdf"){
                
                temp_filename_1<-"Figures/Barplots_selectedfeats.pdf"
                
                #png(temp_filename_1,width=plots.width,height=plots.height,res=plots.res,type=plots.type,units="in")
                #pdf(temp_filename_1,bg="transparent") #, height = 5.5, width = 3)
                pdf(temp_filename_1,width=plots.width,height=plots.height)
              }
              
              
              
              plot(0:10, type = "n", xaxt="n", yaxt="n", bty="n", xlab = "", ylab = "")
              
              
              text(5, 8, "Barplots of selected features using the\n normalized intensities/adundance levels")
              
              par(mfrow=c(1,1),family="sans",cex=cex.plots,pty="s")
              try(get_barplots(feature_table_file,class_labels_file,X=goodfeats_temp,Y=classlabels_orig,parentoutput_dir=output_dir
                               ,newdevice=FALSE,ylabel=ylab_text,cex.val=cex.plots,barplot.col.opt=barplot.col.opt,error.bar=error.bar),silent=TRUE)
              
              ###savelist=ls(),file="getbarplots.Rda")
              if(featselmethod=="limma2way" | featselmethod=="limma2wayrepeat" | featselmethod=="pls2wayrepeat" | featselmethod=="spls2wayrepeat" | featselmethod=="pls2way" | featselmethod=="spls2way" | featselmethod=="lm2wayanova" | featselmethod=="lm2wayanovarepeat")
              {
                #if(ggplot.type1==TRUE){
                  barplot.xaxis="Factor2"
               # }else{
                  
                  
               # }
                
              }
              
              get_barplots(feature_table_file,class_labels_file,X=goodfeats_temp,Y=classlabels_orig,parentoutput_dir=output_dir,
                           newdevice=FALSE,ylabel=plot.ylab_text,cex.plots=cex.plots,barplot.col.opt=barplot.col.opt,error.bar=error.bar,
                           barplot.xaxis=barplot.xaxis,alphabetical.order=alphabetical.order,name=goodfeats_name,study.design=analysistype)
              
              if(output.device.type!="pdf"){
                
                try(dev.off(),silent=TRUE)
              }
              
                      if(FALSE){    
                          if(output.device.type!="pdf"){
                            
                            temp_filename_1<-"Figures/Individual_sample_plots_selectedfeats.pdf"
                            
                            #png(temp_filename_1,width=plots.width,height=plots.height,res=plots.res,type=plots.type,units="in")
                            #pdf(temp_filename_1)
                            
                            pdf(temp_filename_1,width=plots.width,height=plots.height)
                          }
                          
                          #  par(mfrow=c(2,2))
                          par(mfrow=c(1,1),family="sans",cex=cex.plots)
                          #try(get_individualsampleplots(feature_table_file,class_labels_file,X=goodfeats_temp,Y=classlabels_orig,parentoutput_dir=output_dir,newdevice=FALSE,ylabel=ylab_text,cex.val=cex.plots,sample.col.opt=sample.col.opt),silent=TRUE)
                          get_individualsampleplots(feature_table_file,class_labels_file,X=goodfeats_temp,Y=classlabels_orig,parentoutput_dir=output_dir,newdevice=FALSE,ylabel=ylab_text,cex.plots=cex.plots,sample.col.opt=individualsampleplot.col.opt,alphabetical.order=alphabetical.order,name=goodfeats_name)
                          
                          
                          
                          if(output.device.type!="pdf"){
                            
                            try(dev.off(),silent=TRUE)
                          }
                      }
              
            }
            if(globalclustering==TRUE){
              
              print("Performing global clustering using EM")
              
              if(output.device.type!="pdf"){
                
                temp_filename_1<-"Figures/GlobalclusteringEM.png"
                
                png(temp_filename_1,width=plots.width,height=plots.height,res=plots.res,type=plots.type,units="in")
              }
              m1<-Mclust(t(data_m_fc_withfeats[,-c(1:2)]))
              s1<-m1$classification #summary(m1)
              
              EMcluster<-m1$classification
              
              col_vec <- colorRampPalette(brewer.pal(10, "RdBu"))(length(levels(as.factor(classlabels_orig[,2]))))
              #col_vec<-topo.colors(length(levels(as.factor(classlabels_orig[,2])))) #patientcolors #heatmap_cols[1:length(levels(classlabels_orig[,2]))]
              t1<-table(EMcluster,classlabels_orig[,2])
              
              par(mfrow=c(1,1))
              plot(t1,col=col_vec,main="EM cluster labels\n using all features",cex.axis=1,ylab="Class",xlab="Cluster number")
              
              par(xpd=TRUE)
              try(legend("bottomright",legend=levels(classlabels_orig[,2]),text.col=col_vec,pch=13,cex=0.4),silent=TRUE)
              
              par(xpd=FALSE)
              
            #  save(m1,EMcluster,classlabels_orig,file="EMres.Rda")
              
              t1<-cbind(EMcluster,classlabels_orig[,2])
              
              write.table(t1,file="Tables/EM_clustering_labels_using_allfeatures.txt",sep="\t")
              
              if(output.device.type!="pdf"){
                
                try(dev.off(),silent=TRUE)
              }
              
              print("Performing global clustering using HCA")
              
              if(output.device.type!="pdf"){
                
                temp_filename_1<-"Figures/GlobalclusteringHCA.png"
                
                png(temp_filename_1,width=plots.width,height=plots.height,res=plots.res,type=plots.type,units="in")
              }
              
              #if(FALSE)
              {
                
                #p1<-heatmap.2(as.matrix(data_m_fc_withfeats[,-c(1:2)]),scale="row",symkey=FALSE,col=topo.colors(n=256))
                
                if(heatmap.col.opt=="RdBu"){
                  
                  heatmap.col.opt="redblue"
                }
                
                heatmap_cols <- colorRampPalette(brewer.pal(10, "RdBu"))(256)
                heatmap_cols<-rev(heatmap_cols)
                
                if(heatmap.col.opt=="topo"){
                  heatmap_cols<-topo.colors(256)
                  heatmap_cols<-rev(heatmap_cols)
                }else
                {
                  if(heatmap.col.opt=="heat"){
                    heatmap_cols<-heat.colors(256)
                    heatmap_cols<-rev(heatmap_cols)
                  }else{
                    
                    if(heatmap.col.opt=="yellowblue"){
                      
                      heatmap_cols<-colorRampPalette(c("yellow","blue"))(256) #colorRampPalette(c("yellow","white","blue"))(256)
                      #heatmap_cols<-blue2yellow(256) #colorRampPalette(c("yellow","blue"))(256)
                      heatmap_cols<-rev(heatmap_cols)
                    }else{
                      
                      if(heatmap.col.opt=="redblue"){
                        
                        heatmap_cols <- colorRampPalette(brewer.pal(10, "RdBu"))(256)
                        heatmap_cols<-rev(heatmap_cols)
                      }else{
                        
                        #my_palette <- colorRampPalette(c("red", "yellow", "green"))(n = 299)
                        if(heatmap.col.opt=="redyellowgreen"){
                          
                          heatmap_cols <- colorRampPalette(c("red", "yellow", "green"))(n = 299)
                          heatmap_cols<-rev(heatmap_cols)
                        }else{
                          if(heatmap.col.opt=="yellowwhiteblue"){
                            
                            heatmap_cols<-colorRampPalette(c("yellow2","white","blue"))(256) #colorRampPalette(c("yellow","white","blue"))(256)
                            heatmap_cols<-rev(heatmap_cols)
                          }else{
                            
                            if(heatmap.col.opt=="redwhiteblue"){
                              
                              heatmap_cols<-colorRampPalette(c("red","white","blue"))(256) #colorRampPalette(c("yellow","white","blue"))(256)
                              heatmap_cols<-rev(heatmap_cols)
                            }else{
                              
                              
                              
                              heatmap_cols <- colorRampPalette(brewer.pal(10, heatmap.col.opt))(256)
                              heatmap_cols<-rev(heatmap_cols)
                              
                            }
                            
                          }
                          
                        }
                        
                      }
                      
                    }
                  }
                  
                  
                }
                
                
                
                #col_vec<-heatmap_cols[1:length(levels(classlabels_orig[,2]))]
                c1<-WGCNA::cor(as.matrix(data_m_fc_withfeats[,-c(1:2)]),method=cor.method,use="pairwise.complete.obs") #cor(d1[,-c(1:2)])
                d2<-as.dist(1-c1)
                clust1<-hclust(d2)
                
                hr <- try(hclust(as.dist(1-WGCNA::cor(t(data_m_fc_withfeats),method=cor.method,use="pairwise.complete.obs"))),silent=TRUE) #metabolites
                #hc <- try(hclust(as.dist(1-WGCNA::cor(data_m,method=cor.method,use="pairwise.complete.obs"))),silent=TRUE) #samples
                
                
                h73<-heatmap.2(as.matrix(data_m_fc_withfeats[,-c(1:2)]), Rowv=as.dendrogram(hr), Colv=as.dendrogram(clust1),  
                               col=heatmap_cols, scale="row",key=TRUE, symkey=FALSE, density.info="none", trace="none",
                               cexRow=1, cexCol=1,xlab="",ylab="", main="Global clustering\n using all features",
                               ColSideColors=patientcolors,labRow = FALSE, labCol = FALSE)
                
                # par(xpd=TRUE)
                #legend("bottomleft",legend=levels(classlabels_orig[,2]),text.col=unique(patientcolors),pch=13,cex=0.4)
                #par(xpd=FALSE)
                
                clust_res<-cutreeDynamic(distM=as.matrix(d2),dendro=clust1)
                
                #mycl_samples <- cutree(clust1, h=max(clust1$height)/2)
                
                HCAcluster<-clust_res
                
                c2<-cbind(clust1$labels,HCAcluster)
                
                rownames(c2)<-c2[,1]
                
                c2<-as.data.frame(c2)
                
                t1<-table(HCAcluster,classlabels_orig[,2])
                
                plot(t1,col=col_vec,main="HCA (Cutree Dynamic) cluster labels\n using all features",cex.axis=1,ylab="Class",xlab="Cluster number")
                
                par(xpd=TRUE)
                try(legend("bottomright",legend=levels(classlabels_orig[,2]),text.col=col_vec,pch=13,cex=0.4),silent=TRUE)
                par(xpd=FALSE)
                
                t1<-cbind(HCAcluster,classlabels_orig[,2])
                write.table(t1,file="Tables/HCA_clustering_labels_using_allfeatures.txt",sep="\t")
                
          
              }
              
              
              
              if(output.device.type!="pdf"){
                
                try(dev.off(),silent=TRUE)
              }
            }
            
          }
          #dev.off()
        }
        else
        {
          #goodfeats_allfields<-as.data.frame(goodfeats)
          
          goodfeats<-goodfeats[,-c(1:time_ind)]
          
          
        }
        
      }
      
      if(length(goodip)>0){
        
        try(dev.off(),silent=TRUE)
      }
    }
    else{
      try(dev.off(),silent=TRUE)
      break;
    }
    
    
    if(analysismode=="classification" & WGCNAmodules==TRUE){
      classlabels_temp<-classlabels_orig_wgcna #cbind(classlabels_sub[,1],classlabels)
      #print(classlabels_temp)
      data_temp<-data_matrix_beforescaling[,-c(1:2)]
      
      
      cl<-makeCluster(num_nodes)
      
      #clusterExport(cl,"do_rsd")
      #feat_rsds<-parApply(cl,data_temp,1,do_rsd)
      #rm(data_temp)
      #feat_rsds<-abs(feat_rsds) #round(max_rsd,2)
      #print(summary(feat_rsds))
      #if(length(which(feat_rsds>0))>0)
      {
        X<-data_m_fc_withfeats #data_matrix[which(feat_rsds>=wgcnarsdthresh),]
        
        #	print(head(X))
        #		print(dim(X))
        
        
        if(output.device.type!="pdf"){
          
          temp_filename_1<-"Figures/WGCNA_preservation_plot.png"
          
          png(temp_filename_1,width=plots.width,height=plots.height,res=plots.res,type=plots.type,units="in")
        }
        
        #  #save(X,classlabels_temp,data_m_fc_withfeats,goodip,file="wgcna.Rda")
        
        print("Performing WGCNA: generating preservation plot")
        #preservationres<-try(do_wgcna(X=X,Y=classlabels_temp,sigfeats=data_m_fc_withfeats[goodip,c(1:2)]),silent=TRUE)
        #pres<-try(do_wgcna(X=X,Y=classlabels_temp,sigfeats=data_m_fc_withfeats[goodip,c(1:2)]),silent=TRUE)
        pres<-try(do_wgcna(X=X,Y=classlabels_temp,sigfeats=data_m_fc_withfeats[goodip,c(1:2)]),silent=TRUE)
        
        #pres<-do_wgcna(X=X,Y=classlabels_temp,sigfeats=data_m_fc_withfeats[goodip,c(1:2)]) #,silent=TRUE)
        
        if(is(pres,"try-error")){
          
          print("WGCNA could not be performed. Error: ")
          print(pres)
        }
        
        if(output.device.type!="pdf"){
          
          try(dev.off(),silent=TRUE)
        }
      }
    }	
    #print(lf)
    
    #print("next iteration")
    #dev.off()
  }
  
  setwd(parentoutput_dir)
  summary_res<-cbind(log2.fold.change.thresh_list,feat_eval,feat_sigfdrthresh,feat_sigfdrthresh_cv,feat_sigfdrthresh_permut,res_score_vec)
  
  if(fdrmethod=="none"){
    exp_fp<-round(fdrthresh*feat_eval)
  }else{
    exp_fp<-round(fdrthresh*feat_sigfdrthresh)
  }
  rank_num<-order(summary_res[,5],decreasing=TRUE)
  
  ##save(allmetabs_res,file="allmetabs_res.Rda")
  
  
  if(featselmethod=="limma" | featselmethod=="limma2way" | featselmethod=="limma2wayrepeat" | featselmethod=="lmreg" | featselmethod=="logitreg" 
     | featselmethod=="lm2wayanova" | featselmethod=="lm1wayanova" | featselmethod=="lm1wayanovarepeat" | featselmethod=="lm2wayanovarepeat" | featselmethod=="wilcox" | featselmethod=="ttest" | featselmethod=="poissonreg" | featselmethod=="limma1wayrepeat" | featselmethod=="lmregrepeat")
  {
    summary_res<-cbind(summary_res,exp_fp)
    
    #print("HERE13134")
    type.statistic="pvalue"
    if(length(allmetabs_res)>0){	
      
      #stat_val<-(-1)*log10(allmetabs_res[,4])
      
      
      stat_val<-allmetabs_res[,4]
      
    }
    
    
    colnames(summary_res)<-c("RSD.thresh","Number of features left after RSD filtering","Number of features selected",paste(pred.eval.method,"-accuracy",sep=""),paste(pred.eval.method," permuted accuracy",sep=""),"Score","Expected_False_Positives")
  }else{
    #exp_fp<-round(fdrthresh*feat_sigfdrthresh)
    
    #if(featselmethod=="MARS" | featselmethod=="RF" | featselmethod=="pls" | featselmethod=="o1pls" | featselmethod=="o2pls"){
    
    exp_fp<-rep(NA,dim(summary_res)[1])
    #}
  #  print("HERE13135")
    if(length(allmetabs_res)>0){
      stat_val<-(allmetabs_res[,4])
    }
    type.statistic="other"
    summary_res<-cbind(summary_res,exp_fp)
    colnames(summary_res)<-c("RSD.thresh","Number of features left after RSD filtering","Number of features selected",paste(pred.eval.method,"-accuracy",sep=""),paste(pred.eval.method," permuted accuracy",sep=""),"Score","Expected_False_Positives")
    
  }
  
  featselmethod<-parentfeatselmethod
  file_name<-paste(parentoutput_dir,"/Results_summary_",featselmethod,".txt",sep="")
  write.table(summary_res,file=file_name,sep="\t",row.names=FALSE)
  
  
  if(output.device.type=="pdf"){
    
    try(dev.off(),silent=TRUE)
  }
  
  print("##############Level 1: processing complete###########")
  
  
  if(length(best_feats)>1)
  {
    
    mz_index<-best_feats
    
    #par(mfrow=c(1,1),family="sans",cex=cex.plots)
    #                get_boxplots(X=goodfeats_raw,Y=classlabels_orig,parentoutput_dir=output_dir,boxplot.col.opt=boxplot.col.opt,alphacol=0.3,newdevice=FALSE,cex=cex.plots,ylabel="raw Intensity",name=goodfeats_name,add.pvalues=add.pvalues,add.jitter=add.jitter,boxplot.type=boxplot.type)
    
    setwd(output_dir)
    
    ###save(goodfeats,goodfeats_temp,classlabels_orig,classlabels_response_mat,output_dir,xlab_text,ylab_text,goodfeats_name,file="debugscatter.Rda")
    
    if(analysismode=="regression"){
      pdf("Figures/Scatterplots.pdf")
      
      if(is.na(xlab_text)==TRUE){
        
        xlab_text=""
      }
      
    #  save(goodfeats_temp,classlabels_orig,output_dir,ylab_text,xlab_text,goodfeats_name,cex.plots,scatterplot.col.opt,file="scdebug.Rda")
      get_scatter_plots(X=goodfeats_temp,Y=classlabels_orig,parentoutput_dir=output_dir,newdevice=FALSE,ylabel=ylab_text,xlabel=xlab_text,
                        name=goodfeats_name,cex.plots=cex.plots,scatterplot.col.opt=scatterplot.col.opt)
      dev.off()
    }
    setwd(parentoutput_dir)
    
    
    
    
    if(analysismode=="classification"){
      
      
      
      log2.fold.change.thresh=log2.fold.change.thresh_list[best_logfc_ind]
      
      print(paste("Best results found at RSD threshold ", log2.fold.change.thresh,sep=""))
      
      print(best_acc)
      #print(paste(kfold,"-fold CV accuracy ", best_acc,sep=""))
      if(pred.eval.method=="CV"){
        
        print(paste(kfold,"-fold CV accuracy: ", best_acc,sep=""))
      }else{
        if(pred.eval.method=="AUC"){
          
          print(paste("Area under the curve (AUC) is : ", best_acc,sep=""))
        }
      }
      
      
    #  data_m<-parent_data_m
     # data_m_fc<-data_m #[which(abs(mean_groups)>log2.fold.change.thresh),]
      
      data_m_fc_withfeats<-data_matrix[,c(1:2)]
      data_m_fc_withfeats<-cbind(data_m_fc_withfeats,data_m_fc)
      
      
      #when using a feature table generated by apLCMS
      
      rnames<-paste("mzid_",seq(1,dim(data_m_fc)[1]),sep="")
      #print(best_limma_res[1:3,])
      
      goodfeats<-best_limma_res[order(best_limma_res$mz),-c(1:2)]
      
      #goodfeats<-best_limma_res[,-c(1:2)]
      
      goodfeats_all<-goodfeats
      
      goodfeats<-goodfeats_all
      rm(goodfeats_all)
    }
    
    try(unlink("Rplots.pdf"),silent=TRUE)
    
    
    if(globalcor==TRUE){
      
      
      if(length(best_feats)>2){
        if(is.na(abs.cor.thresh)==FALSE){
          #setwd(parentoutput_dir)
          print("##############Level 2: Metabolome wide correlation network analysis of differentially expressed metabolites###########")
          print(paste("Generating metabolome-wide ",cor.method," correlation network using RSD threshold ", log2.fold.change.thresh," results",sep=""))
          print(parentoutput_dir)
          print(output_dir)
          setwd(output_dir)
          data_m_fc_withfeats<-as.data.frame(data_m_fc_withfeats)
          
          goodfeats<-as.data.frame(goodfeats)
          #print(goodfeats[1:4,])
          sigfeats_index<-which(data_m_fc_withfeats$mz%in%goodfeats$mz)
          sigfeats<-sigfeats_index
          if(globalcor==TRUE){
            
            #outloc<-paste(parentoutput_dir,"/Allcornetworksigfeats","log2fcthresh",log2.fold.change.thresh,"/",sep="")
            #outloc<-paste(parentoutput_dir,"/Stage2","/",sep="")
            
            #dir.create(outloc)
            #setwd(outloc)
            
            #dir.create("CorrelationAnalysis")
            #setwd("CorrelationAnalysis")
            if(networktype=="complete"){
              
           #   save(data_matrix,sigfeats_index,output_dir,max.cor.num,net_node_colors,net_legend,cor.method,abs.cor.thresh,cor.fdrthresh,file="r1.Rda")
              mwan_fdr<-try(do_cor(data_matrix,subindex=sigfeats_index,targetindex=NA,outloc=output_dir,networkscope="global",cor.method,abs.cor.thresh,cor.fdrthresh,
                               max.cor.num,net_node_colors,net_legend,newdevice=TRUE),silent=TRUE)
            }else{
              if(networktype=="GGM"){
                mwan_fdr<-try(get_partial_cornet(data_matrix, sigfeats.index=sigfeats_index,targeted.index=NA,networkscope="global",
                                             cor.method,abs.cor.thresh,cor.fdrthresh,outloc=output_dir,net_node_colors,net_legend),silent=TRUE)
              }else{
                print("Invalid option. Please use complete or GGM.")
              }
            }
            
            print("##############Level 2: processing complete###########")
          }else{
            print("##############Skipping Level 2: global correlation analysis###########")
            
          }
          
          
          #temp_data_m<-cbind(allmetabs_res[,c("mz","time")],stat_val)
          
          
          
          
          if(analysismode=="classification"){
            #  classlabels_temp<-cbind(classlabels_sub[,1],classlabels)
            #do_wgcna(X=data_matrix,Y=classlabels,sigfeats.index=sigfeats_index)
          }
          
          print("##############Level 3: processing complete###########")
          print("#########################")
        }
        
        
      }
    }
    else{
      print(paste("Can not perform network analysis. Too few metabolites.",sep=""))
    }
  }
  
  if(FALSE){
    if(length(featselmethod)>1){
      abs.cor.thresh=NA
      globalcor=FALSE
    }
  }
  
  
  ###save(stat_val,allmetabs_res,check_names,metab_annot,kegg_species_code,database,reference_set,type.statistic,file="fcsdebug.Rda")
  
  
  setwd(output_dir)
  
  unlink("fdrtoolB.pdf",force=TRUE)
  
  if(is.na(target.data.annot)==FALSE){
    
    #dir.create("NetworkAnalysis")
    #setwd("NetworkAnalysis")
    
    colnames(target.data.annot)<-c("mz","time","KEGGID")
    if(length(check_names)<1){
      
      allmetabs_res<-cbind(stat_val,allmetabs_res)
      metab_data<-merge(allmetabs_res,target.data.annot,by=c("mz","time"))
      dup.feature.check=TRUE
    }else{
      
      allmetabs_res_withnames<-cbind(stat_val,allmetabs_res_withnames)
      metab_data<-merge(allmetabs_res_withnames,target.data.annot,by=c("Name"))
      dup.feature.check=FALSE
    }
    
    ###save(stat_val,allmetabs_res,check_names,metab_annot,kegg_species_code,database,metab_data,reference_set,type.statistic,file="fcsdebug.Rda")
    
    
    if(length(check_names)<1){
      metab_data<-metab_data[,c("KEGGID","stat_val","mz","time")]
      colnames(metab_data)<-c("KEGGID","Statistic","mz","time")
      
    }else{
      
      metab_data<-metab_data[,c("KEGGID","stat_val")]
      colnames(metab_data)<-c("KEGGID","Statistic")
    }
    
    # ##save(metab_annot,kegg_species_code,database,metab_data,reference_set,type.statistic,file="fcsdebug.Rda")
    
    
    #metab_data: KEGGID, Statistic
    fcs_res<-get_fcs(kegg_species_code=kegg_species_code,database=database,target.data=metab_data,target.data.annot=target.data.annot,reference_set=reference_set,type.statistic=type.statistic,fcs.min.hits=fcs.min.hits)
    
    ###save(fcs_res,file="fcs_res.Rda")
    
    write.table(fcs_res,file="Tables/functional_class_scoring.txt",sep="\t",row.names=TRUE)	
    if(length(fcs_res)>0){
      if(length(which(fcs_res$pvalue<pvalue.thresh))>10){
        
        fcs_res_filt<-fcs_res[which(fcs_res$pvalue<pvalue.thresh)[1:10],]
      }else{
        fcs_res_filt<-fcs_res[which(fcs_res$pvalue<pvalue.thresh),]
      }
      fcs_res_filt<-fcs_res_filt[order(fcs_res_filt$pvalue,decreasing=FALSE),]
      
      fcs_res_filt$Name<-gsub(as.character(fcs_res_filt$Name),pattern=" - Homo sapiens \\(human\\)",replacement="")
      
      
      fcs_res_filt$pvalue=(-1)*log10(fcs_res_filt$pvalue)
      
      fcs_res_filt<-fcs_res_filt[order(fcs_res_filt$pvalue,decreasing=FALSE),]
      
      
      print(Sys.time())
      
      p=ggbarplot(fcs_res_filt,x="Name",y="pvalue",orientation="horiz",ylab="(-)log10pvalue",xlab="",color="orange",fill="orange",title=paste("Functional classes significant at p<",pvalue.thresh," threhsold",sep=""))
      p=p+font("title",size=10)
      p=p+font("x.text",size=10)
      p=p+font("y.text",size=10)
      p=p + geom_hline(yintercept = (-1)*log10(pvalue.thresh), linetype="dotted",size=0.7)
      print(Sys.time())
      
      pdf("Figures/Functional_Class_Scoring.pdf")
      print(p)
      dev.off()
    }
    
    
    
    
    
    print(paste(featselmethod, " processing done.",sep=""))
    
  }
  
  
  setwd(parentoutput_dir)
  
  
  
  
  #print("Note A: Please note that log2 fold-change based filtering is only applicable to two-class comparison. 
  #log2fcthresh of 0 will remove only those features that have exactly sample mean intensities between the two groups.
  #More features will be filtered prior to FDR as log2fcthresh increases.")
  
  #print("Note C: Please make sure all the packages are installed. You can use the command install.packages(packagename) to install packages.")
  #print("Eg: install.packages(\"mixOmics\"),install.packages(\"snow\"), install.packages(\"e1071\"), biocLite(\"limma\"), install.packages(\"gplots\").")
  #print("Eg: install.packages("mixOmics""),install.packages("snow"), install.packages("e1071"), biocLite("limma"), install.packages("gplots").")
  ##############################
  ##############################
  ###############################
  
  
  
  if(length(best_feats)>0){
    
    goodfeats<-as.data.frame(goodfeats)
    #goodfeats<-data_matrix_beforescaling[which(data_matrix_beforescaling$mz%in%goodfeats$mz),]
  }else{
    goodfeats-{}
  }
  
  cur_date<-Sys.time()
  cur_date<-gsub(x=cur_date,pattern="-",replacement="")
  cur_date<-gsub(x=cur_date,pattern=":",replacement="")
  cur_date<-gsub(x=cur_date,pattern=" ",replacement="")
  if(saveRda==TRUE){
    fname<-paste("Analysis_",featselmethod,"_",cur_date,".Rda",sep="")
    ###savelist=ls(),file=fname)
  }
  
  ################################
  
  fname_del<-paste(output_dir,"/Rplots.pdf",sep="")
  try(unlink(fname_del),silent=TRUE)
  
  
  
  if(removeRda==TRUE)
  {
    unlink("*.Rda",force=TRUE,recursive=TRUE)
    #unlink("pairwise_results/*.Rda",force=TRUE,recursive=TRUE)
    
  }
  
  return(list("diffexp_metabs"=goodfeats_allfields,  "mw.an.fdr"=mwan_fdr,"targeted.an.fdr"=targetedan_fdr,"classlabels"=classlabels_orig,"all_metabs"=allmetabs_res_withnames))
  
  
}
