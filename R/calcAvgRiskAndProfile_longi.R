calcAvgRiskAndProfile_longi<-function(clusObj,includeFixedEffects=F,proportionalHazards=F,nSweeps){

  clusObjRunInfoObj=NULL
  directoryPath=NULL
  fileStem=NULL
  xModel=NULL
  nCategories=NULL
  varSelect=NULL
  varSelectType=NULL
  includeResponse=NULL
  nFixedEffects=NULL
  nFixedEffects_clust=NULL
  reportBurnIn=NULL
  nBurn=NULL
  nFilter=NULL
  nSweeps=NULL
  nRandomEffects=NULL
  nClusters=NULL
  clustering=NULL
  nCategoriesY=NULL
  nCovariates=NULL
  nDiscreteCovs=NULL
  nContinuousCovs=NULL
  nSubjects=NULL
  nPredictSubjects=NULL
  yModel=NULL
  wMat=NULL
  yMat=NULL
  weibullFixedShape=NULL
  sampleGPmean=FALSE
  nFixedEffects_clust=0

  #library("PReMiuM")
  for (i in 1:length(clusObj)) assign(names(clusObj)[i],clusObj[[i]])
  for (i in 1:length(clusObjRunInfoObj)) assign(names(clusObjRunInfoObj)[i],clusObjRunInfoObj[[i]])

  assign(names(clusObjRunInfoObj)[4],min(clusObjRunInfoObj[[4]],nSweeps))
  # Construct the number of clusters file name
  nClustersFileName<-file.path(directoryPath,paste(fileStem,'_nClusters.txt',sep=''))
  nClustersFile<-file(nClustersFileName,open="r")

  # Construct the allocation file name
  zFileName <- file.path(directoryPath,paste(fileStem,'_z.txt',sep=''))
  zFile<-file(zFileName,open="r")

  if(xModel=="Discrete"){
    # Construct the allocation file name
    phiFileName <- file.path(directoryPath,paste(fileStem,'_phi.txt',sep=''))
    phiFile<-file(phiFileName,open="r")
    # Get the maximum number of categories
    maxNCategories<-max(nCategories)
    if(varSelect){
      nullPhiFileName <- file.path(directoryPath,paste(fileStem,'_nullPhi.txt',sep=''))
      if(varSelectType=="Continuous"){
        gammaFileName <- file.path(directoryPath,paste(fileStem,'_rho.txt',sep=''))
        gammaFile<-file(gammaFileName,open="r")
      }else{
        gammaFileName <- file.path(directoryPath,paste(fileStem,'_gamma.txt',sep=''))
        gammaFile<-file(gammaFileName,open="r")
      }
    }
  }else if(xModel=="Normal"){
    muFileName <- file.path(directoryPath,paste(fileStem,'_mu.txt',sep=''))
    muFile<-file(muFileName,open="r")
    SigmaFileName <- file.path(directoryPath,paste(fileStem,'_Sigma.txt',sep=''))
    SigmaFile<-file(SigmaFileName,open="r")
    if(varSelect){
      nullMuFileName <- file.path(directoryPath,paste(fileStem,'_nullMu.txt',sep=''))
      if(varSelectType=="Continuous"){
        gammaFileName <- file.path(directoryPath,paste(fileStem,'_rho.txt',sep=''))
        gammaFile<-file(gammaFileName,open="r")
      }else{
        gammaFileName <- file.path(directoryPath,paste(fileStem,'_gamma.txt',sep=''))
        gammaFile<-file(gammaFileName,open="r")
      }
    }
  }else if(xModel=="Mixed"){
    # Construct the allocation file name
    phiFileName <- file.path(directoryPath,paste(fileStem,'_phi.txt',sep=''))
    phiFile<-file(phiFileName,open="r")
    # Get the maximum number of categories
    maxNCategories<-max(nCategories)
    muFileName <- file.path(directoryPath,paste(fileStem,'_mu.txt',sep=''))
    muFile<-file(muFileName,open="r")
    SigmaFileName <- file.path(directoryPath,paste(fileStem,'_Sigma.txt',sep=''))
    SigmaFile<-file(SigmaFileName,open="r")
    if(varSelect){
      nullPhiFileName <- file.path(directoryPath,paste(fileStem,'_nullPhi.txt',sep=''))
      nullMuFileName <- file.path(directoryPath,paste(fileStem,'_nullMu.txt',sep=''))
      if(varSelectType=="Continuous"){
        gammaFileName <- file.path(directoryPath,paste(fileStem,'_rho.txt',sep=''))
        gammaFile<-file(gammaFileName,open="r")
      }else{
        gammaFileName <- file.path(directoryPath,paste(fileStem,'_gamma.txt',sep=''))
        gammaFile<-file(gammaFileName,open="r")
      }
    }
  }

  if(includeResponse){
    # Construct the theta file name
    if(!is.element(yModel, c("Longitudinal","LME" ))){
      thetaFileName <- file.path(directoryPath,paste(fileStem,'_theta.txt',sep=''))
      thetaFile<-file(thetaFileName,open="r")
    }
    if (yModel=="Survival"&&!weibullFixedShape){
      nuFileName <- file.path(directoryPath,paste(fileStem,'_nu.txt',sep=''))
      nuFile<-file(nuFileName,open="r")
    }else if(yModel=="Longitudinal"){
      LFileName <- file.path(directoryPath,paste(fileStem,'_L.txt',sep=''))
      LFile<-file(LFileName,open="r")
      if(sampleGPmean){
        GPmeanFileName <- file.path(directoryPath,paste(fileStem,'_meanGP.txt',sep=''))
        GPmeanFile<-file(GPmeanFileName,open="r")
      }
    }else if(yModel=="MVN"){
      MVNmuFileName <- file.path(directoryPath,paste(fileStem,'_MVNmu.txt',sep=''))
      MVNmuFile<-file(MVNmuFileName,open="r")
      MVNSigmaFileName <- file.path(directoryPath,paste(fileStem,'_MVNSigma.txt',sep=''))
      MVNSigmaFile<-file(MVNSigmaFileName,open="r")
    }else if(yModel=="LME"){
      covREFileName <- file.path(directoryPath,paste(fileStem,'_cov_RandomEffects_LME.txt',sep=''))
      covREFile<-file(covREFileName,open="r")
      SigmaLMEFileName <- file.path(directoryPath,paste(fileStem,'_epsilonLME.txt',sep=''))
      SigmaLMEFile<-file(SigmaLMEFileName,open="r")
      RE_LMEFileName <- file.path(directoryPath,paste(fileStem,'_RandomEffectsLME.txt',sep=''))
      RE_LMEFile<-file(RE_LMEFileName,open="r")
    }
    if(nFixedEffects>0){
      # Construct the fixed effect coefficient file name
      betaFileName <-file.path(directoryPath,paste(fileStem,'_beta.txt',sep=''))
      betaFile<-file(betaFileName,open="r")
    }
    if(nFixedEffects_clust>0){
      # Construct the fixed effect coefficient file name
      betamixFileName <-file.path(directoryPath,paste(fileStem,'_beta_cluster-specific.txt',sep=''))
      betamixFile<-file(betamixFileName,open="r")
    }
    if (yModel=="Survival"&&weibullFixedShape){
      nuFileName<-file.path(directoryPath,paste(fileStem,'_nu.txt',sep=''))
      nuFile<-file(nuFileName,open="r")
      nu<-(read.table(nuFile)[,1])
      close(nuFile)
    }
  }

  # Restrict to sweeps after burn in
  firstLine<-floor(ifelse(reportBurnIn,nBurn/nFilter+2,1)) #AR
  lastLine<-floor((nSweeps+ifelse(reportBurnIn,nBurn+1,0))/nFilter) #AR
  nSamples<-lastLine-firstLine+1

  # Make a list of the subjects in each of the optimal clusters
  optAlloc<-vector("list",nClusters)
  for(c in 1:nClusters){
    optAlloc[[c]]<-which(clustering==c)
  }

  if(includeResponse){
    # Initialise the object for storing the risks
    riskArray<-array(0,dim=c(nSamples,nClusters,nCategoriesY))
    thetaArray<-array(0,dim=c(nSamples,nClusters,nCategoriesY))
    if (yModel=="Survival"&&!weibullFixedShape){
      nuArray<-array(0,dim=c(nSamples,nClusters))
    } else if(yModel=="Longitudinal"){
      LArray<-array(0,dim=c(nSamples,nClusters,3))
      if(sampleGPmean)
        GPmeanArray<-array(0,dim=c(nSamples,nClusters,nTimes_unique))
      riskArray<-array(0,dim=c(nSamples,nClusters))
    } else if(yModel=="MVN"){
      MVNmuArray<-array(0,dim=c(nSamples,nClusters,nOutcomes))
      MVNSigmaArray<-array(0,dim=c(nSamples,nClusters,nOutcomes*(nOutcomes+1)/2))
      riskArray<-array(0,dim=c(nSamples,nClusters,nOutcomes))
    } else if(yModel=="LME"){
      nTimes=dim(longMat)[1]
      longMat<-as.data.frame(longMat)
      covREArray<-array(0,dim=c(nSamples,nClusters,nRandomEffects*(nRandomEffects+1)/2))
      SigmaLMEArray<-array(0,dim=c(nSamples))
      RE_LMEArray<-array(0,dim=c(nSamples, nSubjects,nRandomEffects))

      LME_timepoints <- seq(min(longMat$time),max(longMat$time),length.out=40)
      LMEArray <- array(0,dim=c(nSamples,nClusters,length(LME_timepoints)))
      riskArray<-array(0,dim=c(nSamples,nClusters))
    }
    if(nFixedEffects>0){
      betaArray<-array(0,dim=c(nSamples,nFixedEffects,nCategoriesY))
    }
    if(nFixedEffects_clust>0){
      betamixArray<-array(0,dim=c(nSamples,nClusters, nFixedEffects_clust*nCategoriesY))
    }
  }else{
    riskArray<-NULL
  }

  # Initialise the object for storing the profiles
  if(xModel=='Discrete'){
    phiArray<-array(dim=c(nSamples,nClusters,nCovariates,maxNCategories))
    if(varSelect){
      phiStarArray<-array(dim=c(nSamples,nClusters,nCovariates,maxNCategories))
      tmpCurrNullPhi<-scan(nullPhiFileName,what=double(),quiet=T)
      tmpCurrNullPhi<-array(tmpCurrNullPhi,dim=c(maxNCategories,nCovariates))
      currNullPhi<-array(dim=c(1,maxNCategories,nCovariates))
      currNullPhi[1,,]<-tmpCurrNullPhi
    }else{
      phiStarArray<-NULL
    }
  }else if(xModel=='Normal'){
    muArray<-array(dim=c(nSamples,nClusters,nCovariates))
    if(varSelect){
      muStarArray<-array(dim=c(nSamples,nClusters,nCovariates))
      currNullMu<-scan(nullMuFileName,what=double(),quiet=T)
      currNullMu<-array(currNullMu,dim=c(nCovariates,1))
      currNullMu<-t(currNullMu)
    }else{
      muStarArray<-NULL
    }
    sigmaArray<-array(dim=c(nSamples,nClusters,nCovariates,nCovariates))
  }else if(xModel=='Mixed'){
    phiArray<-array(dim=c(nSamples,nClusters,nDiscreteCovs,maxNCategories))
    if(varSelect){
      phiStarArray<-array(dim=c(nSamples,nClusters,nDiscreteCovs,maxNCategories))
      tmpCurrNullPhi<-scan(nullPhiFileName,what=double(),quiet=T)
      tmpCurrNullPhi<-array(tmpCurrNullPhi,dim=c(maxNCategories,nDiscreteCovs))
      currNullPhi<-array(dim=c(1,maxNCategories,nDiscreteCovs))
      currNullPhi[1,,]<-tmpCurrNullPhi
    }else{
      phiStarArray<-NULL
    }
    muArray<-array(dim=c(nSamples,nClusters,nContinuousCovs))
    if(varSelect){
      muStarArray<-array(dim=c(nSamples,nClusters,nContinuousCovs))
      currNullMu<-scan(nullMuFileName,what=double(),quiet=T)
      currNullMu<-array(currNullMu,dim=c(nContinuousCovs,1))
      currNullMu<-t(currNullMu)
    }else{
      muStarArray<-NULL
    }
    sigmaArray<-array(dim=c(nSamples,nClusters,nContinuousCovs,nContinuousCovs))
  }


  for(sweep in firstLine:lastLine){
    if(sweep==firstLine){
      skipVal<-firstLine-1
    }else{
      skipVal<-0
    }

    if(sweep-firstLine==0||(sweep-firstLine+1)%%1000==0){
      cat(paste("Processing sweep",sweep-firstLine+1,"of ",lastLine-firstLine+1,"\n"))
    }
    currMaxNClusters<-scan(nClustersFile,what=integer(),skip=skipVal,n=1,quiet=T)

    # Get the current allocation data for this sweep
    currZ<-scan(zFile,what=integer(),skip=skipVal,n=nSubjects+nPredictSubjects,quiet=T)
    currZ<-1+currZ

    if(includeResponse){
      # Get the risk data corresponding to this sweep
      if (yModel=="Categorical") {
        currThetaVector<-scan(thetaFile,what=double(),skip=skipVal,n=currMaxNClusters*(nCategoriesY-1),quiet=T)
        currTheta<-matrix(currThetaVector,ncol=(nCategoriesY-1),byrow=T)
        currTheta<-cbind(rep(0,dim(currTheta)[1]),currTheta)
      } else {
        if(!is.element(yModel, c("Longitudinal","LME" ))){
          currThetaVector<-scan(thetaFile,what=double(),skip=skipVal,n=currMaxNClusters*nCategoriesY,quiet=T)
          currTheta<-matrix(currThetaVector,ncol=nCategoriesY,byrow=T)
        }
        if (yModel=="Survival"&&!weibullFixedShape){
          currNuVector<-scan(nuFile,what=double(),skip=skipVal,n=currMaxNClusters,quiet=T)
          currNu<-c(currNuVector)
        } else if (yModel=="Longitudinal"){
          currLVector<-scan(LFile,what=double(),skip=skipVal,n=3*currMaxNClusters,quiet=T)
          currL<-matrix(currLVector,nrow=currMaxNClusters,ncol=3,byrow=T)
          if(sampleGPmean){
            GPmeanVector<-scan(GPmeanFile,what=double(),skip=skipVal,n=nTimes_unique*currMaxNClusters,quiet=T)
            currGPmean<-matrix(GPmeanVector,nrow=currMaxNClusters,ncol=nTimes_unique,byrow=T)
          }
        } else if (yModel=="MVN"){
          currMVNmuVector<-scan(MVNmuFile,what=double(),skip=skipVal,n=nOutcomes*currMaxNClusters,quiet=T)
          currMVNmu<-matrix(currMVNmuVector,nrow=currMaxNClusters,ncol=nOutcomes,byrow=T)
          currMVNSigmaVector<-scan(MVNSigmaFile,what=double(),skip=skipVal,n=(nOutcomes+1)*nOutcomes/2*currMaxNClusters,quiet=T)
          currMVNSigma<-matrix(currMVNSigmaVector,nrow=currMaxNClusters,ncol=(nOutcomes+1)*nOutcomes/2,byrow=T)

        } else if (yModel=="LME"){
          currcovREVector<-scan(covREFile,what=double(),skip=skipVal,n=(nRandomEffects+1)*nRandomEffects/2,quiet=T)
          currcovRE<-matrix(currcovREVector,nrow=1,ncol=(nRandomEffects+1)*nRandomEffects/2,byrow=T)
          currSigmaLME<-scan(SigmaLMEFile,what=double(),skip=skipVal,n=1,quiet=T)
          #currSigmaLME<-matrix(currSigmaLMEVector,nrow=currMaxNClusters,ncol=1,byrow=T)
          currRE_LMEVector<-scan(RE_LMEFile,what=double(),skip=skipVal,n=nSubjects*nRandomEffects,quiet=T)
          currRE_LME<-matrix(currRE_LMEVector,ncol=nRandomEffects,byrow=T)
        }
      }
      if(nFixedEffects>0){
        if (yModel=="Categorical") {
          currBetaVector<-scan(betaFile,what=double(),skip=skipVal,n=nFixedEffects*(nCategoriesY-1),quiet=T)
          currBeta<-matrix(currBetaVector,ncol=(nCategoriesY-1),byrow=T)
          currBeta<-cbind(rep(0,dim(currBeta)[1]),currBeta)
        } else {
          currBetaVector<-scan(betaFile,what=double(),skip=skipVal,n=nFixedEffects*nCategoriesY,quiet=T)
          currBeta<-matrix(currBetaVector,ncol=nCategoriesY,byrow=T)
        }
        betaArray[sweep-firstLine+1,,]<-currBeta
      }

      if(nFixedEffects_clust>0){
        currBetamixVector<-scan(betamixFile,what=double(),skip=skipVal,n=nFixedEffects_clust*nCategoriesY*currMaxNClusters,quiet=T)
        currBetamix<-matrix(currBetamixVector,nrow=currMaxNClusters,ncol=nFixedEffects_clust*nCategoriesY,byrow=T)
      }

      # Calculate the average risk (over subjects) for each cluster

      for(c in 1:nClusters){
        if(!is.element(yModel,c("Longitudinal","LME"))){
          currLambdaVector<-currTheta[currZ[optAlloc[[c]]],]
          currLambda<-matrix(currLambdaVector,ncol=nCategoriesY)
        }
        if(includeFixedEffects&&nFixedEffects>0){
          if (yModel=="Categorical"){
            for (i in 1:length(optAlloc[[c]])){
              for (k in 1:nCategoriesY) currLambda[i,k]<-currLambda[i,k]+
                  t(as.matrix(wMat[optAlloc[[c]][i],]))%*%
                  currBeta[,yMat[optAlloc[[c]][i]]+1]
            }
          } else {
            currLambda<-currLambda+as.matrix(wMat[optAlloc[[c]],])%*%currBeta
          }
        }
        if(yModel=="Poisson"){
          currRisk<-exp(currLambda )
        }else if(yModel=="Bernoulli"||yModel=="Binomial"){
          currRisk<-1.0/(1.0+exp(-currLambda))
        }else if(yModel=="Normal"){
          currRisk<-currLambda
        }else if(yModel=="Categorical"){
          currRisk<-matrix(0,ncol=length(optAlloc[[c]]),nrow=nCategoriesY)
          currRisk<-exp(currLambda)/rowSums(exp(currLambda))
        }else if(yModel=="Survival"){
          if (proportionalHazards){
            currRisk<-exp(currLambda)
          }else{
            if (!weibullFixedShape){
              currNuVector<-currNu[currZ[optAlloc[[c]]]]
            } else {
              currNuVector<-nu[sweep]
            }
            currRisk<-1/((exp(currLambda))^(1/currNuVector))*gamma(1+1/currNuVector)
          }
        }else if(yModel=="Longitudinal"){##//RJ current L
          LArray[sweep-firstLine+1,c,]<-apply(matrix(currL[currZ[optAlloc[[c]]],],ncol=3),2,mean)
          IDs <- unique(longMat$ID)
          currIDs <- IDs[optAlloc[[c]]]
          currRisk<-matrix(longMat$outcome[longMat$ID%in%currIDs] + longMean)
          if(sampleGPmean)
            GPmeanArray[sweep-firstLine+1,c,]<-apply(matrix(currGPmean[currZ[optAlloc[[c]]],],ncol=nTimes_unique),2,mean)
        }else if(yModel=="MVN"){##//RJ current MVN
          MVNmuArray[sweep-firstLine+1,c,]<-apply(matrix(currMVNmu[currZ[optAlloc[[c]]],],ncol=nOutcomes),2,mean)
          MVNSigmaArray[sweep-firstLine+1,c,]<-apply(matrix(currMVNSigma[currZ[optAlloc[[c]]],],ncol=(nOutcomes+1)*nOutcomes/2),2,mean)
          currRisk<-matrix(yMat[optAlloc[[c]],],nrow=length(optAlloc[[c]]))
          if(!all(!is.na(MVNmuArray[sweep-firstLine+1,c,])))
            browser()
        }else if(yModel=="LME"){
          #covREArray[sweep-firstLine+1,c,]<-apply(matrix(currcovRE[currZ[optAlloc[[c]]],],ncol=(nRandomEffects+1)*nRandomEffects/2),2,mean)
          covREArray[sweep-firstLine+1,c,]<-currcovRE
          if(c==1){
            SigmaLMEArray[sweep-firstLine+1]<-currSigmaLME
            RE_LMEArray[sweep-firstLine+1,, ]<-currRE_LME
          }
          #ind_time <- which()
          #LMEArray[sweep-firstLine+1,c,] <- apply(matrix(currBetamix[currZ[optAlloc[[c]]],],ncol=nFixedEffects_clust),2,mean)
          IDs <- unique(longMat$ID)
          currIDs <- IDs[optAlloc[[c]]]
          currRisk<-matrix(longMat$outcome[longMat$ID%in%currIDs] + longMean)
        }

        if(nFixedEffects_clust>0)
          betamixArray[sweep-firstLine+1,c,] <- apply(matrix(currBetamix[currZ[optAlloc[[c]]],],ncol=nFixedEffects_clust),2,mean)
        if(!is.element(yModel,c("Longitudinal","LME"))){
          thetaArray[sweep-firstLine+1,c,]<-apply(as.matrix(currTheta[currZ[optAlloc[[c]]],],ncol=nCategoriesY),2,mean)
          riskArray[sweep-firstLine+1,c,]<-apply(currRisk,2,mean)
        }else{
          riskArray[sweep-firstLine+1,c]<-mean(currRisk)
        }
        if(yModel=="Survival"&&!weibullFixedShape){
          nuArray[sweep-firstLine+1,c]<-mean(currNuVector)
        }
      }
    }

    # Calculate the average profile (over subjects) for each cluster
    if(xModel=='Discrete'){
      currPhi<-scan(phiFile,what=double(),
                    skip=skipVal,n=currMaxNClusters*maxNCategories*nCovariates,quiet=T)
      # This is slightly convoluted, because of the way that R reads in by column
      # I switched the order of categories and covariates in column below, and then
      # take the transpose to correct in the loop
      currPhi<-array(currPhi,dim=c(currMaxNClusters,maxNCategories,nCovariates))
      if(varSelect){
        # We increase dimensions of currGamma using duplicates, to
        # enable easier calculation of phiStar
        if(varSelectType=='BinaryCluster'){
          tmpCurrGamma<-scan(gammaFile,what=double(),
                             skip=skipVal,n=nCovariates*currMaxNClusters,quiet=T)
          tmpCurrGamma<-array(tmpCurrGamma,dim=c(currMaxNClusters,nCovariates))
        }else{
          tmpCurrGamma<-scan(gammaFile,what=double(),skip=skipVal,n=nCovariates,quiet=T)
          tmpCurrGamma<-array(tmpCurrGamma,dim=c(nCovariates,currMaxNClusters))
          tmpCurrGamma<-t(tmpCurrGamma)
        }
        currGamma<-array(dim=c(currMaxNClusters,maxNCategories,nCovariates))
        for(p in 1:maxNCategories){
          currGamma[,p,]<-tmpCurrGamma
        }
      }

      for(c in 1:nClusters){
        phiArray[sweep-firstLine+1,c,,]<-t(apply(array(currPhi[currZ[optAlloc[[c]]],,],
                                                       dim=c(length(optAlloc[[c]]),
                                                       dim(currPhi)[2],dim(currPhi)[3])),2:3,mean))
        if(varSelect){
          phiStarArray[sweep-firstLine+1,c,,]<-t(apply(array(currGamma[currZ[optAlloc[[c]]],,],
                                                             dim=c(length(optAlloc[[c]]),dim(currGamma)[2],
                                                                   dim(currGamma)[3]))*array(currPhi[currZ[optAlloc[[c]]],,],
                                                                                             dim=c(length(optAlloc[[c]]),dim(currPhi)[2],dim(currPhi)[3]))+
                                                         (1-array(currGamma[currZ[optAlloc[[c]]],,],
                                                                  dim=c(length(optAlloc[[c]]),dim(currGamma)[2],dim(currGamma)[3])))*
                                                         array(currNullPhi[rep(1,length(optAlloc[[c]])),,],
                                                               dim=c(length(optAlloc[[c]]),dim(currNullPhi)[2],
                                                                     dim(currNullPhi)[3])),2:3,mean))

        }
      }
    }else if(xModel=='Normal'){
      # mu stored like phi
      currMu<-scan(muFile,what=double(),skip=skipVal,n=currMaxNClusters*nCovariates,quiet=T)
      currMu<-array(currMu,dim=c(currMaxNClusters,nCovariates))
      if(varSelect){
        # We increase dimensions of nullPhi and currGamma using duplicates, to
        # enable easier calculation of phiStar
        # We increase dimensions of currGamma using duplicates, to
        # enable easier calculation of phiStar
        if(varSelectType=='BinaryCluster'){
          tmpCurrGamma<-scan(gammaFile,what=double(),skip=skipVal,n=currMaxNClusters*nCovariates,quiet=T)
          currGamma<-array(tmpCurrGamma,dim=c(currMaxNClusters,nCovariates))
        }else{
          tmpCurrGamma<-scan(gammaFile,what=double(),skip=skipVal,n=nCovariates,quiet=T)
          tmpCurrGamma<-array(tmpCurrGamma,dim=c(nCovariates,currMaxNClusters))
          currGamma<-t(tmpCurrGamma)
        }

      }
      for(c in 1:nClusters){
        muArray[sweep-firstLine+1,c,]<-apply(matrix(currMu[currZ[optAlloc[[c]]],],ncol=nCovariates),2,mean)
        if(varSelect){
          muStarArray[sweep-firstLine+1,c,]<-apply(matrix(currGamma[currZ[optAlloc[[c]]],],
                                                          ncol=nCovariates)*matrix(currMu[currZ[optAlloc[[c]]],],ncol=nCovariates)+
                                                     matrix(1-currGamma[currZ[optAlloc[[c]]],],ncol=nCovariates)*
                                                     matrix(currNullMu[rep(1,length(optAlloc[[c]])),],ncol=nCovariates),2,mean)
        }
      }

      currSigma<-scan(SigmaFile,what=double(),skip=skipVal,n=currMaxNClusters*nCovariates*nCovariates,quiet=T)
      currSigma<-array(currSigma,dim=c(currMaxNClusters,nCovariates,nCovariates))
      for(c in 1:nClusters){
        sigmaArray[sweep-firstLine+1,c,,]<-apply(array(currSigma[currZ[optAlloc[[c]]],,],
                                                       dim=c(length(optAlloc[[c]]),dim(currSigma)[2],dim(currSigma)[3])),2:3,mean)
      }
    }else if(xModel=='Mixed'){
      currPhi<-scan(phiFile,what=double(),
                    skip=skipVal,n=currMaxNClusters*maxNCategories*nDiscreteCovs,quiet=T)
      # This is slightly convoluted, because of the way that R reads in by column
      # I switched the order of categories and covariates in column below, and then
      # take the transpose to correct in the loop
      currPhi<-array(currPhi,dim=c(currMaxNClusters,maxNCategories,nDiscreteCovs))
      # mu stored like phi
      currMu<-scan(muFile,what=double(),skip=skipVal,n=currMaxNClusters*nContinuousCovs,quiet=T)
      currMu<-array(currMu,dim=c(currMaxNClusters,nContinuousCovs))
      if(varSelect){
        # We increase dimensions of currGamma using duplicates, to
        # enable easier calculation of phiStar
        if(varSelectType=='BinaryCluster'){
          tmpCurrGamma<-scan(gammaFile,what=double(),
                             skip=skipVal,n=nCovariates*currMaxNClusters,quiet=T)
          tmpCurrGamma<-array(tmpCurrGamma,dim=c(currMaxNClusters,nCovariates))
        }else{
          tmpCurrGamma<-scan(gammaFile,what=double(),skip=skipVal,n=nCovariates,quiet=T)
          tmpCurrGamma<-array(tmpCurrGamma,dim=c(nCovariates,currMaxNClusters))
          tmpCurrGamma<-t(tmpCurrGamma)
        }
        currGamma<-array(dim=c(currMaxNClusters,maxNCategories,nCovariates))
        for(p in 1:maxNCategories){
          currGamma[,p,]<-tmpCurrGamma
        }
      }
      for(c in 1:nClusters){
        phiArray[sweep-firstLine+1,c,,]<-t(apply(array(currPhi[currZ[optAlloc[[c]]],,],
                                                       dim=c(length(optAlloc[[c]]),
                                                             dim(currPhi)[2],dim(currPhi)[3])),2:3,mean))
        if(varSelect){
          currGammaD<-array(currGamma[,,1:nDiscreteCovs],dim=c(currMaxNClusters,maxNCategories,nDiscreteCovs))
          phiStarArray[sweep-firstLine+1,c,,]<-t(apply(array(currGammaD[currZ[optAlloc[[c]]],,],
                                                             dim=c(length(optAlloc[[c]]),dim(currGammaD)[2],
                                                                   dim(currGammaD)[3]))*array(currPhi[currZ[optAlloc[[c]]],,],
                                                                                              dim=c(length(optAlloc[[c]]),dim(currPhi)[2],dim(currPhi)[3]))+
                                                         (1-array(currGammaD[currZ[optAlloc[[c]]],,],
                                                                  dim=c(length(optAlloc[[c]]),dim(currGammaD)[2],dim(currGammaD)[3])))*
                                                         array(currNullPhi[rep(1,length(optAlloc[[c]])),,],
                                                               dim=c(length(optAlloc[[c]]),dim(currNullPhi)[2],
                                                                     dim(currNullPhi)[3])),2:3,mean))
        }
        muArray[sweep-firstLine+1,c,]<-apply(matrix(currMu[currZ[optAlloc[[c]]],],ncol=nContinuousCovs),2,mean)
        if(varSelect){
          currGammaC<-array(currGamma[,1,(nDiscreteCovs+1):nCovariates],
                            dim=c(currMaxNClusters,nContinuousCovs))
          muStarArray[sweep-firstLine+1,c,]<-apply(matrix(currGammaC[currZ[optAlloc[[c]]],],
                                                          ncol=nContinuousCovs)*matrix(currMu[currZ[optAlloc[[c]]],],ncol=nContinuousCovs)+
                                                     matrix(1-currGammaC[currZ[optAlloc[[c]]],],ncol=nContinuousCovs)*
                                                     matrix(currNullMu[rep(1,length(optAlloc[[c]])),],ncol=nContinuousCovs),2,mean)
        }
      }

      currSigma<-scan(SigmaFile,what=double(),skip=skipVal,n=currMaxNClusters*nContinuousCovs*nContinuousCovs,quiet=T)
      currSigma<-array(currSigma,dim=c(currMaxNClusters,nContinuousCovs,nContinuousCovs))
      for(c in 1:nClusters){
        sigmaArray[sweep-firstLine+1,c,,]<-apply(array(currSigma[currZ[optAlloc[[c]]],,],
                                                       dim=c(length(optAlloc[[c]]),dim(currSigma)[2],dim(currSigma)[3])),2:3,mean)
      }
    }
  }

  # Calculate the empiricals
  empiricals<-rep(0,nClusters)
  if(!is.null(yModel)){
    for(c in 1:nClusters){
      if(yModel=='Bernoulli'||yModel=='Normal'||yModel=='Survival'){
        empiricals[c]<-mean(yMat[optAlloc[[c]],1])
      }else if(yModel=='Binomial'){
        empiricals[c]<-mean(yMat[optAlloc[[c]],1]/yMat[optAlloc[[c]],2])
      }else if(yModel=='Poisson'){
        empiricals[c]<-mean(yMat[optAlloc[[c]],1]/yMat[optAlloc[[c]],2])
        #}else if(yModel=='Categorical'){
        # no empiricals for categorical outcome
      }else if(yModel=="Longitudinal"||yModel=="LME"){##//RJ same as empirical
        empiricals[c] <- mean(riskArray[,c])
      }else if(yModel=="MVN"){
        empiricals[c] <- mean(riskArray[,c,1])
      }
    }
  }

  if(xModel=='Discrete'){
    out<-list('riskProfClusObj'=clusObj,'risk'=riskArray,'profile'=phiArray,'profileStar'=phiStarArray,'empiricals'=empiricals)
  }else if(xModel=='Normal'){
    out<-list('riskProfClusObj'=clusObj,'risk'=riskArray,
              'profile'=muArray,'profileStar'=muStarArray,
              'profileStdDev'=sigmaArray,'empiricals'=empiricals)
  }else if(xModel=='Mixed'){
    out<-list('riskProfClusObj'=clusObj,'risk'=riskArray,
              'profilePhi'=phiArray,'profileStarPhi'=phiStarArray,
              'profileMu'=muArray,'profileStarMu'=muStarArray,
              'profileStdDev'=sigmaArray,'empiricals'=empiricals)
  }

  close(zFile)
  close(nClustersFile)
  if(xModel=="Discrete"){
    close(phiFile)
    if(varSelect){
      close(gammaFile)
    }
  }else if(xModel=="Normal"){
    close(muFile)
    close(SigmaFile)
    if(varSelect){
      close(gammaFile)
    }
  }else if(xModel=="Mixed"){
    close(phiFile)
    close(muFile)
    close(SigmaFile)
    if(varSelect){
      close(gammaFile)
    }
  }


  if(includeResponse){
    if(!is.element(yModel, c("Longitudinal","LME" )))
      close(thetaFile)
    if (yModel=="Survival"&&!weibullFixedShape) {
      out$nuArray<-nuArray
      close(nuFile)
    }
    if (yModel=="Longitudinal") {
      out$LArray<-LArray
      close(LFile)
      if(sampleGPmean){
        out$GPmeanArray<-GPmeanArray
        close(GPmeanFile)
      }
    }
    if (yModel=="MVN") {
      out$MVNmuArray<-MVNmuArray
      close(MVNmuFile)
      out$MVNSigmaArray<-MVNSigmaArray
      close(MVNSigmaFile)
    }
    if (yModel=="LME") {
      out$covREArray<-covREArray
      close(covREFile)
      out$SigmaLMEArray<-SigmaLMEArray
      close(SigmaLMEFile)
      out$RE_LMEArray<-RE_LMEArray
      close(RE_LMEFile)
      out$LMEArray <-LMEArray
    }
    if(nFixedEffects>0){
      out$betaArray <- betaArray
      close(betaFile)
    }
    if(nFixedEffects_clust>0){
      out$betamixArray <- betamixArray
      close(betamixFile)
    }
  }
  return(out)
}
