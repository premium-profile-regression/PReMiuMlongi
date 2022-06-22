plotRiskProfile_longi<-function(riskProfObj,outFile,showRelativeRisk=F,orderBy=NULL,whichClusters=NULL,whichCovariates=NULL,
                             useProfileStar=F,riskLim=NULL,bycol=FALSE, profile_X=NULL, timevar=NULL, double_plot = FALSE){

  riskProfClusObj=NULL
  clusObjRunInfoObj=NULL
  includeResponse=NULL
  yModel=NULL
  profileStar=NULL
  xModel=NULL
  whicCov=NULL
  nCategoriesY=NULL
  cluster=NULL
  prob=NULL
  meanProb=NULL
  fillColor=NULL
  lowerProb=NULL
  upperProb=NULL
  meanRisk=NULL
  lowerRisk=NULL
  upperRisk=NULL
  clusterSize=NULL
  mu=NULL
  meanMu=NULL
  lowerMu=NULL
  upperMu=NULL
  sigma=NULL
  meanSigma=NULL
  lowerSigma=NULL
  upperSigma=NULL
  weibullFixedShape=NULL
  nu=NULL
  meanNu=NULL
  lowerNu=NULL
  upperNu=NULL
  sampleGPmean=FALSE
  nFixedEffects_clust=0

  for (i in 1:length(riskProfObj)) assign(names(riskProfObj)[i],riskProfObj[[i]])
  for (i in 1:length(riskProfClusObj)) assign(names(riskProfClusObj)[i],riskProfClusObj[[i]])
  for (i in 1:length(clusObjRunInfoObj)) assign(names(clusObjRunInfoObj)[i],clusObjRunInfoObj[[i]])

  #assign(names(clusObjRunInfoObj)[4],min(nSweeps,clusObjRunInfoObj[[4]]))

  if(!is.null(profile_X)){
    if(length(which(!names(profile_X)%in%c(fixedEffectsNames,fixedEffectsNames_clust)))>0)
      stop("Error: Names of profile_X do not correspond to fixedEffectsNames nor fixedEffectsNames_clust")
  }

  if(is.null(profile_X) && yModel == "LME" && length(which(!c(fixedEffectsNames,fixedEffectsNames_clust)%in%c("intercept",timevar)))>0){
    stop("Error: profile_X should be defined as a list of covariates values such list(cov=1)")
  }

  if (nClusters==1) stop("Cannot produce plots because only one cluster has been found.")

  plotRiskFlag <- 2
  if(includeResponse){
    if(yModel=="Normal"){
      showRelativeRisk<-F
    }
    if(yModel=="Longitudinal"||yModel=="MVN"||yModel=="LME"){
      plotRiskFlag <- 1
    }
  }
  if(useProfileStar){
    profile<-profileStar
  }
  if(!is.null(whichCovariates)){
    if (!is.numeric(whichCovariates)){
      whichCovariatesTmp<-vector()
      for (k in 1:length(whichCovariates)){
        whichCovariatesTmp[k]<-which(riskProfClusObj$clusObjRunInfoObj$covNames==whichCovariates[k])
      }
      whichCovariates<-whichCovariatesTmp
    }
    if(xModel=='Discrete'){
      profile<-profile[,,whichCovariates,]
      nCategories<-nCategories[whichCovariates]
      covNames<-covNames[whichCovariates]
      nCovariates<-length(whichCovariates)
    }else if(xModel=='Normal'){
      profile<-profile[,,whichCovariates]
      profileStdDev<-profileStdDev[,,whichCovariates,whichCovariates]
      covNames<-covNames[whichCovariates]
      nCovariates<-length(whichCovariates)
    }else if(xModel=='Mixed'){
      nDiscreteCovsAll <- nDiscreteCovs
      nContinuousCovsAll <- nContinuousCovs
      whichDiscreteCovs <- whichCovariates[which(whichCovariates<=nDiscreteCovs)]
      whichContinuousCovs <- whichCovariates[which(whichCovariates>nDiscreteCovs)]
      discreteCovs <- discreteCovs[whichDiscreteCovs]
      nDiscreteCovs <- length(discreteCovs)
      continuousCovs <- continuousCovs[whichContinuousCovs-nDiscreteCovsAll]
      nContinuousCovs <- length(continuousCovs)
      profilePhi<-profilePhi[,,whichDiscreteCovs,, drop=FALSE]
      nCategories<-nCategories[whichDiscreteCovs]
      profileMu<-profileMu[,,whichContinuousCovs-nDiscreteCovsAll, drop=FALSE]
      profileStdDev<-profileStdDev[,,whichContinuousCovs-nDiscreteCovsAll,whichContinuousCovs-nDiscreteCovsAll]
      covNames<-c(discreteCovs,continuousCovs)
      nCovariates<-length(covNames)

    }
  }

  png(outFile,width=1200,height=800)
  orderProvided<-F

  if(!is.null(orderBy)){
    if(!includeResponse){
      if(orderBy!='Empirical'&&orderBy!='ClusterSize'&&!orderBy%in%covNames){
        if(is.numeric(orderBy)){
          if(length(orderBy)==nClusters){
            orderProvided<-T
            meanSortIndex<-orderBy
          }else{
            cat("Order vector provided not of same length as number of clusters. Reverting to default ordering.\n")
            orderBy<-NULL
          }
          orderBy<-NULL
        }
        #orderBy<-NULL
      }
    }else{
      if(orderBy!='Risk'&&orderBy!='Empirical'&&orderBy!='ClusterSize'&&!orderBy%in%covNames){
        if(is.numeric(orderBy)){
          if(length(orderBy)==nClusters){
            orderProvided<-T
            meanSortIndex<-orderBy
          }else{
            cat("Order vector provided not of same length as number of clusters. Reverting to default ordering.\n")
            orderBy<-NULL
          }
          orderBy<-NULL
        }
        #orderBy<-NULL
      }
    }

  }

  # Set up the layout for the plot
  ##//RJ remove 'Risk' from Longitudinal summary
  plotLayout<-grid.layout(ncol = nCovariates+plotRiskFlag, nrow = 6)
  grid.newpage()
  pushViewport(viewport(layout = plotLayout))
  if(!orderProvided){
    if(!is.null(risk)){
      if(is.null(orderBy)){
        # Default is to order by posterior theta risk
        # Compute the means
        orderStat<-apply(risk,2,median)
      }else{
        if(orderBy=='Risk'){
          orderStat<-apply(risk,2,median)
        }else if(orderBy=='Empirical'){
          orderStat<-empiricals
        }else if(orderBy=='ClusterSize'){
          orderStat<-clusterSizes
        }else{
          whichCov<-match(orderBy,covNames)
          if(xModel=='Normal'){
            orderStat<-apply(profile[,,whichCov],2,median)
          }else{
            # This assumes that there is some order to the categories
            # and then uses an expected value
            tmpMat<-profile[,,whichCov,1]
            if(nCategories[whichCov]>1){
              for(k in 2:nCategories[whichCov]){
                tmpMat<-tmpMat+k*profile[,,whichCov,k]
              }
            }
            orderStat<-apply(tmpMat,2,median)
          }
        }
      }
    }else{
      if(is.null(orderBy)){
        # Default is to order by empirical risk
        orderStat<-empiricals
      }else{
        if(orderBy=='Empirical'){
          orderStat<-empiricals
        }else if(orderBy=='ClusterSize'){
          orderStat<-clusterSizes
        }else{
          whichCov<-match(orderBy,covNames)
          if(xModel=='Normal'){
            orderStat<-apply(profile[,,whichCov],2,median)
          }else{
            # This assumes that there is some order to the categories
            # and then uses an expected value
            tmpMat<-profile[,,whichCov,1]
            if(nCategories[whichCov]>1){
              for(k in 2:(nCategories[whichCov])){
                tmpMat<-tmpMat+k*profile[,,whichCov,k]
              }
            }
            orderStat<-apply(tmpMat,2,median)
          }
        }
      }
    }
    # Sort into ascending mean size
    meanSortIndex<-order(orderStat,decreasing=F)
    meanSortIndex<-seq(1:max(meanSortIndex))
  }
  if(includeResponse){
    # Reorder the risk matrix
    riskDim<-dim(risk)
    if (!is.element(yModel,c("Longitudinal","LME"))){
      risk<-array(risk[,meanSortIndex,],dim=riskDim)
      if(showRelativeRisk){
        for(c in nClusters:1){
          risk[,c,]<-risk[,c,]/risk[,1,]
        }
      }
    }else{
      risk<-array(risk[,meanSortIndex],dim=riskDim)
      if(showRelativeRisk){
        for(c in nClusters:1){
          risk[,c]<-risk[,c]/risk[,1]
        }
      }
    }

    # reorder the nu matrix
    if (yModel=="Survival"&&!weibullFixedShape){
      nuDim<-dim(nuArray)
      nuArray<-array(nuArray[,meanSortIndex],dim=nuDim)
    }
    ##//RJ reorder the L matrix
    if (yModel=="Longitudinal"){
      LDim<-dim(LArray)
      LArray<-array(LArray[,meanSortIndex,],dim=LDim)
      if(sampleGPmean){
        GPmeanDim<-dim(GPmeanArray)
        GPmeanArray<-array(GPmeanArray[,meanSortIndex,],dim=GPmeanDim)
      }
    }else if (yModel=="MVN"){
      MVNmuDim<-dim(MVNmuArray)
      MVNSigmaDim<-dim(MVNSigmaArray)
      MVNmuArray<-array(MVNmuArray[,meanSortIndex,],dim=MVNmuDim)
      MVNSigmaArray<-array(MVNSigmaArray[,meanSortIndex,],dim=MVNSigmaDim)
    }else if (yModel=="LME"){
      covREArray<-array(covREArray[,meanSortIndex,],dim=dim(covREArray))
      RE_LMEArray<-array(RE_LMEArray[,meanSortIndex,],dim=dim(RE_LMEArray))
      SigmaLMEArray<-array(RE_LMEArray[,meanSortIndex,],dim=dim(RE_LMEArray))
    }
    if(nFixedEffects_clust>0){
      betamixArray<-array(betamixArray[,meanSortIndex,],dim=dim(betamixArray))
    }
  }

  # Reorder the cluster sizes
  clusterSizes<-clusterSizes[meanSortIndex]
  # Reorder the empiricals
  empiricals<-empiricals[meanSortIndex]
  meanEmpirical<-sum(empiricals*clusterSizes)/sum(clusterSizes)
  if(includeResponse){
    # Recompute the means and now also credible intervals
    riskMeans<-apply(risk,2,mean,trim=0.005)
    riskMean<-sum(riskMeans*clusterSizes)/sum(clusterSizes)
    riskLower<-apply(risk,2,quantile,0.05)
    riskUpper<-apply(risk,2,quantile,0.95)
    # The next line is to avoid outliers spoiling plot scales
    plotMax<-max(riskUpper)

    # Get the plot colors
    riskColor<-ifelse(riskLower>rep(riskMean,nClusters),"high",
                      ifelse(riskUpper<rep(riskMean,nClusters),"low","avg"))
    if (yModel=="Categorical"){
      riskDF<-data.frame("risk"=c(),"category"=c(),"cluster"=c(),"meanRisk"=c(),
                         "lowerRisk"=c(),"upperRisk"=c(),"fillColor"=c())
    } else {
      riskDF<-data.frame("risk"=c(),"cluster"=c(),"meanRisk"=c(),
                         "lowerRisk"=c(),"upperRisk"=c(),"fillColor"=c())
    }
    if (yModel=="Survival"&&!weibullFixedShape){
      # Recompute the means and now also credible intervals
      nuMeans<-apply(nuArray,2,mean,trim=0.005)
      nuMean<-sum(nuMeans*clusterSizes)/sum(clusterSizes)
      nuLower<-apply(nuArray,2,quantile,0.05)
      nuUpper<-apply(nuArray,2,quantile,0.95)

      nuDF<-data.frame("nu"=c(),"cluster"=c(),"meanNu"=c(),
                       "lowerNu"=c(),"upperNu"=c(),"fillColor"=c())
    }
    if (yModel=="Longitudinal"){
      # Recompute the means and now also credible intervals
      LMeans <- matrix(0,ncol=3,nrow=dim(LArray)[2])
      #LMean <- c()
      for(i in 1:3){
        LMeans[,i]<-apply(LArray[,,i],2,mean,trim=0.005)
        #LMean[i]<-sum(LMeans*clusterSizes)/sum(clusterSizes)
      }
      if(sampleGPmean){
        GPmeanMeans <- matrix(0,ncol=nTimes_unique,nrow=dim(GPmeanArray)[2])
        GPmeanMeans_sup <- matrix(0,ncol=nTimes_unique,nrow=dim(GPmeanArray)[2])
        GPmeanMeans_inf <- matrix(0,ncol=nTimes_unique,nrow=dim(GPmeanArray)[2])
        #GPmeanMean <- c()

        for(i in 1:nTimes_unique){
          GPmeanMeans[,i]<-apply(GPmeanArray[,,i],2,mean,trim=0.005)
          #GPmeanMean[i]<-sum(GPmeanMeans*clusterSizes)/sum(clusterSizes)
          GPmeanMeans_sup[,i] <- apply(GPmeanArray[,,i],2,quantile,0.95)
          GPmeanMeans_inf[,i] <- apply(GPmeanArray[,,i],2,quantile,0.05)
        }
      }
    }
    if (yModel=="MVN"){
      # Recompute the means and now also credible intervals
      MVNmuMeans<-apply(MVNmuArray,2,mean,trim=0.005)
      MVNmuMean<-sum(MVNmuMeans*clusterSizes)/sum(clusterSizes)
      MVNmuLower<-apply(MVNmuArray,2,quantile,0.05)
      MVNmuUpper<-apply(MVNmuArray,2,quantile,0.95)
      MVNmuDF<-data.frame("MVNmu"=c(),"cluster"=c(),"meanMVNmu"=c(),
                          "lowerMVNmu"=c(),"upperMVNmu"=c(),"fillColor"=c())
    }
    if (yModel=="LME"){
      # LMeans <- matrix(0,ncol=3,nrow=dim(LArray)[2])
      #
      # covREArray<-array(0,dim=c(nSamples,nClusters,nRandomEffects*(nRandomEffects+1)/2))
      # SigmaLMEArray<-array(0,dim=c(nSamples))
      # RE_LMEArray<-array(0,dim=c(nSamples, nSubjects,nRandomEffects))
      #
      # LMEmuMean
      # LMEmuLower
      # LMEmuUpper
    }

  }else{
    riskColor<-ifelse(empiricals>rep(meanEmpirical,length(empiricals)),"high",
                      ifelse(empiricals<rep(meanEmpirical,nClusters),"low","avg"))
  }

  if(is.null(whichClusters)){
    whichClusters<-1:nClusters
  }
  nClusters<-length(whichClusters)

  empiricalDF<-data.frame("empiricals"=c(),"meanEmpirical"=c(),"cluster"=c(),"fillColor"=c())
  sizeDF<-data.frame("clusterSize"=c(),"cluster"=c(),"fillColor"=c())
  # Restructure the data for plotting
  for(c in whichClusters){
    if(includeResponse){
      if (yModel=="Categorical"){
        plotRisk<-risk[,c,]
        nPoints<-dim(plotRisk)[1]
        for (k in 1:nCategoriesY){
          riskDF<-rbind(riskDF,data.frame("risk"=plotRisk[,k],
                                          "category"=rep(k,nPoints),
                                          "cluster"=rep(c,nPoints),
                                          "meanRisk"=rep(riskMean,nPoints),
                                          "lowerRisk"=rep(riskLower[c],nPoints),
                                          "upperRisk"=rep(riskUpper[c],nPoints),
                                          "fillColor"=rep(riskColor[c],nPoints)))
        }
      } else if (yModel=="Longitudinal"||yModel=="MVN"){##//RJ same as empirical
        if(yModel=="MVN"){
          plotRisk<-mean(risk[,c,1])
        }else{
          plotRisk<-mean(risk[,c])
        }
        nPoints<-length(plotRisk)
        riskDF<-rbind(riskDF,data.frame("risk"=plotRisk,"cluster"=rep(c,nPoints),
                                        "meanRisk"=rep(riskMean,nPoints),
                                        "lowerRisk"=rep(riskLower[c],nPoints),
                                        "upperRisk"=rep(riskUpper[c],nPoints),
                                        "fillColor"=rep(riskColor[c],nPoints)))
      } else {
        if(yModel=="LME"){
          plotRisk<-risk[,c]
        }else{
          plotRisk<-risk[,c,]
        }
        plotRisk<-plotRisk[plotRisk<plotMax]
        nPoints<-length(plotRisk)
        riskDF<-rbind(riskDF,data.frame("risk"=plotRisk,"cluster"=rep(c,nPoints),
                                        "meanRisk"=rep(riskMean,nPoints),
                                        "lowerRisk"=rep(riskLower[c],nPoints),
                                        "upperRisk"=rep(riskUpper[c],nPoints),
                                        "fillColor"=rep(riskColor[c],nPoints)))
      }
      if (yModel=="Survival"&&!weibullFixedShape){
        plotNu<-nuArray[,c]
        nPoints<-length(plotNu)
        nuDF<-rbind(nuDF,data.frame("nu"=plotNu,"cluster"=rep(c,nPoints),
                                    "meanNu"=rep(nuMean,nPoints),
                                    "lowerNu"=rep(nuLower[c],nPoints),
                                    "upperNu"=rep(nuUpper[c],nPoints),
                                    "fillColor"=rep(riskColor[c],nPoints)))
      }
    }
    empiricalDF<-rbind(empiricalDF,
                       data.frame("empiricals"=empiricals[c],
                                  "meanEmpirical"=meanEmpirical,"cluster"=c,"fillColor"=riskColor[c]))
    sizeDF<-rbind(sizeDF,
                  data.frame("clusterSize"=clusterSizes[c],"cluster"=c,"fillColor"=riskColor[c]))
  }
  if(includeResponse){
    if(yModel=='Categorical'){
      riskDF<-data.frame("prob"=c(),"cluster"=c(),"category"=c(),"meanProb"=c(),
                         "lowerProb"=c(),"upperProb"=c(),"fillColor"=c())
      for(k in 1:nCategoriesY){
        probMat<-risk[,,k]
        nPoints<-nrow(probMat)
        probMeans<-apply(probMat,2,mean)
        probMean<-sum(probMeans*clusterSizes)/sum(clusterSizes)
        probLower<-apply(probMat,2,quantile,0.05)
        probUpper<-apply(probMat,2,quantile,0.95)

        # Get the plot colors
        probColor<-ifelse(probLower>rep(probMean,length(probLower)),"high",
                          ifelse(probUpper<rep(probMean,length(probUpper)),"low","avg"))

        for(c in whichClusters){
          riskDF<-rbind(riskDF,data.frame("prob"=probMat[,c],"cluster"=rep(c,nPoints),
                                          "category"=rep(k-1,nPoints),
                                          "meanProb"=rep(probMean,nPoints),
                                          "lowerProb"=rep(probLower[c],nPoints),
                                          "upperProb"=rep(probUpper[c],nPoints),
                                          "fillColor"=rep(probColor[c],nPoints)))
          rownames(riskDF)<-seq(1,nrow(riskDF),1)

        }
      }

      plotObj<-ggplot(riskDF)
      plotObj<-plotObj+facet_wrap(~category,ncol=1,as.table=F,scales="free_y")
      plotObj<-plotObj+geom_hline(aes(x=as.factor(cluster),y=prob,yintercept=meanProb))
      plotObj<-plotObj+geom_boxplot(aes(x=as.factor(cluster),y=prob,fill=as.factor(fillColor)),outlier.size=0.5)
      plotObj<-plotObj+geom_point(aes(x=as.factor(cluster),y=lowerProb,colour=as.factor(fillColor)),size=1.5)
      plotObj<-plotObj+geom_point(aes(x=as.factor(cluster),y=upperProb,colour=as.factor(fillColor)),size=1.5)
      plotObj<-plotObj+
        scale_fill_manual(values = c(high ="#CC0033",low ="#0066CC", avg ="#33CC66"))+
        scale_colour_manual(values = c(high ="#CC0033",low ="#0066CC", avg ="#33CC66"))+
        theme(legend.position="none")+labs(x="Cluster")+theme(axis.title.x=element_text(size=10))
      plotObj<-plotObj+labs(y="Probability")+theme(axis.title.y=element_text(size=10,angle=90))
      plotObj<-plotObj+labs(title=ifelse(showRelativeRisk,'Relative Risk','Risk'),plot.title=element_text(size=10))
      # Margin order is (top,right,bottom,left)
      plotObj<-plotObj+theme(plot.margin=unit(c(0.5,0.15,0.5,0.15),'lines'))+
        theme(plot.margin=unit(c(0,0,0,0),'lines'))
      print(plotObj,vp=viewport(layout.pos.row=1:6,layout.pos.col=2))
      	}else if (yModel=="Survival"&&!weibullFixedShape){
      		rownames(riskDF)<-seq(1,nrow(riskDF),by=1)

      		# Create the risk plot
      		plotObj<-ggplot(riskDF)
      		plotObj<-plotObj+geom_hline(aes(x=as.factor(cluster),y=risk,yintercept=meanRisk))
      		plotObj<-plotObj+geom_boxplot(aes(x=as.factor(cluster),y=risk,fill=as.factor(fillColor)),outlier.size=0.5)
      		plotObj<-plotObj+geom_point(aes(x=as.factor(cluster),y=lowerRisk,colour=as.factor(fillColor)),size=1.5)
      		plotObj<-plotObj+geom_point(aes(x=as.factor(cluster),y=upperRisk,colour=as.factor(fillColor)),size=1.5)
      		plotObj<-plotObj+scale_fill_manual(values = c(high ="#CC0033",low ="#0066CC", avg ="#33CC66"))+
      			scale_colour_manual(values = c(high ="#CC0033",low ="#0066CC", avg ="#33CC66"))+
      			theme(legend.position="none")+
      			labs(x="Cluster",y=ifelse(showRelativeRisk,'RR',
      			ifelse(yModel=="Categorical"||yModel=="Bernoulli"||yModel=="Binomial","Probability","E[Y]")))
      		plotObj<-plotObj+theme(axis.title.y=element_text(size=10,angle=90),axis.title.x=element_text(size=10))
      		plotObj<-plotObj+labs(title=ifelse(showRelativeRisk,'Relative Risk','Risk'),plot.title=element_text(size=10))
      		# Margin order is (top,right,bottom,left)
      		plotObj<-plotObj+theme(plot.margin=unit(c(0,0,0,0),'lines'))+theme(plot.margin=unit(c(0.5,0.15,0.5,0.15),'lines'))
      		print(plotObj,vp=viewport(layout.pos.row=1:3,layout.pos.col=2))

      		rownames(nuDF)<-seq(1,nrow(nuDF),by=1)
      		# Create the nu plot
      		plotObj<-ggplot(nuDF)
      		plotObj<-plotObj+geom_hline(aes(x=as.factor(cluster),y=nu,yintercept=meanNu))
      		plotObj<-plotObj+geom_boxplot(aes(x=as.factor(cluster),y=nu,fill=as.factor(fillColor)),outlier.size=0.5)
      		plotObj<-plotObj+geom_point(aes(x=as.factor(cluster),y=lowerNu,colour=as.factor(fillColor)),size=1.5)
      		plotObj<-plotObj+geom_point(aes(x=as.factor(cluster),y=upperNu,colour=as.factor(fillColor)),size=1.5)
      		plotObj<-plotObj+scale_fill_manual(values = c(high ="#CC0033",low ="#0066CC", avg ="#33CC66"))+
      			scale_colour_manual(values = c(high ="#CC0033",low ="#0066CC", avg ="#33CC66"))+
      			theme(legend.position="none")+
      			labs(x="Cluster",y="Shape Parameter")
      		plotObj<-plotObj+theme(axis.title.y=element_text(size=10,angle=90),axis.title.x=element_text(size=10))
      		plotObj<-plotObj+labs(title="",plot.title=element_text(size=10))
      		# Margin order is (top,right,bottom,left)
      		plotObj<-plotObj+theme(plot.margin=unit(c(0,0,0,0),'lines'))+theme(plot.margin=unit(c(0.5,0.15,0.5,0.15),'lines'))
      		print(plotObj,vp=viewport(layout.pos.row=4:6,layout.pos.col=2))
    } else if (yModel=="MVN"){ ##//RJ yModel!='Longitudinal'||
      rownames(riskDF)<-seq(1,nrow(riskDF),by=1)

      # Create the risk plot
      plotObj<-ggplot(riskDF)
      plotObj<-plotObj+geom_hline(aes(x=as.factor(cluster),y=risk,yintercept=meanRisk))
      plotObj<-plotObj+geom_boxplot(aes(x=as.factor(cluster),y=risk,fill=as.factor(fillColor)),outlier.size=0.5)
      if (!is.null(riskLim)) plotObj<-plotObj+coord_cartesian(ylim = riskLim)
      plotObj<-plotObj+geom_point(aes(x=as.factor(cluster),y=lowerRisk,colour=as.factor(fillColor)),size=1.5)
      plotObj<-plotObj+geom_point(aes(x=as.factor(cluster),y=upperRisk,colour=as.factor(fillColor)),size=1.5)
      plotObj<-plotObj+scale_fill_manual(values = c(high ="#CC0033",low ="#0066CC", avg ="#33CC66"))+
        scale_colour_manual(values = c(high ="#CC0033",low ="#0066CC", avg ="#33CC66"))+
        theme(legend.position="none")+
        labs(x="Cluster",y=ifelse(showRelativeRisk,'RR',
                                  ifelse(yModel=="Categorical"||yModel=="Bernoulli"||yModel=="Binomial","Probability","E[Y]")))
      plotObj<-plotObj+theme(axis.title.y=element_text(size=10,angle=90),axis.title.x=element_text(size=10))
      plotObj<-plotObj+labs(title=ifelse(showRelativeRisk,'Relative Risk','Risk'),plot.title=element_text(size=10))
      # Margin order is (top,right,bottom,left)
      plotObj<-plotObj+theme(plot.margin=unit(c(0,0,0,0),'lines'))+theme(plot.margin=unit(c(0.5,0.15,0.5,0.15),'lines'))
      print(plotObj,vp=viewport(layout.pos.row=1:6,layout.pos.col=2))
    }
  }

  # Create a bar chart of cluster empiricals
  if((!is.null(yModel))){
    if(yModel!="Categorical"){
      plotObj<-ggplot(empiricalDF)
      plotObj<-plotObj+geom_point(aes(x=as.factor(cluster),y=empiricals,colour=as.factor(fillColor)),size=3)
      plotObj<-plotObj+geom_hline(aes(x=as.factor(cluster),y=empiricals,yintercept=meanEmpirical))
      plotObj<-plotObj+scale_colour_manual(values = c(high ="#CC0033",low ="#0066CC", avg ="#33CC66"))+
        theme(legend.position="none")
      plotObj<-plotObj+labs(title='Empirical Data',plot.title=element_text(size=10))
      plotObj<-plotObj+theme(axis.title.x=element_text(size=10),axis.title.y=element_text(size=10,angle=90))
      plotObj<-plotObj+
        labs(y=ifelse(yModel=="Bernoulli","Proportion of cases",
                      ifelse(yModel=="Binomial","Avg Proportion of occurrence",
                             ifelse(yModel=="Poisson","Avg Count",
                                    ifelse(yModel=="Survival","Avg Survival Time",
                                           ifelse(yModel=="Categorical","Avg Proportion of occurrence","Avg Y"))))),x="Cluster")
      plotObj<-plotObj+theme(plot.margin=unit(c(0,0,0,0),'lines'))+theme(plot.margin=unit(c(0.15,0.5,0.5,1),'lines'))
      print(plotObj,vp=viewport(layout.pos.row=1:3,layout.pos.col=1))
    }
  }
  # Create a bar chart of cluster sizes
  plotObj<-ggplot(sizeDF)
  plotObj<-plotObj+geom_point(aes(x=as.factor(cluster),y=clusterSize,colour=as.factor(fillColor)),size=3)
  plotObj<-plotObj+scale_colour_manual(values = c(high ="#CC0033",low ="#0066CC", avg ="#33CC66"))+theme(legend.position="none")
  plotObj<-plotObj+labs(title="Size",plot.title=element_text(size=10))
  plotObj<-plotObj+theme(axis.title.x=element_text(size=10),axis.title.y=element_text(size=10,angle=90))
  plotObj<-plotObj+labs(y="No. of Subjects",x="Cluster")
  plotObj<-plotObj+theme(plot.margin=unit(c(0,0,0,0),'lines'))+theme(plot.margin=unit(c(0.15,0.5,0.5,1),'lines'))
  print(plotObj,vp=viewport(layout.pos.row=4:6,layout.pos.col=1))
  # Loop over the covariates
  for(j in 1:nCovariates){
    if(xModel=='Discrete'){
      profileDF<-data.frame("prob"=c(),"cluster"=c(),"category"=c(),"meanProb"=c(),
                            "lowerProb"=c(),"upperProb"=c(),"fillColor"=c())
      for(k in 1:nCategories[j]){
        probMat<-profile[,meanSortIndex,j,k]
        nPoints<-nrow(probMat)
        probMeans<-apply(probMat,2,mean)
        probMean<-sum(probMeans*clusterSizes)/sum(clusterSizes)
        probLower<-apply(probMat,2,quantile,0.05)
        probUpper<-apply(probMat,2,quantile,0.95)
        # Get the plot colors
        probColor<-ifelse(probLower>rep(probMean,length(probLower)),"high",
                          ifelse(probUpper<rep(probMean,length(probUpper)),"low","avg"))


        for(c in whichClusters){
          profileDF<-rbind(profileDF,data.frame("prob"=probMat[,c],"cluster"=rep(c,nPoints),
                                                "category"=rep(k-1,nPoints),
                                                "meanProb"=rep(probMean,nPoints),
                                                "lowerProb"=rep(probLower[c],nPoints),
                                                "upperProb"=rep(probUpper[c],nPoints),
                                                "fillColor"=rep(probColor[c],nPoints)))
          rownames(profileDF)<-seq(1,nrow(profileDF),1)

        }
      }
      plotObj<-ggplot(profileDF)
      plotObj<-plotObj+facet_wrap(~category,ncol=1,as.table=F,scales="free_y")
      plotObj<-plotObj+geom_hline(aes(x=as.factor(cluster),y=prob,yintercept=meanProb))
      plotObj<-plotObj+geom_boxplot(aes(x=as.factor(cluster),y=prob,fill=as.factor(fillColor)),outlier.size=0.5)
      plotObj<-plotObj+geom_point(aes(x=as.factor(cluster),y=lowerProb,colour=as.factor(fillColor)),size=1.5)
      plotObj<-plotObj+geom_point(aes(x=as.factor(cluster),y=upperProb,colour=as.factor(fillColor)),size=1.5)
      plotObj<-plotObj+
        scale_fill_manual(values = c(high ="#CC0033",low ="#0066CC", avg ="#33CC66"))+
        scale_colour_manual(values = c(high ="#CC0033",low ="#0066CC", avg ="#33CC66"))+
        theme(legend.position="none")+labs(x="Cluster")+theme(axis.title.x=element_text(size=10))
      if(j==1){
        plotObj<-plotObj+labs(y="Probability")+theme(axis.title.y=element_text(size=10,angle=90))
      }else{
        plotObj<-plotObj+theme(axis.title.y=element_blank())
      }
      plotObj<-plotObj+labs(title=covNames[j],plot.title=element_text(size=10))
      plotObj<-plotObj+theme(plot.margin=unit(c(0.5,ifelse(j==nCovariates,1,0),0.5,ifelse(j==1,0.5,0)),'lines'))+
        theme(plot.margin=unit(c(0,0,0,0),'lines'))
      print(plotObj,vp=viewport(layout.pos.row=1:6,layout.pos.col=j+plotRiskFlag))

    }else if(xModel=='Normal'){
      # Plot the means
      profileDF<-data.frame("mu"=c(),"cluster"=c(),"muMean"=c(),
                            "lowerMu"=c(),"upperMu"=c(),"fillColor"=c())
      muMat<-profile[,meanSortIndex,j]
      muMeans<-apply(muMat,2,mean)
      muMean<-sum(muMeans*clusterSizes)/sum(clusterSizes)
      muLower<-apply(muMat,2,quantile,0.05)
      muUpper<-apply(muMat,2,quantile,0.95)
      # The next line is to avoid outliers spoiling plot scales
      plotMax<-max(muUpper)
      plotMin<-min(muLower)

      # Get the plot colors
      muColor<-ifelse(muLower>rep(muMean,length(muLower)),"high",
                      ifelse(muUpper<rep(muMean,length(muUpper)),"low","avg"))
      muColor<-c("high","avg","low")
      for(c in whichClusters){
        plotMu<-muMat[,c]
        plotMu<-plotMu[plotMu<plotMax&plotMu>plotMin]
        nPoints<-length(plotMu)
        profileDF<-rbind(profileDF,data.frame("mu"=plotMu,"cluster"=rep(c,nPoints),
                                              "meanMu"=rep(muMean,nPoints),
                                              "lowerMu"=rep(muLower[c],nPoints),
                                              "upperMu"=rep(muUpper[c],nPoints),
                                              "fillColor"=rep(muColor[c],nPoints)))
      }
      rownames(profileDF)<-seq(1,nrow(profileDF),1)
      plotObj<-ggplot(profileDF)
      plotObj<-plotObj+geom_hline(aes(x=as.factor(cluster),y=mu,yintercept=meanMu))
      plotObj<-plotObj+geom_boxplot(aes(x=as.factor(cluster),y=mu,fill=as.factor(fillColor)),outlier.size=0.5)
      plotObj<-plotObj+geom_point(aes(x=as.factor(cluster),y=lowerMu,colour=as.factor(fillColor)),size=1.5)
      plotObj<-plotObj+geom_point(aes(x=as.factor(cluster),y=upperMu,colour=as.factor(fillColor)),size=1.5)
      plotObj<-plotObj+
        scale_fill_manual(values = c(high ="#CC0033",low ="#0066CC", avg ="#33CC66"))+
        scale_colour_manual(values = c(high ="#CC0033",low ="#0066CC", avg ="#33CC66"))+
        theme(legend.position="none")+labs(x="Cluster")+theme(axis.title.x=element_text(size=10))
      if(j==1){
        plotObj<-plotObj+labs(y="Mean")+theme(axis.title.y=element_text(size=10,angle=90))
      }else{
        plotObj<-plotObj+theme(axis.title.y=element_blank())
      }
      plotObj<-plotObj+labs(title=covNames[j],plot.title=element_text(size=10))
      plotObj<-plotObj+
        theme(plot.margin=unit(c(0.5,ifelse(j==nCovariates,1,0),0.5,ifelse(j==1,0.5,0)),'lines'))+
        theme(plot.margin=unit(c(0,0,0,0),'lines'))

      print(plotObj,vp=viewport(layout.pos.row=1:3,layout.pos.col=j+plotRiskFlag))
      # Plot the variances
      profileDF<-data.frame("sigma"=c(),"cluster"=c(),"sigmaMean"=c(),
                            "lowerSigma"=c(),"upperSigma"=c(),"fillColor"=c())
      sigmaMat<-profileStdDev[,meanSortIndex,j,j]
      sigmaMeans<-apply(sigmaMat,2,mean)
      sigmaMean<-sum(sigmaMeans*clusterSizes)/sum(clusterSizes)
      sigmaLower<-apply(sigmaMat,2,quantile,0.05)
      sigmaUpper<-apply(sigmaMat,2,quantile,0.95)
      # The next line is to avoid outliers spoiling plot scales
      plotMax<-max(sigmaUpper)

      # Get the plot colors
      sigmaColor<-ifelse(sigmaLower>rep(sigmaMean,length(sigmaLower)),"high",
                         ifelse(sigmaUpper<rep(sigmaMean,length(sigmaUpper)),"low","avg"))
      for(c in whichClusters){
        plotSigma<-sigmaMat[,c]
        plotSigma<-plotSigma[plotSigma<plotMax]
        nPoints<-length(plotSigma)
        profileDF<-rbind(profileDF,data.frame("sigma"=plotSigma,"cluster"=rep(c,nPoints),
                                              "meanSigma"=rep(sigmaMean,nPoints),
                                              "lowerSigma"=rep(sigmaLower[c],nPoints),
                                              "upperSigma"=rep(sigmaUpper[c],nPoints),
                                              "fillColor"=rep(sigmaColor[c],nPoints)))
      }
      rownames(profileDF)<-seq(1,nrow(profileDF),1)

      plotObj<-ggplot(profileDF)
      plotObj<-plotObj+geom_hline(aes(x=as.factor(cluster),y=sigma,yintercept=meanSigma))
      plotObj<-plotObj+geom_boxplot(aes(x=as.factor(cluster),y=sigma,fill=as.factor(fillColor)),outlier.size=0.5)
      plotObj<-plotObj+geom_point(aes(x=as.factor(cluster),y=lowerSigma,colour=as.factor(fillColor)),size=1.5)
      plotObj<-plotObj+geom_point(aes(x=as.factor(cluster),y=upperSigma,colour=as.factor(fillColor)),size=1.5)
      plotObj<-plotObj+
        scale_fill_manual(values = c(high ="#CC0033",low ="#0066CC", avg ="#33CC66"))+
        scale_colour_manual(values = c(high ="#CC0033",low ="#0066CC", avg ="#33CC66"))+
        theme(legend.position="none")+labs(x="Cluster")+theme(axis.title.x=element_text(size=10))
      if(j==1){
        plotObj<-plotObj+labs(y="Std Dev")+theme(axis.title.y=element_text(size=10,angle=90))
      }else{
        plotObj<-plotObj+theme(axis.title.y=element_blank())
      }
      plotObj<-plotObj+
        theme(plot.margin=unit(c(0.5,ifelse(j==nCovariates,1,0),0.5,ifelse(j==1,0.5,0)),'lines'))+
        theme(plot.margin=unit(c(0,0,0,0),'lines'))
      print(plotObj,vp=viewport(layout.pos.row=4:6,layout.pos.col=j+plotRiskFlag))
    }else if(xModel=='Mixed'){
      if (j<=nDiscreteCovs){
        profileDF<-data.frame("prob"=c(),"cluster"=c(),"category"=c(),"meanProb"=c(),
                              "lowerProb"=c(),"upperProb"=c(),"fillColor"=c())
        for(k in 1:nCategories[j]){
          if (nDiscreteCovs==1) {
            probMat<-profilePhi[,meanSortIndex,1,k]
          } else {
            probMat<-profilePhi[,meanSortIndex,j,k]
          }
          nPoints<-nrow(probMat)
          probMeans<-apply(probMat,2,mean)
          probMean<-sum(probMeans*clusterSizes)/sum(clusterSizes)
          probLower<-apply(probMat,2,quantile,0.05)
          probUpper<-apply(probMat,2,quantile,0.95)

          # Get the plot colors
          probColor<-ifelse(probLower>rep(probMean,length(probLower)),"high",
                            ifelse(probUpper<rep(probMean,length(probUpper)),"low","avg"))


          for(c in whichClusters){
            profileDF<-rbind(profileDF,data.frame("prob"=probMat[,c],"cluster"=rep(c,nPoints),
                                                  "category"=rep(k-1,nPoints),
                                                  "meanProb"=rep(probMean,nPoints),
                                                  "lowerProb"=rep(probLower[c],nPoints),
                                                  "upperProb"=rep(probUpper[c],nPoints),
                                                  "fillColor"=rep(probColor[c],nPoints)))
            rownames(profileDF)<-seq(1,nrow(profileDF),1)

          }
        }
        plotObj<-ggplot(profileDF)
        plotObj<-plotObj+facet_wrap(~category,ncol=1,as.table=F,scales="free_y")
        plotObj<-plotObj+geom_hline(aes(x=as.factor(cluster),y=prob,yintercept=meanProb))
        plotObj<-plotObj+geom_boxplot(aes(x=as.factor(cluster),y=prob,fill=as.factor(fillColor)),outlier.size=0.5)
        plotObj<-plotObj+geom_point(aes(x=as.factor(cluster),y=lowerProb,colour=as.factor(fillColor)),size=1.5)
        plotObj<-plotObj+geom_point(aes(x=as.factor(cluster),y=upperProb,colour=as.factor(fillColor)),size=1.5)
        plotObj<-plotObj+
          scale_fill_manual(values = c(high ="#CC0033",low ="#0066CC", avg ="#33CC66"))+
          scale_colour_manual(values = c(high ="#CC0033",low ="#0066CC", avg ="#33CC66"))+
          theme(legend.position="none")+labs(x="Cluster")+theme(axis.title.x=element_text(size=10))
        if(j==1){
          plotObj<-plotObj+labs(y="Probability")+theme(axis.title.y=element_text(size=10,angle=90))
        }else{
          plotObj<-plotObj+theme(axis.title.y=element_blank())
        }
        plotObj<-plotObj+labs(title=covNames[j],plot.title=element_text(size=10))
        plotObj<-plotObj+theme(plot.margin=unit(c(0.5,ifelse(j==nCovariates,1,0),0.5,ifelse(j==1,0.5,0)),'lines'))+
          theme(plot.margin=unit(c(0,0,0,0),'lines'))
        print(plotObj,vp=viewport(layout.pos.row=1:6,layout.pos.col=j+plotRiskFlag))

      } else {
        # Plot the means
        profileDF<-data.frame("mu"=c(),"cluster"=c(),"muMean"=c(),
                              "lowerMu"=c(),"upperMu"=c(),"fillColor"=c())
        if (nContinuousCovs==1){
          muMat<-profileMu[,meanSortIndex,1]
        } else {
          muMat<-profileMu[,meanSortIndex,(j-nDiscreteCovs)]
        }
        muMeans<-apply(muMat,2,mean)
        muMean<-sum(muMeans*clusterSizes)/sum(clusterSizes)
        muLower<-apply(muMat,2,quantile,0.05)
        muUpper<-apply(muMat,2,quantile,0.95)
        # The next line is to avoid outliers spoiling plot scales
        plotMax<-max(muUpper)
        plotMin<-min(muLower)

        # Get the plot colors
        muColor<-ifelse(muLower>rep(muMean,length(muLower)),"high",
                        ifelse(muUpper<rep(muMean,length(muUpper)),"low","avg"))
        for(c in whichClusters){
          plotMu<-muMat[,c]
          plotMu<-plotMu[plotMu<plotMax&plotMu>plotMin]
          nPoints<-length(plotMu)
          profileDF<-rbind(profileDF,data.frame("mu"=plotMu,"cluster"=rep(c,nPoints),
                                                "meanMu"=rep(muMean,nPoints),
                                                "lowerMu"=rep(muLower[c],nPoints),
                                                "upperMu"=rep(muUpper[c],nPoints),
                                                "fillColor"=rep(muColor[c],nPoints)))
        }
        rownames(profileDF)<-seq(1,nrow(profileDF),1)

        plotObj<-ggplot(profileDF)
        plotObj<-plotObj+geom_hline(aes(x=as.factor(cluster),y=mu,yintercept=meanMu))
        plotObj<-plotObj+geom_boxplot(aes(x=as.factor(cluster),y=mu,fill=as.factor(fillColor)),outlier.size=0.5)
        plotObj<-plotObj+geom_point(aes(x=as.factor(cluster),y=lowerMu,colour=as.factor(fillColor)),size=1.5)
        plotObj<-plotObj+geom_point(aes(x=as.factor(cluster),y=upperMu,colour=as.factor(fillColor)),size=1.5)
        plotObj<-plotObj+
          scale_fill_manual(values = c(high ="#CC0033",low ="#0066CC", avg ="#33CC66"))+
          scale_colour_manual(values = c(high ="#CC0033",low ="#0066CC", avg ="#33CC66"))+
          theme(legend.position="none")+labs(x="Cluster")+theme(axis.title.x=element_text(size=10))
        if(j==1){
          plotObj<-plotObj+labs(y="Mean")+theme(axis.title.y=element_text(size=10,angle=90))
        }else{
          plotObj<-plotObj+theme(axis.title.y=element_blank())
        }
        plotObj<-plotObj+labs(title=covNames[j],plot.title=element_text(size=10))
        plotObj<-plotObj+
          theme(plot.margin=unit(c(0.5,ifelse(j==nCovariates,1,0),0.5,ifelse(j==1,0.5,0)),'lines'))+
          theme(plot.margin=unit(c(0,0,0,0),'lines'))

        print(plotObj,vp=viewport(layout.pos.row=1:3,layout.pos.col=j+plotRiskFlag))

        # Plot the variances
        profileDF<-data.frame("sigma"=c(),"cluster"=c(),"sigmaMean"=c(),
                              "lowerSigma"=c(),"upperSigma"=c(),"fillColor"=c())
        if (nContinuousCovs==1){
          sigmaMat<-profileStdDev[,meanSortIndex,1,1]
        } else {
          sigmaMat<-profileStdDev[,meanSortIndex,(j-nDiscreteCovs),(j-nDiscreteCovs)]
        }
        sigmaMeans<-apply(sigmaMat,2,mean)
        sigmaMean<-sum(sigmaMeans*clusterSizes)/sum(clusterSizes)
        sigmaLower<-apply(sigmaMat,2,quantile,0.05)
        sigmaUpper<-apply(sigmaMat,2,quantile,0.95)
        # The next line is to avoid outliers spoiling plot scales
        plotMax<-max(sigmaUpper)

        # Get the plot colors
        sigmaColor<-ifelse(sigmaLower>rep(sigmaMean,length(sigmaLower)),"high",
                           ifelse(sigmaUpper<rep(sigmaMean,length(sigmaUpper)),"low","avg"))
        for(c in whichClusters){
          plotSigma<-sigmaMat[,c]
          plotSigma<-plotSigma[plotSigma<plotMax]
          nPoints<-length(plotSigma)
          profileDF<-rbind(profileDF,data.frame("sigma"=plotSigma,"cluster"=rep(c,nPoints),
                                                "meanSigma"=rep(sigmaMean,nPoints),
                                                "lowerSigma"=rep(sigmaLower[c],nPoints),
                                                "upperSigma"=rep(sigmaUpper[c],nPoints),
                                                "fillColor"=rep(sigmaColor[c],nPoints)))
        }
        rownames(profileDF)<-seq(1,nrow(profileDF),1)
        plotObj<-ggplot(profileDF)
        plotObj<-plotObj+geom_hline(aes(x=as.factor(cluster),y=sigma,yintercept=meanSigma))
        plotObj<-plotObj+geom_boxplot(aes(x=as.factor(cluster),y=sigma,fill=as.factor(fillColor)),outlier.size=0.5)
        plotObj<-plotObj+geom_point(aes(x=as.factor(cluster),y=lowerSigma,colour=as.factor(fillColor)),size=1.5)
        plotObj<-plotObj+geom_point(aes(x=as.factor(cluster),y=upperSigma,colour=as.factor(fillColor)),size=1.5)
        plotObj<-plotObj+
          scale_fill_manual(values = c(high ="#CC0033",low ="#0066CC", avg ="#33CC66"))+
          scale_colour_manual(values = c(high ="#CC0033",low ="#0066CC", avg ="#33CC66"))+
          theme(legend.position="none")+labs(x="Cluster")+theme(axis.title.x=element_text(size=10))
        if(j==1){
          plotObj<-plotObj+labs(y="Std Dev")+theme(axis.title.y=element_text(size=10,angle=90))
        }else{
          plotObj<-plotObj+theme(axis.title.y=element_blank())
        }
        plotObj<-plotObj+
          theme(plot.margin=unit(c(0.5,ifelse(j==nCovariates,1,0),0.5,ifelse(j==1,0.5,0)),'lines'))+
          theme(plot.margin=unit(c(0,0,0,0),'lines'))
        print(plotObj,vp=viewport(layout.pos.row=4:6,layout.pos.col=j+plotRiskFlag))
      }
    }
  }
  dev.off()

  if(yModel=='LME'){
    longFile <- paste(strsplit(outFile,"\\.")[[1]][1],'-trajectories-data.png',sep="")
    png(longFile,width=1200,height=800)
    plotLayout<-grid.layout(ncol = 1, nrow = 1)
    grid.newpage()
    pushViewport(viewport(layout = plotLayout))

    GPDF<-data.frame("time"=c(),"mu"=c(),"cluster"=c(),"sigma"=c(),"fillColor"=c())
    times <- longMat$time
    yData <- longMat$outcome
    tTimes <- seq(min(times),max(times),length.out=41)
    palette <- rainbow(max(whichClusters))
    times_c <- list()
    yData_c <- list()
    mu <- rep(0,length(tTimes))

    if(nFixedEffects>0){
      if(nFixedEffects==1){
        betamean <-  mean(betaArray[,1,])
      }else{
        betamean <-  as.vector(colMeans(betaArray))
      }

      if(all(timevar %in% fixedEffectsNames)){
        ind_time <- which(fixedEffectsNames %in% c("intercept",timevar))
        mat_times <- rep(1,length(tTimes))
        for(jj in 1:length(timevar))
          mat_times<- cbind(mat_times, sapply(tTimes, function(x) x^jj))
        mu <- betamean[c(1,ind_time)] %*% mat_times
      }

      if(nFixedEffects==1){
        mu <- mu + rep(mean(betaArray[,1,]) * profile_X[[1]],length(tTimes))
      }else{
        mu <- mu + rep(colMeans(betaArray) %*% profile_X,length(tTimes))
      }
    }

    for(c in whichClusters){
      betamix_mean <- colMeans(betamixArray[,c,])
      jj=1
      if(all(timevar %in% fixedEffectsNames_clust)){
        ind_time <- which(fixedEffectsNames_clust %in% c("intercept",timevar))
        mat_times <- rep(1,length(tTimes))
        for(jjj in 1:length(timevar))
          mat_times<- cbind(mat_times, sapply(tTimes, function(x) x^jjj))
        mu <- mu + mat_times %*%betamix_mean[ind_time]
        jj = jj + length(ind_time)
      }

      if(length(fixedEffectsNames_clust)>=jj){
        for(other_spec_cov in jj:length(fixedEffectsNames_clust))
          mu <- mu + rep(betamix_mean[other_spec_cov] * profile_X[[fixedEffectsNames_clust[jj]]],length(tTimes))
      }

      GPDF <- rbind(GPDF,data.frame("time"=tTimes,"mu"=mu,
                                    "cluster"=rep(c,times=length(tTimes)),
                                    "sigma"=mu,#1.645*sqrt(diag(params$GPSigma)),
                                    "inf"=mu,#params$mu+longMean-1.645*sqrt(diag(params$GPSigma)),
                                    "sup"=mu,#params$mu+longMean+1.645*sqrt(diag(params$GPSigma)),
                                    "fillColor"=riskColor[c]))

      #       dev.off()
      #       par(mfrow=c(1,1))
      #       c=1
      #       i=min(which(clustering==c))
      #       plot(yData[tMat[i,1]:tMat[i,2]]+longMean~times[tMat[i,1]:tMat[i,2]],type='l',ylim=c(min(yData)+longMean,max(yData)+longMean),xlim=c(min(times),max(times)))
      #
      #       for(i in 2:nSubjects){
      #         if(clustering[i]==c)
      #           lines(yData[tMat[i,1]:tMat[i,2]]+longMean~times[tMat[i,1]:tMat[i,2]])
      #       }
      #       lines(params$mu+longMean~tTimes,lwd=3,col=c)
      # c=c+1

    }

    rownames(GPDF)<-seq(1,nrow(GPDF),1)
    plotObj <- ggplot(GPDF)
    for(i in 1:nSubjects){
      df <- (data.frame(x=times[tMat[i,1]:tMat[i,2]],y=yData[tMat[i,1]:tMat[i,2]]+longMean,cluster=rep(clustering[i],tMat[i,2]-tMat[i,1]+1)))
      if(bycol){
        #  if(clustering[i]==2)
        plotObj <- plotObj + geom_line(data=df,aes(x,y,colour=as.factor(cluster)),alpha = 0.3)#colour='azure4')
      }else{
        plotObj <- plotObj + geom_line(data=df,aes(x,y),colour='azure4')
      }
    }
    plotObj <- plotObj + geom_line(aes(x=time,y=mu,group=cluster,colour=as.factor(cluster)),size=2)
    plotObj <- plotObj + geom_line(aes(x=time,y=sup,group=cluster,colour=as.factor(cluster)))
    plotObj <- plotObj + geom_line(aes(x=time,y=inf,group=cluster,colour=as.factor(cluster)))
    plotObj <- plotObj + labs(color="Cluster")
    plotObj <- plotObj + labs(y="Outcome\n")+theme(axis.title.y=element_text(size=30,angle=90))
    plotObj <- plotObj + labs(x="\nTime")+theme(axis.title.x=element_text(size=30))
    plotObj <- plotObj + coord_cartesian(ylim = c(min(GPDF$inf,yData),max(GPDF$sup,yData)))
    plotObj <- plotObj + theme(axis.text.x=element_text(size=30)) + theme(axis.text.y=element_text(size=30))
    plotObj <- plotObj + theme(axis.line.x = element_line(colour="black",size=2),
                               axis.line.y = element_line(colour="black",size=2))
    plotObj <- plotObj + theme(legend.title=element_text(size=30)) + theme(legend.text=element_text(size=25))
    plotObj <- plotObj + theme(plot.title=element_text(size=30))
    #   plotObj <- plotObj + coord_cartesian(ylim = c(0,max(GPDF$mu+GPDF$sigma)))#min(GPDF$mu-GPDF$sigma),max(GPDF$mu+GPDF$sigma)))
    #   plotObj <- plotObj + coord_cartesian(ylim = c(min(GPDF$inf),max(GPDF$sup)))
    print(plotObj,vp=viewport(layout.pos.row=1,layout.pos.col=1))
    dev.off()

    ##//RJ data plot
    longFile <- paste(strsplit(outFile,"\\.")[[1]][1],'-trajectories.png',sep="")
    png(longFile,width=1200,height=800)
    plotLayout<-grid.layout(ncol = 1, nrow = 1)
    grid.newpage()
    pushViewport(viewport(layout = plotLayout))
    plotObj <- ggplot(GPDF)
    plotObj <- plotObj + geom_line(aes(x=time,y=mu,group=cluster,colour=as.factor(cluster)),size=2)
    #plotObj <- plotObj + geom_line(aes(x=time,y=sup,group=cluster,colour=as.factor(cluster)))
    #plotObj <- plotObj + geom_line(aes(x=time,y=inf,group=cluster,colour=as.factor(cluster)))
    plotObj <- plotObj + labs(color="Cluster")
    plotObj <- plotObj + labs(y="Outcome")+theme(axis.title.y=element_text(size=40,angle=90))
    plotObj <- plotObj + labs(x="Time")+theme(axis.title.x=element_text(size=40))
    plotObj <- plotObj + coord_cartesian(ylim = c(min(GPDF$inf),max(GPDF$sup)))
    plotObj <- plotObj + theme(axis.text.x=element_text(size=30)) + theme(axis.text.y=element_text(size=30))
    plotObj <- plotObj + theme(axis.line.x = element_line(colour="black",size=2),
                               axis.line.y = element_line(colour="black",size=2))
    plotObj <- plotObj + theme(legend.title=element_text(size=30)) + theme(legend.text=element_text(size=30))
    plotObj <- plotObj + theme(plot.title=element_text(size=30))
    #plotObj <- plotObj + coord_cartesian(ylim = c(min(GPDF$inf),max(GPDF$sup)))
    print(plotObj,vp=viewport(layout.pos.row=1,layout.pos.col=1))
    dev.off()

    ##//AR data plot
    if(double_plot){
      longFile <- paste(strsplit(outFile,"\\.")[[1]][1],'-all_trajectories.png',sep="")
      png(longFile,width=1200,height=800)
      plotObj <- ggplot(GPDF) + theme_bw()
      plotObj <- plotObj + geom_line(aes(x=time,y=mu,group=cluster,colour=as.factor(cluster)),size=2)
      #plotObj <- plotObj + geom_line(aes(x=time,y=sup,group=cluster,colour=as.factor(cluster)))
      #plotObj <- plotObj + geom_line(aes(x=time,y=inf,group=cluster,colour=as.factor(cluster)))
      plotObj <- plotObj + labs(color="Cluster")
      plotObj <- plotObj + labs(y="Cluster-specific outcome")+theme(axis.title.y=element_text(size=30,angle=90))
      plotObj <- plotObj + labs(x="Time") + theme(axis.title.x=element_text(size=30))
      plotObj <- plotObj + coord_cartesian(ylim = c(min(GPDF$inf),max(GPDF$sup)))
      plotObj <- plotObj + theme(axis.text.x=element_text(size=20, color = "black"), axis.text.y=element_text(size=20, color = "black"))

      plotObj <- plotObj + theme(axis.line.x = element_line(colour="black",size=1),
                                 axis.line.y = element_line(colour="black",size=1))
      #plotObj <- plotObj + theme(legend.title=element_text(size=10)) + theme(legend.text=element_text(size=20))
      #plotObj <- plotObj + theme(plot.title=element_text(size=10))

      p1<- plotObj + #geom_line(aes(x=time,y=mu,group=cluster),size=1) +
        #theme(axis.text.x=element_text(size=10),axis.text.y=element_text(size=10))+
        facet_wrap(~cluster,ncol=1, strip.position="left")+theme(legend.position = "none")
      #theme(axis.title.x=element_text(size=20), axis.title.y=element_text(size=20))
      p2 <- plotProfilesByCluster_longi(riskProfObj, rhoMinimum =0.1, useProfileStar=F)
      p2 <-p2 + theme(axis.text.x = element_text(angle = 90, size=17))
      #p2 + theme(axis.text.x = element_blank())

      plot_p2 <- plot_grid(p1,p2, align= 'h', axis='b', rel_widths = c(1,2))
      print(plot_p2,vp=viewport(layout.pos.row=1,layout.pos.col=1))

      dev.off()

      longFile <- paste(strsplit(outFile,"\\.")[[1]][1],'-all_trajectories.pdf',sep="")
      pdf(longFile, width=12, height=10 )
      p2 <- p2 +theme(legend.text=element_text(size=15),legend.title=element_text(size=20))
      p1 <- p1 + geom_vline(xintercept=c(0,13.33,26.67,40),linetype="dashed")
      p1 <- p1 + geom_vline(xintercept=c(0.7,4.7,7.3,9.3,10.7,12.7),linetype="dashed",col="gray")
      timescale<-function(x){x*5}
      p1 + scale_x_continuous(labels=timescale)
      print(plot_p2,vp=viewport(layout.pos.row=1,layout.pos.col=1))
      dev.off()


      longFile <- paste(strsplit(outFile,"\\.")[[1]][1],'-all_trajectories-data.pdf',sep="")
      pdf(longFile, width=12, height=10 )
      #png(longFile,width=1200,height=800)


      plotObj <- ggplot(GPDF)+theme_bw()
      for(i in 1:nSubjects){
        df <- (data.frame(x=times[tMat[i,1]:tMat[i,2]],y=yData[tMat[i,1]:tMat[i,2]]+longMean,cluster=rep(clustering[i],tMat[i,2]-tMat[i,1]+1)))
        if(bycol){
          #  if(clustering[i]==2)
          plotObj <- plotObj + geom_line(data=df,aes(x,y,colour=as.factor(cluster)),alpha = 0.3)#colour='azure4')
        }else{
          plotObj <- plotObj + geom_line(data=df,aes(x,y),colour='azure4')
        }
      }
      #plotObj <- plotObj + geom_line(aes(x=time,y=mu,group=cluster,colour=as.factor(cluster)),size=2)
      #plotObj <- plotObj + geom_line(aes(x=time,y=sup,group=cluster,colour=as.factor(cluster)))
      #plotObj <- plotObj + geom_line(aes(x=time,y=inf,group=cluster,colour=as.factor(cluster)))
      plotObj <- plotObj + labs(color="Cluster")
      plotObj <- plotObj + labs(y="Cluster-specific outcome")+theme(axis.title.y=element_text(size=30,angle=90))
      plotObj <- plotObj + labs(x="Time")+theme(axis.title.x=element_text(size=30))
      plotObj <- plotObj + coord_cartesian(ylim = c(min(GPDF$inf,yData),max(GPDF$sup,yData)))
      plotObj <- plotObj + theme(axis.text.x=element_text(size=20, color = "black"), axis.text.y=element_text(size=20, color = "black"))
      #plotObj <- plotObj + theme(axis.line.x = element_line(colour="black",size=1),
      #axis.line.y = element_line(colour="black",size=1))
      #plotObj <- plotObj + theme(legend.title=element_text(size=20)) #+ theme(legend.text=element_text(size=25))
      #plotObj <- plotObj + theme(plot.title=element_text(size=20))

      p1b<- plotObj + #geom_line(aes(x=time,y=mu,group=cluster),size=1) +
        #theme(axis.text.x=element_text(size=10),axis.text.y=element_text(size=10))+
        facet_wrap(~cluster,ncol=1, strip.position="left")+
        theme(legend.position = "none")
      plot_p1b <- plot_grid(p1b,p2, align= 'h', axis='b', rel_widths = c(1,2))
      print(plot_p1b,vp=viewport(layout.pos.row=1,layout.pos.col=1))

      dev.off()

      longFile <- paste(strsplit(outFile,"\\.")[[1]][1],'-all_trajectories-data.pdf',sep="")
      pdf(longFile, width=12, height=10 )
      p1b <- p1b + geom_vline(xintercept=c(0,13.33,26.67,40),linetype="dashed")
      p1b <- p1b + geom_vline(xintercept=c(0.7,4.7,7.3,9.3,10.7,12.7),linetype="dashed",col="gray")
      timescale<-function(x){x*5}
      p1b <- p1b + scale_x_continuous(labels=timescale)
      plot_p1b <- plot_grid(p1b,p2, align= 'h', axis='b', rel_widths = c(1,2))
      print(plot_p1b,vp=viewport(layout.pos.row=1,layout.pos.col=1))
      dev.off()
      browser()
    }
  }

  ##//RJ new plot for longitudinal data
  if(yModel=='Longitudinal'){
    longFile <- paste(strsplit(outFile,"\\.")[[1]][1],'-trajectories-data.png',sep="")
    png(longFile,width=1200,height=800)
    plotLayout<-grid.layout(ncol = 1, nrow = 1)
    grid.newpage()
    pushViewport(viewport(layout = plotLayout))

    GPDF<-data.frame("time"=c(),"mu"=c(),"cluster"=c(),"sigma"=c(),"fillColor"=c())
    times <- longMat$time
    yData <- longMat$outcome
    tTimes <- seq(min(times),max(times),length.out=41)
    palette <- rainbow(max(whichClusters))
    times_c <- list()
    yData_c <- list()

    for(c in whichClusters){
      times_c0 <- unlist(lapply(1:nSubjects,function(x)if(clustering[x]==c) times[tMat[x,1]:tMat[x,2]]))
      times_c <- all_times[sapply(times_c0, function(x) which(abs(all_times-x)==min(abs(all_times-x))))]
      yData_c <- unlist(lapply(1:nSubjects,function(x)if(clustering[x]==c) yData[tMat[x,1]:tMat[x,2]]))

      if(sampleGPmean){#sampleGPmean
        tTimes0 <- unique(times)[order(unique(times))]
        tTimes <- unique(all_times[sapply(tTimes0, function(x) which(abs(all_times-x)==min(abs(all_times-x))))])
        mu <- GPmeanMeans[c,]
        sigma <- exp(LMeans[c,3]/2)
        inf <- GPmeanMeans_inf[c,]
        sup <- GPmeanMeans_sup[c,]
        GPDF <- rbind(GPDF,data.frame("time"=tTimes,"mu"=mu+longMean,
                                      "cluster"=rep(c,times=length(tTimes)),
                                      "inf"=inf+longMean,
                                      "sup"=sup+longMean,
                                      "sigma"=rep(sigma,length(tTimes)),
                                      "fillColor"=riskColor[c]))

      }else{
        params <- GPfun(t0=times_c,ts=tTimes,y=yData_c,Lp=LMeans[c,],kernel)
        GPDF <- rbind(GPDF,data.frame("time"=tTimes,"mu"=params$mu+longMean,
                                      "cluster"=rep(c,times=length(tTimes)),
                                      "sigma"=1.645*sqrt(diag(params$GPSigma)),
                                      "inf"=params$mu+longMean-1.645*sqrt(diag(params$GPSigma)),
                                      "sup"=params$mu+longMean+1.645*sqrt(diag(params$GPSigma)),
                                      "fillColor"=riskColor[c]))
      }

      #       dev.off()
      #       par(mfrow=c(1,1))
      #       c=1
      #       i=min(which(clustering==c))
      #       plot(yData[tMat[i,1]:tMat[i,2]]+longMean~times[tMat[i,1]:tMat[i,2]],type='l',ylim=c(min(yData)+longMean,max(yData)+longMean),xlim=c(min(times),max(times)))
      #
      #       for(i in 2:nSubjects){
      #         if(clustering[i]==c)
      #           lines(yData[tMat[i,1]:tMat[i,2]]+longMean~times[tMat[i,1]:tMat[i,2]])
      #       }
      #       lines(params$mu+longMean~tTimes,lwd=3,col=c)
      # c=c+1
    }
    rownames(GPDF)<-seq(1,nrow(GPDF),1)
    plotObj <- ggplot(GPDF)
    for(i in 1:nSubjects){
      df <- (data.frame(x=times[tMat[i,1]:tMat[i,2]],y=yData[tMat[i,1]:tMat[i,2]]+longMean,cluster=rep(clustering[i],tMat[i,2]-tMat[i,1]+1)))
      if(bycol){
        #  if(clustering[i]==2)
        plotObj <- plotObj + geom_line(data=df,aes(x,y,colour=as.factor(cluster)),alpha = 0.3)#colour='azure4')
      }else{
        plotObj <- plotObj + geom_line(data=df,aes(x,y),colour='azure4')
      }
    }
    plotObj <- plotObj + geom_line(aes(x=time,y=mu,group=cluster,colour=as.factor(cluster)),size=2)
    plotObj <- plotObj + geom_line(aes(x=time,y=sup,group=cluster,colour=as.factor(cluster)))
    plotObj <- plotObj + geom_line(aes(x=time,y=inf,group=cluster,colour=as.factor(cluster)))
    plotObj <- plotObj + labs(color="Cluster")
    plotObj <- plotObj + labs(y="Outcome\n")+theme(axis.title.y=element_text(size=30,angle=90))
    plotObj <- plotObj + labs(x="\nTime")+theme(axis.title.x=element_text(size=30))
    plotObj <- plotObj + coord_cartesian(ylim = c(min(GPDF$inf,yData),max(GPDF$sup,yData)))
    plotObj <- plotObj + theme(axis.text.x=element_text(size=30)) + theme(axis.text.y=element_text(size=30))
    plotObj <- plotObj + theme(axis.line.x = element_line(colour="black",size=2),
                               axis.line.y = element_line(colour="black",size=2))
    plotObj <- plotObj + theme(legend.title=element_text(size=30)) + theme(legend.text=element_text(size=25))
    plotObj <- plotObj + theme(plot.title=element_text(size=30))
    #   plotObj <- plotObj + coord_cartesian(ylim = c(0,max(GPDF$mu+GPDF$sigma)))#min(GPDF$mu-GPDF$sigma),max(GPDF$mu+GPDF$sigma)))
    #   plotObj <- plotObj + coord_cartesian(ylim = c(min(GPDF$inf),max(GPDF$sup)))
    print(plotObj,vp=viewport(layout.pos.row=1,layout.pos.col=1))
    dev.off()

    ##//RJ data plot
    longFile <- paste(strsplit(outFile,"\\.")[[1]][1],'-trajectories.png',sep="")
    png(longFile,width=1200,height=800)
    plotLayout<-grid.layout(ncol = 1, nrow = 1)
    grid.newpage()
    pushViewport(viewport(layout = plotLayout))
    plotObj <- ggplot(GPDF)
    plotObj <- plotObj + geom_line(aes(x=time,y=mu,group=cluster,colour=as.factor(cluster)),size=2)
    #plotObj <- plotObj + geom_line(aes(x=time,y=sup,group=cluster,colour=as.factor(cluster)))
    #plotObj <- plotObj + geom_line(aes(x=time,y=inf,group=cluster,colour=as.factor(cluster)))
    plotObj <- plotObj + labs(color="Cluster")
    plotObj <- plotObj + labs(y="Outcome")+theme(axis.title.y=element_text(size=40,angle=90))
    plotObj <- plotObj + labs(x="Time")+theme(axis.title.x=element_text(size=40))
    plotObj <- plotObj + coord_cartesian(ylim = c(min(GPDF$inf),max(GPDF$sup)))
    plotObj <- plotObj + theme(axis.text.x=element_text(size=30)) + theme(axis.text.y=element_text(size=30))
    plotObj <- plotObj + theme(axis.line.x = element_line(colour="black",size=2),
                               axis.line.y = element_line(colour="black",size=2))
    plotObj <- plotObj + theme(legend.title=element_text(size=30)) + theme(legend.text=element_text(size=30))
    plotObj <- plotObj + theme(plot.title=element_text(size=30))
    #plotObj <- plotObj + coord_cartesian(ylim = c(min(GPDF$inf),max(GPDF$sup)))
    print(plotObj,vp=viewport(layout.pos.row=1,layout.pos.col=1))
    dev.off()

    ##//AR data plot
    if(double_plot){
      longFile <- paste(strsplit(outFile,"\\.")[[1]][1],'-all_trajectories.png',sep="")
      png(longFile,width=1200,height=800)
      plotObj <- ggplot(GPDF) + theme_bw()
      plotObj <- plotObj + geom_line(aes(x=time,y=mu,group=cluster,colour=as.factor(cluster)),size=2)
      #plotObj <- plotObj + geom_line(aes(x=time,y=sup,group=cluster,colour=as.factor(cluster)))
      #plotObj <- plotObj + geom_line(aes(x=time,y=inf,group=cluster,colour=as.factor(cluster)))
      plotObj <- plotObj + labs(color="Cluster")
      plotObj <- plotObj + labs(y="Cluster-specific outcome")+theme(axis.title.y=element_text(size=30,angle=90))
      plotObj <- plotObj + labs(x="Time") + theme(axis.title.x=element_text(size=30))
      plotObj <- plotObj + coord_cartesian(ylim = c(min(GPDF$inf),max(GPDF$sup)))
      plotObj <- plotObj + theme(axis.text.x=element_text(size=20, color = "black"), axis.text.y=element_text(size=20, color = "black"))

      plotObj <- plotObj + theme(axis.line.x = element_line(colour="black",size=1),
                                 axis.line.y = element_line(colour="black",size=1))
      #plotObj <- plotObj + theme(legend.title=element_text(size=10)) + theme(legend.text=element_text(size=20))
      #plotObj <- plotObj + theme(plot.title=element_text(size=10))

      p1<- plotObj + #geom_line(aes(x=time,y=mu,group=cluster),size=1) +
        #theme(axis.text.x=element_text(size=10),axis.text.y=element_text(size=10))+
        facet_wrap(~cluster,ncol=1, strip.position="left")+theme(legend.position = "none")
      #theme(axis.title.x=element_text(size=20), axis.title.y=element_text(size=20))
      p2 <- plotProfilesByCluster_longi(riskProfObj, rhoMinimum = 0.1, useProfileStar=F)

      #table(clusObj$clusObjRunInfoObj$xMat$FKH1)


      p2bis <- plotProfilesByCluster(riskProfObj, rhoMinimum = 0.1, useProfileStar=F)
      p2 <-p2 + theme(axis.text.x = element_text(angle = 90, size=17))
      #p2 + theme(axis.text.x = element_blank())

      plot_p2 <- plot_grid(p1,p2, align= 'h', axis='b', rel_widths = c(1,2))
      print(plot_p2,vp=viewport(layout.pos.row=1,layout.pos.col=1))

      dev.off()

      # longFile <- paste(strsplit(outFile,"\\.")[[1]][1],'-all_trajectories.pdf',sep="")
      # pdf(longFile, width=12, height=10 )
      # p2 <- p2 +theme(legend.text=element_text(size=15),legend.title=element_text(size=20))
      # p1 <- p1 + geom_vline(xintercept=c(0,13.33,26.67,40),linetype="dashed")
      # p1 <- p1 + geom_vline(xintercept=c(0.7,4.7,7.3,9.3,10.7,12.7),linetype="dashed",col="gray")
      # timescale<-function(x){x*5}
      # p1 + scale_x_continuous(labels=timescale)
      # plot_p2 <- plot_grid(p1,p2, align= 'h', axis='b', rel_widths = c(1,2))
      # print(plot_p2,vp=viewport(layout.pos.row=1,layout.pos.col=1))
      # dev.off()


      longFile <- paste(strsplit(outFile,"\\.")[[1]][1],'-all_trajectories-data.png',sep="")
      png(longFile,width=1200,height=800)


      plotObj <- ggplot(GPDF)+theme_bw()
      for(i in 1:nSubjects){
        df <- (data.frame(x=times[tMat[i,1]:tMat[i,2]],y=yData[tMat[i,1]:tMat[i,2]]+longMean,cluster=rep(clustering[i],tMat[i,2]-tMat[i,1]+1)))
        if(bycol){
          #  if(clustering[i]==2)
          plotObj <- plotObj + geom_line(data=df,aes(x,y,colour=as.factor(cluster)),alpha = 0.3)#colour='azure4')
        }else{
          plotObj <- plotObj + geom_line(data=df,aes(x,y),colour='azure4')
        }
      }
      #plotObj <- plotObj + geom_line(aes(x=time,y=mu,group=cluster,colour=as.factor(cluster)),size=2)
      #plotObj <- plotObj + geom_line(aes(x=time,y=sup,group=cluster,colour=as.factor(cluster)))
      #plotObj <- plotObj + geom_line(aes(x=time,y=inf,group=cluster,colour=as.factor(cluster)))
      plotObj <- plotObj + labs(color="Cluster")
      plotObj <- plotObj + labs(y="Cluster-specific outcome")+theme(axis.title.y=element_text(size=30,angle=90),axis.title.x=element_text(size=30))
      plotObj <- plotObj + labs(x="Time")
      plotObj <- plotObj + coord_cartesian(ylim = c(min(GPDF$inf,yData),max(GPDF$sup,yData)))
      plotObj <- plotObj + theme(axis.text.x=element_text(size=20, color = "black"), axis.text.y=element_text(color = "black"))
      #plotObj <- plotObj + theme(axis.line.x = element_line(colour="black",size=1),
      #axis.line.y = element_line(colour="black",size=1))
      #plotObj <- plotObj + theme(legend.title=element_text(size=20)) #+ theme(legend.text=element_text(size=25))
      #plotObj <- plotObj + theme(plot.title=element_text(size=20))

      p1b<- plotObj + #geom_line(aes(x=time,y=mu,group=cluster),size=1) +
        #theme(axis.text.x=element_text(size=10),axis.text.y=element_text(size=10))+
        facet_wrap(~cluster,ncol=1, strip.position="left")+
        theme(legend.position = "none")
      plot_p1b <- plot_grid(p1b,p2, align= 'h', axis='b', rel_widths = c(1,2))
      print(plot_p1b,vp=viewport(layout.pos.row=1,layout.pos.col=1))

      dev.off()

      # longFile <- paste(strsplit(outFile,"\\.")[[1]][1],'-all_trajectories-data.pdf',sep="")
      # pdf(longFile, width=12, height=10 )
      # p1b <- p1b + geom_vline(xintercept=c(0,13.33,26.67,40),linetype="dashed")
      # p1b <- p1b + geom_vline(xintercept=c(0.7,4.7,7.3,9.3,10.7,12.7),linetype="dashed",col="gray")
      # timescale<-function(x){x*5}
      # p1b <- p1b + scale_x_continuous(labels=timescale)
      # plot_p1b <- plot_grid(p1b,p2, align= 'h', axis='b', rel_widths = c(1,2))
      # print(plot_p1b,vp=viewport(layout.pos.row=1,layout.pos.col=1))
      # dev.off()

      longFile <- paste(strsplit(outFile,"\\.")[[1]][1],'-all.png',sep="")
      png(longFile,width=1200,height=800)
      p3       <- Plot_XY_graph(clusObj, nameX = "Transcription Factors", nameY = "Genes" ) +
        theme(axis.title.x=element_text(size=30, color = "black"),
              axis.title.y=element_text(size=30, color = "black"),
              legend.title = element_text(size=25),
              legend.text = element_text(size=25))
      plot_p1b <- plot_grid(p1b,p3, align= 'h', axis='b', rel_widths = c(1,2))
      print(plot_p1b,vp=viewport(layout.pos.row=1,layout.pos.col=1))
      dev.off()
    }
  }

  if(yModel=="MVN"){
    png(paste(strsplit(outFile,"\\.")[[1]][1],'-MVN.png',sep=""),width=1200,height=800)
    plotLayout<-grid.layout(ncol = nOutcomes, nrow = 6)
    grid.newpage()
    pushViewport(viewport(layout = plotLayout))

    profileDFmean<-data.frame("time"=c(),"mu"=c(),"cluster"=c(),"muMean"=c(),
                              "lowerMu"=c(),"upperMu"=c(),"fillColor"=c())
    for(j in 1:nOutcomes){
      # Plot the means
      profileDF<-data.frame("time"=c(),"mu"=c(),"cluster"=c(),"muMean"=c(),
                            "lowerMu"=c(),"upperMu"=c(),"fillColor"=c())
      muMat<-MVNmuArray[,,j]
      muMeans<-apply(muMat,2,mean)
      muMean<-sum(muMeans*clusterSizes)/sum(clusterSizes)
      muLower<-apply(muMat,2,quantile,0.05)
      muUpper<-apply(muMat,2,quantile,0.95)
      # The next line is to avoid outliers spoiling plot scales
      plotMax<-max(muUpper)
      plotMin<-min(muLower)

      # Get the plot colors
      muColor<-ifelse(muLower>rep(muMean,length(muLower)),"high",
                      ifelse(muUpper<rep(muMean,length(muUpper)),"low","avg"))
      for(c in whichClusters){
        plotMu<-muMat[,c]
        plotMu<-plotMu[plotMu<plotMax&plotMu>plotMin]
        nPoints<-length(plotMu)
        profileDF<-rbind(profileDF,data.frame("time"=j,"mu"=plotMu,"cluster"=rep(c,nPoints),
                                              "meanMu"=rep(muMean,nPoints),
                                              "lowerMu"=rep(muLower[c],nPoints),
                                              "upperMu"=rep(muUpper[c],nPoints),
                                              "fillColor"=rep(muColor[c],nPoints)))
        profileDFmean <- rbind(profileDFmean,data.frame("time"=j,"mu"=muMeans[c],"cluster"=c,
                                                        "meanMu"=muMeans[c],
                                                        "lowerMu"=muLower[c],
                                                        "upperMu"=muUpper[c],
                                                        "fillColor"=muColor[c]))
      }


      rownames(profileDF)<-seq(1,nrow(profileDF),1)

      plotObj<-ggplot(profileDF)
      plotObj<-plotObj+geom_hline(aes(x=as.factor(cluster),y=mu,yintercept=meanMu))
      plotObj<-plotObj+geom_boxplot(aes(x=as.factor(cluster),y=mu,fill=as.factor(cluster)),outlier.size=0.5)
      plotObj<-plotObj+geom_point(aes(x=as.factor(cluster),y=lowerMu,colour=as.factor(fillColor)),size=1.5)
      plotObj<-plotObj+geom_point(aes(x=as.factor(cluster),y=upperMu,colour=as.factor(fillColor)),size=1.5)
      plotObj<-plotObj+
        #scale_fill_manual(values = c(high ="#CC0033",low ="#0066CC", avg ="#33CC66"))+
        #scale_colour_manual(values = c(high ="#CC0033",low ="#0066CC", avg ="#33CC66"))+
        theme(legend.position="none")+labs(x="Cluster")+theme(axis.title.x=element_text(size=10))
      if(j==1){
        plotObj<-plotObj+labs(y="Mean")+theme(axis.title.y=element_text(size=10,angle=90))
      }else{
        plotObj<-plotObj+theme(axis.title.y=element_blank())
      }
      plotObj<-plotObj+labs(title=outcome[j],plot.title=element_text(size=20))+
        theme(plot.title = element_text(size=20))
      plotObj<-plotObj+
        theme(plot.margin=unit(c(0.5,ifelse(j==nOutcomes,1,0),0.5,ifelse(j==1,0.5,0)),'lines'))
      #+ theme(plot.margin=unit(c(0,0,0,0),'lines'))
      print(plotObj,vp=viewport(layout.pos.row=1:3,layout.pos.col=j))

      # Plot the variances
      profileDF<-data.frame("sigma"=c(),"cluster"=c(),"sigmaMean"=c(),
                            "lowerSigma"=c(),"upperSigma"=c(),"fillColor"=c())
      # MVNSigmaArray is a vector of the lower triangular covariance matrix, ordered by row.
      # The index of the variance is the j-th triangular number.
      sigmaMat<-sqrt(MVNSigmaArray[,,j*(j+1)/2])
      sigmaLower<-apply(sigmaMat,2,quantile,0.05)
      sigmaUpper<-apply(sigmaMat,2,quantile,0.95)
      sigmaMeans<-apply(sigmaMat,2,mean)
      sigmaMean<-sum(sigmaMeans*clusterSizes)/sum(clusterSizes)
      # The next line is to avoid outliers spoiling plot scales
      plotMax<-max(sigmaUpper)

      # Get the plot colors
      sigmaColor<-ifelse(sigmaLower>rep(sigmaMean,length(sigmaLower)),"high",
                         ifelse(sigmaUpper<rep(sigmaMean,length(sigmaUpper)),"low","avg"))
      for(c in whichClusters){
        plotSigma<-sigmaMat[,c]
        plotSigma<-plotSigma[plotSigma<plotMax]
        nPoints<-length(plotSigma)
        profileDF<-rbind(profileDF,data.frame("sigma"=plotSigma,"cluster"=rep(c,nPoints),
                                              "meanSigma"=rep(sigmaMean,nPoints),
                                              "lowerSigma"=rep(sigmaLower[c],nPoints),
                                              "upperSigma"=rep(sigmaUpper[c],nPoints),
                                              "fillColor"=rep(sigmaColor[c],nPoints)))
      }
      rownames(profileDF)<-seq(1,nrow(profileDF),1)
      plotObj<-ggplot(profileDF)
      plotObj<-plotObj+geom_hline(aes(x=as.factor(cluster),y=sigma,yintercept=meanSigma))
      plotObj<-plotObj+geom_boxplot(aes(x=as.factor(cluster),y=sigma,fill=as.factor(fillColor)),outlier.size=0.5)
      plotObj<-plotObj+geom_point(aes(x=as.factor(cluster),y=lowerSigma,colour=as.factor(fillColor)),size=1.5)
      plotObj<-plotObj+geom_point(aes(x=as.factor(cluster),y=upperSigma,colour=as.factor(fillColor)),size=1.5)
      plotObj<-plotObj+
        scale_fill_manual(values = c(high ="#CC0033",low ="#0066CC", avg ="#33CC66"))+
        scale_colour_manual(values = c(high ="#CC0033",low ="#0066CC", avg ="#33CC66"))+
        theme(legend.position="none")+labs(x="Cluster")+theme(axis.title.x=element_text(size=10))
      if(j==1){
        plotObj<-plotObj+labs(y="Std Dev")+theme(axis.title.y=element_text(size=10,angle=90))
      }else{
        plotObj<-plotObj+theme(axis.title.y=element_blank())
      }
      plotObj<-plotObj+
        theme(plot.margin=unit(c(0.5,ifelse(j==nOutcomes,1,0),0.5,ifelse(j==1,0.5,0)),'lines'))+
        theme(plot.margin=unit(c(0,0,0,0),'lines'))
      print(plotObj,vp=viewport(layout.pos.row=4:6,layout.pos.col=j))
    }
    dev.off()

    longFile <- paste(strsplit(outFile,"\\.")[[1]][1],'-all_trajectories.png',sep="")

    png(longFile,width=1200,height=800)
    plotObj <- ggplot(profileDFmean) + theme_bw()
    plotObj <- plotObj + geom_line(aes(x=time,y=mu,group=cluster,colour=as.factor(cluster)),size=2)
    plotObj <- plotObj + labs(color="Cluster")
    plotObj <- plotObj + labs(y="Cluster-specific outcome")+theme(axis.title.y=element_text(size=30,angle=90))
    plotObj <- plotObj + labs(x="Time") + theme(axis.title.x=element_text(size=30))
    plotObj <- plotObj + coord_cartesian(ylim = c(min(profileDFmean$mu),max(profileDFmean$mu)))
    plotObj <- plotObj + theme(axis.text.x=element_text(size=20, color = "black"), axis.text.y=element_text(size=20, color = "black"))

    plotObj <- plotObj + theme(axis.line.x = element_line(colour="black",size=1),
                               axis.line.y = element_line(colour="black",size=1))

    p1<- plotObj + facet_wrap(~cluster,ncol=1, strip.position="left")+theme(legend.position = "none")
    ##!! varSelect
    riskProfObj$riskProfClusObj$clusObjRunInfoObj$varSelect <- T
    p2 <- plotProfilesByCluster_longi(riskProfObj, rhoMinimum = 0.1, useProfileStar=F)

    table(clusObj$clusObjRunInfoObj$xMat$FKH1)


    p2bis <- plotProfilesByCluster_longi(riskProfObj, rhoMinimum = 0.1, useProfileStar=F)
    p2 <-p2 + theme(axis.text.x = element_text(angle = 90, size=17))

    plot_p2 <- plot_grid(p1,p2, align= 'h', axis='b', rel_widths = c(1,2))
    print(plot_p2,vp=viewport(layout.pos.row=1,layout.pos.col=1))

    dev.off()

    longFile <- paste(strsplit(outFile,"\\.")[[1]][1],'-all_trajectories-data.png',sep="")
    png(longFile,width=1200,height=800)
    #times<-(times-1)*5
    plotObj <- ggplot(profileDF)+theme_bw()
    for(i in 1:nSubjects){
      df <- (data.frame(x=times,y=yMat[i,],cluster=rep(clustering[i],nOutcomes)))
      if(bycol){
        plotObj <- plotObj + geom_line(data=df,aes(x,y,colour=as.factor(cluster)),alpha = 0.3)#colour='azure4')
      }else{
        plotObj <- plotObj + geom_line(data=df,aes(x,y),colour='azure4')
      }
    }
    plotObj <- plotObj + labs(color="Cluster")
    plotObj <- plotObj + labs(y="Cluster-specific outcome")+theme(axis.title.y=element_text(size=30,angle=90),axis.title.x=element_text(size=30))
    plotObj <- plotObj + labs(x="Time")
    plotObj <- plotObj + coord_cartesian(ylim = c(min(profileDFmean$mu,yMat),max(profileDFmean$mu,yMat)))
    plotObj <- plotObj + theme(axis.text.x=element_text(size=20, color = "black"), axis.text.y=element_text(color = "black"))

    p1b<- plotObj + facet_wrap(~cluster,ncol=1, strip.position="left")+
      theme(legend.position = "none", axis.text.y=element_text(size=20))
    plot_p1b <- plot_grid(p1b,p2, align= 'h', axis='b', rel_widths = c(1,2))
    print(plot_p1b,vp=viewport(layout.pos.row=1,layout.pos.col=1))

    dev.off()

    longFile <- paste(strsplit(outFile,"\\.")[[1]][1],'-all_trajectories-data.pdf',sep="")
    pdf(longFile, width=12, height=10 )
    p1b <- p1b + geom_vline(xintercept=c(0,13.33,26.67,40),linetype="dashed")
    p1b <- p1b + geom_vline(xintercept=c(0.7,4.7,7.3,9.3,10.7,12.7),linetype="dashed",col="gray")
    timescale<-function(x){x*5}
    p1b <- p1b + scale_x_continuous(labels=timescale)
    plot_p1b <- plot_grid(p1b,p2, align= 'h', axis='b', rel_widths = c(1,2))
    print(plot_p1b,vp=viewport(layout.pos.row=1,layout.pos.col=1))
    dev.off()

    longFile <- paste(strsplit(outFile,"\\.")[[1]][1],'-all.png',sep="")
    png(longFile,width=1200,height=800)
    p3       <- Plot_XY_graph(clusObj, nameX = "Transcription Factors", nameY = "Genes" ) +
      theme(axis.title.x=element_text(size=30, color = "black"),
            axis.title.y=element_text(size=30, color = "black"),
            legend.title = element_text(size=25),
            legend.text = element_text(size=25))
    plot_p1b <- plot_grid(p1b,p3, align= 'h', axis='b', rel_widths = c(1,2))
    print(plot_p1b,vp=viewport(layout.pos.row=1,layout.pos.col=1))
    dev.off()

  }

  return(meanSortIndex)
}
