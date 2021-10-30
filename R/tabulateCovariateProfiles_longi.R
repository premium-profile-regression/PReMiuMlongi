tabulateCovariateProfiles_longi<-function (riskProfObj, whichCovariates = NULL, rhoMinimum = NULL,
          useProfileStar = TRUE)
{

  for (i in 1:length(riskProfObj)) assign(names(riskProfObj)[i],
                                          riskProfObj[[i]])
  for (i in 1:length(riskProfClusObj)) assign(names(riskProfClusObj)[i],
                                              riskProfClusObj[[i]])
  for (i in 1:length(clusObjRunInfoObj)) assign(names(clusObjRunInfoObj)[i],
                                                clusObjRunInfoObj[[i]])
  xModel == "Discrete" || stop("Only implemented for discrete covariates.\n\n            Use plotRiskProfile() for Normal- or Mixed-covariate models.")
  if (useProfileStar) {
    profile <- profileStar
  }

  if(varSelect){
    rhotab <- tabulateVarSelectRho(riskProfObj)
    names(rhotab)[1]<-"var"
    if (is.null(whichCovariates) & !is.null(rhoMinimum)) {
      rhotab <- dplyr::filter(rhotab, rhoMean >= rhoMinimum)
      whichCovariates <- rhotab$var
    }
    if (length(whichCovariates) == 1) {
      rhotab <- dplyr::arrange(rhotab, rhoRank)
      whichCovariates <- rhotab$var[1:whichCovariates]
    }
  }

  if (!is.null(whichCovariates)) {
    if (!is.numeric(whichCovariates)) {
      whichCovariates <- match(whichCovariates, covNames)
    }
    covNames <- covNames[whichCovariates]
    nCovariates <- length(whichCovariates)
    profile <- profile[, , whichCovariates, ]
    nCategories <- nCategories[whichCovariates]
  }
  orderStat <- apply(risk, 2, median)
  meanSortIndex <- order(clusterSizes, orderStat, decreasing = T)
  clusterSizes <- clusterSizes[meanSortIndex]
  profile <- profile[, meanSortIndex, , ]
  profDFlist = list()

  for (j in 1:nCovariates) {
    miniDFlist = list()
    for (k in 1:nCategories[j]) {
      probMat <- profile[, , j, k]
      probMeans <- apply(probMat, 2, mean)
      probMean <- sum(probMeans * clusterSizes)/sum(clusterSizes)
      probLower <- apply(probMat, 2, quantile, 0.05)
      probUpper <- apply(probMat, 2, quantile, 0.95)
      clusterDF <- tibble::tibble(cluster = 1:nClusters,
                                  category = k - 1, mean = probMean, lower = probLower,
                                  upper = probUpper)
      profileDF <- tibble::tibble(cluster = rep(1:nClusters,
                                                each = nrow(probMat)), est = c(probMat)) %>%
        dplyr::left_join(clusterDF, by = "cluster")
      miniDFlist[[k]] <- profileDF
    }
    profileDF <- dplyr::bind_rows(miniDFlist)
    profDFlist[[covNames[j]]] <- profileDF
  }

  varselect=TRUE
  if(varselect){
    dplyr::bind_rows(profDFlist, .id = "covname") %>% dplyr::left_join(rhotab,
                                                                       by = c(covname = "var")) %>% dplyr::mutate(fillColor = ifelse(lower >
                                                                                                                                       mean, "high", ifelse(upper < mean, "low", "avg")), fillColor = as.character(fillColor))

  }else{
    dplyr::bind_rows(profDFlist, .id = "covname") %>% #dplyr::left_join(rhotab,by = c(covname = "var")) %>%
      dplyr::mutate(fillColor = ifelse(lower >mean, "high", ifelse(upper < mean, "low", "avg")), fillColor = as.character(fillColor))
  }
 }
