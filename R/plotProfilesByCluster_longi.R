plotProfilesByCluster_longi<-function (riskProfObj, whichCovariates = NULL, rhoMinimum = NULL,
                                    useProfileStar = TRUE, covariate_info = list(title = "Covariate\ncategory",
                                                                                 levels = NULL, labels = NULL, split = NULL))
{
  profileDF <- tabulateCovariateProfiles_longi(riskProfObj = riskProfObj,
                                            whichCovariates = whichCovariates, rhoMinimum = rhoMinimum,
                                            useProfileStar = useProfileStar)
  if (!is.null(covariate_info$levels) & !is.null(covariate_info$labels)) {
    profileDF <- profileDF %>% dplyr::mutate(category = factor(category,
                                                               levels = covariate_info$levels, labels = covariate_info$labels))
  }
  empirical_proportions <- profileDF %>% dplyr::group_by(category,
                                                         covname) %>% dplyr::summarise(x = mean(mean)) %>% dplyr::group_by(covname) %>%
    dplyr::arrange(dplyr::desc(category)) %>% dplyr::mutate(emp_propn = cumsum(x)) %>%
    dplyr::filter(emp_propn < max(emp_propn)) %>% dplyr::ungroup() %>%
    dplyr::select(-x)

  varselect = TRUE
  if(varselect){
    covtab <- profileDF %>% dplyr::left_join(empirical_proportions,
                                             by = c("covname", "category")) %>% dplyr::group_by(cluster,
                                                                                                category, covname, fillColor, rhoMean, rhoRank, emp_propn) %>%
      dplyr::summarise(prop = mean(est)) %>% dplyr::ungroup()
  }else{
    covtab <- profileDF %>% dplyr::left_join(empirical_proportions,by = c("covname", "category")) %>%
      ##!! varSelect
      dplyr::group_by(cluster,category, covname, fillColor, emp_propn) %>% # , rhoMean, rhoRank
      dplyr::summarise(prop = mean(est)) %>% dplyr::ungroup()
  }

  if (!is.null(covariate_info$split)) {
    if (length(covariate_info$split) < 3) {
      covariate_info$split <- rep(covariate_info$split,
                                  3 - length(covariate_info$split))
      covariate_info$split <- c(covariate_info$split, sprintf("not %s",
                                                              covariate_info$split[2]))
    }
    covtab <- covtab %>% dplyr::mutate(type = stringr::str_detect(covname,
                                                                  covariate_info$split[1]), type = dplyr::recode(as.numeric(type),
                                                                                                                 `1` = covariate_info$split[2], `0` = covariate_info$split[3]))
    facetting_layer <- list(ggplot2::facet_grid(cluster ~
                                                  type, scales = "free_x", space = "free_x"))
  }
  else {
    facetting_layer <- list(ggplot2::facet_grid(cluster ~
                                                  .))
  }

  if(varselect){
    p2<- ggplot2::ggplot(covtab, ggplot2::aes(x = reorder(covname, rhoRank), y = prop, fill = factor(category), alpha = fillColor == "high"))
  }else{
    p2<- ggplot2::ggplot(covtab, ggplot2::aes(x = covname, ##!!reorder(covname, rhoRank),
                                              y = prop, fill = factor(category), alpha = fillColor == "high"))
  }
  p2 <- p2 +
    ggplot2::geom_bar(position = "fill", stat = "identity") +
    ggplot2::geom_point(ggplot2::aes(y = emp_propn, group = category), col = "black", fill = "white", alpha = 1, shape = 18) +
    #ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90,hjust = 1, vjust = 0.5)) +
    ggplot2::labs(x = "Covariates", y = "Cluster-specific proportion"#, title = "Covariate profiles",
                  #subtitle = "Black pips: empirical proportions in each covariate category.\n\n
                  #Dark fill: category more prevalent within cluster than overall."
    ) +
    theme(axis.title.x=element_text(size=30), axis.title.y=element_text(size=30))+
    theme(axis.text.x=element_text(size=20)) + theme(axis.text.y=element_text(size=20))+
    theme(legend.text=element_text(size=25),legend.title=element_text(size=25))+
    ggplot2::scale_fill_discrete(name = covariate_info$title) +
    ggplot2::scale_alpha_discrete(name="Significant", labels=c("no","yes"),range = c(0.25, 1)) +
    guides(
      fill = guide_legend(order = 1),
      alpha = guide_legend(order = 2)
    )+ facetting_layer
  return(p2)
}
