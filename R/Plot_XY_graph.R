Plot_XY_graph <- function (clusObj, nameX = "X variables", nameY = "Y")
{
  nSbj          <- clusObj$clusObjRunInfoObj$nSubjects
  clusters_orig <- clusObj$clustering
  ID_orig       <- 1:nSbj
  
  #Reorder by clusters
  ID       <- rev(order(clusters_orig))
  clusters <- clusters_orig[rev(order(clusters_orig))]
  Matrix   <- clusObj$clusObjRunInfoObj$xMat[ID,]
  first_indiv_percluster <- sapply(unique(clusters),function(x) which(clusters==x)[1])[-1]
  transc   <- c(as.matrix(Matrix))
  N        <- dim(Matrix)[1]
  p        <- dim(Matrix)[2]
  
  #Create ggplot dataframe
  data           <- data.frame("ID"=factor(rep(ID,p), levels=unique(ID)))
  data[,2]       <- as.character(rep(names(Matrix),each=N))
  data[,2]       <- factor(data[,2], levels=unique(data[,2]))
  names(data)[2] <- "varX"
  data$fill      <- as.factor(transc)
  data$Clusters     <- c(as.numeric(as.character(rep(clusters,p))))#as.factor(as.character(rep(clusters, p)))

  p       <- ggplot(data, aes(varX, ID)) + geom_tile(aes(fill = as.factor(Clusters), alpha = fill)) +
    xlab(nameX) + ylab(nameY) +
    scale_alpha_manual(values = c(0, 1), name = "Expressed", labels = c("No", "Yes")) + 
    scale_fill_discrete(name = "Clusters", labels = c("4", "3", "2", "1")) + 
    theme(axis.text.x = element_text(angle=90, hjust = 1, size = 10),
          axis.text.y=element_blank(), 
          axis.ticks.y=element_blank() #,legend.spacing.y = unit(0, 'cm')
          ) +
    guides(alpha = guide_legend(override.aes = list(color = "black"))) +
    scale_y_discrete(position="right") 
  
  #p <- p + geom_hline(yintercept = first_indiv_percluster, linetype = "dashed", color = "grey", size = 0.5)
   
  return(p) 
}
